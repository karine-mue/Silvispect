"""Pure standard-library GeoTIFF reader, writer, and ``.tfw`` world-file parser.

Implements SPEC-GIS.md sections 2, 3, and 4 using only the Python standard
library (``struct``, ``zlib``, ``math``, ``pathlib``, ``re``).
"""

from __future__ import annotations

import math
import re
import struct
import zlib
from pathlib import Path
from typing import Any

from silvispect.grid import Grid, GridError

MAX_CELLS = 50_000_000

_NODATA_DEFAULT = -9999.0
_MAX_EXACT_INT = 1 << 53

# TIFF field types
_BYTE = 1
_ASCII = 2
_SHORT = 3
_LONG = 4
_RATIONAL = 5
_SBYTE = 6
_UNDEFINED = 7
_SSHORT = 8
_SLONG = 9
_SRATIONAL = 10
_FLOAT = 11
_DOUBLE = 12

_TYPE_SIZE: dict[int, int] = {
    _BYTE: 1,
    _ASCII: 1,
    _SHORT: 2,
    _LONG: 4,
    _RATIONAL: 8,
    _SBYTE: 1,
    _UNDEFINED: 1,
    _SSHORT: 2,
    _SLONG: 4,
    _SRATIONAL: 8,
    _FLOAT: 4,
    _DOUBLE: 8,
}

_TYPE_CODE: dict[int, str] = {
    _BYTE: "B",
    _SHORT: "H",
    _LONG: "I",
    _SBYTE: "b",
    _SSHORT: "h",
    _SLONG: "i",
    _FLOAT: "f",
    _DOUBLE: "d",
}

# TIFF tag IDs
_NEW_SUBFILE_TYPE = 254
_IMAGE_WIDTH = 256
_IMAGE_LENGTH = 257
_BITS_PER_SAMPLE = 258
_COMPRESSION = 259
_PHOTOMETRIC = 262
_FILL_ORDER = 266
_STRIP_OFFSETS = 273
_ORIENTATION = 274
_SAMPLES_PER_PIXEL = 277
_ROWS_PER_STRIP = 278
_STRIP_BYTE_COUNTS = 279
_PLANAR_CONFIG = 284
_PREDICTOR = 317
_TILE_WIDTH = 322
_TILE_LENGTH = 323
_TILE_OFFSETS = 324
_TILE_BYTE_COUNTS = 325
_SAMPLE_FORMAT = 339
_MODEL_PIXEL_SCALE = 33550
_MODEL_TIEPOINT = 33922
_MODEL_TRANSFORMATION = 34264
_GEO_KEY_DIRECTORY = 34735
_GDAL_METADATA = 42112
_GDAL_NODATA = 42113

_VALID_SAMPLE_FORMATS: dict[int, set[int]] = {
    1: {8, 16, 32, 64},
    2: {8, 16, 32, 64},
    3: {32, 64},
}

_SAMPLE_STRUCT_CODE: dict[tuple[int, int], str] = {
    (1, 8): "B",
    (1, 16): "H",
    (1, 32): "I",
    (1, 64): "Q",
    (2, 8): "b",
    (2, 16): "h",
    (2, 32): "i",
    (2, 64): "q",
    (3, 32): "f",
    (3, 64): "d",
}

_GEOGRAPHIC_EPSG_CODES: frozenset[int] = frozenset(
    {4326, 4612, 4301, 4269, 4258, 6666, 6667, 6668}
)

_GDAL_METADATA_ITEM_RE = re.compile(
    r'<Item\b[^>]*\bname="(SCALE|OFFSET)"[^>]*>([^<]*)</Item>'
)
_INT_LITERAL_RE = re.compile(r"[+-]?[0-9]+")
_CRS_RE = re.compile(r"EPSG:([1-9][0-9]*)")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def read_geotiff(path: Path | str) -> Grid:
    """Read a GeoTIFF (with optional ``.tfw`` fallback) into a :class:`Grid`."""
    file_path = Path(path)
    data = file_path.read_bytes()
    try:
        return _read_geotiff_bytes(data, file_path)
    except (GridError, OSError):
        raise
    except (
        struct.error,
        UnicodeDecodeError,
        OverflowError,
        zlib.error,
        IndexError,
        ValueError,
    ) as exc:
        raise GridError(f"corrupt TIFF: {exc}") from exc


def write_geotiff(grid: Grid, path: Path | str) -> None:
    """Write a :class:`Grid` to an uncompressed classic GeoTIFF."""
    try:
        data = _encode_geotiff(grid)
    except GridError:
        raise
    except (
        struct.error,
        UnicodeDecodeError,
        OverflowError,
        zlib.error,
        IndexError,
        ValueError,
    ) as exc:
        raise GridError(f"cannot write GeoTIFF: {exc}") from exc

    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_bytes(data)


# ---------------------------------------------------------------------------
# Reader implementation (SPEC-GIS.md §2, §3)
# ---------------------------------------------------------------------------


def _read_geotiff_bytes(data: bytes, file_path: Path) -> Grid:
    bo, tags = _parse_first_ifd(data)

    # NewSubfileType (254)
    if _NEW_SUBFILE_TYPE in tags:
        subfile = _get_single_int(tags, _NEW_SUBFILE_TYPE, bo, (_SHORT, _LONG))
        if subfile & 0x5:  # bit 0 (reduced-resolution) or bit 2 (mask)
            raise GridError(f"unsupported TIFF feature: NewSubfileType {subfile}")

    # Required dimensions
    width = _get_single_int(tags, _IMAGE_WIDTH, bo, (_SHORT, _LONG))
    height = _get_single_int(tags, _IMAGE_LENGTH, bo, (_SHORT, _LONG))
    if width < 1 or height < 1:
        raise GridError(f"corrupt TIFF: invalid image dimensions {width}x{height}")

    # SamplesPerPixel (277)
    samples_per_pixel = (
        _get_single_int(tags, _SAMPLES_PER_PIXEL, bo, (_SHORT,))
        if _SAMPLES_PER_PIXEL in tags
        else 1
    )
    if samples_per_pixel != 1:
        raise GridError(
            f"unsupported TIFF feature: multi-band raster ({samples_per_pixel} bands)"
        )

    # BitsPerSample (258) - required
    bps_vals = _get_ints(tags, _BITS_PER_SAMPLE, bo, (_SHORT,))
    if len(bps_vals) != samples_per_pixel:
        raise GridError("corrupt TIFF: BitsPerSample count does not match SamplesPerPixel")
    bits = bps_vals[0]

    # SampleFormat (339) - default 1
    if _SAMPLE_FORMAT in tags:
        sf_vals = _get_ints(tags, _SAMPLE_FORMAT, bo, (_SHORT,))
        if len(sf_vals) != samples_per_pixel:
            raise GridError(
                "corrupt TIFF: SampleFormat count does not match SamplesPerPixel"
            )
        sample_format = sf_vals[0]
    else:
        sample_format = 1

    if sample_format not in _VALID_SAMPLE_FORMATS or bits not in _VALID_SAMPLE_FORMATS[sample_format]:
        raise GridError(
            f"unsupported TIFF feature: sample format {sample_format} with {bits} bits"
        )

    # Compression (259) - default 1
    compression = (
        _get_single_int(tags, _COMPRESSION, bo, (_SHORT,))
        if _COMPRESSION in tags
        else 1
    )
    if compression not in (1, 32773, 8, 32946):
        if compression == 5:
            cname = "LZW"
        elif compression in (6, 7):
            cname = "JPEG"
        else:
            cname = "unknown"
        raise GridError(
            f"unsupported TIFF feature: compression {compression} ({cname})"
        )

    # PhotometricInterpretation (262) - optional
    if _PHOTOMETRIC in tags:
        photometric = _get_single_int(tags, _PHOTOMETRIC, bo, (_SHORT,))
        if photometric not in (0, 1):
            raise GridError(
                f"unsupported TIFF feature: photometric interpretation {photometric}"
            )

    # PlanarConfiguration (284) - default 1
    if _PLANAR_CONFIG in tags:
        planar = _get_single_int(tags, _PLANAR_CONFIG, bo, (_SHORT,))
        if planar not in (1, 2):
            raise GridError(
                f"unsupported TIFF feature: planar configuration {planar}"
            )

    # FillOrder (266) - default 1
    if _FILL_ORDER in tags:
        fill_order = _get_single_int(tags, _FILL_ORDER, bo, (_SHORT,))
        if fill_order != 1:
            raise GridError(f"unsupported TIFF feature: fill order {fill_order}")

    # Orientation (274) - default 1
    if _ORIENTATION in tags:
        orientation = _get_single_int(tags, _ORIENTATION, bo, (_SHORT,))
        if orientation != 1:
            raise GridError(f"unsupported TIFF feature: orientation {orientation}")

    # Predictor (317) - default 1
    predictor = (
        _get_single_int(tags, _PREDICTOR, bo, (_SHORT,))
        if _PREDICTOR in tags
        else 1
    )
    if predictor == 1:
        pass
    elif predictor == 2:
        if compression not in (8, 32946) or sample_format not in (1, 2):
            raise GridError(
                f"unsupported TIFF feature: predictor {predictor} with "
                f"compression {compression} and sample format {sample_format}"
            )
    elif predictor == 3:
        if compression not in (8, 32946) or sample_format != 3:
            raise GridError(
                f"unsupported TIFF feature: predictor {predictor} with "
                f"compression {compression} and sample format {sample_format}"
            )
    else:
        raise GridError(
            f"unsupported TIFF feature: predictor {predictor} with "
            f"compression {compression} and sample format {sample_format}"
        )

    # Raster size guard
    if width * height > MAX_CELLS:
        raise GridError(
            f"unsupported TIFF feature: raster too large ({width}x{height} > {MAX_CELLS})"
        )

    # GDAL_METADATA (42112) scale/offset check
    _check_gdal_metadata(tags)

    # GeoKeys & CRS (§2.7)
    geokeys = _parse_geokeys(tags, bo)
    crs = _resolve_crs(geokeys)

    # Georeferencing (§2.6 & §3)
    xllcorner, yllcorner, cellsize = _resolve_georeferencing(
        tags, bo, geokeys, height, file_path
    )

    # GDAL_NODATA (§2.8)
    grid_nodata, nodata_spec = _parse_nodata(tags, sample_format, bits)

    # Decode segments (§2.3, §2.4, §2.5)
    raw_samples = _read_samples(
        data, tags, bo, width, height, sample_format, bits, compression, predictor
    )
    values = _convert_samples(raw_samples, sample_format, nodata_spec)

    return Grid(
        ncols=width,
        nrows=height,
        xllcorner=xllcorner,
        yllcorner=yllcorner,
        cellsize=cellsize,
        nodata_value=grid_nodata,
        values=values,
        crs=crs,
    )


def _parse_first_ifd(data: bytes) -> tuple[str, dict[int, tuple[int, int, bytes]]]:
    if len(data) < 8:
        raise GridError("not a TIFF file: file shorter than 8 bytes")

    sig = data[:2]
    if sig == b"II":
        bo = "<"
    elif sig == b"MM":
        bo = ">"
    else:
        raise GridError("not a TIFF file: invalid byte order marker")

    (magic,) = struct.unpack_from(bo + "H", data, 2)
    if magic == 43:
        raise GridError("unsupported TIFF feature: BigTIFF")
    if magic != 42:
        raise GridError(f"not a TIFF file: unexpected magic number {magic}")

    (ifd_offset,) = struct.unpack_from(bo + "I", data, 4)
    if ifd_offset < 8 or ifd_offset + 2 > len(data):
        raise GridError("corrupt TIFF: IFD offset out of range")

    (num_entries,) = struct.unpack_from(bo + "H", data, ifd_offset)
    ifd_end = ifd_offset + 2 + num_entries * 12 + 4
    if ifd_end > len(data):
        raise GridError("corrupt TIFF: IFD extends beyond end of file")

    tags: dict[int, tuple[int, int, bytes]] = {}
    file_len = len(data)
    for i in range(num_entries):
        off = ifd_offset + 2 + i * 12
        tag_id, field_type, count = struct.unpack_from(bo + "HHI", data, off)
        if tag_id in tags:
            raise GridError(f"corrupt TIFF: duplicate tag {tag_id}")

        raw_field = data[off + 8 : off + 12]
        type_size = _TYPE_SIZE.get(field_type)
        if type_size is None:
            # Unknown field type: store raw 4-byte field and ignore unless required.
            tags[tag_id] = (field_type, count, raw_field)
            continue

        byte_len = count * type_size
        if byte_len <= 4:
            val_bytes = raw_field[:byte_len]
        else:
            (val_offset,) = struct.unpack(bo + "I", raw_field)
            if val_offset + byte_len > file_len:
                raise GridError(
                    f"corrupt TIFF: tag {tag_id} value out of range"
                )
            val_bytes = data[val_offset : val_offset + byte_len]
        tags[tag_id] = (field_type, count, val_bytes)

    return bo, tags


def _get_ints(
    tags: dict[int, tuple[int, int, bytes]],
    tag_id: int,
    bo: str,
    allowed_types: tuple[int, ...],
) -> tuple[int, ...]:
    if tag_id not in tags:
        raise GridError(f"corrupt TIFF: missing required tag {tag_id}")
    field_type, count, val_bytes = tags[tag_id]
    if field_type not in allowed_types:
        raise GridError(f"corrupt TIFF: unexpected type {field_type} for tag {tag_id}")
    if count == 0:
        return ()
    return struct.unpack(f"{bo}{count}{_TYPE_CODE[field_type]}", val_bytes)


def _get_single_int(
    tags: dict[int, tuple[int, int, bytes]],
    tag_id: int,
    bo: str,
    allowed_types: tuple[int, ...],
) -> int:
    vals = _get_ints(tags, tag_id, bo, allowed_types)
    if len(vals) != 1:
        raise GridError(f"corrupt TIFF: expected count 1 for tag {tag_id}, got {len(vals)}")
    return vals[0]


def _get_doubles(
    tags: dict[int, tuple[int, int, bytes]],
    tag_id: int,
    bo: str,
) -> tuple[float, ...]:
    if tag_id not in tags:
        raise GridError(f"corrupt TIFF: missing required tag {tag_id}")
    field_type, count, val_bytes = tags[tag_id]
    if field_type != _DOUBLE:
        raise GridError(f"corrupt TIFF: expected DOUBLE for tag {tag_id}")
    if count == 0:
        return ()
    return struct.unpack(f"{bo}{count}d", val_bytes)


def _check_gdal_metadata(tags: dict[int, tuple[int, int, bytes]]) -> None:
    if _GDAL_METADATA not in tags:
        return
    field_type, _count, val_bytes = tags[_GDAL_METADATA]
    if field_type != _ASCII:
        return
    try:
        text = val_bytes.rstrip(b"\x00").decode("ascii")
    except UnicodeDecodeError:
        return
    for match in _GDAL_METADATA_ITEM_RE.finditer(text):
        name = match.group(1)
        raw_val = match.group(2).strip()
        try:
            num = float(raw_val)
        except ValueError as exc:
            raise GridError(
                f"unsupported TIFF feature: scale/offset ({name}={raw_val!r})"
            ) from exc
        if name == "SCALE" and num != 1.0:
            raise GridError(f"unsupported TIFF feature: scale/offset (SCALE={num})")
        if name == "OFFSET" and num != 0.0:
            raise GridError(f"unsupported TIFF feature: scale/offset (OFFSET={num})")


def _parse_geokeys(
    tags: dict[int, tuple[int, int, bytes]], bo: str
) -> dict[int, int]:
    if _GEO_KEY_DIRECTORY not in tags:
        return {}
    vals = _get_ints(tags, _GEO_KEY_DIRECTORY, bo, (_SHORT,))
    if len(vals) < 4:
        raise GridError("corrupt TIFF: GeoKeyDirectoryTag shorter than 4 entries")
    version, _rev, _minor, num_keys = vals[0], vals[1], vals[2], vals[3]
    if version != 1:
        raise GridError(f"corrupt TIFF: unsupported GeoKeyDirectoryVersion {version}")
    if len(vals) < 4 + 4 * num_keys:
        raise GridError("corrupt TIFF: GeoKeyDirectoryTag truncated")

    referenced = {1024, 1025, 2048, 3072, 3076}
    seen: set[int] = set()
    result: dict[int, int] = {}
    for i in range(num_keys):
        base = 4 + 4 * i
        key_id, location, count, value_offset = (
            vals[base],
            vals[base + 1],
            vals[base + 2],
            vals[base + 3],
        )
        if key_id in seen:
            raise GridError(f"corrupt TIFF: duplicate GeoKey {key_id}")
        seen.add(key_id)
        if key_id in referenced:
            if location != 0 or count != 1:
                raise GridError(
                    f"corrupt TIFF: GeoKey {key_id} must have TIFFTagLocation=0 and Count=1"
                )
            result[key_id] = value_offset
    return result


def _resolve_crs(geokeys: dict[int, int]) -> str | None:
    model_type = geokeys.get(1024)
    if model_type in (2, 3):
        raise GridError(
            "unsupported CRS: geographic or geocentric coordinate system (degrees)"
        )
    if 3072 in geokeys:
        code = geokeys[3072]
        if code in (0, 32767):
            raise GridError(
                "unsupported CRS: user-defined or undefined projected CRS"
            )
        if code in _GEOGRAPHIC_EPSG_CODES:
            raise GridError(
                f"unsupported CRS: geographic EPSG code in ProjectedCSTypeGeoKey ({code})"
            )
        if code > 32766:
            # Grid.crs holds codes in 1..32766 (SPEC 5.1); refuse here, in
            # CRS terms, rather than letting the constructor's generic
            # message report it.
            raise GridError(f"unsupported CRS: EPSG code {code} is out of range")
        if 3076 in geokeys and geokeys[3076] != 9001:
            raise GridError(
                f"unsupported CRS: linear unit is not metre (code {geokeys[3076]})"
            )
        return f"EPSG:{code}"
    if 2048 in geokeys:
        raise GridError(
            "unsupported CRS: geographic or geocentric coordinate system (degrees)"
        )
    return None


def _resolve_georeferencing(
    tags: dict[int, tuple[int, int, bytes]],
    bo: str,
    geokeys: dict[int, int],
    height: int,
    file_path: Path,
) -> tuple[float, float, float]:
    has_scale = _MODEL_PIXEL_SCALE in tags
    has_tie = _MODEL_TIEPOINT in tags
    has_trans = _MODEL_TRANSFORMATION in tags

    if has_trans and (has_scale or has_tie):
        raise GridError(
            "unsupported georeferencing: ModelTransformationTag combined with scale/tiepoint"
        )

    if not has_trans and (has_scale != has_tie):
        raise GridError("corrupt TIFF: incomplete georeferencing tags")

    if not has_trans and not has_scale:
        return _read_tfw(file_path, height)

    # GeoTIFF tag path: validate GTRasterTypeGeoKey (1025)
    raster_type = geokeys.get(1025, 1)
    if raster_type not in (1, 2):
        raise GridError(f"unsupported georeferencing: raster type {raster_type}")

    if has_scale and has_tie:
        scale = _get_doubles(tags, _MODEL_PIXEL_SCALE, bo)
        tie = _get_doubles(tags, _MODEL_TIEPOINT, bo)
        if len(scale) != 3:
            raise GridError(
                f"corrupt TIFF: ModelPixelScaleTag must have 3 values, got {len(scale)}"
            )
        if len(tie) == 0 or len(tie) % 6 != 0:
            raise GridError(
                f"corrupt TIFF: ModelTiepointTag count {len(tie)} is not a positive multiple of 6"
            )
        if len(tie) > 6:
            raise GridError("unsupported georeferencing: more than one tiepoint")

        sx, sy, _sz = scale
        ti, tj, _tk, tx, ty, _tz = tie
        if not (math.isfinite(sx) and math.isfinite(sy)):
            raise GridError("unsupported georeferencing: non-finite pixel scale")
        if sy <= 0.0 or sx <= 0.0:
            raise GridError("unsupported georeferencing: non-north-up")
        if not math.isclose(sx, sy, rel_tol=1e-9):
            raise GridError("unsupported georeferencing: non-square cells")

        cellsize = sx
        x_ul = tx - ti * sx
        y_ul = ty + tj * sy
    else:
        m = _get_doubles(tags, _MODEL_TRANSFORMATION, bo)
        if len(m) != 16:
            raise GridError(
                f"corrupt TIFF: ModelTransformationTag must have 16 values, got {len(m)}"
            )
        if not all(math.isfinite(v) for v in m):
            raise GridError("unsupported georeferencing: non-finite transformation matrix")
        if m[1] != 0.0 or m[4] != 0.0:
            raise GridError("unsupported georeferencing: rotated or sheared raster")
        if m[12] != 0.0 or m[13] != 0.0 or m[14] != 0.0 or m[15] != 1.0:
            raise GridError("unsupported georeferencing: non-affine transformation")
        if m[0] <= 0.0 or m[5] >= 0.0:
            raise GridError("unsupported georeferencing: non-north-up")
        if not math.isclose(m[0], -m[5], rel_tol=1e-9):
            raise GridError("unsupported georeferencing: non-square cells")

        cellsize = m[0]
        x_ul = m[3]
        y_ul = m[7]

    if raster_type == 2:
        x_ul -= cellsize / 2.0
        y_ul += cellsize / 2.0

    xllcorner = x_ul
    yllcorner = y_ul - height * cellsize
    if not (
        math.isfinite(cellsize)
        and math.isfinite(x_ul)
        and math.isfinite(y_ul)
        and math.isfinite(yllcorner)
    ):
        raise GridError("unsupported georeferencing: non-finite origin or cellsize")
    return xllcorner, yllcorner, cellsize


def _read_tfw(tif_path: Path, height: int) -> tuple[float, float, float]:
    lower = tif_path.with_suffix(".tfw")
    upper = tif_path.with_suffix(".TFW")
    if lower.is_file():
        sidecar = lower
    elif upper.is_file():
        sidecar = upper
    else:
        raise GridError(
            "missing georeferencing: no GeoTIFF tags and no .tfw sidecar"
        )

    raw = sidecar.read_bytes()
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise GridError("corrupt world file: invalid UTF-8 encoding") from exc

    tokens = text.split()
    if len(tokens) != 6:
        raise GridError(f"corrupt world file: expected 6 values, got {len(tokens)}")

    nums: list[float] = []
    for tok in tokens:
        try:
            v = float(tok)
        except ValueError as exc:
            raise GridError(f"corrupt world file: invalid number {tok!r}") from exc
        if not math.isfinite(v):
            raise GridError(f"corrupt world file: non-finite value {tok!r}")
        nums.append(v)

    a, d, b, e, c, f = nums
    if b != 0.0 or d != 0.0:
        raise GridError("unsupported georeferencing: rotated or sheared world file")
    if e >= 0.0 or a <= 0.0:
        raise GridError("unsupported georeferencing: non-north-up world file")
    if not math.isclose(a, -e, rel_tol=1e-9):
        raise GridError("unsupported georeferencing: non-square cells")

    cellsize = a
    xllcorner = c - a / 2.0
    yllcorner = f - (height - 0.5) * cellsize
    if not (
        math.isfinite(cellsize)
        and math.isfinite(xllcorner)
        and math.isfinite(yllcorner)
    ):
        raise GridError("unsupported georeferencing: non-finite world file result")
    return xllcorner, yllcorner, cellsize


def _parse_nodata(
    tags: dict[int, tuple[int, int, bytes]], sample_format: int, bits: int
) -> tuple[float, tuple[str, Any]]:
    """Return ``(grid_nodata_value, match_spec)``."""
    if _GDAL_NODATA not in tags:
        return _NODATA_DEFAULT, ("finite", _NODATA_DEFAULT)

    field_type, count, val_bytes = tags[_GDAL_NODATA]
    if field_type != _ASCII or count < 1:
        raise GridError("corrupt TIFF: bad GDAL_NODATA tag type or length")

    try:
        raw_text = val_bytes.decode("ascii")
    except UnicodeDecodeError as exc:
        raise GridError("corrupt TIFF: bad GDAL_NODATA non-ASCII bytes") from exc

    text = raw_text.rstrip("\x00").strip()
    if not text:
        raise GridError("corrupt TIFF: bad GDAL_NODATA empty value")

    if sample_format in (1, 2) and _INT_LITERAL_RE.fullmatch(text):
        int_val = int(text)
        return float(int_val), ("finite", int_val)

    try:
        num = float(text)
    except ValueError as exc:
        raise GridError(f"corrupt TIFF: bad GDAL_NODATA {text!r}") from exc

    if math.isfinite(num):
        if sample_format == 3 and bits == 32:
            try:
                f32_val = struct.unpack("<f", struct.pack("<f", num))[0]
            except OverflowError:
                return num, ("none", None)
            return num, ("finite", f32_val)
        return num, ("finite", num)

    if math.isnan(num):
        return _NODATA_DEFAULT, ("nan", None)
    if num > 0:
        return _NODATA_DEFAULT, ("+inf", None)
    return _NODATA_DEFAULT, ("-inf", None)


def _read_samples(
    data: bytes,
    tags: dict[int, tuple[int, int, bytes]],
    bo: str,
    width: int,
    height: int,
    sample_format: int,
    bits: int,
    compression: int,
    predictor: int,
) -> list[int | float]:
    has_strips = any(
        t in tags for t in (_STRIP_OFFSETS, _ROWS_PER_STRIP, _STRIP_BYTE_COUNTS)
    )
    has_tiles = any(
        t in tags
        for t in (_TILE_WIDTH, _TILE_LENGTH, _TILE_OFFSETS, _TILE_BYTE_COUNTS)
    )
    if has_strips == has_tiles:
        raise GridError("corrupt TIFF: raster must define either strips or tiles")

    bps = bits // 8
    file_len = len(data)

    if has_strips:
        offsets = _get_ints(tags, _STRIP_OFFSETS, bo, (_SHORT, _LONG))
        byte_counts = _get_ints(tags, _STRIP_BYTE_COUNTS, bo, (_SHORT, _LONG))
        if _ROWS_PER_STRIP in tags:
            rps = _get_single_int(tags, _ROWS_PER_STRIP, bo, (_SHORT, _LONG))
            if rps == 0:
                raise GridError("corrupt TIFF: RowsPerStrip must be >= 1")
            if rps > height:
                rps = height
        else:
            rps = height

        num_strips = (height + rps - 1) // rps
        if len(offsets) != num_strips or len(byte_counts) != num_strips:
            raise GridError("corrupt TIFF: strip offset/bytecount count mismatch")

        out: list[int | float] = []
        for idx in range(num_strips):
            off = offsets[idx]
            cnt = byte_counts[idx]
            if off == 0 or cnt == 0:
                raise GridError("unsupported TIFF feature: sparse segment")
            if off + cnt > file_len:
                raise GridError("corrupt TIFF: segment out of range")

            rows_in_strip = (
                height - rps * (num_strips - 1) if idx == num_strips - 1 else rps
            )
            expected = rows_in_strip * width * bps
            raw = _decompress_segment(data[off : off + cnt], compression, expected)
            out.extend(
                _decode_block(
                    raw, width, rows_in_strip, sample_format, bits, bo, predictor
                )
            )
        return out

    # Tiles
    tw = _get_single_int(tags, _TILE_WIDTH, bo, (_SHORT, _LONG))
    th = _get_single_int(tags, _TILE_LENGTH, bo, (_SHORT, _LONG))
    if tw < 1 or th < 1:
        raise GridError("corrupt TIFF: tile dimensions must be >= 1")
    if tw * th > MAX_CELLS:
        raise GridError(
            f"unsupported TIFF feature: raster too large (tile {tw}x{th} > {MAX_CELLS})"
        )

    offsets = _get_ints(tags, _TILE_OFFSETS, bo, (_SHORT, _LONG))
    byte_counts = _get_ints(tags, _TILE_BYTE_COUNTS, bo, (_SHORT, _LONG))
    tiles_across = (width + tw - 1) // tw
    tiles_down = (height + th - 1) // th
    num_tiles = tiles_across * tiles_down
    if len(offsets) != num_tiles or len(byte_counts) != num_tiles:
        raise GridError("corrupt TIFF: tile offset/bytecount count mismatch")

    expected = th * tw * bps
    grid_samples: list[int | float] = [0] * (width * height)
    for ty in range(tiles_down):
        row_start = ty * th
        valid_rows = min(th, height - row_start)
        for tx in range(tiles_across):
            idx = ty * tiles_across + tx
            off = offsets[idx]
            cnt = byte_counts[idx]
            if off == 0 or cnt == 0:
                raise GridError("unsupported TIFF feature: sparse segment")
            if off + cnt > file_len:
                raise GridError("corrupt TIFF: segment out of range")

            raw = _decompress_segment(data[off : off + cnt], compression, expected)
            tile_vals = _decode_block(
                raw, tw, th, sample_format, bits, bo, predictor
            )
            col_start = tx * tw
            valid_cols = min(tw, width - col_start)
            for r in range(valid_rows):
                dst = (row_start + r) * width + col_start
                src = r * tw
                grid_samples[dst : dst + valid_cols] = tile_vals[
                    src : src + valid_cols
                ]
    return grid_samples


def _decompress_segment(seg: bytes, compression: int, expected: int) -> bytes:
    if compression == 1:
        if len(seg) < expected:
            raise GridError("corrupt TIFF: segment too short")
        return seg[:expected]

    if compression == 32773:
        return _packbits_decode(seg, expected)

    # Deflate (8 or 32946)
    decomp = zlib.decompressobj()
    raw = decomp.decompress(seg, expected)
    if len(raw) < expected:
        raw += decomp.flush(expected - len(raw))
    if len(raw) < expected:
        raise GridError("corrupt TIFF: segment too short")
    return raw[:expected]


def _packbits_decode(seg: bytes, expected: int) -> bytes:
    out = bytearray()
    i = 0
    n_seg = len(seg)
    while len(out) < expected:
        if i >= n_seg:
            raise GridError("corrupt TIFF: segment too short")
        control = seg[i]
        i += 1
        if control > 127:
            control -= 256
        if 0 <= control <= 127:
            count = control + 1
            if i + count > n_seg:
                raise GridError("corrupt TIFF: PackBits literal truncated")
            out.extend(seg[i : i + count])
            i += count
        elif control != -128:
            if i >= n_seg:
                raise GridError("corrupt TIFF: PackBits run truncated")
            out.extend(bytes([seg[i]]) * (1 - control))
            i += 1
    return bytes(out[:expected])


def _decode_block(
    raw: bytes,
    row_width: int,
    num_rows: int,
    sample_format: int,
    bits: int,
    bo: str,
    predictor: int,
) -> list[int | float]:
    if predictor == 1:
        code = _SAMPLE_STRUCT_CODE[(sample_format, bits)]
        return list(struct.unpack(f"{bo}{row_width * num_rows}{code}", raw))

    bps = bits // 8
    row_bytes = row_width * bps
    out: list[int | float] = []

    if predictor == 2:
        u_code = _SAMPLE_STRUCT_CODE[(1, bits)]
        mask = (1 << bits) - 1
        sign_bit = 1 << (bits - 1)
        modulus = 1 << bits
        is_signed = sample_format == 2
        row_fmt = f"{bo}{row_width}{u_code}"
        for r in range(num_rows):
            row = list(struct.unpack_from(row_fmt, raw, r * row_bytes))
            acc = 0
            for i in range(row_width):
                acc = (acc + row[i]) & mask
                if is_signed and acc >= sign_bit:
                    out.append(acc - modulus)
                else:
                    out.append(acc)
        return out

    # Predictor 3 (floating point)
    f_code = ">f" if bits == 32 else ">d"
    for r in range(num_rows):
        base = r * row_bytes
        row_buf = bytearray(row_bytes)
        acc = 0
        for i in range(row_bytes):
            acc = (acc + raw[base + i]) & 0xFF
            row_buf[i] = acc
        sample_buf = bytearray(bps)
        for i in range(row_width):
            for k in range(bps):
                sample_buf[k] = row_buf[k * row_width + i]
            out.append(struct.unpack(f_code, sample_buf)[0])
    return out


def _convert_samples(
    raw_samples: list[int | float],
    sample_format: int,
    nodata_spec: tuple[str, Any],
) -> list[float | None]:
    kind, nd_val = nodata_spec
    is_float = sample_format == 3
    nonfinite_nodata = kind in ("nan", "+inf", "-inf")
    values: list[float | None] = []

    for r in raw_samples:
        if kind == "finite":
            if r == nd_val:
                values.append(None)
                continue
        elif kind == "nan":
            if is_float and math.isnan(r):
                values.append(None)
                continue
        elif kind == "+inf":
            if is_float and r == float("inf"):
                values.append(None)
                continue
        elif kind == "-inf":
            if is_float and r == float("-inf"):
                values.append(None)
                continue

        if is_float:
            if not math.isfinite(r):
                raise GridError("corrupt TIFF: non-finite cell value")
            val = float(r)
        else:
            if abs(r) > _MAX_EXACT_INT:
                raise GridError(
                    "unsupported TIFF feature: integer sample not exactly representable as float"
                )
            val = float(r)

        if nonfinite_nodata and val == _NODATA_DEFAULT:
            raise GridError(
                "unsupported TIFF feature: nodata sentinel collision"
            )
        values.append(val)

    return values


# ---------------------------------------------------------------------------
# Writer implementation (SPEC-GIS.md §4)
# ---------------------------------------------------------------------------


def _is_exact_float32(v: float) -> bool:
    try:
        return struct.unpack("<f", struct.pack("<f", v))[0] == v
    except OverflowError:
        return False


def _validate_crs_for_write(crs: str | None) -> int | None:
    if crs is None:
        return None
    if not isinstance(crs, str):
        raise GridError(f"cannot write GeoTIFF: invalid CRS {crs!r}")
    match = _CRS_RE.fullmatch(crs)
    if match is None:
        raise GridError(f"cannot write GeoTIFF: invalid CRS {crs!r}")
    code = int(match.group(1))
    if not (1 <= code <= 32766):
        raise GridError(f"cannot write GeoTIFF: EPSG code out of range ({code})")
    return code


def _encode_geotiff(grid: Grid) -> bytes:
    ncols = grid.ncols
    nrows = grid.nrows
    if ncols * nrows > MAX_CELLS:
        raise GridError(
            f"cannot write GeoTIFF: raster too large ({ncols}x{nrows} > {MAX_CELLS})"
        )

    epsg_code = _validate_crs_for_write(grid.crs)

    y_ul = grid.yllcorner + nrows * grid.cellsize
    if not (
        math.isfinite(grid.xllcorner)
        and math.isfinite(grid.yllcorner)
        and math.isfinite(grid.cellsize)
        and math.isfinite(y_ul)
        and math.isfinite(grid.nodata_value)
    ):
        raise GridError("cannot write GeoTIFF: non-finite georeferencing or nodata")

    nd = float(grid.nodata_value)
    use_f32 = _is_exact_float32(nd)
    packed_values: list[float] = []
    for v in grid.values:
        if v is None:
            packed_values.append(nd)
        else:
            fv = float(v)
            if fv == nd:
                raise GridError(
                    f"cannot write GeoTIFF: cell value {fv!r} equals the nodata "
                    "value and would be read back as missing"
                )
            if use_f32 and not _is_exact_float32(fv):
                use_f32 = False
            packed_values.append(fv)

    bits = 32 if use_f32 else 64
    bps = bits // 8
    fmt_char = "f" if use_f32 else "d"
    strip_byte_count = ncols * nrows * bps

    nodata_bytes = f"{nd!r}\x00".encode("ascii")

    # Build IFD tag list in ascending order (SPEC-GIS.md §4.1)
    # Each entry: (tag_id, field_type, count, raw_value_bytes)
    entries: list[tuple[int, int, int, bytes]] = [
        (_IMAGE_WIDTH, _LONG, 1, struct.pack("<I", ncols)),
        (_IMAGE_LENGTH, _LONG, 1, struct.pack("<I", nrows)),
        (_BITS_PER_SAMPLE, _SHORT, 1, struct.pack("<H", bits)),
        (_COMPRESSION, _SHORT, 1, struct.pack("<H", 1)),
        (_PHOTOMETRIC, _SHORT, 1, struct.pack("<H", 1)),
        (_STRIP_OFFSETS, _LONG, 1, b"\x00\x00\x00\x00"),  # placeholder
        (_SAMPLES_PER_PIXEL, _SHORT, 1, struct.pack("<H", 1)),
        (_ROWS_PER_STRIP, _LONG, 1, struct.pack("<I", nrows)),
        (_STRIP_BYTE_COUNTS, _LONG, 1, struct.pack("<I", strip_byte_count)),
        (_PLANAR_CONFIG, _SHORT, 1, struct.pack("<H", 1)),
        (_SAMPLE_FORMAT, _SHORT, 1, struct.pack("<H", 3)),
        (
            _MODEL_PIXEL_SCALE,
            _DOUBLE,
            3,
            struct.pack("<3d", float(grid.cellsize), float(grid.cellsize), 0.0),
        ),
        (
            _MODEL_TIEPOINT,
            _DOUBLE,
            6,
            struct.pack(
                "<6d", 0.0, 0.0, 0.0, float(grid.xllcorner), float(y_ul), 0.0
            ),
        ),
    ]

    if epsg_code is not None:
        geokey_values = (
            1, 1, 0, 3,
            1024, 0, 1, 1,
            1025, 0, 1, 1,
            3072, 0, 1, epsg_code,
        )
        entries.append(
            (_GEO_KEY_DIRECTORY, _SHORT, 16, struct.pack("<16H", *geokey_values))
        )

    entries.append((_GDAL_NODATA, _ASCII, len(nodata_bytes), nodata_bytes))

    # Layout:
    # 1. Header (8 bytes)
    # 2. IFD at offset 8: 2 + 12 * N + 4 bytes
    # 3. External tag values in tag order, each on an even offset
    # 4. Pixel data on an even offset
    ifd_offset = 8
    ifd_size = 2 + 12 * len(entries) + 4
    ext_base = ifd_offset + ifd_size
    ext_buf = bytearray()
    serialized_entries: list[tuple[int, int, int, bytes]] = []

    for tag_id, field_type, count, val_bytes in entries:
        if len(val_bytes) <= 4:
            field_bytes = val_bytes.ljust(4, b"\x00")
        else:
            if (ext_base + len(ext_buf)) % 2:
                ext_buf.append(0)
            offset = ext_base + len(ext_buf)
            ext_buf.extend(val_bytes)
            field_bytes = struct.pack("<I", offset)
        serialized_entries.append((tag_id, field_type, count, field_bytes))

    if (ext_base + len(ext_buf)) % 2:
        ext_buf.append(0)
    pixel_offset = ext_base + len(ext_buf)

    # Patch StripOffsets entry with actual pixel_offset
    final_ifd = bytearray(struct.pack("<H", len(serialized_entries)))
    for tag_id, field_type, count, field_bytes in serialized_entries:
        if tag_id == _STRIP_OFFSETS:
            field_bytes = struct.pack("<I", pixel_offset)
        final_ifd.extend(struct.pack("<HHI", tag_id, field_type, count))
        final_ifd.extend(field_bytes)
    final_ifd.extend(struct.pack("<I", 0))  # next IFD offset = 0

    pixel_data = struct.pack(f"<{len(packed_values)}{fmt_char}", *packed_values)

    header = b"II" + struct.pack("<HI", 42, ifd_offset)
    return header + bytes(final_ifd) + bytes(ext_buf) + pixel_data