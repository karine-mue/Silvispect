"""Tests for the GeoTIFF reader (SPEC-GIS.md section 2 and the read-side of 6.2).

Covers byte order, strips and tiles, sample types, compression and predictors
(including the hand-computed vectors of appendix A), georeferencing (appendix B
numbers), GeoKeys and the CRS rules, nodata, and every read-side error
condition.  Files are made with ``tests/tiff_builder.py``, an independent
``struct``-only assembler; the reader is only reached through the public API
(``Grid.read`` and ``silvispect.cli.main``).

GDAL-made fixtures (SPEC 7.1) are replaced by the builder, as the SPEC allows;
the appendix A vectors are included as hand-written byte strings.
"""

from __future__ import annotations

import io
import math
import random
import struct
import zlib

import pytest

from silvispect.cli import EXIT_ERROR, EXIT_OK, main
from silvispect.grid import Grid, GridError

try:
    from tests import tiff_builder as tb
except ImportError:  # tests/ is not a package: its directory is on sys.path
    import tiff_builder as tb  # type: ignore[no-redef]

BYTEORDERS = ["<", ">"]

# --- helpers ----------------------------------------------------------------

# 1 m cells, upper-left corner at (1000, 2000): xll = 1000, yll = 2000 - nrows.
_GEO = [tb.scale_tag(1.0, 1.0, 0.0), tb.tiepoint_tag((0.0, 0.0, 0.0, 1000.0, 2000.0, 0.0))]


def _rows(height: int, width: int) -> list[list[float]]:
    return [[float(r * width + c + 1) for c in range(width)] for r in range(height)]


def _flat(rows) -> list:
    return [v for row in rows for v in row]


def _tif(rows, *, tags=(), **kwargs) -> bytes:
    """A TIFF with the default georeferencing; ``tags`` may override it."""
    return tb.build_tiff(rows, tags=[*_GEO, *tags], **kwargs)


def _write(tmp_path, data: bytes, name: str = "chm.tif"):
    path = tmp_path / name
    path.write_bytes(data)
    return path


def _read(tmp_path, data: bytes, name: str = "chm.tif") -> Grid:
    return Grid.read(_write(tmp_path, data, name))


def _rejects(tmp_path, data: bytes, match: str, name: str = "chm.tif") -> None:
    path = _write(tmp_path, data, name)
    with pytest.raises(GridError, match=match):
        Grid.read(path)


def _f32(value: float) -> float:
    return struct.unpack("<f", struct.pack("<f", value))[0]


_PROJ_KEYS = [(1024, 0, 1, 1), (1025, 0, 1, 1), (3072, 0, 1, 6677)]


# --- builder self-checks (appendix A defines the stored forms) ----------------


def test_builder_encodes_appendix_a_vectors():
    assert tb.predictor2_encode_row([10, 12, 9], 8) == [10, 2, 253]
    assert tb.predictor2_encode_row([100, 101, 99], 16) == [100, 1, 65534]
    assert tb.predictor3_encode_row([1.0, 2.0], 32) == bytes.fromhex("3F01408000000000")
    stream = bytes.fromhex("FEAA02010203" "80007F")
    assert tb.packbits_decode(stream, 6) == bytes.fromhex("AAAAAA010203")


# --- 2.9: result structure -----------------------------------------------------


def test_read_result_fields(tmp_path):
    rows = _rows(2, 3)
    grid = _read(tmp_path, _tif(rows))
    assert (grid.ncols, grid.nrows) == (3, 2)
    assert grid.xllcorner == 1000.0
    assert grid.yllcorner == 1998.0
    assert grid.cellsize == 1.0
    assert grid.nodata_value == -9999.0
    assert grid.crs is None
    assert grid.values == _flat(rows)
    assert all(isinstance(v, float) for v in grid.values)


def test_row_zero_is_the_north_row(tmp_path):
    grid = _read(tmp_path, _tif([[1.0, 2.0], [3.0, 4.0]]))
    assert grid.get(0, 0) == 1.0
    assert grid.get(1, 1) == 4.0
    assert grid.cell_center(0, 0) == (1000.5, 1999.5)


# --- 6.2 #1 / 1.2 / 2: extension handling and OS errors -------------------------


@pytest.mark.parametrize(
    "name", ["chm.txt", "chm.dat", "chm.gpkg", "chm.las", "chm.tif.bak", "chm", "dir/chm", ".tif", "chm."]
)
def test_unknown_or_missing_extension_is_refused_before_touching_the_file(tmp_path, name):
    with pytest.raises(GridError, match="unsupported raster file extension"):
        Grid.read(tmp_path / name)


def test_extension_message_names_the_suffix(tmp_path):
    with pytest.raises(GridError, match=r"unsupported raster file extension '\.txt'"):
        Grid.read(tmp_path / "a.txt")
    with pytest.raises(GridError, match=r"unsupported raster file extension ''"):
        Grid.read(tmp_path / "noext")


@pytest.mark.parametrize("name", ["CHM.TIF", "a.Tiff", "b.tiff", "c.TIFF"])
def test_tif_extension_is_case_insensitive(tmp_path, name):
    grid = _read(tmp_path, _tif(_rows(2, 2)), name)
    assert grid.values == _flat(_rows(2, 2))


def test_missing_tif_propagates_file_not_found(tmp_path):
    with pytest.raises(FileNotFoundError):
        Grid.read(tmp_path / "absent.tif")


def test_directory_named_like_a_tif_propagates_os_error(tmp_path):
    (tmp_path / "d.tif").mkdir()
    with pytest.raises(OSError):
        Grid.read(tmp_path / "d.tif")


# --- 2.1 / 6.2 #2, #3: container --------------------------------------------------


@pytest.mark.parametrize(
    "data",
    [b"", b"II", b"II*\x00\x08\x00\x00", b"XX*\x00\x08\x00\x00\x00", b"hello world, not a tiff"],
)
def test_bad_header_is_not_a_tiff(tmp_path, data):
    _rejects(tmp_path, data, "not a TIFF file")


@pytest.mark.parametrize("magic", [0, 41, 44, 0xFFFF])
@pytest.mark.parametrize("bo", BYTEORDERS)
def test_other_magic_is_not_a_tiff(tmp_path, bo, magic):
    _rejects(tmp_path, _tif(_rows(2, 2), byteorder=bo, magic=magic), "not a TIFF file")


@pytest.mark.parametrize("bo", BYTEORDERS)
def test_bigtiff_is_unsupported(tmp_path, bo):
    _rejects(tmp_path, _tif(_rows(2, 2), byteorder=bo, magic=43), "unsupported TIFF feature: BigTIFF")


@pytest.mark.parametrize("offset", [0, 4, 7, 10**6, 0xFFFFFFFF])
def test_ifd_offset_below_eight_or_beyond_the_file_is_corrupt(tmp_path, offset):
    _rejects(tmp_path, _tif(_rows(2, 2), ifd_offset=offset), "corrupt TIFF:")


def test_entry_count_beyond_the_file_is_corrupt(tmp_path):
    data = _tif(_rows(2, 2))
    data = tb.patch(data, tb.ifd_offset_of(data), b"\xff\xff")
    _rejects(tmp_path, data, "corrupt TIFF:")


def test_tag_value_offset_beyond_the_file_is_corrupt(tmp_path):
    data = tb.set_entry(_tif(_rows(2, 2)), tb.MODEL_TIEPOINT, field=10**7)
    _rejects(tmp_path, data, "corrupt TIFF:")


def test_huge_tag_count_is_refused_without_allocating(tmp_path):
    data = tb.set_entry(_tif(_rows(2, 2)), tb.MODEL_TIEPOINT, count=0xFFFFFFFF)
    _rejects(tmp_path, data, "corrupt TIFF:")


def test_duplicate_tag_is_corrupt(tmp_path):
    _rejects(
        tmp_path, _tif(_rows(2, 2), extra_tags=[tb.short(tb.PHOTOMETRIC, 1)]), "corrupt TIFF:"
    )


@pytest.mark.parametrize("bo", BYTEORDERS)
def test_tag_order_is_not_significant(tmp_path, bo):
    rows = _rows(2, 3)
    grid = _read(tmp_path, _tif(rows, byteorder=bo, tag_order="reversed"))
    assert grid.values == _flat(rows)
    assert grid.yllcorner == 1998.0


def test_unknown_tag_and_unknown_type_are_ignored(tmp_path):
    rows = _rows(2, 2)
    extra = [tb.short(60000, 7), tb.Tag(60001, 99, b"\x01\x02\x03\x04", count=1)]
    assert _read(tmp_path, _tif(rows, extra_tags=extra)).values == _flat(rows)


def test_only_the_first_ifd_is_read(tmp_path):
    rows = _rows(2, 2)
    second = [tb.long(tb.IMAGE_WIDTH, 99), tb.long(tb.IMAGE_LENGTH, 99), tb.long(tb.NEW_SUBFILE_TYPE, 1)]
    grid = _read(tmp_path, _tif(rows, second_ifd=second))
    assert (grid.ncols, grid.nrows) == (2, 2)
    assert grid.values == _flat(rows)


@pytest.mark.parametrize("subfile", [1, 4, 5])
def test_reduced_image_or_mask_first_ifd_is_refused(tmp_path, subfile):
    _rejects(
        tmp_path,
        _tif(_rows(2, 2), tags=[tb.long(tb.NEW_SUBFILE_TYPE, subfile)]),
        "unsupported TIFF feature: NewSubfileType",
    )


@pytest.mark.parametrize("subfile", [0, 2])
def test_other_subfile_type_bits_are_accepted(tmp_path, subfile):
    rows = _rows(2, 2)
    assert _read(tmp_path, _tif(rows, tags=[tb.long(tb.NEW_SUBFILE_TYPE, subfile)])).values == _flat(rows)


def test_every_proper_prefix_of_a_valid_file_is_refused(tmp_path):
    data = _tif(_rows(2, 2))
    path = tmp_path / "cut.tif"
    for length in range(len(data)):
        path.write_bytes(data[:length])
        with pytest.raises(GridError):
            Grid.read(path)


# --- 2.2: structure tags -------------------------------------------------------------


@pytest.mark.parametrize(
    "tag",
    [
        tb.ascii_tag(tb.IMAGE_WIDTH, "7"),
        tb.Tag(tb.IMAGE_WIDTH, tb.LONG, [2, 2]),
        tb.Tag(tb.IMAGE_LENGTH, tb.DOUBLE, [2.0]),
    ],
)
def test_image_size_tag_of_wrong_type_or_count_is_corrupt(tmp_path, tag):
    _rejects(tmp_path, _tif(_rows(2, 2), tags=[tag]), "corrupt TIFF:")


@pytest.mark.parametrize("missing", [tb.IMAGE_WIDTH, tb.IMAGE_LENGTH, tb.BITS_PER_SAMPLE])
def test_missing_required_tag_is_refused(tmp_path, missing):
    path = _write(tmp_path, _tif(_rows(2, 2), drop_tags=[missing]))
    with pytest.raises(GridError):
        Grid.read(path)


def test_zero_image_size_is_refused(tmp_path):
    path = _write(tmp_path, _tif(_rows(2, 2), tags=[tb.long(tb.IMAGE_WIDTH, 0)]))
    with pytest.raises(GridError):
        Grid.read(path)


def test_defaults_when_optional_tags_are_absent(tmp_path):
    # The data is uint-encoded, so the SPEC 2.2 defaults are what is being
    # tested: Compression 1 (none), SamplesPerPixel 1, SampleFormat 1 (uint).
    # (With float-encoded data and no SampleFormat tag, the default would
    # read the float bit patterns as integers, not the float values.)
    rows = [[1, 2], [3, 4]]
    data = _tif(
        rows,
        sample_format=1,
        photometric=None,
        drop_tags=[tb.COMPRESSION, tb.SAMPLE_FORMAT, tb.SAMPLES_PER_PIXEL],
    )
    assert _read(tmp_path, data).values == [1.0, 2.0, 3.0, 4.0]


def test_planar_configuration_two_is_accepted_for_a_single_band(tmp_path):
    rows = _rows(2, 2)
    assert _read(tmp_path, _tif(rows, tags=[tb.short(tb.PLANAR_CONFIG, 2)])).values == _flat(rows)


@pytest.mark.parametrize("photometric", [0, 1])
def test_photometric_zero_and_one_are_read_without_inversion(tmp_path, photometric):
    rows = _rows(2, 2)
    assert _read(tmp_path, _tif(rows, photometric=photometric)).values == _flat(rows)


@pytest.mark.parametrize("photometric", [2, 3, 5])
def test_other_photometric_interpretation_is_refused(tmp_path, photometric):
    _rejects(
        tmp_path,
        _tif(_rows(2, 2), photometric=photometric),
        "unsupported TIFF feature: photometric interpretation",
    )


def test_multi_band_is_refused(tmp_path):
    tags = [
        tb.short(tb.SAMPLES_PER_PIXEL, 3),
        tb.short(tb.BITS_PER_SAMPLE, 32, 32, 32),
        tb.short(tb.SAMPLE_FORMAT, 3, 3, 3),
    ]
    _rejects(tmp_path, _tif(_rows(2, 2), tags=tags), "unsupported TIFF feature: multi-band")


@pytest.mark.parametrize(
    "tag", [tb.short(tb.BITS_PER_SAMPLE, 32, 32), tb.short(tb.SAMPLE_FORMAT, 3, 3)]
)
def test_sample_tag_count_disagreeing_with_samples_per_pixel_is_corrupt(tmp_path, tag):
    _rejects(tmp_path, _tif(_rows(2, 2), tags=[tag]), "corrupt TIFF:")


@pytest.mark.parametrize(
    ("sample_format", "bits", "stored_as"),
    [
        (3, 16, (1, 16)),  # float16
        (1, 12, (1, 8)),
        (1, 4, (1, 8)),
        (1, 1, (1, 8)),
        (3, 8, (1, 8)),  # 8-bit float
        (4, 32, (3, 32)),  # SampleFormat 4 = undefined
        (5, 64, (3, 64)),  # complex
    ],
)
def test_unsupported_sample_format_combinations(tmp_path, sample_format, bits, stored_as):
    # Data is encoded with a valid layout; only the declared tags lie.
    fmt, width = stored_as
    rows = [[1.0, 2.0]] if fmt == 3 else [[1, 2]]
    tags = [tb.short(tb.BITS_PER_SAMPLE, bits), tb.short(tb.SAMPLE_FORMAT, sample_format)]
    data = _tif(rows, sample_format=fmt, bits=width, tags=tags)
    _rejects(tmp_path, data, "unsupported TIFF feature: sample format")


def test_fill_order_two_is_refused(tmp_path):
    _rejects(
        tmp_path, _tif(_rows(2, 2), tags=[tb.short(tb.FILL_ORDER, 2)]), "unsupported TIFF feature"
    )


def test_orientation_one_is_accepted_and_others_refused(tmp_path):
    rows = _rows(2, 2)
    assert _read(tmp_path, _tif(rows, tags=[tb.short(tb.ORIENTATION, 1)])).values == _flat(rows)
    for orientation in (2, 4, 8):
        _rejects(
            tmp_path,
            _tif(rows, tags=[tb.short(tb.ORIENTATION, orientation)]),
            "unsupported TIFF feature: orientation",
        )


def test_raster_larger_than_max_cells_is_refused_before_allocation(tmp_path):
    tags = [
        tb.long(tb.IMAGE_WIDTH, 10_000),
        tb.long(tb.IMAGE_LENGTH, 10_000),
        tb.long(tb.ROWS_PER_STRIP, 10_000),
    ]
    _rejects(tmp_path, _tif(_rows(2, 2), tags=tags), "unsupported TIFF feature: raster too large")


def _metadata(scale: str | None, offset: str | None) -> tb.Tag:
    items = ""
    if scale is not None:
        items += f'<Item name="SCALE" sample="0" role="scale">{scale}</Item>'
    if offset is not None:
        items += f'<Item name="OFFSET" sample="0" role="offset">{offset}</Item>'
    return tb.ascii_tag(tb.GDAL_METADATA, f"<GDALMetadata>{items}</GDALMetadata>")


@pytest.mark.parametrize(
    ("scale", "offset"), [("0.01", None), ("1", "5"), (None, "-1.5"), ("abc", None), ("1", "n/a")]
)
def test_gdal_metadata_scale_or_offset_is_refused(tmp_path, scale, offset):
    _rejects(
        tmp_path,
        _tif(_rows(2, 2), tags=[_metadata(scale, offset)]),
        "unsupported TIFF feature: scale/offset",
    )


@pytest.mark.parametrize(
    "tag",
    [
        _metadata("1", "0"),
        _metadata("1.0", "0.0"),
        tb.ascii_tag(tb.GDAL_METADATA, '<GDALMetadata><Item name="AREA_OR_POINT">Area</Item></GDALMetadata>'),
        tb.ascii_tag(tb.GDAL_METADATA, "not xml at all"),
        tb.Tag(tb.GDAL_METADATA, tb.ASCII, b"\xff\xfe\x00"),
    ],
)
def test_neutral_or_unreadable_gdal_metadata_is_ignored(tmp_path, tag):
    rows = _rows(2, 2)
    assert _read(tmp_path, _tif(rows, tags=[tag])).values == _flat(rows)


# --- 2.3: strips and tiles --------------------------------------------------------------


@pytest.mark.parametrize("bo", BYTEORDERS)
@pytest.mark.parametrize("code", [1, 32773, 8, 32946])
def test_multiple_strips_with_a_short_last_strip(tmp_path, bo, code):
    rows = _rows(5, 3)
    grid = _read(tmp_path, _tif(rows, byteorder=bo, compression=code, rows_per_strip=2))
    assert grid.values == _flat(rows)


@pytest.mark.parametrize("bo", BYTEORDERS)
@pytest.mark.parametrize("code", [1, 32773, 8, 32946])
def test_tiles_with_right_and_bottom_padding(tmp_path, bo, code):
    rows = _rows(3, 5)
    data = _tif(rows, byteorder=bo, compression=code, layout="tiles", tile_width=3, tile_length=2)
    grid = _read(tmp_path, data)
    assert (grid.ncols, grid.nrows) == (5, 3)
    assert grid.values == _flat(rows)


def test_tiles_numbered_row_major_across_then_down(tmp_path):
    rows = [[float(10 * r + c) for c in range(4)] for r in range(4)]
    data = _tif(rows, layout="tiles", tile_width=2, tile_length=2)
    assert _read(tmp_path, data).values == _flat(rows)


def test_rows_per_strip_larger_than_the_image_is_clamped(tmp_path):
    rows = _rows(3, 2)
    assert _read(tmp_path, _tif(rows, tags=[tb.long(tb.ROWS_PER_STRIP, 1000)])).values == _flat(rows)


def test_rows_per_strip_absent_means_one_strip(tmp_path):
    rows = _rows(3, 2)
    assert _read(tmp_path, _tif(rows, drop_tags=[tb.ROWS_PER_STRIP])).values == _flat(rows)


def test_rows_per_strip_zero_is_corrupt(tmp_path):
    _rejects(
        tmp_path, _tif(_rows(3, 2), tags=[tb.long(tb.ROWS_PER_STRIP, 0)]), "corrupt TIFF:"
    )


def test_strip_count_disagreeing_with_rows_per_strip_is_corrupt(tmp_path):
    data = _tif(_rows(4, 2), rows_per_strip=2, tags=[tb.long(tb.STRIP_OFFSETS, 8)])
    _rejects(tmp_path, data, "corrupt TIFF:")


def test_missing_strip_byte_counts_is_corrupt(tmp_path):
    _rejects(tmp_path, _tif(_rows(2, 2), drop_tags=[tb.STRIP_BYTE_COUNTS]), "corrupt TIFF:")


def test_strips_and_tiles_together_are_corrupt(tmp_path):
    tags = [
        tb.long(tb.TILE_WIDTH, 2),
        tb.long(tb.TILE_LENGTH, 2),
        tb.long(tb.TILE_OFFSETS, 8),
        tb.long(tb.TILE_BYTE_COUNTS, 16),
    ]
    _rejects(tmp_path, _tif(_rows(2, 2), tags=tags), "corrupt TIFF:")


def test_neither_strips_nor_tiles_is_corrupt(tmp_path):
    drop = [tb.STRIP_OFFSETS, tb.ROWS_PER_STRIP, tb.STRIP_BYTE_COUNTS]
    _rejects(tmp_path, _tif(_rows(2, 2), drop_tags=drop), "corrupt TIFF:")


def test_tile_count_disagreeing_with_the_grid_is_corrupt(tmp_path):
    data = _tif(_rows(4, 4), layout="tiles", tile_width=4, tile_length=4, tags=[tb.long(tb.TILE_OFFSETS, 8, 16)])
    _rejects(tmp_path, data, "corrupt TIFF:")


def test_tile_area_beyond_max_cells_is_refused(tmp_path):
    tags = [tb.long(tb.TILE_WIDTH, 10_000), tb.long(tb.TILE_LENGTH, 10_000)]
    data = _tif(_rows(4, 4), layout="tiles", tile_width=4, tile_length=4, tags=tags)
    _rejects(tmp_path, data, "unsupported TIFF feature: raster too large")


@pytest.mark.parametrize("tag", [tb.long(tb.TILE_WIDTH, 0), tb.long(tb.TILE_LENGTH, 0)])
def test_zero_tile_dimension_is_refused(tmp_path, tag):
    data = _tif(_rows(4, 4), layout="tiles", tile_width=4, tile_length=4, tags=[tag])
    with pytest.raises(GridError):
        Grid.read(_write(tmp_path, data))


def test_segment_beyond_the_file_is_corrupt(tmp_path):
    data = _tif(_rows(2, 2), tags=[tb.long(tb.STRIP_OFFSETS, 10_000_000)])
    _rejects(tmp_path, data, "corrupt TIFF: segment out of range")
    data = _tif(_rows(2, 2), tags=[tb.long(tb.STRIP_BYTE_COUNTS, 10_000_000)])
    _rejects(tmp_path, data, "corrupt TIFF: segment out of range")


@pytest.mark.parametrize(
    "tag", [tb.long(tb.STRIP_OFFSETS, 0), tb.long(tb.STRIP_BYTE_COUNTS, 0)]
)
def test_sparse_segment_is_unsupported(tmp_path, tag):
    _rejects(tmp_path, _tif(_rows(2, 2), tags=[tag]), "unsupported TIFF feature: sparse segment")


def test_uncompressed_segment_shorter_than_expected_is_corrupt(tmp_path):
    data = _tif(_rows(2, 2), tags=[tb.long(tb.STRIP_BYTE_COUNTS, 4)])
    _rejects(tmp_path, data, "corrupt TIFF:")


# --- 2.5: sample types ---------------------------------------------------------------------

SAMPLE_CASES = [
    (1, 8, [0, 255, 7, 100]),
    (1, 16, [0, 65535, 300, 1]),
    (1, 32, [0, 4294967295, 70000, 5]),
    (1, 64, [0, 2**53, 12345678901, 3]),
    (2, 8, [-128, 127, -1, 0]),
    (2, 16, [-32768, 32767, -2, 9]),
    (2, 32, [-(2**31), 2**31 - 1, -70000, 0]),
    (2, 64, [-(2**53), 2**53, -5, 0]),
    (3, 32, [1.5, -2.25, 0.0, 1.0e10]),
    (3, 64, [0.1, -1e300, 1e-300, 3.0]),
]


@pytest.mark.parametrize("bo", BYTEORDERS)
@pytest.mark.parametrize(("sample_format", "bits", "values"), SAMPLE_CASES)
def test_sample_types_read_back_as_floats(tmp_path, bo, sample_format, bits, values):
    rows = [values[:2], values[2:]]
    grid = _read(tmp_path, _tif(rows, byteorder=bo, sample_format=sample_format, bits=bits))
    assert grid.values == [float(v) for v in values]
    assert all(isinstance(v, float) for v in grid.values)


def test_float32_is_widened_exactly(tmp_path):
    value = _f32(0.1)
    grid = _read(tmp_path, _tif([[0.1]]))
    assert grid.values == [value]
    assert grid.values[0] != 0.1


def test_negative_zero_is_preserved(tmp_path):
    for bits in (32, 64):
        grid = _read(tmp_path, _tif([[-0.0, 1.0]], bits=bits))
        assert math.copysign(1.0, grid.values[0]) == -1.0


@pytest.mark.parametrize("bits", [32, 64])
@pytest.mark.parametrize("bad", [float("nan"), float("inf"), float("-inf")])
def test_non_finite_float_cell_without_matching_nodata_is_corrupt(tmp_path, bits, bad):
    _rejects(tmp_path, _tif([[1.0, bad]], bits=bits), "corrupt TIFF: non-finite cell value")


@pytest.mark.parametrize(("sample_format", "value"), [(1, 2**53 + 1), (1, 2**64 - 1), (2, -(2**53) - 1), (2, 2**62)])
def test_integer_sample_beyond_two_to_the_53_is_refused(tmp_path, sample_format, value):
    data = _tif([[1, value]], sample_format=sample_format, bits=64)
    _rejects(tmp_path, data, "unsupported TIFF feature: integer sample")


# --- 2.4: compression and predictor ---------------------------------------------------------


@pytest.mark.parametrize(("code", "name"), [(5, "LZW"), (6, "JPEG"), (7, "JPEG"), (50000, "unknown"), (99, "unknown")])
def test_unsupported_compression_is_named(tmp_path, code, name):
    data = _tif(_rows(2, 2), compression=code)
    _rejects(tmp_path, data, rf"unsupported TIFF feature: compression {code} \({name}\)")


def test_none_compression_uses_only_the_first_expected_bytes(tmp_path):
    rows = [[1.0, 2.0]]
    padded = struct.pack("<2f", 1.0, 2.0) + b"\xaa" * 6
    data = _tif(rows, raw_segments=[padded])
    assert _read(tmp_path, data).values == [1.0, 2.0]


def test_packbits_appendix_a_vector(tmp_path):
    # Appendix A: FE AA 02 01 02 03 80 00 7F -> AA AA AA 01 02 03 (7F unread).
    stream = bytes.fromhex("FEAA02010203" "80007F")
    data = _tif([[0] * 6], sample_format=1, bits=8, compression=32773, raw_segments=[stream])
    assert _read(tmp_path, data).values == [170.0, 170.0, 170.0, 1.0, 2.0, 3.0]
    data = _tif([[0] * 7], sample_format=1, bits=8, compression=32773, raw_segments=[stream])
    assert _read(tmp_path, data, "seven.tif").values == [170.0, 170.0, 170.0, 1.0, 2.0, 3.0, 127.0]


@pytest.mark.parametrize("bo", BYTEORDERS)
def test_predictor_2_appendix_a_vectors(tmp_path, bo):
    # uint8 row [10, 12, 9] is stored as [10, 2, 253] (SPEC appendix A,
    # as corrected 2026-11-19: 9 - 12 = -3 = 253 mod 256; 255 would
    # restore to [10, 12, 11]).
    stored = zlib.compress(bytes([10, 2, 253]))
    data = _tif([[0, 0, 0]], byteorder=bo, sample_format=1, bits=8, compression=8, predictor=2, raw_segments=[stored])
    assert _read(tmp_path, data).values == [10.0, 12.0, 9.0]
    # uint16 row [100, 101, 99] is stored as [100, 1, 65534] in the file byte order.
    stored = zlib.compress(struct.pack(bo + "3H", 100, 1, 65534))
    data = _tif([[0, 0, 0]], byteorder=bo, sample_format=1, bits=16, compression=8, predictor=2, raw_segments=[stored])
    assert _read(tmp_path, data, "u16.tif").values == [100.0, 101.0, 99.0]


@pytest.mark.parametrize("bo", BYTEORDERS)
@pytest.mark.parametrize("code", [8, 32946])
def test_predictor_3_appendix_a_vector(tmp_path, bo, code):
    # [1.0, 2.0] stored as 3F 01 40 80 00 00 00 00, whatever the file byte order.
    stored = zlib.compress(bytes.fromhex("3F01408000000000"))
    data = _tif([[0.0, 0.0]], byteorder=bo, compression=code, predictor=3, raw_segments=[stored])
    assert _read(tmp_path, data).values == [1.0, 2.0]


_STRIPS = {"rows_per_strip": 1}
_TILES = {"layout": "tiles", "tile_width": 2, "tile_length": 2}

INT_PREDICTOR_CASES = [
    (1, 8, [[10, 12, 9], [255, 0, 1]]),
    (1, 16, [[100, 101, 99], [65535, 0, 65535]]),
    (1, 32, [[4294967295, 0, 1], [7, 7, 7]]),
    (1, 64, [[2**53, 0, 5], [3, 2, 1]]),
    (2, 8, [[127, -128, 0], [1, -1, 2]]),
    (2, 16, [[-1, 5, -32768], [32767, -32768, 0]]),
    (2, 32, [[2**31 - 1, -(2**31), 0], [1, 2, 3]]),
    (2, 64, [[2**53, -(2**53), 0], [1, -1, 2]]),
]


@pytest.mark.parametrize("bo", BYTEORDERS)
@pytest.mark.parametrize("code", [8, 32946])
@pytest.mark.parametrize("layout", [_STRIPS, _TILES], ids=["strips", "tiles"])
@pytest.mark.parametrize(("sample_format", "bits", "rows"), INT_PREDICTOR_CASES)
def test_predictor_2_restores_integer_rows(tmp_path, bo, code, layout, sample_format, bits, rows):
    data = _tif(rows, byteorder=bo, sample_format=sample_format, bits=bits, compression=code, predictor=2, **layout)
    assert _read(tmp_path, data).values == [float(v) for v in _flat(rows)]


FLOAT_PREDICTOR_CASES = [
    (32, [[1.0, 2.0, -3.5], [0.0, 1.0e10, -0.25]]),
    (64, [[0.1, -2.5e100, 3.0], [1e-300, 0.0, 7.25]]),
]


@pytest.mark.parametrize("bo", BYTEORDERS)
@pytest.mark.parametrize("code", [8, 32946])
@pytest.mark.parametrize("layout", [_STRIPS, _TILES], ids=["strips", "tiles"])
@pytest.mark.parametrize(("bits", "rows"), FLOAT_PREDICTOR_CASES)
def test_predictor_3_restores_float_rows(tmp_path, bo, code, layout, bits, rows):
    data = _tif(rows, byteorder=bo, sample_format=3, bits=bits, compression=code, predictor=3, **layout)
    assert _read(tmp_path, data).values == _flat(rows)


@pytest.mark.parametrize("bo", BYTEORDERS)
def test_predictor_1_deflate_float_roundtrip(tmp_path, bo):
    rows = _rows(3, 4)
    assert _read(tmp_path, _tif(rows, byteorder=bo, compression=8, predictor=1)).values == _flat(rows)


@pytest.mark.parametrize(
    ("code", "sample_format", "bits", "predictor"),
    [
        (1, 1, 8, 2),  # predictor with no compression
        (1, 3, 32, 3),
        (32773, 1, 8, 2),  # predictor with PackBits
        (8, 1, 16, 3),  # float predictor on integers
        (8, 3, 32, 2),  # integer predictor on floats
        (8, 3, 32, 4),  # unknown predictor
        (8, 3, 32, 0),
    ],
)
def test_predictor_combinations_that_are_refused(tmp_path, code, sample_format, bits, predictor):
    rows = [[1, 2]] if sample_format != 3 else [[1.0, 2.0]]
    data = _tif(rows, sample_format=sample_format, bits=bits, compression=code, tags=[tb.short(tb.PREDICTOR, predictor)])
    _rejects(tmp_path, data, "unsupported TIFF feature: predictor")


def test_deflate_output_is_capped_at_the_expected_size(tmp_path):
    bomb = zlib.compress(b"\x00" * 5_000_000)
    data = _tif([[0.0, 0.0]], compression=8, raw_segments=[bomb])
    assert _read(tmp_path, data).values == [0.0, 0.0]


def test_deflate_surplus_bytes_are_ignored(tmp_path):
    payload = struct.pack("<2f", 1.0, 2.0) + b"\x55" * 8
    data = _tif([[0.0, 0.0]], compression=8, raw_segments=[zlib.compress(payload)])
    assert _read(tmp_path, data).values == [1.0, 2.0]


@pytest.mark.parametrize(
    "segment",
    # "empty" is a well-formed deflate stream of zero bytes: a raw b""
    # segment would have bytecount 0, which SPEC 6.2 assigns to the
    # sparse-segment refusal, not to stream corruption.
    [b"\x78\x9cthis is not a zlib stream", zlib.compress(b""), zlib.compress(b"\x00" * 4)],
    ids=["garbage", "empty", "short-stream"],
)
@pytest.mark.parametrize("code", [8, 32946])
def test_bad_deflate_stream_is_corrupt(tmp_path, code, segment):
    data = _tif([[0.0, 0.0]], compression=code, raw_segments=[segment])
    _rejects(tmp_path, data, "corrupt TIFF:")


@pytest.mark.parametrize(
    "segment",
    # No "empty" (b"") case: bytecount 0 is the sparse-segment refusal
    # (SPEC 6.2), and a PackBits stream yielding zero bytes is exactly
    # the "only-noop" case above.
    [bytes([0xFE]), bytes([0x05, 0x01, 0x02]), bytes([0x80]), bytes([0x00, 0x7F])],
    ids=["run-without-byte", "literal-cut-short", "only-noop", "too-few-bytes"],
)
def test_bad_packbits_stream_is_corrupt(tmp_path, segment):
    data = _tif([[0] * 4], sample_format=1, bits=8, compression=32773, raw_segments=[segment])
    _rejects(tmp_path, data, "corrupt TIFF:")


# --- 2.6 / appendix B: georeferencing ----------------------------------------------------------


def _geo(scale, tiepoint) -> list[tb.Tag]:
    return [tb.scale_tag(*scale), tb.tiepoint_tag(tiepoint)]


def _matrix(*, m0=0.5, m1=0.0, m3=612500.0, m4=0.0, m5=-0.5, m7=6642050.0, last=(0.0, 0.0, 0.0, 1.0)):
    return [m0, m1, 0.0, m3, m4, m5, 0.0, m7, 0.0, 0.0, 1.0, 0.0, *last]


def _tall(tags, **kwargs) -> bytes:
    """100 rows x 2 columns, so the appendix B numbers apply (nrows = 100)."""
    return tb.build_tiff(_rows(100, 2), tags=tags, **kwargs)


@pytest.mark.parametrize("bo", BYTEORDERS)
def test_appendix_b_tiepoint_pixel_is_area(tmp_path, bo):
    tags = _geo((0.5, 0.5, 0.0), (0, 0, 0, 612500.0, 6642050.0, 0))
    grid = _read(tmp_path, _tall(tags, byteorder=bo))
    assert (grid.xllcorner, grid.yllcorner, grid.cellsize) == (612500.0, 6642000.0, 0.5)
    assert grid.nrows == 100


def test_appendix_b_pixel_is_point_shifts_half_a_cell(tmp_path):
    tags = [*_geo((0.5, 0.5, 0.0), (0, 0, 0, 612500.25, 6642049.75, 0)), tb.geokey_tag((1025, 0, 1, 2))]
    grid = _read(tmp_path, _tall(tags))
    assert (grid.xllcorner, grid.yllcorner, grid.cellsize) == (612500.0, 6642000.0, 0.5)
    assert grid.crs is None


def test_appendix_b_tiepoint_with_nonzero_pixel_coordinates(tmp_path):
    tags = _geo((0.5, 0.5, 0.0), (2, 3, 0, 612501.0, 6642048.5, 0))
    grid = _read(tmp_path, _tall(tags))
    assert (grid.xllcorner, grid.yllcorner) == (612500.0, 6642000.0)


def test_appendix_b_world_file_example(tmp_path):
    path = _write(tmp_path, _tall([], drop_tags=[tb.MODEL_PIXEL_SCALE, tb.MODEL_TIEPOINT]))
    (tmp_path / "chm.tfw").write_text("0.5\n0\n0\n-0.5\n612500.25\n6642049.75\n", encoding="utf-8")
    grid = Grid.read(path)
    assert (grid.xllcorner, grid.yllcorner, grid.cellsize) == (612500.0, 6642000.0, 0.5)


def test_tiepoint_k_and_scale_z_are_ignored(tmp_path):
    tags = _geo((0.5, 0.5, 9.0), (0, 0, 5, 612500.0, 6642050.0, 123.0))
    grid = _read(tmp_path, _tall(tags))
    assert (grid.xllcorner, grid.yllcorner) == (612500.0, 6642000.0)


def test_scales_equal_within_relative_1e9_are_square(tmp_path):
    tags = _geo((0.5, 0.5 * (1 + 1e-10), 0.0), (0, 0, 0, 10.0, 20.0, 0))
    assert _read(tmp_path, _tall(tags)).cellsize == 0.5


def test_non_square_scale_is_refused(tmp_path):
    tags = _geo((0.5, 0.6, 0.0), (0, 0, 0, 10.0, 20.0, 0))
    _rejects(tmp_path, _tall(tags), "unsupported georeferencing: non-square cells")


@pytest.mark.parametrize("sy", [-0.5, 0.0])
def test_non_positive_y_scale_is_not_north_up(tmp_path, sy):
    tags = _geo((0.5, sy, 0.0), (0, 0, 0, 10.0, 20.0, 0))
    _rejects(tmp_path, _tall(tags), "unsupported georeferencing: non-north-up")


def test_nan_scale_is_refused(tmp_path):
    tags = _geo((float("nan"), float("nan"), 0.0), (0, 0, 0, 10.0, 20.0, 0))
    with pytest.raises(GridError):
        Grid.read(_write(tmp_path, _tall(tags)))


@pytest.mark.parametrize("bad", [float("inf"), float("nan"), 1e308])
def test_non_finite_derived_origin_is_refused(tmp_path, bad):
    tags = _geo((0.5, 0.5, 0.0), (0, 0, 0, bad, 20.0, 0))
    if bad == 1e308:  # finite tiepoint, but I * scale overflows
        tags = _geo((1e308, 1e308, 0.0), (1e308, 0, 0, 10.0, 20.0, 0))
        tags = _geo((1.7e308, 1.7e308, 0.0), (1e300, 1e300, 0, 10.0, 20.0, 0))
    _rejects(tmp_path, _tall(tags), "unsupported georeferencing")


def test_two_tiepoints_are_refused(tmp_path):
    tags = _geo((0.5, 0.5, 0.0), (0, 0, 0, 10.0, 20.0, 0, 1, 1, 0, 10.5, 19.5, 0))
    _rejects(tmp_path, _tall(tags), "unsupported georeferencing: more than one tiepoint")


@pytest.mark.parametrize(
    "tags",
    [
        [tb.scale_tag(0.5, 0.5), tb.double(tb.MODEL_TIEPOINT, 0, 0, 0, 1, 2, 0, 5)],  # count 7
        [tb.double(tb.MODEL_PIXEL_SCALE, 0.5, 0.5), tb.tiepoint_tag((0, 0, 0, 1, 2, 0))],  # count 2
    ],
)
def test_wrong_geotag_counts_are_corrupt(tmp_path, tags):
    _rejects(tmp_path, _tall(tags), "corrupt TIFF:")


@pytest.mark.parametrize("drop", [tb.MODEL_TIEPOINT, tb.MODEL_PIXEL_SCALE])
def test_one_of_scale_and_tiepoint_alone_is_corrupt_and_never_uses_the_sidecar(tmp_path, drop):
    data = _tif(_rows(2, 2), drop_tags=[drop])
    path = _write(tmp_path, data)
    (tmp_path / "chm.tfw").write_text("1\n0\n0\n-1\n1000.5\n1999.5\n", encoding="utf-8")
    with pytest.raises(GridError, match="corrupt TIFF: incomplete georeferencing tags"):
        Grid.read(path)


@pytest.mark.parametrize("bo", BYTEORDERS)
def test_model_transformation_is_accepted(tmp_path, bo):
    data = _tall([tb.transform_tag(_matrix())], byteorder=bo, drop_tags=[tb.MODEL_PIXEL_SCALE, tb.MODEL_TIEPOINT])
    grid = _read(tmp_path, data)
    assert (grid.xllcorner, grid.yllcorner, grid.cellsize) == (612500.0, 6642000.0, 0.5)


def test_model_transformation_pixel_is_point_shifts_half_a_cell(tmp_path):
    tags = [tb.transform_tag(_matrix(m3=612500.25, m7=6642049.75)), tb.geokey_tag((1025, 0, 1, 2))]
    grid = _read(tmp_path, _tall(tags, drop_tags=[tb.MODEL_PIXEL_SCALE, tb.MODEL_TIEPOINT]))
    assert (grid.xllcorner, grid.yllcorner) == (612500.0, 6642000.0)


@pytest.mark.parametrize("override", [{"m1": 0.1}, {"m4": -0.1}])
def test_rotated_or_sheared_transformation_is_refused(tmp_path, override):
    data = _tall([tb.transform_tag(_matrix(**override))], drop_tags=[tb.MODEL_PIXEL_SCALE, tb.MODEL_TIEPOINT])
    _rejects(tmp_path, data, "unsupported georeferencing: rotated or sheared raster")


@pytest.mark.parametrize(
    "override",
    [{"m0": 0.6}, {"m5": 0.5}, {"m0": -0.5}, {"last": (1.0, 0.0, 0.0, 1.0)}, {"last": (0.0, 0.0, 0.0, 2.0)}],
    ids=["non-square", "south-up", "negative-x", "perspective-row", "scale-row"],
)
def test_non_affine_or_non_north_up_transformation_is_refused(tmp_path, override):
    data = _tall([tb.transform_tag(_matrix(**override))], drop_tags=[tb.MODEL_PIXEL_SCALE, tb.MODEL_TIEPOINT])
    _rejects(tmp_path, data, "unsupported georeferencing:")


def test_transformation_with_wrong_count_is_refused(tmp_path):
    tag = tb.double(tb.MODEL_TRANSFORMATION, *_matrix()[:15])
    data = _tall([tag], drop_tags=[tb.MODEL_PIXEL_SCALE, tb.MODEL_TIEPOINT])
    with pytest.raises(GridError):
        Grid.read(_write(tmp_path, data))


@pytest.mark.parametrize(
    "keep",
    [(tb.MODEL_PIXEL_SCALE, tb.MODEL_TIEPOINT), (tb.MODEL_PIXEL_SCALE,), (tb.MODEL_TIEPOINT,)],
    ids=["both", "scale", "tiepoint"],
)
def test_transformation_combined_with_scale_or_tiepoint_is_refused(tmp_path, keep):
    drop = [t for t in (tb.MODEL_PIXEL_SCALE, tb.MODEL_TIEPOINT) if t not in keep]
    data = _tall([*_geo((0.5, 0.5, 0.0), (0, 0, 0, 1.0, 2.0, 0)), tb.transform_tag(_matrix())], drop_tags=drop)
    _rejects(tmp_path, data, "unsupported georeferencing: ModelTransformationTag combined")


# --- 3 (smoke): .tfw fallback conditions that belong to the read path -----------------------------


def test_no_georeferencing_at_all_is_missing(tmp_path):
    data = _tif(_rows(2, 2), drop_tags=[tb.MODEL_PIXEL_SCALE, tb.MODEL_TIEPOINT])
    _rejects(tmp_path, data, "missing georeferencing: no GeoTIFF tags and no .tfw sidecar")


def test_corrupt_world_file_is_reported(tmp_path):
    path = _write(tmp_path, _tif(_rows(2, 2), drop_tags=[tb.MODEL_PIXEL_SCALE, tb.MODEL_TIEPOINT]))
    (tmp_path / "chm.tfw").write_text("1\n0\n0\n-1\n1000.5\n", encoding="utf-8")
    with pytest.raises(GridError, match="corrupt world file"):
        Grid.read(path)


def test_rotated_world_file_is_refused(tmp_path):
    path = _write(tmp_path, _tif(_rows(2, 2), drop_tags=[tb.MODEL_PIXEL_SCALE, tb.MODEL_TIEPOINT]))
    (tmp_path / "chm.tfw").write_text("1\n0\n0.2\n-1\n1000.5\n1999.5\n", encoding="utf-8")
    with pytest.raises(GridError, match="unsupported georeferencing: rotated"):
        Grid.read(path)


def test_geotiff_tags_win_over_an_invalid_world_file(tmp_path):
    rows = _rows(2, 2)
    path = _write(tmp_path, _tif(rows))
    (tmp_path / "chm.tfw").write_text("garbage", encoding="utf-8")
    assert Grid.read(path).values == _flat(rows)


# --- 2.7: GeoKeyDirectory and CRS -----------------------------------------------------------------


def _crs(*keys, **kwargs) -> list[tb.Tag]:
    return [tb.geokey_tag(*keys, **kwargs)]


@pytest.mark.parametrize("bo", BYTEORDERS)
def test_projected_epsg_code_is_kept(tmp_path, bo):
    grid = _read(tmp_path, _tif(_rows(2, 2), byteorder=bo, tags=_crs(*_PROJ_KEYS)))
    assert grid.crs == "EPSG:6677"


@pytest.mark.parametrize("code", [32654, 6669, 2443, 3857])
def test_any_projected_code_is_formatted_as_epsg(tmp_path, code):
    keys = [(1024, 0, 1, 1), (3072, 0, 1, code)]
    assert _read(tmp_path, _tif(_rows(2, 2), tags=_crs(*keys))).crs == f"EPSG:{code}"


def test_projected_and_geographic_keys_together_are_accepted(tmp_path):
    keys = [*_PROJ_KEYS, (2048, 0, 1, 4612)]
    assert _read(tmp_path, _tif(_rows(2, 2), tags=_crs(*keys))).crs == "EPSG:6677"


def test_metre_unit_key_is_accepted(tmp_path):
    keys = [*_PROJ_KEYS, (3076, 0, 1, 9001)]
    assert _read(tmp_path, _tif(_rows(2, 2), tags=_crs(*keys))).crs == "EPSG:6677"


@pytest.mark.parametrize("unit", [9002, 9003, 9005, 0])
def test_non_metre_linear_unit_is_refused(tmp_path, unit):
    keys = [*_PROJ_KEYS, (3076, 0, 1, unit)]
    _rejects(tmp_path, _tif(_rows(2, 2), tags=_crs(*keys)), "unsupported CRS: linear unit")


@pytest.mark.parametrize("model_type", [2, 3])
def test_geographic_or_geocentric_model_type_is_refused(tmp_path, model_type):
    keys = [(1024, 0, 1, model_type), (2048, 0, 1, 4326)]
    _rejects(tmp_path, _tif(_rows(2, 2), tags=_crs(*keys)), "unsupported CRS: geographic")


def test_geographic_model_type_wins_over_a_projected_code(tmp_path):
    keys = [(1024, 0, 1, 2), (3072, 0, 1, 6677)]
    _rejects(tmp_path, _tif(_rows(2, 2), tags=_crs(*keys)), "unsupported CRS: geographic")


def test_geographic_type_alone_is_refused(tmp_path):
    keys = [(1025, 0, 1, 1), (2048, 0, 1, 4612)]
    _rejects(tmp_path, _tif(_rows(2, 2), tags=_crs(*keys)), "unsupported CRS: geographic")


@pytest.mark.parametrize("code", [0, 32767])
def test_user_defined_projected_crs_is_refused(tmp_path, code):
    keys = [(1024, 0, 1, 1), (3072, 0, 1, code)]
    _rejects(
        tmp_path,
        _tif(_rows(2, 2), tags=_crs(*keys)),
        "unsupported CRS: user-defined or undefined projected CRS",
    )


@pytest.mark.parametrize("code", [4326, 4612, 4301, 4269, 4258, 6666, 6667, 6668])
def test_geographic_epsg_in_projected_key_is_refused(tmp_path, code):
    keys = [(1024, 0, 1, 1), (3072, 0, 1, code)]
    _rejects(
        tmp_path,
        _tif(_rows(2, 2), tags=_crs(*keys)),
        "unsupported CRS: geographic EPSG code in ProjectedCSTypeGeoKey",
    )


def test_geokeys_without_crs_keys_give_unknown_crs(tmp_path):
    assert _read(tmp_path, _tif(_rows(2, 2), tags=_crs((1025, 0, 1, 1)))).crs is None
    assert _read(tmp_path, _tif(_rows(2, 2)), "nokeys.tif").crs is None


@pytest.mark.parametrize("raster_type", [0, 32767, 3])
def test_unknown_raster_type_is_refused(tmp_path, raster_type):
    _rejects(
        tmp_path,
        _tif(_rows(2, 2), tags=_crs((1025, 0, 1, raster_type))),
        "unsupported georeferencing: raster type",
    )


def test_unreferenced_keys_and_ascii_params_are_ignored(tmp_path):
    tags = [
        tb.geokey_tag(*_PROJ_KEYS, (1026, 34737, 5, 0)),
        tb.ascii_tag(tb.GEO_ASCII_PARAMS, "abcd|"),
    ]
    assert _read(tmp_path, _tif(_rows(2, 2), tags=tags)).crs == "EPSG:6677"


@pytest.mark.parametrize(
    "tag",
    [
        tb.geokey_tag(*_PROJ_KEYS, number_of_keys=5),  # directory shorter than announced
        tb.geokey_tag(*_PROJ_KEYS, version=(2, 1, 0)),
        tb.geokey_tag((1024, 0, 1, 1), (3072, 34736, 1, 0)),  # value in another tag
        tb.geokey_tag((1024, 0, 1, 1), (3072, 0, 2, 6677)),  # count != 1
        tb.geokey_tag((1024, 0, 1, 1), (1024, 0, 1, 1), (3072, 0, 1, 6677)),  # duplicate
        tb.geokey_tag((3072, 0, 1, 6677), (3072, 0, 1, 6678)),
    ],
    ids=["too-short", "bad-version", "other-location", "bad-count", "duplicate", "duplicate-crs"],
)
def test_corrupt_geokey_directory(tmp_path, tag):
    _rejects(tmp_path, _tif(_rows(2, 2), tags=[tag]), "corrupt TIFF:")


# --- 2.8: GDAL_NODATA ---------------------------------------------------------------------------------


def test_missing_nodata_tag_uses_the_default_sentinel(tmp_path):
    grid = _read(tmp_path, _tif([[1.0, -9999.0], [3.0, 4.0]]))
    assert grid.nodata_value == -9999.0
    assert grid.values == [1.0, None, 3.0, 4.0]


def test_finite_nodata_marks_matching_cells(tmp_path):
    grid = _read(tmp_path, _tif([[1.0, -32768.0]], tags=[tb.nodata_tag("-32768")]))
    assert grid.nodata_value == -32768.0
    assert grid.values == [1.0, None]


def test_default_sentinel_is_a_valid_value_when_another_nodata_is_declared(tmp_path):
    grid = _read(tmp_path, _tif([[-9999.0, 0.0]], tags=[tb.nodata_tag("0")]))
    assert grid.values == [-9999.0, None]


@pytest.mark.parametrize("text", ["-99", "-99\0", "  -99  "])
def test_short_and_padded_nodata_strings(tmp_path, text):
    tag = tb.ascii_tag(tb.GDAL_NODATA, text, nul=False)
    grid = _read(tmp_path, _tif([[-99.0, 1.0]], tags=[tag]))
    assert grid.nodata_value == -99.0
    assert grid.values == [None, 1.0]


def test_integer_nodata_on_integer_samples(tmp_path):
    grid = _read(tmp_path, _tif([[-32768, 5]], sample_format=2, bits=16, tags=[tb.nodata_tag("-32768")]))
    assert grid.values == [None, 5.0]


def test_int64_min_nodata_is_matched_exactly(tmp_path):
    data = _tif([[-(2**63), 5]], sample_format=2, bits=64, tags=[tb.nodata_tag("-9223372036854775808")])
    grid = _read(tmp_path, data)
    assert grid.values == [None, 5.0]
    assert grid.nodata_value == float(-(2**63))


@pytest.mark.parametrize("text", ["0.5", "-1", "300", "1e30"])
def test_integer_samples_with_unmatchable_nodata_have_no_missing_cells(tmp_path, text):
    grid = _read(tmp_path, _tif([[0, 1]], sample_format=1, bits=8, tags=[tb.nodata_tag(text)]))
    assert grid.values == [0.0, 1.0]
    assert grid.nodata_value == float(text)


def test_float32_nodata_is_compared_after_float32_rounding(tmp_path):
    cell = _f32(-9999.1234)
    assert cell != -9999.1234
    grid = _read(tmp_path, _tif([[-9999.1234, 3.5]], tags=[tb.nodata_tag("-9999.1234")]))
    assert grid.nodata_value == -9999.1234
    assert grid.values == [None, 3.5]
    assert cell not in [v for v in grid.values if v is not None]


def test_float32_nodata_outside_float32_range_matches_nothing(tmp_path):
    grid = _read(tmp_path, _tif([[1.0, 2.0]], tags=[tb.nodata_tag("1e39")]))
    assert grid.values == [1.0, 2.0]
    assert grid.nodata_value == 1e39


def test_float64_nodata_is_compared_exactly(tmp_path):
    rounded = _f32(-9999.1234)
    data = _tif([[-9999.1234, -9999.1235, rounded]], bits=64, tags=[tb.nodata_tag("-9999.1234")])
    assert _read(tmp_path, data).values == [None, -9999.1235, rounded]


def test_nan_nodata_marks_nan_cells_and_uses_the_default_sentinel(tmp_path):
    grid = _read(tmp_path, _tif([[float("nan"), 1.5]], tags=[tb.nodata_tag("nan")]))
    assert grid.values == [None, 1.5]
    assert grid.nodata_value == -9999.0


@pytest.mark.parametrize(("text", "marker"), [("inf", float("inf")), ("-inf", float("-inf"))])
def test_infinite_nodata_marks_same_signed_infinities(tmp_path, text, marker):
    grid = _read(tmp_path, _tif([[marker, 2.0]], bits=64, tags=[tb.nodata_tag(text)]))
    assert grid.values == [None, 2.0]
    assert grid.nodata_value == -9999.0


def test_infinite_nodata_does_not_cover_the_opposite_infinity(tmp_path):
    data = _tif([[float("-inf"), 2.0]], tags=[tb.nodata_tag("inf")])
    _rejects(tmp_path, data, "corrupt TIFF: non-finite cell value")


def test_non_finite_nodata_on_integer_samples_matches_nothing(tmp_path):
    grid = _read(tmp_path, _tif([[1, 2]], sample_format=1, bits=8, tags=[tb.nodata_tag("nan")]))
    assert grid.values == [1.0, 2.0]
    assert grid.nodata_value == -9999.0


def test_non_finite_nodata_with_a_default_sentinel_collision_is_refused(tmp_path):
    data = _tif([[float("nan"), -9999.0]], tags=[tb.nodata_tag("nan")])
    _rejects(tmp_path, data, "unsupported TIFF feature: nodata sentinel collision")


@pytest.mark.parametrize(
    "tag",
    [
        tb.Tag(tb.GDAL_NODATA, tb.SHORT, [1]),
        tb.Tag(tb.GDAL_NODATA, tb.BYTE, [45, 49, 0]),
        tb.ascii_tag(tb.GDAL_NODATA, "", nul=False),
        tb.ascii_tag(tb.GDAL_NODATA, ""),
        tb.ascii_tag(tb.GDAL_NODATA, "   "),
        tb.Tag(tb.GDAL_NODATA, tb.ASCII, b"\xff\xfe\x00"),
    ],
    ids=["short", "byte", "empty", "nul-only", "blank", "not-ascii"],
)
def test_malformed_nodata_tag_is_corrupt(tmp_path, tag):
    _rejects(tmp_path, _tif([[1.0, 2.0]], tags=[tag]), "corrupt TIFF:")


def test_unparseable_nodata_text_is_reported_as_bad_gdal_nodata(tmp_path):
    _rejects(tmp_path, _tif([[1.0]], tags=[tb.nodata_tag("abc")]), "corrupt TIFF: bad GDAL_NODATA")


# --- CLI: reading goes through the same dispatcher ------------------------------------------------------


def test_cli_inspect_reads_a_tif(tmp_path):
    path = _write(tmp_path, _tif(_rows(3, 4)))
    code = main(["inspect", "--chm", str(path)], out=io.StringIO())
    assert code == EXIT_OK


def test_cli_exits_with_error_for_unsupported_tiff(tmp_path):
    path = _write(tmp_path, _tif(_rows(3, 4), compression=5))
    code = main(["inspect", "--chm", str(path)], out=io.StringIO())
    assert code == EXIT_ERROR


# --- 6.4: no raw exception may escape (fixed-seed fuzzing) ------------------------------------------------


def _fuzz_bases() -> list[bytes]:
    rows = _rows(3, 4)
    return [
        _tif(rows),
        _tif(rows, byteorder=">", compression=8, predictor=3, rows_per_strip=2),
        _tif(
            [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
            sample_format=1,
            bits=16,
            compression=32773,
            layout="tiles",
            tile_width=2,
            tile_length=2,
        ),
        _tif(rows, tags=[tb.geokey_tag(*_PROJ_KEYS), tb.nodata_tag("-9999")]),
    ]


def test_random_bytes_never_leak_a_raw_exception(tmp_path):
    rng = random.Random(20240601)
    path = tmp_path / "fuzz.tif"
    prefixes = [b"", b"II*\x00\x08\x00\x00\x00", b"MM\x00*\x00\x00\x00\x08"]
    for _ in range(300):
        data = rng.choice(prefixes) + bytes(rng.randrange(256) for _ in range(rng.randrange(0, 200)))
        path.write_bytes(data)
        try:
            Grid.read(path)
        except GridError:
            pass


def test_mutated_valid_files_never_leak_a_raw_exception(tmp_path):
    rng = random.Random(4242)
    path = tmp_path / "mutant.tif"
    for base in _fuzz_bases():
        ifd_start, ifd_end = tb.ifd_span(base)
        for _ in range(150):
            data = bytearray(base)
            for _ in range(rng.randint(1, 6)):
                region = rng.choice(["head", "ifd", "any"])
                if region == "head":
                    pos = rng.randrange(0, 16)
                elif region == "ifd":
                    pos = rng.randrange(ifd_start, ifd_end)
                else:
                    pos = rng.randrange(0, len(data))
                data[pos] = rng.randrange(256)
            path.write_bytes(bytes(data))
            try:
                Grid.read(path)
            except GridError:
                pass
