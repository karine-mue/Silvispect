"""Test-only minimal TIFF assembler (independent of ``silvispect``).

Everything here is built from ``struct`` and ``zlib`` alone so that the
reader under test and the files it is tested with share no code, and therefore
no misunderstanding of the format.  The encodings follow SPEC-GIS.md section
2.4 and appendix A:

* Predictor 2: per row, ``d[0] = s[0]``, ``d[i] = (s[i] - s[i-1]) mod 2**bits``
  on unsigned samples in the file byte order.
* Predictor 3: per row, each sample is serialised big-endian IEEE 754, the
  bytes are regrouped into planes (most significant byte plane first), and
  the whole row is byte-differenced (``b[i] = (b[i] - b[i-1]) mod 256``).

Beyond well-formed files the builder can inject faults: overriding any tag,
adding duplicates, lying about counts, types and offsets, supplying raw
(pre-compressed) segments, and patching bytes of a finished file.
"""

from __future__ import annotations

import struct
import zlib
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

# --- TIFF field types -------------------------------------------------------
BYTE, ASCII, SHORT, LONG, RATIONAL = 1, 2, 3, 4, 5
SBYTE, UNDEFINED, SSHORT, SLONG, SRATIONAL, FLOAT, DOUBLE = 6, 7, 8, 9, 10, 11, 12

TYPE_SIZE = {
    BYTE: 1,
    ASCII: 1,
    SHORT: 2,
    LONG: 4,
    RATIONAL: 8,
    SBYTE: 1,
    UNDEFINED: 1,
    SSHORT: 2,
    SLONG: 4,
    SRATIONAL: 8,
    FLOAT: 4,
    DOUBLE: 8,
}
_TYPE_CODE = {
    BYTE: "B",
    SHORT: "H",
    LONG: "I",
    SBYTE: "b",
    SSHORT: "h",
    SLONG: "i",
    FLOAT: "f",
    DOUBLE: "d",
}

# --- tag ids ---------------------------------------------------------------
NEW_SUBFILE_TYPE = 254
IMAGE_WIDTH = 256
IMAGE_LENGTH = 257
BITS_PER_SAMPLE = 258
COMPRESSION = 259
PHOTOMETRIC = 262
FILL_ORDER = 266
STRIP_OFFSETS = 273
ORIENTATION = 274
SAMPLES_PER_PIXEL = 277
ROWS_PER_STRIP = 278
STRIP_BYTE_COUNTS = 279
PLANAR_CONFIG = 284
PREDICTOR = 317
TILE_WIDTH = 322
TILE_LENGTH = 323
TILE_OFFSETS = 324
TILE_BYTE_COUNTS = 325
SAMPLE_FORMAT = 339
MODEL_PIXEL_SCALE = 33550
MODEL_TIEPOINT = 33922
MODEL_TRANSFORMATION = 34264
GEO_KEY_DIRECTORY = 34735
GEO_DOUBLE_PARAMS = 34736
GEO_ASCII_PARAMS = 34737
GDAL_METADATA = 42112
GDAL_NODATA = 42113


@dataclass
class Tag:
    """One IFD entry.

    ``count`` overrides the count written into the entry (the data written is
    unaffected), and ``field`` overrides the 4-byte value/offset field: an
    ``int`` is written in the file byte order, ``bytes`` verbatim.  Giving
    ``field`` suppresses the external value data.  For a type that is not a
    standard TIFF type, ``values`` is written as raw bytes.
    """

    tag: int
    type: int
    values: Any = ()
    count: int | None = None
    field: int | bytes | None = None


def short(tag: int, *values: int) -> Tag:
    return Tag(tag, SHORT, list(values))


def long(tag: int, *values: int) -> Tag:
    return Tag(tag, LONG, list(values))


def double(tag: int, *values: float) -> Tag:
    return Tag(tag, DOUBLE, list(values))


def ascii_tag(tag: int, text: str | bytes, *, nul: bool = True) -> Tag:
    data = text.encode("ascii") if isinstance(text, str) else bytes(text)
    return Tag(tag, ASCII, data + (b"\0" if nul else b""))


def nodata_tag(text: str | bytes) -> Tag:
    return ascii_tag(GDAL_NODATA, text)


def scale_tag(sx: float, sy: float, sz: float = 0.0) -> Tag:
    return double(MODEL_PIXEL_SCALE, sx, sy, sz)


def tiepoint_tag(values: Sequence[float]) -> Tag:
    return double(MODEL_TIEPOINT, *values)


def transform_tag(matrix: Sequence[float]) -> Tag:
    return double(MODEL_TRANSFORMATION, *matrix)


def geokey_tag(
    *keys: tuple[int, int, int, int],
    version: tuple[int, int, int] = (1, 1, 0),
    number_of_keys: int | None = None,
) -> Tag:
    """GeoKeyDirectoryTag from ``(KeyID, TIFFTagLocation, Count, Value)`` tuples."""
    n = len(keys) if number_of_keys is None else number_of_keys
    flat = [version[0], version[1], version[2], n]
    for key in keys:
        flat.extend(key)
    return short(GEO_KEY_DIRECTORY, *flat)


# --- sample encoding --------------------------------------------------------
def sample_code(sample_format: int, bits: int) -> str:
    table = {
        1: {8: "B", 16: "H", 32: "I", 64: "Q"},
        2: {8: "b", 16: "h", 32: "i", 64: "q"},
        3: {32: "f", 64: "d"},
    }
    return table[sample_format][bits]


def predictor2_encode_row(samples: Sequence[int], bits: int) -> list[int]:
    """Horizontal differencing on unsigned samples (appendix A)."""
    mask = (1 << bits) - 1
    out: list[int] = []
    previous = 0
    for i, sample in enumerate(samples):
        current = sample & mask
        out.append(current if i == 0 else (current - previous) & mask)
        previous = current
    return out


def predictor3_encode_row(samples: Sequence[float], bits: int) -> bytes:
    """Floating-point predictor: big-endian byte planes, then byte deltas."""
    bps = bits // 8
    fmt = ">f" if bits == 32 else ">d"
    width = len(samples)
    chunks = [struct.pack(fmt, s) for s in samples]
    planes = bytes(chunks[i][k] for k in range(bps) for i in range(width))
    return bytes([planes[0]] + [(planes[i] - planes[i - 1]) & 0xFF for i in range(1, len(planes))])


def encode_row(
    samples: Sequence[float], sample_format: int, bits: int, byteorder: str, predictor: int
) -> bytes:
    if predictor == 3:
        return predictor3_encode_row(samples, bits)
    if predictor == 2:
        deltas = predictor2_encode_row(list(samples), bits)
        return struct.pack(f"{byteorder}{len(deltas)}{sample_code(1, bits)}", *deltas)
    code = sample_code(sample_format, bits)
    return struct.pack(f"{byteorder}{len(samples)}{code}", *samples)


# --- compression ------------------------------------------------------------
def packbits_encode(data: bytes) -> bytes:
    out = bytearray()
    i, n = 0, len(data)
    while i < n:
        run = 1
        while i + run < n and run < 128 and data[i + run] == data[i]:
            run += 1
        if run >= 2:
            out.append((257 - run) & 0xFF)  # control byte 1 - run, as unsigned
            out.append(data[i])
            i += run
        else:
            j = i
            while j < n and j - i < 128 and (j + 1 >= n or data[j] != data[j + 1]):
                j += 1
            out.append(j - i - 1)
            out += data[i:j]
            i = j
    return bytes(out)


def packbits_decode(data: bytes, expected: int) -> bytes:
    """Reference PackBits decoder used only to self-check the appendix A vector."""
    out = bytearray()
    i = 0
    while len(out) < expected and i < len(data):
        control = data[i]
        i += 1
        signed = control - 256 if control > 127 else control
        if 0 <= signed <= 127:
            out += data[i : i + signed + 1]
            i += signed + 1
        elif signed != -128:
            out += bytes([data[i]]) * (1 - signed)
            i += 1
    return bytes(out[:expected])


def compress(data: bytes, code: int) -> bytes:
    if code == 32773:
        return packbits_encode(data)
    if code in (8, 32946):
        return zlib.compress(data)
    return data  # none, and unsupported codes: stored as-is


# --- assembling -------------------------------------------------------------
def _encode_values(tag: Tag, bo: str) -> tuple[bytes, int]:
    t, v = tag.type, tag.values
    if t not in TYPE_SIZE:
        data = bytes(v)
        return data, len(data)
    if t in (ASCII, UNDEFINED):
        data = v.encode("ascii") if isinstance(v, str) else bytes(v)
        return data, len(data)
    values = list(v)
    if t in (RATIONAL, SRATIONAL):
        code = "I" if t == RATIONAL else "i"
        return struct.pack(f"{bo}{len(values)}{code}", *values), len(values) // 2
    return struct.pack(f"{bo}{len(values)}{_TYPE_CODE[t]}", *values), len(values)


def _serialize_ifd(entries: list[Tag], pos: int, bo: str, next_ifd: int) -> tuple[bytes, bytes]:
    n = len(entries)
    ext_pos = pos + 2 + 12 * n + 4
    head = bytearray(struct.pack(bo + "H", n))
    ext = bytearray()
    for entry in entries:
        data, count = _encode_values(entry, bo)
        if entry.count is not None:
            count = entry.count
        if entry.field is None:
            if len(data) <= 4:
                fld = data.ljust(4, b"\0")  # left-justified, whatever the byte order
            else:
                if (ext_pos + len(ext)) % 2:
                    ext.append(0)
                fld = struct.pack(bo + "I", ext_pos + len(ext))
                ext += data
        elif isinstance(entry.field, int):
            fld = struct.pack(bo + "I", entry.field)
        else:
            fld = bytes(entry.field).ljust(4, b"\0")[:4]
        head += struct.pack(bo + "HHI", entry.tag, entry.type, count) + fld
    head += struct.pack(bo + "I", next_ifd)
    return bytes(head), bytes(ext)


def build_tiff(
    rows: Sequence[Sequence[float]],
    *,
    byteorder: str = "<",
    sample_format: int = 3,
    bits: int = 32,
    compression: int = 1,
    predictor: int = 1,
    layout: str = "strips",
    rows_per_strip: int | None = None,
    tile_width: int = 16,
    tile_length: int = 16,
    photometric: int | None = 1,
    tags: Sequence[Tag] = (),
    extra_tags: Sequence[Tag] = (),
    drop_tags: Sequence[int] = (),
    tag_order: str = "sorted",
    raw_segments: Sequence[bytes] | None = None,
    magic: int = 42,
    ifd_offset: int | None = None,
    second_ifd: Sequence[Tag] | None = None,
) -> bytes:
    """Assemble a classic TIFF.

    ``rows`` is north row first.  ``tags`` replace (or add) tags by id, later
    entries winning; ``drop_tags`` removes ids after merging; ``extra_tags``
    are appended verbatim (so duplicates can be injected).  ``tag_order`` is
    ``"sorted"``, ``"reversed"`` or ``"given"``.  ``raw_segments`` supplies
    already-stored segment bytes (used for hand-made vectors and corrupt
    streams); the ``predictor`` and ``compression`` arguments then only affect
    the tags written.  File layout: header, segments, IFD, external values.
    """
    bo = byteorder
    height, width = len(rows), len(rows[0])
    zero = 0.0 if sample_format == 3 else 0

    def segment(block: Sequence[Sequence[float]]) -> bytes:
        raw = b"".join(encode_row(r, sample_format, bits, bo, predictor) for r in block)
        return compress(raw, compression)

    rps = rows_per_strip or height
    if raw_segments is not None:
        segments = list(raw_segments)
    elif layout == "strips":
        segments = [segment(rows[r : r + rps]) for r in range(0, height, rps)]
    elif layout == "tiles":
        segments = []
        for ty in range(-(-height // tile_length)):
            for tx in range(-(-width // tile_width)):
                block = []
                for r in range(ty * tile_length, (ty + 1) * tile_length):
                    line = [
                        rows[r][c] if (r < height and c < width) else zero
                        for c in range(tx * tile_width, (tx + 1) * tile_width)
                    ]
                    block.append(line)
                segments.append(segment(block))
    else:
        raise ValueError(layout)

    body = bytearray()
    offsets: list[int] = []
    counts: list[int] = []
    for seg in segments:
        offsets.append(8 + len(body))
        counts.append(len(seg))
        body += seg
        if len(body) % 2:
            body.append(0)
    ifd_pos = 8 + len(body)

    base: dict[int, Tag] = {
        IMAGE_WIDTH: long(IMAGE_WIDTH, width),
        IMAGE_LENGTH: long(IMAGE_LENGTH, height),
        BITS_PER_SAMPLE: short(BITS_PER_SAMPLE, bits),
        COMPRESSION: short(COMPRESSION, compression),
        SAMPLES_PER_PIXEL: short(SAMPLES_PER_PIXEL, 1),
        SAMPLE_FORMAT: short(SAMPLE_FORMAT, sample_format),
    }
    if photometric is not None:
        base[PHOTOMETRIC] = short(PHOTOMETRIC, photometric)
    if predictor != 1:
        base[PREDICTOR] = short(PREDICTOR, predictor)
    if layout == "strips":
        base[STRIP_OFFSETS] = long(STRIP_OFFSETS, *offsets)
        base[ROWS_PER_STRIP] = long(ROWS_PER_STRIP, rps)
        base[STRIP_BYTE_COUNTS] = long(STRIP_BYTE_COUNTS, *counts)
    else:
        base[TILE_WIDTH] = long(TILE_WIDTH, tile_width)
        base[TILE_LENGTH] = long(TILE_LENGTH, tile_length)
        base[TILE_OFFSETS] = long(TILE_OFFSETS, *offsets)
        base[TILE_BYTE_COUNTS] = long(TILE_BYTE_COUNTS, *counts)
    for tag in tags:
        base[tag.tag] = tag
    for tag_id in drop_tags:
        base.pop(tag_id, None)
    entries = list(base.values()) + list(extra_tags)
    if tag_order == "sorted":
        entries.sort(key=lambda t: t.tag)
    elif tag_order == "reversed":
        entries.sort(key=lambda t: t.tag, reverse=True)

    ifd, ext = _serialize_ifd(entries, ifd_pos, bo, 0)
    next_pos = 0
    if second_ifd is not None:
        end = ifd_pos + len(ifd) + len(ext)
        next_pos = end + (end % 2)
        ifd, ext = _serialize_ifd(entries, ifd_pos, bo, next_pos)
    out = bytearray(b"II" if bo == "<" else b"MM")
    out += struct.pack(bo + "H", magic)
    out += struct.pack(bo + "I", ifd_pos if ifd_offset is None else ifd_offset)
    out += body + ifd + ext
    if second_ifd is not None:
        if len(out) % 2:
            out.append(0)
        second, second_ext = _serialize_ifd(list(second_ifd), len(out), bo, 0)
        out += second + second_ext
    return bytes(out)


# --- inspecting and corrupting finished files --------------------------------
def byteorder_of(data: bytes) -> str:
    if data[:2] == b"II":
        return "<"
    if data[:2] == b"MM":
        return ">"
    raise ValueError("not a TIFF")


def ifd_offset_of(data: bytes) -> int:
    return struct.unpack_from(byteorder_of(data) + "I", data, 4)[0]


def ifd_span(data: bytes) -> tuple[int, int]:
    """``(start, end)`` of the first IFD (entry count through next-IFD pointer)."""
    bo = byteorder_of(data)
    start = ifd_offset_of(data)
    (n,) = struct.unpack_from(bo + "H", data, start)
    return start, start + 2 + 12 * n + 4


def entry_offset(data: bytes, tag: int, occurrence: int = 0) -> int:
    bo = byteorder_of(data)
    start = ifd_offset_of(data)
    (n,) = struct.unpack_from(bo + "H", data, start)
    seen = 0
    for i in range(n):
        off = start + 2 + 12 * i
        (found,) = struct.unpack_from(bo + "H", data, off)
        if found == tag:
            if seen == occurrence:
                return off
            seen += 1
    raise KeyError(tag)


def patch(data: bytes, offset: int, new: bytes) -> bytes:
    return data[:offset] + new + data[offset + len(new) :]


def truncate(data: bytes, length: int) -> bytes:
    return data[:length]


def set_entry(
    data: bytes,
    tag: int,
    *,
    type: int | None = None,  # noqa: A002 - mirrors the TIFF field name
    count: int | None = None,
    field: int | bytes | None = None,
) -> bytes:
    """Overwrite parts of the IFD entry for ``tag`` in a finished file."""
    bo = byteorder_of(data)
    off = entry_offset(data, tag)
    if type is not None:
        data = patch(data, off + 2, struct.pack(bo + "H", type))
    if count is not None:
        data = patch(data, off + 4, struct.pack(bo + "I", count))
    if field is not None:
        raw = struct.pack(bo + "I", field) if isinstance(field, int) else bytes(field)
        data = patch(data, off + 8, raw)
    return data
