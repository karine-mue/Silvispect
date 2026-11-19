"""Pinning tests added by the post-green completeness sweep (2026-11-19).

The sweep (Sonnet, after the seat swap) reviewed the integrated GeoTIFF
slice against SPEC-GIS.md and returned ADD verdicts for C1-C4 of the
ChatGPT boundary list, plus a regression test for the B1 adjudication
(overlapping segments are read, not detected).  Each test here pins
behaviour the main test files left implicit.  Files are built with
``tests/tiff_builder.py`` and, where the builder cannot express the
layout (physical segment placement, shared offsets), patched at the
byte level afterwards.
"""

from __future__ import annotations

import struct
import zlib

import pytest

from silvispect.grid import Grid, GridError

try:
    from tests import tiff_builder as tb
except ImportError:  # tests/ is not a package: its directory is on sys.path
    import tiff_builder as tb  # type: ignore[no-redef]

_GEO = [tb.scale_tag(1.0, 1.0, 0.0), tb.tiepoint_tag((0.0, 0.0, 0.0, 1000.0, 2000.0, 0.0))]


def _rows(height: int, width: int) -> list[list[float]]:
    return [[float(r * width + c + 1) for c in range(width)] for r in range(height)]


def _write(tmp_path, data: bytes, name: str = "chm.tif"):
    path = tmp_path / name
    path.write_bytes(data)
    return path


def _patch_segment_offsets(data: bytes, tag_id: int, new_offsets: list[int]) -> bytes:
    """Rewrite a strip/tile offsets tag (LONG array, held externally)."""
    out = bytearray(data)
    ifd_off = struct.unpack_from("<I", out, 4)[0]
    count = struct.unpack_from("<H", out, ifd_off)[0]
    for i in range(count):
        entry = ifd_off + 2 + i * 12
        if struct.unpack_from("<H", out, entry)[0] != tag_id:
            continue
        n = struct.unpack_from("<I", out, entry + 4)[0]
        assert n == len(new_offsets)
        if n == 1:
            struct.pack_into("<I", out, entry + 8, new_offsets[0])
        else:
            values_at = struct.unpack_from("<I", out, entry + 8)[0]
            for j, value in enumerate(new_offsets):
                struct.pack_into("<I", out, values_at + 4 * j, value)
        return bytes(out)
    raise AssertionError(f"tag {tag_id} not found")


def _segment_bytes(data: bytes, offsets_tag: int, counts_tag: int) -> tuple[list[int], list[int]]:
    """Return (offsets, byte_counts) for a strip or tile tag pair."""
    ifd_off = struct.unpack_from("<I", data, 4)[0]
    count = struct.unpack_from("<H", data, ifd_off)[0]
    found: dict[int, list[int]] = {}
    for i in range(count):
        entry = ifd_off + 2 + i * 12
        tag = struct.unpack_from("<H", data, entry)[0]
        if tag not in (offsets_tag, counts_tag):
            continue
        n = struct.unpack_from("<I", data, entry + 4)[0]
        if n == 1:
            found[tag] = [struct.unpack_from("<I", data, entry + 8)[0]]
        else:
            values_at = struct.unpack_from("<I", data, entry + 8)[0]
            found[tag] = [struct.unpack_from("<I", data, values_at + 4 * j)[0] for j in range(n)]
    return found[offsets_tag], found[counts_tag]


# --- B1 adjudication: overlapping segments are read, not detected ------------


def test_two_strips_sharing_one_offset_read_the_same_bytes_twice(tmp_path):
    data = tb.build_tiff(_rows(4, 2), tags=_GEO, rows_per_strip=2)
    offsets, _ = _segment_bytes(data, tb.STRIP_OFFSETS, tb.STRIP_BYTE_COUNTS)
    data = _patch_segment_offsets(data, tb.STRIP_OFFSETS, [offsets[0], offsets[0]])
    grid = Grid.read(_write(tmp_path, data))
    assert grid.values == [1.0, 2.0, 3.0, 4.0, 1.0, 2.0, 3.0, 4.0]


# --- C1 (+C6): physical placement is not logical order -----------------------


def test_strips_physically_reversed_still_read_in_logical_order(tmp_path):
    rows = _rows(4, 2)
    data = tb.build_tiff(rows, tags=_GEO, rows_per_strip=2)
    offsets, counts = _segment_bytes(data, tb.STRIP_OFFSETS, tb.STRIP_BYTE_COUNTS)
    first = data[offsets[0] : offsets[0] + counts[0]]
    second = data[offsets[1] : offsets[1] + counts[1]]
    assert len(first) == len(second)
    out = bytearray(data)
    out[offsets[0] : offsets[0] + counts[0]] = second
    out[offsets[1] : offsets[1] + counts[1]] = first
    out = _patch_segment_offsets(bytes(out), tb.STRIP_OFFSETS, [offsets[1], offsets[0]])
    out += b"\xde\xad\xbe\xef unreferenced trailing bytes"
    grid = Grid.read(_write(tmp_path, out))
    assert grid.values == [v for row in rows for v in row]


def test_tiles_physically_reversed_still_read_in_logical_order(tmp_path):
    rows = _rows(2, 4)
    data = tb.build_tiff(rows, tags=_GEO, layout="tiles", tile_width=2, tile_length=2)
    offsets, counts = _segment_bytes(data, tb.TILE_OFFSETS, tb.TILE_BYTE_COUNTS)
    assert len(offsets) == 2
    first = data[offsets[0] : offsets[0] + counts[0]]
    second = data[offsets[1] : offsets[1] + counts[1]]
    assert len(first) == len(second)
    out = bytearray(data)
    out[offsets[0] : offsets[0] + counts[0]] = second
    out[offsets[1] : offsets[1] + counts[1]] = first
    out = _patch_segment_offsets(bytes(out), tb.TILE_OFFSETS, [offsets[1], offsets[0]])
    grid = Grid.read(_write(tmp_path, out))
    assert grid.values == [v for row in rows for v in row]


# --- C2: deflate stream variants ---------------------------------------------


def test_raw_deflate_stream_is_corrupt(tmp_path):
    compressor = zlib.compressobj(9, zlib.DEFLATED, -15)
    raw = compressor.compress(struct.pack("<2f", 1.0, 2.0)) + compressor.flush()
    data = tb.build_tiff([[0.0, 0.0]], tags=_GEO, compression=8, raw_segments=[raw])
    with pytest.raises(GridError, match="corrupt TIFF"):
        Grid.read(_write(tmp_path, data))


def test_deflate_with_a_bad_adler32_checksum_is_corrupt(tmp_path):
    stream = bytearray(zlib.compress(struct.pack("<f", 1.0)))
    stream[-1] ^= 0xFF
    data = tb.build_tiff([[0.0, 0.0]], tags=_GEO, compression=8, raw_segments=[bytes(stream)])
    with pytest.raises(GridError, match="corrupt TIFF"):
        Grid.read(_write(tmp_path, data))


# --- C3: PackBits control-byte boundaries in one complete stream --------------


def test_packbits_control_byte_boundaries(tmp_path):
    literals = bytes(range(128))
    stream = b"".join(
        [
            bytes([0, 7]),  # control 0: one literal
            bytes([127]) + literals,  # control 127: 128 literals
            bytes([0x81, 200]),  # control -127: 128 repeats
            bytes([0x80, 0x80]),  # control -128: no-ops
            bytes([0xFF, 55]),  # control -1: two repeats
        ]
    )
    data = tb.build_tiff(
        [[0] * 259], tags=_GEO, sample_format=1, bits=8, compression=32773, raw_segments=[stream]
    )
    grid = Grid.read(_write(tmp_path, data))
    assert grid.values == [7.0] + [float(v) for v in range(128)] + [200.0] * 128 + [55.0, 55.0]


# --- C4: big-endian non-finite nodata -----------------------------------------


def test_big_endian_infinite_nodata_marks_matching_cells(tmp_path):
    data = tb.build_tiff(
        [[float("inf"), 2.0]],
        tags=[*_GEO, tb.nodata_tag("inf")],
        byteorder=">",
    )
    grid = Grid.read(_write(tmp_path, data))
    assert grid.values == [None, 2.0]
    assert grid.nodata_value == -9999.0


def test_big_endian_infinite_nodata_with_sentinel_collision_is_refused(tmp_path):
    data = tb.build_tiff(
        [[float("inf"), -9999.0]],
        tags=[*_GEO, tb.nodata_tag("inf")],
        byteorder=">",
    )
    with pytest.raises(GridError, match="nodata sentinel collision"):
        Grid.read(_write(tmp_path, data))


# --- D4: projected CRS codes above the Grid range ------------------------------


def test_projected_crs_code_above_32766_is_refused_as_unsupported_crs(tmp_path):
    keys = [(1024, 0, 1, 1), (1025, 0, 1, 1), (3072, 0, 1, 54009)]
    data = tb.build_tiff([[1.0]], tags=[*_GEO, tb.geokey_tag(*keys)])
    with pytest.raises(GridError, match="unsupported CRS"):
        Grid.read(_write(tmp_path, data))
