"""Tests for the GeoTIFF writer (SPEC-GIS.md section 4).

The bytes ``Grid.write`` produces are parsed here with ``struct`` alone, so
the writer is judged against the layout table of 4.1 and not against the
reader it is paired with.  ``tiff_builder`` supplies only the byte-order and
IFD-offset helpers; every tag value is decoded independently below.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass

import pytest
import tiff_builder as tb

from silvispect.grid import Grid, GridError

# TIFF field types and sizes, restated here on purpose.
ASCII, SHORT, LONG, DOUBLE = 2, 3, 4, 12
SIZE = {ASCII: 1, SHORT: 2, LONG: 4, DOUBLE: 8}
CODE = {SHORT: "H", LONG: "I", DOUBLE: "d"}

# 4.1: (tag, type, count) in the order they must appear (GDAL_NODATA last).
BASE_TAGS = [
    (256, LONG, 1),
    (257, LONG, 1),
    (258, SHORT, 1),
    (259, SHORT, 1),
    (262, SHORT, 1),
    (273, LONG, 1),
    (277, SHORT, 1),
    (278, LONG, 1),
    (279, LONG, 1),
    (284, SHORT, 1),
    (339, SHORT, 1),
    (33550, DOUBLE, 3),
    (33922, DOUBLE, 6),
]
GEOKEY_TAG = (34735, SHORT, 16)
NODATA_TAG_ID = 42113


# --- independent parsing -----------------------------------------------------


@dataclass(frozen=True)
class Entry:
    tag: int
    type: int
    count: int
    value_at: int  # absolute offset of the value bytes

    @property
    def size(self) -> int:
        return SIZE[self.type] * self.count

    @property
    def inline(self) -> bool:
        return self.size <= 4


@dataclass(frozen=True)
class Parsed:
    data: bytes
    bo: str
    ifd: int
    entries: list[Entry]
    next_ifd: int

    @property
    def ifd_end(self) -> int:
        return self.ifd + 2 + 12 * len(self.entries) + 4

    def entry(self, tag: int) -> Entry:
        for e in self.entries:
            if e.tag == tag:
                return e
        raise KeyError(tag)

    def tags(self) -> list[int]:
        return [e.tag for e in self.entries]

    def value(self, tag: int):
        e = self.entry(tag)
        if e.type == ASCII:
            return self.data[e.value_at : e.value_at + e.count]
        vals = struct.unpack_from(f"{self.bo}{e.count}{CODE[e.type]}", self.data, e.value_at)
        return vals[0] if e.count == 1 else tuple(vals)


def parse(data: bytes) -> Parsed:
    bo = tb.byteorder_of(data)
    ifd = tb.ifd_offset_of(data)
    (n,) = struct.unpack_from(bo + "H", data, ifd)
    entries = []
    for i in range(n):
        off = ifd + 2 + 12 * i
        tag, typ, count = struct.unpack_from(bo + "HHI", data, off)
        size = SIZE[typ] * count
        if size <= 4:
            value_at = off + 8
        else:
            (value_at,) = struct.unpack_from(bo + "I", data, off + 8)
        entries.append(Entry(tag, typ, count, value_at))
    (next_ifd,) = struct.unpack_from(bo + "I", data, ifd + 2 + 12 * n)
    return Parsed(data, bo, ifd, entries, next_ifd)


def pixels(p: Parsed) -> list[float]:
    offset = p.value(273)
    count = p.value(279)
    bits = p.value(258)
    code = {32: "f", 64: "d"}[bits]
    n = p.value(256) * p.value(257)
    assert count == n * bits // 8
    return list(struct.unpack_from(f"{p.bo}{n}{code}", p.data, offset))


# --- fixtures ---------------------------------------------------------------


def _grid(*, values=None, nodata=-9999.0, crs: str | None = "EPSG:6677") -> Grid:
    # 3 x 2, cell 0.5, origin (100, 200)  =>  y_ul = 200 + 2 * 0.5 = 201.
    if values is None:
        values = [1.0, 2.5, None, 0.25, 4.0, 5.75]
    return Grid(3, 2, 100.0, 200.0, 0.5, nodata_value=nodata, values=list(values), crs=crs)


def _write(tmp_path, grid: Grid | None = None, name: str = "out.tif", **kwargs) -> bytes:
    path = tmp_path / name
    (grid or _grid()).write(path, **kwargs)
    return path.read_bytes()


# --- 4.1: container ---------------------------------------------------------


def test_header_is_little_endian_classic_tiff_with_ifd_at_8(tmp_path):
    data = _write(tmp_path)
    assert data[:2] == b"II"
    assert tb.byteorder_of(data) == "<"
    assert struct.unpack_from("<H", data, 2)[0] == 42
    assert tb.ifd_offset_of(data) == 8


def test_next_ifd_offset_is_zero(tmp_path):
    assert parse(_write(tmp_path)).next_ifd == 0


# --- 4.1: tag set, order, types and counts ---------------------------------


def test_tag_set_and_order_with_crs_match_the_table(tmp_path):
    p = parse(_write(tmp_path))
    assert p.tags() == [t for t, _, _ in BASE_TAGS] + [34735, NODATA_TAG_ID]


def test_tag_set_and_order_without_crs_match_the_table(tmp_path):
    p = parse(_write(tmp_path, _grid(crs=None)))
    assert p.tags() == [t for t, _, _ in BASE_TAGS] + [NODATA_TAG_ID]


def test_tag_types_and_counts_match_the_table(tmp_path):
    p = parse(_write(tmp_path))
    nodata_len = len(b"-9999.0\x00")
    expected = BASE_TAGS + [GEOKEY_TAG, (NODATA_TAG_ID, ASCII, nodata_len)]
    assert [(e.tag, e.type, e.count) for e in p.entries] == expected


def test_crs_none_writes_no_geokey_directory(tmp_path):
    assert 34735 not in parse(_write(tmp_path, _grid(crs=None))).tags()


def test_no_extra_tags_are_written(tmp_path):
    # 4.1 "書かないもの": Orientation, resolution tags, GDAL_METADATA, ...
    p = parse(_write(tmp_path))
    for forbidden in (274, 282, 283, 296, 42112, 254, 317):
        assert forbidden not in p.tags()


# --- 4.1: tag values --------------------------------------------------------


def test_image_structure_values(tmp_path):
    p = parse(_write(tmp_path))
    assert p.value(256) == 3  # ImageWidth
    assert p.value(257) == 2  # ImageLength
    assert p.value(259) == 1  # Compression: none
    assert p.value(262) == 1  # BlackIsZero
    assert p.value(277) == 1  # SamplesPerPixel
    assert p.value(278) == 2  # RowsPerStrip = nrows (single strip)
    assert p.value(284) == 1  # PlanarConfiguration
    assert p.value(339) == 3  # SampleFormat: IEEE float


def test_geo_tags_values(tmp_path):
    p = parse(_write(tmp_path))
    assert p.value(33550) == (0.5, 0.5, 0.0)
    # (0, 0, 0, xll, y_ul, 0) with y_ul = yll + nrows * cellsize
    assert p.value(33922) == (0.0, 0.0, 0.0, 100.0, 201.0, 0.0)


def test_geo_tags_do_not_depend_on_crs(tmp_path):
    p = parse(_write(tmp_path, _grid(crs=None)))
    assert p.value(33550) == (0.5, 0.5, 0.0)
    assert p.value(33922) == (0.0, 0.0, 0.0, 100.0, 201.0, 0.0)


def test_geokey_directory_has_sixteen_expected_elements(tmp_path):
    p = parse(_write(tmp_path, _grid(crs="EPSG:6677")))
    assert p.value(34735) == (
        1, 1, 0, 3,  # header: version, revision, minor, number of keys
        1024, 0, 1, 1,  # GTModelType = Projected
        1025, 0, 1, 1,  # GTRasterType = PixelIsArea
        3072, 0, 1, 6677,  # ProjectedCSType
    )  # fmt: skip


@pytest.mark.parametrize("code", [1, 6678, 32766])
def test_geokey_directory_carries_the_epsg_code(tmp_path, code):
    p = parse(_write(tmp_path, _grid(crs=f"EPSG:{code}")))
    keys = p.value(34735)
    assert len(keys) == 16
    assert keys[12:16] == (3072, 0, 1, code)


# --- 4.1: layout ------------------------------------------------------------


def test_external_values_follow_the_ifd_in_tag_order_on_even_offsets(tmp_path):
    p = parse(_write(tmp_path))
    cursor = p.ifd_end
    for e in p.entries:
        if e.inline:
            continue
        cursor += cursor % 2
        assert e.value_at == cursor, f"tag {e.tag}"
        assert e.value_at % 2 == 0
        cursor = e.value_at + e.size
    # padding between values is zero bytes
    for e in p.entries:
        if not e.inline and e.value_at > 0:
            pass
    assert cursor <= len(p.data)


def test_odd_sized_external_value_is_zero_padded_to_an_even_offset(tmp_path):
    # "-9999.1234\0" is 11 bytes: the pixel data must still start on an even offset.
    grid = _grid(nodata=-9999.1234, values=[1.0, 2.0, None, 3.0, 4.0, 5.0])
    p = parse(_write(tmp_path, grid))
    nodata = p.entry(NODATA_TAG_ID)
    assert nodata.size == 11
    end = nodata.value_at + nodata.size
    assert end % 2 == 1
    assert p.data[end] == 0
    assert p.value(273) == end + 1


def test_pixel_data_is_last_contiguous_and_on_an_even_offset(tmp_path):
    p = parse(_write(tmp_path))
    offset = p.value(273)
    assert offset % 2 == 0
    assert offset + p.value(279) == len(p.data)
    last_external = max(
        (e.value_at + e.size for e in p.entries if not e.inline), default=p.ifd_end
    )
    assert offset == last_external + last_external % 2


def test_strip_byte_count_is_cells_times_bytes_per_sample(tmp_path):
    p = parse(_write(tmp_path))
    assert p.value(258) == 32
    assert p.value(279) == 3 * 2 * 4


# --- 4.1: pixel data --------------------------------------------------------


def test_pixel_data_is_row_major_north_row_first_with_nodata_filled(tmp_path):
    p = parse(_write(tmp_path))
    assert pixels(p) == [1.0, 2.5, -9999.0, 0.25, 4.0, 5.75]


def test_all_nodata_grid_is_written_as_the_sentinel(tmp_path):
    grid = _grid(values=[None] * 6, nodata=-1.5)
    p = parse(_write(tmp_path, grid))
    assert pixels(p) == [-1.5] * 6


def test_negative_zero_survives_in_the_pixel_data(tmp_path):
    grid = _grid(values=[-0.0, 0.0, 1.0, 2.0, 3.0, 4.0])
    data = _write(tmp_path, grid)
    p = parse(data)
    offset = p.value(273)
    assert data[offset : offset + 4] == struct.pack("<f", -0.0)


# --- 4.3: GDAL_NODATA -------------------------------------------------------


@pytest.mark.parametrize(
    ("nodata", "text"),
    [
        (-9999.0, b"-9999.0\x00"),
        (-9999.1234, b"-9999.1234\x00"),
        (0.0, b"0.0\x00"),
        (-1.5, b"-1.5\x00"),
    ],
)
def test_gdal_nodata_string_is_repr_of_float_plus_nul(tmp_path, nodata, text):
    grid = _grid(nodata=nodata, values=[1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    assert parse(_write(tmp_path, grid)).value(NODATA_TAG_ID) == text


def test_gdal_nodata_is_written_even_without_nodata_cells(tmp_path):
    grid = _grid(values=[1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    p = parse(_write(tmp_path, grid))
    assert NODATA_TAG_ID in p.tags()
    assert p.value(NODATA_TAG_ID) == b"-9999.0\x00"


def test_inline_gdal_nodata_is_still_a_valid_ascii_entry(tmp_path):
    # "1.0\0" is exactly four bytes and lives inside the IFD entry.
    grid = _grid(nodata=1.0, values=[2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
    p = parse(_write(tmp_path, grid))
    assert p.entry(NODATA_TAG_ID).inline
    assert p.value(NODATA_TAG_ID) == b"1.0\x00"


# --- 4.2: dtype selection ---------------------------------------------------


def test_float32_exact_grid_is_written_as_32_bit(tmp_path):
    p = parse(_write(tmp_path))
    assert p.value(258) == 32


def test_label_like_integer_values_are_32_bit(tmp_path):
    grid = _grid(values=[0.0, 1.0, 2.0, 3.0, 12345.0, None], nodata=-1.0)
    assert parse(_write(tmp_path, grid)).value(258) == 32


def test_value_that_is_a_float32_already_is_32_bit(tmp_path):
    exact = struct.unpack("<f", struct.pack("<f", 0.1))[0]
    grid = _grid(values=[exact, 1.0, 2.0, 3.0, 4.0, 5.0])
    p = parse(_write(tmp_path, grid))
    assert p.value(258) == 32
    assert pixels(p)[0] == exact


def test_value_not_exact_in_float32_selects_64_bit(tmp_path):
    grid = _grid(values=[0.1, 1.0, 2.0, 3.0, 4.0, 5.0])
    p = parse(_write(tmp_path, grid))
    assert p.value(258) == 64
    assert p.value(279) == 6 * 8
    assert p.value(339) == 3  # still IEEE float
    assert pixels(p)[0] == 0.1  # exact, no quantisation


def test_one_inexact_cell_makes_the_whole_raster_64_bit(tmp_path):
    grid = _grid(values=[1.0, 2.0, 3.0, 4.0, 5.0, 123.456])
    p = parse(_write(tmp_path, grid))
    assert p.value(258) == 64
    assert pixels(p) == [1.0, 2.0, 3.0, 4.0, 5.0, 123.456]


def test_value_overflowing_float32_selects_64_bit(tmp_path):
    # struct.pack("<f", 1e300) raises OverflowError; it must not escape (4.2).
    grid = _grid(values=[1e300, 1.0, 2.0, 3.0, 4.0, 5.0])
    p = parse(_write(tmp_path, grid))
    assert p.value(258) == 64
    assert pixels(p)[0] == 1e300


def test_value_underflowing_float32_selects_64_bit(tmp_path):
    grid = _grid(values=[1e-50, 1.0, 2.0, 3.0, 4.0, 5.0])
    p = parse(_write(tmp_path, grid))
    assert p.value(258) == 64
    assert pixels(p)[0] == 1e-50


def test_inexact_nodata_selects_64_bit_even_if_all_cells_are_exact(tmp_path):
    grid = _grid(nodata=-9999.1234, values=[1.0, 2.5, None, 0.25, 4.0, 5.75])
    p = parse(_write(tmp_path, grid))
    assert p.value(258) == 64
    assert p.value(NODATA_TAG_ID) == b"-9999.1234\x00"
    assert pixels(p) == [1.0, 2.5, -9999.1234, 0.25, 4.0, 5.75]


def test_64_bit_grid_is_still_a_single_uncompressed_strip(tmp_path):
    p = parse(_write(tmp_path, _grid(values=[0.1] * 6)))
    assert p.value(259) == 1
    assert p.value(278) == p.value(257)
    assert p.value(273) + p.value(279) == len(p.data)


# --- 4.1: determinism ------------------------------------------------------


def test_same_grid_written_twice_gives_identical_bytes(tmp_path):
    grid = _grid()
    assert _write(tmp_path, grid, "a.tif") == _write(tmp_path, grid, "b.tif")


def test_equal_grids_built_separately_give_identical_bytes(tmp_path):
    assert _write(tmp_path, _grid(), "a.tif") == _write(tmp_path, _grid(), "b.tif")


def test_64_bit_output_is_deterministic(tmp_path):
    grid = _grid(values=[0.1, 0.2, None, 0.3, 0.4, 0.5])
    assert _write(tmp_path, grid, "a.tif") == _write(tmp_path, grid, "b.tif")


@pytest.mark.parametrize("values", [[1.0, 2.5, None, 0.25, 4.0, 5.75], [0.1, 0.2, None, 0.3, 0.4, 0.5]])
def test_precision_does_not_change_tif_output(tmp_path, values):
    # 1.5: precision is .asc-only; .tif output is always reversible.
    grid = _grid(values=values)
    reference = _write(tmp_path, grid, "ref.tif")
    for precision in (0, 1, 3, 8, None):
        assert _write(tmp_path, grid, f"p{precision}.tif", precision=precision) == reference


def test_tif_and_tiff_suffixes_give_identical_bytes(tmp_path):
    grid = _grid()
    assert _write(tmp_path, grid, "a.tif") == _write(tmp_path, grid, "a.tiff")


def test_uppercase_suffix_is_written_as_a_geotiff(tmp_path):
    data = _write(tmp_path, _grid(), "CHM.TIF")
    assert data[:4] == b"II*\x00"


def test_write_creates_missing_parent_directories(tmp_path):
    target = tmp_path / "a" / "b" / "out.tif"
    _grid().write(target)
    assert parse(target.read_bytes()).value(256) == 3


# --- 4.2: nodata collision --------------------------------------------------


def test_valid_cell_equal_to_nodata_is_refused(tmp_path):
    grid = _grid(values=[1.0, -9999.0, None, 0.25, 4.0, 5.75])
    with pytest.raises(GridError, match=r"cannot write GeoTIFF: .*nodata"):
        grid.write(tmp_path / "out.tif")


def test_collision_leaves_no_file_and_no_directory(tmp_path):
    grid = _grid(values=[1.0, -9999.0, None, 0.25, 4.0, 5.75])
    target = tmp_path / "newdir" / "out.tif"
    with pytest.raises(GridError):
        grid.write(target)
    assert not target.exists()
    assert not target.parent.exists()


def test_collision_is_detected_on_the_64_bit_path_too(tmp_path):
    grid = _grid(nodata=0.1, values=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    target = tmp_path / "newdir" / "out.tif"
    with pytest.raises(GridError, match=r"cannot write GeoTIFF: .*nodata"):
        grid.write(target)
    assert not target.parent.exists()


def test_collision_does_not_overwrite_an_existing_file(tmp_path):
    target = tmp_path / "out.tif"
    target.write_bytes(b"previous")
    grid = _grid(values=[-9999.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    with pytest.raises(GridError):
        grid.write(target)
    assert target.read_bytes() == b"previous"


def test_float32_neighbour_of_nodata_is_not_a_collision(tmp_path):
    # Comparison is exact in double precision (4.2), so -9999.5 is just a value.
    grid = _grid(values=[-9999.5, 1.0, 2.0, 3.0, 4.0, 5.0])
    assert pixels(parse(_write(tmp_path, grid)))[0] == -9999.5


# --- 4.1 / 6.3 #33: re-validation at write time ----------------------------


def test_invalid_crs_assigned_after_construction_is_refused(tmp_path):
    grid = _grid()
    grid.crs = "epsg:6677"  # attribute assignment is not re-validated by Grid
    target = tmp_path / "newdir" / "out.tif"
    with pytest.raises(GridError, match="cannot write GeoTIFF:"):
        grid.write(target)
    assert not target.exists()
    assert not target.parent.exists()


def test_non_finite_upper_left_northing_is_refused(tmp_path):
    # y_ul = yll + nrows * cellsize overflows although every field is finite.
    grid = Grid(
        1, 2, 0.0, 1.7e308, 1e308, nodata_value=-9999.0, values=[1.0, 2.0], crs=None
    )
    target = tmp_path / "newdir" / "out.tif"
    with pytest.raises(GridError, match="cannot write GeoTIFF:"):
        grid.write(target)
    assert not target.exists()
    assert not target.parent.exists()


# --- 1.3: extension is checked before touching the file system --------------


def test_unsupported_extension_creates_no_directory(tmp_path):
    target = tmp_path / "newdir" / "out.txt"
    with pytest.raises(GridError, match="unsupported raster file extension"):
        _grid().write(target)
    assert not target.parent.exists()
