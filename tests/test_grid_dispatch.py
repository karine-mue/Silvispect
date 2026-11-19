"""Tests for extension dispatch in ``Grid.read`` / ``Grid.write`` (SPEC-GIS.md section 1).

``.asc``, ``.tif`` and ``.tiff`` are routed by the lower-cased last suffix;
anything else is a ``GridError`` raised before the file system is touched.
The ``.asc`` tests also pin today's behaviour, so they are expected to pass
on the current code.
"""

from __future__ import annotations

import re

import pytest

import tiff_builder
from silvispect.grid import Grid, GridError

# A fixed grid with fractional header fields and nodata, and its exact text.
FIXED_VALUES = [1.0, 2.5, None, 0.125, None, 3.0]

FIXED_TEXT_P3 = (
    "ncols 3\n"
    "nrows 2\n"
    "xllcorner 100.5\n"
    "yllcorner 200.25\n"
    "cellsize 0.5\n"
    "NODATA_value -9999.0\n"
    "1 2.5 -9999.0\n"
    "0.125 -9999.0 3\n"
)

FIXED_TEXT_EXACT = (
    "ncols 3\n"
    "nrows 2\n"
    "xllcorner 100.5\n"
    "yllcorner 200.25\n"
    "cellsize 0.5\n"
    "NODATA_value -9999.0\n"
    "1.0 2.5 -9999.0\n"
    "0.125 -9999.0 3.0\n"
)

BAD_NAMES = [
    "chm.txt",
    "chm.dat",
    "chm.gpkg",
    "chm.las",
    "chm.tif.bak",
    "chm.asc.bak",
    "chm",
    "chm.",
    ".tif",
    "sub/chm",
]


def _fixed_grid() -> Grid:
    return Grid(3, 2, 100.5, 200.25, 0.5, -9999.0, list(FIXED_VALUES))


def _tif_bytes() -> bytes:
    return tiff_builder.build_tiff(
        [[1.0, 2.0], [3.0, 4.0]],
        tags=[
            tiff_builder.scale_tag(0.5, 0.5),
            tiff_builder.tiepoint_tag([0, 0, 0, 100.0, 200.0, 0]),
            tiff_builder.nodata_tag("-9999.0"),
        ],
    )


# --- 1.1 / 1.2: reading ----------------------------------------------------


@pytest.mark.parametrize("name", ["chm.asc", "CHM.ASC", "X.Asc", "plot.a1.asc"])
def test_read_dispatches_asc_in_any_case(tmp_path, name):
    path = tmp_path / name
    path.write_text(FIXED_TEXT_P3, encoding="utf-8")  # SPEC 1.1: .asc
    grid = Grid.read(path)
    assert (grid.ncols, grid.nrows) == (3, 2)
    assert grid.values == FIXED_VALUES


@pytest.mark.parametrize("name", ["chm.tif", "chm.tiff", "CHM.TIF", "a.Tiff", "x.TIFF", "v1.2.tif"])
def test_read_dispatches_tiff_in_any_case(tmp_path, name):
    # SPEC 1.1 / 1.2
    path = tmp_path / name
    path.write_bytes(_tif_bytes())
    grid = Grid.read(path)
    assert (grid.ncols, grid.nrows) == (2, 2)
    assert grid.values == [1.0, 2.0, 3.0, 4.0]
    assert grid.xllcorner == 100.0
    assert grid.yllcorner == 199.0
    assert grid.cellsize == 0.5


def test_read_accepts_str_paths_for_tif(tmp_path):
    path = tmp_path / "chm.tif"
    path.write_bytes(_tif_bytes())
    assert Grid.read(str(path)).values == [1.0, 2.0, 3.0, 4.0]


def test_directory_part_of_the_path_is_not_lowercased(tmp_path):
    folder = tmp_path / "MiXeD"
    folder.mkdir()
    path = folder / "Stem.TIF"
    path.write_bytes(_tif_bytes())
    assert Grid.read(path).ncols == 2


# --- 1.1 / 1.2: writing ----------------------------------------------------


@pytest.mark.parametrize("name", ["chm.tif", "chm.tiff", "CHM.TIF", "a.Tiff"])
def test_write_dispatches_tiff_in_any_case(tmp_path, name):
    # SPEC 1.1 / 1.2: the file is a classic little-endian TIFF, under exactly this name.
    path = tmp_path / name
    returned = _fixed_grid().write(path)
    assert returned == path
    assert [p.name for p in tmp_path.iterdir()] == [name]
    assert path.read_bytes()[:4] == b"II\x2a\x00"


@pytest.mark.parametrize("name", ["chm.asc", "X.ASC"])
def test_write_dispatches_asc_in_any_case(tmp_path, name):
    path = tmp_path / name
    _fixed_grid().write(path)
    assert [p.name for p in tmp_path.iterdir()] == [name]
    assert path.read_text(encoding="utf-8") == FIXED_TEXT_P3


def test_write_tif_creates_parent_directories(tmp_path):
    path = tmp_path / "a" / "b" / "chm.tif"
    _fixed_grid().write(path)
    assert path.is_file()


def test_write_asc_creates_parent_directories(tmp_path):
    path = tmp_path / "a" / "b" / "chm.asc"
    _fixed_grid().write(path)
    assert path.is_file()


def test_tif_written_through_dispatch_reads_back(tmp_path):
    path = tmp_path / "chm.tif"
    grid = _fixed_grid()
    grid.write(path)
    back = Grid.read(path)
    assert back.values == grid.values
    assert back.nodata_value == grid.nodata_value
    assert back.cellsize == grid.cellsize


def test_precision_does_not_change_tif_output(tmp_path):
    # SPEC 1.5: precision is .asc-only; .tif output is always lossless.
    grid = Grid(2, 1, 0.0, 0.0, 1.0, values=[0.1, 2.0])
    grid.write(tmp_path / "default.tif")
    grid.write(tmp_path / "p0.tif", precision=0)
    grid.write(tmp_path / "pn.tif", precision=None)
    reference = (tmp_path / "default.tif").read_bytes()
    assert (tmp_path / "p0.tif").read_bytes() == reference
    assert (tmp_path / "pn.tif").read_bytes() == reference
    assert Grid.read(tmp_path / "p0.tif").values == [0.1, 2.0]


def test_precision_still_quantises_asc(tmp_path):
    grid = Grid(2, 1, 0.0, 0.0, 1.0, values=[0.123456, 2.0])
    grid.write(tmp_path / "q.asc", precision=2)
    assert Grid.read(tmp_path / "q.asc").values == [0.12, 2.0]


# --- 1.3: unsupported extensions --------------------------------------------


@pytest.mark.parametrize("name", BAD_NAMES)
def test_read_rejects_unsupported_extension(tmp_path, name):
    with pytest.raises(GridError, match="unsupported raster file extension"):
        Grid.read(tmp_path / name)


@pytest.mark.parametrize("name", BAD_NAMES)
def test_write_rejects_unsupported_extension(tmp_path, name):
    with pytest.raises(GridError, match="unsupported raster file extension"):
        _fixed_grid().write(tmp_path / name)


def test_message_names_the_suffix_and_the_accepted_ones(tmp_path):
    pattern = r"unsupported raster file extension '\.txt': expected \.asc, \.tif or \.tiff"
    with pytest.raises(GridError, match=pattern):
        Grid.read(tmp_path / "chm.txt")
    with pytest.raises(GridError, match=pattern):
        _fixed_grid().write(tmp_path / "chm.txt")


def test_message_for_missing_extension_shows_empty_suffix(tmp_path):
    pattern = re.escape("unsupported raster file extension ''")
    with pytest.raises(GridError, match=pattern):
        Grid.read(tmp_path / "chm")
    with pytest.raises(GridError, match=pattern):
        Grid.read(tmp_path / "dir" / "chm")


def test_only_the_last_suffix_counts(tmp_path):
    with pytest.raises(GridError, match=re.escape("'.bak'")):
        Grid.read(tmp_path / "chm.tif.bak")


def test_hidden_file_named_dot_tif_has_no_extension(tmp_path):
    with pytest.raises(GridError, match=re.escape("unsupported raster file extension ''")):
        Grid.read(tmp_path / ".tif")


def test_unknown_extension_in_upper_case_is_rejected(tmp_path):
    with pytest.raises(GridError, match="unsupported raster file extension"):
        Grid.read(tmp_path / "chm.TXT")


def test_existing_ascii_grid_with_other_extension_is_refused(tmp_path):
    # SPEC 1.5: a valid ASCII Grid named grid.txt used to read; it no longer does.
    path = tmp_path / "grid.txt"
    path.write_text(FIXED_TEXT_P3, encoding="utf-8")
    with pytest.raises(GridError, match="unsupported raster file extension"):
        Grid.read(path)


# --- 1.3: the check precedes any file-system access ---------------------------


def test_read_of_missing_path_with_bad_extension_is_grid_error_not_oserror(tmp_path):
    missing = tmp_path / "nowhere" / "chm.txt"
    with pytest.raises(GridError) as info:
        Grid.read(missing)
    assert not isinstance(info.value, OSError)


def test_read_of_missing_asc_still_raises_file_not_found(tmp_path):
    with pytest.raises(FileNotFoundError):
        Grid.read(tmp_path / "missing.asc")


def test_read_of_missing_tif_raises_file_not_found(tmp_path):
    # SPEC section 2: OS errors propagate untouched.
    with pytest.raises(FileNotFoundError):
        Grid.read(tmp_path / "missing.tif")


def test_read_of_directory_named_like_a_tif_raises_oserror(tmp_path):
    (tmp_path / "d.tif").mkdir()
    with pytest.raises(OSError):
        Grid.read(tmp_path / "d.tif")


@pytest.mark.parametrize("name", BAD_NAMES)
def test_rejected_write_creates_neither_directory_nor_file(tmp_path, name):
    target = tmp_path / "newdir" / name
    with pytest.raises(GridError):
        _fixed_grid().write(target)
    assert not (tmp_path / "newdir").exists()
    assert list(tmp_path.iterdir()) == []


def test_rejected_write_into_existing_directory_leaves_it_empty(tmp_path):
    with pytest.raises(GridError):
        _fixed_grid().write(tmp_path / "chm.txt")
    assert list(tmp_path.iterdir()) == []


# --- 5.3 / 7.1: .asc regression ------------------------------------------------


def test_to_text_matches_the_fixed_literal():
    assert _fixed_grid().to_text() == FIXED_TEXT_P3


def test_to_text_without_quantisation_matches_the_fixed_literal():
    assert _fixed_grid().to_text(precision=None) == FIXED_TEXT_EXACT


def test_asc_file_bytes_match_the_fixed_literal(tmp_path):
    path = tmp_path / "fixed.asc"
    _fixed_grid().write(path)
    assert path.read_bytes() == FIXED_TEXT_P3.encode("utf-8")


def test_parse_of_fixed_text_has_no_crs():
    grid = Grid.parse(FIXED_TEXT_P3)
    assert grid.crs is None


def test_asc_read_has_no_crs(tmp_path):
    path = tmp_path / "fixed.asc"
    path.write_text(FIXED_TEXT_P3, encoding="utf-8")
    assert Grid.read(path).crs is None


def test_asc_round_trip_keeps_geometry_and_values(tmp_path):
    path = tmp_path / "fixed.asc"
    _fixed_grid().write(path, precision=None)
    back = Grid.read(path)
    assert back.values == FIXED_VALUES
    assert (back.xllcorner, back.yllcorner, back.cellsize) == (100.5, 200.25, 0.5)
    assert back.nodata_value == -9999.0
