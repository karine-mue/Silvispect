"""Regression tests — Lane A (grid.py). Issues #6, #11, #13, #23, #27, #28."""

import pytest

from silvispect.grid import Grid, GridError


def make_grid(values, *, ncols=2, nrows=2, cellsize=0.5, nodata=-9999.0):
    return Grid(
        ncols=ncols,
        nrows=nrows,
        xllcorner=100.25,
        yllcorner=200.5,
        cellsize=cellsize,
        nodata_value=nodata,
        values=list(values),
    )


# --- #6: header fields must not be rounded -------------------------------

def test_header_cellsize_not_rounded_precision0():
    grid = make_grid([1.4, 2.6, 3.5, 4.5], cellsize=0.5)
    text = grid.to_text(precision=0)
    assert "cellsize 0.5" in text
    back = Grid.parse(text)
    assert back.cellsize == 0.5
    assert back.xllcorner == 100.25


def test_header_cellsize_0125_survives_precision2():
    grid = make_grid([1.0, 2.0, 3.0, 4.0], cellsize=0.125)
    text = grid.to_text(precision=2)
    assert "cellsize 0.125" in text
    assert Grid.parse(text).cellsize == 0.125


# --- #11: nodata collision on write ---------------------------------------

def test_write_rejects_value_equal_to_nodata():
    grid = make_grid([0.0, 5.0, 5.0, 5.0], nodata=0.0)
    with pytest.raises(GridError):
        grid.to_text()


def test_write_rejects_value_rounding_into_nodata():
    grid = make_grid([0.0004, 5.0, 5.0, 5.0], nodata=0.0)
    with pytest.raises(GridError):
        grid.to_text(precision=3)


def test_write_normal_grid_still_works():
    grid = make_grid([1.0, None, 3.0, 4.0], nodata=0.0)
    back = Grid.parse(grid.to_text())
    assert back.values == [1.0, None, 3.0, 4.0]


# --- #13: finite-value enforcement and values handling ---------------------

def test_construct_rejects_nan():
    with pytest.raises(GridError):
        make_grid([float("nan"), 1.0, 1.0, 1.0])


def test_construct_rejects_inf():
    with pytest.raises(GridError):
        make_grid([float("inf"), 1.0, 1.0, 1.0])


def test_set_rejects_nan():
    grid = make_grid([1.0, 2.0, 3.0, 4.0])
    with pytest.raises(GridError):
        grid.set(0, 0, float("nan"))


def test_replace_nodata_rejects_inf():
    grid = make_grid([None, 2.0, 3.0, 4.0])
    with pytest.raises(GridError):
        grid.replace_nodata(float("inf"))


def test_map_values_rejects_nan_result():
    grid = make_grid([1.0, 2.0, 3.0, 4.0])
    with pytest.raises(GridError):
        grid.map_values(lambda v: float("nan"))


def test_focal_rejects_nonfinite_result():
    grid = make_grid([1.0, 2.0, 3.0, 4.0])
    with pytest.raises(GridError):
        grid.focal(lambda vs: float("inf"))


def test_values_empty_list_is_an_error():
    with pytest.raises(GridError):
        Grid(ncols=2, nrows=2, xllcorner=0.0, yllcorner=0.0,
             cellsize=1.0, nodata_value=-9999.0, values=[])


def test_values_omitted_means_all_nodata():
    grid = Grid(ncols=2, nrows=2, xllcorner=0.0, yllcorner=0.0,
                cellsize=1.0, nodata_value=-9999.0)
    assert grid.values == [None, None, None, None]


# --- #28: parse robustness --------------------------------------------------

def test_parse_accepts_utf8_bom():
    grid = make_grid([1.0, 2.0, 3.0, 4.0])
    back = Grid.parse("\ufeff" + grid.to_text())
    assert back.values == [1.0, 2.0, 3.0, 4.0]


def test_parse_rejects_extra_header_token():
    grid = make_grid([1.0, 2.0, 3.0, 4.0])
    lines = grid.to_text().splitlines()
    lines[0] = lines[0] + " 99"
    with pytest.raises(GridError):
        Grid.parse("\n".join(lines))


# --- #27 (histogram part) ---------------------------------------------------

def test_histogram_rejects_lo_greater_than_hi():
    grid = make_grid([1.0, 2.0, 3.0, 4.0])
    with pytest.raises(GridError):
        grid.histogram(lo=10.0, hi=5.0)


def test_histogram_rejects_bool_bins():
    grid = make_grid([1.0, 2.0, 3.0, 4.0])
    with pytest.raises(GridError):
        grid.histogram(bins=True)


# --- #23 (median part) --------------------------------------------------------

def test_smooth_median_does_not_overflow():
    grid = Grid.from_rows([[1e308, 1e308], [1e308, 1e308]])
    out = grid.smooth_median(radius=1)
    assert all(v == 1e308 for v in out.values)


# --- #14: far-edge coordinates belong to the edge cells -------------------------

def test_cell_of_east_and_south_edges():
    grid = Grid.from_rows([[5.0, 5.0, 5.0]] * 3, cellsize=1.0)
    assert grid.cell_of(3.0, 1.5) == (1, 2)
    assert grid.cell_of(1.5, 0.0) == (2, 1)
    assert grid.cell_of(3.0, 0.0) == (2, 2)
    assert grid.cell_of(3.5, 1.5) is None
    assert grid.cell_of(1.5, -0.5) is None


def test_edge_tree_is_observed_and_not_flagged_outside():
    from silvispect.inspection import inspect_stand
    from silvispect.inventory import Tree

    chm = Grid.from_rows([[5.0, 5.0, 5.0]] * 3, cellsize=1.0)
    edge = Tree("EDGE", 3.0, 1.5, dbh_cm=20.0, height_m=5.0)
    outside = Tree("OUT", 3.5, 1.5, dbh_cm=20.0, height_m=5.0)
    report = inspect_stand(chm=chm, trees=[edge, outside], area_ha=0.0009)
    flagged = {f.subject for f in report.findings if f.code == "SV043"}
    assert "tree:OUT" in flagged
    assert "tree:EDGE" not in flagged
