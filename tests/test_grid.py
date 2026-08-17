"""Tests for the ASCII grid raster."""

from __future__ import annotations

import math

import pytest

from silvispect.grid import Grid, GridError

SAMPLE = """\
ncols 3
nrows 2
xllcorner 100.0
yllcorner 200.0
cellsize 2.0
NODATA_value -9999
1.0 2.0 3.0
4.0 -9999 6.0
"""


def test_parse_round_trip():
    grid = Grid.parse(SAMPLE)
    assert (grid.ncols, grid.nrows) == (3, 2)
    assert grid.cellsize == 2.0
    assert grid.get(0, 0) == 1.0
    assert grid.get(1, 1) is None
    reparsed = Grid.parse(grid.to_text())
    assert reparsed.values == grid.values
    assert reparsed.extent == grid.extent


def test_write_and_read(tmp_path):
    grid = Grid.parse(SAMPLE)
    path = grid.write(tmp_path / "nested" / "grid.asc")
    assert path.exists()
    assert Grid.read(path).values == grid.values


def test_header_is_case_insensitive_and_accepts_centers():
    text = SAMPLE.replace("xllcorner 100.0", "XLLCENTER 101.0").replace(
        "NODATA_value", "nodata_value"
    )
    grid = Grid.parse(text)
    assert grid.xllcorner == pytest.approx(100.0)


def test_parse_rejects_wrong_cell_count():
    with pytest.raises(GridError, match="expected 6 cell values"):
        Grid.parse(SAMPLE.replace("4.0 -9999 6.0", "4.0 -9999"))


def test_parse_rejects_missing_header():
    with pytest.raises(GridError, match="missing required header"):
        Grid.parse("nrows 2\ncellsize 1\n1 2\n3 4\n")


def test_parse_rejects_non_numeric_cell():
    with pytest.raises(GridError, match="non-numeric"):
        Grid.parse(SAMPLE.replace("2.0 3.0", "2.0 tree"))


def test_geometry():
    grid = Grid.parse(SAMPLE)
    # Row 0 is the northern row.
    assert grid.cell_center(0, 0) == pytest.approx((101.0, 203.0))
    assert grid.cell_center(1, 2) == pytest.approx((105.0, 201.0))
    assert grid.cell_of(101.0, 203.0) == (0, 0)
    assert grid.cell_of(105.0, 201.0) == (1, 2)
    assert grid.cell_of(0.0, 0.0) is None
    assert grid.extent == (100.0, 200.0, 106.0, 204.0)
    assert grid.area == pytest.approx(24.0)
    assert grid.area_ha == pytest.approx(0.0024)


def test_cell_center_round_trips_through_cell_of(tiny_grid):
    for row, col in tiny_grid.cells():
        x, y = tiny_grid.cell_center(row, col)
        assert tiny_grid.cell_of(x, y) == (row, col)


def test_index_bounds(tiny_grid):
    with pytest.raises(IndexError):
        tiny_grid.get(4, 0)
    with pytest.raises(IndexError):
        tiny_grid.cell_center(-1, 0)


def test_stats_skips_nodata():
    grid = Grid.parse(SAMPLE)
    stats = grid.stats()
    assert stats.count == 5
    assert stats.nodata_count == 1
    assert stats.minimum == 1.0
    assert stats.maximum == 6.0
    assert stats.mean == pytest.approx(3.2)
    assert stats.stdev == pytest.approx(math.sqrt(3.7), rel=1e-9)
    assert stats.as_dict()["count"] == 5


def test_stats_of_empty_grid():
    grid = Grid.filled(2, 2, None)
    stats = grid.stats()
    assert stats.count == 0
    assert stats.mean is None


def test_combine_propagates_nodata():
    grid = Grid.parse(SAMPLE)
    doubled = grid.combine(grid, lambda a, b: a + b)
    assert doubled.get(0, 0) == 2.0
    assert doubled.get(1, 1) is None


def test_combine_rejects_misaligned():
    grid = Grid.parse(SAMPLE)
    other = Grid.filled(2, 3, 1.0, cellsize=5.0)
    with pytest.raises(GridError, match="cellsize"):
        grid.combine(other, lambda a, b: a + b)
    with pytest.raises(GridError, match="shape"):
        grid.combine(Grid.filled(3, 3, 1.0, cellsize=2.0), lambda a, b: a + b)


def test_map_and_clip():
    grid = Grid.parse(SAMPLE)
    assert grid.map_values(lambda v: v * 2).get(0, 1) == 4.0
    clipped = grid.clip(minimum=2.0, maximum=5.0)
    assert clipped.get(0, 0) == 2.0
    assert clipped.get(1, 2) == 5.0
    assert clipped.get(1, 1) is None


def test_replace_nodata():
    grid = Grid.parse(SAMPLE).replace_nodata(0.0)
    assert grid.get(1, 1) == 0.0
    assert grid.stats().nodata_count == 0


def test_focal_mean_smooths_the_peak(tiny_grid):
    smoothed = tiny_grid.smooth_mean(1)
    assert smoothed.get(1, 1) < tiny_grid.get(1, 1)
    # The kernel is circular, so radius 1 covers the four edge neighbours only.
    assert smoothed.get(1, 1) == pytest.approx((5.0 + 0.0 + 4.0 + 0.0 + 4.0) / 5.0)


def test_focal_median_removes_a_spike():
    grid = Grid.from_rows([[1.0] * 5 for _ in range(5)])
    grid.set(2, 2, 99.0)
    assert grid.smooth_median(1).get(2, 2) == 1.0


def test_focal_max_and_window_shapes(tiny_grid):
    assert tiny_grid.focal_max(1).get(1, 0) == 5.0
    assert tiny_grid.focal_max(1).get(0, 0) == 0.0  # (1, 1) is outside a circular r=1
    square = set(tiny_grid.window(1, 1, 1))
    circular = set(tiny_grid.window(1, 1, 1, circular=True))
    assert circular < square
    assert (0, 0) in square and (0, 0) not in circular


def test_neighbors_connectivity(tiny_grid):
    assert len(list(tiny_grid.neighbors(1, 1, 4))) == 4
    assert len(list(tiny_grid.neighbors(1, 1, 8))) == 8
    assert len(list(tiny_grid.neighbors(0, 0, 8))) == 3
    with pytest.raises(GridError):
        list(tiny_grid.neighbors(0, 0, 6))


def test_histogram():
    grid = Grid.from_rows([[0.0, 1.0], [2.0, 3.0]])
    edges, counts = grid.histogram(bins=3, lo=0.0, hi=3.0)
    assert len(edges) == 4
    assert sum(counts) == 4


def test_from_rows_validates():
    with pytest.raises(GridError, match="same length"):
        Grid.from_rows([[1.0, 2.0], [3.0]])
    with pytest.raises(GridError, match="zero rows"):
        Grid.from_rows([])


def test_construction_validates():
    with pytest.raises(GridError, match="positive number"):
        Grid(0, 1, 0.0, 0.0, 1.0)
    with pytest.raises(GridError, match="cellsize"):
        Grid(1, 1, 0.0, 0.0, 0.0)
    with pytest.raises(GridError, match="expected 4 values"):
        Grid(2, 2, 0.0, 0.0, 1.0, values=[1.0])


def test_copy_is_independent(tiny_grid):
    clone = tiny_grid.copy()
    clone.set(0, 0, 42.0)
    assert tiny_grid.get(0, 0) == 0.0


def test_item_access(tiny_grid):
    tiny_grid[0, 0] = 9.0
    assert tiny_grid[0, 0] == 9.0


def test_valid_cells_matches_finite_values():
    grid = Grid.parse(SAMPLE)
    assert [value for _, _, value in grid.valid_cells()] == list(grid.finite_values())


def test_serialisation_trims_trailing_zeros():
    text = Grid.from_rows([[1.5, 2.0]]).to_text()
    assert "1.5 2" in text
    assert "NODATA_value -9999" in text
