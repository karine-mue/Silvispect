"""Tests for treetop detection and crown delineation."""

from __future__ import annotations

import math

import pytest

from silvispect.detect import (
    DetectionConfig,
    crown_radius_limit,
    detect_trees,
    find_treetops,
    segment_crowns,
    window_radius,
)
from silvispect.grid import Grid, GridError
from silvispect.inventory import match_trees


def cone_grid(peaks, *, size=40, cellsize=0.5):
    """Render paraboloid crowns onto an otherwise empty canopy."""
    grid = Grid.filled(size, size, 0.0, cellsize=cellsize)
    for x, y, height, radius in peaks:
        for row, col in grid.cells():
            cx, cy = grid.cell_center(row, col)
            distance = math.hypot(cx - x, cy - y)
            if distance > radius:
                continue
            value = height * (1.0 - (distance / radius) ** 2)
            if value > (grid.get(row, col) or 0.0):
                grid.set(row, col, value)
    return grid


def test_window_radius_grows_with_height():
    config = DetectionConfig()
    assert window_radius(0.0, config) == pytest.approx(config.window_intercept)
    assert window_radius(20.0, config) > window_radius(10.0, config)
    assert crown_radius_limit(20.0, config) > crown_radius_limit(0.0, config)


def test_config_validation():
    with pytest.raises(GridError):
        DetectionConfig(min_height=-1.0)
    with pytest.raises(GridError):
        DetectionConfig(drop_fraction=1.0)
    with pytest.raises(GridError):
        DetectionConfig(smooth_radius=-1)
    with pytest.raises(GridError):
        DetectionConfig(window_intercept=0.0)
    with pytest.raises(GridError):
        DetectionConfig(min_crown_cells=0)
    assert DetectionConfig().as_dict()["min_height"] == 2.0


def test_single_peak_is_found_once():
    grid = cone_grid([(10.0, 10.0, 20.0, 3.0)])
    tops = find_treetops(grid)
    assert len(tops) == 1
    assert tops[0].height == pytest.approx(20.0, abs=0.5)
    assert tops[0].x == pytest.approx(10.0, abs=0.5)


def test_two_separated_peaks_are_both_found():
    grid = cone_grid([(6.0, 6.0, 18.0, 2.5), (14.0, 14.0, 22.0, 3.0)])
    assert len(find_treetops(grid)) == 2


def test_below_min_height_is_ignored():
    grid = cone_grid([(10.0, 10.0, 1.5, 2.0)])
    assert find_treetops(grid, DetectionConfig(min_height=2.0)) == []


def test_plateau_yields_a_single_top():
    rows = [[0.0] * 9 for _ in range(9)]
    for r in (3, 4, 5):
        for c in (3, 4, 5):
            rows[r][c] = 12.0
    tops = find_treetops(Grid.from_rows(rows, cellsize=1.0))
    assert len(tops) == 1


def test_crowns_are_delineated_and_bounded():
    grid = cone_grid([(10.0, 10.0, 20.0, 3.0)])
    result = detect_trees(grid, DetectionConfig(smooth_radius=0))
    assert len(result.crowns) == 1
    crown = result.crowns[0]
    assert crown.area > 0
    assert crown.radius <= crown_radius_limit(crown.height, result.config) + grid.cellsize
    assert crown.max_extent <= crown_radius_limit(crown.height, result.config)
    assert crown.diameter == pytest.approx(2 * crown.radius)
    assert crown.mean_height <= crown.height
    assert crown.as_dict()["tree_id"] == 1


def test_crowns_do_not_overlap():
    grid = cone_grid([(6.0, 6.0, 18.0, 2.5), (13.0, 13.0, 22.0, 3.0)])
    crowns, labels = segment_crowns(grid, find_treetops(grid))
    owned = [cell for crown in crowns for cell in crown.cells]
    assert len(owned) == len(set(owned))
    assigned = [v for v in labels.values if v is not None]
    assert len(assigned) == len(owned)
    assert set(assigned) == {float(crown.tree_id) for crown in crowns}


def test_crowns_are_ranked_by_height():
    grid = cone_grid([(6.0, 6.0, 18.0, 2.5), (14.0, 14.0, 24.0, 3.0)])
    result = detect_trees(grid, DetectionConfig(smooth_radius=0))
    heights = [crown.height for crown in result.crowns]
    assert heights == sorted(heights, reverse=True)
    assert [crown.tree_id for crown in result.crowns] == [1, 2]


def test_tiny_crowns_are_discarded():
    grid = Grid.filled(20, 20, 0.0, cellsize=0.5)
    grid.set(10, 10, 15.0)  # a single-cell spike
    result = detect_trees(grid, DetectionConfig(smooth_radius=0, min_crown_cells=3))
    assert result.crowns == []


def test_density_per_ha_and_dict():
    grid = cone_grid([(10.0, 10.0, 20.0, 3.0)], size=40, cellsize=0.5)
    result = detect_trees(grid)
    assert result.density_per_ha() == pytest.approx(1 / grid.area_ha)
    assert result.density_per_ha(area_ha=0.5) == pytest.approx(2.0)
    with pytest.raises(GridError):
        result.density_per_ha(area_ha=0.0)
    payload = result.as_dict()
    assert payload["tree_count"] == len(result) == 1
    assert len(payload["trees"]) == 1
    assert next(iter(result)) is result.crowns[0]


def test_detection_recovers_a_sparse_synthetic_stand(sparse_stand):
    result = detect_trees(sparse_stand.chm)
    match = match_trees(result.crowns, sparse_stand.trees, tolerance=2.5)
    assert match.recall > 0.9
    assert match.precision > 0.9
    assert match.height_rmse < 1.5
    assert match.mean_offset < 1.0


def test_detection_is_deterministic(stand):
    first = detect_trees(stand.chm)
    second = detect_trees(stand.chm)
    assert [crown.as_dict() for crown in first] == [crown.as_dict() for crown in second]


def test_detection_work_is_bounded_by_the_raster_not_the_window():
    """An implausible but accepted height must not make a tiny raster slow.

    The search window grows with the height of the cell being tested, so a
    height of 1e5 m asks for a radius of 5,501 cells.  On a one-cell raster
    that is one comparison of real work; walking the whole requested square
    spent sixteen seconds on it, and `inspect` runs detection before it can
    report the implausible height at all.
    """
    import time

    def elapsed(height: float) -> float:
        grid = Grid.from_rows([[height]], cellsize=1.0)
        start = time.perf_counter()
        tops = find_treetops(grid)
        assert len(tops) == 1 and tops[0].height == height
        return time.perf_counter() - start

    small = max(elapsed(10.0), 1e-6)
    assert elapsed(1e5) / small < 50.0

    # Ordinary detection is unchanged: a peak is still found, and it still
    # suppresses every cell out to the window its own height buys — a 30 m
    # tree searches 2 cells at this resolution.
    grid = Grid.filled(9, 9, 5.0, cellsize=1.0)
    grid.set(4, 4, 30.0)
    found = {(top.row, top.col) for top in find_treetops(grid)}
    assert (4, 4) in found
    suppressed = {
        (row, col)
        for row in range(2, 7)
        for col in range(2, 7)
        if (row - 4) ** 2 + (col - 4) ** 2 <= 4 and (row, col) != (4, 4)
    }
    assert found & suppressed == set()
