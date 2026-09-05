"""Tests for treetop detection and crown delineation."""

from __future__ import annotations

import math
import random

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

    # Ordinary detection is unchanged: the peak is still found, and every cell
    # that can see it inside its *own* search window is still suppressed.  A
    # 5 m cell searches one cell at this resolution, so that is the four cells
    # orthogonally touching the peak; the surrounding flat field is one plateau
    # and contributes one apex of its own, as it did before.
    grid = Grid.filled(9, 9, 5.0, cellsize=1.0)
    grid.set(4, 4, 30.0)
    found = {(top.row, top.col) for top in find_treetops(grid)}
    assert (4, 4) in found
    assert len(found) == 2
    dominated = {(3, 4), (4, 3), (4, 5), (5, 4)}
    assert found & dominated == set()


# ----------------------------------------------------------------------
# equivariance under the symmetries of a raster
# ----------------------------------------------------------------------
def _rows_of(grid):
    return [[grid.get(row, col) for col in range(grid.ncols)] for row in range(grid.nrows)]


#: The eight symmetries of a rectangle, each as a way of moving the rows of a
#: raster and the matching way of moving one of its cells.  Naming both halves
#: lets a test state its expectation as the *transformed output* of the
#: original run rather than as a second table of constants.
SYMMETRIES = {
    "fliplr": (
        lambda rows: [row[::-1] for row in rows],
        lambda r, c, nr, nc: (r, nc - 1 - c),
    ),
    "flipud": (
        lambda rows: rows[::-1],
        lambda r, c, nr, nc: (nr - 1 - r, c),
    ),
    "rot180": (
        lambda rows: [row[::-1] for row in rows[::-1]],
        lambda r, c, nr, nc: (nr - 1 - r, nc - 1 - c),
    ),
    "rot90": (
        lambda rows: [list(col) for col in zip(*rows[::-1], strict=True)],
        lambda r, c, nr, nc: (c, nr - 1 - r),
    ),
    "rot270": (
        lambda rows: [list(col) for col in zip(*rows, strict=True)][::-1],
        lambda r, c, nr, nc: (nc - 1 - c, r),
    ),
    "transpose": (
        lambda rows: [list(col) for col in zip(*rows, strict=True)],
        lambda r, c, nr, nc: (c, r),
    ),
    "antitranspose": (
        lambda rows: [list(col) for col in zip(*[row[::-1] for row in rows[::-1]], strict=True)],
        lambda r, c, nr, nc: (nc - 1 - c, nr - 1 - r),
    ),
}


def moved(grid, name):
    """Return ``grid`` under one symmetry, and the map that follows a cell."""
    rows = _rows_of(grid)
    move_rows, move_cell = SYMMETRIES[name]
    nrows, ncols = grid.nrows, grid.ncols
    turned = Grid.from_rows(move_rows(rows), cellsize=grid.cellsize)

    def follow(cell):
        return move_cell(cell[0], cell[1], nrows, ncols)

    # The two halves have to describe the same motion, or the test proves
    # nothing; every value must be found where the cell map says it went.
    for row in range(nrows):
        for col in range(ncols):
            new_row, new_col = follow((row, col))
            assert turned.get(new_row, new_col) == grid.get(row, col)
    return turned, follow


def _self_symmetric(grid):
    """True when some symmetry other than the identity leaves the raster alone.

    Such a raster has two cells that no property of the plot can tell apart:
    ``5 5 5 5`` is its own mirror image, and a transpose-symmetric plot has a
    pair of corners that swap.  Whichever of them is chosen, the mirror of that
    choice is the other one, so an apex — and the crown grown from it — cannot
    be equivariant however the tie is broken.  The count and the crown sizes
    still are, and those are what is asserted for these rasters.
    """
    return any(moved(grid, name)[0].values == grid.values for name in SYMMETRIES)


def test_a_mirrored_raster_is_not_a_different_forest():
    """The reported case: one raster, its mirror image, and different answers.

    ``8 7 4 8`` and its reflection ``8 4 7 8`` are the same four trees seen
    from the other side of the plot.  Suppressing a cell that merely *equalled*
    a neighbour made the left-hand 8 lose to the right-hand one in one
    direction and survive in the other, so the mirror held a crown the
    original did not — and the stocking findings came out reversed.
    """
    config = DetectionConfig(smooth_radius=0, min_height=0.0, min_crown_cells=1)
    grid = Grid.from_rows([[8.0, 7.0, 4.0, 8.0]], cellsize=0.25)
    mirror, follow = moved(grid, "fliplr")

    original = detect_trees(grid, config)
    reflected = detect_trees(mirror, config)
    assert len(original.crowns) == len(reflected.crowns) == 2
    assert {follow(cell) for crown in original.crowns for cell in crown.cells} == {
        cell for crown in reflected.crowns for cell in crown.cells
    }


CURATED = [
    [[8.0, 7.0, 4.0, 8.0]],
    [[8.0, 4.0, 7.0, 8.0]],
    [[5.0, 5.0, 5.0, 5.0]],  # one flat plateau, no apex anywhere in particular
    [[9.0, 3.0, 9.0], [3.0, 3.0, 3.0], [9.0, 3.0, 9.0]],  # four equal maxima
    [[12.0, 12.0, 3.0], [12.0, 12.0, 3.0], [3.0, 3.0, 3.0]],  # a square plateau
    [[20.0, 15.0, 10.0, 5.0]],  # a single descending path
    [[10.0, 9.0, 8.0, 9.0, 10.0]],  # tied maxima with a descent between them
    [[7.0, 7.0], [7.0, 7.0]],  # every cell tied
    [[4.0, 3.0, 3.0], [20.0, 3.0, 8.0], [8.0, 3.0, 8.0], [8.0, 8.0, 8.0]],
    [[6.0, 0.0, 6.0], [0.0, 0.0, 0.0], [6.0, 0.0, 6.0]],  # crowns at the size floor
]


@pytest.mark.parametrize("rows", CURATED, ids=range(len(CURATED)))
@pytest.mark.parametrize("name", sorted(SYMMETRIES))
@pytest.mark.parametrize("cellsize", [0.25, 1.0])
def test_detection_follows_the_raster_through_every_symmetry(rows, name, cellsize):
    """Rotating or reflecting a plot must rotate or reflect the answer.

    Detection has three places where equal values have to be separated — which
    candidate survives, which cell of a plateau is its apex, and which of two
    crowns claims a contested cell — and each of them originally fell back on
    the row and column, which reverse with the raster.  The expectation here is
    the *transformed* result of the original run, so there is no second table
    of constants to keep in step.

    A raster the transform leaves unchanged is the one case where the mapped
    crowns cannot be asked for.  ``5 5 5 5`` is its own mirror image, and its
    plateau has no middle cell to be the apex of; whichever of the two central
    cells is chosen, the mirror image of that choice is the other one.  The
    count and the crown sizes are still fixed, and those are what is asserted
    there.
    """
    config = DetectionConfig(smooth_radius=0, min_height=0.0)
    grid = Grid.from_rows(rows, cellsize=cellsize)
    turned, follow = moved(grid, name)

    original = detect_trees(grid, config)
    after = detect_trees(turned, config)

    assert len(after.crowns) == len(original.crowns)
    assert sorted(crown.height for crown in after.crowns) == sorted(
        crown.height for crown in original.crowns
    )
    assert sorted(len(crown.cells) for crown in after.crowns) == sorted(
        len(crown.cells) for crown in original.crowns
    )
    if turned.values != grid.values:
        assert sorted(
            tuple(sorted(map(follow, crown.cells))) for crown in original.crowns
        ) == sorted(tuple(sorted(crown.cells)) for crown in after.crowns)


@pytest.mark.parametrize("cellsize,smooth", [(0.25, 0), (1.0, 0), (1.0, 1)])
def test_random_rasters_keep_their_tree_count_through_every_symmetry(cellsize, smooth):
    """The count is the claim a report makes, so hold it over many rasters.

    Heights are drawn from a short list on purpose: ties are what the ordering
    rules exist for, and a continuous draw would almost never produce one.
    """
    rng = random.Random(20240607)
    config = DetectionConfig(smooth_radius=smooth)
    levels = (8.0, 7.0, 4.0, 3.0, 20.0)
    for _ in range(120):
        nrows, ncols = rng.randint(1, 4), rng.randint(1, 4)
        rows = [[rng.choice(levels) for _ in range(ncols)] for _ in range(nrows)]
        grid = Grid.from_rows(rows, cellsize=cellsize)
        expected = len(detect_trees(grid, config).crowns)
        for name in SYMMETRIES:
            turned, _ = moved(grid, name)
            assert len(detect_trees(turned, config).crowns) == expected, (rows, name)


def test_equal_neighbours_do_not_delete_each_other():
    """Two cells of the same height are both still treetops.

    The window rejected a cell that any neighbour merely *reached*, so of two
    equal maxima the one scanned second was deleted — which is what made the
    answer depend on the scan direction.  A run of equal cells is one treetop
    because it is one plateau, not because one of them won.
    """
    config = DetectionConfig(smooth_radius=0)
    apart = find_treetops(Grid.from_rows([[9.0, 1.0, 1.0, 1.0, 9.0]], cellsize=1.0), config)
    assert [(top.row, top.col) for top in apart] == [(0, 0), (0, 4)]

    touching = find_treetops(Grid.from_rows([[9.0, 9.0, 1.0]], cellsize=1.0), config)
    assert len(touching) == 1


def test_a_plateau_apex_sits_in_the_middle_of_the_plateau():
    """The representative cell of an equal-height run is its centre.

    A corner is as defensible as any other member until the raster is
    reflected, when the corner moves to the opposite side and the crown grows
    somewhere else.  The centre travels with the run.
    """
    rows = [[0.0] * 7 for _ in range(7)]
    for row in (2, 3, 4):
        for col in (2, 3, 4):
            rows[row][col] = 15.0
    tops = find_treetops(Grid.from_rows(rows, cellsize=1.0))
    assert [(top.row, top.col) for top in tops] == [(3, 3)]


def test_smoothing_may_make_two_maxima_equal_but_not_choose_between_them():
    """The reported case: the mean filter erases the only asymmetry there was.

    ``2 3 / 4 9 / 3 2`` is not its own mirror image, but its smoothed surface
    is — the two maxima it leaves are indistinguishable in every respect.
    Detection was reading its tie-break off that surface, so the flipped plot
    produced the *same* crown rather than the flipped one, and the exported
    positions did not travel with the raster.  Ties are read from the raster as
    measured, which still tells the two apart.
    """
    rows = [[2.0, 3.0], [4.0, 9.0], [3.0, 2.0]]
    grid = Grid.from_rows(rows, cellsize=0.5)
    turned, follow = moved(grid, "flipud")

    # The premise: smoothing really does destroy the difference.
    assert grid.smooth_mean(1).values == turned.smooth_mean(1).values

    original = detect_trees(grid)
    after = detect_trees(turned)
    assert len(original.crowns) == len(after.crowns) == 1
    assert sorted(map(follow, original.crowns[0].cells)) == sorted(after.crowns[0].cells)
    assert follow((original.crowns[0].apex.row, original.crowns[0].apex.col)) == (
        after.crowns[0].apex.row,
        after.crowns[0].apex.col,
    )


@pytest.mark.parametrize("name", sorted(SYMMETRIES))
def test_smoothed_detection_follows_the_raster_through_every_symmetry(name):
    """The default pipeline smooths, so the property has to hold through it."""
    rng = random.Random(20240608)
    levels = (8.0, 7.0, 4.0, 3.0, 20.0)
    for _ in range(150):
        nrows, ncols = rng.randint(2, 5), rng.randint(2, 5)
        rows = [[rng.choice(levels) for _ in range(ncols)] for _ in range(nrows)]
        grid = Grid.from_rows(rows, cellsize=1.0)
        turned, follow = moved(grid, name)
        original = detect_trees(grid)
        after = detect_trees(turned)
        assert len(after.crowns) == len(original.crowns), (rows, name)
        assert sorted(len(c.cells) for c in after.crowns) == sorted(
            len(c.cells) for c in original.crowns
        ), (rows, name)
        if not _self_symmetric(grid):
            assert sorted(
                tuple(sorted(map(follow, crown.cells))) for crown in original.crowns
            ) == sorted(tuple(sorted(crown.cells)) for crown in after.crowns), (rows, name)


def test_a_whole_valued_count_is_stored_as_a_whole_number():
    """``smooth_radius=1.0`` is one cell, so it has to behave as one.

    Validation accepted it as a whole number and then handed the float to
    ``range``, which raised a TypeError from inside detection — the value was
    neither refused nor usable.
    """
    config = DetectionConfig(smooth_radius=1.0, min_crown_cells=3.0)
    assert config.smooth_radius == 1 and isinstance(config.smooth_radius, int)
    assert config.min_crown_cells == 3 and isinstance(config.min_crown_cells, int)
    assert detect_trees(Grid.from_rows([[3.0]], cellsize=1.0), config).crowns == []

    with pytest.raises(GridError, match="whole number"):
        DetectionConfig(smooth_radius=1.5)


@pytest.mark.parametrize("scale", [1.0, 10.0, 1000.0, 0.001])
def test_crown_membership_does_not_depend_on_the_vertical_unit(scale):
    """Measuring the same forest in decimetres must find the same crowns.

    The allowance for a cell that rises by a hair on the way down from an apex
    was a fixed 1e-9 m.  That is generous at one scale and nothing at another,
    so multiplying every height by ten — which cannot change which cell is
    above which — moved a four-cell crown to three.
    """
    config = DetectionConfig(
        smooth_radius=0,
        min_height=0.0,
        min_crown_cells=1,
        window_intercept=10.0,
        window_slope=0.0,
        crown_intercept=10.0,
        crown_slope=0.0,
        drop_fraction=0.0,
    )
    rows = [[3.0 * scale, 2.0 * scale, 1.0 * scale, scale * (1.0 + 5e-10)]]
    crowns = detect_trees(Grid.from_rows(rows, cellsize=1.0), config).crowns
    assert [len(crown.cells) for crown in crowns] == [4]
