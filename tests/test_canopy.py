"""Tests for canopy height models, cover and gap mapping."""

from __future__ import annotations

import math

import pytest

from silvispect.canopy import (
    _opened_gap_mask,
    canopy_cover,
    canopy_height_model,
    distance_to_canopy,
    find_gaps,
    gap_fraction,
    opening_radius_field,
    rugosity,
    vertical_strata,
)
from silvispect.grid import Grid, GridError


def test_chm_is_dsm_minus_dtm():
    dsm = Grid.from_rows([[12.0, 11.0], [10.0, 9.0]])
    dtm = Grid.from_rows([[2.0, 2.0], [3.0, 3.0]])
    chm = canopy_height_model(dsm, dtm)
    assert chm.get(0, 0) == 10.0
    assert chm.get(1, 1) == 6.0


def test_chm_clamps_negative_heights():
    dsm = Grid.from_rows([[1.0]])
    dtm = Grid.from_rows([[3.0]])
    assert canopy_height_model(dsm, dtm).get(0, 0) == 0.0


def test_chm_keeps_nodata():
    dsm = Grid.from_rows([[10.0, None]])
    dtm = Grid.from_rows([[1.0, 1.0]])
    assert canopy_height_model(dsm, dtm).get(0, 1) is None


def test_canopy_cover_and_gap_fraction():
    chm = Grid.from_rows([[10.0, 0.5], [3.0, 1.0]])
    assert canopy_cover(chm, 2.0) == pytest.approx(0.5)
    assert gap_fraction(chm, threshold=2.0) == pytest.approx(0.5)
    assert canopy_cover(Grid.filled(2, 2, None)) == 0.0


def test_find_gaps_labels_one_opening():
    rows = [[10.0] * 6 for _ in range(6)]
    for r in (2, 3):
        for c in (2, 3):
            rows[r][c] = 0.0
    chm = Grid.from_rows(rows, cellsize=2.0)  # each cell is 4 m2
    gaps = find_gaps(chm, threshold=2.0, min_area=4.0)
    assert len(gaps) == 1
    gap = gaps[0]
    assert gap.gap_id == 1
    assert len(gap.cells) == 4
    assert gap.area == pytest.approx(16.0)
    assert gap.touches_edge is False
    assert gap.equivalent_radius == pytest.approx((16.0 / 3.141592653589793) ** 0.5)
    assert gap.centroid_x == pytest.approx(6.0)
    assert gap.as_dict()["cell_count"] == 4


def test_find_gaps_respects_min_area():
    rows = [[10.0] * 4 for _ in range(4)]
    rows[1][1] = 0.0
    chm = Grid.from_rows(rows)
    assert find_gaps(chm, min_area=25.0) == []
    assert len(find_gaps(chm, min_area=0.5)) == 1


def test_find_gaps_flags_and_can_drop_edge_gaps():
    rows = [[10.0] * 4 for _ in range(4)]
    rows[0][0] = 0.0
    rows[0][1] = 0.0
    chm = Grid.from_rows(rows)
    gaps = find_gaps(chm, min_area=1.0)
    assert gaps[0].touches_edge is True
    assert find_gaps(chm, min_area=1.0, include_edge_gaps=False) == []


def test_find_gaps_sorted_by_area_and_renumbered():
    rows = [[10.0] * 9 for _ in range(5)]
    rows[2][1] = 0.0
    for c in (5, 6, 7):
        rows[2][c] = 0.0
    gaps = find_gaps(Grid.from_rows(rows), min_area=0.5)
    assert [len(gap.cells) for gap in gaps] == [3, 1]
    assert [gap.gap_id for gap in gaps] == [1, 2]


def test_find_gaps_rejects_negative_min_area():
    with pytest.raises(GridError):
        find_gaps(Grid.filled(2, 2, 0.0), min_area=-1.0)


def test_rugosity_is_zero_for_a_flat_canopy():
    assert rugosity(Grid.filled(3, 3, 20.0)) == 0.0
    assert rugosity(Grid.from_rows([[10.0, 20.0, 30.0]])) == pytest.approx(10.0)


def test_vertical_strata_sum_to_one():
    chm = Grid.from_rows([[0.1, 1.0, 3.0], [10.0, 20.0, 30.0]])
    strata = vertical_strata(chm)
    assert sum(strata.values()) == pytest.approx(1.0)
    assert strata["ground"] == pytest.approx(1 / 6)
    assert strata["emergent"] == pytest.approx(1 / 6)


def test_vertical_strata_requires_ascending_breaks():
    with pytest.raises(GridError):
        vertical_strata(Grid.filled(2, 2, 1.0), breaks=(5.0, 2.0))


def test_vertical_strata_of_empty_grid():
    strata = vertical_strata(Grid.filled(2, 2, None))
    assert set(strata.values()) == {0.0}


def test_synthetic_stand_has_realistic_cover(stand):
    cover = canopy_cover(stand.chm, 2.0)
    assert 0.4 < cover < 1.0
    assert stand.chm.stats().maximum > 15.0


def test_distance_to_canopy_measures_to_the_nearest_tree():
    """The distance is to the canopy's edge, as the edge of the plot already was.

    Part of the step from an open cell's centre to a canopy cell's centre lies
    inside the canopy cell, so it is not clearance.  Counting it made the two
    halves of this field disagree — the plot edge below is measured to the
    extent, half a cell from the outermost centre — and a one-metre opening
    claimed a two-metre inscribed circle.

    How much of that step is inside depends on which way the tree lies.
    Straight across, the boundary is half a cell away.  Diagonally, the nearest
    point of the cell is its *corner*, ``sqrt(2)/2`` of a cell away — taking
    off a flat half cell there left a diagonal boundary a quarter of a cell
    further away than it is.
    """
    chm = Grid.filled(7, 7, 0.0, cellsize=1.0)
    chm.set(3, 3, 10.0)
    distances = distance_to_canopy(chm, threshold=2.0)
    assert distances.get(3, 3) == 0.0
    assert distances.get(3, 4) == pytest.approx(0.5)
    assert distances.get(2, 2) == pytest.approx(math.sqrt(2.0) / 2.0)


def test_distance_to_canopy_is_bounded_by_the_plot_edge():
    """Nothing is known beyond the extent, so the edge is a boundary too."""
    chm = Grid.filled(7, 7, 0.0, cellsize=1.0)
    chm.set(3, 3, 10.0)
    distances = distance_to_canopy(chm, threshold=2.0)
    # The corner is 4.24 m from the only tree but half a cell from the edge.
    assert distances.get(0, 0) == pytest.approx(0.5)
    assert max(v for v in distances.values if v is not None) <= 3.5


def test_treeless_raster_reports_a_finite_gap_width():
    """With no canopy anywhere the plot edge is the only boundary there is."""
    chm = Grid.filled(10, 10, 0.0, cellsize=1.0)
    gaps = find_gaps(chm, threshold=2.0, min_area=1.0, min_width=5.0)
    assert len(gaps) == 1
    assert math.isfinite(gaps[0].width)
    # The largest circle inside a 10 m x 10 m opening has a 10 m diameter; the
    # chamfer transform measures from cell centres, so it lands just under.
    assert 8.0 <= gaps[0].width <= 10.0
    assert math.isfinite(gaps[0].as_dict()["width_m"])


def test_treeless_raster_distance_field_is_finite():
    distances = distance_to_canopy(Grid.filled(6, 6, 0.0, cellsize=1.0), threshold=2.0)
    assert all(value is not None and math.isfinite(value) for value in distances.values)


def test_gap_width_cannot_exceed_the_raster():
    """An inscribed circle must fit inside the plot that was actually mapped."""
    chm = Grid.filled(10, 10, 0.0, cellsize=1.0)
    chm.set(0, 0, 20.0)  # a single corner tree; the rest is open
    diagonal = math.hypot(10.0, 10.0)
    for min_width in (0.0, 5.0):
        gaps = find_gaps(chm, threshold=2.0, min_area=1.0, min_width=min_width)
        assert gaps, f"expected an opening at min_width={min_width}"
        for gap in gaps:
            assert gap.width <= diagonal
            assert gap.touches_edge is True


def test_interior_gap_width_is_unaffected_by_the_edge_bound():
    """A gap far from the border must still measure to the surrounding canopy."""
    chm = Grid.filled(40, 40, 20.0, cellsize=1.0)
    for row in range(16, 24):  # an 8 m x 8 m opening in the middle
        for col in range(16, 24):
            chm.set(row, col, 0.0)
    gap = find_gaps(chm, threshold=2.0, min_area=1.0, min_width=0.0)[0]
    assert gap.touches_edge is False
    assert 7.0 <= gap.width <= 9.0


def test_opening_is_monotone_in_the_requested_width():
    """A stricter width must never enlarge an opening.

    Erosion followed by dilation mixed two quantisations and oscillated, so a
    tighter `min_width` could grow a gap and raise a finding a looser setting
    did not.  The opening function makes the family nested by construction.
    """
    import random

    rng = random.Random(11)
    for _ in range(60):
        chm = Grid.filled(12, 12, 0.0, cellsize=1.0)
        for row, col in chm.cells():
            chm.set(row, col, 20.0 if rng.random() < 0.35 else 0.0)
        previous: set[int] | None = None
        for width in (1.0, 2.0, 3.0, 4.0, 5.0, 6.0):
            mask, _ = _opened_gap_mask(chm, 2.0, width / 2.0)
            current = {i for i, keep in enumerate(mask) if keep}
            if previous is not None:
                assert current <= previous, f"width {width} added cells"
            previous = current


def test_tightening_the_width_never_raises_a_new_area_finding():
    chm = Grid.filled(10, 10, 0.0, cellsize=1.0)
    areas = [
        (find_gaps(chm, threshold=2.0, min_area=1.0, min_width=w) or [None])[0] for w in (5.9, 6.0)
    ]
    assert areas[1] is None or areas[0] is None or areas[1].area <= areas[0].area


def _reflect_horizontally(grid: Grid) -> Grid:
    out = grid.like()
    for row in range(grid.nrows):
        for col in range(grid.ncols):
            out.set(row, col, grid.get(row, grid.ncols - 1 - col))
    return out


def test_opening_covers_every_disc_that_fits():
    """Four equal discs cover a treeless 4x4 plot; none may be dropped.

    Skipping a disc because its centre was already covered is unsound — a disc
    of the same radius centred elsewhere does not contain it — and it made the
    result depend on which cell was visited first.
    """
    chm = Grid.filled(4, 4, 0.0, cellsize=1.0)
    mask, _ = _opened_gap_mask(chm, 2.0, 1.5)
    assert sum(mask) == 16


def test_opening_commutes_with_reflection():
    """Geometry has no preferred direction, so neither may the opening."""
    import random

    rng = random.Random(7)
    for _ in range(40):
        chm = Grid.filled(10, 10, 0.0, cellsize=1.0)
        for row, col in chm.cells():
            chm.set(row, col, 20.0 if rng.random() < 0.3 else 0.0)
        direct, _ = _opened_gap_mask(chm, 2.0, 1.5)
        as_grid = Grid.from_rows(
            [[1.0 if direct[r * 10 + c] else 0.0 for c in range(10)] for r in range(10)]
        )
        reflected_after = [value == 1.0 for value in _reflect_horizontally(as_grid).values]
        reflected_before, _ = _opened_gap_mask(_reflect_horizontally(chm), 2.0, 1.5)
        assert reflected_after == reflected_before


def test_opening_grows_with_the_canopy_threshold():
    """A higher threshold admits more sub-canopy area, so the opening can only grow."""
    import random

    rng = random.Random(3)
    for _ in range(40):
        chm = Grid.filled(10, 10, 0.0, cellsize=1.0)
        for row, col in chm.cells():
            chm.set(row, col, rng.uniform(0.0, 5.0))
        previous: set[int] | None = None
        for threshold in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0):
            mask, _ = _opened_gap_mask(chm, threshold, 1.5)
            current = {i for i, keep in enumerate(mask) if keep}
            if previous is not None:
                assert previous <= current, f"threshold {threshold} lost cells"
            previous = current


def _rotate_quarter(grid: Grid) -> Grid:
    """Rotate a raster 90 degrees, which is an isometry of the cell lattice."""
    out = Grid(
        ncols=grid.nrows,
        nrows=grid.ncols,
        xllcorner=grid.xllcorner,
        yllcorner=grid.yllcorner,
        cellsize=grid.cellsize,
        nodata_value=grid.nodata_value,
    )
    out.values = [None] * (grid.ncols * grid.nrows)
    for row in range(grid.nrows):
        for col in range(grid.ncols):
            out.set(col, grid.nrows - 1 - row, grid.get(row, col))
    return out


def test_opening_field_is_covariant_with_the_cell_size():
    """Rescaling the plot rescales the answer; it does not change its shape.

    The opening field carries map units, so the only fixed quantity a
    comparison against it may use is a fraction of the cell size.  An absolute
    tolerance is invisible on a metre-wide cell and swamps a millimetre-wide
    one: at a cell size of 1e-9 the four equal discs of a treeless 4x4 plot
    collapsed to two distinct radii, and the width filter admitted every cell
    whatever was asked for.
    """
    import random

    rng = random.Random(31)
    heights = [rng.choice([0.0, 0.0, 0.0, 20.0]) for _ in range(36)]
    reference: list[float] | None = None
    for cellsize in (1e-9, 1e-3, 1.0, 10.0, 1e6):
        chm = Grid.filled(6, 6, 0.0, cellsize=cellsize)
        chm.values = list(heights)
        field = opening_radius_field(chm, threshold=2.0)
        shape = [value / cellsize for value in field]
        if reference is None:
            reference = shape
        assert all(
            math.isclose(a, b, rel_tol=1e-9) for a, b in zip(shape, reference, strict=True)
        ), f"cell size {cellsize:g} changed the shape of the field"

        # The width filter has to scale with it, so the same relative width
        # selects the same cells.
        mask, _ = _opened_gap_mask(chm, 2.0, 1.5 * cellsize)
        expected = [value >= 1.5 * cellsize * (1 - 1e-9) for value in field]
        assert mask == [keep and chm.values[i] < 2.0 for i, keep in enumerate(expected)], (
            f"cell size {cellsize:g} changed which cells the width filter kept"
        )


def test_opening_field_is_equivariant_under_rotation():
    """A quarter turn of the plot must be a quarter turn of the field."""
    import random

    rng = random.Random(77)
    for _ in range(40):
        chm = Grid.filled(9, 7, 0.0, cellsize=1.0)
        for row, col in chm.cells():
            chm.set(row, col, 20.0 if rng.random() < 0.3 else 0.0)
        before = opening_radius_field(chm, threshold=2.0)
        field = Grid.from_rows(
            [[before[r * chm.ncols + c] for c in range(chm.ncols)] for r in range(chm.nrows)]
        )
        turned = _rotate_quarter(chm)
        after = opening_radius_field(turned, threshold=2.0)
        assert _rotate_quarter(field).values == after


def test_opening_field_is_the_union_of_every_disc_that_fits():
    """The field is defined by geometry alone, so a brute-force scan must agree.

    This is the invariant behind the reflection and rotation checks: no disc
    that fits inside the opening may be dropped, whatever order the scan visits
    the cells in and whatever bookkeeping it uses to stay cheap.  The oracle
    shares nothing with the code: it samples circle centres densely, keeps
    those whose clearance (measured from the definition) reaches the cell
    centre, and brackets the true value using only that clearance is
    1-Lipschitz.  Every cell must land in its bracket.
    """
    import random

    rng = random.Random(505)
    grid = 16
    slack = math.sqrt(2.0) / (2.0 * grid)
    for _ in range(40):
        nrows, ncols = rng.randint(1, 8), rng.randint(1, 8)
        cellsize = rng.choice([0.5, 1.0, 3.0])
        chm = Grid.filled(nrows, ncols, 0.0, cellsize=cellsize)
        chm.values = [
            None if rng.random() < 0.08 else rng.uniform(0.0, 6.0) for _ in range(nrows * ncols)
        ]
        blocked = [value is None or value >= 2.0 for value in chm.values]
        samples: list[tuple[float, float, float]] = []
        for index, is_blocked in enumerate(blocked):
            if is_blocked:
                continue
            row, col = divmod(index, ncols)
            for i in range(grid + 1):
                for j in range(grid + 1):
                    x, y = col + j / grid, row + i / grid
                    samples.append((x, y, _clearance(blocked, nrows, ncols, x, y)))
        field = opening_radius_field(chm, threshold=2.0)
        for index, is_blocked in enumerate(blocked):
            if is_blocked:
                assert field[index] == 0.0
                continue
            row, col = divmod(index, ncols)
            px, py = col + 0.5, row + 0.5
            # Any circle that covers the centre has a sample within half a
            # sample diagonal of its own centre; that sample's clearance is
            # at most that much smaller and reaches at most that much less far.
            low = max((d for x, y, d in samples if math.hypot(x - px, y - py) <= d), default=0.0)
            high = (
                max(
                    (d for x, y, d in samples if math.hypot(x - px, y - py) <= d + 2.0 * slack),
                    default=0.0,
                )
                + slack
            )
            value = field[index] / cellsize
            assert low - 1e-9 <= value <= high + 1e-9, (chm.values, index, value, low, high)


@pytest.mark.parametrize("nrows,ncols", [(1, 1), (2, 2), (3, 3), (4, 6), (5, 5), (3, 8), (7, 2)])
def test_opening_field_of_a_treeless_plot_in_closed_form(nrows, ncols):
    """On an empty rectangle the largest circle covering a cell is known exactly.

    It is either the circle filling the short side — which can slide along
    the long axis, so it covers everything within its radius of that run — or
    a circle wedged into the nearest corner and passing through the cell
    centre, whose radius solves ``(r - a)^2 + (r - b)^2 = r^2`` for the cell's
    distances ``a`` and ``b`` to the two walls of that corner.
    """
    chm = Grid.filled(nrows, ncols, 0.0, cellsize=1.0)
    field = opening_radius_field(chm, threshold=2.0)
    half = min(nrows, ncols) / 2.0
    for row in range(nrows):
        for col in range(ncols):
            px, py = col + 0.5, row + 0.5
            # The sliding circle's run of centres.
            run_x = (max(half, px), min(ncols - half, px)) if ncols > nrows else (half, half)
            run_y = (max(half, py), min(nrows - half, py)) if nrows > ncols else (half, half)
            run_x = (min(run_x), max(run_x))
            run_y = (min(run_y), max(run_y))
            near_x = min(max(px, half), ncols - half)
            near_y = min(max(py, half), nrows - half)
            expected = half if math.hypot(px - near_x, py - near_y) <= half + 1e-12 else 0.0
            for a in (px, ncols - px):
                for b in (py, nrows - py):
                    for root in (a + b - math.sqrt(2 * a * b), a + b + math.sqrt(2 * a * b)):
                        if root <= half + 1e-12:
                            expected = max(expected, root)
            assert field[row * ncols + col] == pytest.approx(expected, rel=1e-12), (row, col)


def test_opening_field_cost_grows_with_the_raster_not_the_discs():
    """Painting every maximal disc cell by cell is cubic in the plot size.

    A treeless plot is the worst case: every cell carries a disc and the discs
    are as large as the plot.  Dropping discs contained in a neighbour's and
    settling each surviving cell once keeps the work proportional to the number
    of cells, so sixteen times the cells must not cost sixteen times over: the
    superlinear version ran 61x here where this one runs 17x.
    """
    import time

    def elapsed(side: int) -> float:
        chm = Grid.filled(side, side, 0.0, cellsize=1.0)
        start = time.perf_counter()
        opening_radius_field(chm, threshold=2.0)
        return time.perf_counter() - start

    small = max(elapsed(32), 1e-4)
    large = elapsed(128)
    assert large / small < 32.0, f"sixteen times the cells cost {large / small:.1f}x"


def test_edge_contact_follows_the_requested_connectivity():
    """Under four-connectivity a corner touch does not join two gaps.

    Edge contact is inherited from the untrimmed sub-canopy component, so that
    component has to be built with the connectivity the caller asked for.  Built
    with eight regardless, a diagonal contact carried the edge flag inward and
    ``include_edge_gaps=False`` discarded an interior opening.
    """
    chm = Grid.from_rows([[0.0, 20.0, 20.0], [20.0, 0.0, 20.0], [20.0, 20.0, 20.0]])
    four = find_gaps(chm, threshold=2.0, min_area=0.0, min_width=0.0, connectivity=4)
    assert [gap.touches_edge for gap in four] == [True, False]
    assert (
        len(
            find_gaps(
                chm,
                threshold=2.0,
                min_area=0.0,
                min_width=0.0,
                connectivity=4,
                include_edge_gaps=False,
            )
        )
        == 1
    )

    eight = find_gaps(chm, threshold=2.0, min_area=0.0, min_width=0.0, connectivity=8)
    assert [gap.touches_edge for gap in eight] == [True]


def test_edge_flag_agrees_with_the_component_it_labels():
    """Whatever the connectivity, a gap is an edge gap exactly when it has an edge cell.

    With ``min_width=0`` nothing is trimmed away, so the reported flag and the
    cells of the gap have to tell the same story — for every connectivity, on
    every raster.
    """
    import random

    rng = random.Random(88)
    for _ in range(60):
        nrows, ncols = rng.randint(2, 8), rng.randint(2, 8)
        chm = Grid.filled(nrows, ncols, 0.0, cellsize=1.0)
        for row, col in chm.cells():
            chm.set(row, col, 20.0 if rng.random() < 0.45 else 0.0)
        for connectivity in (4, 8):
            gaps = find_gaps(
                chm,
                threshold=2.0,
                min_area=0.0,
                min_width=0.0,
                connectivity=connectivity,
            )
            for gap in gaps:
                on_border = any(
                    row in (0, nrows - 1) or col in (0, ncols - 1) for row, col in gap.cells
                )
                assert gap.touches_edge == on_border, (
                    f"connectivity {connectivity}: gap {gap.gap_id} flagged "
                    f"{gap.touches_edge} with border cells {on_border}"
                )
            kept = find_gaps(
                chm,
                threshold=2.0,
                min_area=0.0,
                min_width=0.0,
                connectivity=connectivity,
                include_edge_gaps=False,
            )
            assert len(kept) == sum(not gap.touches_edge for gap in gaps)


def test_gap_fraction_is_exact_at_a_round_share():
    """Three open cells in ten is three tenths, not a hair more.

    Counted directly the answer is exactly representable; taken as one minus
    the cover it is ``0.30000000000000004``, and a rule phrased as *above* a
    30 % limit then fires on a plot sitting exactly on it.
    """
    chm = Grid.from_rows([[0.0] * 3 + [20.0] * 7])
    assert gap_fraction(chm, threshold=2.0) == 0.3
    assert gap_fraction(chm, threshold=2.0) + canopy_cover(chm, 2.0) == 1.0

    for open_cells in range(0, 11):
        rows = [[0.0] * open_cells + [20.0] * (10 - open_cells)]
        assert gap_fraction(Grid.from_rows(rows), threshold=2.0) == open_cells / 10

    assert gap_fraction(Grid.filled(2, 2, None), threshold=2.0) == 0.0


def test_rugosity_stays_finite_for_a_raster_it_accepts():
    """The deviations are taken in units of the largest height present.

    Summing the squares first overflowed for any canopy above about 1.3e154 m
    — an absurd forest, but a raster the reader accepts, writes back out and
    describes everywhere else without complaint.
    """
    import sys

    peak = sys.float_info.max
    flat = Grid.from_rows([[peak, peak, peak]], cellsize=1.0)
    assert rugosity(flat, threshold=0.0) == 0.0

    spread = Grid.from_rows([[peak / 4, -peak / 4, peak / 4, peak / 4]], cellsize=1.0)
    value = rugosity(spread, threshold=-peak)
    assert math.isfinite(value) and value > 0.0


def test_gap_geometry_stays_finite_for_a_raster_it_accepts():
    """A gap centroid is an average, and averages were summed before dividing.

    The cell centres here are each comfortably finite; only their running
    total is not, which is exactly the case a plain ``sum`` gets wrong.  The
    centroid it produced was ``Infinity``, and that went straight into a JSON
    document.
    """
    grid = Grid.from_rows([[0.0] * 12 for _ in range(12)], cellsize=1e306)
    gaps = find_gaps(grid, threshold=1.0, min_area=0.0, min_width=0.0)
    assert gaps
    for gap in gaps:
        assert math.isfinite(gap.centroid_x)
        assert math.isfinite(gap.centroid_y)
        assert math.isfinite(gap.mean_height)


def test_gap_width_is_the_circle_that_fits_inside_the_opening():
    """One open cell on a one-metre raster is one metre wide, not two.

    The distance field walked centre to centre, so the open cell counted a
    whole cell of clearance when half of that step lies inside the tree beside
    it.  The opening then survived a ``min_width`` of 1.5 m that no circle
    inside it could have met.
    """
    plot = Grid.from_rows([[9.0, 9.0, 9.0], [9.0, 0.0, 9.0], [9.0, 9.0, 9.0]], cellsize=1.0)
    unfiltered = find_gaps(plot, threshold=2.0, min_area=0.0, min_width=0.0)
    assert [gap.width for gap in unfiltered] == [1.0]
    assert find_gaps(plot, threshold=2.0, min_area=0.0, min_width=1.5) == []
    assert len(find_gaps(plot, threshold=2.0, min_area=0.0, min_width=1.0)) == 1


@pytest.mark.parametrize("cellsize", [0.25, 1.0, 4.0])
def test_a_square_opening_reports_the_side_it_actually_has(cellsize):
    """A k-by-k open square has an inscribed circle of diameter k cells."""
    for span in (1, 3, 5):
        size = span + 4
        rows = [[9.0] * size for _ in range(size)]
        edge = (size - span) // 2
        for row in range(edge, edge + span):
            for col in range(edge, edge + span):
                rows[row][col] = 0.0
        gaps = find_gaps(
            Grid.from_rows(rows, cellsize=cellsize),
            threshold=2.0,
            min_area=0.0,
            min_width=0.0,
        )
        assert len(gaps) == 1
        assert gaps[0].width == pytest.approx(span * cellsize)


def test_a_canopy_height_that_cannot_be_written_is_refused_before_writing():
    """The reader forbids non-finite cells, so the writer must not produce one.

    A surface at the largest finite float over a terrain at its negative has a
    height no float can hold.  The subtraction returned infinity, the file was
    written with ``inf`` in it, the command reported success, and the model it
    had just written would not parse.
    """
    import sys

    peak = sys.float_info.max
    with pytest.raises(GridError, match="not a representable number"):
        canopy_height_model(Grid.from_rows([[peak]]), Grid.from_rows([[-peak]]))

    # Whatever the model does produce reads back as itself.
    derived = canopy_height_model(
        Grid.from_rows([[20.0, 5.0], [1e200, -1e200]]),
        Grid.from_rows([[1.0, 6.0], [1.0, 1.0]]),
    )
    assert Grid.parse(derived.to_text(precision=None)).values == derived.values


def test_a_diagonal_boundary_is_a_corner_away_not_a_side_away():
    """The reported case: a plus-shaped opening, walled in on the diagonals.

    The centre of a plus has open cells straight across from it and canopy on
    every diagonal, so the nearest boundary is a *corner* — ``sqrt(2)/2`` of a
    cell, not the half cell a flat correction assumed.  It reported 1.83 m of
    clearance where the largest circle that fits is 1.41 m across, and survived
    a 1.6 m minimum it should have failed.
    """
    rows = [[9.0] * 5 for _ in range(5)]
    for row, col in ((1, 2), (2, 1), (2, 2), (2, 3), (3, 2)):
        rows[row][col] = 0.0
    plot = Grid.from_rows(rows, cellsize=1.0)

    gaps = find_gaps(plot, threshold=2.0, min_area=0.0, min_width=0.0)
    assert [gap.width for gap in gaps] == [pytest.approx(math.sqrt(2.0))]
    assert find_gaps(plot, threshold=2.0, min_area=0.0, min_width=1.6) == []
    assert len(find_gaps(plot, threshold=2.0, min_area=0.0, min_width=1.4)) == 1


@pytest.mark.parametrize("cellsize", [0.25, 1.0, 4.0])
def test_openings_report_the_circle_that_fits_in_them(cellsize):
    """One shape per adjacency: orthogonal, diagonal, corner, edge, and wider.

    The expected widths are worked out from the geometry rather than read back
    from the implementation: the largest circle inside an opening touches the
    nearest canopy *boundary*, which is half a cell away straight across and
    ``sqrt(2)/2`` of a cell away on the diagonal.
    """
    half, corner = 0.5, math.sqrt(2.0) / 2.0

    def plot(rows):
        return Grid.from_rows([[9.0 if v else 0.0 for v in row] for row in rows], cellsize=cellsize)

    def width_of(rows):
        gaps = find_gaps(plot(rows), threshold=2.0, min_area=0.0, min_width=0.0)
        assert len(gaps) == 1
        return gaps[0].width / cellsize

    # A single cell: four orthogonal walls, so half a cell of clearance.
    assert width_of([[1, 1, 1], [1, 0, 1], [1, 1, 1]]) == pytest.approx(2 * half)
    # A plus: the walls it has are all diagonal.
    assert width_of(
        [[1, 1, 1, 1, 1], [1, 1, 0, 1, 1], [1, 0, 0, 0, 1], [1, 1, 0, 1, 1], [1, 1, 1, 1, 1]]
    ) == pytest.approx(2 * corner)
    # A three-by-three block: the centre is a cell and a half from the wall.
    assert width_of(
        [[1, 1, 1, 1, 1], [1, 0, 0, 0, 1], [1, 0, 0, 0, 1], [1, 0, 0, 0, 1], [1, 1, 1, 1, 1]]
    ) == pytest.approx(3.0)
    # Against the plot edge, which bounds the circle the same way.
    assert width_of([[0, 1, 1], [1, 1, 1], [1, 1, 1]]) == pytest.approx(2 * half)


def test_gap_width_transforms_with_the_plot():
    """A gap is a shape on the ground, so reflecting the plot cannot resize it."""
    rows = [[9.0] * 6 for _ in range(6)]
    for row, col in ((1, 1), (1, 2), (2, 2), (3, 3), (3, 4), (4, 4)):
        rows[row][col] = 0.0
    forward = find_gaps(
        Grid.from_rows(rows, cellsize=0.5), threshold=2.0, min_area=0.0, min_width=0.0
    )
    mirrored = find_gaps(
        Grid.from_rows([row[::-1] for row in rows], cellsize=0.5),
        threshold=2.0,
        min_area=0.0,
        min_width=0.0,
    )
    assert sorted(gap.width for gap in forward) == sorted(gap.width for gap in mirrored)


@pytest.mark.parametrize("floor", [math.inf, -math.inf, math.nan])
def test_a_clamp_that_is_not_a_number_is_refused(floor):
    """The subtraction was checked; the clamp applied afterwards was not.

    A floor of ``inf`` replaced every measured height with infinity, and the
    command wrote that model out and reported success — after which its own
    reader refused the file.
    """
    with pytest.raises(GridError, match="finite"):
        canopy_height_model(Grid.from_rows([[5.0]]), Grid.from_rows([[1.0]]), floor=floor)
    with pytest.raises(GridError, match="finite"):
        Grid.from_rows([[5.0]]).clip(maximum=floor)
    assert canopy_height_model(Grid.from_rows([[5.0]]), Grid.from_rows([[1.0]])).values == [4.0]


def test_a_raster_cannot_hold_a_value_a_file_cannot():
    """``parse`` has always refused non-finite cells; construction now agrees.

    Accepting them in memory left a raster that could be built and analysed but
    not written down and read back, which is the shape every one of these
    defects has taken.
    """
    for bad in (math.inf, -math.inf, math.nan):
        with pytest.raises(GridError, match="finite"):
            Grid.from_rows([[bad]])
        with pytest.raises(GridError, match="finite"):
            Grid.filled(1, 1, bad)
    assert Grid.from_rows([[1.0, None]]).values == [1.0, None]


# ----------------------------------------------------------------------
# the width is the largest circle that fits, exactly
# ----------------------------------------------------------------------
def _clearance(blocked, nrows, ncols, x, y):
    """Independent oracle: Euclidean distance to the nearest blocked square or edge.

    Written from the definition — the nearest point of each canopy square, and
    the plot boundary — with no transform, walk or envelope in it.
    """
    best = min(x, ncols - x, y, nrows - y)
    for row in range(nrows):
        for col in range(ncols):
            if not blocked[row * ncols + col]:
                continue
            dx = max(0.0, col - x, x - (col + 1))
            dy = max(0.0, row - y, y - (row + 1))
            best = min(best, math.hypot(dx, dy))
    return best


def _sampled_radius(blocked, nrows, ncols, members, grid=24):
    """A lower bound on the largest inscribed circle, and how far it can be off.

    The clearance is sampled on a fine lattice over the member cells.  It is a
    1-Lipschitz function, so the true maximum lies within half a lattice
    diagonal of the best sample: the pair returned brackets it.
    """
    best = 0.0
    for index in members:
        row, col = divmod(index, ncols)
        for i in range(grid + 1):
            for j in range(grid + 1):
                best = max(best, _clearance(blocked, nrows, ncols, col + j / grid, row + i / grid))
    return best, best + math.sqrt(2.0) / (2.0 * grid)


def _plot(rows, cellsize=1.0):
    return Grid.from_rows([[9.0 if v else 0.0 for v in row] for row in rows], cellsize=cellsize)


def _box(height, width, pad=1):
    rows = [[1] * (width + 2 * pad) for _ in range(height + 2 * pad)]
    for row in range(pad, pad + height):
        for col in range(pad, pad + width):
            rows[row][col] = 0
    return rows


def _plus(arm):
    size = 2 * arm + 1 + 2
    rows = [[1] * size for _ in range(size)]
    middle = size // 2
    for offset in range(-arm, arm + 1):
        rows[middle][middle + offset] = 0
        rows[middle + offset][middle] = 0
    return rows


def _cross_of_squares():
    """A 3x3 block with one cell added to the middle of each side."""
    rows = [[1] * 7 for _ in range(7)]
    for row in range(2, 5):
        for col in range(2, 5):
            rows[row][col] = 0
    for row, col in ((1, 3), (5, 3), (3, 1), (3, 5)):
        rows[row][col] = 0
    return rows


#: Openings whose largest inscribed circle is known in closed form, in cells.
ANALYTIC = [
    ("1x1", _box(1, 1), 1.0),
    ("2x2", _box(2, 2), 2.0),
    ("3x3", _box(3, 3), 3.0),
    ("4x4", _box(4, 4), 4.0),
    ("5x5", _box(5, 5), 5.0),
    ("2x5", _box(2, 5), 2.0),
    ("3x6", _box(3, 6), 3.0),
    ("4x7", _box(4, 7), 4.0),
    ("1x7 corridor", _box(1, 7), 1.0),
    ("2x9 corridor", _box(2, 9), 2.0),
    ("plus of single cells", _plus(1), math.sqrt(2.0)),
    ("plus of arms of two", _plus(2), math.sqrt(2.0)),
    ("3x3 with a cell on each side", _cross_of_squares(), math.sqrt(10.0)),
    (
        "L of three cells",
        [[1, 1, 1, 1], [1, 0, 0, 1], [1, 0, 1, 1], [1, 1, 1, 1]],
        4.0 - 2.0 * math.sqrt(2.0),
    ),
    ("single cell in the plot corner", [[0, 1, 1], [1, 1, 1], [1, 1, 1]], 1.0),
    ("2x2 in the plot corner", [[0, 0, 1], [0, 0, 1], [1, 1, 1]], 2.0),
    ("2x2 against the plot edge", [[1, 0, 0, 1], [1, 0, 0, 1], [1, 1, 1, 1]], 2.0),
    ("two cells touching at a corner", [[0, 1], [1, 0]], 1.0),
]


@pytest.mark.parametrize("cellsize", [0.25, 1.0, 4.0])
@pytest.mark.parametrize("name,rows,diameter", ANALYTIC, ids=[case[0] for case in ANALYTIC])
def test_gap_width_is_the_largest_inscribed_circle_exactly(name, rows, diameter, cellsize):
    """Squares, rectangles, corridors, crosses and corners, by closed form.

    Two things about that circle are easy to get wrong on a grid, and both
    were.  Its centre need not be a cell centre — the circle that fills a 2x2
    block sits on the corner where the four cells meet, and read at cell
    centres the block was one cell wide.  And its radius is a Euclidean
    distance, not a walk — a chamfer measures a knight's move as
    ``1 + sqrt(2)`` where the crow flies ``sqrt(5)``, so a 3x3 block with a
    cell on each side claimed a 3.41-cell circle where 3.16 is the largest that
    fits.  The expected values here are geometry, not the implementation.
    """
    plot = _plot(rows, cellsize)
    gaps = find_gaps(plot, threshold=2.0, min_area=0.0, min_width=0.0)
    assert len(gaps) == 1, name
    assert gaps[0].width == pytest.approx(diameter * cellsize, rel=1e-9)

    # The filter decides on the same circle: just under keeps, just over drops.
    below = find_gaps(plot, threshold=2.0, min_area=0.0, min_width=diameter * cellsize * (1 - 1e-6))
    above = find_gaps(plot, threshold=2.0, min_area=0.0, min_width=diameter * cellsize * (1 + 1e-6))
    assert len(below) == 1, name
    assert above == [], name


def test_an_even_sided_opening_keeps_its_whole_width():
    """The reported case: a 2x2 block is a two-metre square.

    Its largest circle is centred where the four cells meet, as far from every
    cell centre as a point can be; the field read at cell centres alone gave
    it half its width and a 1.5 m filter rejected it.
    """
    plot = _plot(_box(2, 2))
    assert [gap.width for gap in find_gaps(plot, threshold=2.0, min_area=0.0, min_width=0.0)] == [
        pytest.approx(2.0)
    ]
    kept = find_gaps(plot, threshold=2.0, min_area=0.0, min_width=1.5)
    assert len(kept) == 1 and len(kept[0].cells) == 4

    # A corridor of even width is retained in full at exactly its width: the
    # circle slides along the midline, and every position is painted.
    corridor = _plot(_box(2, 9))
    kept = find_gaps(corridor, threshold=2.0, min_area=0.0, min_width=2.0)
    assert len(kept) == 1 and len(kept[0].cells) == 18


def test_an_oblique_boundary_is_measured_as_the_crow_flies():
    """The reported case: a 3x3 block with one cell added to each side.

    From the centre the nearest walls are the corners of the cells diagonal to
    the arms, a knight's move away: ``sqrt(2.5)`` cells, for a circle of
    diameter ``sqrt(10)``.  The chamfer walk reached those corners in a
    straight step and a diagonal one, ``1 + sqrt(2)/2``, and reported 3.41 m.
    """
    plot = _plot(_cross_of_squares())
    (gap,) = find_gaps(plot, threshold=2.0, min_area=0.0, min_width=0.0)
    assert gap.width == pytest.approx(math.sqrt(10.0))
    assert find_gaps(plot, threshold=2.0, min_area=0.0, min_width=3.3) == []
    assert len(find_gaps(plot, threshold=2.0, min_area=0.0, min_width=3.1)) == 1


def test_the_distance_field_is_exact_euclidean_at_every_cell():
    """No octile approximation anywhere: every cell agrees with the definition."""
    import random

    rng = random.Random(1201)
    for _ in range(40):
        nrows, ncols = rng.randint(1, 9), rng.randint(1, 9)
        cellsize = rng.choice([0.5, 1.0, 2.0])
        chm = Grid.filled(nrows, ncols, 0.0, cellsize=cellsize)
        chm.values = [
            None if rng.random() < 0.05 else rng.choice([0.0, 0.0, 9.0])
            for _ in range(nrows * ncols)
        ]
        blocked = [value is None or value >= 2.0 for value in chm.values]
        field = distance_to_canopy(chm, threshold=2.0)
        for index in range(len(chm.values)):
            if blocked[index]:
                continue
            row, col = divmod(index, ncols)
            expected = _clearance(blocked, nrows, ncols, col + 0.5, row + 0.5) * cellsize
            assert field.values[index] == pytest.approx(expected, rel=1e-12, abs=1e-12)


def test_gap_width_on_irregular_openings_is_certified_and_bounded():
    """Random unions of cells, against an oracle that shares nothing with the code.

    The circle the implementation reports is checked two ways.  Its centre is
    handed to the definition, which must find exactly that much clearance
    there — so the circle really fits.  And the clearance is sampled densely
    over the gap; being 1-Lipschitz, the true maximum lies within half a
    sample diagonal of the best sample, so the report must fall in that
    bracket — so no larger circle was missed.
    """
    import random

    from silvispect.geometry import Walls, distance_field, inscribed_circle

    rng = random.Random(4711)
    checked = 0
    for _ in range(50):
        nrows, ncols = rng.randint(2, 7), rng.randint(2, 7)
        rows = [[1 if rng.random() < 0.35 else 0 for _ in range(ncols)] for _ in range(nrows)]
        blocked = [bool(v) for row in rows for v in row]
        if all(blocked):
            continue
        plot = _plot(rows)
        cells = distance_field(blocked, nrows, ncols)
        walls = Walls(blocked, nrows, ncols)
        for gap in find_gaps(plot, threshold=2.0, min_area=0.0, min_width=0.0):
            members = [r * ncols + c for r, c in gap.cells]
            disc = inscribed_circle(walls, cells, members)
            x, y = (disc.x1 + disc.x2) / 2.0, (disc.y1 + disc.y2) / 2.0
            assert _clearance(blocked, nrows, ncols, x, y) == pytest.approx(disc.radius, abs=1e-8)
            assert gap.width == pytest.approx(2.0 * disc.radius)
            low, high = _sampled_radius(blocked, nrows, ncols, members)
            assert low - 1e-9 <= disc.radius <= high + 1e-9, (rows, gap.cells)
            checked += 1
    assert checked > 60


def test_gap_widths_transform_with_the_plot():
    """A gap is a shape on the ground; turning the plot cannot resize it."""
    import random

    rng = random.Random(2718)
    for _ in range(30):
        nrows, ncols = rng.randint(2, 7), rng.randint(2, 7)
        rows = [[rng.choice([0.0, 0.0, 9.0]) for _ in range(ncols)] for _ in range(nrows)]
        images = [
            rows,
            [row[::-1] for row in rows],
            rows[::-1],
            [list(col) for col in zip(*rows, strict=True)],
            [list(col) for col in zip(*rows[::-1], strict=True)],
        ]
        widths = [
            sorted(
                round(gap.width, 9)
                for gap in find_gaps(
                    Grid.from_rows(image, cellsize=0.5), threshold=2.0, min_area=0.0, min_width=0.0
                )
            )
            for image in images
        ]
        assert all(w == widths[0] for w in widths), rows
        # And the filter's decisions, on the same geometry.
        kept = [
            len(
                find_gaps(
                    Grid.from_rows(image, cellsize=0.5), threshold=2.0, min_area=0.0, min_width=1.0
                )
            )
            for image in images
        ]
        assert all(k == kept[0] for k in kept), rows
