"""Tests for canopy height models, cover and gap mapping."""

from __future__ import annotations

import math

import pytest

from silvispect.canopy import (
    _canopy_distance_field,
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
    chm = Grid.filled(7, 7, 0.0, cellsize=1.0)
    chm.set(3, 3, 10.0)
    distances = distance_to_canopy(chm, threshold=2.0)
    assert distances.get(3, 3) == 0.0
    assert distances.get(3, 4) == pytest.approx(1.0)
    assert distances.get(2, 2) == pytest.approx(math.sqrt(2.0))


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
        field = Grid.from_rows(
            [
                [
                    opening_radius_field(chm, threshold=2.0)[r * chm.ncols + c]
                    for c in range(chm.ncols)
                ]
                for r in range(chm.nrows)
            ]
        )
        turned = _rotate_quarter(chm)
        after = opening_radius_field(turned, threshold=2.0)
        assert _rotate_quarter(field).values == after


def test_opening_field_is_the_union_of_every_disc_that_fits():
    """The field is defined by geometry alone, so a brute-force scan must agree.

    This is the invariant behind the reflection and rotation checks: no disc
    that fits inside the opening may be dropped, whatever order the scan visits
    the cells in and whatever bookkeeping it uses to stay cheap.
    """
    import random

    rng = random.Random(505)
    for _ in range(60):
        nrows, ncols = rng.randint(1, 9), rng.randint(1, 9)
        chm = Grid.filled(nrows, ncols, 0.0, cellsize=rng.choice([0.5, 1.0, 3.0]))
        chm.values = [
            None if rng.random() < 0.08 else rng.uniform(0.0, 6.0) for _ in range(nrows * ncols)
        ]
        distances = _canopy_distance_field(chm, 2.0)
        expected = [0.0] * len(chm.values)
        for centre, value in enumerate(chm.values):
            if value is None or value >= 2.0 or distances[centre] <= 0.0:
                continue
            crow, ccol = divmod(centre, ncols)
            for position in range(len(chm.values)):
                prow, pcol = divmod(position, ncols)
                offset = math.hypot((prow - crow) * chm.cellsize, (pcol - ccol) * chm.cellsize)
                if offset <= distances[centre] and expected[position] < distances[centre]:
                    expected[position] = distances[centre]
        assert opening_radius_field(chm, threshold=2.0) == expected


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
