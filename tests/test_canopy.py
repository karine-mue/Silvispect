"""Tests for canopy height models, cover and gap mapping."""

from __future__ import annotations

import math

import pytest

from silvispect.canopy import (
    canopy_cover,
    canopy_height_model,
    distance_to_canopy,
    find_gaps,
    gap_fraction,
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
