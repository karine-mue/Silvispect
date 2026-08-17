"""Tests for canopy height models, cover and gap mapping."""

from __future__ import annotations

import pytest

from silvispect.canopy import (
    canopy_cover,
    canopy_height_model,
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
