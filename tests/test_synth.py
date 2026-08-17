"""Tests for the synthetic stand generator.

The generator is the ground truth the detector is measured against, so it is
tested as carefully as the code it exercises.
"""

from __future__ import annotations

import pytest

from silvispect.canopy import canopy_height_model
from silvispect.synth import SynthConfig, synthesize


def test_config_validation():
    with pytest.raises(ValueError, match="extent"):
        SynthConfig(width=0.0)
    with pytest.raises(ValueError, match="cellsize"):
        SynthConfig(cellsize=0.0)
    with pytest.raises(ValueError, match="same length"):
        SynthConfig(species=("PIAB",), species_weights=(0.5, 0.5))
    with pytest.raises(ValueError, match="stems_per_ha"):
        SynthConfig(stems_per_ha=-1.0)


def test_area_and_raster_shape():
    stand = synthesize(SynthConfig(width=50.0, height=25.0, cellsize=0.5, seed=1))
    assert stand.area_ha == pytest.approx(0.125)
    assert stand.chm.ncols == 100
    assert stand.chm.nrows == 50
    assert stand.chm.cellsize == 0.5


def test_same_seed_gives_the_same_stand():
    first = synthesize(SynthConfig(width=40.0, height=40.0, seed=99))
    second = synthesize(SynthConfig(width=40.0, height=40.0, seed=99))
    assert first.chm.values == second.chm.values
    assert [t.as_dict() for t in first.trees] == [t.as_dict() for t in second.trees]


def test_different_seeds_give_different_stands():
    first = synthesize(SynthConfig(width=40.0, height=40.0, seed=1))
    second = synthesize(SynthConfig(width=40.0, height=40.0, seed=2))
    assert first.chm.values != second.chm.values


def test_overrides_apply():
    stand = synthesize(SynthConfig(seed=4), width=30.0, height=30.0, gap_count=0)
    assert stand.config.width == 30.0
    assert stand.gaps == []
    assert synthesize(width=30.0, height=30.0).config.width == 30.0


def test_trees_are_inside_the_extent_and_measured():
    stand = synthesize(SynthConfig(width=40.0, height=40.0, seed=8))
    xmin, ymin, xmax, ymax = stand.chm.extent
    for tree in stand.trees:
        assert xmin <= tree.x <= xmax
        assert ymin <= tree.y <= ymax
        assert tree.dbh_cm and tree.dbh_cm > 0
        assert tree.height_m and tree.height_m > 0
        assert tree.crown_diameter_m and tree.crown_diameter_m > 0
        assert tree.species in stand.config.species


def test_dimensions_stay_inside_the_configured_bounds():
    config = SynthConfig(width=60.0, height=60.0, seed=13, dbh_min=10.0, dbh_max=50.0)
    stand = synthesize(config)
    assert all(10.0 <= (t.dbh_cm or 0) <= 50.0 for t in stand.trees)


def test_density_is_close_to_the_target():
    config = SynthConfig(width=100.0, height=100.0, stems_per_ha=300.0, seed=21, gap_count=0)
    stand = synthesize(config)
    assert 0.85 * 300 <= len(stand.trees) / stand.area_ha <= 300


def test_chm_matches_dsm_minus_dtm():
    stand = synthesize(SynthConfig(width=30.0, height=30.0, seed=6))
    recomputed = canopy_height_model(stand.dsm, stand.dtm)
    for expected, actual in zip(stand.chm.values, recomputed.values, strict=True):
        assert expected == pytest.approx(actual, abs=1e-9)


def test_terrain_is_not_flat():
    stand = synthesize(SynthConfig(width=60.0, height=60.0, seed=6))
    assert stand.dtm.stats().stdev > 0.1


def test_canopy_reaches_the_tallest_tree():
    stand = synthesize(SynthConfig(width=50.0, height=50.0, seed=14, surface_noise=0.0))
    tallest = max(tree.height_m or 0 for tree in stand.trees)
    assert stand.chm.stats().maximum == pytest.approx(tallest, rel=0.05)


def test_gaps_are_cut_free_of_trees():
    config = SynthConfig(width=80.0, height=80.0, seed=17, gap_count=2, gap_radius=(8.0, 10.0))
    stand = synthesize(config)
    assert len(stand.gaps) == 2
    for gx, gy, radius in stand.gaps:
        for tree in stand.trees:
            assert ((tree.x - gx) ** 2 + (tree.y - gy) ** 2) ** 0.5 >= radius


def test_zero_density_produces_bare_ground():
    stand = synthesize(SynthConfig(width=20.0, height=20.0, stems_per_ha=0.0, seed=2))
    assert stand.trees == []
    assert stand.chm.stats().maximum == pytest.approx(0.0, abs=0.5)


def test_georeferencing_is_honoured():
    stand = synthesize(
        SynthConfig(width=20.0, height=20.0, seed=3, xllcorner=500_000.0, yllcorner=6_000_000.0)
    )
    assert stand.chm.extent[0] == pytest.approx(500_000.0)
    assert all(tree.x >= 500_000.0 for tree in stand.trees)


def test_extent_rounds_to_whole_cells():
    config = SynthConfig(width=14.0, height=100.0, cellsize=10.0)
    assert config.shape == (10, 1)
    assert config.extent == (10.0, 100.0)
    assert config.area_ha == pytest.approx(0.1)


def test_every_generated_tree_is_inside_the_raster():
    """The returned trees are exactly the stems rendered into the CHM."""
    stand = synthesize(
        SynthConfig(
            width=14.0, height=100.0, cellsize=10.0, stems_per_ha=300.0, seed=4, gap_count=0
        )
    )
    xmin, ymin, xmax, ymax = stand.chm.extent
    assert stand.trees
    for tree in stand.trees:
        assert xmin <= tree.x <= xmax
        assert ymin <= tree.y <= ymax
    assert stand.area_ha == pytest.approx(stand.chm.area_ha)


def test_reported_area_follows_the_raster_not_the_request():
    stand = synthesize(SynthConfig(width=14.0, height=100.0, cellsize=10.0, seed=1))
    assert stand.area_ha == pytest.approx(0.1)


def test_every_truth_tree_paints_a_cell_even_on_coarse_rasters():
    """Every stem reaches its own cell, even when its crown is sub-cell sized."""
    from silvispect.grid import Grid
    from silvispect.synth import _render_crown

    stand = synthesize(
        SynthConfig(
            width=14.0,
            height=100.0,
            cellsize=10.0,
            stems_per_ha=300.0,
            seed=4,
            gap_count=0,
            surface_noise=0.0,
        )
    )
    assert stand.trees
    for tree in stand.trees:
        probe = Grid.filled(
            stand.chm.nrows,
            stand.chm.ncols,
            0.0,
            cellsize=stand.chm.cellsize,
            xllcorner=stand.chm.xllcorner,
            yllcorner=stand.chm.yllcorner,
        )
        _render_crown(probe, tree, stand.config)
        assert any(value for value in probe.values), f"{tree.tree_id} painted nothing"


def test_taller_crowns_occlude_shorter_neighbours_in_the_composite():
    """The surface keeps the maximum, so a stem can be painted yet invisible.

    This is not a generator defect — a canopy height model records the top of
    the canopy, and a suppressed stem leaves no trace in one.  The test pins the
    behaviour so the guarantee is not over-read as "every tree is visible".
    """
    stand = synthesize(
        SynthConfig(
            width=100.0,
            height=100.0,
            cellsize=10.0,
            stems_per_ha=300.0,
            seed=4,
            gap_count=0,
            surface_noise=0.0,
        )
    )
    occupied = sum(1 for value in stand.chm.values if value and value > 0)
    assert occupied < len(stand.trees), "coarse cells should hide some stems"
    # At a resolution that can actually resolve crowns, nothing is lost.
    fine = synthesize(
        SynthConfig(
            width=60.0,
            height=60.0,
            cellsize=0.5,
            stems_per_ha=300.0,
            seed=4,
            gap_count=0,
            surface_noise=0.0,
        )
    )
    assert sum(1 for value in fine.chm.values if value and value > 0) > len(fine.trees)
