"""The sample plot in ``data/`` is part of the documented interface.

The README, the docs and the CI smoke test all quote figures from it, so these
tests guard against it drifting or becoming unreadable.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from silvispect.canopy import canopy_height_model
from silvispect.grid import Grid
from silvispect.inspection import InspectionConfig, inspect_stand
from silvispect.inventory import read_trees

DATA = Path(__file__).resolve().parent.parent / "data"
PLOT_AREA_HA = 0.25


@pytest.fixture(scope="module")
def sample_chm() -> Grid:
    return Grid.read(DATA / "plot-a1-chm.asc")


@pytest.fixture(scope="module")
def sample_trees():
    return read_trees(DATA / "plot-a1-trees.csv")


def test_sample_files_are_present():
    for name in (
        "plot-a1-chm.asc",
        "plot-a1-dsm.asc",
        "plot-a1-dtm.asc",
        "plot-a1-trees.csv",
        "thresholds-strict.json",
    ):
        assert (DATA / name).is_file(), f"missing sample file {name}"


def test_sample_raster_geometry(sample_chm):
    assert (sample_chm.nrows, sample_chm.ncols) == (100, 100)
    assert sample_chm.cellsize == 0.5
    assert sample_chm.area_ha == pytest.approx(PLOT_AREA_HA)
    assert sample_chm.xllcorner == pytest.approx(612_500.0)
    assert sample_chm.stats().nodata_count == 0


def test_sample_chm_matches_the_shipped_surfaces(sample_chm):
    derived = canopy_height_model(
        Grid.read(DATA / "plot-a1-dsm.asc"), Grid.read(DATA / "plot-a1-dtm.asc")
    )
    for expected, actual in zip(sample_chm.values, derived.values, strict=True):
        assert expected == pytest.approx(actual, abs=0.02)


def test_sample_inventory(sample_trees):
    assert len(sample_trees) == 85
    assert all(tree.is_live for tree in sample_trees)
    assert all(tree.dbh_cm and tree.height_m for tree in sample_trees)
    assert {tree.species for tree in sample_trees} == {"PIAB", "FASY", "BEPE"}


def test_sample_stems_fall_inside_the_raster(sample_chm, sample_trees):
    xmin, ymin, xmax, ymax = sample_chm.extent
    for tree in sample_trees:
        assert xmin <= tree.x <= xmax
        assert ymin <= tree.y <= ymax


def test_sample_inspects_cleanly(sample_chm, sample_trees):
    """The documented sample must stay a well-behaved example."""
    report = inspect_stand(chm=sample_chm, trees=sample_trees, area_ha=PLOT_AREA_HA)
    assert report.metrics.tree_count == 85
    assert report.metrics.stems_per_ha == pytest.approx(340.0)
    assert report.match["recall"] > 0.95
    assert report.match["precision"] > 0.95
    assert report.max_severity.label == "notice"


def test_strict_profile_loads_and_bites(sample_chm, sample_trees):
    config = InspectionConfig.load(DATA / "thresholds-strict.json")
    assert config.min_canopy_cover == 0.70
    strict = inspect_stand(chm=sample_chm, trees=sample_trees, area_ha=PLOT_AREA_HA, config=config)
    default = inspect_stand(chm=sample_chm, trees=sample_trees, area_ha=PLOT_AREA_HA)
    assert len(strict.findings) > len(default.findings)
