"""Regression tests for the Gemini 4 Neon evaluation findings (#31-#34).

Each test pins one execution-verified finding from the 2026-11-19 bench of
the fixed code: a round trip that corrupted absence, a label raster the
writer refused, and validation that stopped one parameter short of the
non-finite values it was added to refuse.
"""

from __future__ import annotations

import io

import pytest

from silvispect.allometry import default_model
from silvispect.canopy import find_gaps
from silvispect.cli import EXIT_OK, main
from silvispect.detect import TreeTop, segment_crowns
from silvispect.grid import Grid, GridError
from silvispect.inventory import Tree
from silvispect.metrics import MetricsError, above_ground_biomass
from silvispect.report import ascii_map
from silvispect.synth import SynthConfig

# --- #31: to_text wrote absence at cell precision -------------------------


def test_nodata_cell_round_trips_at_any_precision():
    grid = Grid(
        2, 1, 0.0, 0.0, 1.0, nodata_value=-9999.1234, values=[5.0, None]
    )
    for precision in (0, 2, 3, None):
        back = Grid.parse(grid.to_text(precision=precision))
        assert back.values[0] == 5.0
        assert back.values[1] is None


def test_valid_cell_matching_nodata_after_formatting_still_refused():
    # The #11 guard is untouched: a measurement that would be written as
    # the sentinel is refused, it does not silently become absence.
    grid = Grid(1, 1, 0.0, 0.0, 1.0, nodata_value=-9999.0, values=[-9999.04])
    with pytest.raises(GridError):
        grid.to_text(precision=1)


# --- #32: label grid inherited the CHM sentinel ---------------------------


def _collision_chm() -> Grid:
    return Grid(
        5,
        5,
        0.0,
        0.0,
        1.0,
        nodata_value=1.0,
        values=[
            5.0, 5.0, 5.0, 5.0, 5.0,
            5.0, 12.0, 15.0, 12.0, 5.0,
            5.0, 15.0, 20.0, 15.0, 5.0,
            5.0, 12.0, 15.0, 12.0, 5.0,
            5.0, 5.0, 5.0, 5.0, 5.0,
        ],
    )


def test_label_grid_writes_when_chm_nodata_is_a_tree_id():
    chm = _collision_chm()
    top = TreeTop(tree_id=1, row=2, col=2, x=2.0, y=2.0, height=20.0)
    crowns, labels = segment_crowns(chm, [top])
    assert len(crowns) == 1
    back = Grid.parse(labels.to_text())
    assert back.values == labels.values
    assert back.values[2 * 5 + 2] == 1.0


def test_cli_detect_labels_survives_chm_nodata_collision(tmp_path):
    chm_path = tmp_path / "chm.asc"
    labels_path = tmp_path / "labels.asc"
    _collision_chm().write(chm_path)
    buffer = io.StringIO()
    code = main(["detect", str(chm_path), "--labels", str(labels_path)], out=buffer)
    assert code == EXIT_OK
    written = Grid.read(labels_path)
    assert written.stats().maximum == 1.0


# --- #33: robustness of the numeric entry points --------------------------


def test_invert_of_non_finite_height_is_none():
    model = default_model("")
    assert model.invert(float("nan")) is None
    assert model.invert(float("inf")) is None
    assert model.invert(float("-inf")) is None


def test_biomass_rejects_impossible_wood_density():
    tree = Tree(tree_id="1", x=0.0, y=0.0, dbh_cm=30.0, height_m=20.0, species="SUGI")
    for density in (-0.5, 0.0, float("nan"), float("inf")):
        with pytest.raises(MetricsError):
            above_ground_biomass(tree, wood_density=density)


def test_ascii_map_rejects_non_finite_bounds():
    grid = Grid(2, 1, 0.0, 0.0, 1.0, values=[1.0, 2.0])
    with pytest.raises(ValueError, match="finite"):
        ascii_map(grid, minimum=float("nan"))
    with pytest.raises(ValueError, match="finite"):
        ascii_map(grid, maximum=float("inf"))


# --- #34: finiteness validation completeness ------------------------------


def _gap_chm() -> Grid:
    values = [20.0] * 100
    for row in range(4, 7):
        for col in range(4, 7):
            values[row * 10 + col] = 0.0
    return Grid(10, 10, 0.0, 0.0, 1.0, values=values)


def test_find_gaps_rejects_non_finite_filters():
    chm = _gap_chm()
    with pytest.raises(GridError):
        find_gaps(chm, min_area=float("nan"))
    with pytest.raises(GridError):
        find_gaps(chm, min_area=0.0, min_width=float("nan"))
    with pytest.raises(GridError):
        find_gaps(chm, min_area=float("inf"))


def test_histogram_rejects_near_reversed_range_and_non_finite_bounds():
    grid = Grid(2, 1, 0.0, 0.0, 1.0, values=[0.5, 1.5])
    with pytest.raises(GridError):
        grid.histogram(lo=1.0, hi=1.0 - 1e-12)
    with pytest.raises(GridError):
        grid.histogram(lo=float("nan"))
    with pytest.raises(GridError):
        grid.histogram(hi=float("inf"))


def test_histogram_equal_bounds_still_widen():
    grid = Grid(2, 1, 0.0, 0.0, 1.0, values=[2.0, 2.0])
    edges, counts = grid.histogram()
    assert sum(counts) == 2
    edges, counts = grid.histogram(lo=2.0, hi=2.0)
    assert sum(counts) == 2


def test_synth_config_rejects_non_finite_dimensions():
    for overrides in (
        {"width": float("nan")},
        {"height": float("inf")},
        {"cellsize": float("nan")},
        {"stems_per_ha": float("inf")},
        {"dbh_max": float("nan")},
    ):
        with pytest.raises(ValueError):
            SynthConfig(**overrides)
