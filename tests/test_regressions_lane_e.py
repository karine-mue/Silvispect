"""Regression tests — Lane E (detect.py / synth.py). Issues #25, #27."""

import pytest

from silvispect.detect import (
    DetectionConfig,
    TreeTop,
    find_treetops,
    segment_crowns,
)
from silvispect.grid import Grid, GridError
from silvispect.synth import SynthConfig


# --- #25: config validation ---------------------------------------------------

def test_detection_config_rejects_negative_window_slope():
    with pytest.raises(GridError):
        DetectionConfig(window_slope=-0.1)


def test_detection_config_rejects_negative_crown_slope():
    with pytest.raises(GridError):
        DetectionConfig(crown_slope=-0.5)


def test_detection_config_defaults_still_valid():
    DetectionConfig()


def test_synth_config_rejects_inverted_dbh_range():
    with pytest.raises(ValueError):
        SynthConfig(dbh_min=90.0, dbh_max=10.0)


def test_synth_config_rejects_all_zero_weights():
    with pytest.raises(ValueError):
        SynthConfig(species_weights=(0.0, 0.0, 0.0))


def test_synth_config_rejects_negative_weight():
    with pytest.raises(ValueError):
        SynthConfig(species_weights=(0.6, -0.1, 0.5))


def test_synth_config_rejects_negative_surface_noise():
    with pytest.raises(ValueError):
        SynthConfig(surface_noise=-0.1)


def test_synth_config_rejects_float_gap_count():
    with pytest.raises(ValueError):
        SynthConfig(gap_count=1.0)


def test_synth_config_rejects_inverted_gap_radius():
    with pytest.raises(ValueError):
        SynthConfig(gap_radius=(9.0, 5.0))


def test_synth_config_defaults_still_valid():
    SynthConfig()


# --- #27: metric search window ---------------------------------------------------

def test_treetop_window_uses_metric_distance():
    chm = Grid.from_rows([
        [10.5, 5.0, 5.0],
        [5.0, 10.0, 5.0],
        [5.0, 5.0, 5.0],
    ], cellsize=1.0)
    config = DetectionConfig(smooth_radius=0, window_intercept=1.49,
                             window_slope=0.0)
    tops = find_treetops(chm, config)
    assert {(t.row, t.col) for t in tops} == {(0, 0)}


# --- #27: crowns do not climb ------------------------------------------------------

def test_crown_does_not_climb_micro_rises():
    step = 5e-9
    chm = Grid.from_rows(
        [[10.0, 10.0 + step, 10.0 + 2 * step, 10.0 + 3 * step, 10.0 + 4 * step]],
        cellsize=1.0,
    )
    top = TreeTop(tree_id=1, row=0, col=0, x=0.0, y=0.0, height=10.0)
    config = DetectionConfig(min_crown_cells=1, crown_intercept=100.0,
                             crown_slope=0.0, drop_fraction=0.0)
    _, labels = segment_crowns(chm, [top], config)
    assert labels.values[0] == 1
    assert all(v is None for v in labels.values[1:])
