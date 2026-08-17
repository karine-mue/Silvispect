"""Shared fixtures."""

from __future__ import annotations

import pytest

from silvispect.grid import Grid
from silvispect.synth import SynthConfig, synthesize


@pytest.fixture(scope="session")
def stand():
    """A small deterministic synthetic stand shared by the whole session."""
    return synthesize(SynthConfig(width=60.0, height=60.0, stems_per_ha=280.0, seed=7, gap_count=1))


@pytest.fixture(scope="session")
def sparse_stand():
    """A widely spaced stand where detection should be near perfect."""
    return synthesize(
        SynthConfig(
            width=60.0,
            height=60.0,
            stems_per_ha=140.0,
            seed=11,
            gap_count=0,
            surface_noise=0.0,
        )
    )


@pytest.fixture
def tiny_grid() -> Grid:
    """A 4x4 grid with one obvious peak, cellsize 1 m."""
    return Grid.from_rows(
        [
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 5.0, 4.0, 0.0],
            [0.0, 4.0, 3.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ],
        cellsize=1.0,
    )
