"""Silvispect — inspect forest canopy models and inventories.

The package is organised as a short pipeline:

``grid`` -> ``canopy`` -> ``detect`` -> ``metrics`` -> ``inspection`` -> ``report``

with :mod:`silvispect.inventory` and :mod:`silvispect.allometry` supplying the
field-data side and :mod:`silvispect.synth` generating reproducible test
stands.  Everything runs on the Python standard library.

Typical use::

    from silvispect import Grid, inspect_stand, render_text

    chm = Grid.read("chm.asc")
    report = inspect_stand(chm=chm)
    print(render_text(report))
"""

from __future__ import annotations

__version__ = "0.1.1"

from .allometry import FittedModel, fit_height_diameter
from .canopy import canopy_cover, canopy_height_model, find_gaps
from .detect import DetectionConfig, detect_trees
from .grid import Grid, GridError
from .inspection import (
    Finding,
    InspectionConfig,
    InspectionReport,
    Severity,
    inspect_stand,
)
from .inventory import Tree, match_trees, read_trees, write_trees
from .metrics import StandMetrics, stand_metrics
from .report import ascii_map, render_json, render_markdown, render_text
from .synth import SynthConfig, synthesize

__all__ = [
    "DetectionConfig",
    "Finding",
    "FittedModel",
    "Grid",
    "GridError",
    "InspectionConfig",
    "InspectionReport",
    "Severity",
    "StandMetrics",
    "SynthConfig",
    "Tree",
    "__version__",
    "ascii_map",
    "canopy_cover",
    "canopy_height_model",
    "detect_trees",
    "find_gaps",
    "fit_height_diameter",
    "inspect_stand",
    "match_trees",
    "read_trees",
    "render_json",
    "render_markdown",
    "render_text",
    "stand_metrics",
    "synthesize",
    "write_trees",
]
