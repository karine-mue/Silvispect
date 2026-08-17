#!/usr/bin/env python3
"""Regenerate the sample plot committed in ``data/``.

The sample is a 0.25 ha synthetic stand of Norway spruce, beech and birch with
one cut gap, georeferenced somewhere plausible in a projected metric system.
Because the generator is seeded, running this script reproduces the committed
files byte for byte — it exists so the data can be regenerated deliberately,
not so it drifts.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:  # allow running from a checkout without installing
    sys.path.insert(0, str(ROOT))

from silvispect.inspection import InspectionConfig  # noqa: E402
from silvispect.inventory import write_trees  # noqa: E402
from silvispect.synth import SynthConfig, synthesize  # noqa: E402

DATA = ROOT / "data"

PLOT = SynthConfig(
    width=50.0,
    height=50.0,
    cellsize=0.5,
    stems_per_ha=340.0,
    seed=1848,
    gap_count=1,
    gap_radius=(6.0, 8.0),
    xllcorner=612_500.0,
    yllcorner=6_642_000.0,
)

STRICT = InspectionConfig(
    min_canopy_cover=0.70,
    max_gap_area=250.0,
    min_gini=0.25,
    height_residual_z=2.0,
    min_recall=0.80,
    min_precision=0.80,
)


def main() -> int:
    DATA.mkdir(parents=True, exist_ok=True)
    stand = synthesize(PLOT)
    stand.chm.write(DATA / "plot-a1-chm.asc", precision=2)
    stand.dsm.write(DATA / "plot-a1-dsm.asc", precision=2)
    stand.dtm.write(DATA / "plot-a1-dtm.asc", precision=2)
    write_trees(stand.trees, DATA / "plot-a1-trees.csv", precision=2)

    profile = DATA / "thresholds-strict.json"
    profile.write_text(
        json.dumps(STRICT.as_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    print(
        f"wrote {len(stand.trees)} trees over {stand.area_ha:.2f} ha "
        f"({len(stand.trees) / stand.area_ha:.0f} stems/ha) to {DATA}/"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
