# Silvispect

**Inspect a forest the way you inspect a codebase.**

Silvispect reads a canopy height model and/or a field inventory, detects
individual trees, computes stand mensuration metrics, and runs a fixed set of
rules that emit **coded findings** with severities — so a stand can be checked
in CI and fail a build when it breaks a threshold.

The name is *silva* (forest) + *spectare* (to inspect), and the tool takes both
halves literally.

- **No dependencies.** Pure Python standard library, Python 3.10+.
- **Plain text in, plain text out.** ESRI ASCII Grid rasters and CSV
  inventories, Markdown / JSON / terminal reports.
- **Reproducible.** A seeded synthetic stand generator provides ground truth,
  so detection accuracy is measured, not asserted.

```console
$ silvispect inspect --chm data/plot-a1-chm.asc --inventory data/plot-a1-trees.csv
Silvispect inspection
=====================
area              0.250 ha (field data)
trees             85
stems/ha          340
yield basis       85 stems
basal area/ha     25.11 m2
qmd               30.67 cm
lorey height      23.05 m
dominant height   24.38 m
volume/ha         289.38 m3
canopy cover      63.6%
gaps mapped       3
detection         recall 98.8%, precision 100.0%

findings (2 notice):
  [-] SV012 stand         High open-canopy fraction: 36.4% of the plot is below 2.0 m, above the 30% limit.
  [-] SV030 tree:T0072    Height-diameter outlier: Tree T0072 is 1.9 m taller than the FASY naslund curve predicts (z = +2.6); check for crown breakage or a transcription error.
```

## Why

Forest inventory data goes wrong in specific, repetitive ways: a transposed
digit in a diameter, a stem measured twice under two identifiers, a plot whose
GPS drifted, a canopy model built on the wrong terrain surface, a detector
tuned until it splits every crown in two. These are all *checkable*. Silvispect
turns that checking into a fixed rule set with stable codes, so the same
question — "is this stand plausible?" — gets the same answer every time, in a
form a script can act on.

## Install

```console
$ git clone https://github.com/karine-mue/Silvispect.git
$ cd Silvispect
$ pip install -e .
```

There is nothing to compile and nothing to download. You can also run it
straight from a checkout without installing:

```console
$ python -m silvispect demo
```

## Quick tour

Every subcommand works on the sample plot in [`data/`](data/), a 0.25 ha
synthetic stand of Norway spruce, beech and birch with one cut gap.

```console
# See it
$ silvispect render data/plot-a1-chm.asc --width 50 --tops

# Detect trees and write them as an inventory
$ silvispect detect data/plot-a1-chm.asc -o /tmp/detected.csv
detected 84 trees (336 stems/ha over 0.250 ha)
height  min 13.7 m  mean 21.5 m  max 27.5 m
crown   min 2.8 m  mean 4.1 m  max 5.5 m

# Map the canopy openings
$ silvispect gaps data/plot-a1-chm.asc
3 gaps at or above 25 m2 and 5.0 m wide; 36.4% of the plot is below 2.0 m
  gap 1        170 m2  width  12.3 m  centre (612525.4, 6642009.3)  [edge]
  gap 2         68 m2  width   5.8 m  centre (612545.6, 6642035.5)  [edge]
  gap 3         51 m2  width   5.8 m  centre (612508.6, 6642018.5)  [edge]

# Compare the detections against the field record
$ silvispect match data/plot-a1-trees.csv --chm data/plot-a1-chm.asc
matched          84
omissions        1
commissions      0
recall           0.988
precision        1
height_bias_m    -0.688

# Fit the stand's own height-diameter curve
$ silvispect allometry data/plot-a1-trees.csv

# Summarise the inventory
$ silvispect metrics data/plot-a1-trees.csv --area-ha 0.25

# Derive a canopy height model from a surface and a terrain model
$ silvispect chm --dsm data/plot-a1-dsm.asc --dtm data/plot-a1-dtm.asc -o /tmp/chm.asc

# Generate a fresh stand with known ground truth
$ silvispect synth --out /tmp/stand --seed 7 --stems-per-ha 450

# The whole pipeline on a generated stand
$ silvispect demo --format markdown
```

```
-*^0*-:==-.+0^o=*o*-*00o+.=o^o+=:-*^o+-++=:=**+:-:
.=++++*o^*o0oo*+o^*-o0^0*-.:+0^0+:*0o+*0^o+o0^0*^o
  +0^#ooo*o^0*:.--:.+*oo00o-=o0o=  . .+o0o+*00o*oo
  +0##o-.:+++-:.*00o- +0^0o-*^*-: -++-  .=^*= -o^#
*00*--:   -*0#0**0^o--*0^*o0o*0#0*-+^=   .--. :*00
o^0*:   :--o#^0*:--. -*ooo#^0o#^#o++- .=+=:
-++- .::+^*=*o*=*oo+: .::+oo*+*o*=++- =o^0*:    =+
  .--*^o=++:  .+o^0*=**+-*oo+^*=   .---o0o+.   -*^
+=o##0o*o00*:--=*oo+*0^o=0^0*--.   =*^+:--.    :+*
^o0#^0+*0^#o0##0*+0^#*o+:+*o*: -==.-+*+*^0*.
o*=*o^o=ooo+0@^#*=0##*:::=o^*::+^+-   :*0o+.  -*o*
  low 0.0 m [ .:-=+*o0#@] 27.9 m high
```

## Use it as a check

`silvispect inspect` exits non-zero when a finding reaches the severity you
name, which makes it a gate rather than a report:

```yaml
- name: Inspect inventory
  run: |
    silvispect inspect --chm survey/chm.asc \
                       --inventory survey/trees.csv \
                       --config profiles/beech-target.json \
                       --fail-on warning
```

Exit codes: `0` clean, `1` findings at or above `--fail-on`, `2` bad input.

## Findings

Rules are grouped by the number range of their code. Full descriptions and the
reasoning behind each threshold are in [`docs/findings.md`](docs/findings.md).

| Code | Severity | Finding |
| --- | --- | --- |
| SV001 / SV002 | warning / notice | Understocked / overstocked stem density |
| SV003 / SV004 | notice / warning | Stand density index below / above the managed range |
| SV010 | warning | Canopy cover below target |
| SV011 | critical | Canopy gap larger than allowed |
| SV012 | notice | Too much of the plot is open canopy |
| SV020 | notice | Low structural diversity (Gini of basal area) |
| SV021 / SV022 | notice | Low species diversity / single-species dominance |
| SV030 | warning | Height-diameter outlier against the stand's own curve |
| SV031 | notice | Diameter outlier against the plot distribution |
| SV040 | notice | Live stem with a missing measurement |
| SV041 | critical | Duplicate tree identifier |
| SV042 | critical | Non-positive or implausible dimension |
| SV043 | warning | Stem coordinate outside the raster extent |
| SV050 / SV051 | warning | Detection recall / precision below target |
| SV052 | warning | Systematic height bias between detection and field |

Every threshold is a documented field of `InspectionConfig` and can be
overridden from a JSON file — see [`data/thresholds-strict.json`](data/thresholds-strict.json).

## Use it as a library

```python
from silvispect import Grid, detect_trees, inspect_stand, render_markdown
from silvispect.inventory import read_trees

chm = Grid.read("data/plot-a1-chm.asc")
trees = read_trees("data/plot-a1-trees.csv")

detection = detect_trees(chm)
print(f"{len(detection)} crowns, {detection.density_per_ha():.0f} stems/ha")

report = inspect_stand(chm=chm, trees=trees, area_ha=0.25)
for finding in report.findings:
    print(finding.severity.label, finding.code, finding.detail)

open("report.md", "w").write(render_markdown(report, chm=chm))
```

The package is a short pipeline of focused modules:

| Module | Responsibility |
| --- | --- |
| `grid` | ESRI ASCII Grid raster: I/O, geometry, focal operations |
| `canopy` | Canopy height model, cover, morphological gap mapping, strata |
| `detect` | Variable-window treetop detection, marker-controlled crown growing |
| `inventory` | Tree records, CSV interchange, detection/field matching |
| `allometry` | Height-diameter model families, fitting, residual screening |
| `optimize` | A self-contained Nelder-Mead simplex search |
| `metrics` | Basal area, QMD, Lorey height, SDI, Gini, volume, biomass |
| `inspection` | The rule registry and the report it produces |
| `report` | ASCII canopy maps, Markdown / JSON / text renderers |
| `synth` | Deterministic synthetic stands with ground truth |
| `cli` | The `silvispect` command |

## How it works

**Tree detection** uses a search window whose radius grows with the height of
the cell being tested, because crown width scales with tree height — a fixed
window either merges small neighbours or splits large crowns. Each maximum
seeds a crown, and cells are absorbed tallest-first while the surface keeps
descending away from the apex, placing the boundary in the valley between
adjacent crowns.

**Gap mapping** does not simply threshold the canopy. On a fine raster the
ground visible between crowns forms one connected web that a naive
connected-components pass reports as a single enormous opening. Silvispect
morphologically opens the sub-canopy mask with a disc of the minimum gap width
(the Brokaw criterion) before labelling, which trims that web — including the
parts attached to a genuine opening — and leaves the real gaps intact.

**Height-diameter curves** are fitted per species, because a single pooled
curve on a mixed stand systematically flags the minority species as outliers.
Näslund's equation is the default: it passes near the origin and saturates at
`1.3 + 1/b²`, so a plausible asymptotic height is readable straight off the
parameters. Curtis, a power law and Chapman-Richards are also available.

More detail — including the mensuration formulae and their assumptions — is in
[`docs/concepts.md`](docs/concepts.md).

## Data formats

Rasters are [ESRI ASCII Grid](docs/data-formats.md) (`.asc`). Inventories are
CSV with required `x`, `y` columns and optional `tree_id`, `species`, `dbh_cm`,
`height_m`, `crown_diameter_m`, `status`; common aliases (`dbh`, `height`,
`easting`, `northing`, …) are accepted case-insensitively. Full specification
in [`docs/data-formats.md`](docs/data-formats.md).

## Development

```console
$ pip install -e ".[dev]"
$ pytest              # 250 tests
$ ruff check .
$ ruff format --check .
```

Detection accuracy is a *test*, not a claim: the suite generates a synthetic
stand, runs the detector against the known stems, and asserts recall,
precision, height RMSE and positional offset. Fitting is checked by recovering
the parameters the generator used.

## Limitations

Read these before pointing it at real data.

- **Canopy models cannot see suppressed trees.** Recall against a full field
  inventory is bounded by what reaches the canopy surface; SV050 is often a
  true statement about the sensor, not about the detector.
- **The thresholds are placeholders for your silviculture.** The defaults are
  plausible for a temperate managed forest and wrong for anything else. Set
  them deliberately.
- **Species defaults are illustrative.** The tabulated Näslund parameters are
  representative, not regional yield-table values. Fit your own where you can;
  Silvispect prefers a fitted curve automatically whenever the data supports it.
- **Volume uses a form factor**, not a taper function, and biomass uses a
  pantropical equation. Both are order-of-magnitude tools.
- **Everything assumes metric map units** and a projected coordinate system.

## Licence

MIT — see [LICENSE](LICENSE).
