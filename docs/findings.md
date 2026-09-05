# Findings reference

Every rule in Silvispect emits findings with a stable code, a severity, a
subject and a numeric value against the threshold that produced it. Codes never
change meaning; new rules take new codes.

## Severities

| Severity | Meaning | Typical response |
| --- | --- | --- |
| `info` | Contextual, not a problem | None |
| `notice` | Worth knowing; may be entirely normal for the site | Read it |
| `warning` | Likely a real problem with the data or the stand | Investigate |
| `critical` | The data cannot be right, or the stand is in a state that demands action | Fix before use |

`silvispect inspect --fail-on <severity>` exits `1` when any finding reaches
that level, which is what makes the tool usable as a gate.

## Which rules run

Rules skip silently when their inputs are absent, so partial data is not an
error:

| Input supplied | Rules that run |
| --- | --- |
| Canopy model only | Stocking (from detections), canopy, gaps, structure |
| Canopy model with `--no-detection` | Canopy and gaps only |
| Inventory only | Stocking, structure, composition, all data-quality rules |
| Both | Everything, including detection agreement |

With a canopy model alone, stem metrics come from the detected crowns with
diameters back-calculated from a default height-diameter curve; the report
labels this `metrics_source: detected` so the numbers are never mistaken for
measurements.

`--no-detection` leaves no stem information at all. The report then says
`metrics_source: none` and the stem-based rules (SV001–SV004, SV020) stay
silent: a stand nobody counted has *unknown* stocking, not zero stocking. A
detector that did run and found nothing is a different matter — that is
evidence of an empty stand, and SV001 fires.

Coverage is judged **cell by cell**. A field stem standing where the raster
holds no measurement is not a stem the detector missed — there was nothing
there to miss it with — so it is left out of the detection-agreement rules
(SV050–SV052) rather than counted as an omission. Stems outside the extent are
excluded for the same reason and reported separately by SV043. Everything else
about the stem, including every data-quality and stocking rule, is unaffected.

A canopy model whose cells are *all* nodata is the first case, not the second.
There is nothing to look at, so detection, gap mapping and the canopy summary
do not run at all: `metrics_source` is `none`, `tree_count` is `null`, the
`detection`, `match` and `gaps` blocks are absent, `canopy` is empty, and the
rules that read them — SV010, SV012, SV011 and the detection-agreement rules
SV050–SV052 — stay silent. In particular no inventory is matched against an
empty detection, which used to report every field stem as missed by a sensor
that had never seen anything and could fail a `--fail-on warning` run. One
valid cell is enough to make the raster an observation.

---

## SV00x — stocking and density

### SV001 Understocked stand — `warning`
Stem density is below `min_stems_per_ha` (default 200). Either the stand is
genuinely poorly stocked, the plot area passed to Silvispect is too large, or —
when the source is detection — the canopy is hiding stems from the sensor.

### SV002 Overstocked stand — `notice`
Density above `max_stems_per_ha` (default 1400). Competition-driven mortality
is likely; this is a management signal, not a data error, hence `notice`.

### SV003 Low stand density index — `notice`
Reineke SDI below `min_sdi` (default 150). The site is not fully occupying the
available growing space. SDI is a better cross-stand comparison than raw stem
count because it accounts for tree size.

### SV004 High stand density index — `warning`
SDI above `max_sdi` (default 900). The stand is approaching the zone of
imminent competition mortality.

---

## SV01x — canopy condition

### SV010 Canopy cover below target — `warning`
Cover at the `canopy_threshold` (default 2 m) is below `min_canopy_cover`
(default 0.60). Check the terrain model before concluding the canopy is open:
an inflated DTM lowers every CHM cell at once.

### SV011 Large canopy gap — `critical`, or `warning` at the plot edge
A single opening exceeds `max_gap_area` (default 400 m²). Gaps are found by
morphological opening with a `min_gap_width` disc (default 5 m), so slivers of
inter-crown ground are already excluded — a finding here is a real opening.

Raising `min_gap_width` can only ever shrink the openings reported, never grow
them, so tightening the criterion cannot conjure this finding into existence.
Gaps whose *untrimmed* sub-canopy area reaches the plot boundary are truncated
by the extent, so their area is a lower bound and their severity is reduced to
`warning`.

`min_gap_area` is the reporting floor and `max_gap_area` the limit this rule
enforces, so a profile setting the floor above the limit is refused rather than
applied: it would delete exactly the openings the rule exists to report, and a
plot with a 100 m² gap against a 50 m² limit came back clean.

### SV012 High open-canopy fraction — `notice`
More than `max_gap_fraction` (default 0.30) of the plot is below the canopy
threshold. Unlike SV011 this counts *all* sub-canopy cells, including the
inter-crown web, so an open-grown stand can trip it without having any single
large opening.

---

## SV02x — structure and composition

### SV020 Low structural diversity — `notice`
The Gini coefficient of basal area is below `min_gini` (default 0.15). A
structurally uniform stand is vulnerable to a single disturbance agent, since
every tree is at the same stage and the same size.

### SV021 Low species diversity — `notice`
Shannon index below `min_shannon` (default 0.35). Requires field data — a
canopy model carries no species information.

### SV022 Single-species dominance — `notice`
One species accounts for more than `max_species_share` of stems (default 0.90).

---

## SV03x — measurement plausibility

### SV030 Height-diameter outlier — `warning` (`notice` when marginal)
A stem departs from its species' fitted height-diameter curve by more than
`height_residual_z` standard deviations of the residuals (default 2.5).

Common causes, in rough order of frequency: a transcription error in one of the
two dimensions; a broken or dead top, which makes the tree genuinely short for
its diameter; a forked stem measured below the fork; a species miscoded.

Curves are fitted **per species** when a species has at least
`min_allometry_points` paired measurements, precisely so this rule does not
simply rediscover the species mixture. The finding names the curve used.

### SV031 Diameter outlier — `notice`
A diameter is more than `dbh_outlier_z` standard deviations from the plot mean
(default 3.0). Uninformative in a strongly uneven-aged stand, where large
outliers are the structure rather than an error — read it alongside SV020.

---

## SV04x — record integrity

These are the rules that catch broken data rather than an unusual forest.

### SV040 Missing measurement — `notice`
A live stem has no diameter, no height, or neither. Dead stems are exempt:
snags are routinely recorded without a full set of measurements.

### SV041 Duplicate tree identifier — `critical`
The same `tree_id` appears more than once. Records cannot be linked across
re-measurements, so any growth analysis built on this dataset is unsound.

### SV042 Implausible dimension — `critical`
A non-positive diameter or height, or one above `max_plausible_dbh_cm` (250) /
`max_plausible_height_m` (120). A zero or negative dimension is never a
measurement; it is a placeholder that escaped.

### SV043 Coordinate outside raster extent — `warning`
A field stem falls outside the canopy model envelope. Usually a coordinate
system mismatch, a transposed easting and northing, or a plot origin offset —
all of which invalidate the detection comparison rather than just one row.

---

## SV05x — remote sensing agreement

These require both a canopy model and a field inventory.

### SV050 Low detection recall — `warning`
Fewer than `min_recall` (default 0.70) of field stems were matched to a
detection. Before tuning the detector, ask whether the missing stems are
suppressed: a canopy height model cannot see a tree that never reaches the
canopy surface, and no parameter change will make it.

### SV051 Low detection precision — `warning`
Fewer than `min_precision` (default 0.70) of detections matched a field stem.
Usually crown splitting: broad or irregular crowns produce several local maxima.
Increase `window_intercept` / `window_slope`, or smooth more aggressively.

### SV052 Systematic height bias — `warning`
Detected heights differ from field heights by more than `max_height_bias_m`
(default 1.5 m) on average. A *systematic* bias points at the terrain model,
not at the detector — an error in the DTM shifts every CHM cell by the same
amount. Note that pre-smoothing produces a small, expected negative bias of a
few decimetres.

---

## Configuring thresholds

Every threshold above is a field of `InspectionConfig`. Override any subset
from JSON:

```json
{
  "min_canopy_cover": 0.7,
  "max_gap_area": 250.0,
  "height_residual_z": 2.0
}
```

```console
$ silvispect inspect --chm chm.asc --inventory trees.csv --config profile.json
```

A profile is validated when it loads, and anything unusable is rejected with
exit status `2` rather than surfacing later as a traceback:

- unknown keys, so a typo fails loudly instead of leaving a default in place;
- non-numeric or non-object payloads, and non-finite numbers;
- fractional counts, rather than truncating `1.9` to `1` without a word;
- negative thresholds, and proportions outside `[0, 1]`;
- `min_allometry_points` below 1;
- contradictory pairs such as `min_stems_per_ha` above `max_stems_per_ha`,
  which no stand could ever satisfy.

A worked example ships as
[`data/thresholds-strict.json`](../data/thresholds-strict.json).

The defaults are plausible for a temperate managed forest and wrong for
anything else. They are a starting point for your own silvicultural
expectations, not a standard.

## Adding a rule

Rules are plain functions registered with a decorator:

```python
from silvispect.inspection import Finding, InspectionContext, Severity, rule


@rule
def check_veteran_trees(ctx: InspectionContext):
    """SV060 — no large-diameter habitat trees retained."""
    if not ctx.has_field_data:
        return
    veterans = [t for t in ctx.field_trees if (t.dbh_cm or 0) >= 70]
    if not veterans:
        yield Finding(
            "SV060",
            Severity.NOTICE,
            "No veteran trees",
            "No stem reaches 70 cm DBH; habitat continuity is not provided.",
        )
```

A rule receives the whole context — metrics, crowns, gaps, the fitted models,
the raster — and yields zero or more findings. It must not raise on missing
inputs; check `ctx.has_field_data` / `ctx.has_raster` and return quietly.
