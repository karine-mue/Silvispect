# Command line reference

```console
$ silvispect --help
usage: silvispect [-h] [--version]
                  {chm,detect,metrics,gaps,render,match,allometry,inspect,synth,demo}
```

Run without installing with `python -m silvispect ...`.

**Exit codes** are uniform across subcommands:

| Code | Meaning |
| --- | --- |
| `0` | Success |
| `1` | `inspect` found something at or above `--fail-on` |
| `2` | Bad input: unreadable file, malformed raster or CSV, invalid config |

---

## `inspect`

Run the rule set and produce a report. The heart of the tool.

```console
$ silvispect inspect [--chm CHM] [--inventory CSV] [--area-ha HA]
                     [--config JSON] [--format text|markdown|json]
                     [-o OUTPUT] [--no-detection] [--with-map]
                     [--fail-on none|notice|warning|critical]
```

| Option | Default | Notes |
| --- | --- | --- |
| `--chm` | — | Canopy height model. Enables the canopy, gap and detection rules |
| `--inventory` | — | Field tree list. Enables the data-quality rules |
| `--area-ha` | Raster footprint | Required when no raster is given |
| `--config` | Built-in defaults | JSON threshold overrides; unknown keys are rejected |
| `--format` | `text` | `markdown` for a full report, `json` for machines |
| `--no-detection` | off | Skip tree detection; canopy and gap rules still run |
| `--with-map` | off | Append an ASCII canopy map |
| `--fail-on` | `none` | Exit `1` when a finding reaches this severity |

At least one of `--chm` and `--inventory` is required.

```console
# Human report to a file
$ silvispect inspect --chm chm.asc --inventory trees.csv \
      --format markdown --with-map -o report.md

# Machine-readable, piped into jq
$ silvispect inspect --chm chm.asc --format json | jq '.findings[].code'

# As a gate in CI
$ silvispect inspect --chm chm.asc --inventory trees.csv --fail-on warning
```

---

## `detect`

Detect treetops and delineate crowns from a canopy height model.

```console
$ silvispect detect CHM [--min-height M] [--smooth-radius CELLS]
                        [--window-intercept M] [--window-slope M_PER_M]
                        [-o TREES.CSV] [--labels LABELS.ASC] [--json]
```

| Option | Default | Notes |
| --- | --- | --- |
| `--min-height` | `2.0` | Cells below this cannot belong to a tree |
| `--smooth-radius` | `1` | Mean-filter radius in cells; `0` keeps raw apex heights |
| `--window-intercept` | `0.8` | Constant term of the search-window radius (m) |
| `--window-slope` | `0.055` | Growth of that radius per metre of height |
| `-o` | — | Write detections as an inventory CSV |
| `--labels` | — | Write a raster of crown IDs |

Raise the window parameters if crowns are being split; lower them if
neighbouring trees are being merged.

---

## `gaps`

Map canopy openings.

```console
$ silvispect gaps CHM [--threshold M] [--min-area M2] [--min-width M] [--json]
```

| Option | Default | Notes |
| --- | --- | --- |
| `--threshold` | `2.0` | Height at or above which a cell counts as canopy |
| `--min-area` | `25.0` | Smallest reported opening (m²) |
| `--min-width` | `5.0` | Brokaw criterion; `0` reports the raw threshold components |

`--min-width 0` is instructive rather than useful: it shows the connected
inter-crown web that the morphological opening exists to remove. See
[`concepts.md`](concepts.md#2-canopy-cover-gaps-and-strata).

---

## `match`

Compare detections against a field inventory.

```console
$ silvispect match INVENTORY.CSV --chm CHM [--tolerance M] [--json]
```

Reports matched / omitted / spurious counts, recall, precision, F1, height
bias, height RMSE and mean positional offset. `--tolerance` (default 2.5 m) is
the matching radius — set it to roughly the mean crown radius of the stand.

---

## `metrics`

Summarise an inventory.

```console
$ silvispect metrics INVENTORY.CSV --area-ha HA [--include-dead] [--json]
```

`--area-ha` is required: Silvispect will not guess an expansion factor. Dead
stems are excluded unless `--include-dead` is given.

---

## `allometry`

Fit a height-diameter curve to a whole inventory.

```console
$ silvispect allometry INVENTORY.CSV [--model naslund|curtis|power|chapman] [--json]
```

Reports the fitted parameters, RMSE and R². This fits one pooled curve; the
`inspect` command additionally fits per species.

---

## `chm`

Derive a canopy height model.

```console
$ silvispect chm --dsm DSM.ASC --dtm DTM.ASC [-o CHM.ASC]
```

The two rasters must share shape, cell size and origin. Negative heights are
clamped to zero; nodata in either input stays nodata.

---

## `render`

Draw a raster as an ASCII map.

```console
$ silvispect render CHM [--width CHARS] [--tops]
```

Downsampling is by block **maximum**, so treetops survive the reduction rather
than being averaged away. `--tops` overlays detected apexes as `^`.

---

## `synth`

Generate a synthetic stand with known ground truth.

```console
$ silvispect synth [--out DIR] [--seed N] [--width M] [--height M]
                   [--cellsize M] [--stems-per-ha N] [--gaps N]
```

Writes `chm.asc`, `dsm.asc`, `dtm.asc` and `trees.csv`. The same seed always
produces the same stand, byte for byte.

---

## `demo`

Generate a stand and inspect it end to end — the fastest way to see what the
tool does.

```console
$ silvispect demo [--seed N] [--format text|markdown|json]
```
