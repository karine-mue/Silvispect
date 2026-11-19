# Data formats

Silvispect reads and writes rasters as ESRI ASCII Grid (`.asc`) and GeoTIFF
(`.tif`, `.tiff`); inventories and profiles are plain text (CSV, JSON). All
are widely supported interchange formats, so data can move in and out of GIS
and statistical software without a conversion step in the middle.

## Rasters — ESRI ASCII Grid (`.asc`)

```
ncols 100
nrows 100
xllcorner 612500.0
yllcorner 6642000.0
cellsize 0.5
NODATA_value -9999
0 0 1.2 3.4 ...
0 2.1 5.6 7.8 ...
...
```

### Header

| Key | Required | Meaning |
| --- | --- | --- |
| `ncols`, `nrows` | yes | Grid dimensions in cells |
| `cellsize` | yes | Edge length of a cell in map units |
| `xllcorner`, `yllcorner` | no | Lower-left corner of the lower-left cell (default 0) |
| `xllcenter`, `yllcenter` | no | Alternative: centre of that cell; converted on read |
| `NODATA_value` | no | Sentinel for missing cells (default −9999) |

Keys are matched case-insensitively. Cell values follow the header as
whitespace-separated numbers in row-major order, and line breaks within the
body are irrelevant — the parser reads tokens, not lines.

### Conventions

- **Row 0 is the northern row**, column 0 the western column, as the format
  specification requires. `Grid.cell_center(row, col)` returns map coordinates
  accordingly.
- **Cells are square.** Non-square cells are not representable in this format.
- **Nodata becomes `None` in memory**, never a number. Statistics, focal
  operations and cell-wise arithmetic all propagate the absence rather than
  averaging in a `-9999`.
- **Every number must be finite.** `inf` and `nan` parse as floats but are not
  measurements, and they travel unnoticed into comparisons that are always
  false and into JSON no strict parser will read back, so they are rejected on
  input.
- **Map units are assumed metric.** Areas are reported in m² and hectares.
- A projected coordinate system is assumed; degrees will produce meaningless
  areas and distances. The format carries no CRS, so nothing can be checked:
  a grid read from `.asc` has `crs = None`, and writing a grid that carries a
  CRS to `.asc` drops it — no `.prj` sidecar is written. GeoTIFF (below) is
  the format that preserves and validates the CRS.

### Reading and writing

```python
from silvispect.grid import Grid

grid = Grid.read("chm.asc")
print(grid.nrows, grid.ncols, grid.cellsize, grid.extent)
print(grid.stats().as_dict())

grid.smooth_median(1).write("smoothed.asc", precision=2)
```

`Grid` also offers `from_rows`, `filled`, `map_values`, `combine`, `clip`,
`focal`, `window`, `neighbors`, `histogram` and `cell_of` / `cell_center` for
coordinate conversion. Trailing zeros are trimmed on write to keep files small.

Files produced by Silvispect are readable by GDAL (`AAIGrid` driver), QGIS,
GRASS and R's `terra`.

## Rasters — GeoTIFF (`.tif`, `.tiff`)

`Grid.read` and `Grid.write` choose the format from the file extension,
case-insensitively: `.asc` is ASCII Grid, `.tif` and `.tiff` are GeoTIFF.
Any other extension — or none — is refused with a `GridError` before the
file is touched, so an ASCII grid saved under a name like `chm.txt` is no
longer read as a raster.

### Reading

Classic TIFF only, single band, stored as strips or tiles. Compression may
be none, PackBits, or Deflate (codes 8 and 32946), with predictors 1, 2 and
3. Sample formats are unsigned and signed integers of 8, 16, 32 or 64 bits,
and IEEE floats of 32 or 64 bits. Anything outside that scope — LZW or JPEG
compression, BigTIFF, multiple bands, a rotated or sheared raster, GDAL
scale/offset metadata, sparse segments — is refused with a `GridError`
naming the unsupported feature, never read approximately. Rasters over
50,000,000 cells are refused before any memory is allocated.

Georeferencing comes from ModelPixelScale + ModelTiepoint or from a
ModelTransformationTag, and must be north-up with square cells. If a file
carries no GeoTIFF georeferencing tags at all, a `.tfw` world file beside it
is read instead; the world file must be axis-aligned (its rotation terms
zero) and is interpreted on the pixel-centre convention. A file with only
some of the GeoTIFF tags is corrupt, not a candidate for the fallback.

### The CRS

The CRS is **preserved and validated, never reprojected**. A projected CRS
is kept on the grid as `Grid.crs`, the string `"EPSG:<code>"`, written back
to the GeoKey directory on write, and carried through derived grids; two
grids whose CRS values are both known and different refuse to combine. A
geographic or geocentric coordinate system (degrees) is refused, as is a
user-defined CRS (code 32767), which cannot be kept as an EPSG code. When a
file states its linear unit it must be metres. Beyond that stated unit and
a small list of known geographic codes, whether an EPSG code is a metric
projected system is *not* verified — Silvispect carries no EPSG database,
so an unusual code is trusted as written.

### Writing

Output is an uncompressed classic TIFF and is always lossless: float32 is
written only when every cell value and the nodata sentinel are exactly
representable in float32, and float64 otherwise, so a grid read back is the
grid that was written, cell for cell. Values originating in a three-decimal
`.asc` file are usually not float32-exact and are therefore written as
float64 — still a valid GeoTIFF, though not the 32-bit form the Forestry
Agency standard specifies; no quantising write option is offered, because
one would break the rule that the same data analyses identically in either
format. The nodata sentinel is always written as a `GDAL_NODATA` tag. The
`precision` argument applies to `.asc` output only and is ignored for
GeoTIFF. Files produced by Silvispect are readable by GDAL (the `GTiff`
driver), QGIS, GRASS and R's `terra`.

## Inventories — CSV

```csv
tree_id,x,y,species,dbh_cm,height_m,crown_diameter_m,status
T0001,612512.34,6642031.87,PIAB,32.4,23.1,4.12,live
T0002,612518.02,6642029.55,FASY,27.9,21.6,3.88,live
T0003,612521.77,6642040.10,BEPE,,,,dead
```

### Columns

| Column | Required | Meaning |
| --- | --- | --- |
| `x`, `y` | **yes** | Stem position in map units, same system as the raster |
| `tree_id` | no | Identifier; generated sequentially when absent |
| `species` | no | Species code, uppercased on use |
| `dbh_cm` | no | Diameter at breast height (1.3 m), centimetres |
| `height_m` | no | Total height, metres |
| `crown_diameter_m` | no | Crown diameter, metres |
| `status` | no | `live` (default), or any other value for a non-live stem |

Only `x` and `y` are required. Anything missing is carried as `None` and
*reported* by the inspection rules rather than imputed — SV040 exists precisely
so a gap in the data is visible instead of invented.

### Aliases

Column names are normalised (lowercased, spaces and hyphens to underscores) and
these aliases are accepted:

| Canonical | Also accepted |
| --- | --- |
| `tree_id` | `id`, `treeid`, `tree`, `stem_id` |
| `x` | `easting`, `x_m` |
| `y` | `northing`, `y_m` |
| `species` | `sp`, `species_code` |
| `dbh_cm` | `dbh`, `d`, `diameter`, `diameter_cm` |
| `height_m` | `height`, `h`, `ht` |
| `crown_diameter_m` | `crown_diameter`, `cd`, `crown_width_m` |
| `status` | `state` |

### Missing values

Empty cells and the literals `NA`, `NaN`, `NULL` and `-` (case-insensitive) all
parse as missing. A non-numeric value anywhere else is an **error**, not a
silent `None`: `row 7: dbh_cm value 'approx 30' is not numeric`. A numeric but
non-finite value is an error too — `inf` is not a diameter. Blank rows are
skipped.

A file with a header and no rows is a valid, *empty* inventory: somebody counted
and found nothing. That is different from supplying no inventory at all, and the
report says so — `metrics_source: field` with a count of zero, rather than
`none`.

### Live status

`live`, `living`, `alive`, `l` and `1` count as living. Everything else —
`dead`, `snag`, `stump`, `windthrow` — is a non-live stem, excluded from
metrics by default and exempt from the missing-measurement rule.

### Reading and writing

```python
from silvispect.inventory import read_trees, write_trees

trees = read_trees("trees.csv")
live = [t for t in trees if t.is_live]
write_trees(live, "live-only.csv", precision=2)
```

Output always uses the canonical column order, so a read/write round trip
normalises a messy file.

**Inventory CSV is a quantised record, not an exact copy.** Coordinates and
measurements are written to `precision` decimals — three by default, a
millimetre, which is finer than any field instrument reports and is the
resolution a stem position is actually known to. Writing and reading back is
therefore *not* an identity on the numbers: a stem at `x = 2.5004` is written
as `2.5`, so a match tolerance of exactly `2.5 m` excludes it before the round
trip and admits it afterwards. Raise `precision` if a downstream comparison
needs a finer record than a millimetre.

That resolution bounds what a diameter *class* can mean. `class_width` is
refused below **0.001 cm**, a hundredth of a millimetre: a band narrower than
the record it is cut from separates stems the file cannot tell apart, and the
decimal arithmetic behind the binning is sized from that floor rather than left
open-ended. Ordinary widths — 0.1, 1, 2, 5 cm — are unaffected.

This is the opposite of the choice `silvispect chm` makes, and deliberately so.
A derived canopy height model is an *intermediate* that the next command
analyses, so it defaults to text that reads back as the same float; an
inventory is a *record of measurements*, and writing more decimals than the
measurement has claims a precision the field data does not carry.

## Threshold profiles — JSON

A flat object of `InspectionConfig` field names to numbers. Any subset may be
given; the rest keep their defaults.

```json
{
  "canopy_threshold": 2.0,
  "min_canopy_cover": 0.7,
  "max_gap_area": 250.0,
  "min_gap_width": 5.0,
  "height_residual_z": 2.0,
  "min_recall": 0.8
}
```

Unknown keys, non-numeric values, negative thresholds, proportions outside
`[0, 1]` and self-contradictory pairs are all rejected when the profile loads,
so a broken profile fails loudly instead of quietly leaving a default in place.
See
[`data/thresholds-strict.json`](../data/thresholds-strict.json) for a complete
profile and [`findings.md`](findings.md) for what each threshold controls.

## Report JSON

`silvispect inspect --format json` emits a single object:

```jsonc
{
  "area_ha": 0.25,
  "metrics_source": "field",          // "detected", or "none" when nothing counted stems
  "summary": { "finding_count": 2, "max_severity": "notice", "counts": {...} },
  "metrics": { "stems_per_ha": 340.0, "basal_area_per_ha_m2": 25.11,
               "yield_basis_count": 85, ... },
  "canopy": { "cover": 0.636, "rugosity_m": 4.2, "strata": {...}, ... },
  "findings": [ { "code": "SV012", "severity": "notice", "title": ..., ... } ],
  "detection": { "tree_count": 84, "trees": [...] },
  "match": { "recall": 0.988, "precision": 1.0, ... },
  "allometry": { "model": "naslund", "params": {...}, "by_species": {...} },
  "gaps": [ { "gap_id": 1, "area_m2": 164.0, "width_m": 12.3, ... } ],
  "config": { ... }
}
```

`detection`, `match`, `allometry` and `gaps` are present only when the
corresponding input was supplied. The output is strict JSON: absent statistics
are `null`, never `NaN` or `Infinity`, so any parser can read it. That holds
for every command that writes JSON, and it is enforced rather than assumed —
a number with no finite value stops the command instead of being printed.

`recall`, `precision` and `f1` count only the field stems standing on a cell
the canopy model actually measured; a stem in a nodata hole is not evidence
about the detector either way.

`recall`, `precision` and `f1` are `null` when their own denominator is empty:
with no field stems there is no recall to report and with no detections there
is no precision, and neither is zero. The agreement rules read them the same
way and stay silent.

`yield_basis_count` is the number of stems behind `volume_per_ha_m3` and
`biomass_per_ha_t`. Both need a height as well as a diameter, so the totals
cover a subset of `measured_count` and are `null` when that subset is empty —
a stand with diameters but no heights has an unknown volume, not a zero one.
