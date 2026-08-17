# Data formats

Silvispect reads and writes plain text only. Both formats are widely supported
interchange formats, so data can move in and out of GIS and statistical
software without a conversion step in the middle.

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
- **Map units are assumed metric.** Areas are reported in m² and hectares.
- A projected coordinate system is assumed; degrees will produce meaningless
  areas and distances.

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
silent `None`: `row 7: dbh_cm value 'approx 30' is not numeric`. Blank rows are
skipped.

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
are `null`, never `NaN` or `Infinity`, so any parser can read it.

`yield_basis_count` is the number of stems behind `volume_per_ha_m3` and
`biomass_per_ha_t`. Both need a height as well as a diameter, so the totals
cover a subset of `measured_count` and are `null` when that subset is empty —
a stand with diameters but no heights has an unknown volume, not a zero one.
