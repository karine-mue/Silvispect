# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[semantic versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] — 2026-08-17

First release. Everything below is new.

### Rasters and canopy

- `Grid`: ESRI ASCII Grid reader and writer with nodata carried as `None`,
  cell/map coordinate conversion, focal operations (mean, median, max), cell
  algebra, statistics and histograms.
- `canopy_height_model` derives a CHM from a surface and a terrain model,
  clamping normalisation artefacts rather than discarding measured cells.
- Canopy cover, open fraction, rugosity and vertical strata.
- Gap mapping with the **Brokaw minimum-width criterion**, implemented as a
  morphological opening of the sub-canopy mask over a chamfer distance
  transform. Without it, the ground visible between crowns on a fine raster
  labels as one enormous spurious opening.

### Detection

- Variable-window local maxima, with the search radius a linear function of
  cell height.
- Marker-controlled crown region growing that stops where the surface rises
  again, placing boundaries in the valley between adjacent crowns.
- Crown geometry, a crown label raster, and deterministic tallest-first tree
  numbering.
- Greedy one-to-one matching of detections against field stems, reporting
  recall, precision, F1, height bias, height RMSE and positional offset.

### Inventory and allometry

- `Tree` records with CSV interchange, column aliases, tolerant missing-value
  handling and strict errors on non-numeric data.
- Four height-diameter families (Näslund, Curtis, power, Chapman-Richards),
  fitted by a self-contained Nelder-Mead simplex search, with inversion by
  bisection and residual z-score screening.
- **Per-species fitting**, because a pooled curve on a mixed stand
  systematically flags the minority species as outliers.
- Tabulated Näslund defaults for six European species.

### Metrics

- Basal area, stems per hectare, quadratic mean diameter, Lorey's mean height,
  dominant height, Reineke SDI, Gini coefficient of basal area, Shannon and
  Simpson diversity, diameter distribution, form-factor volume and
  Chave-form above-ground biomass.

### Inspection and reporting

- A rule registry producing coded findings (SV001–SV052) at four severities,
  with every threshold a documented, JSON-overridable field of
  `InspectionConfig`.
- Markdown, JSON and terminal report renderers, plus ASCII canopy maps that
  downsample by block maximum so treetops survive.
- Report JSON is strict: absent statistics are `null`, never `NaN`.

### Tooling

- `silvispect` CLI with ten subcommands, uniform exit codes, and
  `inspect --fail-on <severity>` for use as a CI gate.
- Seeded synthetic stand generator providing ground truth for the test suite
  and for `silvispect demo`.
- 227 tests, ruff lint and format configuration, GitHub Actions CI across
  Python 3.10–3.13, and a committed sample plot whose reproducibility CI
  verifies.
- PEP 561 `py.typed` marker, and a source distribution that carries the test
  fixtures the shipped tests need.

### Corrections before release

An independent review of the first pushed commit found five behavioural
defects, all fixed here with regression tests:

- Volume and biomass totals reported `0` for a stand with diameters but no
  heights. They are now absent, and `yield_basis_count` exposes how many stems
  the totals actually rest on.
- `inspect --no-detection` invented a zero-tree "detected" stand and raised
  SV001, failing a `--fail-on warning` gate for a stand nobody had counted.
  Metrics now report `metrics_source: none` and the stem-based rules stay quiet.
- A raster with no canopy at all produced an infinite gap width, which is not
  representable in JSON. With no canopy to measure from, distances are now
  measured from the plot edge.
- A malformed threshold profile (a JSON `null`, a list value, an out-of-range
  or self-contradictory threshold) surfaced as a traceback. Profiles are now
  validated on load and rejected with exit status `2`.
- The synthetic generator drew stems across the *requested* extent while
  rendering the *rounded* raster, so a non-divisible extent returned ground-truth
  trees that were never rendered. Stems are now drawn inside the realised
  raster.

[0.1.0]: https://github.com/karine-mue/Silvispect/releases/tag/v0.1.0
