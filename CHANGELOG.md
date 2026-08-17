# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[semantic versioning](https://semver.org/spec/v2.0.0.html).

## 0.1.0 — 2026-08-17

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
- 236 tests, ruff lint and format configuration, mypy, GitHub Actions CI across
  Python 3.10–3.13, and a committed sample plot whose reproducibility CI
  verifies.
- CI builds both distributions, installs the wheel, checks the `py.typed`
  marker ships, and runs the whole suite from the unpacked source
  distribution — the packaging claims are enforced, not asserted.
- PEP 561 `py.typed` marker, and a source distribution that carries the test
  fixtures the shipped tests need.

### Corrections before release

Two rounds of independent review found nine behavioural defects before release.
Each is fixed with a regression test that fails against the code as it stood.

First round — absent data was being turned into confident numbers:

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

Second round — three of those fixes were incomplete, and the re-review found
the remaining halves:

- The gap-width fix only covered rasters with *no* canopy. Wherever any canopy
  existed the plot edge was still ignored, so an opening running off the side
  could claim an inscribed circle larger than the raster: a 10 m x 10 m plot
  with one corner tree reported a 25 m width. The distance field is now bounded
  by the extent in every case, which subsumes the no-canopy special case.
- Because the morphological opening insets the mask from the border, an
  edge-touching gap lost its `touches_edge` flag to the erosion and was rated
  `critical` instead of `warning`. Edge contact is now judged by proximity to
  the border, using the same radius the opening removed.
- Drawing stems inside the realised raster was not enough to make every
  ground-truth tree visible: a crown narrower than its own cell fell between
  cell centres and painted nothing, so 23 of 30 truth trees left no trace at a
  10 m cell size. The painting radius now has a half-diagonal floor.
- `metrics_source: none` still published `tree_count: 0` and
  `stems_per_ha: 0.0`, which is the same "unknown becomes zero" mistake one
  layer further out. Those figures are now `null`, via `StandMetrics.unknown`.
- Threshold validation rejected ordinary bad values but not infinities
  (`1e999` raised an uncaught `OverflowError`) or fractional counts (`1.9` was
  silently truncated to `1`). Both are now rejected with exit status `2`.
- `yield_basis_count` reached the JSON but not the text or Markdown reports, so
  the human-facing output could show a partial volume without its sample size.

The changelog previously claimed the first five were "all fixed with regression
tests". Three of them were only partly fixed, and the tests covered the
reproductions rather than the general cases — hence this second round.
