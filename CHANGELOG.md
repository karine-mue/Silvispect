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
- 250 tests, ruff lint and format configuration, mypy, GitHub Actions CI across
  Python 3.10–3.13, and a committed sample plot whose reproducibility CI
  verifies.
- CI builds both distributions, installs the wheel, checks the `py.typed`
  marker ships, and runs the whole suite from the unpacked source
  distribution — the packaging claims are enforced, not asserted.
- PEP 561 `py.typed` marker, and a source distribution that carries the test
  fixtures the shipped tests need.

### Corrections before release

Four rounds of independent review found twenty behavioural defects before
release. Each is fixed with a regression test that fails against the code as it
stood; where the defect violated a general property — monotonicity, spatial
symmetry, permutation invariance — the property itself is now the test.

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

Third round — a property check over random rasters found the gap filter was
not monotone, plus four narrower issues:

- **`min_gap_width` was not monotone.** Erosion followed by dilation mixes two
  quantisations, and the opened area oscillated (64, 96, 64, 88, 60 cells on a
  treeless 10 m plot as the width rose a metre at a time), so *tightening* the
  criterion could enlarge a gap and raise SV011 where a looser setting did not.
  Replaced by an opening function — one scalar field per cell, thresholded —
  which is nested by construction. Verified against 800 random rasters.
- Judging edge contact by proximity to the border, introduced in round two,
  flagged interior openings a single tree-row from the edge. Edge contact is
  now inherited from the untrimmed sub-canopy component, which is exact.
- An inventory with a header and no rows was treated as no inventory at all by
  `inspect` and as a counted empty stand by `metrics`. It is data either way:
  somebody counted and found nothing. All paths now agree on
  `metrics_source: field` with a count of zero.
- Non-finite numbers were rejected in threshold profiles but accepted from
  rasters, inventories and CLI options, where they produced tracebacks, exit 1,
  or `NaN` in JSON. Every input surface now rejects them with exit 2 — as does
  `inspect` with no arguments, which previously exited 1, the code reserved for
  findings.
- The synthetic generator's contract claimed the returned trees were exactly
  the stems rendered into the canopy model. The surface keeps the maximum over
  overlapping crowns, so a shorter neighbour can be painted and then hidden —
  which is what a real canopy model does. The claim is corrected rather than
  the code, and a test pins the behaviour at both coarse and fine resolutions.
- The CI check that `py.typed` ships was importing the checkout rather than the
  installed wheel, so it would have passed with the marker missing from the
  wheel. It now runs from a neutral directory and asserts the import resolved
  to site-packages.

Fourth round — a property campaign over 914 randomized stands applied
transformations (translation, reflection, rotation, permutation) and checked
what should be invariant under them:

- **The opening function was not spatially equivariant.** Discs were skipped
  when something had already covered their centre, but covering a centre is not
  containment — that needs `|p - q| + r <= R`. Equal-radius discs therefore
  dropped each other in row-major order: on a treeless 4 m x 4 m plot the four
  equal discs collapsed to the north-west 3x3 block, so reflecting the input
  changed the answer (199 of 1,200 comparisons), and raising `canopy_threshold`
  — which can only admit more open ground — lost cells in 1,302 of 1,920
  transitions. The containment margin is now tracked properly. Round three's
  monotonicity in `min_gap_width` was real, but it thresholded a field whose
  construction was itself incomplete.
- **Matching depended on CSV row order at equal distances.** Two crowns and two
  stems all exactly one metre apart matched one pair in one row order and two in
  the other, moving recall from 0.5 to 1.0 and RMSE from 1.0 m to 16.2 m on
  identical geometry. Ties now break on position and identifier.
- Dominant height had the same defect at the 100-stems-per-hectare cutoff: two
  30 cm stems of 10 m and 30 m gave whichever was listed first. Ranking is now
  diameter, then height, then identifier.
- `silvispect chm` wrote the derived model at three decimals with no way to
  change it, turning near-equal heights into exact plateaus and changing which
  cells the detector called apexes — 267 of 300 randomized stands differed from
  the same analysis run on the unrounded model. It now writes six decimals by
  default and takes `--precision`.
- Fitted allometry parameters drifted in the sixth decimal under permutation,
  because floating-point summation order changes the optimiser's path.
  Observations are sorted before fitting.

The changelog previously claimed the first five defects were "all fixed with
regression tests". Three were only partly fixed, and the tests pinned the
reproductions rather than the general cases — which is how each later round
happened. Where a general property exists, it is now tested as one: opening
monotonicity and spatial equivariance, matching and fitting under permutation,
and threshold inclusion are all checked over randomized inputs rather than a
single hand-built example.
