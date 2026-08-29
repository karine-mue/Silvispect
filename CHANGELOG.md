# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[semantic versioning](https://semver.org/spec/v2.0.0.html).

Sections below the newest one describe a release as it was published. Defects
found *after* a release are recorded under the version that fixes them, never
folded back into the version that shipped with them.

## 0.1.1 — unreleased

Behavioural fixes on top of the released 0.1.0. Nothing here was present in
0.1.0; each item describes a defect that shipped in it.

Numbered as a patch release: nothing is removed or renamed, and every entry is
a correction to behaviour that was wrong. Two of them are visible in output —
`silvispect chm` writes heights that round-trip instead of rounding to three
decimals, and gains an optional `--precision` to round on request — but neither
is a new capability. The rounding was the defect.

Three rounds of review after the release found these. As before, each fix
carries a test that fails against the code as it stood — and where the defect
violated a general property, that property is what the test asserts, over
randomized inputs rather than the one raster that exposed it.

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
  transitions. 0.1.0's monotonicity in `min_gap_width` was real, but it
  thresholded a field whose construction was itself incomplete.
- **Matching depended on CSV row order at equal distances.** Two crowns and two
  stems all exactly one metre apart matched one pair in one row order and two in
  the other, moving recall from 0.5 to 1.0 and RMSE from 1.0 m to 16.2 m on
  identical geometry.
- Dominant height had the same defect at the 100-stems-per-hectare cutoff: two
  30 cm stems of 10 m and 30 m gave whichever was listed first. Ranking is now
  diameter, then height, then identifier.
- `silvispect chm` wrote the derived model at three decimals with no way to
  change it, turning near-equal heights into exact plateaus and changing which
  cells the detector called apexes — 267 of 300 randomized stands differed from
  the same analysis run on the unrounded model.
- Fitted allometry parameters drifted in the sixth decimal under permutation,
  because floating-point summation order changes the optimiser's path.
  Observations are sorted before fitting.

Fifth round — a review of those repairs found that three of them had replaced
one dependence with another, and that a fourth had made the opening function
too slow to use:

- **Matching had traded row order for map orientation.** Breaking equal-distance
  ties on coordinates makes the answer depend on which way the plot faces:
  rotating the same two crowns and two stems by 180 degrees moved recall from
  0.5 to 1.0 and RMSE from 1.0 m to 16.2 m — the exact defect the tie-break was
  meant to remove. Equal-distance candidates are now resolved together, by a
  maximum bipartite matching over that one distance, so the number of pairs a
  tie yields is a fact about the geometry. Nearer pairs are settled first and
  are never given up. Checked under six isometries, including translation,
  rotation, reflection and transposition.
- **Edge classification ignored the requested connectivity.** `find_gaps`
  labelled components with the caller's connectivity but always inherited edge
  contact from eight-connected components, so under `connectivity=4` a diagonal
  contact carried the edge flag inward: on a 3x3 raster the interior cell was
  reported as touching the border and `include_edge_gaps=False` discarded it.
  Both passes now use the same connectivity, and a randomized check asserts that
  a gap is flagged as an edge gap exactly when it contains a border cell.
- **The opening tolerance was a fixed number of metres.** Comparisons against
  the radius field absorbed a `1e-9` slack, which is invisible on a metre-wide
  cell and total on a nanometre-wide one: at that cell size the four equal discs
  of a treeless 4x4 plot resolved to two different radii and the width filter
  admitted every cell whatever width was asked for. The tolerance is now a
  fraction of the cell size, so the field is covariant with the plot's scale —
  checked across cell sizes from 1e-9 to 1e6.
- **`silvispect chm` had moved the rounding boundary, not removed it.** Writing
  six decimals instead of three still merges heights that differ below the sixth
  into a plateau, and a plateau has no apex: heights differing at 1e-8 m again
  gave a different tree count than the same analysis on the unrounded surface.
  The command now writes the shortest text that reads back as the same float, so
  the derived model *is* the model the analysis would have used, at any
  magnitude. `--precision` still rounds on request.
- **The corrected opening function was superlinear.** Painting every maximal
  disc cell by cell cost 8.6 s on a treeless 384 x 384 raster, against 0.16 s
  for the 0.1.0 code — which was fast only because it was dropping discs it
  should have kept. Discs contained in a neighbour's are now dropped by the
  honest test `|p - q| + d(p) <= d(q)`, which reads only the distance field and
  so cannot depend on scan order, and the survivors settle each cell once
  through a per-row disjoint-set walk. The same raster takes 0.59 s and the cost
  is now linear in the number of cells. A brute-force comparison over 400
  randomized rasters confirms the field is unchanged.

Sixth round — an exhaustive campaign found that the matching repair was still
incomplete, and turned up five defects that have been there since 0.1.0:

- **Matching still let identifiers decide how many pairs it found.** Settling
  one distance at a time and never looking further is not enough to honour the
  documented objective. Two stems a metre either side of a crown are an
  arbitrary choice on their own, but taking the wrong one strands a second
  crown two metres from the left-hand stem: relabelling the two stems, with the
  geometry untouched, moved the match count from one to two. Matching now
  reaches the stated optimum exactly — as many pairs as possible at the
  shortest distance, then as many as possible at the next without giving any of
  those up — by encoding it as an integer weight per distance and solving a
  minimum-cost assignment over each connected cluster of candidates. Checked
  against brute-force enumeration of every possible pairing on randomized
  geometries. Only the *count* is fixed by geometry; where the geometry is
  symmetric several pairings are equally optimal, and identifiers still choose
  between them, which is now stated rather than incidental.
- **Materially misaligned rasters were combined without complaint.** The origin
  check kept a relative tolerance, so the allowance grew with the coordinates:
  at an easting of 1e9 two rasters half a cell apart were subtracted cell for
  cell into a quietly wrong canopy model. Origins must now agree to within a
  millionth of a cell, which means the same thing at every magnitude.
- **A one-cell raster could take sixteen seconds to detect one tree.** The
  search window grows with the height of the cell being tested, so an accepted
  height of 1e5 m asks for a radius of 5,501 cells, and the window walked the
  whole requested square before discarding everything outside the raster. The
  offsets are now clipped first, so the cost is the number of cells actually
  looked at. The cells yielded, and their order, are unchanged.
- **A plot exactly 30 % open was reported as above a 30 % limit.** The open
  share was taken as one minus the cover, and `1.0 - 0.7` is not `0.3`. It is
  now counted from the open cells directly, so a rule phrased as *above* a
  limit reads strictly, and both numbers no longer print as "30%" beside a
  finding that says one exceeds the other.
- **Values close to the nodata sentinel were read back as absent.** Absence was
  matched approximately, so with a sentinel of `0` the smallest positive number
  a raster can hold vanished on the way back in, and the band of swallowed
  values widened with the magnitude of the sentinel. The sentinel is written
  exactly and is now matched exactly.
- **`chm` wrote its output and then failed describing it.** A raster holding
  `0` and `1e200` is finite and is accepted everywhere, but squaring the
  deviations overflowed while computing display statistics, so the command
  reported failure over a file it had already written correctly. Statistics are
  now scaled before squaring and summed exactly.
- Fractional raster dimensions — `ncols 1.9` — were truncated to a smaller
  raster instead of being rejected as malformed.
- Lorey's mean height depended on row order for extreme inputs, because a
  running total loses the small stems once a large one has been added. Summed
  exactly instead.

Default inventory CSV serialisation was examined and **left unchanged**. It
rounds to three decimals, so a stem four ten-thousandths outside a match
tolerance is inside it after a round trip. That is the intended contract rather
than a defect: `precision` is a documented parameter, an inventory is a record
of measurements rather than an intermediate to be re-analysed, and a millimetre
is finer than any field instrument reports. `docs/data-formats.md` now says so
outright, and a test pins it.

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
- 242 tests, ruff lint and format configuration, mypy, GitHub Actions CI across
  Python 3.10–3.13, and a committed sample plot whose reproducibility CI
  verifies.
- CI builds both distributions, installs the wheel, checks the `py.typed`
  marker ships, and runs the whole suite from the unpacked source
  distribution — the packaging claims are enforced, not asserted.
- PEP 561 `py.typed` marker, and a source distribution that carries the test
  fixtures the shipped tests need.

### Corrections before release

Three rounds of independent review found fifteen behavioural defects before
release. Each is fixed with a regression test that fails against the code as it
stood.

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

The changelog previously claimed the first five defects were "all fixed with
regression tests". Three were only partly fixed, and the tests pinned the
reproductions rather than the general cases — which is how rounds two and three
happened. Where a general property exists, it is now tested as one.
