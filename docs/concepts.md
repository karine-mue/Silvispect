# Concepts

What Silvispect computes, how it computes it, and what each number assumes.
This page is the reference for the science; [`findings.md`](findings.md) covers
the rules built on top of it and [`cli.md`](cli.md) the commands.

## 1. Surfaces

Three rasters describe a stand vertically.

| Surface | Meaning |
| --- | --- |
| **DTM** — digital terrain model | Elevation of the bare ground |
| **DSM** — digital surface model | Elevation of the topmost reflecting surface |
| **CHM** — canopy height model | `DSM − DTM`: vegetation height above ground |

```python
from silvispect.canopy import canopy_height_model

chm = canopy_height_model(dsm, dtm)
```

Negative results are clamped to zero rather than discarded. A cell that
normalises to −0.2 m is a *measured, treeless* cell, and throwing it away would
inflate every cover statistic computed afterwards. Cells with no data in either
input stay nodata in the result: Silvispect never lets a `-9999` sentinel enter
an average.

**Everything downstream reads the CHM.** If the terrain model is wrong, every
height is wrong by the same amount, which is exactly what finding SV052
(systematic height bias) is there to catch.

## 2. Canopy cover, gaps and strata

**Canopy cover** is the fraction of valid cells at or above a height threshold,
2 m by default — the conventional break between regeneration and canopy.

**Gap mapping** is where a naive implementation goes wrong. Thresholding the
CHM and labelling connected components seems obvious, but on a 0.5 m raster the
ground visible *between* neighbouring crowns forms one connected web across the
whole plot. Labelling it reports a single enormous "gap" that is really an
artefact of resolution:

```console
$ silvispect gaps data/plot-a1-chm.asc --min-width 0
3 gaps at or above 25 m2 and 0.0 m wide; 36.4% of the plot is below 2.0 m
  gap 1        511 m2  width  11.8 m  ...   <- one opening fused to the web
```

Silvispect applies the **Brokaw criterion**: a gap must be wide enough to
contain a circle of a minimum diameter (5 m by default). It is implemented as a
morphological opening of the sub-canopy mask, which removes the web *including
the tentacles attached to a genuine opening* — something a simple "discard small
components" filter cannot do. Labelling then runs on the opened mask:

```console
$ silvispect gaps data/plot-a1-chm.asc
3 gaps at or above 25 m2 and 5.0 m wide; 36.4% of the plot is below 2.0 m
  gap 1        161 m2  width  11.8 m  ...   <- the opening alone
```

Both report the *same* inscribed width: the difference is entirely tentacles of
inter-crown ground, which is why a "discard narrow components" filter cannot
remove it and an opening can.

Each gap reports its area, its **width** (the diameter of the largest inscribed
circle), its equivalent radius, its centroid, and whether it touches the plot
edge — an edge gap is truncated by the extent, so its area is a lower bound
and Silvispect downgrades its severity accordingly.

The opening is computed as an **opening function**: one pass records, for every
cell, the radius of the largest disc that fits in the sub-canopy area and still
covers that cell, and the mask is a threshold on that field.

That shape is deliberate, and two earlier attempts at it were wrong in
instructive ways.

*Erode by the radius, then dilate by it* is the textbook construction, but it
mixes two quantisations — a cut on the distance field, then a fresh chamfer walk
out of the eroded core — and the result is not monotone. On a treeless
10 m × 10 m plot the reported opening ran 64, 96, 64, 88, 60 cells as the
requested width rose a metre at a time, so *tightening* `min_width` could
enlarge a gap and raise a finding a looser setting did not. A threshold on a
single scalar field cannot behave that way: a larger radius always selects a
subset.

Painting the field needs care too. A disc is dropped only when it is genuinely
*contained* in another — `|p − q| + r ≤ R`, not merely "something covered p's
centre". Skipping on coverage alone drops discs that reach cells no other disc
does, which left the field asymmetric: on a treeless 4 m × 4 m plot the four
equal discs collapsed to the north-west 3×3 block, so reflecting the input
changed the answer, and raising `canopy_threshold` — which can only admit *more*
open ground — could lose cells.

Containment is tested against the eight adjacent cells, reading the distance
field alone: `d(q) ≥ d(p) + |p − q|` for a neighbour `q` means the disc at `p`
is inside the disc at `q`. Because the test never consults the work done so
far, it cannot depend on the order cells are visited in — the failure mode the
coverage test had. On open ground it removes nearly every disc, since almost
all of them sit on a chamfer-shortest path towards the boundary. What survives
is painted largest-first, so the first disc to reach a cell settles it; already
settled cells are stepped over through a per-row disjoint-set walk rather than
rewritten. Painting the survivors cell by cell instead cost 8.6 s on a treeless
384 × 384 raster where this costs 0.59 s, and the cost is now linear in the
number of cells.

One tolerance is unavoidable — the distance field accumulates step lengths, so
exact comparisons against it need slack — and it is expressed as a fraction of
the cell size rather than as a fixed number of metres. A fixed `1e-9` is
invisible on a metre-wide cell and total on a nanometre-wide one, where it made
the field lose its scale covariance and the width filter admit everything.

The **connectivity** asked for is used consistently. Gaps are labelled with it,
and edge contact — inherited from the untrimmed sub-canopy component, since the
opening insets the mask away from the border — is judged on components built
with the same connectivity. Using eight-connectivity there regardless let a
diagonal contact carry the edge flag into an opening that four-connectivity had
deliberately kept separate.

Distances are measured to the canopy's **boundary**, not to the centre of the
cell beyond it. Half of the step from an open cell's centre to a canopy cell's
centre lies inside the canopy, so it is not clearance: counting it let a
one-metre square opening claim a two-metre inscribed circle and survive a
`min_width` of 1.5 m that it should have failed.

Two consequences of "as mapped" are worth stating plainly. The **plot edge is a
boundary** for the inscribed circle too, measured the same way — half a cell
out from the outermost centre: nothing is known beyond the extent, so a
gap running off the side is measured only within it. Without that bound a
10 m x 10 m plot with a single corner tree reported a 25 m gap width — larger
than the raster's own diagonal. And because the opening insets the mask
from the border, edge contact is inherited from the *untrimmed* sub-canopy
component the gap was carved from: the trimmed remnant keeps no cell on the
border, and proximity is not a usable stand-in, since an interior opening one
tree-row from the edge is just as close.

**Vertical strata** bin the canopy cells into ground / regeneration /
understorey / midstorey / canopy / emergent layers, and **rugosity** is the
standard deviation of canopy heights: a compact proxy for vertical complexity
that is low in even-aged plantations and high in multi-layered stands.

## 3. Individual tree detection

Two classical steps, both dependency-free.

### Variable-window local maxima

Crown width grows with tree height, so a fixed search window is wrong at both
ends of the size distribution: too large and it merges small neighbouring
trees, too small and it finds several maxima inside one large crown. The window
radius is therefore a linear function of the height of the cell being tested:

```
r(h) = window_intercept + window_slope · h        [metres]
```

with defaults `0.8 + 0.055·h`. A cell is a treetop if a **strictly** taller
cell lies within `r(h)`; equal cells never suppress one another, because which
of two equal maxima is scanned first is a fact about the array and not about
the forest — it reverses when the plot is mirrored.

A connected run of exactly equal height is one treetop, not one per cell. Its
apex is the member nearest the run's centre, measured in whole cells so that
the answer mirrors exactly; where several are equally central, the one whose
surroundings are denser wins, judged by the sorted list of (distance, height)
pairs seen from each — a quantity every rotation and reflection of a raster
preserves. Only members that even that cannot tell apart fall back on row and
column order.

The same ordering settles which crown claims a cell two equally tall apexes
both reach. **Ties are read from the raster as measured, not from the smoothed
surface**: the mean filter can make two candidates equal — that is what it is
for — but it must not also be what decides between them. On `2 3 / 4 9 / 3 2`
smoothing produces a surface that is its own mirror image, and nothing in it
distinguishes the two maxima that the raster distinguishes easily.

Detection is therefore equivariant: rotating or reflecting a plot rotates or
reflects the trees found in it. The exception is a plot that really is
symmetric — a flat run with no middle cell has no apex that mirrors onto
itself — where the count and the crown sizes are fixed but the choice between
two interchangeable cells cannot be.

The CHM is mean-filtered first (radius 1 cell by default). This suppresses the
single-cell spikes that would otherwise each become a tree, at the cost of a
small negative height bias, visible as `height_bias_m ≈ −0.7` in the sample
plot. Set `--smooth-radius 0` to keep raw apex heights.

### Marker-controlled region growing

Each maximum seeds a crown. Cells are absorbed from a priority queue,
tallest-first, and a candidate joins a crown only if:

- it is at or above the minimum tree height,
- it stays above `drop_fraction` of the apex height (0.45 by default), which
  stops crowns bleeding down into the understorey,
- it is not *taller* than the cell that recruited it — the surface must keep
  descending away from the apex, and
- it is within the maximum plausible crown radius for a tree of that height,
  `crown_intercept + crown_slope · h`.

The "must keep descending" rule is what places the boundary in the valley
between two adjacent crowns. Crowns smaller than `min_crown_cells` are
discarded as noise, and the survivors are renumbered tallest-first, so tree
IDs are stable across runs.

Reported crown geometry: area (cell count × cell area), equivalent radius and
diameter, mean height, and maximum extent from the apex.

## 4. Matching detections to field trees

`match_trees` pairs crowns with field stems within a tolerance radius under one
objective:

> Honour as many pairs as possible at the shortest distance present; then,
> without giving up any of those, as many as possible at the next distance; and
> so on.

That is a function of the geometry alone. It is deliberately *not* "as many
pairs as possible": one pair at half a metre outranks two at two metres,
because a match is a claim that the detector found a particular stem, and a
nearer claim is the better-evidenced one.

**How many pairs are reported is fixed by the geometry.** Only that. Where the
geometry is symmetric — two crowns and two stems mutually equidistant — several
pairings are equally optimal, they differ in which crown is credited with which
stem, and something outside the geometry has to choose. Identifiers do, because
they travel with the records through reordering and through any rigid motion of
the plot. Height bias and RMSE can therefore depend on identifiers in a
perfectly symmetric plot; the counts cannot.

Three narrower rules were tried first and all three failed the same way, by
letting something that is not geometry decide the count. Taking the earlier row
made the answer depend on how the CSV happened to be written: two crowns and
two stems all exactly one metre apart matched one pair in one row order and two
in the other, moving recall from 0.5 to 1.0 and RMSE from 1.0 m to 16.2 m on
identical geometry. Ranking ties by *position* fixed that and introduced the
same defect one step out — the answer then depended on which way the plot
faced, and rotating the same four points by 180° moved recall from 0.5 back to
1.0. Resolving each distance by a maximum matching *restricted to that
distance* fixed both and still fell short: a tie can be broken two ways that
both honour one pair now while only one of them leaves a partner free later, so
two stems a metre either side of a crown stranded a second crown two metres
away, and relabelling the stems moved the count from one to two.

Reaching the objective needs the whole cluster considered at once. Each
distance gets an integer weight — a base larger than any achievable pair count,
raised to a power that falls with distance, so no number of farther pairs can
outweigh a single nearer one — and the exact optimum is a minimum-cost
assignment. Candidates that share no stem cannot affect each other, so the
graph is split into connected clusters first; on real inventories those are
almost all a single crown and a single stem. The suite checks the result
against brute-force enumeration of every possible pairing, and checks the
counts under row permutation, identifier permutation, translation, rotation,
reflection and transposition.

Dominant height has an ordering rule of its own — diameter, then height, then
identifier — so a tie at the 100-stems-per-hectare cutoff does not depend on
which stem was typed first.

From the pairing come the accuracy statistics:

| Statistic | Meaning |
| --- | --- |
| **Recall** | Share of field stems that were detected |
| **Precision** | Share of detections with a field counterpart |
| **F1** | Harmonic mean of the two |
| **Height bias** | Mean of `detected − field` height |
| **Height RMSE** | Root mean square of the same errors |
| **Mean offset** | Mean planimetric distance between apex and stem |

Recall against a *full* field inventory is bounded above by what reaches the
canopy surface: a suppressed tree under a closed canopy is invisible to a CHM
at any resolution. A low recall is often a true statement about the sensor
rather than about the detector.

## 5. Height-diameter allometry

The relationship between diameter at breast height (DBH, 1.3 m) and total
height is used twice: to give a remotely detected tree a diameter it has no
stem for, and to screen field measurements against the stand's own curve.

Four families are available:

| Name | Equation |
| --- | --- |
| `naslund` (default) | `H = 1.3 + D² / (a + b·D)²` |
| `curtis` | `H = 1.3 + a·(D / (1 + D))^b` |
| `power` | `H = a·D^b` |
| `chapman` | `H = 1.3 + a·(1 − e^(−b·D))^c` |

Näslund's equation is the default because it behaves correctly at both ends:
it passes near the origin and saturates at `1.3 + 1/b²`, so the asymptotic
height is readable straight off the parameters. (Curtis with a small exponent
saturates so fast that a 7 cm stem is predicted at full canopy height — a real
trap when choosing starting values.)

Fitting is least squares by Nelder-Mead simplex search
([`optimize.py`](../silvispect/optimize.py)), restarted once from the incumbent
to tighten the simplex. No derivatives, no SciPy. Observations are sorted before
fitting: the objective is permutation-invariant in real arithmetic but not in
floating point, and without a canonical order the same stand reported parameters
that differed in the sixth decimal depending on row order.

**Curves are fitted per species.** A single pooled curve on a mixed stand
systematically flags the minority species as outliers — in a test stand of
spruce, beech and birch, the pooled fit produced twelve spurious height
outliers that per-species fits reduced to four genuine tail observations.
Species with too few paired measurements fall back to the pooled curve, and a
stand with too little data altogether falls back to tabulated defaults.

**Residual screening** z-scores the residuals against their own standard
deviation, so it measures departure from the stand's *own* curve rather than
from an external expectation. That is the right question for data quality: a
stand may be unusual, but a single stem should not be unusual within it.

## 6. Stand metrics

All per-hectare values are expansions onto an area the caller supplies.
Silvispect never guesses the plot area — an unstated expansion factor is the
classic way an inventory summary becomes quietly wrong.

| Metric | Definition | Note |
| --- | --- | --- |
| Basal area | `π·(D/200)²` per stem, summed | Cross-section at breast height, m² |
| Stems per hectare | `n / area_ha` | Live stems by default |
| Quadratic mean diameter | `√(Σ D² / n)` | Diameter of the tree of mean basal area; always ≥ the arithmetic mean |
| Lorey's mean height | `Σ(g·h) / Σg` | Basal-area weighted; robust to the many small stems that drag an arithmetic mean down |
| Dominant height | Mean height of the 100 thickest stems per hectare | "Top height"; the usual site-index basis |
| Reineke SDI | `N · (QMD/25)^1.605` | Density referenced to a 25 cm QMD |
| Gini coefficient | Of the basal-area distribution | 0 = perfectly uniform; > 0.5 = strongly size-differentiated, often uneven-aged |
| Shannon / Simpson | Of the species mixture | Natural log; Simpson is the Gini-Simpson form |
| Volume | `form_factor · g · h` | Form factor 0.5, **not** a taper function |
| Biomass | `0.0673·(ρ·D²·H)^0.976` kg | Chave et al. (2014) form, wood density by species code |

Volume and biomass are order-of-magnitude tools. Anyone needing merchantable
volume should substitute a regional taper equation.

Both need a height as well as a diameter, so they are computed over a subset of
the measured stems. That subset is reported as `yield_basis_count`, and the
totals are absent rather than zero when it is empty: a stand with diameters but
no heights has an unknown volume, and reporting `0` would be a fabrication.

## 7. Synthetic stands

The detector is tested against ground truth, which requires a stand whose trees
are known exactly. `silvispect.synth` generates one from a seed:

1. Cut circular gaps at random positions.
2. Reject-sample stem positions at a minimum spacing, skipping the gaps.
3. Draw diameters from a lognormal distribution, clamped to a plausible range.
4. Derive heights from the species' Näslund curve plus relative noise.
5. Paint paraboloid crowns onto the canopy surface, keeping the per-cell
   maximum, with crown radius scaled to height.
6. Add per-cell measurement noise, lay the canopy over a sloped and undulating
   terrain model, and subtract to obtain the CHM.

A raster holds a whole number of cells, so a requested extent that is not a
multiple of the cell size is rounded; stems are then drawn inside the *realised*
extent, and `SyntheticStand.area_ha` reports the raster's own area. The painting
radius is also floored just above a cell's half-diagonal, so a crown narrower
than its own cell still reaches that cell's centre instead of falling between
sample points. Together these keep the guarantee the generator exists for: the
returned trees are exactly the stems rendered into the canopy model. At normal
resolutions neither adjustment binds.

The same seed always produces the same stand, byte for byte, which is what
makes the accuracy assertions in the test suite meaningful rather than
decorative. Detection on a well-spaced synthetic stand recovers > 90 % of stems
at > 90 % precision with sub-metre positional error — and the test suite
asserts exactly that.

## References

The methods here are standard forest mensuration and remote sensing; the
implementations are original. For background:

- Brokaw, N. (1982). The definition of treefall gap and its effect on measures
  of forest dynamics. *Biotropica* 14(2).
- Chave, J. et al. (2014). Improved allometric models to estimate the
  aboveground biomass of tropical trees. *Global Change Biology* 20(10).
- Näslund, M. (1937). Skogsförsöksanstaltens gallringsförsök i tallskog.
  *Meddelanden från Statens Skogsförsöksanstalt* 29.
- Popescu, S. & Wynne, R. (2004). Seeing the trees in the forest: using lidar
  and multispectral data fusion with local filtering and variable window size
  for estimating tree height. *Photogrammetric Engineering & Remote Sensing*
  70(5).
- Reineke, L. H. (1933). Perfecting a stand-density index for even-aged
  forests. *Journal of Agricultural Research* 46.
