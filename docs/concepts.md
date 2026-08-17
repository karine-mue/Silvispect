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
  gap 1        511 m2  width  12.3 m  ...   <- one opening fused to the web
```

Silvispect applies the **Brokaw criterion**: a gap must be wide enough to
contain a circle of a minimum diameter (5 m by default). It is implemented as a
morphological opening of the sub-canopy mask — an erosion by the disc radius
followed by a dilation of the same size — computed from a chamfer distance
transform. The opening removes the web *including the tentacles attached to a
genuine opening*, which a simple "discard small components" filter cannot do.
Labelling then runs on the opened mask:

```console
$ silvispect gaps data/plot-a1-chm.asc
4 gaps at or above 25 m2 and 5.0 m wide; 36.4% of the plot is below 2.0 m
  gap 1        164 m2  width  12.3 m  ...   <- the opening alone
```

The 511 m² blob and the 164 m² opening have the *same* inscribed width: the
extra 347 m² is entirely tentacles of inter-crown ground, which is why a
"discard narrow components" filter cannot remove it and an opening can.

Each gap reports its area, its **width** (the diameter of the largest inscribed
circle), its equivalent radius, its centroid, and whether it touches the plot
edge — an edge gap is truncated by the extent, so its area is a lower bound
and Silvispect downgrades its severity accordingly.

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

with defaults `0.8 + 0.055·h`. A cell is a treetop if no cell within `r(h)` is
taller. Plateaus of exactly equal height resolve to a single apex — the first
in row-major order — so results are deterministic.

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

`match_trees` pairs crowns with field stems greedily: all pairs within a
tolerance radius are sorted by distance and accepted while both partners are
still unpaired. This yields a stable one-to-one assignment without a full
Hungarian solve, which the problem size does not justify.

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
to tighten the simplex. No derivatives, no SciPy.

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
