"""Canopy height model derivation and canopy-level descriptors.

A canopy height model (CHM) is the height of vegetation above the ground,
normally obtained by subtracting a digital terrain model (DTM) from a digital
surface model (DSM).  Everything downstream in Silvispect — tree detection,
gap analysis, structural metrics — reads a CHM.
"""

from __future__ import annotations

import heapq
import math
from collections import deque
from collections.abc import Iterable
from dataclasses import dataclass

from .grid import Grid, GridError

__all__ = [
    "STRATA_BREAKS",
    "CanopyGap",
    "canopy_cover",
    "canopy_height_model",
    "distance_to_canopy",
    "find_gaps",
    "gap_fraction",
    "rugosity",
    "vertical_strata",
]

#: Default vertical strata breaks in metres: ground, regeneration, understorey,
#: mid-storey, canopy, emergent.
STRATA_BREAKS: tuple[float, ...] = (0.5, 2.0, 5.0, 15.0, 25.0)

STRATA_NAMES: tuple[str, ...] = (
    "ground",
    "regeneration",
    "understorey",
    "midstorey",
    "canopy",
    "emergent",
)


def canopy_height_model(dsm: Grid, dtm: Grid, *, floor: float = 0.0) -> Grid:
    """Return ``dsm - dtm`` clamped at ``floor``.

    Negative heights are a normalisation artefact rather than a physical
    quantity, so they are clamped instead of being discarded: a cell that
    measures slightly below ground is still a measured, treeless cell.
    """
    chm = dsm.combine(dtm, lambda surface, terrain: surface - terrain)
    return chm.clip(minimum=floor)


def canopy_cover(chm: Grid, threshold: float = 2.0) -> float:
    """Fraction of valid cells whose height is at or above ``threshold``.

    Returns ``0.0`` for an entirely empty grid.
    """
    total = 0
    covered = 0
    for _, _, height in chm.valid_cells():
        total += 1
        if height >= threshold:
            covered += 1
    if total == 0:
        return 0.0
    return covered / total


@dataclass(frozen=True)
class CanopyGap:
    """A connected patch of canopy below the cover threshold."""

    gap_id: int
    cells: tuple[tuple[int, int], ...]
    area: float
    centroid_x: float
    centroid_y: float
    mean_height: float
    max_height: float
    touches_edge: bool
    inscribed_radius: float = 0.0

    @property
    def area_ha(self) -> float:
        return self.area / 10_000.0

    @property
    def equivalent_radius(self) -> float:
        """Radius of a circle with the same area, in map units."""
        return math.sqrt(self.area / math.pi)

    @property
    def width(self) -> float:
        """Diameter of the largest circle that fits inside the gap."""
        return 2.0 * self.inscribed_radius

    def as_dict(self) -> dict[str, float | int | bool]:
        return {
            "gap_id": self.gap_id,
            "cell_count": len(self.cells),
            "area_m2": round(self.area, 2),
            "equivalent_radius_m": round(self.equivalent_radius, 2),
            "width_m": round(self.width, 2),
            "centroid_x": round(self.centroid_x, 2),
            "centroid_y": round(self.centroid_y, 2),
            "mean_height_m": round(self.mean_height, 2),
            "max_height_m": round(self.max_height, 2),
            "touches_edge": self.touches_edge,
        }


def _chamfer(chm: Grid, seeds: Iterable[int], *, ceiling: float | None = None) -> list[float]:
    """Distance in map units from every cell to the nearest seed cell.

    Orthogonal steps cost one cell width and diagonal steps ``sqrt(2)`` cell
    widths, which approximates the Euclidean distance closely enough for
    morphology on a canopy raster.  ``ceiling`` stops the search early.
    """
    distances: list[float] = [math.inf] * len(chm.values)
    frontier: list[tuple[float, int]] = []
    for position in seeds:
        distances[position] = 0.0
        frontier.append((0.0, position))
    heapq.heapify(frontier)

    straight = chm.cellsize
    diagonal = chm.cellsize * math.sqrt(2.0)
    while frontier:
        distance, position = heapq.heappop(frontier)
        if distance > distances[position]:
            continue
        if ceiling is not None and distance > ceiling:
            continue
        row, col = divmod(position, chm.ncols)
        for nrow, ncol in chm.neighbors(row, col, connectivity=8):
            step = straight if (nrow == row or ncol == col) else diagonal
            candidate = distance + step
            index = nrow * chm.ncols + ncol
            if candidate < distances[index]:
                distances[index] = candidate
                heapq.heappush(frontier, (candidate, index))
    return distances


def _distance_to_outside(chm: Grid) -> list[float]:
    """Distance from each cell centre to just beyond the nearest grid edge."""
    half = chm.cellsize / 2.0
    return [
        min(col, chm.ncols - 1 - col, row, chm.nrows - 1 - row) * chm.cellsize + half
        for row in range(chm.nrows)
        for col in range(chm.ncols)
    ]


def _canopy_distance_field(chm: Grid, threshold: float) -> list[float]:
    """Distance from every cell to the nearest canopy cell or plot edge.

    **The plot edge counts as a boundary.**  Measuring only to canopy would let
    an opening that runs off the side of the raster claim an inscribed circle
    larger than the raster itself — a 10 m x 10 m plot with one corner tree
    reported a 25 m gap width, which no circle inside it could have.  Nothing is
    known about the forest beyond the extent, so the honest bound is the extent:
    the reported width describes the opening *as mapped*, and ``touches_edge``
    marks it as possibly continuing further.

    Taking the minimum against the edge distance also covers a raster with no
    canopy at all, where the chamfer search alone would leave every cell at
    infinity — and JSON has no way to write that down.
    """
    seeds = [
        position for position, value in enumerate(chm.values) if value is None or value >= threshold
    ]
    to_canopy = _chamfer(chm, seeds)
    to_outside = _distance_to_outside(chm)
    return [min(canopy, outside) for canopy, outside in zip(to_canopy, to_outside, strict=True)]


def distance_to_canopy(chm: Grid, *, threshold: float = 2.0) -> Grid:
    """Distance from every cell to the nearest canopy cell or plot edge.

    Cells at or above ``threshold`` are the sources and get distance zero.
    Cells with no data are treated as unknown rather than as openings, so they
    are sources too.  The plot edge bounds the result — see
    :func:`_canopy_distance_field` for why.
    """
    out = chm.like()
    for position, distance in enumerate(_canopy_distance_field(chm, threshold)):
        out.values[position] = None if math.isinf(distance) else distance
    return out


def _opened_gap_mask(chm: Grid, threshold: float, radius: float) -> tuple[list[bool], list[float]]:
    """Morphologically open the sub-canopy mask with a disc of ``radius``.

    The opening is the union of every disc of that radius that fits entirely
    inside the sub-canopy area, computed as an erosion (cells further than
    ``radius`` from any canopy cell) followed by a dilation of the same size.
    It removes the narrow web of ground visible between neighbouring crowns
    without eroding genuine openings.

    Returns:
        The opened mask and the distance-to-canopy values behind it.
    """
    to_canopy = _canopy_distance_field(chm, threshold)
    eroded = [position for position, distance in enumerate(to_canopy) if distance > radius]
    mask = [False] * len(chm.values)
    if not eroded:
        return mask, to_canopy
    from_core = _chamfer(chm, eroded, ceiling=radius)
    for position, value in enumerate(chm.values):
        if value is None or value >= threshold:
            continue
        if from_core[position] <= radius + 1e-9:
            mask[position] = True
    return mask, to_canopy


def find_gaps(
    chm: Grid,
    *,
    threshold: float = 2.0,
    min_area: float = 25.0,
    min_width: float = 0.0,
    connectivity: int = 8,
    include_edge_gaps: bool = True,
) -> list[CanopyGap]:
    """Label canopy openings as connected components below ``threshold``.

    A pure threshold rule is not enough on a fine raster: the ground visible
    between neighbouring crowns forms a connected web that would be reported as
    one enormous opening.  ``min_width`` applies the Brokaw criterion by
    morphologically opening the sub-canopy mask with a disc of that diameter,
    which trims the inter-crown web — including the parts of it attached to a
    genuine opening — before the components are labelled.

    Args:
        chm: Canopy height model.
        threshold: Height at or above which a cell counts as canopy.
        min_area: Discard gaps smaller than this area in square map units.
        min_width: Narrowest opening to report, in map units.  ``0`` disables
            the morphological filter and reports the raw threshold components.
        connectivity: ``4`` or ``8`` cell connectivity.
        include_edge_gaps: Keep gaps that touch the grid border.  Such gaps are
            truncated by the extent, so their area is a lower bound.

    Returns:
        Gaps sorted by descending area.
    """
    if min_area < 0:
        raise GridError("min_area must be non-negative")
    if min_width < 0:
        raise GridError("min_width must be non-negative")

    if min_width > 0:
        mask, to_canopy = _opened_gap_mask(chm, threshold, min_width / 2.0)
    else:
        mask = [value is not None and value < threshold for value in chm.values]
        to_canopy = _canopy_distance_field(chm, threshold)

    # Morphological opening insets the mask from the border, so a gap that
    # genuinely runs off the plot no longer has a cell *on* the border.  Judge
    # edge contact by proximity to the outside instead, using the same radius
    # the opening ate away; with no opening this reduces to "is a border cell".
    to_outside = _distance_to_outside(chm)
    edge_reach = (min_width / 2.0) + chm.cellsize / 2.0

    seen = [False] * len(chm.values)
    gaps: list[CanopyGap] = []
    next_id = 1
    for row, col, _ in chm.valid_cells():
        start = row * chm.ncols + col
        if seen[start] or not mask[start]:
            continue
        queue: deque[tuple[int, int]] = deque([(row, col)])
        seen[start] = True
        members: list[tuple[int, int]] = []
        heights: list[float] = []
        while queue:
            r, c = queue.popleft()
            members.append((r, c))
            value = chm.values[r * chm.ncols + c]
            if value is not None:  # the mask only covers valid cells
                heights.append(value)
            for nr, nc in chm.neighbors(r, c, connectivity):
                idx = nr * chm.ncols + nc
                if seen[idx] or not mask[idx]:
                    continue
                seen[idx] = True
                queue.append((nr, nc))
        area = len(members) * chm.cell_area
        if area < min_area:
            continue
        touches_edge = min(to_outside[r * chm.ncols + c] for r, c in members) <= edge_reach
        if touches_edge and not include_edge_gaps:
            continue
        inscribed = max(to_canopy[r * chm.ncols + c] for r, c in members)
        xs, ys = zip(*(chm.cell_center(r, c) for r, c in members), strict=True)
        gaps.append(
            CanopyGap(
                gap_id=next_id,
                cells=tuple(members),
                area=area,
                centroid_x=sum(xs) / len(xs),
                centroid_y=sum(ys) / len(ys),
                mean_height=sum(heights) / len(heights),
                max_height=max(heights),
                touches_edge=touches_edge,
                inscribed_radius=inscribed,
            )
        )
        next_id += 1
    gaps.sort(key=lambda gap: gap.area, reverse=True)
    return [
        CanopyGap(
            gap_id=position,
            cells=gap.cells,
            area=gap.area,
            centroid_x=gap.centroid_x,
            centroid_y=gap.centroid_y,
            mean_height=gap.mean_height,
            max_height=gap.max_height,
            touches_edge=gap.touches_edge,
            inscribed_radius=gap.inscribed_radius,
        )
        for position, gap in enumerate(gaps, start=1)
    ]


def gap_fraction(chm: Grid, *, threshold: float = 2.0) -> float:
    """Fraction of the canopy below ``threshold`` — the complement of cover."""
    return 1.0 - canopy_cover(chm, threshold)


def rugosity(chm: Grid, *, threshold: float = 0.0) -> float:
    """Standard deviation of canopy heights above ``threshold``.

    Canopy rugosity is a compact proxy for vertical structural complexity:
    even-aged plantations score low, multi-layered stands score high.
    """
    heights = [h for _, _, h in chm.valid_cells() if h > threshold]
    if len(heights) < 2:
        return 0.0
    mean = sum(heights) / len(heights)
    return math.sqrt(sum((h - mean) ** 2 for h in heights) / (len(heights) - 1))


def vertical_strata(chm: Grid, breaks: tuple[float, ...] = STRATA_BREAKS) -> dict[str, float]:
    """Return the proportion of canopy cells falling in each vertical stratum."""
    if list(breaks) != sorted(breaks):
        raise GridError("strata breaks must be ascending")
    names = _strata_names(len(breaks) + 1)
    counts = dict.fromkeys(names, 0)
    total = 0
    for _, _, height in chm.valid_cells():
        total += 1
        slot = len(breaks)
        for i, edge in enumerate(breaks):
            if height < edge:
                slot = i
                break
        counts[names[slot]] += 1
    if total == 0:
        return dict.fromkeys(names, 0.0)
    return {name: counts[name] / total for name in names}


def _strata_names(count: int) -> tuple[str, ...]:
    if count == len(STRATA_NAMES):
        return STRATA_NAMES
    return tuple(f"stratum_{i}" for i in range(count))
