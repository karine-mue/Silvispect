"""Deterministic synthetic stands.

A synthetic forest with known trees is the only way to test a detector without
shipping proprietary lidar.  Given a seed, :func:`synthesize` produces a
terrain model, a surface model, the derived canopy height model, and the exact
tree list that generated them, so detection accuracy can be measured against
ground truth in the test suite and demonstrated from the command line.
"""

from __future__ import annotations

import math
import random
from collections.abc import Sequence
from dataclasses import dataclass, field

from .allometry import DEFAULT_MODEL, MODELS, SPECIES_DEFAULTS
from .canopy import canopy_height_model
from .grid import Grid
from .inventory import Tree

__all__ = ["SynthConfig", "SyntheticStand", "synthesize"]


@dataclass(frozen=True)
class SynthConfig:
    """Parameters of a synthetic stand.

    Attributes:
        width: East-west extent in metres.
        height: North-south extent in metres.
        cellsize: Raster resolution in metres.
        stems_per_ha: Target stem density before gaps are cut.
        seed: Random seed; the same seed always yields the same stand.
        species: Species codes to draw from.
        species_weights: Relative frequencies, aligned with ``species``.
        dbh_log_mean: Mean of ``log(DBH)`` in centimetres.
        dbh_log_sigma: Standard deviation of ``log(DBH)``.
        dbh_min / dbh_max: Clamp for the drawn diameters.
        height_noise: Relative noise added to the allometric height.
        crown_intercept / crown_slope: Crown radius in metres per tree height.
        surface_noise: Standard deviation of the per-cell measurement noise.
        gap_count: Number of circular canopy openings to cut.
        gap_radius: Range of gap radii in metres.
        terrain_slope: Terrain rise per metre of easting.
        terrain_relief: Amplitude of the sinusoidal terrain undulation.
    """

    width: float = 100.0
    height: float = 100.0
    cellsize: float = 0.5
    stems_per_ha: float = 320.0
    seed: int = 20240917
    species: tuple[str, ...] = ("PIAB", "FASY", "BEPE")
    species_weights: tuple[float, ...] = (0.6, 0.3, 0.1)
    dbh_log_mean: float = math.log(28.0)
    dbh_log_sigma: float = 0.32
    dbh_min: float = 8.0
    dbh_max: float = 85.0
    height_noise: float = 0.04
    crown_intercept: float = 0.9
    crown_slope: float = 0.085
    surface_noise: float = 0.05
    gap_count: int = 2
    gap_radius: tuple[float, float] = (5.0, 9.0)
    terrain_slope: float = 0.03
    terrain_relief: float = 1.5
    xllcorner: float = 0.0
    yllcorner: float = 0.0

    def __post_init__(self) -> None:
        if self.width <= 0 or self.height <= 0:
            raise ValueError("stand extent must be positive")
        if self.cellsize <= 0:
            raise ValueError("cellsize must be positive")
        if len(self.species) != len(self.species_weights):
            raise ValueError("species and species_weights must be the same length")
        if not self.species:
            raise ValueError("at least one species is required")
        if self.stems_per_ha < 0:
            raise ValueError("stems_per_ha must be non-negative")

    @property
    def shape(self) -> tuple[int, int]:
        """Raster dimensions in cells as ``(nrows, ncols)``."""
        return (
            max(1, round(self.height / self.cellsize)),
            max(1, round(self.width / self.cellsize)),
        )

    @property
    def extent(self) -> tuple[float, float]:
        """Realised ``(width, height)`` in metres, rounded to whole cells.

        A raster holds a whole number of cells, so a requested extent that is
        not a multiple of the cell size is rounded.  Stems are drawn inside
        *this* extent rather than the requested one, so that no generated tree
        lands outside the raster describing it.
        """
        nrows, ncols = self.shape
        return (ncols * self.cellsize, nrows * self.cellsize)

    @property
    def area_ha(self) -> float:
        """Ground area of the raster that will be generated, in hectares."""
        width, height = self.extent
        return width * height / 10_000.0


@dataclass
class SyntheticStand:
    """A generated stand together with its ground truth."""

    chm: Grid
    dsm: Grid
    dtm: Grid
    trees: list[Tree]
    gaps: list[tuple[float, float, float]]
    config: SynthConfig = field(repr=False)

    @property
    def area_ha(self) -> float:
        """Ground area of the generated raster, in hectares."""
        return self.chm.area_ha


def synthesize(config: SynthConfig | None = None, **overrides) -> SyntheticStand:
    """Generate a synthetic stand.

    Args:
        config: Base configuration; a default stand is used when omitted.
        **overrides: Field overrides applied on top of ``config``.

    Returns:
        A :class:`SyntheticStand` whose ``trees`` are the stems painted onto the
        canopy surface.  Each one reaches at least its own cell, but the surface
        keeps the *maximum* over overlapping crowns, so a shorter neighbour can
        end up completely hidden by a taller one — heavily so on coarse rasters,
        where many stems share a cell.  That is what a real canopy height model
        does too, and it is why detection recall is bounded by the sensor rather
        than by the detector.
    """
    if config is None:
        config = SynthConfig(**overrides)
    elif overrides:
        config = SynthConfig(**{**config.__dict__, **overrides})

    rng = random.Random(config.seed)
    nrows, ncols = config.shape

    dtm = Grid.filled(
        nrows,
        ncols,
        0.0,
        cellsize=config.cellsize,
        xllcorner=config.xllcorner,
        yllcorner=config.yllcorner,
    )
    for row, col in dtm.cells():
        x, y = dtm.cell_center(row, col)
        local_x = x - config.xllcorner
        local_y = y - config.yllcorner
        elevation = config.terrain_slope * local_x + config.terrain_relief * math.sin(
            local_x / 37.0
        ) * math.cos(local_y / 29.0)
        dtm.set(row, col, elevation)

    gaps = _draw_gaps(rng, config)
    positions = _draw_positions(rng, config, gaps)
    trees = _draw_trees(rng, config, positions)

    canopy = Grid.filled(
        nrows,
        ncols,
        0.0,
        cellsize=config.cellsize,
        xllcorner=config.xllcorner,
        yllcorner=config.yllcorner,
    )
    for tree in trees:
        _render_crown(canopy, tree, config)

    if config.surface_noise > 0:
        for row, col in canopy.cells():
            value = canopy.get(row, col) or 0.0
            canopy.set(row, col, max(0.0, value + rng.gauss(0.0, config.surface_noise)))

    dsm = canopy.combine(dtm, lambda crown, ground: crown + ground)
    chm = canopy_height_model(dsm, dtm)
    return SyntheticStand(chm=chm, dsm=dsm, dtm=dtm, trees=trees, gaps=gaps, config=config)


def _draw_gaps(rng: random.Random, config: SynthConfig) -> list[tuple[float, float, float]]:
    gaps: list[tuple[float, float, float]] = []
    low, high = config.gap_radius
    width, height = config.extent
    for _ in range(max(0, config.gap_count)):
        radius = rng.uniform(low, high)
        x = config.xllcorner + rng.uniform(radius, max(radius, width - radius))
        y = config.yllcorner + rng.uniform(radius, max(radius, height - radius))
        gaps.append((x, y, radius))
    return gaps


def _draw_positions(
    rng: random.Random,
    config: SynthConfig,
    gaps: Sequence[tuple[float, float, float]],
) -> list[tuple[float, float]]:
    target = round(config.stems_per_ha * config.area_ha)
    if target <= 0:
        return []
    spacing = 0.62 * math.sqrt(10_000.0 / max(config.stems_per_ha, 1.0))
    width, height = config.extent
    accepted: list[tuple[float, float]] = []
    attempts = 0
    max_attempts = target * 80
    while len(accepted) < target and attempts < max_attempts:
        attempts += 1
        x = config.xllcorner + rng.uniform(0.0, width)
        y = config.yllcorner + rng.uniform(0.0, height)
        if any(math.hypot(x - gx, y - gy) < gr for gx, gy, gr in gaps):
            continue
        if any(math.hypot(x - ax, y - ay) < spacing for ax, ay in accepted):
            continue
        accepted.append((x, y))
    return accepted


def _draw_trees(
    rng: random.Random, config: SynthConfig, positions: Sequence[tuple[float, float]]
) -> list[Tree]:
    curve = MODELS[DEFAULT_MODEL]
    trees: list[Tree] = []
    for index, (x, y) in enumerate(positions, start=1):
        species = _weighted_choice(rng, config.species, config.species_weights)
        dbh = math.exp(rng.gauss(config.dbh_log_mean, config.dbh_log_sigma))
        dbh = min(config.dbh_max, max(config.dbh_min, dbh))
        params = SPECIES_DEFAULTS.get(species.upper(), SPECIES_DEFAULTS[""])
        height = curve.evaluate(dbh, params) * (1.0 + rng.gauss(0.0, config.height_noise))
        height = max(2.5, height)
        crown_radius = (config.crown_intercept + config.crown_slope * height) * rng.uniform(
            0.88, 1.12
        )
        trees.append(
            Tree(
                tree_id=f"T{index:04d}",
                x=x,
                y=y,
                species=species,
                dbh_cm=round(dbh, 1),
                height_m=round(height, 2),
                crown_diameter_m=round(2.0 * crown_radius, 2),
                status="live",
            )
        )
    return trees


def _render_crown(canopy: Grid, tree: Tree, config: SynthConfig) -> None:
    """Paint a paraboloid crown onto the canopy surface, keeping the maximum."""
    radius = (tree.crown_diameter_m or 2.0) / 2.0
    height = tree.height_m or 0.0
    if radius <= 0 or height <= 0:
        return
    # A crown narrower than the cell it stands in can fall entirely between cell
    # centres and paint nothing, which would leave a ground-truth stem with no
    # trace in the raster it is supposed to describe.  Floor the painting radius
    # just above the half-diagonal of a cell so every tree reaches its own cell
    # centre.  At the resolutions this generator is normally used at, every
    # crown is already wider than this and the floor never binds.
    radius = max(radius, 0.75 * canopy.cellsize)
    cell = canopy.cell_of(tree.x, tree.y)
    if cell is None:
        return
    reach = math.ceil(radius / canopy.cellsize) + 1
    row0, col0 = cell
    for row, col in canopy.window(row0, col0, reach):
        x, y = canopy.cell_center(row, col)
        distance = math.hypot(x - tree.x, y - tree.y)
        if distance > radius:
            continue
        ratio = distance / radius
        value = height * (1.0 - ratio * ratio)
        current = canopy.get(row, col) or 0.0
        if value > current:
            canopy.set(row, col, value)


def _weighted_choice(rng: random.Random, options: Sequence[str], weights: Sequence[float]) -> str:
    total = sum(weights)
    if total <= 0:
        return options[0]
    threshold = rng.uniform(0.0, total)
    cumulative = 0.0
    for option, weight in zip(options, weights, strict=True):
        cumulative += weight
        if threshold <= cumulative:
            return option
    return options[-1]
