"""The inspection engine: turn a stand into a list of coded findings.

Silvispect's premise is that a forest dataset should be *inspectable* the way a
codebase is: run a fixed set of rules over it, get back stable codes with
severities, and fail a pipeline when something crosses a threshold.  Rules are
plain functions registered in :data:`RULES`; each one receives an
:class:`InspectionContext` and yields zero or more :class:`Finding` objects.
"""

from __future__ import annotations

import json
import math
from collections.abc import Callable, Iterable, Sequence
from dataclasses import asdict, dataclass, field, fields
from enum import IntEnum
from pathlib import Path

from .allometry import (
    DEFAULT_MODEL,
    FittedModel,
    default_model,
    fit_by_species,
    residual_scores,
)
from .canopy import CanopyGap, canopy_cover, find_gaps, gap_fraction, rugosity, vertical_strata
from .detect import Crown, DetectionConfig, DetectionResult, detect_trees
from .grid import Grid
from .inventory import MatchResult, Tree, match_trees, trees_from_crowns
from .metrics import StandMetrics, stand_metrics

__all__ = [
    "RULES",
    "Finding",
    "InspectionConfig",
    "InspectionContext",
    "InspectionReport",
    "Severity",
    "inspect_stand",
    "rule",
]


class Severity(IntEnum):
    """How much attention a finding deserves."""

    INFO = 0
    NOTICE = 1
    WARNING = 2
    CRITICAL = 3

    @property
    def label(self) -> str:
        return self.name.lower()

    @classmethod
    def parse(cls, text: str) -> Severity:
        try:
            return cls[text.strip().upper()]
        except KeyError as exc:
            raise ValueError(
                f"unknown severity {text!r}; choose from {[s.label for s in cls]}"
            ) from exc


@dataclass(frozen=True)
class Finding:
    """One rule outcome."""

    code: str
    severity: Severity
    title: str
    detail: str
    subject: str = "stand"
    value: float | None = None
    threshold: float | None = None

    def as_dict(self) -> dict[str, object]:
        return {
            "code": self.code,
            "severity": self.severity.label,
            "title": self.title,
            "detail": self.detail,
            "subject": self.subject,
            "value": _round(self.value),
            "threshold": _round(self.threshold),
        }


@dataclass
class InspectionConfig:
    """Thresholds driving the rules.

    Every field is a documented, editable number: an inspection is only as
    meaningful as the silvicultural expectations behind it, and those differ by
    region, species and management goal.
    """

    canopy_threshold: float = 2.0
    min_stems_per_ha: float = 200.0
    max_stems_per_ha: float = 1400.0
    min_sdi: float = 150.0
    max_sdi: float = 900.0
    min_canopy_cover: float = 0.60
    max_gap_area: float = 400.0
    max_gap_fraction: float = 0.30
    min_gap_area: float = 25.0
    min_gap_width: float = 5.0
    min_gini: float = 0.15
    min_shannon: float = 0.35
    max_species_share: float = 0.90
    height_residual_z: float = 2.5
    dbh_outlier_z: float = 3.0
    max_plausible_dbh_cm: float = 250.0
    max_plausible_height_m: float = 120.0
    min_recall: float = 0.70
    min_precision: float = 0.70
    max_height_bias_m: float = 1.50
    match_tolerance_m: float = 2.50
    min_allometry_points: int = 8

    #: Fields constrained to a proportion in [0, 1].
    _FRACTIONS = (
        "min_canopy_cover",
        "max_gap_fraction",
        "min_recall",
        "min_precision",
        "max_species_share",
    )

    def __post_init__(self) -> None:
        """Reject configurations that would make a rule meaningless or crash.

        A profile is user input, so it is validated at the boundary rather than
        being allowed to surface later as a traceback from deep inside a fit.
        """
        if isinstance(self.min_allometry_points, bool) or (
            self.min_allometry_points != int(self.min_allometry_points)
        ):
            raise ValueError(
                f"min_allometry_points must be a whole number, got {self.min_allometry_points!r}"
            )
        if self.min_allometry_points < 1:
            raise ValueError("min_allometry_points must be at least 1")
        for name, value in asdict(self).items():
            # Infinity survives every range check below and only fails later, at
            # JSON serialisation or inside a fit, so it is rejected here.
            if isinstance(value, float) and not math.isfinite(value):
                raise ValueError(f"{name} must be a finite number, got {value}")
        for name in self._FRACTIONS:
            value = getattr(self, name)
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must be a proportion between 0 and 1, got {value}")
        for name, value in asdict(self).items():
            if isinstance(value, (int, float)) and value < 0:
                raise ValueError(f"{name} must not be negative, got {value}")
        for low, high in (
            ("min_stems_per_ha", "max_stems_per_ha"),
            ("min_sdi", "max_sdi"),
        ):
            if getattr(self, low) > getattr(self, high):
                raise ValueError(f"{low} must not exceed {high}")

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> InspectionConfig:
        if not isinstance(data, dict):
            raise ValueError(f"configuration must be a JSON object, got {type(data).__name__}")
        known = {f.name for f in fields(cls)}
        unknown = set(data) - known
        if unknown:
            raise ValueError(f"unknown configuration keys: {sorted(unknown)}")
        converted: dict[str, float | int] = {}
        for key, value in data.items():
            if isinstance(value, bool) or not isinstance(value, (int, float, str)):
                raise ValueError(f"{key} must be a number, got {type(value).__name__}")
            number = float(value)
            if not math.isfinite(number):
                raise ValueError(f"{key} must be a finite number, got {value}")
            if key == "min_allometry_points":
                # Truncating 1.9 to 1 would silently apply a threshold the
                # profile never asked for.
                if number != int(number):
                    raise ValueError(f"{key} must be a whole number, got {value}")
                converted[key] = int(number)
            else:
                converted[key] = number
        return cls(**converted)  # type: ignore[arg-type]

    @classmethod
    def load(cls, path: str | Path) -> InspectionConfig:
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))

    def as_dict(self) -> dict[str, float]:
        return asdict(self)


@dataclass
class InspectionContext:
    """Everything the rules are allowed to look at."""

    config: InspectionConfig
    area_ha: float
    trees: list[Tree]
    metrics: StandMetrics
    metrics_source: str
    field_trees: list[Tree] | None = None
    chm: Grid | None = None
    detection: DetectionResult | None = None
    gaps: list[CanopyGap] = field(default_factory=list)
    match: MatchResult | None = None
    model: FittedModel | None = None
    species_models: dict[str, FittedModel] = field(default_factory=dict)
    fitted: bool = False

    @property
    def crowns(self) -> list[Crown]:
        return [] if self.detection is None else self.detection.crowns

    @property
    def has_field_data(self) -> bool:
        """Whether a field inventory was supplied, even an empty one."""
        return self.field_trees is not None

    @property
    def has_stem_source(self) -> bool:
        """Whether any process actually looked for stems.

        False when detection was skipped and no inventory was supplied: the
        stand then has *unknown* stocking, not zero stocking, so the stem-based
        rules must stay quiet instead of reporting an empty forest.
        """
        return self.detection is not None or self.field_trees is not None

    @property
    def has_raster(self) -> bool:
        return self.chm is not None


Rule = Callable[[InspectionContext], Iterable[Finding]]

#: Registered rules, in registration order.
RULES: list[Rule] = []


def rule(fn: Rule) -> Rule:
    """Register a rule function."""
    RULES.append(fn)
    return fn


# ----------------------------------------------------------------------
# stocking and density
# ----------------------------------------------------------------------
@rule
def check_stocking(ctx: InspectionContext) -> Iterable[Finding]:
    """SV001/SV002 — stem density outside the acceptable band."""
    if not ctx.has_stem_source:
        return
    density = ctx.metrics.stems_per_ha
    if density is None:  # nothing counted stems, so there is nothing to judge
        return
    cfg = ctx.config
    if density < cfg.min_stems_per_ha:
        yield Finding(
            "SV001",
            Severity.WARNING,
            "Understocked stand",
            f"{density:.0f} stems/ha is below the minimum of "
            f"{cfg.min_stems_per_ha:.0f} stems/ha ({ctx.metrics_source} data).",
            value=density,
            threshold=cfg.min_stems_per_ha,
        )
    elif density > cfg.max_stems_per_ha:
        yield Finding(
            "SV002",
            Severity.NOTICE,
            "Overstocked stand",
            f"{density:.0f} stems/ha exceeds the maximum of "
            f"{cfg.max_stems_per_ha:.0f} stems/ha; competition-driven mortality "
            "is likely.",
            value=density,
            threshold=cfg.max_stems_per_ha,
        )


@rule
def check_density_index(ctx: InspectionContext) -> Iterable[Finding]:
    """SV003/SV004 — Reineke stand density index outside the managed range."""
    if not ctx.has_stem_source:
        return
    sdi = ctx.metrics.sdi
    if sdi is None:
        return
    cfg = ctx.config
    if sdi < cfg.min_sdi:
        yield Finding(
            "SV003",
            Severity.NOTICE,
            "Low stand density index",
            f"SDI {sdi:.0f} is below {cfg.min_sdi:.0f}; the site is not fully "
            "occupying the growing space.",
            value=sdi,
            threshold=cfg.min_sdi,
        )
    elif sdi > cfg.max_sdi:
        yield Finding(
            "SV004",
            Severity.WARNING,
            "High stand density index",
            f"SDI {sdi:.0f} exceeds {cfg.max_sdi:.0f}; the stand is approaching "
            "the zone of imminent competition mortality.",
            value=sdi,
            threshold=cfg.max_sdi,
        )


# ----------------------------------------------------------------------
# canopy condition
# ----------------------------------------------------------------------
@rule
def check_canopy_cover(ctx: InspectionContext) -> Iterable[Finding]:
    """SV010 — canopy cover below the silvicultural target."""
    if ctx.chm is None:
        return
    cover = canopy_cover(ctx.chm, ctx.config.canopy_threshold)
    if cover < ctx.config.min_canopy_cover:
        yield Finding(
            "SV010",
            Severity.WARNING,
            "Canopy cover below target",
            f"Cover of {cover:.1%} at the {ctx.config.canopy_threshold:.1f} m "
            f"threshold is under the target of {ctx.config.min_canopy_cover:.0%}.",
            value=cover,
            threshold=ctx.config.min_canopy_cover,
        )


@rule
def check_large_gaps(ctx: InspectionContext) -> Iterable[Finding]:
    """SV011 — individual canopy openings larger than the allowed size."""
    for gap in ctx.gaps:
        if gap.area <= ctx.config.max_gap_area:
            continue
        yield Finding(
            "SV011",
            Severity.WARNING if gap.touches_edge else Severity.CRITICAL,
            "Large canopy gap",
            f"Gap {gap.gap_id} covers {gap.area:.0f} m2 "
            f"(width {gap.width:.1f} m) centred at "
            f"({gap.centroid_x:.1f}, {gap.centroid_y:.1f})"
            + (", truncated by the plot edge" if gap.touches_edge else "")
            + ".",
            subject=f"gap:{gap.gap_id}",
            value=gap.area,
            threshold=ctx.config.max_gap_area,
        )


@rule
def check_gap_fraction(ctx: InspectionContext) -> Iterable[Finding]:
    """SV012 — too much of the plot is open canopy."""
    if ctx.chm is None:
        return
    fraction = gap_fraction(ctx.chm, threshold=ctx.config.canopy_threshold)
    if fraction > ctx.config.max_gap_fraction:
        yield Finding(
            "SV012",
            Severity.NOTICE,
            "High open-canopy fraction",
            f"{fraction:.1%} of the plot is below "
            f"{ctx.config.canopy_threshold:.1f} m, above the "
            f"{ctx.config.max_gap_fraction:.0%} limit.",
            value=fraction,
            threshold=ctx.config.max_gap_fraction,
        )


# ----------------------------------------------------------------------
# structure and composition
# ----------------------------------------------------------------------
@rule
def check_structural_diversity(ctx: InspectionContext) -> Iterable[Finding]:
    """SV020 — size distribution is nearly uniform."""
    if not ctx.has_stem_source:
        return
    gini = ctx.metrics.gini_basal_area
    if gini is None:
        return
    if gini < ctx.config.min_gini:
        yield Finding(
            "SV020",
            Severity.NOTICE,
            "Low structural diversity",
            f"Gini coefficient of basal area is {gini:.2f}, below "
            f"{ctx.config.min_gini:.2f}: the stand is structurally uniform and "
            "vulnerable to a single disturbance agent.",
            value=gini,
            threshold=ctx.config.min_gini,
        )


@rule
def check_species_diversity(ctx: InspectionContext) -> Iterable[Finding]:
    """SV021/SV022 — species mixture too poor or too dominated."""
    if not ctx.has_field_data:
        return
    shannon = ctx.metrics.shannon
    if shannon is not None and shannon < ctx.config.min_shannon:
        yield Finding(
            "SV021",
            Severity.NOTICE,
            "Low species diversity",
            f"Shannon index {shannon:.2f} is below {ctx.config.min_shannon:.2f}.",
            value=shannon,
            threshold=ctx.config.min_shannon,
        )
    shares = ctx.metrics.species_shares
    if shares:
        species, share = next(iter(shares.items()))
        if share > ctx.config.max_species_share:
            yield Finding(
                "SV022",
                Severity.NOTICE,
                "Single-species dominance",
                f"{species} accounts for {share:.0%} of stems, above the "
                f"{ctx.config.max_species_share:.0%} limit.",
                subject=f"species:{species}",
                value=share,
                threshold=ctx.config.max_species_share,
            )


# ----------------------------------------------------------------------
# per-tree data quality
# ----------------------------------------------------------------------
@rule
def check_height_anomalies(ctx: InspectionContext) -> Iterable[Finding]:
    """SV030 — stems far from the stand's own height-diameter curve."""
    if ctx.model is None or not ctx.has_field_data:
        return
    groups: dict[str, list[Tree]] = {}
    for tree in ctx.field_trees or []:
        code = tree.species.strip().upper()
        key = code if code in ctx.species_models else ""
        groups.setdefault(key, []).append(tree)

    for key, members in sorted(groups.items()):
        model = ctx.species_models.get(key, ctx.model)
        scope = f"the {key} " if key else "the pooled "
        for tree, residual, z in residual_scores(model, members):
            if abs(z) < ctx.config.height_residual_z:
                continue
            direction = "taller" if residual > 0 else "shorter"
            yield Finding(
                "SV030",
                Severity.WARNING if abs(z) >= ctx.config.height_residual_z + 1 else Severity.NOTICE,
                "Height-diameter outlier",
                f"Tree {tree.tree_id} is {abs(residual):.1f} m {direction} than "
                f"{scope}{model.name} curve predicts (z = {z:+.1f}); check for "
                "crown breakage or a transcription error.",
                subject=f"tree:{tree.tree_id}",
                value=abs(z),
                threshold=ctx.config.height_residual_z,
            )


@rule
def check_dbh_outliers(ctx: InspectionContext) -> Iterable[Finding]:
    """SV031 — diameters far outside the plot's own distribution."""
    if not ctx.has_field_data:
        return
    diameters = [tree.dbh_cm for tree in ctx.field_trees or [] if tree.dbh_cm and tree.dbh_cm > 0]
    if len(diameters) < 5:
        return
    mean = sum(diameters) / len(diameters)
    variance = sum((d - mean) ** 2 for d in diameters) / (len(diameters) - 1)
    stdev = math.sqrt(variance)
    if stdev == 0:
        return
    for tree in ctx.field_trees or []:
        if not tree.dbh_cm or tree.dbh_cm <= 0:
            continue
        z = (tree.dbh_cm - mean) / stdev
        if abs(z) < ctx.config.dbh_outlier_z:
            continue
        yield Finding(
            "SV031",
            Severity.NOTICE,
            "Diameter outlier",
            f"Tree {tree.tree_id} has DBH {tree.dbh_cm:.1f} cm, {abs(z):.1f} "
            f"standard deviations from the plot mean of {mean:.1f} cm.",
            subject=f"tree:{tree.tree_id}",
            value=abs(z),
            threshold=ctx.config.dbh_outlier_z,
        )


@rule
def check_missing_measurements(ctx: InspectionContext) -> Iterable[Finding]:
    """SV040 — live stems lacking a diameter or a height."""
    if not ctx.has_field_data:
        return
    for tree in ctx.field_trees or []:
        if not tree.is_live:
            continue
        missing = [
            name
            for name, value in (("dbh_cm", tree.dbh_cm), ("height_m", tree.height_m))
            if value is None
        ]
        if not missing:
            continue
        yield Finding(
            "SV040",
            Severity.NOTICE,
            "Missing measurement",
            f"Tree {tree.tree_id} has no {' and no '.join(missing)}.",
            subject=f"tree:{tree.tree_id}",
        )


@rule
def check_duplicate_ids(ctx: InspectionContext) -> Iterable[Finding]:
    """SV041 — the same tree identifier used more than once."""
    if not ctx.has_field_data:
        return
    seen: dict[str, int] = {}
    for tree in ctx.field_trees or []:
        seen[tree.tree_id] = seen.get(tree.tree_id, 0) + 1
    for tree_id, count in sorted(seen.items()):
        if count > 1:
            yield Finding(
                "SV041",
                Severity.CRITICAL,
                "Duplicate tree identifier",
                f"Identifier {tree_id!r} appears {count} times; records cannot be "
                "linked across re-measurements.",
                subject=f"tree:{tree_id}",
                value=float(count),
                threshold=1.0,
            )


@rule
def check_implausible_dimensions(ctx: InspectionContext) -> Iterable[Finding]:
    """SV042 — non-positive or physically implausible dimensions."""
    if not ctx.has_field_data:
        return
    cfg = ctx.config
    for tree in ctx.field_trees or []:
        problems: list[str] = []
        if tree.dbh_cm is not None and tree.dbh_cm <= 0:
            problems.append(f"DBH {tree.dbh_cm:.1f} cm is not positive")
        if tree.height_m is not None and tree.height_m <= 0:
            problems.append(f"height {tree.height_m:.1f} m is not positive")
        if tree.dbh_cm is not None and tree.dbh_cm > cfg.max_plausible_dbh_cm:
            problems.append(f"DBH {tree.dbh_cm:.1f} cm exceeds {cfg.max_plausible_dbh_cm:.0f} cm")
        if tree.height_m is not None and tree.height_m > cfg.max_plausible_height_m:
            problems.append(
                f"height {tree.height_m:.1f} m exceeds {cfg.max_plausible_height_m:.0f} m"
            )
        if problems:
            yield Finding(
                "SV042",
                Severity.CRITICAL,
                "Implausible dimension",
                f"Tree {tree.tree_id}: " + "; ".join(problems) + ".",
                subject=f"tree:{tree.tree_id}",
            )


@rule
def check_coordinates_in_extent(ctx: InspectionContext) -> Iterable[Finding]:
    """SV043 — field stems that fall outside the raster extent."""
    if ctx.chm is None or not ctx.has_field_data:
        return
    xmin, ymin, xmax, ymax = ctx.chm.extent
    for tree in ctx.field_trees or []:
        if xmin <= tree.x <= xmax and ymin <= tree.y <= ymax:
            continue
        yield Finding(
            "SV043",
            Severity.WARNING,
            "Coordinate outside raster extent",
            f"Tree {tree.tree_id} at ({tree.x:.1f}, {tree.y:.1f}) lies outside "
            f"the canopy model extent ({xmin:.0f}, {ymin:.0f}) - "
            f"({xmax:.0f}, {ymax:.0f}).",
            subject=f"tree:{tree.tree_id}",
        )


# ----------------------------------------------------------------------
# remote sensing agreement
# ----------------------------------------------------------------------
@rule
def check_detection_agreement(ctx: InspectionContext) -> Iterable[Finding]:
    """SV050/SV051/SV052 — detection quality against the field reference."""
    match = ctx.match
    if match is None:
        return
    cfg = ctx.config
    if match.recall < cfg.min_recall:
        yield Finding(
            "SV050",
            Severity.WARNING,
            "Low detection recall",
            f"Only {match.recall:.0%} of field stems were detected "
            f"({match.matched} of {match.matched + len(match.omissions)}); "
            "suppressed and understorey trees are invisible to a canopy model.",
            subject="detection",
            value=match.recall,
            threshold=cfg.min_recall,
        )
    if match.precision < cfg.min_precision:
        yield Finding(
            "SV051",
            Severity.WARNING,
            "Low detection precision",
            f"{len(match.commissions)} of {match.matched + len(match.commissions)} "
            "detections have no field counterpart; the crown segmentation may be "
            "splitting single crowns.",
            subject="detection",
            value=match.precision,
            threshold=cfg.min_precision,
        )
    bias = match.height_bias
    if bias is not None and abs(bias) > cfg.max_height_bias_m:
        direction = "over" if bias > 0 else "under"
        yield Finding(
            "SV052",
            Severity.WARNING,
            "Systematic height bias",
            f"Detected heights {direction}estimate field heights by "
            f"{abs(bias):.2f} m on average (RMSE "
            f"{match.height_rmse:.2f} m); check the terrain model.",
            subject="detection",
            value=abs(bias),
            threshold=cfg.max_height_bias_m,
        )


@dataclass
class InspectionReport:
    """The result of an inspection run."""

    findings: list[Finding]
    metrics: StandMetrics
    metrics_source: str
    config: InspectionConfig
    area_ha: float
    canopy: dict[str, object] = field(default_factory=dict)
    detection: dict[str, object] | None = None
    match: dict[str, object] | None = None
    allometry: dict[str, object] | None = None
    gaps: list[CanopyGap] = field(default_factory=list)

    @property
    def max_severity(self) -> Severity:
        if not self.findings:
            return Severity.INFO
        return max(finding.severity for finding in self.findings)

    def counts(self) -> dict[str, int]:
        counts = {severity.label: 0 for severity in Severity}
        for finding in self.findings:
            counts[finding.severity.label] += 1
        return counts

    def findings_at_or_above(self, severity: Severity) -> list[Finding]:
        return [f for f in self.findings if f.severity >= severity]

    def as_dict(self) -> dict[str, object]:
        payload: dict[str, object] = {
            "area_ha": round(self.area_ha, 4),
            "metrics_source": self.metrics_source,
            "summary": {
                "finding_count": len(self.findings),
                "max_severity": self.max_severity.label,
                "counts": self.counts(),
            },
            "metrics": self.metrics.as_dict(),
            "canopy": self.canopy,
            "findings": [finding.as_dict() for finding in self.findings],
            "config": self.config.as_dict(),
        }
        if self.detection is not None:
            payload["detection"] = self.detection
        if self.match is not None:
            payload["match"] = self.match
        if self.allometry is not None:
            payload["allometry"] = self.allometry
        if self.gaps:
            payload["gaps"] = [gap.as_dict() for gap in self.gaps]
        return payload


def inspect_stand(
    *,
    chm: Grid | None = None,
    trees: Sequence[Tree] | None = None,
    area_ha: float | None = None,
    config: InspectionConfig | None = None,
    detection_config: DetectionConfig | None = None,
    run_detection: bool = True,
    allometry_model: str = DEFAULT_MODEL,
) -> InspectionReport:
    """Inspect a stand described by a canopy model, a tree list, or both.

    Args:
        chm: Canopy height model. Enables the canopy, gap and detection rules.
        trees: Field inventory records. Enables the data-quality rules.
        area_ha: Plot area; defaults to the raster footprint when a CHM is given.
        config: Threshold configuration.
        detection_config: Detection parameters used when ``run_detection``.
        run_detection: Set to ``False`` to skip tree detection on the raster.
        allometry_model: Height-diameter family fitted to the field data.

    Raises:
        ValueError: If neither a raster nor a tree list is supplied, or if the
            plot area cannot be determined.
    """
    if chm is None and trees is None:
        raise ValueError("inspect_stand needs a canopy height model, trees, or both")
    config = config or InspectionConfig()

    if area_ha is None:
        area_ha = chm.area_ha if chm is not None else None
    if area_ha is None or area_ha <= 0:
        raise ValueError("area_ha must be supplied and positive when no raster is given")

    detection: DetectionResult | None = None
    gaps: list[CanopyGap] = []
    if chm is not None:
        if run_detection:
            detection = detect_trees(chm, detection_config or DetectionConfig())
        gaps = find_gaps(
            chm,
            threshold=config.canopy_threshold,
            min_area=config.min_gap_area,
            min_width=config.min_gap_width,
        )

    # An inventory with no rows is still an inventory: somebody counted and
    # found nothing.  Only a missing inventory is missing data.
    field_trees = None if trees is None else list(trees)
    model: FittedModel | None = None
    species_models: dict[str, FittedModel] = {}
    fitted = False
    if field_trees:  # fitting needs actual measurements, not just a header row
        species_models = fit_by_species(
            field_trees, model=allometry_model, min_points=config.min_allometry_points
        )
        model = species_models[""]
        fitted = model.n > 0
    elif detection is not None:
        model = default_model("")

    if field_trees is not None:
        analysis_trees = field_trees
        source = "field"
    elif detection is not None:
        analysis_trees = trees_from_crowns(detection.crowns, dbh_model=model)
        source = "detected"
    else:
        # Only a raster, with detection switched off: nothing looked for stems,
        # so the summary below records that the stocking is unknown.
        analysis_trees = []
        source = "none"

    metrics = (
        StandMetrics.unknown(area_ha)
        if source == "none"
        else stand_metrics(analysis_trees, area_ha)
    )
    match: MatchResult | None = None
    if detection is not None and field_trees is not None:
        match = match_trees(detection.crowns, field_trees, tolerance=config.match_tolerance_m)

    ctx = InspectionContext(
        config=config,
        area_ha=area_ha,
        trees=analysis_trees,
        metrics=metrics,
        metrics_source=source,
        field_trees=field_trees,
        chm=chm,
        detection=detection,
        gaps=gaps,
        match=match,
        model=model,
        species_models=species_models,
        fitted=fitted,
    )

    findings: list[Finding] = []
    for check in RULES:
        findings.extend(check(ctx))
    findings.sort(key=lambda f: (-int(f.severity), f.code, f.subject))

    canopy_summary: dict[str, object] = {}
    if chm is not None:
        stats = chm.stats()
        canopy_summary = {
            "cover": round(canopy_cover(chm, config.canopy_threshold), 4),
            "gap_fraction": round(gap_fraction(chm, threshold=config.canopy_threshold), 4),
            "rugosity_m": round(rugosity(chm, threshold=config.canopy_threshold), 3),
            "mean_height_m": _round(stats.mean, 3),
            "max_height_m": _round(stats.maximum, 3),
            "strata": {k: round(v, 4) for k, v in vertical_strata(chm).items()},
            "gap_count": len(gaps),
        }

    return InspectionReport(
        findings=findings,
        metrics=metrics,
        metrics_source=source,
        config=config,
        area_ha=area_ha,
        canopy=canopy_summary,
        detection=detection.as_dict() if detection is not None else None,
        match=match.as_dict() if match is not None else None,
        allometry=(
            model.as_dict()
            | {
                "fitted": fitted,
                "by_species": {
                    code: sub.as_dict() for code, sub in sorted(species_models.items()) if code
                },
            }
        )
        if model is not None
        else None,
        gaps=gaps,
    )


def _round(value: float | None, digits: int = 4) -> float | None:
    if value is None:
        return None
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    return round(value, digits)
