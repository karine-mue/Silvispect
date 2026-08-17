"""Height-diameter allometry: model families, fitting, and residual screening.

The relationship between diameter at breast height (DBH) and total height is
the workhorse of forest mensuration.  Silvispect uses it twice:

* to estimate a diameter for a remotely detected tree, which has a height but
  no measurable stem, and
* to screen field measurements — a stem whose height departs strongly from the
  stand's own fitted curve is either damaged, mis-measured, or mis-linked.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass

from .inventory import Tree
from .optimize import nelder_mead

__all__ = [
    "DEFAULT_MODEL",
    "MODELS",
    "SPECIES_DEFAULTS",
    "AllometryError",
    "FittedModel",
    "ModelSpec",
    "default_model",
    "fit_by_species",
    "fit_height_diameter",
    "paired_measurements",
    "residual_scores",
]

BREAST_HEIGHT = 1.3


class AllometryError(ValueError):
    """Raised when a model cannot be built or fitted."""


@dataclass(frozen=True)
class ModelSpec:
    """A height-diameter model family."""

    name: str
    description: str
    param_names: tuple[str, ...]
    predict: Callable[[float, Sequence[float]], float]
    initial: tuple[float, ...]

    def evaluate(self, dbh_cm: float, params: Sequence[float]) -> float:
        return self.predict(dbh_cm, params)


def _curtis(dbh: float, p: Sequence[float]) -> float:
    if dbh <= 0:
        return BREAST_HEIGHT
    return BREAST_HEIGHT + p[0] * (dbh / (1.0 + dbh)) ** p[1]


def _naslund(dbh: float, p: Sequence[float]) -> float:
    if dbh <= 0:
        return BREAST_HEIGHT
    denominator = p[0] + p[1] * dbh
    if denominator <= 0:
        raise ValueError("invalid Naslund parameters")
    return BREAST_HEIGHT + dbh * dbh / (denominator * denominator)


def _power(dbh: float, p: Sequence[float]) -> float:
    if dbh <= 0:
        return 0.0
    if p[0] <= 0:
        raise ValueError("invalid power-model scale")
    return p[0] * dbh ** p[1]


def _chapman_richards(dbh: float, p: Sequence[float]) -> float:
    if dbh <= 0:
        return BREAST_HEIGHT
    if p[0] <= 0 or p[1] <= 0:
        raise ValueError("invalid Chapman-Richards parameters")
    base = 1.0 - math.exp(-p[1] * dbh)
    return BREAST_HEIGHT + p[0] * base ** p[2]


#: Available model families keyed by name.
MODELS: dict[str, ModelSpec] = {
    "curtis": ModelSpec(
        name="curtis",
        description="H = 1.3 + a * (D / (1 + D))^b",
        param_names=("a", "b"),
        predict=_curtis,
        initial=(30.0, 3.0),
    ),
    "naslund": ModelSpec(
        name="naslund",
        description="H = 1.3 + D^2 / (a + b * D)^2",
        param_names=("a", "b"),
        predict=_naslund,
        initial=(1.45, 0.178),
    ),
    "power": ModelSpec(
        name="power",
        description="H = a * D^b",
        param_names=("a", "b"),
        predict=_power,
        initial=(2.0, 0.6),
    ),
    "chapman": ModelSpec(
        name="chapman",
        description="H = 1.3 + a * (1 - exp(-b * D))^c",
        param_names=("a", "b", "c"),
        predict=_chapman_richards,
        initial=(30.0, 0.04, 1.0),
    ),
}

#: Default family used when nothing else is specified.  Naslund's equation is
#: preferred because it passes near the origin and saturates at ``1.3 + 1/b^2``,
#: so a plausible asymptotic height is readable straight off the parameters.
DEFAULT_MODEL = "naslund"

#: Fallback Naslund parameters for stands with too few paired measurements.
#: Codes follow the four-letter genus/species convention; the comment gives the
#: implied asymptotic height.
SPECIES_DEFAULTS: dict[str, tuple[float, float]] = {
    "PIAB": (1.300, 0.168),  # Picea abies, Norway spruce, ~36.7 m
    "PISY": (1.550, 0.185),  # Pinus sylvestris, Scots pine, ~30.5 m
    "FASY": (1.400, 0.170),  # Fagus sylvatica, European beech, ~35.9 m
    "QURO": (1.700, 0.195),  # Quercus robur, pedunculate oak, ~27.6 m
    "PSME": (1.200, 0.152),  # Pseudotsuga menziesii, Douglas fir, ~44.6 m
    "BEPE": (1.600, 0.200),  # Betula pendula, silver birch, ~26.3 m
    "": (1.450, 0.178),  # generic broadleaf/conifer mix, ~32.9 m
}


@dataclass(frozen=True)
class FittedModel:
    """A model family bound to concrete parameters plus fit diagnostics."""

    spec: ModelSpec
    params: tuple[float, ...]
    n: int
    rmse: float
    r_squared: float
    converged: bool

    @property
    def name(self) -> str:
        return self.spec.name

    def predict(self, dbh_cm: float) -> float:
        """Predict total height in metres for a diameter in centimetres."""
        return self.spec.evaluate(dbh_cm, self.params)

    def invert(self, height_m: float, *, lower: float = 0.1, upper: float = 400.0) -> float | None:
        """Back out the diameter implied by a height, or ``None`` if outside range.

        The supported families are monotonically increasing in diameter, so a
        bisection is both sufficient and robust.
        """
        try:
            low_value = self.predict(lower)
            high_value = self.predict(upper)
        except ValueError:  # pragma: no cover - guarded by fit validation
            return None
        if height_m <= low_value or height_m >= high_value:
            return None
        low, high = lower, upper
        for _ in range(200):
            mid = (low + high) / 2.0
            if self.predict(mid) < height_m:
                low = mid
            else:
                high = mid
            if high - low < 1e-6:
                break
        return (low + high) / 2.0

    def residuals(self, pairs: Sequence[tuple[float, float]]) -> list[float]:
        """Return ``observed - predicted`` heights for ``(dbh, height)`` pairs."""
        return [height - self.predict(dbh) for dbh, height in pairs]

    def as_dict(self) -> dict[str, object]:
        return {
            "model": self.spec.name,
            "equation": self.spec.description,
            "params": {
                name: round(value, 6)
                for name, value in zip(self.spec.param_names, self.params, strict=True)
            },
            "n": self.n,
            # A default curve has no fit statistics; JSON has no NaN literal,
            # so report the absence as null rather than emitting invalid JSON.
            "rmse_m": _finite(self.rmse),
            "r_squared": _finite(self.r_squared),
            "converged": self.converged,
        }


def paired_measurements(
    trees: Sequence[Tree], *, live_only: bool = True
) -> list[tuple[float, float]]:
    """Extract usable ``(dbh_cm, height_m)`` pairs from inventory records."""
    pairs: list[tuple[float, float]] = []
    for tree in trees:
        if live_only and not tree.is_live:
            continue
        if tree.dbh_cm is None or tree.height_m is None:
            continue
        if tree.dbh_cm <= 0 or tree.height_m <= 0:
            continue
        pairs.append((tree.dbh_cm, tree.height_m))
    return pairs


def fit_height_diameter(
    data: Sequence[Tree] | Sequence[tuple[float, float]],
    *,
    model: str = DEFAULT_MODEL,
    min_points: int = 6,
) -> FittedModel:
    """Fit a height-diameter curve by least squares.

    Args:
        data: Either inventory records or ready-made ``(dbh, height)`` pairs.
        model: Key into :data:`MODELS`.
        min_points: Refuse to fit fewer paired observations than this.

    Raises:
        AllometryError: If the family is unknown or there is too little data.
    """
    spec = MODELS.get(model)
    if spec is None:
        raise AllometryError(f"unknown model {model!r}; choose from {sorted(MODELS)}")
    pairs = _as_pairs(data)
    if len(pairs) < min_points:
        raise AllometryError(f"need at least {min_points} paired measurements, got {len(pairs)}")

    def cost(params: Sequence[float]) -> float:
        total = 0.0
        for dbh, height in pairs:
            total += (height - spec.evaluate(dbh, params)) ** 2
        return total

    result = nelder_mead(cost, spec.initial, step=0.25, max_iterations=4000)
    # A second restart from the incumbent tightens the simplex around the optimum.
    refined = nelder_mead(cost, result.x, step=0.02, max_iterations=4000)
    best = refined if refined.fun <= result.fun else result

    n = len(pairs)
    sse = best.fun
    mean_height = sum(height for _, height in pairs) / n
    sst = sum((height - mean_height) ** 2 for _, height in pairs)
    r_squared = 1.0 - sse / sst if sst > 0 else 0.0
    return FittedModel(
        spec=spec,
        params=tuple(best.x),
        n=n,
        rmse=math.sqrt(sse / n),
        r_squared=r_squared,
        converged=best.converged,
    )


def fit_by_species(
    trees: Sequence[Tree],
    *,
    model: str = DEFAULT_MODEL,
    min_points: int = 8,
) -> dict[str, FittedModel]:
    """Fit one curve per species plus a pooled fallback.

    A single pooled curve on a mixed stand systematically flags the minority
    species as outliers, so Silvispect fits every species with enough paired
    measurements separately.  The pooled curve is returned under the empty
    string and is used for species that are too rare to fit.

    Returns:
        Mapping of species code to fitted model, always containing ``""``.
    """
    groups: dict[str, list[Tree]] = {}
    for tree in trees:
        groups.setdefault(tree.species.strip().upper(), []).append(tree)

    models: dict[str, FittedModel] = {}
    try:
        models[""] = fit_height_diameter(trees, model=model, min_points=min_points)
    except AllometryError:
        dominant = max(sorted(groups), key=lambda code: len(groups[code])) if groups else ""
        models[""] = default_model(dominant)

    for code, members in groups.items():
        if not code:
            continue
        try:
            models[code] = fit_height_diameter(members, model=model, min_points=min_points)
        except AllometryError:
            continue
    return models


def default_model(species: str = "", *, model: str = DEFAULT_MODEL) -> FittedModel:
    """Return a literature-default curve for a species code.

    Only the Naslund family has tabulated defaults; other families fall back to
    their generic starting parameters.  A default curve is flagged by ``n = 0``.
    """
    spec = MODELS.get(model)
    if spec is None:
        raise AllometryError(f"unknown model {model!r}")
    code = species.strip().upper()
    if model == DEFAULT_MODEL:
        params = SPECIES_DEFAULTS.get(code, SPECIES_DEFAULTS[""])
    else:
        params = spec.initial
    return FittedModel(
        spec=spec,
        params=tuple(params),
        n=0,
        rmse=float("nan"),
        r_squared=float("nan"),
        converged=True,
    )


def residual_scores(model: FittedModel, trees: Sequence[Tree]) -> list[tuple[Tree, float, float]]:
    """Return ``(tree, residual, z_score)`` for every tree with both dimensions.

    The z-score uses the standard deviation of the residuals themselves, so it
    measures departure from the stand's own curve rather than from an external
    expectation.
    """
    usable = [
        tree
        for tree in trees
        if tree.dbh_cm and tree.height_m and tree.dbh_cm > 0 and tree.height_m > 0
    ]
    if not usable:
        return []
    residuals = [tree.height_m - model.predict(tree.dbh_cm) for tree in usable]  # type: ignore[arg-type]
    n = len(residuals)
    mean = sum(residuals) / n
    if n > 1:
        variance = sum((r - mean) ** 2 for r in residuals) / (n - 1)
        stdev = math.sqrt(variance)
    else:
        stdev = 0.0
    scored: list[tuple[Tree, float, float]] = []
    for tree, residual in zip(usable, residuals, strict=True):
        z = 0.0 if stdev == 0 else (residual - mean) / stdev
        scored.append((tree, residual, z))
    return scored


def _finite(value: float, digits: int = 4) -> float | None:
    """Round a float, mapping NaN and infinity onto ``None``."""
    if value is None or math.isnan(value) or math.isinf(value):
        return None
    return round(value, digits)


def _as_pairs(data: Sequence[Tree] | Sequence[tuple[float, float]]) -> list[tuple[float, float]]:
    if not data:
        return []
    first = data[0]
    if isinstance(first, Tree):
        return paired_measurements(data)  # type: ignore[arg-type]
    pairs: list[tuple[float, float]] = []
    for item in data:
        dbh, height = item  # type: ignore[misc]
        if dbh > 0 and height > 0:
            pairs.append((float(dbh), float(height)))
    return pairs
