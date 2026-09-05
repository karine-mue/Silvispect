"""Stand-level mensuration metrics computed from a list of trees.

All per-hectare quantities are expansions of the sampled trees onto the plot
area supplied by the caller; Silvispect never guesses the area, because an
unstated expansion factor is the classic way an inventory summary becomes
quietly wrong.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass, field
from decimal import Decimal, localcontext

from .grid import mean_of, rms_of, sum_of
from .inventory import Tree

__all__ = [
    "DEFAULT_FORM_FACTOR",
    "DEFAULT_WOOD_DENSITY",
    "MIN_CLASS_WIDTH_CM",
    "MetricsError",
    "StandMetrics",
    "above_ground_biomass",
    "basal_area",
    "diameter_distribution",
    "dominant_height",
    "gini_coefficient",
    "lorey_height",
    "quadratic_mean_diameter",
    "reineke_sdi",
    "shannon_index",
    "simpson_index",
    "species_composition",
    "stand_metrics",
    "stem_volume",
]

#: Cylindrical form factor used when converting basal area x height to volume.
DEFAULT_FORM_FACTOR = 0.5

#: Wood density in g/cm3 used by the biomass estimator when species is unknown.
DEFAULT_WOOD_DENSITY = 0.5

#: Wood densities in g/cm3 for the species codes Silvispect ships defaults for.
WOOD_DENSITY: dict[str, float] = {
    "PIAB": 0.40,
    "PISY": 0.42,
    "FASY": 0.58,
    "QURO": 0.56,
    "PSME": 0.45,
    "BEPE": 0.52,
}


class MetricsError(ValueError):
    """Raised when metrics are requested for impossible inputs."""


def _measured(trees: Sequence[Tree]) -> list[Tree]:
    return [tree for tree in trees if tree.dbh_cm is not None and tree.dbh_cm > 0]


def basal_area(trees: Sequence[Tree]) -> float:
    """Total cross-sectional area at breast height in square metres."""
    return sum_of([tree.basal_area_m2 or 0.0 for tree in trees])


def quadratic_mean_diameter(trees: Sequence[Tree]) -> float | None:
    """Diameter of the tree of mean basal area, in centimetres."""
    diameters = [tree.dbh_cm for tree in _measured(trees) if tree.dbh_cm is not None]
    if not diameters:
        return None
    return rms_of(diameters)


def _basal_area_parts(dbh_cm: float) -> tuple[float, int]:
    """Basal area in square metres as ``(mantissa, exponent)`` — never a float.

    ``pi * (dbh / 200)^2`` with the diameter split into digits and a power of
    two first: the digits are squared safely inside ``[0.25, 1)`` and the power
    of two is doubled as an integer, so a diameter of any accepted size has a
    weight the mean can use, where the float form underflows below about
    5e-160 cm and overflows above about 1.5e156 cm.
    """
    digits, scale = math.frexp(dbh_cm / 200.0)
    digits, extra = math.frexp(math.pi * digits * digits)
    return digits, 2 * scale + extra


def _shifted(digits: float, exponent: int) -> float:
    """``digits * 2**exponent``, reading a shift past the range as zero.

    A term more than a thousand binary orders below the largest one in a sum
    cannot change it, so its disappearance is the right answer rather than an
    error to raise.
    """
    if exponent < -2100:
        return 0.0
    return math.ldexp(digits, exponent)


def lorey_height(trees: Sequence[Tree]) -> float | None:
    """Basal-area weighted mean height, in metres.

    Weighting by basal area makes the statistic robust to the many small stems
    that dominate a raw arithmetic mean.
    """
    numerators: list[tuple[float, int]] = []
    denominators: list[tuple[float, int]] = []
    for tree in trees:
        if tree.dbh_cm is None or tree.dbh_cm <= 0:
            continue
        if tree.height_m is None or tree.height_m <= 0:
            continue
        # The weight is never formed as a float.  Basal area is ``pi (D/200)^2``,
        # and squaring a diameter of 2e-160 cm underflows to zero before any
        # careful summation can see it — which not only lost that stem's
        # weight but, because its enormous height still set the scale of the
        # numerator, shifted the other stem's real contribution out of range
        # and returned 0 for a mean of 7e-12.  Splitting the diameter with
        # ``frexp`` first, the square is a mantissa product and a doubled
        # exponent, and neither can leave the range.
        weight_digits, weight_scale = _basal_area_parts(tree.dbh_cm)
        height_digits, height_scale = math.frexp(tree.height_m)
        numerators.append((weight_digits * height_digits, weight_scale + height_scale))
        denominators.append((weight_digits, weight_scale))
    if not numerators:
        return None
    # Each sum is taken in units of its own largest term, so overflow is
    # impossible and only terms too small to change the sum are lost; the two
    # scales are put back with a single ``ldexp`` at the end — the one place
    # the answer is allowed to leave the range, and only if it truly does.

    top_scale = max(scale for _, scale in numerators)
    bottom_scale = max(scale for _, scale in denominators)
    # Summed exactly, so the answer is the same however the rows are ordered.
    # Running totals lose the small stems once a large one has been added, and
    # the loss depends on when it arrived: reversing a list of one 1e16 m stem
    # among twenty 1 m stems moved the result by a whole metre.
    top = math.fsum(_shifted(digits, scale - top_scale) for digits, scale in numerators)
    bottom = math.fsum(_shifted(digits, scale - bottom_scale) for digits, scale in denominators)
    if bottom == 0.0:
        return None
    return math.ldexp(top / bottom, top_scale - bottom_scale)


def dominant_height(trees: Sequence[Tree], area_ha: float) -> float | None:
    """Mean height of the 100 thickest stems per hectare ("top height")."""
    if area_ha <= 0:
        raise MetricsError("area_ha must be positive")
    candidates = [
        tree for tree in _measured(trees) if tree.height_m is not None and tree.height_m > 0
    ]
    if not candidates:
        return None
    count = max(1, min(len(candidates), round(100 * area_ha)))
    # Ties at the cutoff must not be settled by inventory row order: two 30 cm
    # stems of 10 m and 30 m gave whichever height happened to be listed first.
    # Rank by diameter, then height, then identifier — all properties of the
    # tree rather than of the file.
    ranked = sorted(
        candidates,
        key=lambda tree: (-(tree.dbh_cm or 0.0), -(tree.height_m or 0.0), tree.tree_id),
    )[:count]
    return mean_of([tree.height_m for tree in ranked])  # type: ignore[misc]


def reineke_sdi(stems_per_ha: float, qmd_cm: float) -> float:
    """Reineke's stand density index, referenced to a 25 cm quadratic mean diameter."""
    if stems_per_ha < 0 or qmd_cm <= 0:
        raise MetricsError("stems_per_ha must be non-negative and qmd positive")
    try:
        index = stems_per_ha * (qmd_cm / 25.0) ** 1.605
    except OverflowError:
        index = math.inf
    if not math.isfinite(index):
        raise MetricsError(
            f"stand density index is not representable for qmd {qmd_cm!r} cm "
            f"at {stems_per_ha!r} stems/ha"
        )
    return index


def gini_coefficient(values: Sequence[float]) -> float | None:
    """Gini coefficient of a non-negative distribution.

    Applied to tree basal areas it summarises structural heterogeneity: ``0``
    is a perfectly uniform even-aged stand, values above ~0.5 indicate a
    strongly size-differentiated, often uneven-aged structure.
    """
    positives = [v for v in values if v is not None and v >= 0]
    n = len(positives)
    if n < 2:
        return None
    ordered = sorted(positives)
    scale = max(ordered)
    if scale == 0.0:
        return 0.0
    # Ranks multiply the values, so the products are formed in units of the
    # largest value; the scale cancels between the two sums.
    shares = [value / scale for value in ordered]
    total = math.fsum(shares)
    if total == 0:
        return 0.0
    weighted = math.fsum((i + 1) * share for i, share in enumerate(shares))
    coefficient = (2.0 * weighted) / (n * total) - (n + 1.0) / n
    # A perfectly uniform distribution can land a few ulps below zero.
    return min(1.0, max(0.0, coefficient))


def species_composition(trees: Sequence[Tree]) -> dict[str, float]:
    """Proportion of stems per species code, descending by share."""
    labelled = [tree.species.strip().upper() or "UNKNOWN" for tree in trees]
    if not labelled:
        return {}
    counts = Counter(labelled)
    total = sum(counts.values())
    ordered = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    return {species: count / total for species, count in ordered}


def shannon_index(trees: Sequence[Tree]) -> float | None:
    """Shannon entropy of the species mixture (natural log)."""
    shares = species_composition(trees)
    if not shares:
        return None
    return -sum(p * math.log(p) for p in shares.values() if p > 0)


def simpson_index(trees: Sequence[Tree]) -> float | None:
    """Gini-Simpson diversity: the chance two random stems differ in species."""
    shares = species_composition(trees)
    if not shares:
        return None
    return 1.0 - sum(p * p for p in shares.values())


#: Narrowest diameter class Silvispect will bin into, in centimetres.
#:
#: An inventory is written to three decimals by default (see
#: ``docs/data-formats.md``), so a hundredth of a millimetre is already finer
#: than the record a class is cut from; a band below it separates stems the
#: file cannot tell apart.  Fixing the floor also bounds the arithmetic: the
#: widest band number is the largest accepted diameter over this, so the
#: decimal context below is a stated size rather than an open-ended one.
MIN_CLASS_WIDTH_CM = 1e-3

#: Digits the class arithmetic runs at.  ``MAX_DBH_CM / MIN_CLASS_WIDTH_CM``
#: is about 1.5e159, so 160 digits carry the largest band number this domain
#: can produce; the rest is headroom for the boundary labels computed from it.
_CLASS_DIGITS = 200


def _class_bound(value: Decimal) -> str:
    """Render a class boundary as the decimal it is, without a trailing zero."""
    text = format(value.normalize(), "f")
    return text if text else "0"


def _as_decimal(value: float) -> Decimal:
    """Read a float as the decimal a person would have written for it.

    ``Decimal(0.1)`` is the binary value, which is a shade above a tenth;
    ``Decimal("0.1")`` is a tenth.  Class widths and diameters are written by
    people in decimal, and ``repr`` gives back the shortest text that reads as
    the same float — so going through it recovers what was meant.
    """
    return Decimal(repr(value))


def diameter_distribution(trees: Sequence[Tree], *, class_width: float = 5.0) -> dict[str, int]:
    """Stem counts per diameter class, keyed by the class lower bound.

    The class a stem falls in is its band *number*, not a rounded-off
    centimetre value.  Truncating the lower bound to an integer collapsed every
    band narrower than a centimetre onto its neighbours — at a class width of
    0.5 cm, 0.2 cm and 0.7 cm were counted together in a class labelled
    ``0-0`` — so the width is kept as given and the boundaries are computed
    from the band number.

    The band number is worked out in decimal.  Binary floats put 0.3 / 0.1 a
    hair below 3, so a stem sitting exactly on a class boundary fell into the
    band below it and took a label like ``0.2-0.30000000000000004`` with it.
    Widths and diameters are written in decimal by the people who measure
    them, and dividing them as decimals puts the stem where they meant it.

    Raises:
        MetricsError: If ``class_width`` is below
            :data:`MIN_CLASS_WIDTH_CM`.  A band finer than the inventory
            records is not a measurement question, and the default decimal
            context could not answer it either: at 3e-28 cm the division ran
            out of digits and emitted a class whose own interval excluded the
            stem in it, and a shade narrower it raised
            ``decimal.InvalidOperation`` from inside the library.
    """
    if class_width < MIN_CLASS_WIDTH_CM:
        raise MetricsError(
            f"class_width must be at least {MIN_CLASS_WIDTH_CM} cm, got {class_width!r}"
        )
    width = _as_decimal(class_width)
    counts: Counter[int] = Counter()
    with localcontext() as context:
        context.prec = _CLASS_DIGITS
        for tree in _measured(trees):
            assert tree.dbh_cm is not None
            counts[int(_as_decimal(tree.dbh_cm) // width)] += 1
        return {
            f"{_class_bound(index * width)}-{_class_bound((index + 1) * width)}": counts[index]
            for index in sorted(counts)
        }


def stem_volume(tree: Tree, *, form_factor: float = DEFAULT_FORM_FACTOR) -> float | None:
    """Merchantable stem volume in cubic metres via the form-factor method."""
    area = tree.basal_area_m2
    if area is None or tree.height_m is None or tree.height_m <= 0:
        return None
    return form_factor * area * tree.height_m


def above_ground_biomass(tree: Tree, *, wood_density: float | None = None) -> float | None:
    """Above-ground biomass in kilograms (pantropical Chave et al. 2014 form).

    The equation ``AGB = 0.0673 * (rho * D^2 * H)^0.976`` is used with a wood
    density looked up from the species code when one is not supplied.
    """
    if tree.dbh_cm is None or tree.dbh_cm <= 0:
        return None
    if tree.height_m is None or tree.height_m <= 0:
        return None
    if wood_density is None:
        wood_density = WOOD_DENSITY.get(tree.species.strip().upper(), DEFAULT_WOOD_DENSITY)
    # ``0.0673 * (rho * D^2 * H) ** 0.976``, evaluated so that neither end of
    # the range is lost.  Forming the bracket directly overflows above about
    # 1.3e154 cm of diameter; distributing the exponent as
    # ``(rho * H) ** 0.976 * D ** 1.952`` moved the failure to the other end,
    # where a height of 5e-324 m sent ``rho * H`` to zero and the whole biomass
    # with it.  Both were answers an ordinary float can hold.
    #
    # Splitting the bracket into a mantissa and a power of two keeps every
    # intermediate in range: the mantissas multiply safely, the exponents add
    # as integers, and raising to the power is then
    # ``mantissa**0.976 * 2**(0.976 * exponent)`` — whose own integer part is
    # applied last, by the one operation that is allowed to leave the range.
    digits, exponent = 1.0, 0
    for value, times in ((wood_density, 1), (tree.dbh_cm, 2), (tree.height_m, 1)):
        mantissa, scale = math.frexp(value)
        digits *= mantissa**times
        exponent += scale * times
    power = 0.976 * exponent
    whole = math.floor(power)
    return math.ldexp(0.0673 * digits**0.976 * 2.0 ** (power - whole), whole)


@dataclass(frozen=True)
class StandMetrics:
    """A complete stand summary for one plot."""

    area_ha: float
    #: Stems considered.  ``None`` means nothing counted stems at all, which is
    #: different from having counted and found none.
    tree_count: int | None
    measured_count: int | None
    #: Stems carrying both a diameter and a height, the basis of the volume and
    #: biomass totals.  Lower than ``measured_count`` when heights are missing.
    yield_basis_count: int | None
    stems_per_ha: float | None
    basal_area_per_ha: float | None
    quadratic_mean_diameter_cm: float | None
    mean_diameter_cm: float | None
    mean_height_m: float | None
    lorey_height_m: float | None
    dominant_height_m: float | None
    max_height_m: float | None
    sdi: float | None
    gini_basal_area: float | None
    shannon: float | None
    simpson: float | None
    volume_per_ha_m3: float | None
    biomass_per_ha_t: float | None
    species_shares: dict[str, float] = field(default_factory=dict)
    diameter_classes: dict[str, int] = field(default_factory=dict)

    def as_dict(self) -> dict[str, object]:
        return {
            "area_ha": round(self.area_ha, 4),
            "tree_count": self.tree_count,
            "measured_count": self.measured_count,
            "yield_basis_count": self.yield_basis_count,
            "stems_per_ha": _round(self.stems_per_ha, 2),
            "basal_area_per_ha_m2": _round(self.basal_area_per_ha, 3),
            "qmd_cm": _round(self.quadratic_mean_diameter_cm, 2),
            "mean_dbh_cm": _round(self.mean_diameter_cm, 2),
            "mean_height_m": _round(self.mean_height_m, 2),
            "lorey_height_m": _round(self.lorey_height_m, 2),
            "dominant_height_m": _round(self.dominant_height_m, 2),
            "max_height_m": _round(self.max_height_m, 2),
            "sdi": _round(self.sdi, 1),
            "gini_basal_area": _round(self.gini_basal_area, 4),
            "shannon": _round(self.shannon, 4),
            "simpson": _round(self.simpson, 4),
            "volume_per_ha_m3": _round(self.volume_per_ha_m3, 2),
            "biomass_per_ha_t": _round(self.biomass_per_ha_t, 2),
            "species_shares": {k: round(v, 4) for k, v in self.species_shares.items()},
            "diameter_classes": dict(self.diameter_classes),
        }

    @classmethod
    def unknown(cls, area_ha: float) -> StandMetrics:
        """A summary for a plot where nothing counted stems.

        Every stem-derived figure is ``None`` rather than ``0``: an uncounted
        stand has an unknown stocking, and a zero would be indistinguishable
        from a stand that was counted and found empty.
        """
        if area_ha <= 0:
            raise MetricsError("area_ha must be positive")
        return cls(
            area_ha=area_ha,
            tree_count=None,
            measured_count=None,
            yield_basis_count=None,
            stems_per_ha=None,
            basal_area_per_ha=None,
            quadratic_mean_diameter_cm=None,
            mean_diameter_cm=None,
            mean_height_m=None,
            lorey_height_m=None,
            dominant_height_m=None,
            max_height_m=None,
            sdi=None,
            gini_basal_area=None,
            shannon=None,
            simpson=None,
            volume_per_ha_m3=None,
            biomass_per_ha_t=None,
        )


def stand_metrics(
    trees: Sequence[Tree],
    area_ha: float,
    *,
    live_only: bool = True,
    form_factor: float = DEFAULT_FORM_FACTOR,
    class_width: float = 5.0,
) -> StandMetrics:
    """Summarise a stand.

    Args:
        trees: Inventory records for the plot.
        area_ha: Ground area the records represent, in hectares.
        live_only: Exclude stems whose status is not a living state.
        form_factor: Form factor for the volume estimate.
        class_width: Width of the diameter classes in the distribution.

    Raises:
        MetricsError: If ``area_ha`` is not positive.
    """
    if area_ha <= 0:
        raise MetricsError("area_ha must be positive")
    selected = [tree for tree in trees if tree.is_live] if live_only else list(trees)
    measured = _measured(selected)

    diameters = [tree.dbh_cm for tree in measured if tree.dbh_cm is not None]
    heights = [tree.height_m for tree in selected if tree.height_m and tree.height_m > 0]
    qmd = quadratic_mean_diameter(measured)
    stems_ha = len(selected) / area_ha

    # Volume and biomass both need a height as well as a diameter, so they are
    # computed over a subset of the measured stems.  That subset is reported as
    # `yield_basis_count` and the totals stay None when it is empty: a stand
    # with diameters but no heights has an *unknown* volume, not a zero one.
    volumes = [
        volume
        for volume in (stem_volume(tree, form_factor=form_factor) for tree in measured)
        if volume is not None
    ]
    biomasses = [
        biomass
        for biomass in (above_ground_biomass(tree) for tree in measured)
        if biomass is not None
    ]

    return StandMetrics(
        area_ha=area_ha,
        tree_count=len(selected),
        measured_count=len(measured),
        stems_per_ha=stems_ha,
        basal_area_per_ha=basal_area(measured) / area_ha if measured else None,
        quadratic_mean_diameter_cm=qmd,
        mean_diameter_cm=mean_of(diameters) if diameters else None,
        mean_height_m=mean_of(heights) if heights else None,
        lorey_height_m=lorey_height(selected),
        dominant_height_m=dominant_height(selected, area_ha),
        max_height_m=max(heights) if heights else None,
        sdi=reineke_sdi(len(measured) / area_ha, qmd) if qmd else None,
        gini_basal_area=gini_coefficient([tree.basal_area_m2 or 0.0 for tree in measured]),
        shannon=shannon_index(selected),
        simpson=simpson_index(selected),
        yield_basis_count=len(volumes),
        volume_per_ha_m3=sum_of(volumes) / area_ha if volumes else None,
        biomass_per_ha_t=(sum_of(biomasses) / 1000.0) / area_ha if biomasses else None,
        species_shares=species_composition(selected),
        diameter_classes=diameter_distribution(measured, class_width=class_width),
    )


def _round(value: float | None, digits: int) -> float | None:
    return None if value is None else round(value, digits)
