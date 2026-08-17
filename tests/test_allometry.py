"""Tests for height-diameter allometry and the optimiser behind it."""

from __future__ import annotations

import math

import pytest

from silvispect.allometry import (
    MODELS,
    SPECIES_DEFAULTS,
    AllometryError,
    default_model,
    fit_by_species,
    fit_height_diameter,
    paired_measurements,
    residual_scores,
)
from silvispect.inventory import Tree
from silvispect.optimize import nelder_mead

DIAMETERS = tuple(float(d) for d in range(8, 60, 2))


def naslund_pairs(a: float, b: float, diameters=DIAMETERS):
    spec = MODELS["naslund"]
    return [(float(d), spec.evaluate(float(d), (a, b))) for d in diameters]


def test_all_models_are_monotonic_and_sane():
    for name, spec in MODELS.items():
        heights = [spec.evaluate(d, spec.initial) for d in (5.0, 15.0, 30.0, 60.0)]
        assert heights == sorted(heights), f"{name} is not increasing in diameter"
        assert 0.0 < heights[-1] < 100.0, f"{name} predicts an implausible height"


def test_species_defaults_are_plausible():
    for code in SPECIES_DEFAULTS:
        model = default_model(code)
        assert 5.0 < model.predict(10.0) < 20.0
        assert 15.0 < model.predict(30.0) < 40.0
        assert model.predict(60.0) > model.predict(30.0)


def test_fit_recovers_known_parameters():
    model = fit_height_diameter(naslund_pairs(1.3, 0.17), model="naslund")
    assert model.params[0] == pytest.approx(1.3, abs=0.01)
    assert model.params[1] == pytest.approx(0.17, abs=0.005)
    assert model.rmse < 1e-3
    assert model.r_squared > 0.999


def test_fit_recovers_curtis_parameters():
    spec = MODELS["curtis"]
    pairs = [(float(d), spec.evaluate(float(d), (30.0, 6.0))) for d in range(8, 60, 2)]
    model = fit_height_diameter(pairs, model="curtis")
    assert model.predict(25.0) == pytest.approx(spec.evaluate(25.0, (30.0, 6.0)), abs=0.05)


def test_fit_tolerates_noise():
    import random

    rng = random.Random(3)
    pairs = [(d, h + rng.gauss(0.0, 0.5)) for d, h in naslund_pairs(1.3, 0.17)]
    model = fit_height_diameter(pairs, model="naslund")
    assert model.predict(30.0) == pytest.approx(
        MODELS["naslund"].evaluate(30.0, (1.3, 0.17)), abs=1.0
    )
    assert model.r_squared > 0.9


def test_fit_from_tree_records():
    trees = [
        Tree(f"t{i}", 0, 0, dbh_cm=d, height_m=h)
        for i, (d, h) in enumerate(naslund_pairs(1.4, 0.18))
    ]
    model = fit_height_diameter(trees)
    assert model.n == len(trees)
    assert model.as_dict()["model"] == "naslund"
    assert set(model.as_dict()["params"]) == {"a", "b"}


def test_fit_requires_enough_points():
    with pytest.raises(AllometryError, match="at least"):
        fit_height_diameter(naslund_pairs(1.3, 0.17, diameters=(10.0, 12.0)))


def test_fit_rejects_unknown_model():
    with pytest.raises(AllometryError, match="unknown model"):
        fit_height_diameter(naslund_pairs(1.3, 0.17), model="magic")


def test_paired_measurements_filters_unusable():
    trees = [
        Tree("a", 0, 0, dbh_cm=20.0, height_m=15.0),
        Tree("b", 0, 0, dbh_cm=None, height_m=15.0),
        Tree("c", 0, 0, dbh_cm=20.0, height_m=None),
        Tree("d", 0, 0, dbh_cm=-1.0, height_m=15.0),
        Tree("e", 0, 0, dbh_cm=20.0, height_m=15.0, status="dead"),
    ]
    assert paired_measurements(trees) == [(20.0, 15.0)]
    assert len(paired_measurements(trees, live_only=False)) == 2


def test_invert_round_trips():
    model = default_model("PIAB")
    for dbh in (12.0, 25.0, 45.0):
        height = model.predict(dbh)
        assert model.invert(height) == pytest.approx(dbh, abs=0.01)


def test_invert_outside_range_returns_none():
    model = default_model("PIAB")
    assert model.invert(0.5) is None
    assert model.invert(500.0) is None


def test_residual_scores_flag_the_odd_tree():
    trees = [
        Tree(f"t{i}", 0, 0, dbh_cm=d, height_m=h)
        for i, (d, h) in enumerate(naslund_pairs(1.3, 0.17))
    ]
    model = fit_height_diameter(trees)
    trees.append(Tree("broken", 0, 0, dbh_cm=40.0, height_m=8.0))
    scored = residual_scores(model, trees)
    worst = max(scored, key=lambda item: abs(item[2]))
    assert worst[0].tree_id == "broken"
    assert worst[1] < 0
    assert abs(worst[2]) > 2.0


def test_residual_scores_without_usable_trees():
    assert residual_scores(default_model(), [Tree("a", 0, 0)]) == []


def test_residual_scores_with_a_single_tree():
    scored = residual_scores(default_model(), [Tree("a", 0, 0, dbh_cm=20.0, height_m=15.0)])
    assert len(scored) == 1
    assert scored[0][2] == 0.0


def test_fit_by_species_separates_curves():
    trees = []
    for code, (a, b) in (("PIAB", (1.30, 0.168)), ("BEPE", (1.60, 0.200))):
        for i, (d, h) in enumerate(naslund_pairs(a, b)):
            trees.append(Tree(f"{code}{i}", 0, 0, species=code, dbh_cm=d, height_m=h))
    models = fit_by_species(trees)
    assert set(models) == {"", "PIAB", "BEPE"}
    assert models["PIAB"].predict(30.0) > models["BEPE"].predict(30.0)
    assert models["PIAB"].rmse < models[""].rmse


def test_fit_by_species_falls_back_when_data_is_thin():
    models = fit_by_species([Tree("a", 0, 0, species="PIAB", dbh_cm=20.0, height_m=15.0)])
    assert set(models) == {""}
    assert models[""].n == 0  # a default curve, not a fit


def test_default_model_for_other_families():
    model = default_model("PIAB", model="power")
    assert model.name == "power"
    assert math.isnan(model.rmse)
    with pytest.raises(AllometryError):
        default_model("PIAB", model="nope")


def test_nelder_mead_finds_a_quadratic_minimum():
    result = nelder_mead(lambda p: (p[0] - 3.0) ** 2 + (p[1] + 1.0) ** 2, [0.0, 0.0])
    assert result.x[0] == pytest.approx(3.0, abs=1e-4)
    assert result.x[1] == pytest.approx(-1.0, abs=1e-4)
    assert result.fun < 1e-8
    assert result.converged
    assert result.evaluations > 0


def test_nelder_mead_survives_invalid_regions():
    def objective(p):
        if p[0] < 0:
            return float("inf")
        return (p[0] - 2.0) ** 2

    assert nelder_mead(objective, [1.0]).x[0] == pytest.approx(2.0, abs=1e-3)


def test_nelder_mead_handles_exceptions_and_nan():
    assert nelder_mead(lambda p: 1 / 0 if p[0] < 0 else p[0] ** 2, [1.0]).fun < 1e-6
    assert nelder_mead(lambda p: float("nan") if p[0] < 0 else p[0] ** 2, [1.0]).fun < 1e-6


def test_nelder_mead_requires_a_start():
    with pytest.raises(ValueError):
        nelder_mead(lambda p: 0.0, [])


def test_nelder_mead_respects_the_iteration_cap():
    result = nelder_mead(lambda p: (p[0] - 5.0) ** 2, [0.0], max_iterations=2)
    assert result.iterations <= 2
    assert not result.converged
