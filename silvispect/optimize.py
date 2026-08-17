"""A small derivative-free optimiser.

Silvispect fits non-linear height-diameter curves without SciPy, so it carries
its own Nelder-Mead simplex search.  The algorithm needs only function
evaluations, which keeps model definitions to a single lambda each.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

__all__ = ["OptimizeResult", "nelder_mead"]

Objective = Callable[[Sequence[float]], float]


@dataclass(frozen=True)
class OptimizeResult:
    """Outcome of a minimisation."""

    x: tuple[float, ...]
    fun: float
    iterations: int
    evaluations: int
    converged: bool


def nelder_mead(
    objective: Objective,
    x0: Sequence[float],
    *,
    step: float = 0.1,
    max_iterations: int = 2000,
    ftol: float = 1e-10,
    xtol: float = 1e-10,
    alpha: float = 1.0,
    gamma: float = 2.0,
    rho: float = 0.5,
    sigma: float = 0.5,
) -> OptimizeResult:
    """Minimise ``objective`` starting from ``x0``.

    Args:
        objective: Function of a parameter vector returning a scalar cost.
            Non-finite costs are treated as very large, so models may reject
            invalid parameters by returning ``float('inf')``.
        x0: Starting point.
        step: Relative size of the initial simplex.
        max_iterations: Hard iteration cap.
        ftol: Convergence tolerance on the spread of function values.
        xtol: Convergence tolerance on the spread of vertex coordinates.

    Returns:
        An :class:`OptimizeResult` with the best vertex found.
    """
    dimension = len(x0)
    if dimension == 0:
        raise ValueError("x0 must not be empty")

    evaluations = 0

    def evaluate(point: Sequence[float]) -> float:
        nonlocal evaluations
        evaluations += 1
        try:
            value = float(objective(point))
        except (ValueError, ZeroDivisionError, OverflowError):
            return float("inf")
        if value != value:  # NaN
            return float("inf")
        return value

    simplex: list[list[float]] = [list(x0)]
    for i in range(dimension):
        vertex = list(x0)
        delta = step * abs(vertex[i]) if vertex[i] != 0 else step
        vertex[i] += delta
        simplex.append(vertex)
    scores = [evaluate(vertex) for vertex in simplex]

    converged = False
    iterations = 0
    while iterations < max_iterations:
        iterations += 1
        order = sorted(range(len(simplex)), key=lambda i: scores[i])
        simplex = [simplex[i] for i in order]
        scores = [scores[i] for i in order]

        if _spread(scores) <= ftol and _vertex_spread(simplex) <= xtol:
            converged = True
            break

        centroid = [sum(vertex[i] for vertex in simplex[:-1]) / dimension for i in range(dimension)]
        worst = simplex[-1]

        reflected = [centroid[i] + alpha * (centroid[i] - worst[i]) for i in range(dimension)]
        reflected_score = evaluate(reflected)

        if reflected_score < scores[0]:
            expanded = [
                centroid[i] + gamma * (reflected[i] - centroid[i]) for i in range(dimension)
            ]
            expanded_score = evaluate(expanded)
            if expanded_score < reflected_score:
                simplex[-1], scores[-1] = expanded, expanded_score
            else:
                simplex[-1], scores[-1] = reflected, reflected_score
            continue

        if reflected_score < scores[-2]:
            simplex[-1], scores[-1] = reflected, reflected_score
            continue

        if reflected_score < scores[-1]:
            contracted = [
                centroid[i] + rho * (reflected[i] - centroid[i]) for i in range(dimension)
            ]
        else:
            contracted = [centroid[i] + rho * (worst[i] - centroid[i]) for i in range(dimension)]
        contracted_score = evaluate(contracted)
        if contracted_score < min(reflected_score, scores[-1]):
            simplex[-1], scores[-1] = contracted, contracted_score
            continue

        best = simplex[0]
        for i in range(1, len(simplex)):
            simplex[i] = [best[j] + sigma * (simplex[i][j] - best[j]) for j in range(dimension)]
            scores[i] = evaluate(simplex[i])

    order = sorted(range(len(simplex)), key=lambda i: scores[i])
    best_index = order[0]
    return OptimizeResult(
        x=tuple(simplex[best_index]),
        fun=scores[best_index],
        iterations=iterations,
        evaluations=evaluations,
        converged=converged,
    )


def _spread(scores: Sequence[float]) -> float:
    finite = [s for s in scores if s != float("inf")]
    if len(finite) < 2:
        return float("inf")
    return max(finite) - min(finite)


def _vertex_spread(simplex: Sequence[Sequence[float]]) -> float:
    dimension = len(simplex[0])
    return max(
        max(vertex[i] for vertex in simplex) - min(vertex[i] for vertex in simplex)
        for i in range(dimension)
    )
