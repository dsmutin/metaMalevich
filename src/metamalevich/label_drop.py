"""Drop taxon labels that neighbour agreement does not support.

``logistic_label_drop`` fits on the same graph it is applied to. The targets
are pseudo-labels from that graph's neighbour agreement, not ground-truth
taxa and not a held-out split. The coefficients are not an unbiased
classifier and must not be treated as one.
"""

from __future__ import annotations

import math

import numpy as np

from metamalevich.aggregate import normalize
from metamalevich.resolve import argmax_taxon


def false_positive_report(
    predicted: dict[str, dict[int, float]],
    truth: dict[str, int],
) -> dict[str, int | float]:
    """Score argmax guesses against truth labels.

    ``n`` counts nodes in ``predicted``. Each guess is ``argmax_taxon``, so
    ties break toward the smaller taxon id. A false positive is a guess that
    differs from that node's truth. ``false_positive_rate`` is that count
    divided by ``n``, and ``micro_precision`` is the matching fraction.
    ``n_false_negative`` counts mismatches whose truth taxon is not 0.

    Every predicted node must be present in ``truth``. Nodes that appear only
    in ``truth`` are not scored. An empty prediction map has rate 0.
    """
    n_false_positive = 0
    n_false_negative = 0
    n_correct = 0
    for node_id, weights in predicted.items():
        if node_id not in truth:
            raise ValueError(f"missing truth label for node {node_id}")
        guess, _probability = argmax_taxon(weights)
        actual = truth[node_id]
        if guess != actual:
            n_false_positive += 1
            if actual != 0:
                n_false_negative += 1
        else:
            n_correct += 1
    n = len(predicted)
    if n == 0:
        return {
            "n": 0,
            "n_false_positive": 0,
            "false_positive_rate": 0.0,
            "n_false_negative": 0,
            "micro_precision": 0.0,
        }
    return {
        "n": n,
        "n_false_positive": n_false_positive,
        "false_positive_rate": n_false_positive / n,
        "n_false_negative": n_false_negative,
        "micro_precision": n_correct / n,
    }


def drop_unsupported_labels(
    distributions: dict[str, dict[int, float]],
    edges: list[dict],
    *,
    min_own: float = 0.05,
    min_neighbour: float = 0.05,
) -> dict[str, dict[int, float]]:
    """Drop taxa that are neither confident locally nor supported by neighbours.

    Neighbours of a node are the targets of edges with that node as source.
    A missing edge weight defaults to 1. Non-positive weights are ignored.
    Support for a taxon is the weight-averaged probability of that taxon on
    the node's neighbours, or 0 when the node has no such neighbours.

    A taxon is kept when it is the argmax, its own probability is at least
    ``min_own``, or its neighbour support is at least ``min_neighbour``.
    Ties in the argmax break toward the smaller taxon id. The kept masses are
    renormalized. If that removes every taxon, the argmax is returned as a
    one-hot. An empty input distribution stays empty. Inputs are not modified.
    """
    neighbours = _neighbours(edges)
    resolved: dict[str, dict[int, float]] = {}
    for node_id, weights in distributions.items():
        resolved[node_id] = _drop_node(node_id, weights, neighbours, distributions, min_own, min_neighbour)
    return resolved


def logistic_label_drop(
    distributions: dict[str, dict[int, float]],
    edges: list[dict],
    *,
    iterations: int = 200,
    step: float = 0.5,
    l2: float = 1.0,
    threshold: float = 0.5,
) -> dict[str, dict[int, float]]:
    """Drop non-argmax taxa whose pseudo-label model predicts no support.

    The model is fit on this graph alone. A training row is a node-taxon pair
    with positive own probability. Taxon 0 is never a positive example.
    Features, each expected in ``[0, 1]``, are the own probability, the
    weight-averaged neighbour probability, and the fraction of positive-weight
    neighbours whose argmax is this taxon (0 when there are no neighbours).
    The pseudo-label is 1 only when the taxon is the node's argmax and that
    fraction is at least 0.5. It is 0 when the taxon is not the argmax, the
    fraction is 0, and the own probability is below 0.2. Other rows are
    ambiguous and are not used. Ground-truth labels are not used. Because the
    fit uses pseudo-labels from the same graph it filters, it is not an
    unbiased classifier.

    With at least four rows and both classes, binary logistic regression is
    fit by gradient ascent on the log-likelihood for ``iterations`` steps of
    size ``step``. The feature vector includes a leading bias of 1. L2
    penalizes the non-bias weights by subtracting ``0.5 * l2 * ||w||^2`` from
    the log-likelihood. Otherwise the filter falls back to
    ``drop_unsupported_labels``.

    After the fit, each non-argmax taxon is dropped when its predicted
    probability of the positive class is below ``threshold``. The argmax is
    always kept. The result is renormalized, with the same empty-distribution
    and last-taxon fallback as ``drop_unsupported_labels``. ``threshold`` must
    lie in ``[0, 1]``. Inputs are not modified.
    """
    if not 0.0 <= threshold <= 1.0:
        raise ValueError("threshold must be in [0, 1]")
    neighbours = _neighbours(edges)
    features, labels = _pseudo_rows(distributions, neighbours)
    if len(labels) < 4 or len(set(labels)) < 2:
        return drop_unsupported_labels(distributions, edges)
    coefficients = _fit_logistic(features, labels, iterations=iterations, step=step, l2=l2)
    resolved: dict[str, dict[int, float]] = {}
    for node_id, weights in distributions.items():
        resolved[node_id] = _apply_logistic(node_id, weights, neighbours, distributions, coefficients, threshold)
    return resolved


def _drop_node(
    node_id: str,
    weights: dict[int, float],
    neighbours: dict[str, list[tuple[str, float]]],
    distributions: dict[str, dict[int, float]],
    min_own: float,
    min_neighbour: float,
) -> dict[int, float]:
    """Keep the argmax and taxa with enough own mass or neighbour support."""
    if not weights:
        return {}
    support, _guesses = _neighbour_context(node_id, neighbours, distributions)
    best, _probability = argmax_taxon(weights)
    kept: dict[int, float] = {}
    for taxon_id, probability in weights.items():
        if taxon_id == best or probability >= min_own or support.get(taxon_id, 0.0) >= min_neighbour:
            kept[taxon_id] = probability
    return _finish(kept, best)


def _apply_logistic(
    node_id: str,
    weights: dict[int, float],
    neighbours: dict[str, list[tuple[str, float]]],
    distributions: dict[str, dict[int, float]],
    coefficients: np.ndarray,
    threshold: float,
) -> dict[int, float]:
    """Drop non-argmax taxa whose positive-class probability is below ``threshold``."""
    if not weights:
        return {}
    support, guesses = _neighbour_context(node_id, neighbours, distributions)
    best, _probability = argmax_taxon(weights)
    kept: dict[int, float] = {}
    for taxon_id, probability in weights.items():
        if taxon_id == best:
            if probability > 0:
                kept[taxon_id] = probability
            continue
        if probability <= 0:
            continue
        agreement = _agreement(guesses, taxon_id)
        score = _predict(coefficients, probability, support.get(taxon_id, 0.0), agreement)
        if score < threshold:
            continue
        kept[taxon_id] = probability
    return _finish(kept, best)


def _finish(kept: dict[int, float], best: int) -> dict[int, float]:
    """Renormalize, or keep a one-hot argmax when nothing positive remains."""
    normalized = normalize(kept)
    if normalized:
        return normalized
    return {best: 1.0}


def _pseudo_rows(
    distributions: dict[str, dict[int, float]],
    neighbours: dict[str, list[tuple[str, float]]],
) -> tuple[list[tuple[float, float, float]], list[float]]:
    """Build pseudo-labelled rows. Ambiguous pairs and taxon-0 positives are omitted."""
    features: list[tuple[float, float, float]] = []
    labels: list[float] = []
    for node_id, weights in distributions.items():
        if not weights:
            continue
        support, guesses = _neighbour_context(node_id, neighbours, distributions)
        best, _probability = argmax_taxon(weights)
        for taxon_id, probability in weights.items():
            if probability <= 0:
                continue
            agreement = _agreement(guesses, taxon_id)
            if taxon_id == best and agreement >= 0.5:
                if taxon_id == 0:
                    continue
                label = 1.0
            elif taxon_id != best and agreement == 0.0 and probability < 0.2:
                label = 0.0
            else:
                continue
            features.append((float(probability), float(support.get(taxon_id, 0.0)), agreement))
            labels.append(label)
    return features, labels


def _fit_logistic(
    features: list[tuple[float, float, float]],
    labels: list[float],
    *,
    iterations: int,
    step: float,
    l2: float,
) -> np.ndarray:
    """Gradient ascent on the Bernoulli log-likelihood with L2 on non-bias weights.

    Coefficients start at 0. Column 0 is a bias of 1. The penalty is
    ``0.5 * l2 * ||w[1:]||^2``, so its gradient contribution is ``l2 * w[1:]``.
    """
    design = np.ones((len(features), 4), dtype=float)
    for row_index, (own, support, agreement) in enumerate(features):
        design[row_index, 1] = own
        design[row_index, 2] = support
        design[row_index, 3] = agreement
    response = np.asarray(labels, dtype=float)
    coefficients = np.zeros(4, dtype=float)
    for _ in range(iterations):
        residual = response - _sigmoid(design @ coefficients)
        gradient = design.T @ residual
        gradient[1:] -= l2 * coefficients[1:]
        coefficients = coefficients + step * gradient
    return coefficients


def _predict(coefficients: np.ndarray, own: float, support: float, agreement: float) -> float:
    """Predicted probability that the pseudo-label is positive."""
    logit = float(
        coefficients[0] + coefficients[1] * own + coefficients[2] * support + coefficients[3] * agreement
    )
    if logit >= 0.0:
        if logit > 500.0:
            return 1.0
        return 1.0 / (1.0 + math.exp(-logit))
    if logit < -500.0:
        return 0.0
    exp_logit = math.exp(logit)
    return exp_logit / (1.0 + exp_logit)


def _sigmoid(logits: np.ndarray) -> np.ndarray:
    """Logistic function. Logits outside [-500, 500] saturate without overflow."""
    logits = np.asarray(logits, dtype=float)
    out = np.empty(logits.shape, dtype=float)
    positive = logits >= 0.0
    clipped_pos = np.minimum(logits[positive], 500.0)
    out[positive] = 1.0 / (1.0 + np.exp(-clipped_pos))
    clipped_neg = np.maximum(logits[~positive], -500.0)
    exp_logit = np.exp(clipped_neg)
    out[~positive] = exp_logit / (1.0 + exp_logit)
    return out


def _neighbours(edges: list[dict]) -> dict[str, list[tuple[str, float]]]:
    """Outgoing targets with positive weight. A missing weight defaults to 1."""
    neighbours: dict[str, list[tuple[str, float]]] = {}
    for edge in edges:
        weight = float(edge.get("weight", 1.0))
        if weight <= 0:
            continue
        neighbours.setdefault(edge["source"], []).append((edge["target"], weight))
    return neighbours


def _neighbour_context(
    node_id: str,
    neighbours: dict[str, list[tuple[str, float]]],
    distributions: dict[str, dict[int, float]],
) -> tuple[dict[int, float], list[int]]:
    """Neighbour support and each neighbour's argmax, in edge order.

    Support is the weight-averaged probability. Neighbours with an empty
    distribution are left out of that average. They still count in the argmax
    list, where an empty distribution has argmax 0. No neighbours yields an
    empty support map and an empty argmax list.
    """
    links = neighbours.get(node_id, [])
    total = 0.0
    accumulated: dict[int, float] = {}
    guesses: list[int] = []
    for other, weight in links:
        other_weights = distributions.get(other, {})
        guesses.append(argmax_taxon(other_weights)[0])
        if not other_weights:
            continue
        total += weight
        for taxon_id, probability in other_weights.items():
            accumulated[taxon_id] = accumulated.get(taxon_id, 0.0) + weight * probability
    if total <= 0:
        return {}, guesses
    support = {taxon_id: value / total for taxon_id, value in accumulated.items() if value > 0}
    return support, guesses


def _agreement(guesses: list[int], taxon_id: int) -> float:
    """Fraction of neighbours whose argmax is ``taxon_id``. Zero if there are none."""
    if not guesses:
        return 0.0
    hits = sum(guess == taxon_id for guess in guesses)
    if hits == 0:
        return 0.0
    return hits / len(guesses)
