"""Scores against a closed community with known node taxa and abundances."""

from __future__ import annotations

import math

from metamalevich.resolve import argmax_taxon


def _l1(left: dict[int, float], right: dict[int, float]) -> float:
    keys = set(left) | set(right)
    return sum(abs(left.get(key, 0.0) - right.get(key, 0.0)) for key in keys)


def _bray_curtis(left: dict[int, float], right: dict[int, float]) -> float:
    keys = set(left) | set(right)
    union = sum(left.get(key, 0.0) + right.get(key, 0.0) for key in keys)
    if union <= 0:
        return 0.0
    shared = sum(min(left.get(key, 0.0), right.get(key, 0.0)) for key in keys)
    return 1.0 - (2.0 * shared / union)


def _pearson(xs: list[float], ys: list[float]) -> float | None:
    n = len(xs)
    if n < 2:
        return None
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    num = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    den_x = math.sqrt(sum((x - mean_x) ** 2 for x in xs))
    den_y = math.sqrt(sum((y - mean_y) ** 2 for y in ys))
    if den_x == 0 or den_y == 0:
        return None
    return num / (den_x * den_y)


def _rank(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda index: values[index])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start
        while end + 1 < len(order) and values[order[end + 1]] == values[order[start]]:
            end += 1
        average = (start + end) / 2.0 + 1.0
        for index in order[start : end + 1]:
            ranks[index] = average
        start = end + 1
    return ranks


def classification_scores(
    predicted: dict[str, dict[int, float]],
    truth: dict[str, int],
) -> dict[str, float | None]:
    """Node-level precision, recall, F1, and accuracy at the supplied taxon ids.

    A node is a true positive for its truth taxon when the argmax matches.
    Macro scores average taxa that appear in the truth. A truth taxon that
    receives no predictions contributes precision, recall, and F1 of 0.
    """
    labels = sorted({taxon_id for taxon_id in truth.values() if taxon_id != 0})
    per_taxon = []
    correct = 0
    used = 0
    for node_id, taxon_id in truth.items():
        if node_id not in predicted:
            continue
        used += 1
        guess, _prob = argmax_taxon(predicted[node_id])
        if guess == taxon_id and taxon_id != 0:
            correct += 1
    for taxon_id in labels:
        tp = fp = fn = 0
        for node_id, true_taxon in truth.items():
            if node_id not in predicted:
                continue
            guess, _prob = argmax_taxon(predicted[node_id])
            if guess == taxon_id and true_taxon == taxon_id:
                tp += 1
            elif guess == taxon_id and true_taxon != taxon_id:
                fp += 1
            elif true_taxon == taxon_id and guess != taxon_id:
                fn += 1
        if (tp + fp) == 0 and (tp + fn) == 0:
            continue
        # A truth species with no predicted nodes scores 0. Dropping it would
        # inflate the macro average toward the species that were predicted.
        precision = tp / (tp + fp) if (tp + fp) else 0.0
        recall = tp / (tp + fn) if (tp + fn) else 0.0
        if precision + recall == 0:
            f1 = 0.0
        else:
            f1 = 2 * precision * recall / (precision + recall)
        per_taxon.append((precision, recall, f1))

    def _mean(index: int) -> float | None:
        if not per_taxon:
            return None
        return sum(row[index] for row in per_taxon) / len(per_taxon)

    return {
        "accuracy": (correct / used) if used else None,
        "precision": _mean(0),
        "recall": _mean(1),
        "f1": _mean(2),
        "n_nodes": float(used),
    }


def abundance_scores(predicted: dict[int, float], truth: dict[int, float]) -> dict[str, float | None]:
    """L1, Bray-Curtis, Pearson, and Spearman on relative abundances.

    ``truth`` is rescaled to sum to 1. ``predicted`` is used as given, so an
    unclassified fraction that was left out of the taxon rows must already be
    included (taxon 0) if it should affect L1.
    """
    truth_sum = sum(truth.values())
    truth_rel = {key: value / truth_sum for key, value in truth.items()} if truth_sum else {}
    pred_rel = {key: float(value) for key, value in predicted.items() if value > 0}
    keys = sorted(set(truth_rel) | set(pred_rel))
    xs = [truth_rel.get(key, 0.0) for key in keys]
    ys = [pred_rel.get(key, 0.0) for key in keys]
    return {
        "l1": _l1(pred_rel, truth_rel),
        "bray_curtis": _bray_curtis(pred_rel, truth_rel),
        "pearson": _pearson(xs, ys),
        "spearman": _pearson(_rank(xs), _rank(ys)),
    }


def neighbour_agreement(edges: list[dict], predicted: dict[str, dict[int, float]]) -> float | None:
    """Fraction of edges whose endpoint argmax taxa are equal and not unclassified.

    These edges are composition neighbours. A smoother that mixes a node with
    its neighbours raises this number by construction, so it is not an
    independent measure of taxonomic correctness.
    """
    if not edges:
        return None
    same = 0
    used = 0
    for edge in edges:
        if edge["source"] not in predicted or edge["target"] not in predicted:
            continue
        used += 1
        left, _ = argmax_taxon(predicted[edge["source"]])
        right, _ = argmax_taxon(predicted[edge["target"]])
        if left != 0 and left == right:
            same += 1
    if used == 0:
        return None
    return same / used
