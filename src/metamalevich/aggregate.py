"""Explicit aggregations from raw taxonomic counts to colour weights.

No function in this module silently picks a majority. Callers name the method.
"""

from __future__ import annotations

from metamalevich.evidence import AGGREGATIONS
from metamalevich.taxonomy import Taxonomy


def aggregate(counts: dict[int, float], method: str, taxonomy: Taxonomy | None = None, *, min_fraction: float = 0.0) -> dict[int, float]:
    """Return weights for one graph element.

    ``count`` and ``weighted_count`` keep the numeric counts.
    ``probability_sum`` divides by the total.
    ``majority`` gives equal weight to every taxon tied for the maximum count.
    ``lca`` places weight 1 on the LCA of taxa at or above ``min_fraction``.
    """
    if method not in AGGREGATIONS:
        raise ValueError(f"unknown aggregation method: {method}")
    cleaned = {int(taxon_id): float(value) for taxon_id, value in counts.items() if float(value) > 0}
    if method in {"count", "weighted_count"}:
        return cleaned
    total = sum(cleaned.values())
    if total <= 0:
        return {}
    if method == "probability_sum":
        return {taxon_id: value / total for taxon_id, value in cleaned.items()}
    if method == "majority":
        best = max(cleaned.values())
        tied = [taxon_id for taxon_id, value in cleaned.items() if value == best]
        share = 1.0 / len(tied)
        return {taxon_id: share for taxon_id in tied}
    if taxonomy is None:
        raise ValueError("lca aggregation requires a taxonomy")
    supported = [taxon_id for taxon_id, value in cleaned.items() if taxon_id != 0 and value / total >= min_fraction]
    if not supported:
        return {0: 1.0}
    return {taxonomy.lca(supported): 1.0}


def roll_counts(counts: dict[int, float], taxonomy: Taxonomy, rank: str = "S") -> dict[int, float]:
    """Sum counts onto ``rank``. Mass with no ancestor at that rank becomes taxon 0."""
    rolled: dict[int, float] = {}
    for taxon_id, value in counts.items():
        if value <= 0:
            continue
        if int(taxon_id) == 0:
            destination: int | None = 0
        else:
            destination = taxonomy.ancestor_at_rank(int(taxon_id), rank)
        key = 0 if destination is None else destination
        rolled[key] = rolled.get(key, 0.0) + float(value)
    return {taxon_id: value for taxon_id, value in rolled.items() if value > 0}


def normalize(weights: dict[int, float]) -> dict[int, float]:
    """Divide by the sum. An empty or zero map stays empty."""
    total = sum(value for value in weights.values() if value > 0)
    if total <= 0:
        return {}
    return {taxon_id: value / total for taxon_id, value in weights.items() if value > 0}
