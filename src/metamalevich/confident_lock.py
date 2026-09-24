"""Lock confident calls and adopt a same-genus neighbour vote.

A node at or above the confidence cutoff keeps its distribution. An uncertain
node may be replaced by the weight-averaged distribution of its neighbours
when that vote is strong and stays inside the same genus.
"""

from __future__ import annotations

from metamalevich.aggregate import normalize
from metamalevich.resolve import argmax_taxon
from metamalevich.taxonomy import Taxonomy


def confident_lock(
    distributions: dict[str, dict[int, float]],
    edges: list[dict],
    taxonomy: Taxonomy,
    *,
    confident: float = 0.8,
    vote_min: float = 0.75,
) -> dict[str, dict[int, float]]:
    """Keep confident calls and replace uncertain ones only inside the same genus.

    A node is copied unchanged when its argmax probability is at least
    ``confident``, or when it has no positive-weight outgoing neighbours.
    Otherwise the neighbours' distributions are averaged by edge weight, the
    same idea as ``resolve._vote``. That average replaces the node only when
    its argmax is not taxon 0, the vote probability is at least ``vote_min``,
    and both argmax taxa share a genus (``ancestor_at_rank`` at rank ``G``).
    A missing genus leaves the node unchanged. ``confident`` and ``vote_min``
    must lie in ``[0, 1]``. Inputs are not modified.
    """
    if not 0.0 <= confident <= 1.0:
        raise ValueError("confident must be in [0, 1]")
    if not 0.0 <= vote_min <= 1.0:
        raise ValueError("vote_min must be in [0, 1]")

    neighbours = _neighbours(edges)
    resolved: dict[str, dict[int, float]] = {}
    for node_id, weights in distributions.items():
        own_taxon, own_prob = argmax_taxon(weights)
        links = neighbours.get(node_id, [])
        if own_prob >= confident or not links:
            resolved[node_id] = dict(weights)
            continue
        vote = normalize(_weighted_vote(links, distributions))
        vote_taxon, vote_prob = argmax_taxon(vote)
        if vote_taxon == 0 or vote_prob < vote_min:
            resolved[node_id] = dict(weights)
            continue
        own_genus = taxonomy.ancestor_at_rank(own_taxon, "G")
        vote_genus = taxonomy.ancestor_at_rank(vote_taxon, "G")
        if own_genus is None or vote_genus is None or own_genus != vote_genus:
            resolved[node_id] = dict(weights)
            continue
        resolved[node_id] = vote
    return resolved


def _neighbours(edges: list[dict]) -> dict[str, list[tuple[str, float]]]:
    """Outgoing targets with positive weight. A missing weight defaults to 1."""
    neighbours: dict[str, list[tuple[str, float]]] = {}
    for edge in edges:
        weight = float(edge.get("weight", 1.0))
        if weight <= 0:
            continue
        neighbours.setdefault(edge["source"], []).append((edge["target"], weight))
    return neighbours


def _weighted_vote(
    links: list[tuple[str, float]],
    distributions: dict[str, dict[int, float]],
) -> dict[int, float]:
    """Weight-averaged neighbour distribution. Empty neighbours yield an empty map."""
    total = 0.0
    accumulated: dict[int, float] = {}
    for other, weight in links:
        other_weights = distributions.get(other, {})
        if not other_weights:
            continue
        total += weight
        for taxon_id, probability in other_weights.items():
            accumulated[taxon_id] = accumulated.get(taxon_id, 0.0) + weight * probability
    if total <= 0:
        return {}
    return {taxon_id: value / total for taxon_id, value in accumulated.items()}
