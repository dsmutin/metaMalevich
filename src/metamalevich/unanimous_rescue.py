"""Rescue an uncertain node when every neighbour names one taxon.

Gated smoothing can overwrite a correct node. This resolver copies every node
and replaces one only when that node is uncertain and every positive-weight
outgoing neighbour already carries the same taxon in a compatible genus.
"""

from __future__ import annotations

from metamalevich.aggregate import normalize
from metamalevich.resolve import argmax_taxon
from metamalevich.taxonomy import Taxonomy


def unanimous_rescue(
    distributions: dict[str, dict[int, float]],
    edges: list[dict],
    taxonomy: Taxonomy,
    *,
    uncertain: float = 0.55,
) -> dict[str, dict[int, float]]:
    """Copy each node, adopting a neighbour taxon only when the vote is unanimous.

    A node is kept as a copy when its argmax probability is at least
    ``uncertain``, when it has no positive-weight outgoing neighbour, when
    those neighbours do not share one argmax taxon, or when that taxon is 0.
    Otherwise the shared neighbour taxon is adopted only if
    ``ancestor_at_rank(..., "G")`` is the same non-null genus as the node's
    argmax, or the node's argmax is 0 and the neighbour genus is not None.
    The replacement is ``normalize(0.5 * own + 0.5 * one-hot neighbour taxon)``.

    A missing edge weight counts as 1. ``uncertain`` must lie in ``[0, 1]``.
    Inputs are not modified.
    """
    if not 0.0 <= uncertain <= 1.0:
        raise ValueError("uncertain must be in [0, 1]")
    outgoing = _outgoing(edges)
    resolved: dict[str, dict[int, float]] = {}
    for node, weights in distributions.items():
        own_taxon, own_probability = argmax_taxon(weights)
        links = outgoing.get(node, [])
        if own_probability >= uncertain or not links:
            resolved[node] = dict(weights)
            continue
        agreed = _shared_argmax(links, distributions)
        if agreed is None or agreed == 0 or not _genus_allows(taxonomy, own_taxon, agreed):
            resolved[node] = dict(weights)
            continue
        resolved[node] = _half_with_one_hot(weights, agreed)
    return resolved


def _outgoing(edges: list[dict]) -> dict[str, list[tuple[str, float]]]:
    """Positive-weight outgoing neighbours. A missing weight counts as 1."""
    neighbours: dict[str, list[tuple[str, float]]] = {}
    for edge in edges:
        weight = float(edge.get("weight", 1.0))
        if weight <= 0:
            continue
        neighbours.setdefault(edge["source"], []).append((edge["target"], weight))
    return neighbours


def _shared_argmax(
    links: list[tuple[str, float]],
    distributions: dict[str, dict[int, float]],
) -> int | None:
    """Return the shared neighbour argmax, or None when neighbours differ."""
    agreed: int | None = None
    for other, _weight in links:
        taxon, _probability = argmax_taxon(distributions.get(other, {}))
        if agreed is None:
            agreed = taxon
        elif taxon != agreed:
            return None
    return agreed


def _genus_allows(taxonomy: Taxonomy, own_taxon: int, neighbour_taxon: int) -> bool:
    """True when the neighbour taxon is in the node's genus, or the node is unclassified."""
    neighbour_genus = taxonomy.ancestor_at_rank(neighbour_taxon, "G")
    if own_taxon == 0:
        return neighbour_genus is not None
    own_genus = taxonomy.ancestor_at_rank(own_taxon, "G")
    return own_genus is not None and own_genus == neighbour_genus


def _half_with_one_hot(weights: dict[int, float], taxon_id: int) -> dict[int, float]:
    """Return ``normalize(0.5 * own + 0.5 * one-hot taxon)``."""
    mixed = {key: 0.5 * value for key, value in weights.items()}
    mixed[taxon_id] = mixed.get(taxon_id, 0.0) + 0.5
    return normalize(mixed)
