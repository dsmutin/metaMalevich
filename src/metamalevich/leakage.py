"""Decaying label leakage along graph neighbours.

Each step mixes a node with its neighbours, then restarts on that node's own
evidence. The restart term is what makes leaked mass decay with hop distance.
The same update applies to any graph of this shape, whether the edges are
contig links or a later flipped edge graph.
"""

from __future__ import annotations

from metamalevich.aggregate import normalize


def decaying_leakage(
    distributions: dict[str, dict[int, float]],
    edges: list[dict],
    *,
    decay: float = 0.5,
    iterations: int = 4,
) -> dict[str, dict[int, float]]:
    """Spread label mass along neighbours without replacing a node's own evidence.

    ``evidence`` is a renormalized copy of ``distributions`` and stays fixed.
    ``current`` starts as that evidence. Each iteration replaces every node with
    ``normalize((1 - decay) * evidence + decay * neighbour_vote)``, where
    ``neighbour_vote`` is the weight-normalized mixture of the current
    distributions on positive-weight neighbours. A node with no such neighbours
    keeps its evidence. The restart term is what makes leaked mass decay with
    hop distance: influence from a neighbour is scaled by ``decay`` at every hop
    and cannot erase the original labels.

    ``decay`` must lie in ``[0, 1]`` and ``iterations`` must be at least 1.
    ``decay`` of 0 leaves the renormalized evidence unchanged. Input maps are
    not modified. Every input node is present in the result.
    """
    if not 0.0 <= decay <= 1.0:
        raise ValueError("decay must be in [0, 1]")
    if iterations < 1:
        raise ValueError("iterations must be >= 1")

    evidence = {node: normalize(weights) for node, weights in distributions.items()}
    current = {node: dict(weights) for node, weights in evidence.items()}
    neighbours = _neighbours(edges)
    for _ in range(iterations):
        updated: dict[str, dict[int, float]] = {}
        for node, evidence_weights in evidence.items():
            if node not in neighbours:
                updated[node] = dict(evidence_weights)
                continue
            vote = _neighbour_vote(neighbours[node], current)
            if not vote:
                updated[node] = dict(evidence_weights)
                continue
            mixed = {
                taxon_id: (1.0 - decay) * evidence_weights.get(taxon_id, 0.0) + decay * vote.get(taxon_id, 0.0)
                for taxon_id in set(evidence_weights) | set(vote)
            }
            updated[node] = normalize(mixed)
        current = updated
    return current


def _neighbours(edges: list[dict]) -> dict[str, list[tuple[str, float]]]:
    """Outgoing targets with positive weight. A missing weight defaults to 1."""
    neighbours: dict[str, list[tuple[str, float]]] = {}
    for edge in edges:
        weight = float(edge.get("weight", 1.0))
        if weight <= 0:
            continue
        neighbours.setdefault(edge["source"], []).append((edge["target"], weight))
    return neighbours


def _neighbour_vote(
    links: list[tuple[str, float]],
    current: dict[str, dict[int, float]],
) -> dict[int, float]:
    """Weight-normalized mean of the current distributions on ``links``."""
    total = 0.0
    accumulated: dict[int, float] = {}
    for other, weight in links:
        other_weights = current.get(other)
        if not other_weights:
            continue
        total += weight
        for taxon_id, probability in other_weights.items():
            accumulated[taxon_id] = accumulated.get(taxon_id, 0.0) + weight * probability
    if total <= 0:
        return {}
    return {taxon_id: value / total for taxon_id, value in accumulated.items() if value > 0}
