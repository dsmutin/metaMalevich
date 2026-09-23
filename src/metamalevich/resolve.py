"""Graph-aware resolvers.

Each resolver returns a full distribution per node. Hard labels are only
produced by methods whose name says so. Neighbour and edge methods read the
supplied distributions and do not edit the input maps.
"""

from __future__ import annotations

from metamalevich.aggregate import normalize


def argmax_taxon(weights: dict[int, float]) -> tuple[int, float]:
    """Highest weight. Ties break toward the smaller taxon id."""
    if not weights:
        return (0, 0.0)
    taxon_id = min(weights, key=lambda key: (-weights[key], key))
    return taxon_id, weights[taxon_id]


def hard_assignment(distributions: dict[str, dict[int, float]]) -> dict[str, dict[int, float]]:
    """One-hot at the argmax. Tied argmax values share the one-hot mass equally."""
    resolved = {}
    for element_id, weights in distributions.items():
        if not weights:
            resolved[element_id] = {}
            continue
        best = max(weights.values())
        tied = [taxon_id for taxon_id, value in weights.items() if value == best]
        share = 1.0 / len(tied)
        resolved[element_id] = {taxon_id: share for taxon_id in tied}
    return resolved


def blend(left: dict[int, float], right: dict[int, float], alpha: float) -> dict[int, float]:
    """``(1-alpha) * left + alpha * right``, then renormalize."""
    if alpha < 0 or alpha > 1:
        raise ValueError("blend alpha must be in [0, 1]")
    keys = set(left) | set(right)
    mixed = {
        taxon_id: (1.0 - alpha) * left.get(taxon_id, 0.0) + alpha * right.get(taxon_id, 0.0)
        for taxon_id in keys
    }
    return normalize(mixed)


def _adjacency(edges: list[dict]) -> dict[str, list[tuple[str, float]]]:
    neighbours: dict[str, list[tuple[str, float]]] = {}
    for edge in edges:
        weight = float(edge.get("weight", 1.0))
        if weight <= 0:
            continue
        neighbours.setdefault(edge["source"], []).append((edge["target"], weight))
    return neighbours


def _vote(element_id: str, neighbours: dict[str, list[tuple[str, float]]], distributions: dict[str, dict[int, float]]) -> dict[int, float]:
    total = 0.0
    acc: dict[int, float] = {}
    for other, weight in neighbours.get(element_id, []):
        other_dist = distributions.get(other, {})
        if not other_dist:
            continue
        total += weight
        for taxon_id, prob in other_dist.items():
            acc[taxon_id] = acc.get(taxon_id, 0.0) + weight * prob
    if total <= 0:
        return {}
    return {taxon_id: value / total for taxon_id, value in acc.items()}


def edge_support(left: dict[int, float], right: dict[int, float]) -> dict[int, float]:
    """Shared mass on an edge: the minimum weight of each taxon, renormalized.

    An edge can therefore carry a taxon set smaller than either endpoint.
    """
    raw = {}
    for taxon_id in set(left) | set(right):
        shared = min(left.get(taxon_id, 0.0), right.get(taxon_id, 0.0))
        if shared > 0:
            raw[taxon_id] = shared
    return normalize(raw)


def edge_distributions(
    edges: list[dict],
    node_distributions: dict[str, dict[int, float]],
) -> dict[str, dict[int, float]]:
    """Colour every edge from its endpoints' evidence, independently of later resolution."""
    coloured = {}
    for edge in edges:
        coloured[edge["edge_id"]] = edge_support(
            node_distributions.get(edge["source"], {}),
            node_distributions.get(edge["target"], {}),
        )
    return coloured


def gated_neighbour(
    distributions: dict[str, dict[int, float]],
    edges: list[dict],
    *,
    iterations: int = 3,
    confident: float = 0.8,
    consensus: float = 0.75,
    uncertain_mix: float = 0.7,
    conflict_mix: float = 0.5,
) -> dict[str, dict[int, float]]:
    """Mix a node with its neighbours when it is uncertain or conflicts with them.

    A node whose maximum probability is at least ``confident`` keeps its own
    distribution unless the neighbour vote reaches ``consensus`` on a different
    taxon. Uncertain nodes take ``uncertain_mix`` of the neighbour vote.
    """
    if iterations < 1:
        raise ValueError("iterations must be >= 1")
    neighbours = _adjacency(edges)
    current = {key: normalize(value) for key, value in distributions.items()}
    for _ in range(iterations):
        updated: dict[str, dict[int, float]] = {}
        for element_id, dist in current.items():
            vote = _vote(element_id, neighbours, current)
            if not vote or not dist:
                updated[element_id] = dict(dist)
                continue
            own_taxon, own_prob = argmax_taxon(dist)
            vote_taxon, vote_prob = argmax_taxon(vote)
            if own_prob < confident:
                updated[element_id] = blend(dist, vote, uncertain_mix)
            elif vote_taxon != own_taxon and vote_prob >= consensus:
                updated[element_id] = blend(dist, vote, conflict_mix)
            else:
                updated[element_id] = dict(dist)
        current = updated
    return current


def bayesian_edge(
    distributions: dict[str, dict[int, float]],
    edges: list[dict],
    edge_dist: dict[str, dict[int, float]],
    *,
    iterations: int = 3,
    neighbour_weight: float = 1.0,
    edge_weight: float = 1.0,
    floor: float = 1e-6,
) -> dict[str, dict[int, float]]:
    """Posterior proportional to evidence, neighbour vote, and incident edge colours.

    ``log p(t) = log evidence + neighbour_weight * log vote + edge_weight * log edge``.
    Missing taxa use ``floor`` so a zero in one view does not erase the others
    before renormalization. The returned map is still a distribution.
    """
    if iterations < 1:
        raise ValueError("iterations must be >= 1")
    neighbours = _adjacency(edges)
    incident: dict[str, list[str]] = {}
    for edge in edges:
        incident.setdefault(edge["source"], []).append(edge["edge_id"])
    evidence = {key: normalize(value) for key, value in distributions.items()}
    current = {key: dict(value) for key, value in evidence.items()}
    for _ in range(iterations):
        updated: dict[str, dict[int, float]] = {}
        for element_id, _dist in current.items():
            vote = _vote(element_id, neighbours, current)
            edge_ids = incident.get(element_id, [])
            edge_acc: dict[int, float] = {}
            if edge_ids:
                for edge_id in edge_ids:
                    for taxon_id, prob in edge_dist.get(edge_id, {}).items():
                        edge_acc[taxon_id] = edge_acc.get(taxon_id, 0.0) + prob
                edge_view = normalize(edge_acc)
            else:
                edge_view = {}
            own = evidence[element_id]
            keys = set(own) | set(vote) | set(edge_view)
            raw: dict[int, float] = {}
            for taxon_id in keys:
                score = max(own.get(taxon_id, 0.0), floor)
                if vote:
                    score *= max(vote.get(taxon_id, 0.0), floor) ** neighbour_weight
                if edge_view:
                    score *= max(edge_view.get(taxon_id, 0.0), floor) ** edge_weight
                raw[taxon_id] = score
            updated[element_id] = normalize(raw)
        current = updated
    return current


def resolve(
    distributions: dict[str, dict[int, float]],
    edges: list[dict],
    method: str,
    edge_dist: dict[str, dict[int, float]] | None = None,
) -> dict[str, dict[int, float]]:
    """Dispatch a named resolver. ``probability_sum`` returns the input distributions."""
    if method == "hard":
        return hard_assignment(distributions)
    if method == "probability_sum":
        return {key: normalize(value) for key, value in distributions.items()}
    if method == "gated_neighbour":
        return gated_neighbour(distributions, edges)
    if method == "bayesian_edge":
        colours = edge_dist if edge_dist is not None else edge_distributions(edges, distributions)
        return bayesian_edge(distributions, edges, colours)
    raise ValueError(f"unknown resolver: {method}")
