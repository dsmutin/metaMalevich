"""Graph-aware resolvers.

Each resolver returns a full distribution per node. Hard labels are only
produced by methods whose name says so. Neighbour and edge methods read the
supplied distributions and do not edit the input maps.
"""

from __future__ import annotations

from metamalevich.aggregate import normalize


RESOLVER_PARAMETERS = {
    "initial_colouring": {"source": "kraken2 classified taxid", "rollup": "S", "assignment": "hard"},
    "probability_sum": {"aggregation": "probability_sum", "rollup": "S", "graph": False},
    "lca": {"aggregation": "lca", "min_fraction": 0.05, "rollup": "S"},
    "gated_neighbour": {
        "iterations": 3,
        "confident": 0.8,
        "consensus": 0.75,
        "uncertain_mix": 0.7,
        "conflict_mix": 0.5,
    },
    "bayesian_edge": {
        "iterations": 3,
        "neighbour_weight": 1.0,
        "edge_weight": 1.0,
        "floor": 1e-6,
        "import_min_probability": 0.75,
    },
}


def resolver_parameters(method: str) -> dict:
    """Copy of the parameters recorded for one hypothesis."""
    return dict(RESOLVER_PARAMETERS.get(method, {}))


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
    iterations: int | None = None,
    confident: float | None = None,
    consensus: float | None = None,
    uncertain_mix: float | None = None,
    conflict_mix: float | None = None,
) -> dict[str, dict[int, float]]:
    """Mix a node with its neighbours when it is uncertain or conflicts with them.

    A node whose maximum probability is at least ``confident`` keeps its own
    distribution unless the neighbour vote reaches ``consensus`` on a different
    taxon. Uncertain nodes take ``uncertain_mix`` of the neighbour vote.
    """
    cfg = RESOLVER_PARAMETERS["gated_neighbour"]
    iterations = cfg["iterations"] if iterations is None else iterations
    confident = cfg["confident"] if confident is None else confident
    consensus = cfg["consensus"] if consensus is None else consensus
    uncertain_mix = cfg["uncertain_mix"] if uncertain_mix is None else uncertain_mix
    conflict_mix = cfg["conflict_mix"] if conflict_mix is None else conflict_mix
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
    iterations: int | None = None,
    neighbour_weight: float | None = None,
    edge_weight: float | None = None,
    floor: float | None = None,
    import_min_probability: float | None = None,
) -> dict[str, dict[int, float]]:
    """Posterior proportional to evidence, neighbour vote, and incident edge colours.

    ``log p(t) = log evidence + neighbour_weight * log vote + edge_weight * log edge``.
    A taxon absent from the node's own evidence is added only when the neighbour
    vote or the edge view reaches ``import_min_probability`` on that taxon.
    Edge colours are computed once from the original node evidence and are not
    updated between iterations. Each edge colour is the minimum of its endpoints,
    so multiplying by it counts endpoint evidence a second time.
    Missing taxa use ``floor`` so a zero in one view does not erase the others
    before renormalization. The returned map is still a distribution.
    """
    cfg = RESOLVER_PARAMETERS["bayesian_edge"]
    iterations = cfg["iterations"] if iterations is None else iterations
    neighbour_weight = cfg["neighbour_weight"] if neighbour_weight is None else neighbour_weight
    edge_weight = cfg["edge_weight"] if edge_weight is None else edge_weight
    floor = cfg["floor"] if floor is None else floor
    import_min_probability = (
        cfg["import_min_probability"] if import_min_probability is None else import_min_probability
    )
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
            keys = set(own)
            if vote:
                vote_taxon, vote_prob = argmax_taxon(vote)
                if vote_prob >= import_min_probability and vote_taxon not in keys:
                    keys.add(vote_taxon)
            if edge_view:
                edge_taxon, edge_prob = argmax_taxon(edge_view)
                if edge_prob >= import_min_probability and edge_taxon not in keys:
                    keys.add(edge_taxon)
            if not keys:
                updated[element_id] = {}
                continue
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
