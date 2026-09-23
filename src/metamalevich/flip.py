"""Flip a contig graph into its line graph and project labels back.

Each original edge becomes a node. That node's taxon distribution is the
renormalized sum of the distributions on its two endpoints. Flipped edges
join original edges that share a vertex. Inference runs on the line graph
and is projected back onto the original contigs.
"""

from __future__ import annotations

from collections.abc import Callable

from metamalevich.aggregate import normalize

Distribution = dict[int, float]
NodeDistributions = dict[str, Distribution]
Edge = dict


def flip_graph(
    edges: list[Edge],
    node_distributions: NodeDistributions,
) -> tuple[NodeDistributions, list[Edge]]:
    """Turn each contig edge into a node and link edges that share a vertex.

    The new node's distribution is the renormalized sum of the two endpoint
    distributions: every taxon present on either neighbour is kept. Non-positive
    taxon weights are dropped. Taxon 0 is kept. Two empty endpoints yield an
    empty distribution.

    Each unordered pair of original edges that share an endpoint becomes one
    flipped edge. Pairs are ordered by ``source < target`` using string order,
    so the edge list does not depend on input order. The new ``edge_id`` is
    ``line:{left}:{right}``, the orientation is ``+``, and the weight is the
    minimum of the two original edge weights. A missing weight counts as 1.
    A self-loop is not emitted.
    """
    for edge in edges:
        _require_edge(edge)
    flipped_distributions = {
        edge["edge_id"]: _sum_endpoints(edge, node_distributions) for edge in edges
    }
    return flipped_distributions, _line_edges(edges)


def project_to_nodes(
    edges: list[Edge],
    flipped_distributions: NodeDistributions,
    fallback: NodeDistributions,
) -> NodeDistributions:
    """Sum incident flipped-node distributions back onto each original node.

    The sum is renormalized. A node with no incident flipped mass keeps a copy
    of ``fallback[node]``. Every key in ``fallback`` is returned. Inputs are
    not modified.
    """
    for edge in edges:
        _require_edge(edge)
    incident = _incident_edge_ids(edges)
    projected: NodeDistributions = {}
    for node, prior in fallback.items():
        acc: Distribution = {}
        for edge_id in incident.get(node, []):
            _add_positive(acc, flipped_distributions.get(edge_id, {}))
        if sum(acc.values()) <= 0:
            projected[node] = dict(prior)
        else:
            projected[node] = normalize(acc)
    return projected


def infer_on_flipped_graph(
    edges: list[Edge],
    node_distributions: NodeDistributions,
    resolver: Callable[[NodeDistributions, list[Edge]], NodeDistributions],
) -> NodeDistributions:
    """Run ``resolver`` once on the flipped graph and project onto contigs.

    ``resolver`` receives the flipped distributions and flipped edges. It must
    return a distribution for every flipped node id. Those distributions are
    projected onto the original nodes. ``node_distributions`` is the fallback
    where a node has no incident flipped mass.
    """
    flipped_distributions, flipped_edges = flip_graph(edges, node_distributions)
    resolved = resolver(flipped_distributions, flipped_edges)
    if not isinstance(resolved, dict):
        raise TypeError("resolver must return a dict of distributions keyed by flipped node id")
    missing = [edge_id for edge_id in flipped_distributions if edge_id not in resolved]
    if missing:
        listed = ", ".join(str(edge_id) for edge_id in missing)
        raise ValueError(f"resolver omitted flipped nodes: {listed}")
    return project_to_nodes(edges, resolved, node_distributions)


def _require_edge(edge: Edge) -> None:
    """Fail if an edge lacks ``edge_id``, ``source``, or ``target``."""
    missing = [key for key in ("edge_id", "source", "target") if key not in edge]
    if missing:
        listed = ", ".join(missing)
        raise ValueError(f"edge is missing required keys: {listed}")


def _edge_weight(edge: Edge) -> float:
    """Return the edge weight, or 1 when the key is absent."""
    if "weight" not in edge:
        return 1.0
    return float(edge["weight"])


def _add_positive(acc: Distribution, distribution: Distribution) -> None:
    """Add strictly positive taxon weights from ``distribution`` into ``acc``."""
    for taxon_id, value in distribution.items():
        weight = float(value)
        if weight <= 0:
            continue
        key = int(taxon_id)
        acc[key] = acc.get(key, 0.0) + weight


def _sum_endpoints(edge: Edge, node_distributions: NodeDistributions) -> Distribution:
    """Renormalize the sum of the distributions on the two endpoints."""
    acc: Distribution = {}
    _add_positive(acc, node_distributions.get(edge["source"], {}))
    _add_positive(acc, node_distributions.get(edge["target"], {}))
    return normalize(acc)


def _incident_edge_ids(edges: list[Edge]) -> dict[str, list[str]]:
    """Map each vertex to the original edge ids that touch it, once each."""
    incident: dict[str, list[str]] = {}
    for edge in edges:
        edge_id = edge["edge_id"]
        for endpoint in (edge["source"], edge["target"]):
            bucket = incident.setdefault(endpoint, [])
            if edge_id not in bucket:
                bucket.append(edge_id)
    return incident


def _line_edges(edges: list[Edge]) -> list[Edge]:
    """Build deterministic line-graph edges for pairs that share a vertex."""
    weights: dict[str, float] = {}
    for edge in edges:
        weights[edge["edge_id"]] = _edge_weight(edge)
    pairs: set[tuple[str, str]] = set()
    for edge_ids in _incident_edge_ids(edges).values():
        for index, left_id in enumerate(edge_ids):
            for right_id in edge_ids[index + 1 :]:
                if str(left_id) == str(right_id):
                    continue
                if str(left_id) < str(right_id):
                    pair = (left_id, right_id)
                else:
                    pair = (right_id, left_id)
                pairs.add(pair)
    flipped: list[Edge] = []
    for left_id, right_id in sorted(pairs, key=lambda item: (str(item[0]), str(item[1]))):
        flipped.append(
            {
                "edge_id": f"line:{left_id}:{right_id}",
                "source": left_id,
                "target": right_id,
                "orientation": "+",
                "weight": min(weights[left_id], weights[right_id]),
            }
        )
    return flipped
