"""Keep two close same-genus labels and drop every other taxon.

The argmax is preserved: the leading taxon stays the leader after the mixture.
``edges`` are ignored.
"""

from __future__ import annotations

from metamalevich.aggregate import normalize
from metamalevich.taxonomy import Taxonomy


def margin_mixture(
    distributions: dict[str, dict[int, float]],
    edges: list[dict],
    taxonomy: Taxonomy,
    *,
    margin: float = 0.4,
) -> dict[str, dict[int, float]]:
    """Keep a close same-genus pair, or copy the node distribution unchanged.

    ``edges`` are ignored. For each node with at least two taxa of positive
    mass, the top two are ordered by the argmax key: higher probability, then
    smaller taxon id. When those two share a genus, taxon 0 is not either of
    them, and the probability gap is strictly below ``margin``, the result is
    the renormalized pair and every other taxon is dropped. Otherwise the full
    distribution is copied. The argmax is preserved because the leader is
    always kept and its weight stays at least as large as the runner-up.
    ``margin`` must lie in ``[0, 1]``. Inputs are not modified.
    """
    del edges
    if not 0.0 <= margin <= 1.0:
        raise ValueError("margin must be in [0, 1]")
    resolved: dict[str, dict[int, float]] = {}
    for node_id, weights in distributions.items():
        resolved[node_id] = _mixture_node(weights, taxonomy, margin)
    return resolved


def _mixture_node(
    weights: dict[int, float],
    taxonomy: Taxonomy,
    margin: float,
) -> dict[int, float]:
    """Return the two-taxon mixture, or a copy of ``weights``."""
    ranked = sorted(
        ((int(taxon_id), float(probability)) for taxon_id, probability in weights.items() if float(probability) > 0),
        key=lambda item: (-item[1], item[0]),
    )
    if len(ranked) < 2:
        return dict(weights)
    (first, first_p), (second, second_p) = ranked[0], ranked[1]
    if first == 0 or second == 0:
        return dict(weights)
    genus = taxonomy.ancestor_at_rank(first, "G")
    if genus is None or genus != taxonomy.ancestor_at_rank(second, "G"):
        return dict(weights)
    if first_p - second_p >= margin:
        return dict(weights)
    return normalize({first: first_p, second: second_p})
