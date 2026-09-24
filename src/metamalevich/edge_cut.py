"""Drop confident disagreeing edges, then leak labels on what remains.

Some kNN edges join nodes whose true species differ. Leakage across those
edges pulls a correct call toward the neighbour. This resolver removes an
edge when both ends are already confident and disagree, then runs decaying
leakage on the kept edges.
"""

from __future__ import annotations

from metamalevich.leakage import decaying_leakage
from metamalevich.resolve import argmax_taxon


def cut_then_leak(
    distributions: dict[str, dict[int, float]],
    edges: list[dict],
    taxonomy=None,
    *,
    confident: float = 0.8,
    decay: float = 0.5,
    iterations: int = 4,
) -> dict[str, dict[int, float]]:
    """Drop confident disagreements, then spread labels with decaying leakage.

    An edge is dropped when both endpoints are present in ``distributions``,
    both argmax probabilities are at least ``confident``, and the argmax taxa
    differ. Ties break toward the smaller taxon id. An endpoint that is
    missing, or either end below ``confident``, keeps the edge. ``taxonomy``
    is accepted so callers can pass it and is not used.

    The kept edges are passed to ``decaying_leakage``. ``confident`` and
    ``decay`` must lie in ``[0, 1]`` and ``iterations`` must be at least 1.
    ``distributions`` and ``edges`` are not modified.
    """
    del taxonomy
    if not 0.0 <= confident <= 1.0:
        raise ValueError("confident must be in [0, 1]")
    if not 0.0 <= decay <= 1.0:
        raise ValueError("decay must be in [0, 1]")
    if iterations < 1:
        raise ValueError("iterations must be >= 1")

    kept = [edge for edge in edges if not _confident_disagreement(edge, distributions, confident)]
    return decaying_leakage(distributions, kept, decay=decay, iterations=iterations)


def _confident_disagreement(
    edge: dict,
    distributions: dict[str, dict[int, float]],
    confident: float,
) -> bool:
    """True when both ends exist, are confident, and their argmax taxa differ."""
    source = edge.get("source")
    target = edge.get("target")
    if source not in distributions or target not in distributions:
        return False
    source_taxon, source_probability = argmax_taxon(distributions[source])
    target_taxon, target_probability = argmax_taxon(distributions[target])
    return (
        source_probability >= confident
        and target_probability >= confident
        and source_taxon != target_taxon
    )
