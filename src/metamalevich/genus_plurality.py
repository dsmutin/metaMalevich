"""Fill uncertain nodes with the plurality species of the same genus.

Edges are ignored. This is a genus prior, not a graph model: the label comes
from confident nodes that share a genus, not from graph neighbours.
"""

from __future__ import annotations

from metamalevich.resolve import argmax_taxon
from metamalevich.taxonomy import Taxonomy


def genus_plurality(
    distributions: dict[str, dict[int, float]],
    edges: list[dict],
    taxonomy: Taxonomy,
    calls: dict[str, int] | None = None,
    lengths: dict[str, int] | None = None,
    *,
    confident: float = 0.9,
    uncertain: float = 0.5,
) -> dict[str, dict[int, float]]:
    """Rewrite uncertain nodes to the plurality species of their genus.

    Edges are ignored. This is a genus prior, not a graph model.

    The genus of a node is ``taxonomy.ancestor_at_rank(argmax, "G")`` when
    the argmax taxon is not 0. Ties in the argmax break toward the smaller
    taxon id. When the argmax is 0 and ``calls`` contains that node, the
    genus is ``ancestor_at_rank(calls[node], "G")``, which may itself be a
    genus. A node with no genus is copied unchanged.

    A donor has an argmax other than 0 whose probability is at least
    ``confident``. The plurality species of a genus is the donor argmax with
    the largest summed length. A missing length counts as 1. An equal sum
    takes the smaller taxon id.

    A node whose genus has a donor is replaced by a one-hot of that
    plurality when its argmax probability is below ``uncertain`` or its
    argmax is 0. Every other node is copied. Inputs are not modified.

    ``confident`` and ``uncertain`` must each lie in ``[0, 1]``.
    """
    del edges
    _require_unit_interval("confident", confident)
    _require_unit_interval("uncertain", uncertain)
    _require_lengths(lengths)
    plurality = _plurality_by_genus(distributions, taxonomy, lengths, confident)
    resolved: dict[str, dict[int, float]] = {}
    for node_id, weights in distributions.items():
        taxon_id, probability = argmax_taxon(weights)
        genus_id = _genus_of(node_id, taxon_id, taxonomy, calls)
        species_id = plurality[genus_id] if genus_id is not None and genus_id in plurality else None
        if species_id is not None and (probability < uncertain or taxon_id == 0):
            resolved[node_id] = {species_id: 1.0}
        else:
            resolved[node_id] = dict(weights)
    return resolved


def _require_unit_interval(name: str, value: float) -> None:
    """Raise when ``value`` is not a real number in ``[0, 1]``."""
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not 0.0 <= float(value) <= 1.0:
        raise ValueError(f"{name} must be in [0, 1]")


def _require_lengths(lengths: dict[str, int] | None) -> None:
    """Raise when a provided length is not a positive int."""
    if lengths is None:
        return
    for node_id, length in lengths.items():
        if isinstance(length, bool) or not isinstance(length, int) or length < 1:
            raise ValueError(f"length for node {node_id} must be a positive int")


def _length(node_id: str, lengths: dict[str, int] | None) -> int:
    """Node length used in the plurality sum. Missing counts as 1."""
    if lengths is None or node_id not in lengths:
        return 1
    return lengths[node_id]


def _genus_of(
    node_id: str,
    taxon_id: int,
    taxonomy: Taxonomy,
    calls: dict[str, int] | None,
) -> int | None:
    """Genus ancestor of the argmax, or of the raw call when the argmax is 0."""
    if taxon_id != 0:
        return taxonomy.ancestor_at_rank(taxon_id, "G")
    if calls is None or node_id not in calls:
        return None
    return taxonomy.ancestor_at_rank(calls[node_id], "G")


def _plurality_by_genus(
    distributions: dict[str, dict[int, float]],
    taxonomy: Taxonomy,
    lengths: dict[str, int] | None,
    confident: float,
) -> dict[int, int]:
    """Map each genus that has a donor to its plurality species."""
    sums: dict[int, dict[int, int]] = {}
    for node_id, weights in distributions.items():
        taxon_id, probability = argmax_taxon(weights)
        if taxon_id == 0 or probability < confident:
            continue
        genus_id = taxonomy.ancestor_at_rank(taxon_id, "G")
        if genus_id is None:
            continue
        genus_sums = sums.setdefault(genus_id, {})
        genus_sums[taxon_id] = genus_sums.get(taxon_id, 0) + _length(node_id, lengths)
    plurality: dict[int, int] = {}
    for genus_id, species_lengths in sums.items():
        plurality[genus_id] = min(
            species_lengths,
            key=lambda species_id: (-species_lengths[species_id], species_id),
        )
    return plurality
