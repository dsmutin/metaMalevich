"""Length-weighted taxonomic profile from node posteriors.

Shared nodes contribute to every taxon in proportion to the posterior.
The profile does not replace that posterior with a hard label unless the
caller already passed a hard distribution.
"""

from __future__ import annotations

from metamalevich.resolve import argmax_taxon
from metamalevich.taxonomy import Taxonomy


def _partition(distribution: dict[int, float], *, shared_floor: float = 0.1, ambiguous_ceiling: float = 0.5) -> str:
    if not distribution:
        return "unclassified"
    _taxon, top = argmax_taxon(distribution)
    ordered = sorted(distribution.values(), reverse=True)
    second = ordered[1] if len(ordered) > 1 else 0.0
    if top < ambiguous_ceiling:
        return "ambiguous"
    if second >= shared_floor:
        return "shared"
    return "unique"


def profile_nodes(
    distributions: dict[str, dict[int, float]],
    lengths: dict[str, int],
    taxonomy: Taxonomy,
) -> tuple[list[dict], dict]:
    """Build the abundance table and a small unclassified summary.

    ``estimated_abundance`` is assigned bases. ``relative_abundance`` divides
    by the sum of node lengths, so unclassified bases remain in the total.
    """
    bases: dict[int, float] = {}
    unique: dict[int, float] = {}
    shared: dict[int, float] = {}
    ambiguous: dict[int, float] = {}
    confidence_mass: dict[int, float] = {}
    confidence_bases: dict[int, float] = {}
    total_length = 0
    unclassified_bases = 0.0
    for node_id, distribution in distributions.items():
        length = int(lengths.get(node_id, 0))
        if length < 0:
            raise ValueError(f"negative length for {node_id}")
        total_length += length
        kind = _partition(distribution)
        if kind == "unclassified" or not distribution:
            unclassified_bases += length
            continue
        assigned = 0.0
        for taxon_id, prob in distribution.items():
            if prob <= 0 or taxon_id == 0:
                continue
            portion = prob * length
            assigned += prob
            bases[taxon_id] = bases.get(taxon_id, 0.0) + portion
            confidence_mass[taxon_id] = confidence_mass.get(taxon_id, 0.0) + portion
            confidence_bases[taxon_id] = confidence_bases.get(taxon_id, 0.0) + length
            if kind == "unique":
                unique[taxon_id] = unique.get(taxon_id, 0.0) + portion
            elif kind == "shared":
                shared[taxon_id] = shared.get(taxon_id, 0.0) + portion
            else:
                ambiguous[taxon_id] = ambiguous.get(taxon_id, 0.0) + portion
        unclassified_bases += (1.0 - assigned) * length
    denominator = float(total_length) if total_length else 1.0
    rows = []
    for taxon_id in sorted(bases):
        assigned_bases = bases[taxon_id]
        conf_den = confidence_bases.get(taxon_id, 0.0)
        rows.append(
            {
                "taxon_id": taxon_id,
                "rank": taxonomy.rank(taxon_id),
                "name": taxonomy.name(taxon_id),
                "estimated_abundance": assigned_bases,
                "relative_abundance": assigned_bases / denominator,
                "assigned_bases": assigned_bases,
                "assigned_reads": 0,
                "unique_bases": unique.get(taxon_id, 0.0),
                "shared_bases": shared.get(taxon_id, 0.0),
                "ambiguous_bases": ambiguous.get(taxon_id, 0.0),
                "confidence": (confidence_mass[taxon_id] / conf_den) if conf_den else 0.0,
            }
        )
    summary = {
        "total_bases": total_length,
        "unclassified_bases": unclassified_bases,
        "unclassified_reads": 0,
        "unclassified_fraction": unclassified_bases / denominator if total_length else 0.0,
        "read_weight": "unused",
        "weight": "node_length",
    }
    return rows, summary
