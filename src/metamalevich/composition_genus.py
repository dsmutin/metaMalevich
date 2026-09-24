"""Pass a genus label across a canonical 4-mer graph.

Assembly-overlap graphs on the low-coverage communities are almost edgeless,
so an unclassified contig has no labelled neighbour there. Contigs from one
genome still share tetranucleotide composition. This module links each
unlabelled node to the nearest labelled node in that feature space and copies
the label. It does not invent a species.
"""

from __future__ import annotations

import numpy as np

_BASE = {"A": 0, "C": 1, "G": 2, "T": 3}
_BINS = 256
_K = 4


def canonical_4mer(sequence: str) -> np.ndarray:
    """L2-normalised canonical 4-mer frequencies.

    The coding matches ``cpp/kmer_knn.cpp``: A,C,G,T are 0..3, any other base
    breaks the k-mer, and each k-mer is stored as the minimum of the word and
    its reverse complement. An empty count vector stays the zero vector.
    """
    counts = np.zeros(_BINS, dtype=np.float64)
    code = 0
    filled = 0
    for base_char in sequence.upper():
        base = _BASE.get(base_char)
        if base is None:
            code = 0
            filled = 0
            continue
        code = ((code << 2) | base) & (_BINS - 1)
        filled += 1
        if filled < _K:
            continue
        reverse = 0
        remaining = code
        for _ in range(_K):
            reverse = (reverse << 2) | ((remaining & 3) ^ 3)
            remaining >>= 2
        counts[min(code, reverse)] += 1.0
    norm = float(np.linalg.norm(counts))
    if norm > 0.0:
        counts /= norm
    return counts


def composition_genus_graph(
    sequences: dict[str, str],
    labels: dict[str, int],
    *,
    min_sim: float = 0.0,
) -> dict[str, int]:
    """Copy the nearest positive label onto every unlabelled node.

    A label of 0 is unlabelled. Nodes absent from ``sequences`` keep their
    input label. A node whose 4-mer vector is zero, or whose best cosine is
    below ``min_sim``, stays unlabelled. Labelled nodes are not changed.
    ``min_sim`` must lie in ``[-1, 1]``. Inputs are not modified.
    """
    if min_sim < -1.0 or min_sim > 1.0:
        raise ValueError("min_sim must lie in [-1, 1]")
    result = {node_id: int(label) for node_id, label in labels.items()}
    donors = [node_id for node_id, label in result.items() if label > 0 and node_id in sequences]
    queries = [node_id for node_id in sequences if result.get(node_id, 0) <= 0]
    if not donors or not queries:
        for node_id in queries:
            result.setdefault(node_id, 0)
        return result
    donor_matrix = np.vstack([canonical_4mer(sequences[node_id]) for node_id in donors])
    donor_labels = [result[node_id] for node_id in donors]
    keep = np.linalg.norm(donor_matrix, axis=1) > 0.0
    if not bool(np.any(keep)):
        for node_id in queries:
            result.setdefault(node_id, 0)
        return result
    donor_matrix = donor_matrix[keep]
    donor_labels = [label for label, flag in zip(donor_labels, keep) if flag]
    for node_id in queries:
        query = canonical_4mer(sequences[node_id])
        if float(np.linalg.norm(query)) == 0.0:
            result[node_id] = 0
            continue
        scores = donor_matrix @ query
        best = int(np.argmax(scores))
        if float(scores[best]) < min_sim:
            result[node_id] = 0
            continue
        result[node_id] = int(donor_labels[best])
    return result
