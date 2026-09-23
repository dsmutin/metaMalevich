"""Tiny closed community used by the mandatory tests.

Node n2 is classified as species B, while both neighbours are species A and
the truth label of n2 is A. Gated neighbour resolution must lower the L1
abundance error relative to that hard colouring.
"""

from __future__ import annotations

from metamalevich.aggregate import normalize
from metamalevich.evaluate import abundance_scores, classification_scores
from metamalevich.reprofile import profile_nodes
from metamalevich.resolve import resolve
from metamalevich.taxonomy import parse_kraken_report

REPORT = """\
100.00\t3\t0\tR\t1\troot
100.00\t3\t0\tD\t2\t  Bacteria
66.67\t2\t2\tS\t10\t    speciesA
33.33\t1\t1\tS\t20\t    speciesB
"""


def run_toy() -> dict:
    """Return contract keys plus the toy L1 comparison."""
    taxonomy = parse_kraken_report(REPORT, source="toy", version="toy-report")
    evidence = {
        "n1": normalize({10: 10}),
        "n2": normalize({20: 6, 10: 4}),
        "n3": normalize({10: 10}),
    }
    edges = [
        {"edge_id": "e1", "source": "n1", "target": "n2", "orientation": "++", "weight": 1.0},
        {"edge_id": "e2", "source": "n2", "target": "n1", "orientation": "++", "weight": 1.0},
        {"edge_id": "e3", "source": "n2", "target": "n3", "orientation": "++", "weight": 1.0},
        {"edge_id": "e4", "source": "n3", "target": "n2", "orientation": "++", "weight": 1.0},
    ]
    lengths = {"n1": 100, "n2": 100, "n3": 100}
    truth = {"n1": 10, "n2": 10, "n3": 10}
    abundance = {10: 1.0}
    initial = {"n1": {10: 1.0}, "n2": {20: 1.0}, "n3": {10: 1.0}}
    resolved = resolve(evidence, edges, "gated_neighbour")
    initial_l1 = _l1(initial, lengths, taxonomy, abundance)
    resolved_l1 = _l1(resolved, lengths, taxonomy, abundance)
    scores = classification_scores(resolved, truth)
    return {
        "status": "ok",
        "ok": True,
        "input_path": None,
        "method": "gated_neighbour",
        "initial_l1": initial_l1,
        "resolved_l1": resolved_l1,
        "beats_initial": resolved_l1 < initial_l1,
        "accuracy": scores["accuracy"],
    }


def _l1(distributions, lengths, taxonomy, abundance) -> float:
    rows, summary = profile_nodes(distributions, lengths, taxonomy)
    predicted = {int(row["taxon_id"]): float(row["relative_abundance"]) for row in rows}
    fraction = float(summary["unclassified_fraction"])
    if fraction > 0:
        predicted[0] = predicted.get(0, 0.0) + fraction
    return float(abundance_scores(predicted, abundance)["l1"])
