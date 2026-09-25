"""Gate 4-mer genus or family transfer by cosine, then by the assembly graph.

Unconstrained nearest-neighbour transfer labels every contig. This script
keeps that label only when the cosine is at least a threshold. A second rule
uses the MEGAHIT neighbour genus when 80% of neighbour length agrees, and
otherwise the 4-mer label at the same threshold.
"""

from __future__ import annotations

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "examples" / "heldout_genera"))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "examples"))
from ncbi_taxonomy import ncbi_parser  # noqa: E402

from community import assembly_graph_from_fastg  # noqa: E402
from metamalevich.composition_genus import canonical_4mer  # noqa: E402
from metamalevich.tables import write_tsv  # noqa: E402

THRESHOLDS = (0.0, 0.5, 0.7, 0.85, 0.95)


def _parser():
    return ncbi_parser()


def _rank(parser, taxon_id: int, rank: str) -> int:
    if taxon_id <= 0:
        return 0
    return parser.get_ancestor_by_rank(taxon_id, rank) or 0


def _fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    name = ""
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            if name:
                sequences[name] = "".join(chunks)
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if name:
        sequences[name] = "".join(chunks)
    return sequences


def _best(sequences: dict[str, str], labels: dict[str, int]) -> dict[str, tuple[int, float]]:
    donors = [node_id for node_id, label in labels.items() if label > 0 and node_id in sequences]
    queries = [node_id for node_id, label in labels.items() if label <= 0 and node_id in sequences]
    donor_matrix = np.vstack([canonical_4mer(sequences[node_id]) for node_id in donors])
    donor_labels = [labels[node_id] for node_id in donors]
    found: dict[str, tuple[int, float]] = {}
    for start in range(0, len(queries), 2000):
        batch = queries[start : start + 2000]
        query_matrix = np.vstack([canonical_4mer(sequences[node_id]) for node_id in batch])
        scores = query_matrix @ donor_matrix.T
        best = scores.argmax(axis=1)
        for index, node_id in enumerate(batch):
            if float(np.linalg.norm(query_matrix[index])) == 0.0:
                continue
            found[node_id] = (donor_labels[int(best[index])], float(scores[index, int(best[index])]))
    return found


def _neighbours(labels: dict[str, int], edges: list[dict], lengths: dict[str, int]) -> dict[str, int]:
    neighbours: dict[str, set[str]] = defaultdict(set)
    for edge in edges:
        neighbours[edge["source"]].add(edge["target"])
        neighbours[edge["target"]].add(edge["source"])
    chosen = {}
    for node_id, label in labels.items():
        if label:
            continue
        votes: Counter = Counter()
        for other in neighbours.get(node_id, ()):
            other_label = labels.get(other, 0)
            if other_label:
                votes[other_label] += lengths.get(other, 1)
        if not votes:
            continue
        best, mass = votes.most_common(1)[0]
        if mass / sum(votes.values()) >= 0.8:
            chosen[node_id] = best
    return chosen


def _placement(path: Path) -> dict[str, str]:
    best: dict[str, tuple[int, str]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split("\t")
        if len(parts) < 10:
            continue
        read_id = parts[0].split("/")[0]
        matches = int(parts[9])
        previous = best.get(read_id)
        if previous is None or matches > previous[0]:
            best[read_id] = (matches, parts[5])
    return {read_id: target for read_id, (_matches, target) in best.items()}


def _classifier(path: Path) -> dict[str, int]:
    calls = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split("\t")
        calls[parts[1].split("/")[0]] = 0 if parts[0] == "U" else int(parts[2])
    return calls


def _truth(example: Path, column: str) -> dict[int, float]:
    accessions = {
        row["accession"]: int(row[column])
        for row in csv.DictReader((example / "accessions.tsv").open(encoding="utf-8"), delimiter="\t")
    }
    counts: Counter = Counter()
    for line in (example / "work" / "abundance.csv").read_text(encoding="utf-8").splitlines()[1:]:
        if not line.strip():
            continue
        accession, raw = line.split(",")
        counts[accessions[accession]] += int(raw)
    total = sum(counts.values())
    return {taxon_id: count / total for taxon_id, count in counts.items()}


def _l1(reads, placement, node_labels, parser, rank, truth) -> tuple[float, float]:
    counts: Counter = Counter()
    for read_id, taxon_id in reads.items():
        label = _rank(parser, taxon_id, rank)
        if label == 0:
            contig_id = placement.get(read_id)
            if contig_id is not None and node_labels is not None:
                label = node_labels.get(contig_id, 0)
        counts[label] += 1
    total = sum(counts.values())
    predicted = {taxon_id: count / total for taxon_id, count in counts.items()}
    keys = set(predicted) | set(truth)
    return sum(abs(predicted.get(key, 0.0) - truth.get(key, 0.0)) for key in keys), predicted.get(0, 0.0)


def score(example_name: str, rank: str, column: str, parser) -> list[dict]:
    example = ROOT / "examples" / example_name
    work = example / "work"
    edges, sequences = assembly_graph_from_fastg((work / "megahit" / "assembly.fastg").read_text(encoding="utf-8", errors="replace"))
    fasta = _fasta(work / "megahit" / "graph_nodes.fa")
    calls = {}
    for line in (work / "reprofile" / "contigs.kraken").read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split("\t")
        calls[parts[1]] = 0 if parts[0] == "U" else int(parts[2])
    labels = {node_id: _rank(parser, calls.get(node_id, 0), rank) for node_id in sequences}
    lengths = {node_id: len(sequence) for node_id, sequence in sequences.items()}
    nearest = _best(fasta, labels)
    graph_fill = _neighbours(labels, edges, lengths)
    placement = _placement(work / "reprofile" / "reads_to_contigs.paf")
    truth = _truth(example, column)
    rows = []
    for tool, path in (("kraken2", work / "classify" / "kraken2.output"), ("kaiju", work / "classify" / "kaiju.output")):
        reads = _classifier(path)
        for threshold in THRESHOLDS:
            gated = dict(labels)
            for node_id, (label, similarity) in nearest.items():
                if similarity >= threshold:
                    gated[node_id] = label
            mixed = dict(gated)
            mixed.update({node_id: label for node_id, label in graph_fill.items()})
            for method, table in ((f"{tool}_4mer_{threshold}", gated), (f"{tool}_graph_then_4mer_{threshold}", mixed)):
                l1, unclassified = _l1(reads, placement, table, parser, rank, truth)
                row = {
                    "dataset": example_name,
                    "method": method,
                    "rank": rank,
                    "l1": f"{l1:.6g}",
                    "unclassified": f"{unclassified:.6g}",
                }
                rows.append(row)
                print(f"{example_name}\t{method}\t{row['l1']}", flush=True)
    return rows


def main() -> None:
    parser = _parser()
    rows = []
    rows.extend(score("low75_enriched", "genus", "genus_taxid", parser))
    rows.extend(score("high100_enriched", "family", "family_taxid", parser))
    destination = ROOT / "examples" / "low75_enriched" / "ds" / "threshold_metrics.tsv"
    write_tsv(destination, rows, list(rows[0].keys()))
    (ROOT / "examples" / "high100_enriched" / "ds" / "threshold_metrics.tsv").write_text(
        destination.read_text(encoding="utf-8"), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
