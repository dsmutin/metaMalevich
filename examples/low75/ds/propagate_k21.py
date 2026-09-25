"""Propagate genus labels along similar edges of the k21 MEGAHIT graph.

One hop cannot cross the edges that join two unlabelled contigs. Each round
copies a genus onto an unlabelled node when labelled neighbours that also
share canonical 4-mer cosine of at least 0.5 agree on 80% of their length.
Unclassified reads mapped to a node then inherit that genus.
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

MIN_SIM = 0.5
AGREE = 0.8
ROUNDS = 8


def _parser():
    return ncbi_parser()


def _genus(parser, taxon_id: int) -> int:
    if taxon_id <= 0:
        return 0
    return parser.get_ancestor_by_rank(taxon_id, "genus") or 0


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


def _similar_neighbours(edges, labels, vectors) -> dict[str, list[tuple[str, int]]]:
    neighbours: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for edge in edges:
        left, right = edge["source"], edge["target"]
        left_vector, right_vector = vectors.get(left), vectors.get(right)
        if left_vector is None or right_vector is None:
            continue
        if float(np.dot(left_vector, right_vector)) < MIN_SIM:
            continue
        neighbours[left].append(right)
        neighbours[right].append(left)
    return neighbours


def _propagate(labels, neighbours, lengths) -> tuple[dict[str, int], int]:
    current = dict(labels)
    changed_total = 0
    for _round in range(ROUNDS):
        updates = {}
        for node_id, genus_id in current.items():
            if genus_id:
                continue
            votes: Counter = Counter()
            for other in neighbours.get(node_id, ()):
                other_genus = current.get(other, 0)
                if other_genus:
                    votes[other_genus] += lengths.get(other, 1)
            if not votes:
                continue
            best, mass = votes.most_common(1)[0]
            if mass / sum(votes.values()) >= AGREE:
                updates[node_id] = best
        if not updates:
            break
        current.update(updates)
        changed_total += len(updates)
    return current, changed_total


def _classifier(path: Path) -> dict[str, int]:
    calls = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split("\t")
        calls[parts[1].split("/")[0]] = 0 if parts[0] == "U" else int(parts[2])
    return calls


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


def _truth(example: Path) -> dict[int, float]:
    accessions = {
        row["accession"]: int(row["genus_taxid"])
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


def _node_truth(example: Path, paf: Path) -> dict[str, int]:
    genus_of = {
        row["accession"].split(".")[0]: int(row["genus_taxid"])
        for row in csv.DictReader((example / "accessions.tsv").open(encoding="utf-8"), delimiter="\t")
        if row["role"] == "sim"
    }
    best: dict[str, tuple[int, str]] = {}
    for line in paf.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split("\t")
        if len(parts) < 10:
            continue
        matches = int(parts[9])
        previous = best.get(parts[0])
        if previous is None or matches > previous[0]:
            best[parts[0]] = (matches, parts[5].split()[0])
    return {node_id: genus_of[target] for node_id, (_matches, target) in best.items() if target in genus_of}


def score(example_name: str, parser) -> dict:
    example = ROOT / "examples" / example_name
    work = example / "work"
    edges, sequences = assembly_graph_from_fastg((work / "megahit" / "initial_k21.fastg").read_text(encoding="utf-8", errors="replace"))
    fasta = _fasta(work / "megahit" / "initial_k21_nodes.fa")
    calls = {}
    for line in (work / "reprofile" / "initial_k21.kraken").read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split("\t")
        calls[parts[1]] = 0 if parts[0] == "U" else int(parts[2])
    labels = {node_id: _genus(parser, calls.get(node_id, 0)) for node_id in sequences}
    lengths = {node_id: len(sequence) for node_id, sequence in sequences.items()}
    vectors = {node_id: canonical_4mer(fasta[node_id]) for node_id in sequences if node_id in fasta}
    neighbours = _similar_neighbours(edges, labels, vectors)
    filled, changed = _propagate(labels, neighbours, lengths)
    truth_nodes = _node_truth(example, work / "reprofile" / "initial_k21.paf")
    bases = Counter()
    for node_id, genus_truth in truth_nodes.items():
        predicted = filled.get(node_id, 0)
        length = lengths.get(node_id, 0)
        if predicted and predicted == genus_truth:
            bases["correct"] += length
        elif predicted == 0:
            bases["unclassified"] += length
        else:
            bases["wrong_genus"] += length
    total_bases = sum(bases.values()) or 1
    placement = _placement(work / "reprofile" / "initial_k21_reads.paf")
    truth = _truth(example)

    def profile(path: Path) -> tuple[float, float]:
        counts: Counter = Counter()
        for read_id, taxon_id in _classifier(path).items():
            genus_id = _genus(parser, taxon_id)
            if genus_id == 0:
                contig_id = placement.get(read_id)
                if contig_id is not None:
                    genus_id = filled.get(contig_id, 0)
            counts[genus_id] += 1
        total = sum(counts.values())
        predicted = {taxon_id: count / total for taxon_id, count in counts.items()}
        keys = set(predicted) | set(truth)
        return sum(abs(predicted.get(key, 0.0) - truth.get(key, 0.0)) for key in keys), predicted.get(0, 0.0)

    kraken_l1, kraken_u = profile(work / "classify" / "kraken2.output")
    kaiju_l1, kaiju_u = profile(work / "classify" / "kaiju.output")
    return {
        "dataset": example_name,
        "nodes_relabelled": changed,
        "correct_base_share": f"{bases['correct'] / total_bases:.6g}",
        "wrong_base_share": f"{bases['wrong_genus'] / total_bases:.6g}",
        "unclassified_base_share": f"{bases['unclassified'] / total_bases:.6g}",
        "kraken2_l1": f"{kraken_l1:.6g}",
        "kraken2_unclassified": f"{kraken_u:.6g}",
        "kaiju_l1": f"{kaiju_l1:.6g}",
        "kaiju_unclassified": f"{kaiju_u:.6g}",
    }


def main() -> None:
    parser = _parser()
    rows = [score(name, parser) for name in ("low75", "low75half")]
    for row in rows:
        print(row, flush=True)
    destination = ROOT / "examples" / "low75" / "ds" / "propagate_k21.tsv"
    write_tsv(destination, rows, list(rows[0].keys()))
    (ROOT / "examples" / "low75half" / "ds" / "propagate_k21.tsv").write_text(destination.read_text(encoding="utf-8"), encoding="utf-8")


if __name__ == "__main__":
    main()
