"""L1, Pearson R2, and presence F1 for every read-level method.

Baselines are Kraken2 and Kaiju alone. The other rows fill an unclassified
read from a MEGAHIT neighbour genus or family, or from the nearest labelled
contig in canonical 4-mer cosine. Outputs land in ``figures/inference/``.
"""

from __future__ import annotations

import csv
import importlib.util
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "examples" / "heldout_genera"))

from community import assembly_graph_from_fastg, presence_f1, r_squared  # noqa: E402
from metamalevich.composition_genus import canonical_4mer  # noqa: E402
from metamalevich.evaluate import abundance_scores  # noqa: E402

OUT = ROOT / "figures" / "inference"
TAXDUMP = Path("/mnt/tank/scratch/partition-metagenomics/databases/taxdump/nodes.dmp")
ENGINE = Path("/mnt/tank/scratch/dsmutin/tools/my/samovar/samovar/src/samovar/taxonomy_engine.py")
EXAMPLES = [
    ("low75", "genus", "genus_taxid"),
    ("low75half", "genus", "genus_taxid"),
    ("low75_enriched", "genus", "genus_taxid"),
    ("low75half_enriched", "genus", "genus_taxid"),
    ("high100", "family", "family_taxid"),
    ("high100_enriched", "family", "family_taxid"),
    ("half100half_enriched", "family", "family_taxid"),
]


def _parser():
    spec = importlib.util.spec_from_file_location("taxonomy_engine", ENGINE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.NCBITaxonomyParser(str(TAXDUMP))


def _at_rank(parser, taxon_id: int, rank: str) -> int:
    if taxon_id <= 0:
        return 0
    return parser.get_ancestor_by_rank(taxon_id, rank) or 0


def _classifier(path: Path) -> dict[str, int]:
    calls = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split("\t")
        calls[parts[1].split("/")[0]] = 0 if parts[0] == "U" else int(parts[2])
    return calls


def _kraken_calls(path: Path) -> dict[str, int]:
    calls = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split("\t")
        calls[parts[1]] = 0 if parts[0] == "U" else int(parts[2])
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


def _reads(example: Path) -> Path:
    candidates = [
        path
        for path in (example / "work" / "iss").rglob("*_full_R1.fastq")
        if path.is_file() and path.stat().st_size > 0 and ".iss_full" not in path.parts
    ]
    if len(candidates) != 1:
        raise SystemExit(f"expected one R1 in {example}, found {candidates}")
    return candidates[0]


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


def _fill(labels: dict[str, int], edges: list[dict], lengths: dict[str, int]) -> dict[str, int]:
    neighbours: dict[str, set[str]] = defaultdict(set)
    for edge in edges:
        neighbours[edge["source"]].add(edge["target"])
        neighbours[edge["target"]].add(edge["source"])
    filled = dict(labels)
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
            filled[node_id] = best
    return filled


def _nearest(sequences: dict[str, str], labels: dict[str, int]) -> dict[str, int]:
    donors = [node_id for node_id, label in labels.items() if label > 0 and node_id in sequences]
    queries = [node_id for node_id, label in labels.items() if label <= 0 and node_id in sequences]
    result = dict(labels)
    if not donors or not queries:
        return result
    donor_matrix = np.vstack([canonical_4mer(sequences[node_id]) for node_id in donors])
    donor_labels = np.array([labels[node_id] for node_id in donors])
    for start in range(0, len(queries), 2000):
        batch = queries[start : start + 2000]
        query_matrix = np.vstack([canonical_4mer(sequences[node_id]) for node_id in batch])
        scores = query_matrix @ donor_matrix.T
        best = scores.argmax(axis=1)
        for index, node_id in enumerate(batch):
            if float(np.linalg.norm(query_matrix[index])) == 0.0:
                continue
            result[node_id] = int(donor_labels[int(best[index])])
    return result


def _fasta_from_graph(sequences: dict[str, str]) -> dict[str, str]:
    return sequences


def _profile(calls: dict[str, int], placement: dict[str, str], node_labels: dict[str, int] | None, parser, rank: str) -> dict[int, float]:
    counts: Counter = Counter()
    for read_id, taxon_id in calls.items():
        label = _at_rank(parser, taxon_id, rank)
        if node_labels is not None and label == 0:
            contig_id = placement.get(read_id)
            if contig_id is not None:
                label = node_labels.get(contig_id, 0)
        counts[label] += 1
    total = sum(counts.values())
    return {taxon_id: count / total for taxon_id, count in counts.items()}


def _metric_row(example: str, rank: str, group: str, method: str, predicted: dict[int, float], truth: dict[int, float]) -> dict[str, str]:
    abundance = abundance_scores(predicted, truth)
    pearson = abundance["pearson"]
    presence = presence_f1(predicted, truth)
    return {
        "example": example,
        "rank": rank,
        "group": group,
        "method": method,
        "l1": f"{abundance['l1']:.6g}",
        "r2": "" if pearson is None else f"{r_squared(pearson):.6g}",
        "f1": f"{presence['f1']:.6g}",
        "unclassified": f"{predicted.get(0, 0.0):.6g}",
    }


def _graph(path: Path):
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f"missing graph: {path}")
    return assembly_graph_from_fastg(path.read_text(encoding="utf-8", errors="replace"))


def score_example(name: str, rank: str, column: str, parser) -> list[dict[str, str]]:
    example = ROOT / "examples" / name
    work = example / "work"
    print(f"score {name}", flush=True)
    final_edges, final_sequences = _graph(work / "megahit" / "assembly.fastg")
    initial_edges, initial_sequences = _graph(work / "megahit" / "initial_k21.fastg")
    final_lengths = {node_id: len(sequence) for node_id, sequence in final_sequences.items()}
    initial_lengths = {node_id: len(sequence) for node_id, sequence in initial_sequences.items()}
    final_calls = _kraken_calls(work / "reprofile" / "contigs.kraken")
    initial_calls = _kraken_calls(work / "reprofile" / "initial_k21.kraken")
    final_labels = {node_id: _at_rank(parser, final_calls.get(node_id, 0), rank) for node_id in final_sequences}
    initial_labels = {node_id: _at_rank(parser, initial_calls.get(node_id, 0), rank) for node_id in initial_sequences}
    final_filled = _fill(final_labels, final_edges, final_lengths)
    initial_filled = _fill(initial_labels, initial_edges, initial_lengths)
    print(f"  4mer {name}", flush=True)
    composed = _nearest(_fasta_from_graph(final_sequences), final_labels)
    final_place = _placement(work / "reprofile" / "reads_to_contigs.paf")
    initial_place = _placement(work / "reprofile" / "initial_k21_reads.paf")
    truth = _truth(example, column)
    rows = []
    for tool, path in (("Kraken2", work / "classify" / "kraken2.output"), ("Kaiju", work / "classify" / "kaiju.output")):
        calls = _classifier(path)
        rows.append(_metric_row(name, rank, "baseline", tool, _profile(calls, {}, None, parser, rank), truth))
        rows.append(
            _metric_row(
                name,
                rank,
                "method",
                f"{tool} + k141",
                _profile(calls, final_place, final_filled, parser, rank),
                truth,
            )
        )
        rows.append(
            _metric_row(
                name,
                rank,
                "method",
                f"{tool} + k21",
                _profile(calls, initial_place, initial_filled, parser, rank),
                truth,
            )
        )
        rows.append(
            _metric_row(
                name,
                rank,
                "method",
                f"{tool} + 4-mer",
                _profile(calls, final_place, composed, parser, rank),
                truth,
            )
        )
        print(f"  {name} {tool} baseline L1 {rows[-4]['l1']} F1 {rows[-4]['f1']} R2 {rows[-4]['r2']}", flush=True)
    return rows


def main() -> None:
    parser = _parser()
    rows: list[dict[str, str]] = []
    for name, rank, column in EXAMPLES:
        rows.extend(score_example(name, rank, column, parser))
    path = OUT / "read_metrics.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(path)


if __name__ == "__main__":
    main()
