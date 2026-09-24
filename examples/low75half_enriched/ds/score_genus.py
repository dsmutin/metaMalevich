"""Genus scores for the 10x-depth low75half assembly.

The final graph is ``contig2fastg`` on ``k141.contigs.fa``. The initial graph
is the same tool on ``k21.contigs.fa``. Read baselines are Kraken2 and Kaiju.
Unclassified reads inherit a contig genus from neighbour agreement on the
final graph, from the same rule on the initial graph, or from the nearest
Kraken-labelled contig in canonical 4-mer cosine.
"""

from __future__ import annotations

import csv
import importlib.util
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "examples" / "heldout_genera"))
sys.path.insert(0, str(ROOT / "src"))

from community import assembly_graph_from_fastg  # noqa: E402
from metamalevich.composition_genus import canonical_4mer  # noqa: E402
from metamalevich.tables import write_tsv  # noqa: E402

EXAMPLE = ROOT / "examples" / "low75half_enriched"
WORK = EXAMPLE / "work"
TAXDUMP = Path("/mnt/tank/scratch/partition-metagenomics/databases/taxdump/nodes.dmp")
ENGINE = Path("/mnt/tank/scratch/dsmutin/tools/my/samovar/samovar/src/samovar/taxonomy_engine.py")
KRAKEN = Path("/nfs/home/dsmutin/miniconda3/bin/kraken2")


def _parser():
    spec = importlib.util.spec_from_file_location("taxonomy_engine", ENGINE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.NCBITaxonomyParser(str(TAXDUMP))


def _genus(parser, taxon_id: int) -> int:
    if taxon_id <= 0:
        return 0
    return parser.get_ancestor_by_rank(taxon_id, "genus") or 0


def _write_fasta(sequences: dict[str, str], path: Path) -> None:
    if path.is_file() and path.stat().st_size > 0:
        return
    with path.open("w", encoding="utf-8") as handle:
        for node_id, sequence in sequences.items():
            handle.write(f">{node_id}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def _kraken(fasta: Path, output: Path) -> dict[str, int]:
    if not output.is_file() or output.stat().st_size == 0:
        completed = subprocess.run(
            [str(KRAKEN), "--db", str(WORK / "kraken_db"), "--threads", "8", "--output", str(output), str(fasta)],
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode != 0:
            raise SystemExit(completed.stderr.strip() or "kraken2 failed")
    calls = {}
    for line in output.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split("\t")
        calls[parts[1]] = 0 if parts[0] == "U" else int(parts[2])
    if not calls:
        raise SystemExit(f"empty kraken output: {output}")
    return calls


def _fill(labels: dict[str, int], edges: list[dict], lengths: dict[str, int]) -> dict[str, int]:
    neighbours: dict[str, set[str]] = defaultdict(set)
    for edge in edges:
        neighbours[edge["source"]].add(edge["target"])
        neighbours[edge["target"]].add(edge["source"])
    filled = dict(labels)
    for node_id, genus_id in labels.items():
        if genus_id:
            continue
        votes: Counter = Counter()
        for other in neighbours.get(node_id, ()):
            other_genus = labels.get(other, 0)
            if other_genus:
                votes[other_genus] += lengths.get(other, 1)
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


def _map_reads(fasta: Path, fastq: Path, destination: Path) -> dict[str, str]:
    if not destination.is_file() or destination.stat().st_size == 0:
        minimap = shutil.which("minimap2")
        completed = subprocess.run(
            [minimap, "-x", "sr", "-t", "8", str(fasta), str(fastq)],
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode != 0:
            raise SystemExit(completed.stderr.strip() or "minimap2 failed")
        destination.write_text(completed.stdout, encoding="utf-8")
    best: dict[str, tuple[int, str]] = {}
    for line in destination.read_text(encoding="utf-8", errors="replace").splitlines():
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


def _truth() -> dict[int, float]:
    accessions = {
        row["accession"]: int(row["genus_taxid"])
        for row in csv.DictReader((EXAMPLE / "accessions.tsv").open(encoding="utf-8"), delimiter="\t")
    }
    counts: Counter = Counter()
    for line in (WORK / "abundance.csv").read_text(encoding="utf-8").splitlines()[1:]:
        if not line.strip():
            continue
        accession, raw = line.split(",")
        counts[accessions[accession]] += int(raw)
    total = sum(counts.values())
    return {taxon_id: count / total for taxon_id, count in counts.items()}


def _l1(predicted: dict[int, float], truth: dict[int, float]) -> float:
    keys = set(predicted) | set(truth)
    return float(sum(abs(predicted.get(key, 0.0) - truth.get(key, 0.0)) for key in keys))


def _profile(reads: dict[str, int], placement: dict[str, str], labels: dict[str, int] | None, parser) -> dict[int, float]:
    counts: Counter = Counter()
    for read_id, taxon_id in reads.items():
        genus_id = _genus(parser, taxon_id)
        if labels is not None and genus_id == 0:
            contig_id = placement.get(read_id)
            if contig_id is not None:
                genus_id = labels.get(contig_id, 0)
        counts[genus_id] += 1
    total = sum(counts.values())
    return {taxon_id: count / total for taxon_id, count in counts.items()}


def _graph(path: Path) -> tuple[list[dict], dict[str, str]]:
    return assembly_graph_from_fastg(path.read_text(encoding="utf-8", errors="replace"))


def main() -> None:
    parser = _parser()
    final_edges, final_sequences = _graph(WORK / "megahit" / "assembly.fastg")
    initial_edges, initial_sequences = _graph(WORK / "megahit" / "initial_k21.fastg")
    final_fasta = WORK / "megahit" / "graph_nodes.fa"
    initial_fasta = WORK / "megahit" / "initial_k21_nodes.fa"
    _write_fasta(final_sequences, final_fasta)
    _write_fasta(initial_sequences, initial_fasta)
    reprofile = WORK / "reprofile"
    reprofile.mkdir(parents=True, exist_ok=True)
    final_calls = _kraken(final_fasta, reprofile / "contigs.kraken")
    initial_calls = _kraken(initial_fasta, reprofile / "initial_k21.kraken")
    final_labels = {node_id: _genus(parser, final_calls.get(node_id, 0)) for node_id in final_sequences}
    initial_labels = {node_id: _genus(parser, initial_calls.get(node_id, 0)) for node_id in initial_sequences}
    final_lengths = {node_id: len(sequence) for node_id, sequence in final_sequences.items()}
    initial_lengths = {node_id: len(sequence) for node_id, sequence in initial_sequences.items()}
    final_filled = _fill(final_labels, final_edges, final_lengths)
    initial_filled = _fill(initial_labels, initial_edges, initial_lengths)
    print("4mer start", flush=True)
    composed = _nearest(final_sequences, final_labels)
    print("4mer done", flush=True)
    candidates = [
        path
        for path in (WORK / "iss").rglob("*_full_R1.fastq")
        if path.is_file() and path.stat().st_size > 0 and ".iss_full" not in path.parts
    ]
    if len(candidates) != 1:
        raise SystemExit(f"expected one enriched R1, found {[str(path) for path in candidates]}")
    fastq = candidates[0]
    placement = _map_reads(final_fasta, fastq, reprofile / "reads_to_contigs.paf")
    initial_placement = _map_reads(initial_fasta, fastq, reprofile / "initial_k21_reads.paf")
    truth = _truth()
    kraken = _classifier(WORK / "classify" / "kraken2.output")
    kaiju = _classifier(WORK / "classify" / "kaiju.output")
    rows = []
    methods = [
        ("kraken2", _profile(kraken, {}, None, parser)),
        ("kaiju", _profile(kaiju, {}, None, parser)),
        ("kraken2_k141_neighbour", _profile(kraken, placement, final_filled, parser)),
        ("kaiju_k141_neighbour", _profile(kaiju, placement, final_filled, parser)),
        ("kraken2_k21_neighbour", _profile(kraken, initial_placement, initial_filled, parser)),
        ("kaiju_k21_neighbour", _profile(kaiju, initial_placement, initial_filled, parser)),
        ("kraken2_4mer_graph", _profile(kraken, placement, composed, parser)),
        ("kaiju_4mer_graph", _profile(kaiju, placement, composed, parser)),
    ]
    for method, predicted in methods:
        row = {
            "dataset": "low75half_enriched",
            "method": method,
            "rank": "genus",
            "l1": f"{_l1(predicted, truth):.6g}",
            "unclassified": f"{predicted.get(0, 0.0):.6g}",
            "k141_nodes": len(final_sequences),
            "k141_edges": len(final_edges),
            "k21_nodes": len(initial_sequences),
            "k21_edges": len(initial_edges),
        }
        rows.append(row)
        print(f"{method}\t{row['l1']}\t{row['unclassified']}", flush=True)
    destination = EXAMPLE / "ds"
    destination.mkdir(parents=True, exist_ok=True)
    write_tsv(destination / "genus_metrics.tsv", rows, list(rows[0].keys()))


if __name__ == "__main__":
    main()
