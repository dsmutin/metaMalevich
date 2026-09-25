"""Score the first MEGAHIT graph against the final one on low75.

``initial_k21.fastg`` is ``megahit_toolkit contig2fastg 21`` on
``k21.contigs.fa``. ``assembly.fastg`` is the same tool at k=141. Neither
file is a mock. An unclassified node copies the genus shared by at least 80%
of the length of its labelled neighbours. Unclassified reads that map to a
node inherit that genus.
"""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "examples" / "heldout_genera"))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "examples"))
from ncbi_taxonomy import kraken2_bin, ncbi_parser  # noqa: E402

from community import assembly_graph_from_fastg  # noqa: E402
from metamalevich.tables import write_tsv  # noqa: E402

AGREE = 0.8


def _parser():
    return ncbi_parser()


def _genus(parser, taxon_id: int) -> int:
    if taxon_id <= 0:
        return 0
    return parser.get_ancestor_by_rank(taxon_id, "genus") or 0


def _write_fasta(sequences: dict[str, str], path: Path) -> None:
    if path.is_file() and path.stat().st_size > 0:
        return
    chunks = []
    for node_id, sequence in sequences.items():
        chunks.append(f">{node_id}\n")
        for start in range(0, len(sequence), 80):
            chunks.append(sequence[start : start + 80] + "\n")
    path.write_text("".join(chunks), encoding="utf-8")


def _kraken(fasta: Path, database: Path, output: Path) -> dict[str, int]:
    if not output.is_file() or output.stat().st_size == 0:
        report = output.with_suffix(".report")
        completed = subprocess.run(
            [str(kraken2_bin()), "--db", str(database), "--threads", "8", "--report", str(report), "--output", str(output), str(fasta)],
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode != 0:
            raise SystemExit(completed.stderr.strip() or "kraken2 failed")
    calls = {}
    for line in output.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            raise SystemExit(f"kraken row has fewer than 3 columns: {output}")
        calls[parts[1]] = 0 if parts[0] == "U" else int(parts[2])
    if not calls:
        raise SystemExit(f"kraken output is empty: {output}")
    return calls


def _truth_nodes(fasta: Path, reference: Path, paf_path: Path, genus_of_target: dict[str, int]) -> dict[str, int]:
    minimap = shutil.which("minimap2")
    if not minimap:
        raise SystemExit("required program is not on PATH: minimap2")
    if not paf_path.is_file() or paf_path.stat().st_size == 0:
        completed = subprocess.run(
            [minimap, "-x", "asm20", "-t", "8", str(reference), str(fasta)],
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode != 0:
            raise SystemExit(completed.stderr.strip() or "minimap2 failed")
        paf_path.write_text(completed.stdout, encoding="utf-8")
    best: dict[str, tuple[int, str]] = {}
    for line in paf_path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split("\t")
        if len(parts) < 10:
            continue
        matches = int(parts[9])
        previous = best.get(parts[0])
        if previous is None or matches > previous[0]:
            best[parts[0]] = (matches, parts[5])
    truth = {}
    for node_id, (_matches, target) in best.items():
        genus_id = genus_of_target.get(target.split()[0])
        if genus_id:
            truth[node_id] = genus_id
    return truth


def _genus_of_sim(example: Path) -> dict[str, int]:
    table = {}
    for row in csv.DictReader((example / "accessions.tsv").open(encoding="utf-8"), delimiter="\t"):
        if row["role"] != "sim":
            continue
        table[row["accession"].split(".")[0]] = int(row["genus_taxid"])
    if not table:
        raise SystemExit(f"no simulated accessions in {example}")
    return table


def _read_truth(example: Path) -> dict[int, float]:
    accessions = {
        row["accession"]: int(row["genus_taxid"])
        for row in csv.DictReader((example / "accessions.tsv").open(encoding="utf-8"), delimiter="\t")
    }
    counts: Counter = Counter()
    for line in (example / "work" / "abundance.csv").read_text(encoding="utf-8").splitlines()[1:]:
        if not line.strip():
            continue
        accession, raw = line.split(",")
        genus_id = accessions.get(accession)
        if genus_id is None:
            raise SystemExit(f"abundance accession is not in the pin table: {accession}")
        counts[genus_id] += int(raw)
    total = sum(counts.values())
    return {taxon_id: count / total for taxon_id, count in counts.items()}


def _fill(labels: dict[str, int], edges: list[dict], lengths: dict[str, int]) -> dict[str, int]:
    """Copy a genus onto an unlabelled node when neighbours agree on 80% of length."""
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
        if mass / sum(votes.values()) >= AGREE:
            filled[node_id] = best
    return filled


def _edge_rows(dataset: str, graph: str, labels: dict[str, int], truth: dict[str, int], edges: list[dict]) -> list[dict]:
    counts: Counter = Counter()
    rescue_right = rescue_wrong = 0
    for edge in edges:
        source, target = edge["source"], edge["target"]
        left, right = labels.get(source, 0), labels.get(target, 0)
        if left and right:
            kind = "both_labelled_same" if left == right else "both_labelled_different"
        elif left or right:
            kind = "one_labelled"
            donor = left or right
            blank = target if left else source
            if blank in truth:
                if donor == truth[blank]:
                    rescue_right += 1
                else:
                    rescue_wrong += 1
        else:
            kind = "both_unlabelled"
        counts[kind] += 1
    rows = []
    total = sum(counts.values()) or 1
    for kind in ("both_labelled_same", "both_labelled_different", "one_labelled", "both_unlabelled"):
        rows.append(
            {
                "dataset": dataset,
                "graph": graph,
                "kind": kind,
                "edges": counts[kind],
                "edge_share": f"{counts[kind] / total:.6g}",
                "rescue_matches_truth": rescue_right if kind == "one_labelled" else "",
                "rescue_wrong_genus": rescue_wrong if kind == "one_labelled" else "",
            }
        )
    return rows


def _base_share(dataset: str, graph: str, method: str, labels: dict[str, int], truth: dict[str, int], lengths: dict[str, int]) -> list[dict]:
    counts: Counter = Counter()
    for node_id, genus_truth in truth.items():
        predicted = labels.get(node_id, 0)
        length = lengths.get(node_id, 0)
        if predicted and predicted == genus_truth:
            kind = "correct"
        elif predicted == 0:
            kind = "unclassified"
        else:
            kind = "wrong_genus"
        counts[kind] += length
    total = sum(counts.values()) or 1
    return [
        {
            "dataset": dataset,
            "graph": graph,
            "method": method,
            "kind": kind,
            "bases": counts[kind],
            "base_share": f"{counts[kind] / total:.6g}",
        }
        for kind in ("correct", "unclassified", "wrong_genus")
    ]


def _map_reads(fasta: Path, fastq: Path, destination: Path) -> dict[str, str]:
    minimap = shutil.which("minimap2")
    if not destination.is_file() or destination.stat().st_size == 0:
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
        if len(parts) < 3:
            raise SystemExit(f"classifier row has fewer than 3 columns: {path}")
        calls[parts[1].split("/")[0]] = 0 if parts[0] == "U" else int(parts[2])
    if not calls:
        raise SystemExit(f"classifier output is empty: {path}")
    return calls


def _l1(predicted: dict[int, float], truth: dict[int, float]) -> float:
    keys = set(predicted) | set(truth)
    return float(sum(abs(predicted.get(key, 0.0) - truth.get(key, 0.0)) for key in keys))


def _fraction(counts: Counter) -> dict[int, float]:
    total = sum(counts.values())
    if total <= 0:
        raise SystemExit("profile has no mass")
    return {taxon_id: count / total for taxon_id, count in counts.items()}


def score_example(example_name: str, parser) -> tuple[list[dict], list[dict], list[dict]]:
    example = ROOT / "examples" / example_name
    work = example / "work"
    megahit = work / "megahit"
    fastg = megahit / "initial_k21.fastg"
    reference = work / "sim_reference.fna"
    fastq = work / "iss" / "initial" / "sample_full_R1.fastq"
    for path in (fastg, reference, fastq, work / "kraken_db" / "hash.k2d"):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"missing input: {path}")
    edges, sequences = assembly_graph_from_fastg(fastg.read_text(encoding="utf-8", errors="replace"))
    lengths = {node_id: len(sequence) for node_id, sequence in sequences.items()}
    fasta = megahit / "initial_k21_nodes.fa"
    _write_fasta(sequences, fasta)
    calls = _kraken(fasta, work / "kraken_db", work / "reprofile" / "initial_k21.kraken")
    labels = {node_id: _genus(parser, calls.get(node_id, 0)) for node_id in sequences}
    truth = _truth_nodes(fasta, reference, work / "reprofile" / "initial_k21.paf", _genus_of_sim(example))
    filled = _fill(labels, edges, lengths)
    edge_rows = _edge_rows(example_name, "k21", labels, truth, edges)
    final_fastg = megahit / "assembly.fastg"
    final_edges, _final_sequences = assembly_graph_from_fastg(final_fastg.read_text(encoding="utf-8", errors="replace"))
    edge_rows.append(
        {
            "dataset": example_name,
            "graph": "k141",
            "kind": "all",
            "edges": len(final_edges),
            "edge_share": "1",
            "rescue_matches_truth": "",
            "rescue_wrong_genus": "",
        }
    )
    error_rows = _base_share(example_name, "k21", "kraken_contig", labels, truth, lengths)
    error_rows.extend(_base_share(example_name, "k21", "neighbour_agree", filled, truth, lengths))
    placement = _map_reads(fasta, fastq, work / "reprofile" / "initial_k21_reads.paf")
    read_truth = _read_truth(example)
    kraken = _classifier(work / "classify" / "kraken2.output")
    kaiju = _classifier(work / "classify" / "kaiju.output")

    def profile(calls_by_read: dict[str, int], node_labels: dict[str, int] | None) -> dict[int, float]:
        counts: Counter = Counter()
        for read_id, taxon_id in calls_by_read.items():
            genus_id = _genus(parser, taxon_id)
            if node_labels is not None and genus_id == 0:
                contig_id = placement.get(read_id)
                if contig_id is not None:
                    genus_id = node_labels.get(contig_id, 0)
            counts[genus_id] += 1
        return _fraction(counts)

    metrics = []
    for method, predicted in (
        ("kraken2", profile(kraken, None)),
        ("kaiju", profile(kaiju, None)),
        ("kraken2_k21_neighbour", profile(kraken, filled)),
        ("kaiju_k21_neighbour", profile(kaiju, filled)),
    ):
        metrics.append(
            {
                "dataset": example_name,
                "method": method,
                "rank": "genus",
                "weight": "read_count",
                "l1": f"{_l1(predicted, read_truth):.6g}",
                "unclassified": f"{predicted.get(0, 0.0):.6g}",
                "k21_nodes": len(sequences),
                "k21_edges": len(edges),
                "k141_edges": len(final_edges),
            }
        )
    return metrics, edge_rows, error_rows


def main() -> None:
    parser = _parser()
    metrics: list[dict] = []
    edges: list[dict] = []
    errors: list[dict] = []
    for name in ("low75", "low75half"):
        metric_rows, edge_rows, error_rows = score_example(name, parser)
        metrics.extend(metric_rows)
        edges.extend(edge_rows)
        errors.extend(error_rows)
        for row in metric_rows:
            print(f"{row['dataset']}\t{row['method']}\t{row['l1']}\t{row['unclassified']}", flush=True)
    destination = ROOT / "examples" / "low75" / "ds"
    write_tsv(destination / "initial_graph_metrics.tsv", metrics, list(metrics[0].keys()))
    write_tsv(destination / "initial_graph_edges.tsv", edges, list(edges[0].keys()))
    write_tsv(destination / "initial_graph_errors.tsv", errors, list(errors[0].keys()))
    half = ROOT / "examples" / "low75half" / "ds"
    for name in ("initial_graph_metrics.tsv", "initial_graph_edges.tsv", "initial_graph_errors.tsv"):
        (half / name).write_text((destination / name).read_text(encoding="utf-8"), encoding="utf-8")


if __name__ == "__main__":
    main()
