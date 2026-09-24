"""Choose a neighbour genus on the k21 graph from node and edge features.

Plurality copies the wrong genus more often than the right one. A logistic
model is fit on one community and applied to the other. Features are 4-mer
cosine, contig length, MEGAHIT coverage, and the donor's Kraken k-mer
support. An unclassified read inherits the chosen genus only when the model
probability is at least 0.5.
"""

from __future__ import annotations

import csv
import importlib.util
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "examples" / "heldout_genera"))
sys.path.insert(0, str(ROOT / "src"))

from community import assembly_graph_from_fastg, megahit_coverage  # noqa: E402
from metamalevich.composition_genus import canonical_4mer  # noqa: E402
from metamalevich.tables import write_tsv  # noqa: E402

TAXDUMP = Path("/mnt/tank/scratch/partition-metagenomics/databases/taxdump/nodes.dmp")
ENGINE = Path("/mnt/tank/scratch/dsmutin/tools/my/samovar/samovar/src/samovar/taxonomy_engine.py")
COV = re.compile(r"_cov_([0-9]+(?:\.[0-9]+)?)")
SEED = 0


def _parser():
    spec = importlib.util.spec_from_file_location("taxonomy_engine", ENGINE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.NCBITaxonomyParser(str(TAXDUMP))


def _genus(parser, taxon_id: int) -> int:
    if taxon_id <= 0:
        return 0
    return parser.get_ancestor_by_rank(taxon_id, "genus") or 0


def _coverage(node_id: str) -> float:
    found = megahit_coverage(node_id)
    if found is not None:
        return found
    match = COV.search(node_id)
    return float(match.group(1)) if match else 1.0


def _kraken(path: Path) -> tuple[dict[str, int], dict[str, float]]:
    calls = {}
    support = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split("\t")
        if len(parts) < 5:
            raise SystemExit(f"kraken row has fewer than 5 columns: {path}")
        taxon_id = 0 if parts[0] == "U" else int(parts[2])
        calls[parts[1]] = taxon_id
        total = 0
        called = 0
        for piece in parts[4].split():
            if ":" not in piece:
                continue
            raw_tax, raw_count = piece.split(":", 1)
            count = int(raw_count)
            total += count
            if taxon_id and int(raw_tax) == taxon_id:
                called += count
        support[parts[1]] = called / total if total else 0.0
    if not calls:
        raise SystemExit(f"kraken output is empty: {path}")
    return calls, support


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
        counts[accessions[accession]] += int(raw)
    total = sum(counts.values())
    return {taxon_id: count / total for taxon_id, count in counts.items()}


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


def load(example_name: str, parser) -> dict:
    example = ROOT / "examples" / example_name
    work = example / "work"
    edges, sequences = assembly_graph_from_fastg((work / "megahit" / "initial_k21.fastg").read_text(encoding="utf-8", errors="replace"))
    fasta = _fasta(work / "megahit" / "initial_k21_nodes.fa")
    calls, support = _kraken(work / "reprofile" / "initial_k21.kraken")
    labels = {node_id: _genus(parser, calls.get(node_id, 0)) for node_id in sequences}
    degree: Counter = Counter()
    for edge in edges:
        degree[edge["source"]] += 1
        degree[edge["target"]] += 1
    vectors = {node_id: canonical_4mer(fasta[node_id]) for node_id in sequences if node_id in fasta}
    return {
        "example": example,
        "edges": edges,
        "labels": labels,
        "support": support,
        "lengths": {node_id: len(sequence) for node_id, sequence in sequences.items()},
        "coverage": {node_id: _coverage(node_id) for node_id in sequences},
        "degree": degree,
        "vectors": vectors,
        "truth_nodes": _node_truth(example, work / "reprofile" / "initial_k21.paf"),
        "placement": _placement(work / "reprofile" / "initial_k21_reads.paf"),
        "kraken_reads": _classifier(work / "classify" / "kraken2.output"),
        "kaiju_reads": _classifier(work / "classify" / "kaiju.output"),
        "read_truth": _read_truth(example),
    }


def _features(graph, query: str, donor: str) -> np.ndarray | None:
    query_vector = graph["vectors"].get(query)
    donor_vector = graph["vectors"].get(donor)
    if query_vector is None or donor_vector is None:
        return None
    query_cov = graph["coverage"].get(query, 1.0)
    donor_cov = graph["coverage"].get(donor, 1.0)
    return np.array(
        [
            float(np.dot(query_vector, donor_vector)),
            np.log1p(graph["lengths"].get(query, 1)),
            np.log1p(graph["lengths"].get(donor, 1)),
            np.log1p(query_cov),
            np.log1p(donor_cov),
            np.log((donor_cov + 0.1) / (query_cov + 0.1)),
            graph["support"].get(donor, 0.0),
            np.log1p(graph["degree"].get(donor, 0)),
        ],
        dtype=float,
    )


def _pairs(graph) -> tuple[np.ndarray, np.ndarray, list[tuple[str, str, int]]]:
    rows = []
    labels = []
    meta = []
    seen = set()
    for edge in graph["edges"]:
        for query, donor in ((edge["source"], edge["target"]), (edge["target"], edge["source"])):
            donor_genus = graph["labels"].get(donor, 0)
            truth = graph["truth_nodes"].get(query)
            if not donor_genus or truth is None:
                continue
            key = (query, donor)
            if key in seen:
                continue
            seen.add(key)
            features = _features(graph, query, donor)
            if features is None:
                continue
            rows.append(features)
            labels.append(1.0 if donor_genus == truth else 0.0)
            meta.append((query, donor, donor_genus))
    if not rows:
        raise SystemExit("no labelled neighbour pairs")
    return np.vstack(rows), np.array(labels), meta


def _fit(features: np.ndarray, labels: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = features.mean(axis=0)
    scale = features.std(axis=0)
    scale[scale == 0] = 1.0
    scaled = (features - mean) / scale
    design = np.column_stack([np.ones(len(scaled)), scaled])
    rng = np.random.default_rng(SEED)
    weight = rng.normal(0, 0.01, size=design.shape[1])
    for _step in range(400):
        probability = 1.0 / (1.0 + np.exp(-np.clip(design @ weight, -30, 30)))
        gradient = design.T @ (probability - labels) / len(labels)
        weight -= 0.5 * gradient
    return weight, mean, scale


def _predict(features: np.ndarray, weight: np.ndarray, mean: np.ndarray, scale: np.ndarray) -> np.ndarray:
    scaled = (features - mean) / scale
    design = np.column_stack([np.ones(len(scaled)), scaled])
    return 1.0 / (1.0 + np.exp(-np.clip(design @ weight, -30, 30)))


def _apply(graph, weight, mean, scale) -> dict[str, int]:
    chosen = dict(graph["labels"])
    grouped: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for edge in graph["edges"]:
        for query, donor in ((edge["source"], edge["target"]), (edge["target"], edge["source"])):
            if graph["labels"].get(query, 0):
                continue
            donor_genus = graph["labels"].get(donor, 0)
            if not donor_genus:
                continue
            features = _features(graph, query, donor)
            if features is None:
                continue
            grouped[query].append((donor, donor_genus, features))
    for query, candidates in grouped.items():
        matrix = np.vstack([item[2] for item in candidates])
        probability = _predict(matrix, weight, mean, scale)
        best = int(np.argmax(probability))
        if probability[best] >= 0.5:
            chosen[query] = candidates[best][1]
    return chosen


def _l1(graph, filled, reads, parser) -> tuple[float, float]:
    counts: Counter = Counter()
    for read_id, taxon_id in reads.items():
        genus_id = _genus(parser, taxon_id)
        if genus_id == 0:
            contig_id = graph["placement"].get(read_id)
            if contig_id is not None:
                genus_id = filled.get(contig_id, 0)
        counts[genus_id] += 1
    total = sum(counts.values())
    predicted = {taxon_id: count / total for taxon_id, count in counts.items()}
    truth = graph["read_truth"]
    keys = set(predicted) | set(truth)
    return sum(abs(predicted.get(key, 0.0) - truth.get(key, 0.0)) for key in keys), predicted.get(0, 0.0)


def _base_share(graph, filled) -> tuple[float, float, float]:
    counts: Counter = Counter()
    for node_id, genus_truth in graph["truth_nodes"].items():
        predicted = filled.get(node_id, 0)
        length = graph["lengths"].get(node_id, 0)
        if predicted and predicted == genus_truth:
            counts["correct"] += length
        elif predicted == 0:
            counts["unclassified"] += length
        else:
            counts["wrong_genus"] += length
    total = sum(counts.values()) or 1
    return counts["correct"] / total, counts["wrong_genus"] / total, counts["unclassified"] / total


def main() -> None:
    parser = _parser()
    graphs = {name: load(name, parser) for name in ("low75", "low75half")}
    rows = []
    for train_name, test_name in (("low75", "low75half"), ("low75half", "low75")):
        features, labels, _meta = _pairs(graphs[train_name])
        weight, mean, scale = _fit(features, labels)
        filled = _apply(graphs[test_name], weight, mean, scale)
        correct, wrong, unclassified = _base_share(graphs[test_name], filled)
        kraken_l1, kraken_u = _l1(graphs[test_name], filled, graphs[test_name]["kraken_reads"], parser)
        kaiju_l1, kaiju_u = _l1(graphs[test_name], filled, graphs[test_name]["kaiju_reads"], parser)
        train_probability = _predict(features, weight, mean, scale)
        train_accuracy = float(((train_probability >= 0.5) == (labels >= 0.5)).mean())
        row = {
            "train": train_name,
            "test": test_name,
            "train_pairs": len(labels),
            "train_positive_rate": f"{labels.mean():.6g}",
            "train_accuracy": f"{train_accuracy:.6g}",
            "correct_base_share": f"{correct:.6g}",
            "wrong_base_share": f"{wrong:.6g}",
            "unclassified_base_share": f"{unclassified:.6g}",
            "kraken2_l1": f"{kraken_l1:.6g}",
            "kraken2_unclassified": f"{kraken_u:.6g}",
            "kaiju_l1": f"{kaiju_l1:.6g}",
            "kaiju_unclassified": f"{kaiju_u:.6g}",
        }
        rows.append(row)
        print(row, flush=True)
    destination = ROOT / "examples" / "low75" / "ds" / "colour_features.tsv"
    write_tsv(destination, rows, list(rows[0].keys()))
    (ROOT / "examples" / "low75half" / "ds" / "colour_features.tsv").write_text(
        destination.read_text(encoding="utf-8"), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
