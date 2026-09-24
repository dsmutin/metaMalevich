"""Keep a k21 overlap only when the contigs also share a 4-mer genus.

Raw neighbour votes on the initial MEGAHIT graph copy the wrong genus more
often than the right one. This script accepts a labelled neighbour only when
canonical 4-mer cosine is at least 0.5, then fills unclassified reads.
"""

from __future__ import annotations

import csv
import importlib.util
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

TAXDUMP = Path("/mnt/tank/scratch/partition-metagenomics/databases/taxdump/nodes.dmp")
ENGINE = Path("/mnt/tank/scratch/dsmutin/tools/my/samovar/samovar/src/samovar/taxonomy_engine.py")
MIN_SIM = 0.5


def _parser():
    spec = importlib.util.spec_from_file_location("taxonomy_engine", ENGINE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.NCBITaxonomyParser(str(TAXDUMP))


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


def _fill(labels, edges, vectors, lengths) -> tuple[dict[str, int], int, int]:
    neighbours: dict[str, set[str]] = defaultdict(set)
    for edge in edges:
        neighbours[edge["source"]].add(edge["target"])
        neighbours[edge["target"]].add(edge["source"])
    filled = dict(labels)
    kept = dropped = 0
    for node_id, genus_id in labels.items():
        if genus_id:
            continue
        query = vectors.get(node_id)
        if query is None:
            continue
        votes: Counter = Counter()
        for other in neighbours.get(node_id, ()):
            donor = labels.get(other, 0)
            donor_vector = vectors.get(other)
            if not donor or donor_vector is None:
                continue
            similarity = float(np.dot(query, donor_vector))
            if similarity < MIN_SIM:
                dropped += 1
                continue
            kept += 1
            votes[donor] += lengths.get(other, 1)
        if not votes:
            continue
        best, mass = votes.most_common(1)[0]
        if mass / sum(votes.values()) >= 0.8:
            filled[node_id] = best
    return filled, kept, dropped


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


def score(example_name: str, parser) -> dict:
    example = ROOT / "examples" / example_name
    work = example / "work"
    fastg = work / "megahit" / "initial_k21.fastg"
    edges, sequences = assembly_graph_from_fastg(fastg.read_text(encoding="utf-8", errors="replace"))
    fasta = _fasta(work / "megahit" / "initial_k21_nodes.fa")
    calls = {}
    for line in (work / "reprofile" / "initial_k21.kraken").read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split("\t")
        calls[parts[1]] = 0 if parts[0] == "U" else int(parts[2])
    labels = {node_id: _genus(parser, calls.get(node_id, 0)) for node_id in sequences}
    lengths = {node_id: len(sequence) for node_id, sequence in sequences.items()}
    vectors = {node_id: canonical_4mer(sequence) for node_id, sequence in fasta.items()}
    filled, kept, dropped = _fill(labels, edges, vectors, lengths)
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
        l1 = sum(abs(predicted.get(key, 0.0) - truth.get(key, 0.0)) for key in keys)
        return l1, predicted.get(0, 0.0)

    kraken_l1, kraken_u = profile(work / "classify" / "kraken2.output")
    kaiju_l1, kaiju_u = profile(work / "classify" / "kaiju.output")
    return {
        "dataset": example_name,
        "min_sim": MIN_SIM,
        "edges_kept": kept,
        "edges_dropped": dropped,
        "kraken2_k21_4mer_l1": f"{kraken_l1:.6g}",
        "kraken2_unclassified": f"{kraken_u:.6g}",
        "kaiju_k21_4mer_l1": f"{kaiju_l1:.6g}",
        "kaiju_unclassified": f"{kaiju_u:.6g}",
    }


def main() -> None:
    parser = _parser()
    rows = [score(name, parser) for name in ("low75", "low75half")]
    for row in rows:
        print(row, flush=True)
    destination = ROOT / "examples" / "low75" / "ds" / "edge_4mer.tsv"
    write_tsv(destination, rows, list(rows[0].keys()))
    (ROOT / "examples" / "low75half" / "ds" / "edge_4mer.tsv").write_text(destination.read_text(encoding="utf-8"), encoding="utf-8")


if __name__ == "__main__":
    main()
