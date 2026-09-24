"""Where k141 overlap edges agree with the simulated family on high100_enriched.

Kraken family labels come from ``work/reprofile/contigs.kraken``. Truth is the
best minimap2 alignment of each contig to the simulated genomes.
"""

from __future__ import annotations

import csv
import importlib.util
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "examples" / "heldout_genera"))
sys.path.insert(0, str(ROOT / "src"))

from community import assembly_graph_from_fastg  # noqa: E402
from metamalevich.tables import write_tsv  # noqa: E402

EXAMPLE = ROOT / "examples" / "high100_enriched"
WORK = EXAMPLE / "work"
TAXDUMP = Path("/mnt/tank/scratch/partition-metagenomics/databases/taxdump/nodes.dmp")
ENGINE = Path("/mnt/tank/scratch/dsmutin/tools/my/samovar/samovar/src/samovar/taxonomy_engine.py")


def _parser():
    spec = importlib.util.spec_from_file_location("taxonomy_engine", ENGINE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.NCBITaxonomyParser(str(TAXDUMP))


def _genus(parser, taxon_id: int) -> int:
    if taxon_id <= 0:
        return 0
    return parser.get_ancestor_by_rank(taxon_id, "family") or 0


def _reference() -> Path:
    destination = WORK / "sim_reference.fna"
    if destination.is_file() and destination.stat().st_size > 0:
        return destination
    with destination.open("w", encoding="utf-8") as handle:
        for path in sorted((WORK / "sim").glob("*.fna")):
            sequence = []
            for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
                if not line.startswith(">"):
                    sequence.append(line.strip())
            handle.write(f">{path.stem}\n{''.join(sequence)}\n")
    if destination.stat().st_size == 0:
        raise SystemExit("simulated reference is empty")
    return destination


def _truth(fasta: Path, paf_path: Path) -> dict[str, int]:
    genus_of = {
        row["accession"].split(".")[0]: int(row["family_taxid"])
        for row in csv.DictReader((EXAMPLE / "accessions.tsv").open(encoding="utf-8"), delimiter="\t")
        if row["role"] == "sim"
    }
    if not paf_path.is_file() or paf_path.stat().st_size == 0:
        minimap = shutil.which("minimap2")
        with paf_path.open("w", encoding="utf-8") as handle:
            completed = subprocess.run(
                [minimap, "-x", "asm20", "-t", "4", str(_reference()), str(fasta)],
                check=False,
                stdout=handle,
                stderr=subprocess.PIPE,
                text=True,
            )
        if completed.returncode != 0:
            raise SystemExit(completed.stderr.strip() or "minimap2 failed")
    best: dict[str, tuple[int, str]] = {}
    for line in paf_path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split("\t")
        if len(parts) < 10:
            continue
        matches = int(parts[9])
        previous = best.get(parts[0])
        if previous is None or matches > previous[0]:
            best[parts[0]] = (matches, parts[5].split()[0])
    return {node_id: genus_of[target] for node_id, (_matches, target) in best.items() if target in genus_of}


def main() -> None:
    parser = _parser()
    edges, sequences = assembly_graph_from_fastg((WORK / "megahit" / "assembly.fastg").read_text(encoding="utf-8", errors="replace"))
    calls = {}
    for line in (WORK / "reprofile" / "contigs.kraken").read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split("\t")
        calls[parts[1]] = 0 if parts[0] == "U" else int(parts[2])
    labels = {node_id: _genus(parser, calls.get(node_id, 0)) for node_id in sequences}
    truth = _truth(WORK / "megahit" / "graph_nodes.fa", WORK / "reprofile" / "contigs_vs_sim.paf")
    counts: Counter = Counter()
    rescue_right = rescue_wrong = 0
    for edge in edges:
        left = labels.get(edge["source"], 0)
        right = labels.get(edge["target"], 0)
        if left and right:
            counts["both_labelled_same" if left == right else "both_labelled_different"] += 1
        elif left or right:
            counts["one_labelled"] += 1
            donor = left or right
            blank = edge["target"] if left else edge["source"]
            if blank in truth:
                if donor == truth[blank]:
                    rescue_right += 1
                else:
                    rescue_wrong += 1
        else:
            counts["both_unlabelled"] += 1
        source_truth = truth.get(edge["source"])
        target_truth = truth.get(edge["target"])
        if source_truth and target_truth:
            counts["truth_same" if source_truth == target_truth else "truth_different"] += 1
    total = len(edges) or 1
    rows = []
    for kind in ("both_labelled_same", "both_labelled_different", "one_labelled", "both_unlabelled", "truth_same", "truth_different"):
        rows.append(
            {
                "dataset": "high100_enriched",
                "graph": "k141",
                "kind": kind,
                "edges": counts[kind],
                "edge_share": f"{counts[kind] / total:.6g}",
                "rescue_matches_truth": rescue_right if kind == "one_labelled" else "",
                "rescue_wrong_genus": rescue_wrong if kind == "one_labelled" else "",
            }
        )
    destination = EXAMPLE / "ds" / "k141_edges.tsv"
    write_tsv(destination, rows, list(rows[0].keys()))
    for row in rows:
        print(row, flush=True)


if __name__ == "__main__":
    main()
