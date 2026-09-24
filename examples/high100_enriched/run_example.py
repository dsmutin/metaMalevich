"""high100 at ten times the read depth: 1000000 paired fragments.

Families and accessions are the high100 pin. Score at family. Indexes are
reused once the high100 Kraken2 and Kaiju builds exist.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "examples" / "high100"
HELD = ROOT / "examples" / "heldout_genera"
sys.path.insert(0, str(HELD))

import run_example as heldout  # noqa: E402
from community import load_pairs  # noqa: E402

heldout.EXAMPLE = HERE
heldout.WORK = HERE / "work"
heldout.GRAPH_ID = "high100_enriched"
heldout.SCORE_RANK = "family"
heldout.TOTAL_READS = 1_000_000
heldout.REPORT_PATH = ROOT / "data" / "raw" / "high100_assembly_data_report.jsonl"


def stage_download() -> None:
    """Require the high100 FASTA pin. This example does not download again."""
    pinned = load_pairs(HERE / "accessions.tsv")
    parent = load_pairs(PARENT / "accessions.tsv")
    if pinned != parent:
        raise SystemExit("examples/high100_enriched/accessions.tsv is not the high100 pin")
    if len({row["pair_id"] for row in pinned}) != 100:
        raise SystemExit("high100_enriched accession table must contain 100 families")
    missing = [
        row["accession"]
        for row in pinned
        if not (ROOT / "data" / "raw" / "fasta" / f"{row['accession']}.fna").is_file()
    ]
    if missing:
        raise SystemExit("missing FASTA for high100_enriched: " + ", ".join(missing))


def stage_databases() -> None:
    """Reuse the high100 indexes when they exist. Do not build a second copy."""
    for name, marker in (("kraken_db", "hash.k2d"), ("kaiju_db", "kaiju_db.fmi")):
        dest = HERE / "work" / name
        source = PARENT / "work" / name
        if dest.exists():
            continue
        ready = source.is_dir() and (source / marker).is_file()
        if not ready:
            raise SystemExit(f"parent index is not ready: {source / marker}")
        dest.symlink_to(source)
    heldout_databases()


heldout_databases = heldout.stage_databases
heldout.stage_download = stage_download
heldout.stage_databases = stage_databases


if __name__ == "__main__":
    heldout.main()
