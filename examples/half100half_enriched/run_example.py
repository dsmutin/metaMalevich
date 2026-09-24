"""half100half at ten times the read depth: 1000000 paired fragments.

The accession table is every other high100 family. Score at family.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "examples" / "half100half"
HELD = ROOT / "examples" / "heldout_genera"
sys.path.insert(0, str(HELD))

import run_example as heldout  # noqa: E402
from community import load_pairs  # noqa: E402

heldout.EXAMPLE = HERE
heldout.WORK = HERE / "work"
heldout.GRAPH_ID = "half100half_enriched"
heldout.SCORE_RANK = "family"
heldout.TOTAL_READS = 1_000_000
heldout.REPORT_PATH = ROOT / "data" / "raw" / "high100_assembly_data_report.jsonl"


def stage_download() -> None:
    """Require the half100half FASTA pin. This example does not download again."""
    pinned = load_pairs(HERE / "accessions.tsv")
    parent = load_pairs(PARENT / "accessions.tsv")
    if pinned != parent:
        raise SystemExit("examples/half100half_enriched/accessions.tsv is not the half100half pin")
    missing = [
        row["accession"]
        for row in pinned
        if not (ROOT / "data" / "raw" / "fasta" / f"{row['accession']}.fna").is_file()
    ]
    if missing:
        raise SystemExit("missing FASTA for half100half_enriched: " + ", ".join(missing))


def stage_databases() -> None:
    """Reuse the half100half indexes when they exist."""
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
