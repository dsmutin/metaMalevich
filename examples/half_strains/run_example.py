"""Run the held-out pipeline on every other downloaded strain.

The accession table is the even positions of ``examples/heldout_genera/accessions.tsv``
(ten strains, both assemblies of each). FASTA files come from the download
already stored under ``data/raw/fasta``. This script does not download again
and does not rewrite that download's manifest.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
HELD = ROOT / "examples" / "heldout_genera"
sys.path.insert(0, str(HELD))

import run_example as heldout  # noqa: E402
from community import half_strains, load_pairs  # noqa: E402

heldout.EXAMPLE = HERE
heldout.WORK = HERE / "work"
heldout.GRAPH_ID = "half_strains"


def stage_download() -> None:
    """Require the selected FASTA files from the full held-out download."""
    pinned = load_pairs(HERE / "accessions.tsv")
    expected = half_strains(load_pairs(HELD / "accessions.tsv"))
    if pinned != expected:
        raise SystemExit("examples/half_strains/accessions.tsv is not every other held-out strain")
    missing = []
    for row in pinned:
        path = ROOT / "data" / "raw" / "fasta" / f"{row['accession']}.fna"
        if not path.is_file() or path.stat().st_size == 0:
            missing.append(row["accession"])
    if missing:
        raise SystemExit(
            "missing FASTA for the half-strain example; run examples/heldout_genera/run_example.py --stage download first: "
            + ", ".join(missing)
        )


heldout.stage_download = stage_download


if __name__ == "__main__":
    heldout.main()
