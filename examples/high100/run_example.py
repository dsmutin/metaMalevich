"""100 families scored at family rank, one genome simulated and one held out.

Bacteria, archaea, viruses, and small eukaryotes contribute 25 families each.
Thirteen families in each domain use two species of one genus. Twelve use two
genera. FASTA files stay in ``data/raw`` and are not committed. The download
zip is ``data/raw/ncbi_high100.zip``.
"""

from __future__ import annotations

import json
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
HELD = ROOT / "examples" / "heldout_genera"
sys.path.insert(0, str(HELD))

import run_example as heldout  # noqa: E402
from community import load_pairs, sequence_digest  # noqa: E402

heldout.EXAMPLE = HERE
heldout.WORK = HERE / "work"
heldout.GRAPH_ID = "high100"
heldout.SCORE_RANK = "family"
heldout.REPORT_PATH = ROOT / "data" / "raw" / "high100_assembly_data_report.jsonl"
ZIP_PATH = ROOT / "data" / "raw" / "ncbi_high100.zip"


def stage_download() -> None:
    """Download pinned assemblies if needed, then validate this table only."""
    rows = load_pairs(HERE / "accessions.tsv")
    if len({row["pair_id"] for row in rows}) != 100:
        raise SystemExit("high100 accession table must contain 100 families")
    raw = ROOT / "data" / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    log_path = ROOT / "data" / "logs" / "download.log"
    if not ZIP_PATH.is_file() or not zipfile.is_zipfile(ZIP_PATH):
        datasets = heldout._tool("datasets")
        listing = raw / "high100_accessions.txt"
        listing.write_text("".join(f"{row['accession']}\n" for row in rows), encoding="utf-8")
        heldout._run(
            [
                datasets,
                "download",
                "genome",
                "accession",
                "--inputfile",
                str(listing),
                "--include",
                "genome",
                "--filename",
                str(ZIP_PATH),
            ],
            log_path,
        )
    fasta_dir = raw / "fasta"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(ZIP_PATH) as archive:
        report_name = next((name for name in archive.namelist() if name.endswith("assembly_data_report.jsonl")), "")
        if not report_name:
            raise SystemExit("NCBI zip has no assembly_data_report.jsonl")
        heldout.REPORT_PATH.write_bytes(archive.read(report_name))
        for row in rows:
            dest = fasta_dir / f"{row['accession']}.fna"
            if dest.is_file() and dest.stat().st_size > 0:
                continue
            genomic = [name for name in archive.namelist() if row["accession"] in name and name.endswith(".fna")]
            member = next((name for name in genomic if name.endswith("_genomic.fna")), genomic[0] if genomic else "")
            if not member:
                raise SystemExit(f"NCBI zip has no FASTA for {row['accession']}")
            dest.write_bytes(archive.read(member))
    _validate(rows, fasta_dir)


def _validate(rows: list[dict[str, str]], fasta_dir: Path) -> None:
    reports = {}
    for line in heldout.REPORT_PATH.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        record = json.loads(line)
        accession, _name, _tax = heldout._report_identity(record)
        reports[accession] = record
    manifest = [
        "accession\trepository\tdownload_url\tlocal_path\tfilename_bytes\tsha256\t"
        "download_status\tvalidation_status\tbiological_status\tnotes\tdownload_date"
    ]
    checksums = [f"{heldout._sha256(ZIP_PATH)}  {ZIP_PATH.relative_to(ROOT)}"]
    digests: dict[str, dict[str, str]] = {}
    species: dict[str, dict[str, str]] = {}
    genera: dict[str, dict[str, str]] = {}
    for row in rows:
        record = reports.get(row["accession"])
        if record is None:
            raise SystemExit(f"assembly report has no record for {row['accession']}")
        _accession, name, tax_id = heldout._report_identity(record)
        if not name.startswith(row["genus"]):
            raise SystemExit(f"{row['accession']} organism {name!r} is not genus {row['genus']}")
        opposite = {"sim": "db", "db": "sim"}[row["role"]]
        if row["species_taxid"] == species.setdefault(row["pair_id"], {}).get(opposite):
            raise SystemExit(f"{row['pair_id']} uses one species for both roles")
        species[row["pair_id"]][row["role"]] = row["species_taxid"]
        genera.setdefault(row["pair_id"], {})[row["role"]] = row["genus_taxid"]
        path = fasta_dir / f"{row['accession']}.fna"
        digest = sequence_digest(path)
        digests.setdefault(row["pair_id"], {})[row["role"]] = digest
        file_sha = heldout._sha256(path)
        checksums.append(f"{file_sha}  {path.relative_to(ROOT)}")
        manifest.append(
            "\t".join(
                [
                    row["accession"],
                    "NCBI Assembly",
                    "",
                    str(path.relative_to(ROOT)),
                    str(path.stat().st_size),
                    file_sha,
                    "success",
                    "pass",
                    "pass",
                    f"organism {name}; tax_id {tax_id}; family {row['family']}; {row['pair_mode']}",
                    heldout.date.today().isoformat(),
                ]
            )
        )
    for pair_id, pair_digests in digests.items():
        if pair_digests.get("sim") == pair_digests.get("db"):
            raise SystemExit(f"{pair_id} sim and db assemblies have the same sequence digest")
        same = genera[pair_id]["sim"] == genera[pair_id]["db"]
        mode = next(row["pair_mode"] for row in rows if row["pair_id"] == pair_id)
        if same and mode != "same_genus":
            raise SystemExit(f"{pair_id} shares a genus but is marked {mode}")
        if not same and mode != "different_genera":
            raise SystemExit(f"{pair_id} uses two genera but is marked {mode}")
    manifest_path = ROOT / "data" / "manifests" / "high100_download_manifest.tsv"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text("\n".join(manifest) + "\n", encoding="utf-8")
    checksum_path = ROOT / "data" / "checksums" / "high100_checksums.txt"
    checksum_path.parent.mkdir(parents=True, exist_ok=True)
    checksum_path.write_text("\n".join(checksums) + "\n", encoding="utf-8")


heldout.stage_download = stage_download


if __name__ == "__main__":
    heldout.main()
