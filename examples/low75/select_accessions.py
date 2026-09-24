"""Pin two RefSeq assemblies from each of 25 genera in three domains.

Reads an NCBI RefSeq ``assembly_summary.txt`` and an NCBI taxdump. A genus
qualifies when at least two species have a latest full assembly at Complete
Genome or Chromosome level that is not excluded from RefSeq. Genera are
ranked by how many such species they have, then by genus taxid. The first 25
genera are kept. Within a genus the two lowest species taxids are kept. The
lower species taxid is the metagenome assembly and the higher one is the
database assembly. One assembly is kept per species: reference genome, then
representative genome, then the longest remaining assembly, then the lowest
accession.

The written table is the pin. Re-running this script against a newer summary
can select different accessions.
"""

from __future__ import annotations

import argparse
from pathlib import Path

LEVELS = {"Complete Genome", "Chromosome"}
DOMAINS = ("bacteria", "archaea", "viral")
HEADER = (
    "pair_id\trole\tdomain\tgenus\tgenus_taxid\tspecies_taxid\taccession\t"
    "organism_name\ttax_id\tstrain\trefseq_category\ttotal_sequence_length\t"
    "assembly_level\tnote"
)


def load_taxdump(nodes_path: Path, names_path: Path) -> tuple[dict[int, int], dict[int, str], dict[int, str]]:
    """Return parent, rank, and scientific name for every taxid in the dump."""
    parents: dict[int, int] = {}
    ranks: dict[int, str] = {}
    for line in nodes_path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = [part.strip() for part in line.split("|")]
        if len(parts) < 3:
            continue
        tax_id = int(parts[0])
        parents[tax_id] = int(parts[1])
        ranks[tax_id] = parts[2]
    names: dict[int, str] = {}
    for line in names_path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = [part.strip() for part in line.split("|")]
        if len(parts) >= 4 and parts[3] == "scientific name":
            names[int(parts[0])] = parts[1]
    return parents, ranks, names


def genus_of(tax_id: int, parents: dict[int, int], ranks: dict[int, str]) -> int | None:
    """Walk parents until a genus node. Missing or cyclic paths return None."""
    seen: set[int] = set()
    current = tax_id
    while current and current not in seen:
        seen.add(current)
        if ranks.get(current) == "genus":
            return current
        parent = parents.get(current)
        if parent is None or parent == current:
            return None
        current = parent
    return None


def _better(candidate: dict[str, str], current: dict[str, str]) -> bool:
    """True when ``candidate`` should replace ``current`` for the same species."""
    rank = {"reference genome": 0, "representative genome": 1}
    left = rank.get(candidate["refseq_category"], 2)
    right = rank.get(current["refseq_category"], 2)
    if left != right:
        return left < right
    if int(candidate["genome_size"]) != int(current["genome_size"]):
        return int(candidate["genome_size"]) > int(current["genome_size"])
    return candidate["accession"] < current["accession"]


def choose_domain(
    summary: Path,
    domain: str,
    parents: dict[int, int],
    ranks: dict[int, str],
    names: dict[int, str],
    n_genera: int,
) -> list[str]:
    """Return accession-table lines for ``n_genera`` genera in ``domain``."""
    by_species: dict[int, dict[str, str]] = {}
    with summary.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 26:
                continue
            record = {
                "accession": parts[0],
                "refseq_category": parts[4],
                "tax_id": parts[5],
                "species_taxid": parts[6],
                "organism_name": parts[7],
                "strain": parts[8],
                "version_status": parts[10],
                "assembly_level": parts[11],
                "genome_rep": parts[13],
                "excluded": parts[20],
                "group": parts[24],
                "genome_size": parts[25] if parts[25].isdigit() else "0",
            }
            if record["group"] != domain:
                continue
            if record["version_status"] != "latest" or record["genome_rep"] != "Full":
                continue
            if record["assembly_level"] not in LEVELS or record["excluded"] != "na":
                continue
            if not record["accession"].startswith("GCF_"):
                continue
            species = int(record["species_taxid"])
            current = by_species.get(species)
            if current is None or _better(record, current):
                by_species[species] = record
    genera: dict[int, list[int]] = {}
    for species in by_species:
        genus = genus_of(species, parents, ranks)
        if genus is None:
            continue
        genus_name = names.get(genus, "")
        organism = by_species[species]["organism_name"]
        if not genus_name or genus_name.lower().startswith("unclassified"):
            continue
        if not organism.startswith(genus_name + " "):
            continue
        genera.setdefault(genus, []).append(species)
    ranked = sorted(
        ((genus, species_ids) for genus, species_ids in genera.items() if len(species_ids) >= 2),
        key=lambda item: (-len(item[1]), item[0]),
    )[:n_genera]
    if len(ranked) < n_genera:
        raise SystemExit(f"{domain}: found {len(ranked)} genera with two species, need {n_genera}")
    lines = []
    for genus, species_ids in ranked:
        chosen = sorted(species_ids)[:2]
        genus_name = names[genus]
        for role, species in (("sim", chosen[0]), ("db", chosen[1])):
            record = by_species[species]
            pair_id = genus_name.lower().replace(" ", "_").replace("-", "_")
            strain = record["strain"].removeprefix("strain=") if record["strain"] != "na" else ""
            note = "Lower species taxid is simulated; higher species taxid is the database"
            lines.append(
                "\t".join(
                    [
                        pair_id,
                        role,
                        domain,
                        genus_name,
                        str(genus),
                        str(species),
                        record["accession"],
                        record["organism_name"],
                        record["tax_id"],
                        strain,
                        record["refseq_category"],
                        record["genome_size"],
                        record["assembly_level"],
                        note,
                    ]
                )
            )
    return lines


def main() -> None:
    parser = argparse.ArgumentParser(description="Pin low75 RefSeq accessions from assembly summaries")
    parser.add_argument("--nodes", type=Path, required=True)
    parser.add_argument("--names", type=Path, required=True)
    parser.add_argument("--bacteria", type=Path, required=True)
    parser.add_argument("--archaea", type=Path, required=True)
    parser.add_argument("--viral", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--genera", type=int, default=25)
    args = parser.parse_args()
    parents, ranks, names = load_taxdump(args.nodes, args.names)
    summaries = {"bacteria": args.bacteria, "archaea": args.archaea, "viral": args.viral}
    lines = [HEADER]
    for domain in DOMAINS:
        lines.extend(choose_domain(summaries[domain], domain, parents, ranks, names, args.genera))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {len(lines) - 1} rows to {args.output}")


if __name__ == "__main__":
    main()
