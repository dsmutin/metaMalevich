"""Pin two RefSeq assemblies for each of 100 families.

Twenty-five families come from bacteria, archaea, viruses, and small
eukaryotes (fungi and protozoa). Thirteen families in each domain use two
species of one genus. Twelve use one species from each of two genera. The
lower species taxid is simulated. The higher species taxid is the database
genome. Eukaryotes use the smallest complete or chromosome assembly of each
species, and the smallest such families are kept.

The written table is the pin. A newer assembly summary can change it.
"""

from __future__ import annotations

import argparse
from pathlib import Path

LEVELS = {"Complete Genome", "Chromosome"}
HEADER = (
    "pair_id\trole\tdomain\tfamily\tfamily_taxid\tgenus\tgenus_taxid\t"
    "species_taxid\taccession\torganism_name\ttax_id\tstrain\trefseq_category\t"
    "total_sequence_length\tassembly_level\tpair_mode\tnote"
)
SAME_N = 13
DIFF_N = 12


def load_taxdump(nodes_path: Path, names_path: Path):
    parents, ranks, names = {}, {}, {}
    for line in nodes_path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = [part.strip() for part in line.split("|")]
        if len(parts) < 3:
            continue
        tax_id = int(parts[0])
        parents[tax_id] = int(parts[1])
        ranks[tax_id] = parts[2]
    for line in names_path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = [part.strip() for part in line.split("|")]
        if len(parts) >= 4 and parts[3] == "scientific name":
            names[int(parts[0])] = parts[1]
    return parents, ranks, names


def ancestor_rank(tax_id: int, parents: dict, ranks: dict, target: str) -> int | None:
    seen: set[int] = set()
    current = tax_id
    while current and current not in seen:
        seen.add(current)
        if ranks.get(current) == target:
            return current
        parent = parents.get(current)
        if parent is None or parent == current:
            return None
        current = parent
    return None


def _prefer(candidate: dict, current: dict, *, smallest: bool) -> bool:
    rank = {"reference genome": 0, "representative genome": 1}
    if smallest:
        if int(candidate["genome_size"]) != int(current["genome_size"]):
            return int(candidate["genome_size"]) < int(current["genome_size"])
        return candidate["accession"] < current["accession"]
    left = rank.get(candidate["refseq_category"], 2)
    right = rank.get(current["refseq_category"], 2)
    if left != right:
        return left < right
    if int(candidate["genome_size"]) != int(current["genome_size"]):
        return int(candidate["genome_size"]) > int(current["genome_size"])
    return candidate["accession"] < current["accession"]


def index_species(paths: list[Path], groups: set[str], *, smallest: bool) -> dict[int, dict]:
    by_species: dict[int, dict] = {}
    for path in paths:
        with path.open(encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 26 or parts[24] not in groups:
                    continue
                if parts[10] != "latest" or parts[13] != "Full" or parts[11] not in LEVELS or parts[20] != "na":
                    continue
                if not parts[0].startswith("GCF_") or not parts[25].isdigit():
                    continue
                record = {
                    "accession": parts[0],
                    "refseq_category": parts[4],
                    "tax_id": parts[5],
                    "species_taxid": int(parts[6]),
                    "organism_name": parts[7],
                    "strain": parts[8],
                    "assembly_level": parts[11],
                    "genome_size": parts[25],
                }
                current = by_species.get(record["species_taxid"])
                if current is None or _prefer(record, current, smallest=smallest):
                    by_species[record["species_taxid"]] = record
    return by_species


def family_table(by_species, parents, ranks, names, *, smallest: bool, domain: str) -> list[str]:
    families: dict[int, dict[int, list[int]]] = {}
    for species, record in by_species.items():
        genus = ancestor_rank(species, parents, ranks, "genus")
        family = ancestor_rank(species, parents, ranks, "family")
        if genus is None or family is None:
            continue
        genus_name = names.get(genus, "")
        family_name = names.get(family, "")
        if not genus_name or not family_name or family_name.lower().startswith("unclassified"):
            continue
        if not record["organism_name"].startswith(genus_name + " "):
            continue
        families.setdefault(family, {}).setdefault(genus, []).append(species)

    def same_pair(genera: dict[int, list[int]]) -> tuple[int, int, int] | None:
        options = [(genus, sorted(species)) for genus, species in genera.items() if len(species) >= 2]
        if not options:
            return None
        genus, species = sorted(options, key=lambda item: (-len(item[1]), item[0]))[0]
        return genus, species[0], species[1]

    def diff_pair(genera: dict[int, list[int]]) -> tuple[int, int, int, int] | None:
        if len(genera) < 2:
            return None
        first, second = sorted(genera)[:2]
        return first, min(genera[first]), second, min(genera[second])

    def size_of(species_ids: list[int]) -> int:
        return sum(int(by_species[species]["genome_size"]) for species in species_ids)

    ranked = []
    for family, genera in families.items():
        same = same_pair(genera)
        diff = diff_pair(genera)
        if same is None and diff is None:
            continue
        probe = []
        if same:
            probe.extend(same[1:])
        if diff:
            probe.extend((diff[1], diff[3]))
        key = (size_of(probe), family) if smallest else (-sum(len(v) for v in genera.values()), family)
        ranked.append((key, family, genera, same, diff))
    ranked.sort(key=lambda item: item[0])

    chosen = []
    same_left, diff_left = SAME_N, DIFF_N
    for _key, family, genera, same, diff in ranked:
        if same_left and same is not None:
            genus, left, right = same
            chosen.append((family, genera, "same_genus", [(genus, left), (genus, right)]))
            same_left -= 1
        elif diff_left and diff is not None:
            genus_a, species_a, genus_b, species_b = diff
            chosen.append((family, genera, "different_genera", [(genus_a, species_a), (genus_b, species_b)]))
            diff_left -= 1
        if same_left == 0 and diff_left == 0:
            break
    if same_left or diff_left:
        raise SystemExit(f"{domain}: short {same_left} same-genus and {diff_left} cross-genus families")
    lines = []
    for family, _genera, mode, endpoints in chosen:
        family_name = names[family]
        ordered = sorted(endpoints, key=lambda item: item[1])
        for role, (genus, species) in zip(("sim", "db"), ordered):
            record = by_species[species]
            pair_id = domain + "_" + family_name.lower().replace(" ", "_").replace("-", "_")
            strain = record["strain"].removeprefix("strain=") if record["strain"] != "na" else ""
            lines.append(
                "\t".join(
                    [
                        pair_id,
                        role,
                        domain,
                        family_name,
                        str(family),
                        names[genus],
                        str(genus),
                        str(species),
                        record["accession"],
                        record["organism_name"],
                        record["tax_id"],
                        strain,
                        record["refseq_category"],
                        record["genome_size"],
                        record["assembly_level"],
                        mode,
                        "Lower species taxid is simulated; higher species taxid is the database",
                    ]
                )
            )
    return lines


def main() -> None:
    parser = argparse.ArgumentParser(description="Pin high100 RefSeq families")
    parser.add_argument("--nodes", type=Path, required=True)
    parser.add_argument("--names", type=Path, required=True)
    parser.add_argument("--bacteria", type=Path, required=True)
    parser.add_argument("--archaea", type=Path, required=True)
    parser.add_argument("--viral", type=Path, required=True)
    parser.add_argument("--fungi", type=Path, required=True)
    parser.add_argument("--protozoa", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    parents, ranks, names = load_taxdump(args.nodes, args.names)
    domains = [
        ("bacteria", [args.bacteria], {"bacteria"}, False),
        ("archaea", [args.archaea], {"archaea"}, False),
        ("viral", [args.viral], {"viral"}, False),
        ("eukaryota", [args.fungi, args.protozoa], {"fungi", "protozoa"}, True),
    ]
    lines = [HEADER]
    for domain, paths, groups, smallest in domains:
        species = index_species(paths, groups, smallest=smallest)
        lines.extend(family_table(species, parents, ranks, names, smallest=smallest, domain=domain))
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {len(lines) - 1} rows")


if __name__ == "__main__":
    main()
