"""Mandatory checks for the held-out community helpers."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

pytestmark = pytest.mark.mandatory


def _community(repo_root: Path):
    path = repo_root / "examples" / "heldout_genera" / "community.py"
    spec = importlib.util.spec_from_file_location("heldout_community", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_accession_table_has_forty_distinct_genomes(repo_root: Path) -> None:
    """Four genera, five pairs, and sim keys that Samovar will not collapse."""
    community = _community(repo_root)
    rows = community.load_pairs(repo_root / "examples" / "heldout_genera" / "accessions.tsv")
    assert len(rows) == 40
    genera = {row["genus"] for row in rows}
    assert genera == {"Escherichia", "Shigella", "Pseudomonas", "Streptococcus"}
    sim = [row for row in rows if row["role"] == "sim"]
    db = [row for row in rows if row["role"] == "db"]
    assert len(sim) == 20 and len(db) == 20
    assert community.iss_key("GCF_000005845.2") == "GCF_000005845"
    non_coli = {row["pair_id"] for row in sim if row["genus"] == "Escherichia" and "coli" not in row["organism_name"].lower()}
    assert non_coli == {"eco_albertii", "eco_fergusonii", "eco_marmotae", "eco_ruysiae"}


def test_lognormal_counts_are_uneven_and_cover_every_genome() -> None:
    """One hundred thousand reads stay positive, sum exactly, and are not uniform."""
    community = _community(Path(__file__).resolve().parents[1])
    counts = community.lognormal_read_counts(20, 100_000, 42, mu=0.0, sigma=1.5)
    assert sum(counts) == 100_000
    assert min(counts) >= 1
    assert max(counts) > 2 * min(counts)
    assert counts == community.lognormal_read_counts(20, 100_000, 42, mu=0.0, sigma=1.5)
    with pytest.raises(ValueError):
        community.lognormal_read_counts(5, 4, 1)


def test_fastg_links_and_megahit_coverage() -> None:
    """FASTG neighbour lists become directed edges, including reverse links."""
    community = _community(Path(__file__).resolve().parents[1])
    text = ">NODE_1:NODE_2,NODE_3';\nACGT\n>NODE_2:NODE_1;\nTT\n"
    edges, sequences = community.parse_fastg(text)
    assert sequences == {"NODE_1": "ACGT", "NODE_2": "TT"}
    assert edges[0]["target"] == "NODE_2" and edges[0]["orientation"] == "++"
    assert edges[1]["target"] == "NODE_3" and edges[1]["orientation"] == "+-"
    assert community.megahit_coverage(">k141_1 flag=1 multi=2.5 len=9") == pytest.approx(2.5)
    assert community.megahit_coverage(">NODE_1_length_342_cov_2.0000_ID_1;") == pytest.approx(2.0)
    assert community.megahit_coverage(">plain") is None
    names = ["k99.final.contigs.fa", "k141.contigs.fa", "k21.contigs.fa", "k99.contigs.fa"]
    assert community.highest_megahit_contigs(names).endswith("k141.contigs.fa")
    edges, sequences = community.assembly_graph_from_fastg(
        ">NODE_1_length_4_cov_2.0_ID_1:NODE_2_length_4_cov_1.0_ID_2';\nACGT\n"
        ">NODE_1_length_4_cov_2.0_ID_1';\nACGT\n"
        ">NODE_2_length_4_cov_1.0_ID_2;\nTTTT\n"
    )
    assert set(sequences) == {"NODE_1_length_4_cov_2.0_ID_1", "NODE_2_length_4_cov_1.0_ID_2"}
    assert len(edges) == 1
    assert edges[0]["target"] == "NODE_2_length_4_cov_1.0_ID_2"


def test_species_walk_skips_strain_ranks() -> None:
    """A strain taxid rolls to species; a genus does not."""
    community = _community(Path(__file__).resolve().parents[1])
    parents = {511145: 562, 562: 561, 561: 543}
    ranks = {511145: "no rank", 562: "species", 561: "genus", 543: "family"}
    assert community.species_taxon(511145, parents, ranks) == 562
    assert community.species_taxon(562, parents, ranks) == 562
    assert community.species_taxon(561, parents, ranks) is None


def test_half_strains_keeps_every_other_pair(repo_root: Path) -> None:
    """Ten strains, both assemblies each, taken from the pinned table order."""
    community = _community(repo_root)
    full = community.load_pairs(repo_root / "examples" / "heldout_genera" / "accessions.tsv")
    half = community.half_strains(full)
    pinned = community.load_pairs(repo_root / "examples" / "half_strains" / "accessions.tsv")
    assert pinned == half
    assert len(half) == 20
    pair_ids = []
    for row in half:
        if row["pair_id"] not in pair_ids:
            pair_ids.append(row["pair_id"])
    assert pair_ids == [
        "eco_albertii",
        "eco_marmotae",
        "eco_coli",
        "shi_sonnei",
        "shi_dysenteriae",
        "pse_aeruginosa",
        "pse_fluorescens",
        "pse_protegens",
        "str_pyogenes",
        "str_mutans",
    ]
    assert {row["genus"] for row in half} == {"Escherichia", "Shigella", "Pseudomonas", "Streptococcus"}
    with pytest.raises(ValueError):
        community.half_strains([{"pair_id": "only", "role": "sim"}])


def test_low75_pins_two_species_in_each_genus(repo_root: Path) -> None:
    """25 genera in each domain, and low75half is every other genus."""
    community = _community(repo_root)
    rows = community.load_pairs(repo_root / "examples" / "low75" / "accessions.tsv")
    assert len(rows) == 150
    domains = {}
    species: dict[str, set[str]] = {}
    for row in rows:
        domains.setdefault(row["domain"], set()).add(row["pair_id"])
        species.setdefault(row["pair_id"], set()).add(row["species_taxid"])
    assert {name: len(pairs) for name, pairs in domains.items()} == {"bacteria": 25, "archaea": 25, "viral": 25}
    assert all(len(ids) == 2 for ids in species.values())
    half = community.load_pairs(repo_root / "examples" / "low75half" / "accessions.tsv")
    assert half == community.every_other_strain(rows)
    assert len({row["pair_id"] for row in half}) == 38


def test_presence_f1_and_r_squared() -> None:
    """A missed truth taxon lowers F1, and R² is the squared Pearson value."""
    community = _community(Path(__file__).resolve().parents[1])
    scores = community.presence_f1({1: 0.5, 3: 0.1}, {1: 0.4, 2: 0.6})
    assert scores["tp"] == 1
    assert scores["fp"] == 1
    assert scores["fn"] == 1
    assert scores["f1"] == pytest.approx(0.5)
    assert community.r_squared(0.5) == pytest.approx(0.25)
    assert community.r_squared(None) is None
    assert community.presence_f1({}, {1: 1.0})["f1"] == 0.0
