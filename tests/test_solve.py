"""Mandatory: solve reads a MetaMetro benchbuild graph and does not score truth."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from metamalevich.cli import main
from metamalevich.solve import solve_benches
from metamalevich.taxonomy import parse_taxdump

pytestmark = pytest.mark.mandatory

# Tiny NCBI-shaped taxonomy: root > Bacteria > clade > Enterobacteriaceae >
# Escherichia > E. coli > K-12. Only for these tests.
_NODES = [
    (1, 1, "no rank"),
    (2, 1, "superkingdom"),
    (1236, 2, "clade"),
    (543, 1236, "family"),
    (561, 543, "genus"),
    (562, 561, "species"),
    (83333, 562, "strain"),
]
_NAMES = [
    (1, "root", "scientific name"),
    (2, "Bacteria", "scientific name"),
    (2, "eubacteria", "genbank common name"),
    (1236, "Gammaproteobacteria-like clade", "scientific name"),
    (543, "Enterobacteriaceae", "scientific name"),
    (561, "Escherichia", "scientific name"),
    (562, "Escherichia coli", "scientific name"),
    (83333, "Escherichia coli K-12", "scientific name"),
]
# bubble_strain_2 CFA node ids are S, A, B, T. T stays unclassified.
_BUBBLE_KRAKEN = "C\tS\t562\t100\t562:20\nC\tA\t562\t100\t562:20\nC\tB\t561\t100\t561:20\nU\tT\t0\t100\t0:20\n"


@pytest.fixture
def taxdump(tmp_path: Path) -> Path:
    directory = tmp_path / "taxdump_test"
    directory.mkdir()
    (directory / "nodes.dmp").write_text(
        "".join(f"{taxon}\t|\t{parent}\t|\t{rank}\t|\t\t|\n" for taxon, parent, rank in _NODES),
        encoding="utf-8",
    )
    (directory / "names.dmp").write_text(
        "".join(f"{taxon}\t|\t{name}\t|\t\t|\t{kind}\t|\n" for taxon, name, kind in _NAMES),
        encoding="utf-8",
    )
    return directory


@pytest.fixture
def bubble_kraken(tmp_path: Path) -> Path:
    path = tmp_path / "bubble_kraken.tsv"
    path.write_text(_BUBBLE_KRAKEN, encoding="utf-8")
    return path


def test_parse_taxdump_links_ranks_and_names(taxdump: Path) -> None:
    """nodes.dmp gives parents and Kraken rank codes; names.dmp gives scientific names."""
    taxonomy = parse_taxdump(taxdump)
    assert taxonomy.version == "taxdump_test"
    assert taxonomy.parent(1) is None
    assert taxonomy.parent(562) == 561
    assert taxonomy.name(2) == "Bacteria"
    assert taxonomy.name(562) == "Escherichia coli"
    assert [taxonomy.rank(t) for t in (1, 2, 1236, 543, 561, 562, 83333)] == [
        "R", "D", "D1", "F", "G", "S", "S1"
    ]
    assert taxonomy.ancestor_at_rank(83333, "S") == 562
    assert taxonomy.ancestor_at_rank(83333, "G") == 561
    assert taxonomy.lca([83333, 561]) == 561


def test_parse_taxdump_stops_without_names(tmp_path: Path, taxdump: Path) -> None:
    (taxdump / "names.dmp").unlink()
    with pytest.raises(FileNotFoundError, match="names.dmp"):
        parse_taxdump(taxdump)


def test_solve_bubble_reads_taxa_from_kraken2(tmp_path: Path, taxdump: Path, bubble_kraken: Path) -> None:
    """The default layers include kraken2; only its taxon ids become resolver taxa."""
    summary = solve_benches(
        ["bubble_strain_2"], tmp_path / "work", kraken_output=bubble_kraken, taxdump=taxdump
    )
    assert summary["ok"] is True
    assert summary["status"] == "ok"
    assert summary["resolver"] == "gated_neighbour"
    assert summary["colourings"] == ["kraken2", "composition_kmeans", "decaying"]
    assert summary["taxdump"] == "taxdump_test"
    row = summary["benches"][0]
    assert row["solved"] is True
    assert row["canonical_name"] == "bubble_strain_2"
    assert "kraken2" in row["namespaces"]
    assert row["n_taxa"] == 2
    assert row["n_coloured_unitigs"] > 0
    assert "ground_truth" not in row


def test_solve_without_a_taxonomic_layer_fails(tmp_path: Path, taxdump: Path, bubble_kraken: Path) -> None:
    with pytest.raises(ValueError, match="taxonomic colouring kraken2"):
        solve_benches(
            ["bubble_strain_2"],
            tmp_path / "work",
            ["composition_kmeans", "decaying"],
            kraken_output=bubble_kraken,
            taxdump=taxdump,
        )


def test_solve_with_kraken_needs_taxdump(
    tmp_path: Path, bubble_kraken: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.delenv("TAXDUMP", raising=False)
    with pytest.raises(ValueError, match="TAXDUMP"):
        solve_benches(["bubble_strain_2"], tmp_path / "work", kraken_output=bubble_kraken)


def test_solve_stops_on_a_taxid_outside_taxdump(tmp_path: Path, taxdump: Path) -> None:
    kraken = tmp_path / "unknown.tsv"
    kraken.write_text("C\tS\t999999\t100\t999999:20\n", encoding="utf-8")
    with pytest.raises(ValueError, match="not in TAXDUMP"):
        solve_benches(["bubble_strain_2"], tmp_path / "work", kraken_output=kraken, taxdump=taxdump)


def test_solve_skips_a_bench_without_kraken_input(tmp_path: Path) -> None:
    """Until pinned kraken2 outputs exist, a bench without one is skipped with a reason."""
    summary = solve_benches(["bubble_strain_2", "synthetic_reads_2"], tmp_path / "work")
    assert summary["status"] == "skipped"
    assert summary["ok"] is True
    for row in summary["benches"]:
        assert row["build_status"] == "skipped"
        assert row["solved"] is False
        assert "kraken2" in row["skipped"]
    assert not (tmp_path / "work" / "bubble_strain_2").exists()


def test_cli_solve_skip_exits_zero(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    code = main(["solve", "--bench", "bubble_strain_2", "--workdir", str(tmp_path / "work")])
    summary = json.loads(capsys.readouterr().out)
    assert code == 0
    assert summary["status"] == "skipped"


def test_phage_10_contract_does_not_count_as_a_solve(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """phage_10 is phage_species_5_x10. Without an assembly the batch is incomplete."""
    monkeypatch.setattr("metamalevich.solve._phage_ready", lambda _name, _root: False)
    summary = solve_benches(["phage_10"], tmp_path / "work")
    assert summary["ok"] is False
    row = summary["benches"][0]
    assert row["canonical_name"] == "phage_species_5_x10"
    assert row["build_status"] == "contract"
    assert row["solved"] is False


def test_inprocess_solve_stands_beside_a_phage_contract(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, taxdump: Path, bubble_kraken: Path
) -> None:
    """A built in-process graph can be solved while phage_10 stays a contract."""
    monkeypatch.setattr("metamalevich.solve._phage_ready", lambda _name, _root: False)
    summary = solve_benches(
        ["bubble_strain_2", "phage_10"],
        tmp_path / "work",
        kraken_output=bubble_kraken,
        taxdump=taxdump,
    )
    assert summary["ok"] is True
    built, phage = summary["benches"]
    assert built["solved"] is True
    assert phage["canonical_name"] == "phage_species_5_x10"
    assert phage["solved"] is False
