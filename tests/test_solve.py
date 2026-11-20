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
    assert summary["resolver"] == "genus_plurality"
    assert summary["colourings"] == ["kraken2", "composition_kmeans", "decaying"]
    assert summary["taxdump"] == "taxdump_test"
    assert summary["taxdump_source"] == "cli"
    row = summary["benches"][0]
    assert row["kraken2_input"] == "given"
    assert row["taxdump_source"] == "cli"
    assert "kraken_pin_index" not in row
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
        assert "MetaMetro pins no Kraken2 output" in row["skipped"]
        assert "#1163" not in row["skipped"]
    assert not (tmp_path / "work" / "bubble_strain_2").exists()


def test_cli_solve_skip_exits_zero(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    code = main(["solve", "--bench", "bubble_strain_2", "--workdir", str(tmp_path / "work")])
    summary = json.loads(capsys.readouterr().out)
    assert code == 0
    assert summary["status"] == "skipped"


def _phage_not_ready(monkeypatch: pytest.MonkeyPatch, genome_errors: list[str], program_errors: list[str]) -> None:
    """Make MetaMetro's phage preflight checks raise these errors. No network, no programs."""
    import importlib

    from metametro.errors import ContractError

    genomes = importlib.import_module("metametro.bench.genomes")
    catalog = importlib.import_module("metametro.bench.data.universal.catalog")

    def check_genomes(_bench, _genome_dir):
        if genome_errors:
            raise ContractError(genome_errors)

    def require_programs(_spec, _names):
        if program_errors:
            raise ContractError(program_errors)

    monkeypatch.setattr(genomes, "check_genomes", check_genomes)
    monkeypatch.setattr(catalog, "require_programs", require_programs)


_GENOME_ERRORS = [
    "phage_species_5_x10 is missing pinned genome FASTA under data/raw/genomes: T1 (NC_005833.1)",
    "fetch them with: metametro benchfetch phage_species_5_x10",
]
_PROGRAM_ERRORS = ["phage_species_5_x10 needs these programs on PATH before a build can start: megahit"]


def test_phage_10_contract_does_not_count_as_a_solve(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """phage_10 is phage_species_5_x10. Without an assembly the batch is incomplete."""
    _phage_not_ready(monkeypatch, _GENOME_ERRORS, _PROGRAM_ERRORS)
    summary = solve_benches(["phage_10"], tmp_path / "work")
    assert summary["ok"] is False
    assert summary["status"] == "incomplete"
    row = summary["benches"][0]
    assert row["canonical_name"] == "phage_species_5_x10"
    assert row["build_status"] == "contract"
    assert row["solved"] is False


def test_phage_missing_comes_from_metametro_errors(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The missing list is what MetaMetro raised: genomes with the fetch command, then programs."""
    _phage_not_ready(monkeypatch, _GENOME_ERRORS, _PROGRAM_ERRORS)
    row = solve_benches(["phage_10"], tmp_path / "work")["benches"][0]
    assert row["missing"] == _GENOME_ERRORS + _PROGRAM_ERRORS


def test_phage_missing_genomes_name_benchfetch(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """MetaMetro's real genome check on an empty checkout names the benchfetch command."""
    import importlib

    paths = importlib.import_module("metametro.bench.paths")
    monkeypatch.setattr(paths, "repo_root", lambda: tmp_path / "empty_checkout")
    row = solve_benches(["phage_10"], tmp_path / "work")["benches"][0]
    assert row["build_status"] == "contract"
    assert any("metametro benchfetch phage_species_5_x10" in error for error in row["missing"])


def _kraken_pin(monkeypatch: pytest.MonkeyPatch, present: bool) -> None:
    """Make MetaMetro report a Kraken2 pin (or none) for every bench."""
    import importlib

    module = importlib.import_module("metametro.bench.kraken_pin")
    monkeypatch.setattr(
        module, "kraken_pin", lambda _bench, **_kw: (Path("nodes.kraken"), {}) if present else None
    )


def test_phage_10_has_a_metametro_kraken_pin() -> None:
    """MetaMetro v0.22.0 pins Kraken2 output for phage_10 only among these benches."""
    from metametro.bench.kraken_pin import kraken_pin

    assert kraken_pin("phage_10") is not None
    assert kraken_pin("bubble_strain_2") is None
    assert kraken_pin("synthetic_reads_2") is None


def test_pinned_phage_contract_is_not_skipped_for_kraken(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """With the pin, phage_10 is not a missing-kraken skip; unready, it stays a contract."""
    _phage_not_ready(monkeypatch, [], _PROGRAM_ERRORS)
    monkeypatch.delenv("TAXDUMP", raising=False)
    row = solve_benches(["phage_10"], tmp_path / "work")["benches"][0]
    assert row["build_status"] == "contract"
    assert row["kraken2_input"] == "metametro_pin"
    assert "skipped" not in row


def _stop_at_build(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """Replace MetaMetro build() so a test stops right before any assembly."""
    import importlib

    called: list[str] = []
    bench_build = importlib.import_module("metametro.bench.build")

    def build(name, **_kwargs):
        called.append(name)
        raise RuntimeError("stop before build")

    monkeypatch.setattr(bench_build, "build", build)
    return called


def test_pinned_taxdump_is_the_phage_10_taxonomy() -> None:
    """MetaMetro v0.22.1 pins a trimmed TAXDUMP beside the phage_10 Kraken2 output."""
    from metamalevich.solve import _pinned_taxonomy

    taxonomy = _pinned_taxonomy("phage_10")
    assert taxonomy.version == "k2_viral_20261022-phage_species_5_x10"
    assert taxonomy.name(10665) == "Tequatrovirus T4"
    assert taxonomy.rank(10665) == "S"


def test_ready_pinned_phage_uses_the_pinned_taxdump(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Without --taxdump or TAXDUMP, a ready phage_10 loads the pinned TAXDUMP before build."""
    _phage_not_ready(monkeypatch, [], [])
    monkeypatch.delenv("TAXDUMP", raising=False)
    loaded: list[str] = []
    real = __import__("metamalevich.solve", fromlist=["_pinned_taxonomy"])._pinned_taxonomy
    monkeypatch.setattr(
        "metamalevich.solve._pinned_taxonomy", lambda bench: loaded.append(bench) or real(bench)
    )
    called = _stop_at_build(monkeypatch)
    with pytest.raises(RuntimeError, match="stop before build"):
        solve_benches(["phage_10"], tmp_path / "work")
    assert loaded == ["phage_species_5_x10"]
    assert called == ["phage_10"]


def test_explicit_taxdump_wins_over_the_pinned_one(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, taxdump: Path
) -> None:
    _phage_not_ready(monkeypatch, [], [])

    def no_pin(_bench):
        raise AssertionError("pinned TAXDUMP must not be read when --taxdump is given")

    monkeypatch.setattr("metamalevich.solve._pinned_taxonomy", no_pin)
    called = _stop_at_build(monkeypatch)
    with pytest.raises(RuntimeError, match="stop before build"):
        solve_benches(["phage_10"], tmp_path / "work", taxdump=taxdump)
    assert called == ["phage_10"]


def _contract_build(monkeypatch: pytest.MonkeyPatch) -> None:
    """Replace MetaMetro build() by one that returns a contract without assembling."""
    import importlib
    from types import SimpleNamespace

    from metametro.bench.registry import resolve

    bench_build = importlib.import_module("metametro.bench.build")
    monkeypatch.setattr(
        bench_build,
        "build",
        lambda name, **_kw: SimpleNamespace(spec=resolve(name), status="contract", colourings=(), outdir=None),
    )


def test_explicit_kraken_disables_the_pin_without_taxdump(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, bubble_kraken: Path
) -> None:
    monkeypatch.delenv("TAXDUMP", raising=False)
    called = _stop_at_build(monkeypatch)
    with pytest.raises(
        ValueError,
        match=r"explicit --kraken-output/--kraken-db disables the MetaMetro pin for phage_species_5_x10; "
        r"pass --taxdump for your own index",
    ):
        solve_benches(["phage_10"], tmp_path / "work", kraken_output=bubble_kraken)
    assert called == []


def test_env_taxdump_beats_the_pin_and_is_recorded(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, taxdump: Path
) -> None:
    """$TAXDUMP wins over the pinned TAXDUMP; the row shows the pin index and the TAXDUMP source."""
    _phage_not_ready(monkeypatch, [], [])
    monkeypatch.setenv("TAXDUMP", str(taxdump))

    def no_pin(_bench):
        raise AssertionError("pinned TAXDUMP must not be read when TAXDUMP is set")

    monkeypatch.setattr("metamalevich.solve._pinned_taxonomy", no_pin)
    _contract_build(monkeypatch)
    summary = solve_benches(["phage_10"], tmp_path / "work")
    row = summary["benches"][0]
    assert row["kraken2_input"] == "metametro_pin"
    assert row["kraken_pin_index"] == "k2_viral_20261022"
    assert row["taxdump"] == "taxdump_test"
    assert row["taxdump_source"] == "env"
    assert summary["taxdump_source"] == "env"


def test_pinned_taxdump_is_recorded_on_the_row(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _phage_not_ready(monkeypatch, [], [])
    monkeypatch.delenv("TAXDUMP", raising=False)
    _contract_build(monkeypatch)
    row = solve_benches(["phage_10"], tmp_path / "work")["benches"][0]
    assert row["kraken_pin_index"] == "k2_viral_20261022"
    assert row["taxdump"] == "k2_viral_20261022-phage_species_5_x10"
    assert row["taxdump_source"] == "pin"


def test_pin_without_a_taxdump_is_a_row_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """If MetaMetro pinned no TAXDUMP and none is given, that bench errors before build."""
    import importlib

    _phage_not_ready(monkeypatch, [], [])
    monkeypatch.delenv("TAXDUMP", raising=False)
    module = importlib.import_module("metametro.bench.kraken_pin")
    monkeypatch.setattr(module, "pinned_taxdump", lambda _bench, **_kw: None)
    called = _stop_at_build(monkeypatch)
    summary = solve_benches(["phage_10"], tmp_path / "work")
    row = summary["benches"][0]
    assert row["build_status"] == "error"
    assert "pins no TAXDUMP" in row["error"]
    assert summary["ok"] is False
    assert called == []


def test_broken_pin_fails_only_its_bench(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, taxdump: Path, bubble_kraken: Path
) -> None:
    """A pin that MetaMetro rejects gives an error row; the other benches keep their rows."""
    import importlib

    from metametro.errors import ContractError

    module = importlib.import_module("metametro.bench.kraken_pin")
    real = module.kraken_pin

    def kraken_pin(bench, **kw):
        if bench in ("phage_10", "phage_species_5_x10"):
            raise ContractError(["nodes.kraken does not have the SHA-256 pinned in kraken2.yaml"])
        return real(bench, **kw)

    monkeypatch.setattr(module, "kraken_pin", kraken_pin)
    summary = solve_benches(["bubble_strain_2", "phage_10", "synthetic_reads_2"], tmp_path / "work")
    bubble, phage, synthetic = summary["benches"]
    assert phage["build_status"] == "error"
    assert "SHA-256" in phage["error"]
    assert bubble["build_status"] == "skipped"
    assert synthetic["build_status"] == "skipped"
    assert summary["ok"] is False
    assert summary["status"] == "incomplete"


def test_ready_phage_without_any_kraken_input_is_skipped(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Without a pin or a given input, a ready phage bench is skipped, not assembled."""
    _phage_not_ready(monkeypatch, [], [])
    _kraken_pin(monkeypatch, False)
    row = solve_benches(["phage_10"], tmp_path / "work")["benches"][0]
    assert row["build_status"] == "skipped"
    assert "MetaMetro pins no Kraken2 output" in row["skipped"]
    assert "missing" not in row
    assert not (tmp_path / "work" / "phage_species_5_x10").exists()


def test_inprocess_bench_with_a_pin_is_built(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, taxdump: Path
) -> None:
    """A pin is a kraken2 input: the bench is built instead of skipped."""
    _kraken_pin(monkeypatch, True)
    built = []
    import importlib

    bench_build = importlib.import_module("metametro.bench.build")

    def build(name, **kwargs):
        built.append(name)
        raise RuntimeError("stop after the skip decision")

    monkeypatch.setattr(bench_build, "build", build)
    with pytest.raises(RuntimeError, match="stop after the skip decision"):
        solve_benches(["bubble_strain_2"], tmp_path / "work", taxdump=taxdump)
    assert built == ["bubble_strain_2"]


def test_inprocess_solve_stands_beside_a_phage_contract(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, taxdump: Path, bubble_kraken: Path
) -> None:
    """A built in-process graph can be solved while phage_10 stays a contract."""
    _phage_not_ready(monkeypatch, _GENOME_ERRORS, [])
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
    assert phage["missing"] == _GENOME_ERRORS


def _patch_gold(monkeypatch: pytest.MonkeyPatch, nodes: dict[str, str] | None) -> None:
    """After MetaMetro builds, replace (or drop) the bench gold. Test-only gold."""
    import importlib

    bench_build = importlib.import_module("metametro.bench.build")

    real = bench_build.build

    def build(name, **kwargs):
        result = real(name, **kwargs)
        gold = Path(result.outdir) / "ground_truth"
        if nodes is None:
            for path in gold.glob("taxon_*.tsv"):
                path.unlink()
        else:
            (gold / "taxon_nodes.tsv").write_text(
                "node_id\ttaxon\n" + "".join(f"{node}\t{value}\n" for node, value in nodes.items()),
                encoding="utf-8",
            )
            (gold / "taxon_edges.tsv").write_text("edge_id\ttaxon\n", encoding="utf-8")
        return result

    monkeypatch.setattr(bench_build, "build", build)


@pytest.mark.parametrize("resolver", ["genus_plurality", "gated_neighbour", "probability_sum", "hard"])
def test_solve_runs_the_chosen_resolver(
    tmp_path: Path, taxdump: Path, bubble_kraken: Path, resolver: str
) -> None:
    summary = solve_benches(
        ["bubble_strain_2"], tmp_path / "work", kraken_output=bubble_kraken, taxdump=taxdump, resolver=resolver
    )
    assert summary["resolver"] == resolver
    assert summary["benches"][0]["solved"] is True


def test_solve_rejects_an_unknown_resolver(tmp_path: Path, taxdump: Path, bubble_kraken: Path) -> None:
    with pytest.raises(ValueError, match="unknown resolver"):
        solve_benches(
            ["bubble_strain_2"], tmp_path / "work", kraken_output=bubble_kraken, taxdump=taxdump, resolver="lca"
        )


def test_solve_writes_assignments_and_a_length_weighted_profile(
    tmp_path: Path, taxdump: Path, bubble_kraken: Path
) -> None:
    """bubble_strain_2: S (8 bp) and A (6 bp) are 562, B (6 bp) is 561, T (8 bp) is unclassified."""
    summary = solve_benches(["bubble_strain_2"], tmp_path / "work", kraken_output=bubble_kraken, taxdump=taxdump)
    row = summary["benches"][0]
    out = tmp_path / "work" / "solve" / "bubble_strain_2"
    assert row["assignments"] == str(out / "assignments.tsv")
    lines = (out / "assignments.tsv").read_text(encoding="utf-8").splitlines()
    assert lines[0].split("\t") == ["unitig_id", "length", "taxon_id", "name", "rank", "probability", "distribution"]
    by_taxon = sorted(line.split("\t")[2] for line in lines[1:])
    assert by_taxon == ["0", "561", "562", "562"]
    profile = {
        line.split("\t")[0]: float(line.split("\t")[3])
        for line in (out / "profile.tsv").read_text(encoding="utf-8").splitlines()[1:]
    }
    assert profile == pytest.approx({"562": 14 / 28, "561": 6 / 28, "0": 8 / 28}, rel=1e-5)
    assert row["profile_fractions"] == pytest.approx({"562": 0.5, "561": 6 / 28, "0": 8 / 28})


def test_bubble_gold_names_give_null_metrics_with_a_reason(
    tmp_path: Path, taxdump: Path, bubble_kraken: Path
) -> None:
    """bubble gold is taxon_1/taxon_2, not taxon ids, so nothing is mapped or scored."""
    summary = solve_benches(["bubble_strain_2"], tmp_path / "work", kraken_output=bubble_kraken, taxdump=taxdump)
    row = summary["benches"][0]
    assert row["metrics"] is None
    assert "not NCBI taxon ids" in row["metrics_reason"]


def test_missing_gold_gives_null_metrics_with_a_reason(
    tmp_path: Path, taxdump: Path, bubble_kraken: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _patch_gold(monkeypatch, None)
    summary = solve_benches(["bubble_strain_2"], tmp_path / "work", kraken_output=bubble_kraken, taxdump=taxdump)
    row = summary["benches"][0]
    assert row["metrics"] is None
    assert "no gold taxon labels" in row["metrics_reason"]


def test_taxid_gold_is_scored(
    tmp_path: Path, taxdump: Path, bubble_kraken: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """With taxon-id gold, metrics compare the resolved argmax with each unitig's gold set."""
    _patch_gold(monkeypatch, {"S": "562", "A": "562", "B": "561", "T": "561,562"})
    summary = solve_benches(["bubble_strain_2"], tmp_path / "work", kraken_output=bubble_kraken, taxdump=taxdump)
    metrics = summary["benches"][0]["metrics"]
    assert metrics["n_labelled_unitigs"] == 4
    assert metrics["n_single_label_unitigs"] == 3
    assert metrics["accuracy_any"] == pytest.approx(0.75)
    assert metrics["accuracy"] == pytest.approx(1.0)
    assert metrics["f1"] == pytest.approx(1.0)
