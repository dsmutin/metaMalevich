"""Mandatory: solve reads a MetaMetro benchbuild graph and does not score truth."""

from __future__ import annotations

from pathlib import Path

import pytest

from metamalevich.solve import solve_benches

pytestmark = pytest.mark.mandatory


def test_solve_bubble_uses_metametro_colourings(tmp_path: Path) -> None:
    """composition_kmeans and decaying come from the bench, then gated_neighbour runs."""
    summary = solve_benches(
        ["bubble_strain_2"],
        tmp_path / "work",
        ["composition_kmeans", "decaying"],
    )
    assert summary["ok"] is True
    assert summary["resolver"] == "gated_neighbour"
    row = summary["benches"][0]
    assert row["solved"] is True
    assert row["canonical_name"] == "bubble_strain_2"
    assert row["namespaces"] == ["composition", "decaying"]
    assert row["n_unitigs"] > 0
    assert row["n_links"] > 0
    assert row["n_coloured_unitigs"] > 0
    assert "ground_truth" not in row


def test_phage_10_contract_does_not_count_as_a_solve(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """phage_10 is phage_species_5_x10. Without an assembly the batch is incomplete."""
    monkeypatch.setattr("metamalevich.solve._phage_ready", lambda _name, _root: False)
    summary = solve_benches(["phage_10"], tmp_path / "work", ["composition_kmeans"])
    assert summary["ok"] is False
    row = summary["benches"][0]
    assert row["canonical_name"] == "phage_species_5_x10"
    assert row["build_status"] == "contract"
    assert row["solved"] is False


def test_inprocess_solve_stands_beside_a_phage_contract(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A built in-process graph can be solved while phage_10 stays a contract."""
    monkeypatch.setattr("metamalevich.solve._phage_ready", lambda _name, _root: False)
    summary = solve_benches(
        ["bubble_strain_2", "phage_10"],
        tmp_path / "work",
        ["composition_kmeans"],
    )
    assert summary["ok"] is True
    built, phage = summary["benches"]
    assert built["solved"] is True
    assert phage["canonical_name"] == "phage_species_5_x10"
    assert phage["solved"] is False
