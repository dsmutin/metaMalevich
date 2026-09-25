"""Mandatory: ToCUMG colour layers are selected from MetaMetro, not recoloured here."""

from __future__ import annotations

from pathlib import Path

import pytest

from metamalevich.tocumg import load_selected

pytestmark = pytest.mark.mandatory


def test_load_selected_keeps_one_namespace(tmp_path: Path) -> None:
    """composition_kmeans is a MetaMetro colouring, not a metamalevich palette."""
    from metametro.bench.build import build

    result = build("bubble_strain_2", outdir=tmp_path / "bench")
    graph = load_selected(result.outdir, ["composition_kmeans"])
    assert {row["namespace"] for row in graph.colors or []} == {"composition"}
    assert graph.unitigs
