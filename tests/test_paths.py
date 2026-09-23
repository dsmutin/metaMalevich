"""Mandatory checks that installs do not depend on this machine."""

from __future__ import annotations

from pathlib import Path

import pytest

from metamalevich.bridge import metametro_src
from metamalevich.native import tool_source

pytestmark = pytest.mark.mandatory


def test_tool_source_finds_the_repository_sources() -> None:
    """C++ sources ship with the repository, not only with one checkout path."""
    source = tool_source("kraken_count")
    assert source.name == "kraken_count.cpp"
    assert source.is_file()


def test_metametro_search_rejects_a_missing_checkout(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A missing MetaMetro checkout fails with setup instructions and no machine path."""
    monkeypatch.delenv("METAMETRO_SRC", raising=False)
    monkeypatch.chdir(tmp_path)
    with pytest.raises(FileNotFoundError, match="external/MetaMetro") as caught:
        metametro_src(tmp_path)
    assert "/mnt/" not in str(caught.value)
