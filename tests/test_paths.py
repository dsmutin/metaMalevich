"""Mandatory checks that installs do not depend on this machine."""

from __future__ import annotations

import re
from importlib import metadata
from pathlib import Path

import pytest

from metamalevich.native import tool_source

pytestmark = pytest.mark.mandatory


def test_default_dataset_is_samovar10_only() -> None:
    """samovar10_ont1b repeats samovar10 and is not scored again."""
    from metamalevich.bench import DATASETS

    assert list(DATASETS) == ["samovar10"]


def test_tool_source_finds_the_repository_sources() -> None:
    """C++ sources ship with the repository, not only with one checkout path."""
    source = tool_source("kraken_count")
    assert source.name == "kraken_count.cpp"
    assert source.is_file()


def test_installed_metametro_matches_the_pin() -> None:
    """MetaMetro is a declared dependency with an exact version, and the installed one matches it."""
    requirements = [req for req in metadata.requires("metamalevich") or [] if req.startswith("metametro")]
    assert len(requirements) == 1
    pinned = re.fullmatch(r"metametro==(\S+)", requirements[0].replace(" ", ""))
    assert pinned, requirements[0]
    assert metadata.version("metametro") == pinned.group(1)
