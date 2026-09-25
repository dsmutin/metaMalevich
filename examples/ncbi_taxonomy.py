"""NCBI taxdump and Kraken2 locations for example scripts.

Paths come from the environment or from ``PATH``. A missing file stops the
script. Nothing here names a machine directory.
"""

from __future__ import annotations

import importlib.util
import os
import shutil
from pathlib import Path


def nodes_dmp() -> Path:
    """Return ``nodes.dmp``. ``TAXDUMP`` is that file or the directory that holds it."""
    raw = os.environ.get("TAXDUMP", "")
    if not raw:
        raise SystemExit("TAXDUMP is unset. Point it at nodes.dmp or the directory that contains nodes.dmp and names.dmp.")
    path = Path(raw)
    if path.is_dir():
        path = path / "nodes.dmp"
    names = path.parent / "names.dmp"
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f"missing or empty nodes.dmp: {path}")
    if not names.is_file() or names.stat().st_size == 0:
        raise SystemExit(f"missing or empty names.dmp: {names}")
    return path


def taxonomy_engine() -> Path:
    """Return Samovar ``taxonomy_engine.py``. ``SAMOVAR_SRC`` is the checkout or its src tree."""
    raw = os.environ.get("SAMOVAR_SRC", "")
    if not raw:
        raise SystemExit("SAMOVAR_SRC is unset. Point it at the Samovar checkout or at the directory that contains taxonomy_engine.py.")
    root = Path(raw)
    candidates = (
        root / "samovar" / "taxonomy_engine.py",
        root / "src" / "samovar" / "taxonomy_engine.py",
        root / "taxonomy_engine.py",
    )
    for path in candidates:
        if path.is_file():
            return path
    raise SystemExit(f"taxonomy_engine.py was not found under SAMOVAR_SRC={root}")


def kraken2_bin() -> Path:
    """Return the ``kraken2`` executable on ``PATH``."""
    found = shutil.which("kraken2")
    if not found:
        raise SystemExit("kraken2 is not on PATH")
    return Path(found)


def ncbi_parser():
    """Samovar NCBI walker. A missing rank returns None."""
    engine = taxonomy_engine()
    nodes = nodes_dmp()
    spec = importlib.util.spec_from_file_location("taxonomy_engine", engine)
    if spec is None or spec.loader is None:
        raise SystemExit(f"could not load {engine}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.NCBITaxonomyParser(str(nodes))
