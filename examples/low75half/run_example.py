"""Build this community in MetaMetro.

The accession pin and the generator live in MetaMetro. This script forwards
``--stage`` to that generator. Set ``METAMETRO_SRC`` to the MetaMetro ``src``
directory when the package is not already importable.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path


def _import_metametro() -> None:
    try:
        import metametro.bench.legacy  # noqa: F401
        return
    except ImportError:
        pass
    env = os.environ.get("METAMETRO_SRC", "")
    candidates = []
    if env:
        candidates.append(Path(env))
    here = Path(__file__).resolve()
    candidates.extend(parent / "metametro" / "src" for parent in here.parents)
    for candidate in candidates:
        if (candidate / "metametro" / "bench").is_dir():
            sys.path.insert(0, str(candidate))
            return
    raise SystemExit("metametro is not importable and METAMETRO_SRC is unset")


if __name__ == "__main__":
    _import_metametro()
    from metametro.bench.legacy import main_for

    raise SystemExit(main_for("3domain_genus_75_half"))
