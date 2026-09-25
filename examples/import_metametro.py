"""Import MetaMetro from the pinned checkout or from the environment.

The pin is the submodule ``external/MetaMetro`` (https://github.com/dsmutin/MetaMetro).
``METAMETRO_SRC`` overrides that checkout when it points at a directory that
contains ``metametro/__init__.py``.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path


def ensure() -> None:
    """Put MetaMetro on ``sys.path`` when it is not already importable."""
    try:
        import metametro.bench.legacy  # noqa: F401

        return
    except ImportError:
        pass
    candidates: list[Path] = []
    env = os.environ.get("METAMETRO_SRC", "")
    if env:
        candidates.append(Path(env))
    here = Path(__file__).resolve()
    for parent in here.parents:
        candidates.append(parent / "external" / "MetaMetro" / "src")
    for candidate in candidates:
        if (candidate / "metametro" / "__init__.py").is_file():
            src = str(candidate)
            if src not in sys.path:
                sys.path.insert(0, src)
            return
    raise SystemExit(
        "MetaMetro is not importable. Clone https://github.com/dsmutin/MetaMetro "
        "as external/MetaMetro (the submodule) or set METAMETRO_SRC to its src directory."
    )
