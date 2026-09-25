"""Build this community in MetaMetro.

The accession pin and the generator live in MetaMetro. This script forwards
``--stage`` to that generator. The package comes from the ``external/MetaMetro``
submodule. ``METAMETRO_SRC`` overrides it.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from import_metametro import ensure  # noqa: E402


if __name__ == "__main__":
    ensure()
    from metametro.bench.legacy import main_for

    raise SystemExit(main_for('3domain_genus_75_half'))
