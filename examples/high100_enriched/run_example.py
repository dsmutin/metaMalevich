"""Build this community in MetaMetro.

The accession pin and the generator live in MetaMetro. This script forwards
``--stage`` to that generator. The ``metametro`` package is a pinned dependency of metamalevich.
"""

from __future__ import annotations


if __name__ == "__main__":
    from metametro.bench.legacy import main_for

    raise SystemExit(main_for('4domain_family_100_x10'))
