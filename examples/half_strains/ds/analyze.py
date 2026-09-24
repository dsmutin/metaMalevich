"""Write the half-strain error charts.

The shared analysis lives in ``examples/heldout_genera/ds/analyze.py``.
This wrapper keeps a runnable script inside this example's ``ds/`` folder.
"""

from __future__ import annotations

import runpy
from pathlib import Path

TARGET = Path(__file__).resolve().parents[2] / "heldout_genera" / "ds" / "analyze.py"


def main() -> None:
    """Run the shared analysis. Charts for this example land in this folder."""
    runpy.run_path(str(TARGET), run_name="__main__")


if __name__ == "__main__":
    main()
