"""Execute tracked notebooks into gitignored ``docs/_rendered/`` copies.

Tracked notebooks are stored output-free (implementation plan 2.3 / D2);
this produces browsable executed HTML copies outside git:

    python scripts/render_notebooks.py                       # all notebooks
    python scripts/render_notebooks.py canonical_00 layer    # name filters

Requires jupyter/nbconvert (``pip install nbconvert jupyter``); notebooks
that need local data or MF6 binaries will fail individually and are
reported without stopping the batch.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "_rendered"


def main() -> int:
    filters = [arg for arg in sys.argv[1:] if not arg.startswith("-")]
    out = subprocess.run(
        ["git", "ls-files", "*.ipynb"], capture_output=True, text=True, check=True, cwd=ROOT
    )
    notebooks = [
        ROOT / line
        for line in out.stdout.splitlines()
        if line and (not filters or any(f in line for f in filters))
    ]
    OUT.mkdir(parents=True, exist_ok=True)
    failures = []
    for nb in notebooks:
        print(f"rendering {nb.relative_to(ROOT)} ...", flush=True)
        result = subprocess.run(
            [
                sys.executable, "-m", "nbconvert", "--to", "html", "--execute",
                "--output-dir", str(OUT), str(nb),
            ],
            cwd=nb.parent,
        )
        if result.returncode != 0:
            failures.append(str(nb.relative_to(ROOT)))
    print(f"rendered {len(notebooks) - len(failures)}/{len(notebooks)} into {OUT}")
    if failures:
        print("failed:", *failures, sep="\n  ")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
