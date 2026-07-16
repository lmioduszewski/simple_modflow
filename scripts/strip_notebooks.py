"""Strip outputs from every tracked notebook (fallback for the pre-commit hook).

Tracked notebooks must never carry outputs (implementation plan 2.3 / D2);
rendered copies belong in the gitignored ``docs/_rendered/`` via
``scripts/render_notebooks.py``. Run this before committing if you don't use
pre-commit:

    python scripts/strip_notebooks.py          # strip in place
    python scripts/strip_notebooks.py --check  # dry-run (non-zero on dirty)
"""

from __future__ import annotations

import subprocess
import sys


def tracked_notebooks() -> list[str]:
    out = subprocess.run(
        ["git", "ls-files", "*.ipynb"], capture_output=True, text=True, check=True
    )
    return [line for line in out.stdout.splitlines() if line]


def main() -> int:
    notebooks = tracked_notebooks()
    if not notebooks:
        print("no tracked notebooks")
        return 0
    args = [sys.executable, "-m", "nbstripout"]
    if "--check" in sys.argv[1:]:
        args.append("--dry-run")
    result = subprocess.run(args + notebooks)
    print(f"{'checked' if '--check' in sys.argv[1:] else 'stripped'} {len(notebooks)} notebooks")
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
