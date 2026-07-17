"""Sync the vendored figs snapshot at ``src/myflopy/_vendor/figs``.

Copies the minimal module closure that provides myflopy's figs surface
(``Fig``, ``Subplot``, ``Template``, ``create_hover`` and
``figs.mpl.REPORT/Theme/get_mplfig``) from a local figs checkout, applying
the transforms that make the subtree self-contained:

- absolute ``figs.`` imports become relative;
- the three package ``__init__`` files are trimmed so the vendored copy does
  NOT drag in figs' aq-test/scaling/cross-section stacks (which would add
  bokeh, cairosvg, reportlab, svglib, svgpathtools as dependencies — the
  trimmed closure needs only myflopy's existing core deps);
- the AESI logo PNG (loaded at import time by ``layout_template.py``) is
  embedded as base64 in a generated ``_logo_data.py`` so no binary data file
  has to ship in the wheel;
- ``_vendor/README.md`` is stamped with the figs commit hash and sync date.

Idempotent: rerunning against the same figs commit reproduces the same tree.

Usage:
    python scripts/sync_vendored_figs.py [--figs-repo PATH]

The figs checkout is found via --figs-repo, the FIGS_REPO env var, or the
known per-machine default paths (do not hardcode new ones elsewhere).
"""

from __future__ import annotations

import argparse
import base64
import datetime as _dt
import os
import re
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VENDOR_DIR = ROOT / "src" / "myflopy" / "_vendor"
TARGET = VENDOR_DIR / "figs"

DEFAULT_FIGS_REPOS = [
    Path("/home/lukem/python/Figures%20-%20Templates"),  # Linux dev box
    Path(r"C:\Users\lukem\Python\Projects\figs"),  # Windows dev box
]

# Files copied verbatim, then passed through the import-rewrite table.
# The whole mpl/ subpackage ships (myflopy uses plot_cross_section at
# runtime and its deps are all myflopy core deps); the plotly side is
# trimmed to core.py only.
COPY_FILES = [
    "_datatypes.py",
    "_fig.py",
    "layout_template.py",
    "_layout_presets.py",
    "plotly/core.py",
    "mpl/__init__.py",
    "mpl/mpl.py",
    "mpl/theme.py",
    "mpl/cross_section.py",
    "mpl/cross_section_data.py",
]

# Absolute-import rewrites making the subtree self-contained. Applied to
# every copied file; each pattern must remain anchored to full statements.
IMPORT_REWRITES = [
    (re.compile(r"^import figs\.layout_template as (\w+)", re.M),
     r"from . import layout_template as \1"),
    (re.compile(r"^from figs\._layout_presets import", re.M),
     "from ._layout_presets import"),
    (re.compile(r"^from figs\.figure_transforms import", re.M),
     "from .figure_transforms import"),
]

# layout_template.py loads the logo PNG at import time via Path(__file__);
# replace the whole read block with an import of the generated module.
LOGO_BLOCK = re.compile(
    r"# Enocdes the AESI logo image.*?aesi_logo = base64\.b64encode\(file\.read\(\)\)\.decode\(\)",
    re.S,
)
LOGO_REPLACEMENT = (
    "# AESI logo, embedded at sync time (see scripts/sync_vendored_figs.py)\n"
    "from ._logo_data import AESI_LOGO_B64 as aesi_logo"
)

PLOTLY_INIT = '''"""Trimmed vendored copy: only the core Fig/Subplot/Template surface.

figs' real ``plotly/__init__`` also exposes aq-test/layout/preset/scaling
helpers whose imports require bokeh/cairosvg/reportlab/svglib/svgpathtools;
myflopy does not use them, so the vendored copy imports only ``core``.
"""

from .core import Fig, PlotlyExpressProxy, Subplot, Template, snsfig

__all__ = ["Fig", "PlotlyExpressProxy", "Subplot", "Template", "snsfig"]
'''

VENDOR_INIT = '"""Vendored third-party snapshots. Do not edit by hand."""\n'


def find_figs_repo(cli_value: str | None) -> Path:
    candidates = []
    if cli_value:
        candidates.append(Path(cli_value))
    if os.environ.get("FIGS_REPO"):
        candidates.append(Path(os.environ["FIGS_REPO"]))
    candidates.extend(DEFAULT_FIGS_REPOS)
    for candidate in candidates:
        if (candidate / "src" / "figs" / "_fig.py").exists():
            return candidate
    raise SystemExit(
        "figs checkout not found. Pass --figs-repo or set FIGS_REPO; tried: "
        + ", ".join(str(c) for c in candidates)
    )


def figs_commit(repo: Path) -> str:
    try:
        out = subprocess.run(
            ["git", "-C", str(repo), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, check=True,
        )
        return out.stdout.strip()
    except Exception:
        return "unknown"


def rewrite_imports(text: str) -> str:
    for pattern, replacement in IMPORT_REWRITES:
        text = pattern.sub(replacement, text)
    return text


def build_root_init(source: Path) -> str:
    """Derive the vendored root __init__ from figs' by dropping aq_test."""

    text = source.read_text(encoding="utf-8")
    text = text.replace("from .plotly.aq_test import AquiferTestFigure\n", "")
    text = text.replace('    "AquiferTestFigure",\n', "")
    header = (
        '"""Vendored snapshot of figs (see myflopy/_vendor/README.md).\n\n'
        "Trimmed: the aq-test/scaling/cross-section stacks are omitted.\n"
        '"""\n\n'
    )
    # Replace the original module docstring with the vendored header.
    text = re.sub(r'^""".*?"""\n', header, text, count=1, flags=re.S)
    return rewrite_imports(text)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figs-repo", default=None)
    args = parser.parse_args()

    repo = find_figs_repo(args.figs_repo)
    src = repo / "src" / "figs"
    commit = figs_commit(repo)
    today = _dt.date.today().isoformat()

    if TARGET.exists():
        shutil.rmtree(TARGET)
    (TARGET / "plotly").mkdir(parents=True)
    (TARGET / "mpl").mkdir(parents=True)

    for rel in COPY_FILES:
        text = (src / rel).read_text(encoding="utf-8")
        text = rewrite_imports(text)
        if rel == "layout_template.py":
            new_text, count = LOGO_BLOCK.subn(LOGO_REPLACEMENT, text)
            if count != 1:
                raise SystemExit(
                    "layout_template.py logo block not found — figs changed shape; "
                    "update LOGO_BLOCK in scripts/sync_vendored_figs.py"
                )
            text = new_text
        (TARGET / rel).write_text(text, encoding="utf-8")

    logo_b64 = base64.b64encode(
        (src / "examples" / "data" / "aesi_logo.png").read_bytes()
    ).decode()
    (TARGET / "_logo_data.py").write_text(
        '"""AESI logo (base64), generated by scripts/sync_vendored_figs.py."""\n\n'
        f'AESI_LOGO_B64 = "{logo_b64}"\n',
        encoding="utf-8",
    )

    (TARGET / "__init__.py").write_text(build_root_init(src / "__init__.py"), encoding="utf-8")
    (TARGET / "plotly" / "__init__.py").write_text(PLOTLY_INIT, encoding="utf-8")
    (VENDOR_DIR / "__init__.py").write_text(VENDOR_INIT, encoding="utf-8")
    (VENDOR_DIR / "README.md").write_text(
        f"""# Vendored dependencies

## figs

Vendored snapshot of the local `figs` project @ commit `{commit}`, synced
{today}. Do not edit by hand; rerun `python scripts/sync_vendored_figs.py`
(against the figs checkout) after figs changes and commit the result.

This is a TRIMMED closure: the full `mpl/` subpackage (myflopy uses
`plot_cross_section` at runtime) plus the modules providing `Fig`,
`Subplot`, `Template`, `create_hover`. The plotly-side aq-test and
scaling/export stacks are deliberately omitted — they would add
bokeh/cairosvg/reportlab/svglib/svgpathtools as runtime deps. The trimmed
closure requires only myflopy's existing core dependencies (plotly, pandas,
numpy, matplotlib, seaborn, geopandas, shapely). The AESI logo is embedded
as base64 in `figs/_logo_data.py` instead of shipping a PNG.

`myflopy.viz` imports the REAL figs first and falls back to this snapshot,
so the author's machine exercises live figs while installed environments
(and CI) exercise the vendored copy.
""",
        encoding="utf-8",
    )
    print(f"vendored figs @ {commit} ({today}) -> {TARGET}")


if __name__ == "__main__":
    main()
