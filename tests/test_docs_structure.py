"""Structural checks on the markdown docs: fences close, local links resolve.

Neither of these is about prose quality -- they are the two ways a doc breaks
*mechanically*, silently, and only for the reader:

* **An unclosed code fence** swallows the entire rest of the file into a code
  block on GitHub. The author, reading their own diff, sees nothing wrong. A
  June 2026 review found exactly this in `preferred_api.md`.
* **A local link to a file that moved** is the standard decay product of a
  refactor. This repo has moved a lot of modules (plans 4.x, and Phase 8 will
  move the plotting ones), and every move is a chance for a doc to start
  pointing at nothing.

Deliberately NOT checked: external URLs (they need the network and rot for
reasons outside this repo), and anchors within a file (`#section`), which would
mean parsing every heading and slugifying it the way GitHub does -- a lot of
machinery for a much rarer failure.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

#: Directories whose markdown is ours to keep correct.
SCANNED = ("docs", "examples", ".")

#: Link targets that are allowed not to exist on disk.
#:
#: Generated artifacts are the whole category: a doc SHOULD tell you the
#: pamphlet exists, and 7.4 untracked it precisely because it is generated.
GENERATED = {
    "docs/myflopy_api_pamphlet.pdf",
}

#: Manual chapters that are planned but NOT YET WRITTEN.
#:
#: `docs/manual/README.md` is a full table of contents for a manual whose text
#: is, so far, chapters 3 and 4. Its links to the other 21 are promises, not
#: rot, so failing on them would just mean deleting the test. Listing them here
#: instead makes the promise auditable: this set IS the manual's remaining work,
#: and an entry disappears the moment its chapter lands (if the file exists, the
#: allowlist is not consulted). A typo'd chapter name still fails, because it
#: will not be in this list.
UNWRITTEN_MANUAL_CHAPTERS = {
    f"docs/manual/{name}"
    for name in (
        "01_introduction.md", "02_installation.md", "05_specs.md", "06_grids.md",
        "07_layers.md", "08_flow_packages.md", "09_boundary_conditions.md",
        "10_advanced_packages.md", "11_multiphysics.md", "12_projects_runs.md",
        "13_results.md", "14_visualization.md", "15_particle_tracking.md",
        "16_parallel.md", "17_pest.md", "18_serialization.md",
        "19_troubleshooting.md", "20_api_reference.md", "21_comparison.md",
        "A_canonical_model.md", "B_glossary.md",
    )
}

FENCE = re.compile(r"^\s*(```|~~~)")
#: `[text](target)` -- target captured up to a space (which starts a title) or `)`.
LINK = re.compile(r"\[[^\]]*\]\(([^)\s]+)")


def _markdown_files() -> list[Path]:
    seen: dict[Path, None] = {}
    for folder in SCANNED:
        base = ROOT / folder
        if not base.exists():
            continue
        for path in base.rglob("*.md") if folder != "." else base.glob("*.md"):
            if any(part in {".git", "node_modules", "_vendor", ".claude"} for part in path.parts):
                continue
            seen.setdefault(path, None)
    return sorted(seen)


IDS = [p.relative_to(ROOT).as_posix() for p in _markdown_files()]


@pytest.mark.parametrize("relative", IDS)
def test_code_fences_are_balanced(relative):
    """An odd number of fences means everything after the last one renders as
    code -- the rest of the document simply disappears as prose."""

    lines = (ROOT / relative).read_text(encoding="utf-8").splitlines()
    open_line = None
    marker = None
    for number, line in enumerate(lines, 1):
        match = FENCE.match(line)
        if not match:
            continue
        if open_line is None:
            open_line, marker = number, match.group(1)
        elif line.strip().startswith(marker):
            open_line, marker = None, None

    assert open_line is None, (
        f"{relative}: code fence opened at line {open_line} is never closed, so "
        "everything below it renders as a code block."
    )


@pytest.mark.parametrize("relative", IDS)
def test_local_links_resolve(relative):
    """A relative link to a file that moved is the standard decay product of a
    refactor, and nothing else in CI would notice."""

    path = ROOT / relative
    broken = []
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        for target in LINK.findall(line):
            if target.startswith(("http://", "https://", "mailto:", "#")):
                continue
            local = target.split("#", 1)[0]
            # `path.py:120` is this repo's own convention for a pointer into the
            # source tree (see the manual's conventions table); the line number
            # is not part of the filename.
            local = re.sub(r":\d+$", "", local)
            if not local:
                continue
            resolved = (path.parent / local).resolve()
            try:
                key = resolved.relative_to(ROOT).as_posix()
            except ValueError:
                continue                      # points outside the repo; not ours
            if key in GENERATED or key in UNWRITTEN_MANUAL_CHAPTERS:
                continue
            if resolved.exists():
                continue
            broken.append(f"line {number}: {target}")

    assert not broken, f"{relative} links to files that do not exist: {broken}"


def test_the_unwritten_chapter_list_does_not_outlive_the_chapters():
    """The allowlist must shrink as the manual is written. If a chapter now
    exists and is still listed, the list has stopped describing reality -- and
    the next person reading it would be told work remains that is done."""

    stale = sorted(
        name for name in UNWRITTEN_MANUAL_CHAPTERS if (ROOT / name).exists()
    )
    assert not stale, (
        f"{stale} now exist(s). Delete them from UNWRITTEN_MANUAL_CHAPTERS -- "
        "that set is the manual's remaining work, not a permanent exemption."
    )
