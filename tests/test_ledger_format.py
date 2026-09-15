"""The compromise ledger is a document the code cites, so its shape is testable.

265 entries over ~9,700 lines, cited 1,527 times from `src/`, `tests/`, the
notebooks and CLAUDE.md. Nothing checked that a cited number existed, that two
entries had not claimed the same one, or that the generated index still matched
the entries -- so the failure mode was a reader following a citation to nothing,
which is worse than no citation at all.

The entry shape these pin is the one the entries actually have, measured, not the
one the file's own preamble used to describe.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "docs" / "compromises_and_deferrals.md"

sys.path.insert(0, str(ROOT / "scripts"))
from derive_ledger_index import entries as _entries  # noqa: E402

# A ledger citation in prose or a comment, in every form the repo actually writes:
# `ledger 250`, the plural `ledgers 249-257`, and a comma list. The first version
# matched only the singular with a `(?!\d|-\d)` tail borrowed from the date guard,
# which also swallowed every RANGE -- 17 `ledgers N` and 8 `ledger N-M` citations
# went unchecked, including entry 264's own (ledger 265).
#
# Known and deliberate blind spot: a FOUR-digit number is invisible, so `ledger 1300`
# is not reported. Catching it means matching `\d+`, which then matches the year in
# `ledger 2026-09-14` and in every dated title -- a false alarm on real text is worse
# than missing a typo that names an entry 1,000 past the highest one.
CITATION = re.compile(
    r"\bledgers?\s+\d{1,3}(?:\s*(?:[-\u2013,]|and)\s*\d{1,3})*(?!\d)", re.IGNORECASE
)
NUMBER = re.compile(r"\d{1,3}")

SEARCHED = {".py", ".md", ".ipynb"}

# Dates became the convention partway through; every entry from here on has one.
DATED_FROM = 100

# A title may wrap, but not far: the longest in the file spans three lines.
TITLE_LINES = 4


@pytest.fixture(scope="module")
def ledger():
    """Every ledger entry as ``(number, title, date, resolved, line)``."""

    return _entries(LEDGER.read_text().splitlines())


def test_every_entry_opens_with_a_bold_title(ledger):
    """A numbered line at column 0 is an entry, and every entry titles itself."""

    numbered = [
        line for line in LEDGER.read_text().splitlines() if re.match(r"^\d+\. ", line)
    ]
    assert len(numbered) == len(ledger), (
        "a numbered line at column 0 that is not a bold or struck-through entry -- "
        "an entry opens `N. **Title (date).**`, and a continuation line is indented"
    )


def test_every_entry_closes_its_bold_title():
    """An unclosed title runs its bold into the body, and the parser cannot tell.

    The title is read forward to its closing ``**`` because several titles wrap. Drop
    that closing marker and the read simply continues into the entry body until it
    meets some other ``**``, so the index renders a plausible title and all the other
    checks here stay green -- while the document itself has an unterminated bold span
    (ledger 265). Bound it: a title closes within its own opening lines.
    """

    lines = LEDGER.read_text().splitlines()
    unclosed = []
    for index, line in enumerate(lines):
        match = re.match(r"^\d+\. (~~)?\*\*", line)
        if match is None:
            continue
        window = " ".join(lines[index:index + TITLE_LINES])[match.end():]
        if "**" not in window:
            unclosed.append((int(line.split(".", 1)[0]), index + 1))
    assert not unclosed, f"entries whose bold title does not close: {unclosed}"


def test_entry_numbers_run_from_one_without_a_gap(ledger):
    """A gap means a citation points at nothing; the numbers are the addresses."""

    numbers = sorted({number for number, *_ in ledger})
    assert numbers == list(range(1, numbers[-1] + 1)), (
        f"missing entries: {sorted(set(range(1, numbers[-1] + 1)) - set(numbers))}"
    )


def test_a_repeated_number_is_a_resolution_restatement(ledger):
    """A number may appear twice only when the second is the struck resolution.

    Entry 42 does this: the original stands, and a `~~struck~~` restatement below
    records that it was resolved. Any other repeat is two entries answering to one
    citation.
    """

    rows: dict[int, list[tuple[bool, int]]] = {}
    for number, _title, _date, resolved, line in ledger:
        rows.setdefault(number, []).append((resolved, line))
    for number, claims in rows.items():
        if len(claims) == 1:
            continue
        assert len(claims) == 2, f"entry {number} is claimed {len(claims)} times"
        (first_struck, _), (second_struck, second_line) = claims
        assert not first_struck, f"entry {number}'s original claim is struck through"
        assert second_struck, (
            f"entry {number} is claimed twice and the repeat at line {second_line} is "
            "not a struck-through resolution"
        )


def test_entries_since_100_carry_a_date(ledger):
    """The date is what makes an entry readable next to the commit that made it."""

    undated = [n for n, _t, date, _r, _l in ledger if date is None and n >= DATED_FROM]
    assert not undated, f"entries with no date: {undated}"


def test_the_generated_index_is_current():
    """The index is derived, so an entry added without it is a stale table."""

    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "derive_ledger_index.py"), "--check"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def _tracked() -> list[Path]:
    """The files git tracks, which is the only reproducible definition of "the repo".

    An `rglob` walk also reads `.claude/worktrees/`, which holds complete stale copies
    of this repository -- measured, 1,067 of 1,732 citation occurrences came from
    there. That makes the count unreproducible on another machine and in CI, and lets
    a worktree drafting an entry past the main tree's highest number fail this module
    for the main tree (ledger 265).
    """

    listed = subprocess.run(
        ["git", "ls-files"], cwd=ROOT, capture_output=True, text=True, check=True
    )
    return [ROOT / name for name in listed.stdout.split("\n") if name]


def _cited() -> dict[int, list[str]]:
    """Every ledger citation in the tracked repo, mapped to where it was written."""

    found: dict[int, list[str]] = {}
    for path in _tracked():
        if path.suffix not in SEARCHED or not path.is_file():
            continue
        try:
            text = path.read_text(errors="ignore")
        except OSError:  # pragma: no cover - unreadable file in a working tree
            continue
        for match in CITATION.finditer(text):
            for number in NUMBER.findall(match.group(0)):
                found.setdefault(int(number), []).append(str(path.relative_to(ROOT)))
    return found


def test_every_ledger_citation_resolves(ledger):
    """A citation is a promise that the entry is there to be read."""

    numbers = {number for number, *_ in ledger}
    dangling = {n: sorted(set(w)) for n, w in _cited().items() if n not in numbers}
    assert not dangling, f"citations naming no entry: {dangling}"


def test_the_sweep_reaches_every_area_that_cites_the_ledger(ledger):
    """A sweep that stopped reaching the CODE would leave the check above vacuous.

    Counting distinct entries is not enough: `docs/` alone cites 130 of them and
    contains the newest, so a sweep narrowed to documentation passed both of the
    assertions this test used to make while silently no longer checking `src/` or
    `tests/` (ledger 265). Each area is asserted for itself.
    """

    where = {path for paths in _cited().values() for path in paths}
    for area in ("src/", "tests/", "docs/", "CLAUDE.md"):
        assert any(path.startswith(area) for path in where), (
            f"no ledger citation found under {area}; the sweep no longer reaches it"
        )
    assert len(_cited()) > 100, "far fewer entries cited than expected; the sweep is broken"
