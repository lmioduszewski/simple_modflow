"""Write the compromise ledger's index of entries, from the entries themselves.

The ledger is the largest document in the repo -- 265 entries over ~9,700 lines,
cited 1,527 times from `src/`, `tests/`, the notebooks and CLAUDE.md -- and it has
no table of contents, so finding an entry means grepping for a number you already
know. This writes one, between the fences in the file, so the index cannot drift
from the entries the way a hand-kept list would.

    python scripts/derive_ledger_index.py           # rewrite the index in place
    python scripts/derive_ledger_index.py --check   # report drift, change nothing

`tests/test_ledger_format.py` runs `--check`, so an entry added without
regenerating fails there rather than silently leaving the index a line short.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

LEDGER = Path(__file__).resolve().parents[1] / "docs" / "compromises_and_deferrals.md"

BEGIN = "<!-- generated index: python scripts/derive_ledger_index.py -->"
END = "<!-- end generated index -->"

# An entry opens a markdown list item at column 0 with its number, and titles its
# own subject in bold. A RESOLVED entry restates its number with the title struck
# through, which is why the strike marker is optional here rather than a separate
# shape.
ENTRY = re.compile(r"^(\d+)\. (~~)?\*\*")
# Most titles date themselves in parentheses; a few do it after a comma instead,
# and reading only the parenthesised form reported those as undated.
DATE = re.compile(r"[(,]\s*(\d{4}-\d{2}-\d{2})[^)]*\)?")


def entries(lines: list[str]) -> list[tuple[int, str, str | None, bool, int]]:
    """Return ``(number, title, date, resolved, line)`` for every ledger entry.

    A title may wrap across lines -- several do -- so it is read to its closing
    ``**`` rather than to the end of the first line, which would truncate them.
    """

    found: list[tuple[int, str, str | None, bool, int]] = []
    for index, line in enumerate(lines):
        match = ENTRY.match(line)
        if match is None:
            continue
        number = int(match.group(1))
        resolved = match.group(2) is not None
        title = line[match.end():]
        cursor = index
        while "**" not in title and cursor + 1 < len(lines):
            cursor += 1
            title = f"{title.rstrip()} {lines[cursor].strip()}"
        title = _title_before_close(title)
        date = DATE.search(title)
        if date is not None:
            title = title[: date.start()].strip()
        title = title.rstrip(".,").replace("~~", "").strip()
        found.append((number, title, date.group(1) if date else None, resolved, index + 1))
    return found



def _title_before_close(text: str) -> str:
    """Return the bold title in ``text``, ignoring ``**`` inside a code span.

    Two titles document `**kwargs`, and splitting on the first ``**`` truncated
    both mid-word -- the kind of silent half-answer the index exists to avoid.
    """

    masked = re.sub(r"`[^`]*`", lambda m: "\u0000" * len(m.group(0)), text)
    close = masked.find("**")
    return (text if close < 0 else text[:close]).strip()


def _shorten(title: str, limit: int = 110) -> str:
    """Trim ``title`` for the index without cutting an inline markup span in half.

    A blind slice left three rows with an unterminated backtick or an unbalanced
    ``[``, because a handful of titles open with a bracketed CLOSED-by note or carry
    a code span past the limit. Cutting back to the last safe boundary costs a few
    characters and renders.
    """

    if len(title) <= limit:
        return title
    cut = title[: limit - 3]
    for opener, closer in (("`", "`"), ("[", "]")):
        if cut.count(opener) % 2 if opener == closer else cut.count(opener) > cut.count(closer):
            cut = cut[: cut.rfind(opener)]
    return f"{cut.rstrip()}..."


def render(found: list[tuple[int, str, str | None, bool, int]]) -> list[str]:
    """Render the index block, newest first -- the end of the file is the live end."""

    out = [
        BEGIN,
        "",
        f"<!-- {len({number for number, *_ in found})} entries "
        f"({len(found)} rows; a resolved entry is listed twice). "
        "Search this file for `N. **` to open entry N. -->",
        "",
        "| # | Entry | Date |",
        "| --: | --- | --- |",
    ]
    # Sorted by (number, line), never by the whole tuple: two rows with the same
    # number AND the same title -- the documented resolution shape -- put a `None`
    # date against a string and raise `TypeError` out of `sorted`.
    for number, title, date, resolved, _line in sorted(found, key=lambda e: (e[0], e[4]), reverse=True):
        # A pipe would end the table cell, and a handful of titles carry one in a
        # code span; a few others open with a bracketed CLOSED-by note and run past
        # 200 characters, which an index meant for finding an entry does not need.
        shown = _shorten(title).replace("|", "\\|")
        if resolved:
            shown = f"~~{shown}~~"
        out.append(f"| {number} | {shown} | {date or ''} |")
    out.extend(["", END])
    return out


def rewrite(text: str) -> str:
    """Return ``text`` with the index block replaced, inserting it if absent."""

    lines = text.splitlines()
    block = render(entries(lines))
    try:
        start = lines.index(BEGIN)
        stop = lines.index(END)
    except ValueError:
        # No fences yet: the index goes after the standing-rule preamble, above
        # the first heading, so the entries themselves stay where every citation
        # and every reader's scroll position expects them.
        first = next(i for i, line in enumerate(lines) if line.startswith("## "))
        return "\n".join([*lines[:first], *block, "", *lines[first:]]) + "\n"
    return "\n".join([*lines[:start], *block, *lines[stop + 1:]]) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="report drift without writing")
    args = parser.parse_args()

    current = LEDGER.read_text()
    updated = rewrite(current)
    if current == updated:
        print(f"ledger index is current ({len({n for n, *_ in entries(current.splitlines())})} entries)")
        return 0
    if args.check:
        print("the ledger index is stale; run `python scripts/derive_ledger_index.py`")
        return 1
    LEDGER.write_text(updated)
    print(f"wrote the ledger index ({len({n for n, *_ in entries(updated.splitlines())})} entries)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
