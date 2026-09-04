"""What MF6 actually did, read back from a run's listing file (design step 1).

A run manifest records ``status`` -- ``"created"``, ``"built"``, ``"written"``,
``"completed"`` -- which is **myflopy's own verb**: what the library did to the
run. Nothing recorded what *MF6* did afterwards, so every run executed outside
``run.execute()`` sat at ``status="written", success=None`` forever. Measured on
a real project: four runs, all reading ``written``/``None``; three had reached
``Normal termination`` and one had not.

    ``status`` is myflopy's verb. ``outcome`` is MF6's verdict.
    Two axes; never collapse them.

This module is deliberately stdlib-only and imports nothing from myflopy, so it
sits at the bottom of the import graph and adds no deferred import to
``workspace``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

__all__ = ["ListingOutcome", "read_listing_outcome"]

_LISTING_NAME = "mfsim.lst"

_NORMAL = "Normal termination of simulation."
_PREMATURE = "Premature termination of simulation."

# " Run end date and time (yyyy/mm/dd hh:mm:ss): 2026/09/03 12:50:17"
_END_TIME = re.compile(
    r"Run end date and time \(yyyy/mm/dd hh:mm:ss\):\s*(?P<stamp>[\d/]+\s+[\d:]+)"
)
# " Elapsed run time:  3 Hours, 53 Minutes,  1 Seconds"
_ELAPSED = re.compile(r"Elapsed run time:\s*(?P<elapsed>.+?)\s*$", re.MULTILINE)

# Only the tail is read: a listing can be very large (a long transient run writes
# one block per time step), and everything decided here is in the last few lines.
_TAIL_BYTES = 8192


@dataclass(frozen=True, slots=True)
class ListingOutcome:
    """What ``mfsim.lst`` says MF6 did, and whether it still describes the inputs.

    Attributes
    ----------
    outcome
        ``"normal"`` -- MF6 reported normal termination.
        ``"failed"`` -- MF6 reported premature termination.
        ``"incomplete"`` -- a listing exists but reports neither, so the run was
        killed, crashed, or is still going. This is what a real interrupted run
        looks like: the file simply stops mid-iteration.
        ``"absent"`` -- no listing file, so MF6 never started here.
    finished_at
        MF6's own run-end stamp, verbatim, or ``None``.
    elapsed
        MF6's own elapsed-time line, verbatim (e.g. ``"3 Hours, 53 Minutes,  1
        Seconds"``), or ``None``.
    stale
        Whether any MF6 input file is NEWER than the listing. A listing that
        predates its inputs describes a superseded run -- its verdict is real but
        no longer about what is on disk now.
    """

    outcome: str
    finished_at: str | None = None
    elapsed: str | None = None
    stale: bool = False

    @property
    def succeeded(self) -> bool | None:
        """``True``/``False`` for a run MF6 finished, ``None`` when it did not.

        ``None`` is not "no": an interrupted or never-started run has no verdict,
        and reporting it as failure would be the same class of lie as reporting
        ``status="written"`` as success.
        """

        if self.outcome == "normal":
            return True
        if self.outcome == "failed":
            return False
        return None

    def __repr__(self) -> str:
        """Compact: the verdict, plus staleness when it matters."""

        pieces = [f"outcome={self.outcome!r}"]
        if self.elapsed:
            pieces.append(f"elapsed={self.elapsed!r}")
        if self.stale:
            pieces.append("stale=True")
        return f"ListingOutcome({', '.join(pieces)})"


def _is_stale(listing: Path) -> bool:
    """Whether any MF6 input beside ``listing`` is newer than it.

    Compares against the simulation name file and the package inputs rather than
    every file in the directory: outputs (``.hds``, ``.cbc``) are written DURING
    the run and are legitimately newer, so including them would mark every run
    stale.
    """

    try:
        listing_time = listing.stat().st_mtime
    except OSError:
        return False
    workspace = listing.parent
    outputs = {".hds", ".cbc", ".ucn", ".bud", ".lst", ".csv", ".json"}
    for candidate in workspace.iterdir():
        if not candidate.is_file() or candidate.suffix.lower() in outputs:
            continue
        try:
            if candidate.stat().st_mtime > listing_time:
                return True
        except OSError:  # vanished between listing and stat -- not evidence
            continue
    return False


def read_listing_outcome(workspace: Path | str) -> ListingOutcome:
    """Read what MF6 reported in ``workspace/mfsim.lst``.

    Never raises: an unreadable or truncated listing is reported as
    ``"incomplete"``, because a run you cannot classify must still be listable
    beside the ones you can.
    """

    listing = Path(workspace) / _LISTING_NAME
    try:
        size = listing.stat().st_size
        with listing.open("rb") as handle:
            if size > _TAIL_BYTES:
                handle.seek(size - _TAIL_BYTES)
            tail = handle.read().decode("utf-8", errors="replace")
    except FileNotFoundError:
        return ListingOutcome(outcome="absent")
    except OSError:
        return ListingOutcome(outcome="incomplete")

    if _NORMAL in tail:
        outcome = "normal"
    elif _PREMATURE in tail:
        outcome = "failed"
    else:
        outcome = "incomplete"

    end = _END_TIME.search(tail)
    elapsed = _ELAPSED.search(tail)
    return ListingOutcome(
        outcome=outcome,
        finished_at=end.group("stamp") if end else None,
        elapsed=elapsed.group("elapsed") if elapsed else None,
        stale=_is_stale(listing),
    )


def utc_now() -> str:
    """Return an ISO-8601 UTC timestamp, matching the run manifest's format."""

    return datetime.now(timezone.utc).isoformat(timespec="seconds")
