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

__all__ = ["ListingOutcome", "read_listing_outcome", "pest_tree_root"]

# PEST calibration markers, DUPLICATED from myflopy.modflow.mf6.pest.runs rather
# than imported. That module is import layer 0 and could be imported from here
# without breaking layering -- but importing it executes the pest package
# __init__, measured at **1.43 s** against an `import myflopy` that costs 0.02 s.
# A 70x import regression to reach three string constants is not a trade worth
# making, so they are pinned equal by a test instead (the repo's generate-then-pin
# pattern, as for api_snapshot.json and import_layers.json).
PEST_ROOT_SUFFIX = ".pest"
LEGACY_PEST_DIRNAME = "pest"
PEST_METADATA = "myflopy_pest_metadata.json"

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
    outcome : str
        One of five words. ``"normal"`` -- MF6 reported normal termination.
        ``"failed"`` -- MF6 reported premature termination.
        ``"incomplete"`` -- a listing exists but reports neither, so the run was
        killed, crashed, or is still going. This is what a real interrupted run
        looks like: the file simply stops mid-iteration.
        ``"absent"`` -- no listing file, so MF6 never started here.
        ``"superseded"`` -- a listing exists but PREDATES this run's
        preparation, so its verdict belongs to a previous run in the same
        directory. ``prepare_run(overwrite=True)`` leaves the old listing and
        results in place, so without this a re-prepared, never-executed run
        reported the previous execution's outcome and elapsed time.
    finished_at : str or None, optional
        MF6's own run-end stamp (``yyyy/mm/dd hh:mm:ss``), verbatim, or ``None``
        (the default) when the listing's tail carries no ``Run end date and time``
        line -- always so for ``"absent"``, ``"superseded"`` and an unreadable file.
    elapsed : str or None, optional
        MF6's own elapsed-time line, verbatim (e.g. ``"3 Hours, 53 Minutes,  1
        Seconds"``), or ``None`` (the default) when the tail carries none.
    stale : bool, optional
        Whether any MF6 input file is NEWER than the listing. A listing that
        predates its inputs describes a superseded run -- its verdict is real but
        no longer about what is on disk now. The inputs are the files the
        workspace's own name files declare, so with no ``mfsim.nam`` this reads
        ``False`` ("cannot tell"). Default ``False``.
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
        return None  # incomplete, absent, superseded -- no verdict for THIS run

    def __repr__(self) -> str:
        """Compact: the verdict, plus staleness when it matters."""

        pieces = [f"outcome={self.outcome!r}"]
        if self.elapsed:
            pieces.append(f"elapsed={self.elapsed!r}")
        if self.stale:
            pieces.append("stale=True")
        return f"ListingOutcome({', '.join(pieces)})"


def pest_tree_root(workspace: Path, root: Path) -> Path | None:
    """The calibration directory ``workspace`` belongs to, or ``None``.

    A calibration builds into ``<model workspace>.pest/`` -- a SIBLING of the run
    directory, so inside the project root -- holding the PstFrom template,
    ``<name>_ies_master``, ``<name>_prior_master`` and every pyEMU worker. Each holds
    a copy of the model's MF6 inputs (the whole model workspace until ledger
    383), so each carries an ``mfsim.nam`` and
    would otherwise be listed as a run: one calibration with 8 workers becomes
    ~11 rows in a ledger whose job is telling you what you ran.

    Two markers, either sufficient. An ancestor named ``*.pest`` (or the legacy
    ``pest``) catches the whole tree at once, since pyEMU is launched with
    ``worker_root`` inside it. An ancestor holding
    ``myflopy_pest_metadata.json`` covers a calibration built at a custom
    ``workspace=``.
    """

    outermost: Path | None = None
    current = workspace
    while True:
        name = current.name
        if (
            name.endswith(PEST_ROOT_SUFFIX)
            or name == LEGACY_PEST_DIRNAME
            or (current / PEST_METADATA).is_file()
        ):
            # Keep walking: the OUTERMOST match is the calibration. The template
            # directory carries the metadata marker and sits inside the `.pest`
            # root, so returning the first match split one calibration into two
            # -- the tree, and its own template counted separately.
            outermost = current
        if current == root or current.parent == current:
            return outermost
        current = current.parent


def _input_files(workspace: Path) -> set[Path]:
    """The MF6 input files of a workspace, as its own name files declare them.

    ``mfsim.nam`` names the TDIS file, each model's name file and the solution
    files; a model name file then names its package files. Every whitespace
    token in those that resolves to a real file IS an input -- MF6 says so
    itself. Output paths (``.hds``, ``.cbc``) are named INSIDE package files, not
    in the name files, so they never appear here.

    Returns an empty set when there is no ``mfsim.nam``, which callers read as
    "cannot tell" rather than "nothing changed".
    """

    root = workspace / "mfsim.nam"
    if not root.is_file():
        return set()
    inputs: set[Path] = set()
    pending = [root]
    while pending:
        current = pending.pop()
        if current in inputs:
            continue
        inputs.add(current)
        try:
            text = current.read_text(encoding="utf-8", errors="replace")
        except OSError:  # unreadable -- what we already collected still stands
            continue
        for token in text.split():
            candidate = workspace / token
            if candidate.is_file() and candidate not in inputs:
                # Only a name file can name further inputs; a package file's
                # tokens are data, and following them would sweep in outputs.
                if candidate.suffix.lower() == ".nam":
                    pending.append(candidate)
                else:
                    inputs.add(candidate)
    return inputs


def _is_stale(listing: Path) -> bool:
    """Whether any MF6 INPUT is newer than the listing.

    A listing older than its inputs describes a superseded run -- the verdict is
    real, but no longer about what is on disk now.

    The input set comes from the workspace's own name files rather than from a
    guess about file extensions. An earlier version deny-listed seven output
    suffixes and treated everything else as an input, so saving a figure or a
    note into a run directory marked the run stale -- and it missed MF6's own
    ``.cbb`` budget file, which is an output it did not list. Measured: adding
    ``heads.png``, ``notes.txt`` or ``flow.cbb`` to a clean finished run all
    flipped it to stale.
    """

    try:
        listing_time = listing.stat().st_mtime
    except OSError:
        return False
    for candidate in _input_files(listing.parent):
        try:
            if candidate.stat().st_mtime > listing_time:
                return True
        except OSError:  # vanished between listing and stat -- not evidence
            continue
    return False


def read_listing_outcome(
    workspace: Path | str, *, since: float | None = None
) -> ListingOutcome:
    """Read what MF6 reported in ``workspace/mfsim.lst``.

    Parameters
    ----------
    workspace : str or pathlib.Path
        The MF6 run directory; the listing read is ``workspace / "mfsim.lst"``
        (the name is fixed, not taken from ``mfsim.nam``). Only the last 8192
        bytes of the listing are read.
    since : float or None, optional
        A POSIX timestamp (seconds since the epoch, as ``os.stat().st_mtime``) for
        when the run in this directory was prepared. A listing whose modification
        time is older than that describes a PREVIOUS run --
        ``prepare_run(overwrite=True)`` reuses the directory and leaves the old
        listing and results behind -- so its verdict is reported as
        ``"superseded"`` rather than attributed to a run that has not executed.
        ``None`` (the default) skips that comparison.

    Returns
    -------
    ListingOutcome

    Notes
    -----
    It does not raise on the listing: a missing file is ``"absent"``, and an
    unreadable or truncated one is ``"incomplete"``, because a run you cannot
    classify must still be listable beside the ones you can.
    """

    listing = Path(workspace) / _LISTING_NAME
    try:
        stat = listing.stat()
        if since is not None and stat.st_mtime < since:
            return ListingOutcome(outcome="superseded", stale=_is_stale(listing))
        size = stat.st_size
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


#: MF6's line for a time step its solution did not converge on, one-based. Written to
#: ``mfsim.lst`` whether the run then stops (the default) or carries on (IMS ``CONTINUE``),
#: so it names every unconverged step either way.
_NOT_CONVERGED = re.compile(
    r"Solution Group \d+ did not converge for stress period (\d+) and time step (\d+)"
)


def nonconverged_steps(workspace: Path | str) -> set[tuple[int, int]]:
    """Every time step MF6 reports its solution did NOT converge on, from ``mfsim.lst``.

    A run that fails stops after writing the output of the step that failed -- at the time
    that step would have ended -- and a run told to ``CONTINUE`` past a failure carries on
    with such steps inside it. Either way the numbers saved for those steps are not a
    solution, and MF6 says which steps they were, one line each.

    Parameters
    ----------
    workspace : str or pathlib.Path
        The MF6 run directory; the listing read is ``workspace / "mfsim.lst"``.

    Returns
    -------
    set of tuple of (int, int)
        ``(kstp, kper)``, zero-based. Empty when every step converged, and when there is no
        listing to read (a run with no listing says nothing either way).
    """

    listing = Path(workspace) / _LISTING_NAME
    found: set[tuple[int, int]] = set()
    try:
        with listing.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if "did not converge" not in line:
                    continue
                match = _NOT_CONVERGED.search(line)
                if match:
                    found.add((int(match.group(2)) - 1, int(match.group(1)) - 1))
    except OSError:
        return set()
    return found


def utc_now() -> str:
    """Return an ISO-8601 UTC timestamp, matching the run manifest's format."""

    return datetime.now(timezone.utc).isoformat(timespec="seconds")
