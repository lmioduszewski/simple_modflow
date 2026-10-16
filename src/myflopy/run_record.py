"""What MF6 actually did, read back from a run's listing file (design step 1).

A run manifest records ``status`` -- ``"created"``, ``"built"``, ``"written"``,
``"completed"`` -- which is **myflopy's own verb**: what the library did to the
run. Nothing recorded what *MF6* did afterwards, so every run executed outside
``run.execute()`` sat at ``status="written", success=None`` forever. Measured on
a real project: four runs, all reading ``written``/``None``; three had reached
``Normal termination`` and one had not.

    ``status`` is myflopy's verb. ``outcome`` is MF6's verdict.
    Two axes; never collapse them.

It also holds the two rules a run's RECORD is checked against its FILES by (ledger
611): which files MF6 reads (:func:`_scan_simulation_inputs`, the walk a PEST
build's staging copy and ``stale`` share) and how two written simulations are compared
(:func:`written_differences`), plus the verdict ``run.inputs`` returns
(:class:`InputsVerdict`).

This module is deliberately stdlib-only and imports nothing from myflopy, so it
sits at the bottom of the import graph and adds no deferred import to
``workspace``.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

__all__ = [
    "InputsVerdict",
    "ListingOutcome",
    "pest_tree_root",
    "read_listing_outcome",
    "written_differences",
]

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
        workspace's own name files lead MF6 to read -- their entries and the
        ``FILEIN`` and ``OPEN/CLOSE`` files those name (:func:`_scan_simulation_inputs`,
        ledger 611) -- so with no ``mfsim.nam`` this reads ``False`` ("cannot
        tell"). It compares times, not contents: whether the files still say what the
        run's record says is ``run.inputs``. Default ``False``.
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


#: Name-file blocks whose entries name an INPUT file in their second field: the
#: simulation name file's ``timing``/``models``/``exchanges``/``solutiongroup`` and a
#: model name file's ``packages``. Never ``options``, whose ``LIST`` names an OUTPUT
#: (read as one by the output scan, ledger 453).
_NAME_FILE_BLOCKS = frozenset({"timing", "models", "exchanges", "solutiongroup", "packages"})

#: The file type opening a name-file entry -- ``TDIS6``, ``gwf6``, ``GWF6-GWT6``,
#: ``ims6``, ``DISV6``, ``OBS6`` ... -- which a solution group's ``MXITER 2`` is not.
_MF6_FILE_TYPE = re.compile(r"(?i)^[a-z]+6(-[a-z]+6)?$")

#: One field of an MF6 input line: a quoted name (spaces allowed) or a bare word.
_MF6_FIELD = re.compile(r"'[^']*'|\"[^\"]*\"|\S+")

#: A file an MF6 input names to READ: ``FILEIN <name>`` (TS6, OBS6, TAS6, TAB6, ATS6,
#: SPC6, HPC6, FMI's flows ...) or an external array or list, ``OPEN/CLOSE <name>``.
#: Never ``FILEOUT``, which names an output.
_MF6_INPUT_REFERENCE = re.compile(
    r"(?i)(?<![\w/])(FILEIN|OPEN/CLOSE)\s+('[^']*'|\"[^\"]*\"|\S+)"
)

#: A file an MF6 input names to WRITE: ``FILEOUT <name>`` (OC's head and budget
#: files, an OBS6 continuous block, a budget CSV, LAK/SFR stage ...).
_MF6_OUTPUT_REFERENCE = re.compile(r"(?i)(?<![\w/])FILEOUT\s+('[^']*'|\"[^\"]*\"|\S+)")


def _mf6_field(token: str) -> str:
    """One MF6 field with its quotes removed."""

    if len(token) >= 2 and token[0] == token[-1] and token[0] in "'\"":
        return token[1:-1]
    return token


def _mf6_text(path: Path) -> str | None:
    """``path`` as text with MF6's comments removed, or ``None`` for a binary file."""

    return _mf6_text_of(path.read_bytes())


def _mf6_text_of(data: bytes) -> str | None:
    """:func:`_mf6_text` of a file's bytes already read."""

    if b"\0" in data[:8192]:
        return None
    # Cut only the lines that carry a comment: a regex split per line cost ~1 s over the
    # 200,000 lines of a 20,323-cell model's arrays, where the scan is read on every
    # `run.outcome` and `run.inputs` (ledger 611).
    return "\n".join(
        line.split("#", 1)[0].split("!", 1)[0] if ("#" in line or "!" in line) else line
        for line in data.decode("utf-8", errors="replace").splitlines()
    )




def _scan_simulation_inputs(
    sim_ws: str | Path,
) -> tuple[
    list[Path], list[tuple[str, str]], list[Path], list[tuple[str, str]], dict[Path, Path]
]:
    """Every file MODFLOW 6 reads in ``sim_ws``, followed from ``mfsim.nam`` as MF6 opens them.

    The ONE rule for "the files MF6 reads" (ledger 611): a PEST build's staging
    copy (``pest.project._simulation_input_files``, ledger 383), a listing's ``stale``
    flag and a run's record-against-files check (``Run.inputs``) all walk it. Each
    name-file entry (TDIS, each model's name file and its packages, exchanges,
    solutions), each ``FILEIN`` a scanned file names (TS6, OBS6, TAS6, TAB6, ATS6,
    SPC6, HPC6, FMI's flow files ...) and each ``OPEN/CLOSE`` external array or list --
    names MF6 resolves against the directory it runs in. A ``FILEIN`` file is scanned
    in turn; an ``OPEN/CLOSE`` file holds data and is not.

    Parameters
    ----------
    sim_ws : str or pathlib.Path
        The simulation directory, holding ``mfsim.nam``.

    Returns
    -------
    tuple
        Five items: the input files (relative to ``sim_ws``, ``mfsim.nam`` first, each
        once); the RELATIVE input names leading outside ``sim_ws``, each with what names
        it; the outputs (``FILEOUT``, a model's ``LIST``) written inside, relative; those
        written outside, with what names them; and each input file's OWNER -- the
        name-file entry that leads to it (a package file for its ``FILEIN``/``OPEN/CLOSE``
        files, every name-file entry for itself).
    """

    root = Path(sim_ws)
    found: dict[Path, Path] = {}
    outside: list[tuple[str, str]] = []
    written: dict[Path, None] = {}
    written_outside: list[tuple[str, str]] = []
    # (relative name, how the file is read, what it is): "simulation" / "model"
    # name files are read for their blocks, "input" files scanned for
    # FILEIN/OPEN/CLOSE, "data" files (OPEN/CLOSE) copied unread.
    # The fourth item is the owner: the name-file entry the file belongs to, ``None``
    # for a name-file entry itself (it owns itself).
    pending: list[tuple[str, str, str, Path | None]] = [
        ("mfsim.nam", "simulation", "mfsim.nam", None)
    ]

    def _record_output(name_written: str, owner: str) -> None:
        """File an output name as written inside ``sim_ws`` or outside it."""

        normal = Path(os.path.normpath(name_written))
        if normal.is_absolute() or normal.parts[:1] == ("..",):
            written_outside.append((name_written, owner))
        else:
            written[normal] = None

    while pending:
        name, kind, owner, anchor = pending.pop(0)
        relative = Path(os.path.normpath(name))
        if relative.is_absolute() or relative in found:
            continue
        if relative.parts[:1] == ("..",):
            outside.append((name, owner))
            continue
        path = root / relative
        if not path.is_file():
            continue
        anchor = relative if anchor is None else anchor
        found[relative] = anchor
        if kind == "data":
            continue
        data = path.read_bytes()
        # A package file naming no file at all -- every array and list file of a large
        # model -- is not decoded: its text is only ever searched for these three
        # keywords, and decoding and cutting comments cost 0.2 s per read over a
        # 20,323-cell model's 40 MB of inputs, paid by `outcome`, `inputs` and
        # `execute()` alike (ledger 612). A keyword inside a comment only
        # sends the file down the full path below, which drops it.
        if kind == "input" and not any(
            keyword in data.lower() for keyword in (b"filein", b"open/close", b"fileout")
        ):
            continue
        text = _mf6_text_of(data)
        if text is None:
            continue
        if kind in ("simulation", "model"):
            block = None
            for line in text.splitlines():
                fields = [_mf6_field(token) for token in _MF6_FIELD.findall(line)]
                if not fields:
                    continue
                word = fields[0].lower()
                if word == "begin":
                    block = fields[1].lower() if len(fields) > 1 else None
                elif word == "end":
                    block = None
                elif kind == "model" and block == "options" and word == "list" and len(fields) > 1:
                    # The model's listing file, the one output named outside a FILEOUT
                    # (ledger 453): every realization writes it like any other.
                    _record_output(fields[1], f"LIST in {relative}")
                elif (
                    block in _NAME_FILE_BLOCKS
                    and len(fields) > 1
                    and _MF6_FILE_TYPE.match(fields[0])
                ):
                    entry = fields[0].upper() + (f" {fields[2]!r}" if len(fields) > 2 else "")
                    pending.append((
                        fields[1],
                        "model" if block == "models" else "input",
                        f"the {entry} entry of {relative}",
                        None,
                    ))
        named_by = "" if kind in ("simulation", "model") else f", named by {owner}"
        # The two regexes cost ~1 s over a large model's array files, which name no
        # file at all: look for the keywords first (ledger 611).
        lowered = text.lower()
        names_input = "filein" in lowered or "open/close" in lowered
        for keyword, token in (_MF6_INPUT_REFERENCE.findall(text) if names_input else ()):
            pending.append((
                _mf6_field(token),
                "input" if keyword.upper() == "FILEIN" else "data",
                f"{keyword.upper()} in {relative}{named_by}",
                anchor,
            ))
        for token in (_MF6_OUTPUT_REFERENCE.findall(text) if "fileout" in lowered else ()):
            _record_output(_mf6_field(token), f"FILEOUT in {relative}{named_by}")
    return list(found), outside, list(written), written_outside, dict(found)




def _is_stale(listing: Path) -> bool:
    """Whether any MF6 INPUT is newer than the listing.

    A listing older than its inputs describes a superseded run -- the verdict is
    real, but no longer about what is on disk now.

    The input set is :func:`_scan_simulation_inputs`' -- the files the workspace's own
    name files lead MF6 to, ``FILEIN`` and ``OPEN/CLOSE`` files included (it read the
    name files' entries alone until ledger 611) -- rather than a guess about
    file extensions. An earlier version deny-listed seven output
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
    for name in _scan_simulation_inputs(listing.parent)[0]:
        candidate = listing.parent / name
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


# -- a run's record against its input files (ledger 611) ------------------


@dataclass(frozen=True, slots=True)
class InputsVerdict:
    """Whether a run's MF6 input files still say what its record says (``run.inputs``).

    A run folder holds two descriptions of one model: the RECORD -- ``run.spec``, its
    frozen snapshot (ledger 251), its own library copies (ledger 230) and digests --
    which a rebuild, ``project.runs.drift()`` and the write-backs read; and the input
    FILES, which are what MF6 runs (ledger 611).

    Attributes
    ----------
    verdict : str
        One of five words. ``"ok"`` -- every input file is what the record writes and
        none changed since myflopy last wrote it. ``"edited"`` -- at least one departs
        (``departs`` says which, and how). ``"unfingerprinted"`` -- the files were
        written before runs recorded what they wrote, or rewritten from such files:
        agreement is UNKNOWN, never assumed. ``"unwritten"`` -- the run has not been
        written since it was prepared or rebuilt, so any files present are an earlier
        write's. ``"no record"`` -- the workspace has no ``run.json``.
    departs : dict of str to str, optional
        Input file (relative to the run folder) to how it departs from the record:
        ``"edited through FloPy after run.build()"``, ``"edited through FloPy after
        run.open()"``, ``"changed outside myflopy after its last write"``, ``"added"``
        (an input the last write did not write) or ``"missing"`` (one it wrote that is
        gone). Empty unless ``verdict`` is ``"edited"``. Default empty.
    written_at : str or None, optional
        ISO-8601 UTC time myflopy last wrote the inputs; ``None`` (the default) when it
        recorded none.
    """

    verdict: str
    departs: dict[str, str] = field(default_factory=dict)
    written_at: str | None = None

    @property
    def matches(self) -> bool | None:
        """``True`` for ``"ok"``, ``False`` for ``"edited"``, ``None`` when unknown.

        ``None`` is not "yes": an unfingerprinted run may or may not agree, and an
        unwritten one has nothing on disk that is its own.
        """

        if self.verdict == "ok":
            return True
        if self.verdict == "edited":
            return False
        return None

    def __repr__(self) -> str:
        """Compact: the verdict, plus the departing files when there are any."""

        pieces = [f"verdict={self.verdict!r}"]
        if self.departs:
            pieces.append(f"departs={self.departs!r}")
        return f"InputsVerdict({', '.join(pieces)})"


#: A number as MF6 reads one: optional sign, digits, a decimal point, and an exponent
#: in E or (Fortran's) D.
_MF6_NUMBER = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eEdD][+-]?\d+)?$")


def _number_and_unit(token: str) -> tuple[float, float] | None:
    """``token`` as ``(value, unit of its last written digit)``, or ``None`` for a word."""

    match = _MF6_NUMBER.match(token)
    if match is None:
        return None
    mantissa, exponent = match.group(1), match.group(2)
    decimals = len(mantissa.split(".", 1)[1]) if "." in mantissa else 0
    power = int(exponent[1:]) if exponent else 0
    value = float(token.replace("d", "e").replace("D", "e"))
    return value, 10.0 ** (power - decimals)


def _same_token(files: str, record: str) -> bool:
    """Whether two MF6 tokens say the same thing (ledger 611's token rule).

    Words compare without regard to case (FloPy writes a loaded keyword in the case it
    read, and lower-cases a boundname); numbers by VALUE, equal when they differ by no
    more than half a unit of the coarser token's last written digit -- ``1`` and
    ``1.00000000`` are one number, ``1`` and ``31.00000000`` are not.
    """

    if files == record:
        return True
    left, right = _number_and_unit(files), _number_and_unit(record)
    if left is None or right is None:
        return files.casefold() == record.casefold()
    tolerance = 0.5 * max(left[1], right[1])
    return abs(left[0] - right[0]) <= tolerance * (1 + 1e-9)


def _blocks(text: str) -> list[tuple[str, list[str]]]:
    """``text`` as ``(block, tokens)`` pairs, ``BEGIN <block>`` naming each (``""`` outside)."""

    blocks: list[tuple[str, list[str]]] = [("", [])]
    for line in text.splitlines():
        fields = [_mf6_field(token) for token in _MF6_FIELD.findall(line)]
        if not fields:
            continue
        word = fields[0].casefold()
        if word == "begin":
            blocks.append((" ".join(fields[1:]).upper(), []))
        elif word == "end":
            blocks.append(("", []))
        else:
            blocks[-1][1].extend(fields)
    return [(name, tokens) for name, tokens in blocks if tokens or name]


def _differences_in(name: str, files: str, record: str) -> list[tuple[str, str, str]]:
    """``(file, block, what)`` for each block of one text file where the two differ."""

    left, right = _blocks(files), _blocks(record)
    if [block for block, _ in left] != [block for block, _ in right]:
        count = sum(len(tokens) for _, tokens in left)
        other = sum(len(tokens) for _, tokens in right)
        names, others = [block for block, _ in left], [block for block, _ in right]
        longer = names if len(names) > len(others) else others
        first = next(
            (a for a, b in zip(names, others, strict=False) if a != b),
            longer[min(len(names), len(others))],
        )
        return [(name, first, f"written differently: {count} tokens (files) against "
                              f"{other} (record), in blocks that do not pair up")]
    rows: list[tuple[str, str, str]] = []
    for (block, a), (_, b) in zip(left, right, strict=True):
        if len(a) != len(b):
            rows.append((name, block, f"written differently: {len(a)} tokens (files) "
                                      f"against {len(b)} (record)"))
            continue
        differing = [(x, y) for x, y in zip(a, b, strict=True) if not _same_token(x, y)]
        if not differing:
            continue
        shown = "; ".join(f"{x} (files) vs {y} (record)" for x, y in differing[:3])
        what = f"{len(differing)} of {len(a)} values differ: {shown}"
        if len(differing) > 3:
            gaps = [
                abs(left_[0] - right_[0])
                for x, y in differing
                if (left_ := _number_and_unit(x)) and (right_ := _number_and_unit(y))
            ]
            if gaps:
                what += f"; largest absolute difference {max(gaps):.6g}"
        rows.append((name, block, what))
    return rows


def written_differences(files: str | Path, record: str | Path) -> list[tuple[str, str, str]]:
    """Every difference between two written simulations, file by file (ledger 611).

    The ONE record-against-files comparison -- ``Run._record_differences``, the full
    check behind the write-backs, calls it, and a library import of a run must too --
    with ``files`` a run's own folder and
    ``record`` its record built and written in a scratch directory. Not in memory --
    FloPy writes 8 decimals, so a value read back from a file is not the value built --
    and not text for text: comments are dropped (FloPy stamps the time it wrote), words
    compare without regard to case and numbers by value to the coarser precision the
    two were written at (:func:`_same_token`). A binary file compares by bytes.

    Parameters
    ----------
    files : str or pathlib.Path
        The simulation directory holding the files MF6 ran (``mfsim.nam`` and what it
        leads to, through :func:`_scan_simulation_inputs`).
    record : str or pathlib.Path
        A directory holding the record written by FloPy, walked the same way.

    Returns
    -------
    list of tuple of (str, str, str)
        ``(file, block, what)``, the file relative to both directories, the block as
        ``BEGIN`` names it (``"PERIODDATA"``, ``"PERIOD 1"``; ``""`` for a whole file),
        and what differs -- the first three differing values, ``"(files)"`` beside
        ``"(record)"``, or ``"only in the files"`` / ``"only in the record"`` for a file
        one side does not read. Empty when the two agree.
    """

    files, record = Path(files), Path(record)
    ours = _scan_simulation_inputs(files)[0]
    theirs = _scan_simulation_inputs(record)[0]
    rows: list[tuple[str, str, str]] = []
    for name in list(dict.fromkeys([*ours, *theirs])):
        label = name.as_posix()
        if name not in theirs:
            rows.append((label, "", "only in the files"))
            continue
        if name not in ours:
            rows.append((label, "", "only in the record"))
            continue
        left, right = _mf6_text(files / name), _mf6_text(record / name)
        if left is None or right is None:
            if (files / name).read_bytes() != (record / name).read_bytes():
                rows.append((label, "", "binary file differs"))
            continue
        if left.split() == right.split():
            continue
        rows.extend(_differences_in(label, left, right))
    return rows
