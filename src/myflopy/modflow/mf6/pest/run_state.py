"""The state of a PESTPP-IES run directory: whether a run is going in it, and resuming one.

Two halves, both used by :meth:`PestProject.run_ies`, :meth:`PestProject.prior` and
:meth:`PestProject.resume_ies` (ledger 532 and, for the resume, ledger 533):

* **Is a run going there?** A launch writes :data:`RUN_LOCK` into every directory it
  runs in -- the template for a serial run; the master and its worker root
  (``<master>.workers``) for a parallel one -- naming the host, the process and when
  that process started, and removes it when the run returns. :func:`run_activity`
  reads it back: a lock whose process is still running is a LIVE run, one whose
  process is gone is STALE (the kernel was killed), and a lock written on another
  host cannot be checked and counts as live. On Linux it also looks for PEST++
  processes whose working directory is inside the directory, which catches a run
  whose Python process died while PESTPP-IES kept going.

* **Resuming an interrupted IES run.** PESTPP-IES 5.2.16 refuses the PEST ``/r``
  restart switch ("/r option not supported by pestpp-ies", measured); it restarts
  from ensembles instead -- ``ies_par_en`` (the prior, for regularization),
  ``ies_restart_par_en`` / ``ies_restart_obs_en`` (the last completed iteration's
  parameters and their simulated values, which it does not rerun) and ``ies_obs_en``
  (the noise ensemble). It then numbers its iterations from 0 again, so
  :func:`archive_for_resume` moves the interrupted run's outputs aside first and
  :func:`splice_resumed_run` puts them back with the resumed iterations renumbered
  after them, the phi files joined -- one run directory, iterations ``0..N``, as an
  uninterrupted run leaves it.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import socket
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from myflopy._logging import get_logger
from myflopy.modflow.mf6.pest.runs import pestpp_ies_outputs

logger = get_logger(__name__)

#: The file a launch writes into each directory a PESTPP-IES run is going in, and
#: removes when the run returns (ledger 532).
RUN_LOCK = "myflopy_pest_run.lock"

#: The ``written_by`` value of a :data:`RUN_LOCK`.
_LOCK_WRITER = "myflopy PESTPP-IES run"

_LOCK_NOTE = (
    "Written by a myflopy PESTPP-IES launch while the run is going, and removed when "
    "it returns. While the process it names is running, myflopy refuses to delete or "
    "reuse this directory. If that run is over and this file was left behind, delete it."
)

#: Appended to a parallel run's master directory's NAME for its worker root: the
#: workers are ``<master>.workers/worker_0`` .. ``worker_<n-1>`` (ledger 532).
WORKERS_SUFFIX = ".workers"

#: The record :func:`archive_for_resume` writes into the archive directory it makes.
RESUME_RECORD = "myflopy_resume.json"

#: The two forms PESTPP-IES writes an ensemble in (ledger 512).
_ENSEMBLE_SUFFIXES = (".csv", ".jcb")


def worker_root(master) -> Path:
    """The directory a parallel run's workers are made in: ``<master>.workers``.

    Parameters
    ----------
    master : str or pathlib.Path
        The run's master directory.

    Returns
    -------
    pathlib.Path
        ``master`` with :data:`WORKERS_SUFFIX` appended to its name, beside it.
    """

    master = Path(master)
    return master.with_name(master.name + WORKERS_SUFFIX)


# -- the lock -----------------------------------------------------------------------


def _process_start(pid: int) -> str | None:
    """When process ``pid`` started, as the kernel records it (Linux), or ``None``.

    Field 22 of ``/proc/<pid>/stat``, in clock ticks since boot: with the pid it tells a
    process from a later one the system gave the same number.
    """

    try:
        stat = Path(f"/proc/{pid}/stat").read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    # The command name (field 2) is in parentheses and may hold spaces or ")".
    fields = stat[stat.rfind(")") + 2 :].split()
    return fields[19] if len(fields) > 19 else None


def write_run_lock(directory, *, method: str, calibration: str) -> Path:
    """Write :data:`RUN_LOCK` into ``directory``, naming this process.

    Parameters
    ----------
    directory : str or pathlib.Path
        The directory the run is going in; it must exist.
    method : str
        The launching method, ``"run_ies"``, ``"prior"`` or ``"resume_ies"``, named in
        a refusal.
    calibration : str
        The calibration's name, named in a refusal.

    Returns
    -------
    pathlib.Path
        The lock file written.
    """

    directory = Path(directory).resolve()
    pid = os.getpid()
    path = directory / RUN_LOCK
    path.write_text(
        json.dumps(
            {
                "written_by": _LOCK_WRITER,
                "directory": str(directory),
                "host": socket.gethostname(),
                "pid": pid,
                "process_start": _process_start(pid),
                "started_at": time.strftime("%Y-%m-%d %H:%M:%S"),
                "method": method,
                "calibration": calibration,
                "note": _LOCK_NOTE,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return path


def release_run_lock(directory) -> None:
    """Remove the :data:`RUN_LOCK` THIS process wrote in ``directory``, if it is there.

    A lock another process wrote (a run started since in another kernel) is left.

    Parameters
    ----------
    directory : str or pathlib.Path
        The directory the run went in.
    """

    lock = _read_lock(Path(directory).resolve())
    if lock is not None and lock.get("pid") == os.getpid() and lock.get("host") == socket.gethostname():
        try:
            (Path(directory) / RUN_LOCK).unlink()
        except FileNotFoundError:
            pass


def _read_lock(directory: Path) -> dict | None:
    """The :data:`RUN_LOCK` naming ``directory`` (resolved), or ``None``.

    A lock naming ANOTHER directory is a copy -- pyEMU copies the template, a serial
    run's lock and all, into each worker -- and is not this directory's.
    """

    try:
        lock = json.loads((directory / RUN_LOCK).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        # Absent, unreadable, not UTF-8 or not JSON: no lock.
        return None
    if (
        not isinstance(lock, dict)
        or lock.get("written_by") != _LOCK_WRITER
        or lock.get("directory") != str(directory)
    ):
        return None
    return lock


def _process_is_running(pid: int, start: str | None) -> bool | None:
    """Whether process ``pid`` -- the one that started at ``start`` -- still runs here.

    ``None`` when this platform gives no way to ask without harm: on Windows
    ``os.kill(pid, 0)`` TERMINATES the process.
    """

    if sys.platform.startswith("win"):
        return None
    try:
        pid = int(pid)
    except (TypeError, ValueError):
        return False
    if pid <= 0:
        return False   # kill(0 or -n, 0) would ask about a whole process group
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        pass   # it exists, owned by someone else
    except OverflowError:
        return False
    if start is not None:
        now = _process_start(int(pid))
        if now is not None and now != start:
            return False   # the pid was given to a later process
    return True


def pestpp_processes_in(directory) -> tuple[tuple[int, str], ...]:
    """The PEST++ processes whose working directory is ``directory`` or inside it.

    Read from ``/proc`` (Linux); empty elsewhere, and for a process whose working
    directory this user may not read. Only processes whose command name starts
    ``pestpp`` count -- a shell someone left in the folder is not a run.

    Parameters
    ----------
    directory : str or pathlib.Path
        The run directory (a master, a worker root, a template).

    Returns
    -------
    tuple of (int, str)
        ``(pid, command name)`` of each, sorted by pid.
    """

    directory = Path(directory).resolve()
    proc = Path("/proc")
    if not proc.is_dir():
        return ()
    found = []
    try:
        entries = list(proc.iterdir())
    except OSError:
        return ()
    for entry in entries:
        if not entry.name.isdigit() or int(entry.name) == os.getpid():
            continue
        try:
            name = (entry / "comm").read_text(encoding="utf-8", errors="replace").strip()
            if not name.lower().startswith("pestpp"):
                continue
            cwd = Path(os.readlink(entry / "cwd"))
        except OSError:
            # Gone since the listing, or another user's: not readable, not counted.
            continue
        if cwd == directory or cwd.is_relative_to(directory):
            found.append((int(entry.name), name))
    return tuple(sorted(found))


@dataclass(frozen=True)
class RunActivity:
    """Whether a PESTPP-IES run is going in a directory, as :func:`run_activity` reads it.

    Attributes
    ----------
    directory : pathlib.Path
        The directory asked about (resolved).
    state : str
        ``"free"`` (no lock, no PEST++ process), ``"live"`` (the lock's process is
        running, or a PEST++ process works there), ``"unknown"`` (a lock this machine
        cannot check: written on another host, or on Windows -- treated as live) or
        ``"stale"`` (a lock whose process is gone, and no PEST++ process there).
    lock : dict or None
        The :data:`RUN_LOCK` naming the directory, or ``None``.
    processes : tuple of (int, str)
        The PEST++ processes working in it (:func:`pestpp_processes_in`).
    lock_running : bool or None
        Whether the process the lock names still runs: ``None`` without a lock, or
        when this machine cannot tell.
    """

    directory: Path
    state: str
    lock: dict | None
    processes: tuple[tuple[int, str], ...] = ()
    lock_running: bool | None = None

    @property
    def blocks(self) -> bool:
        """Whether the directory must not be deleted or reused: ``live`` or ``unknown``."""

        return self.state in ("live", "unknown")

    def describe(self) -> str:
        """One sentence saying what holds the directory, for a refusal or a warning."""

        lock = self.lock or {}
        started = (
            f"{lock.get('method', 'a run')}() of calibration {lock.get('calibration')!r}, "
            f"started {lock.get('started_at', '?')} by process {lock.get('pid')} on host "
            f"{lock.get('host')!r}"
        )
        processes = ", ".join(f"{name} (pid {pid})" for pid, name in self.processes)
        if self.state == "free":
            return "no run holds it."
        if self.state == "stale":
            return (
                f"its {RUN_LOCK} names {started}, which is no longer running: the lock is "
                "STALE (that run ended without removing it -- its kernel was killed, say)."
            )
        if self.state == "unknown":
            return (
                f"its {RUN_LOCK} names {started}, which this machine cannot check, so it "
                "counts as a run still going."
            )
        if self.lock_running:
            return f"a run is going in it: {started}" + (
                f"; PEST++ working there: {processes}" if processes else ""
            ) + "."
        return f"PEST++ is running in it: {processes}" + (
            f" (its {RUN_LOCK} names {started}, which has exited)" if self.lock else ""
        ) + "."


def run_activity(directory) -> RunActivity:
    """Whether a PESTPP-IES run is going in ``directory`` (ledger 532).

    Parameters
    ----------
    directory : str or pathlib.Path
        A template, master or worker root. One that does not exist is ``"free"``.

    Returns
    -------
    RunActivity
        Its :attr:`~RunActivity.state` and what decided it.
    """

    directory = Path(directory).resolve()
    if not directory.is_dir():
        return RunActivity(directory, "free", None)
    lock = _read_lock(directory)
    processes = pestpp_processes_in(directory)
    if lock is None:
        return RunActivity(directory, "live" if processes else "free", None, processes)
    running = None
    if lock.get("host") == socket.gethostname():
        running = _process_is_running(lock.get("pid", -1), lock.get("process_start"))
    if running or processes:
        state = "live"
    elif running is None:
        state = "unknown"
    else:
        state = "stale"
    return RunActivity(directory, state, lock, processes, running)


# -- resuming an interrupted IES run --------------------------------------------------


def _iteration_files(directory: Path, case: str) -> dict[int, list[Path]]:
    """Each iteration's files ``<case>.<n>.*`` in ``directory``, by ``n``."""

    pattern = re.compile(re.escape(case) + r"\.(\d+)\..+")
    found: dict[int, list[Path]] = {}
    for path in pestpp_ies_outputs(directory, case):
        match = pattern.fullmatch(path.name)
        if match:
            found.setdefault(int(match.group(1)), []).append(path)
    return found


def ensemble_file(directory, stem: str, *, binary_first: bool = False) -> Path | None:
    """``<stem>.csv`` or ``<stem>.jcb`` in ``directory``, or ``None``.

    Parameters
    ----------
    directory : str or pathlib.Path
        Where to look.
    stem : str
        The ensemble's name without its suffix, e.g. ``cal.2.par``.
    binary_first : bool, optional
        Prefer ``.jcb`` when both exist (the run set ``ies_save_binary``). Default
        ``False``.

    Returns
    -------
    pathlib.Path or None
        The file found, or ``None``.
    """

    suffixes = tuple(reversed(_ENSEMBLE_SUFFIXES)) if binary_first else _ENSEMBLE_SUFFIXES
    for suffix in suffixes:
        path = Path(directory) / f"{stem}{suffix}"
        if path.is_file():
            return path
    return None


def recorded_iterations(directory, case: str) -> list[int]:
    """The iterations the run's ``<case>.phi.actual.csv`` records, sorted; ``[]`` without one.

    Parameters
    ----------
    directory : str or pathlib.Path
        A run directory.
    case : str
        The control-file stem the run used.

    Returns
    -------
    list of int
        The ``iteration`` column's values.
    """

    path = Path(directory) / f"{case}.phi.actual.csv"
    try:
        column = pd.read_csv(path, usecols=["iteration"])["iteration"]
    except (OSError, ValueError) as error:
        # Absent, empty, or a phi file with no `iteration` column: nothing recorded.
        logger.debug("no iteration record in %s: %s", path, error)
        return []
    return sorted({int(value) for value in pd.to_numeric(column, errors="coerce").dropna()})


def last_completed_iteration(directory, case: str, *, binary_first: bool = False) -> int | None:
    """The last iteration a run in ``directory`` completed, or ``None`` when none did.

    An iteration counts when the run's phi file records it AND its parameter and
    observation ensembles (``<case>.<n>.par`` / ``.obs``, CSV or ``.jcb``) are on disk --
    what a resume restarts from.

    Parameters
    ----------
    directory : str or pathlib.Path
        The run directory (a master, or the template for a serial run).
    case : str
        The control-file stem the run used.
    binary_first : bool, optional
        As for :func:`ensemble_file`. Default ``False``.

    Returns
    -------
    int or None
        The iteration number, ``0`` being the evaluated prior.
    """

    complete = [
        n
        for n in recorded_iterations(directory, case)
        if ensemble_file(directory, f"{case}.{n}.par", binary_first=binary_first)
        and ensemble_file(directory, f"{case}.{n}.obs", binary_first=binary_first)
    ]
    return max(complete) if complete else None


_LAMBDA = re.compile(r"(?:current lambda:|updating lambda to)\s*([-+0-9.eE]+)")


def lambda_in_force(record, iteration: int) -> float | None:
    """The Marquardt lambda PESTPP-IES would have started iteration ``iteration + 1`` with.

    Read from the run's ``.rec``: the last ``current lambda:`` or ``updating lambda to``
    it printed before ``starting solve for iteration: <iteration + 1>`` (the whole record
    when that iteration never started). A resume that does not carry it starts from
    PESTPP-IES's own initial lambda, and its iterations follow another path.

    Parameters
    ----------
    record : str or pathlib.Path
        The interrupted run's ``<case>.rec``.
    iteration : int
        The last completed iteration.

    Returns
    -------
    float or None
        The lambda, or ``None`` when the record is unreadable or names none.
    """

    try:
        text = Path(record).read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    cut = re.search(rf"starting solve for iteration:\s*{int(iteration) + 1}\b", text)
    found = _LAMBDA.findall(text[: cut.start()] if cut else text)
    try:
        return float(found[-1]) if found else None
    except ValueError:
        return None


def pending_resume(directory, case: str) -> Path | None:
    """An archive a resume of ``case`` made in ``directory`` and never spliced back, or ``None``.

    What a resume killed too hard for its ``finally`` to run leaves behind.

    Parameters
    ----------
    directory : str or pathlib.Path
        The run directory.
    case : str
        The control-file stem.

    Returns
    -------
    pathlib.Path or None
        The archive directory.
    """

    for archive in sorted(Path(directory).glob(f"{case}.resumed_from_*")):
        try:
            record = json.loads((archive / RESUME_RECORD).read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if isinstance(record, dict) and not record.get("spliced"):
            return archive
    return None


def archive_for_resume(directory, case: str, iteration: int) -> Path:
    """Move a run's outputs aside before PESTPP-IES restarts from ``iteration``.

    Every file :func:`~myflopy.modflow.mf6.pest.runs.pestpp_ies_outputs` lists, and the
    run's ``<case>.rec`` and ``<case>.log``, go into ``<case>.resumed_from_<n>/`` (a
    numbered ``_2``, ``_3`` .. when a resume from the same iteration was made before),
    with a :data:`RESUME_RECORD` naming the iteration -- so the restarted run, which
    numbers its iterations from 0 again, overwrites none of them.

    Parameters
    ----------
    directory : str or pathlib.Path
        The run directory.
    case : str
        The control-file stem.
    iteration : int
        The last completed iteration, which the restart starts from.

    Returns
    -------
    pathlib.Path
        The archive directory.
    """

    directory = Path(directory)
    base = directory / f"{case}.resumed_from_{iteration}"
    archive, number = base, 1
    while archive.exists():
        number += 1
        archive = base.with_name(f"{base.name}_{number}")
    archive.mkdir()
    moved = list(pestpp_ies_outputs(directory, case))
    moved += [directory / f"{case}.{suffix}" for suffix in ("rec", "log")
              if (directory / f"{case}.{suffix}").is_file()]
    for path in moved:
        shutil.move(str(path), str(archive / path.name))
    (archive / RESUME_RECORD).write_text(
        json.dumps({"case": case, "from_iteration": int(iteration), "spliced": False}, indent=2),
        encoding="utf-8",
    )
    return archive


def _join_phi(archived: Path, resumed: Path, out: Path, iteration: int) -> None:
    """Write the interrupted run's phi rows up to ``iteration``, then the resumed run's after it.

    The resumed run's iteration ``k >= 1`` is iteration ``iteration + k``, and its
    ``total_runs`` (cumulative) continue from the interrupted run's at ``iteration``.
    Its iteration 0 restates the restart ensembles, which the interrupted run's own
    row ``iteration`` already holds, and is dropped.
    """

    before = pd.read_csv(archived)
    before = before[pd.to_numeric(before["iteration"], errors="coerce") <= iteration]
    if resumed.is_file():
        after = pd.read_csv(resumed)
        after = after[pd.to_numeric(after["iteration"], errors="coerce") >= 1].copy()
        after["iteration"] = after["iteration"].astype(int) + int(iteration)
        if "total_runs" in after and "total_runs" in before and len(before):
            done = pd.to_numeric(before["total_runs"], errors="coerce").max()
            after["total_runs"] = after["total_runs"] + (0 if pd.isna(done) else int(done))
        joined = pd.concat([before, after], ignore_index=True, sort=False)
    else:
        joined = before
    joined.to_csv(out, index=False)


def splice_resumed_run(directory, case: str, archive) -> list[int]:
    """Join a resumed run's outputs to the interrupted run's, as one run of iterations ``0..N``.

    PESTPP-IES numbered the resumed run's iterations from 0, its iteration 0 being the
    restart ensembles (the interrupted run's last completed iteration ``n``, read from
    the archive). This deletes the resumed iteration 0 (a copy of ``n``), renumbers each
    resumed iteration ``k >= 1`` to ``n + k``, moves the interrupted run's iterations
    ``0..n`` back, and joins each phi file (:func:`_join_phi`). The noise ensemble and
    the final ``.base.par``/``.base.rei`` are the resumed run's when it wrote them, else
    the archived ones come back. What stays in the archive -- the interrupted run's
    ``.rec``/``.log``, and any file of an iteration past ``n`` it had begun -- is kept,
    with its record marked spliced. Safe to call on a resume that wrote nothing.

    Parameters
    ----------
    directory : str or pathlib.Path
        The run directory.
    case : str
        The control-file stem.
    archive : str or pathlib.Path
        The directory :func:`archive_for_resume` made.

    Returns
    -------
    list of int
        The iterations the joined phi file records.
    """

    directory, archive = Path(directory), Path(archive)
    record_path = archive / RESUME_RECORD
    record = json.loads(record_path.read_text(encoding="utf-8"))
    n = int(record["from_iteration"])

    resumed = _iteration_files(directory, case)
    for path in resumed.pop(0, []):
        path.unlink()
    for k in sorted(resumed, reverse=True):   # descending: n + k never lands on a k' > k
        for path in resumed[k]:
            rest = path.name[len(f"{case}.{k}.") :]
            path.rename(directory / f"{case}.{n + k}.{rest}")
    for k, paths in _iteration_files(archive, case).items():
        if k <= n:
            for path in paths:
                shutil.move(str(path), str(directory / path.name))

    for archived in sorted(archive.glob(f"{case}.phi.*.csv")):
        resumed_phi = directory / archived.name
        _join_phi(archived, resumed_phi, resumed_phi, n)
        archived.unlink()
    for archived in pestpp_ies_outputs(archive, case):
        if re.fullmatch(re.escape(case) + r"\.\d+\..+", archived.name):
            continue   # an iteration past n the interrupted run had begun: kept here
        target = directory / archived.name
        if target.exists():
            archived.unlink()   # the resumed run's own (noise ensemble, final base)
        else:
            shutil.move(str(archived), str(target))

    record["spliced"] = True
    record_path.write_text(json.dumps(record, indent=2), encoding="utf-8")
    return recorded_iterations(directory, case)
