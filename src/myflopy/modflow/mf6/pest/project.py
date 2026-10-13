"""Project orchestration for the first ``myflopy`` PEST slice."""

from __future__ import annotations

import contextlib
import io
import json
import os
import re
import shutil
import sys
import time
import warnings
from collections.abc import Mapping, Sequence
from pathlib import Path

import flopy
from flopy.mf6.mfbase import FlopyException, MFDataException

from myflopy._logging import get_logger
from myflopy._optional import require
from myflopy._pythonpath import pythonpath_for_children
from myflopy.modflow.mf6.lak_implicit import keep_implicit
from myflopy.modflow.mf6.observations import (
    ConcDifferences,
    ConcTargets,
    DrnFlowTargets,
    HeadDifferences,
    HeadTargets,
    LakeStageTargets,
    RivFlowTargets,
    SfrFlowTargets,
    SfrStageTargets,
)
from myflopy.modflow.mf6.observations.registry import bare_targets
from myflopy.modflow.mf6.package_registry import model_packages_of_type
from myflopy.modflow.mf6.period_end import period_lengths
from myflopy.modflow.mf6.pest.geostats import build_geostruct, checked_nugget, unit_sill
from myflopy.modflow.mf6.pest.native_parameters import (
    DEFAULT_PP_NEAREST,
    NativeParameterSpec,
    add_native_parameter,
    flopy_model_for,
    list_columns,
    model_name_for,
    pyemu_spatial_reference,
    refuse_an_undeclared_array,
    relayer_array_target,
    resolve_features,
    resolve_target,
    select_packages,
    staged_model,
)
from myflopy.modflow.mf6.pest.observations import (
    finalize_observations,
    prepare_conc_observations,
    prepare_difference_observations,
    prepare_drn_flow_observations,
    prepare_head_target_observations,
    prepare_lake_stage_observations,
    prepare_riv_flow_observations,
    prepare_sfr_flow_observations,
    prepare_sfr_stage_observations,
)
from myflopy.modflow.mf6.pest.pilot_points import (
    add_pilot_point_parameter,
    emit_pilot_point_commands,
    register_pilot_point_parameters,
)
from myflopy.modflow.mf6.pest.run_state import (
    RESUME_RECORD,
    RUN_LOCK,
    archive_for_resume,
    ensemble_file,
    lambda_in_force,
    last_completed_iteration,
    pending_resume,
    recorded_iterations,
    release_run_lock,
    run_activity,
    splice_resumed_run,
    worker_root,
    write_run_lock,
)
from myflopy.modflow.mf6.pest.runs import (
    BUILD_MARKER,
    default_pest_root,
    pestpp_ies_outputs,
)
from myflopy.modflow.mf6.pest.specs import (
    ConcDifferenceObservationSpec,
    ConcObservationSpec,
    DrnFlowObservationSpec,
    ExpGeoStruct,
    HeadDifferenceObservationSpec,
    HeadTargetObservationSpec,
    LakeStageObservationSpec,
    RivFlowObservationSpec,
    SfrFlowObservationSpec,
    SfrStageObservationSpec,
)
from myflopy.modflow.mf6.pest.summary import PestSettings

METADATA_FILENAME = "myflopy_pest_metadata.json"

#: The log ``PstFrom`` opens in the working directory; the build moves it into the template.
PYEMU_LOG = "PstFrom.log"

#: Appended to the template directory's NAME for the build's staging copy of the
#: run (ledger 356): a sibling of the template, removed when the build ends.
STAGING_SUFFIX = ".staging"

# BUILD_MARKER (imported from `runs`, where run discovery reads it too) is the file a
# build writes into its template, and into its staging copy, naming the directory it
# was written for and that directory's role: the one sign by which a later build may
# DELETE that directory (ledger 360). A build deletes a directory only when it is
# absent, empty, or carries this file naming it -- a folder the build did not write is
# refused, whatever else it holds. VISIBLE, not a dot-file (ledger 361): a hidden
# marker survived ``rm <template>/*`` and Finder's select-all delete, so a folder
# emptied that way and refilled with the user's files still carried it, and the next
# build deleted them.

#: The ``written_by`` value of a :data:`BUILD_MARKER`.
_MARKER_WRITER = "myflopy PEST build"

#: What a :data:`BUILD_MARKER` says about itself, for the person who finds one.
_MARKER_NOTE = (
    "Written by a myflopy PEST build. While this file names this directory, a later "
    "build may DELETE the directory and everything in it (ledger 360). Delete this "
    "file to keep the directory from being deleted."
)

logger = get_logger(__name__)

# Pre-built observation specs accepted by observe()/forecast() as-is.
_OBSERVATION_SPEC_TYPES = (
    HeadTargetObservationSpec,
    ConcObservationSpec,
    HeadDifferenceObservationSpec,
    ConcDifferenceObservationSpec,
    LakeStageObservationSpec,
    SfrStageObservationSpec,
    SfrFlowObservationSpec,
    DrnFlowObservationSpec,
    RivFlowObservationSpec,
)

# High-level target sets and the observation spec each is wrapped in, so
# observe()/forecast() accept head, lake-stage, SFR stage/flow, DRN seepage and RIV
# exchange
# targets uniformly (each spec carries its own default prefix).
_TARGET_SPEC_TYPES = (
    # ConcTargets subclasses HeadTargets, so it MUST come first -- an isinstance
    # ladder would otherwise wrap every concentration target as a head spec.
    (ConcTargets, ConcObservationSpec),
    (HeadTargets, HeadTargetObservationSpec),
    # Differenced sets (ledger 565); ConcDifferences subclasses
    # HeadDifferences, so it comes first for the same reason.
    (ConcDifferences, ConcDifferenceObservationSpec),
    (HeadDifferences, HeadDifferenceObservationSpec),
    (LakeStageTargets, LakeStageObservationSpec),
    (SfrStageTargets, SfrStageObservationSpec),
    (SfrFlowTargets, SfrFlowObservationSpec),
    (DrnFlowTargets, DrnFlowObservationSpec),
    # Neither zone family subclasses the other (both subclass a private base).
    (RivFlowTargets, RivFlowObservationSpec),
)


#: How many standard deviations a parameter's bounds span in every prior draw: a
#: log-transformed parameter of bounds ``(lb, ub)`` is drawn with a log10 standard
#: deviation of ``log10(ub / lb) / PRIOR_SIGMA_RANGE``. ``6`` is ``PstFrom.draw``'s
#: own default, which the correlated grid draw used all along; the pilot-point draw
#: now uses it too (ledger 424).
PRIOR_SIGMA_RANGE = 6

#: The ``overdue_giveup_fac`` :meth:`PestProject.run_ies` and :meth:`PestProject.prior`
#: hand PESTPP-IES unless told otherwise: a parallel model run taking longer than this
#: many times the average run is killed and its realization counted failed. PESTPP-IES
#: resets an unset factor to 2.0, which on a shared machine kills a slow NEWTON run
#: and drops its realization -- selectively, the hard-to-converge ones (ledger
#: 510).
DEFAULT_OVERDUE_GIVEUP_FAC = 5.0

#: ``PstFrom``'s ``chunk_len``: how many multiplier files one process of a forward
#: run's ``apply_list_and_array_pars`` takes. pyEMU's default, 50, splits 50 or more
#: files into chunks and applies them in a ``spawn`` pool of ``min(cpu_count, 60)``
#: processes -- in EVERY realization, beside every other worker's (measured at 61
#: files: 21 child processes, ~3.4 GB PSS and 4.4 s against 1.3 s in one process;
#: ledger 511). Parallelism belongs to the IES workers, so one chunk.
PYEMU_APPLY_CHUNK_LEN = 10**6


def _pest_run_slug(name: str) -> str:
    """Filesystem-safe directory name for a PEST run's default workspace."""

    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(name)).strip("_")
    return slug or "pest_run"


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

    data = path.read_bytes()
    if b"\0" in data[:8192]:
        return None
    lines = data.decode("utf-8", errors="replace").splitlines()
    return "\n".join(re.split(r"[#!]", line, maxsplit=1)[0] for line in lines)


def _simulation_input_files(sim_ws: str | Path) -> list[Path]:
    """Every file MODFLOW 6 READS when it runs the simulation in ``sim_ws``, relative to it.

    Followed from ``mfsim.nam`` the way MF6 opens them: each name-file entry (TDIS,
    each model's name file and its packages, exchanges, solutions), each ``FILEIN``
    a scanned file names (TS6, OBS6, TAS6, TAB6, ATS6, SPC6, HPC6, FMI's flow files
    ...) and each ``OPEN/CLOSE`` external array or list -- names MF6 resolves against
    the directory it runs in. A ``FILEIN`` file is scanned in turn; an ``OPEN/CLOSE``
    file holds data and is not. Outputs -- the listing files, every ``FILEOUT``, and
    any file nothing names -- are not inputs, and neither is a name that is absolute
    (read in place from anywhere), outside ``sim_ws`` (see
    :func:`_inputs_outside_the_run`), or not a file there.

    This is what a PEST build's staging copy holds (ledger 383): the template
    needs the simulation's inputs and nothing else of the run.

    Parameters
    ----------
    sim_ws : str or pathlib.Path
        The simulation directory, holding ``mfsim.nam``.

    Returns
    -------
    list of pathlib.Path
        The input files, relative to ``sim_ws``, ``mfsim.nam`` first, each once.
    """

    return _scan_simulation_inputs(sim_ws)[0]


def _inputs_outside_the_run(sim_ws: str | Path) -> list[tuple[str, str]]:
    """Each RELATIVE input name that leads outside ``sim_ws``, with what names it.

    ``../gwf/gwf.hds`` from a GWT model's FMI, an ``OPEN/CLOSE ../arrays/k.txt``, a
    model name file one directory up: MF6 resolves such a name against the
    directory it runs in, and every PEST realization runs in its own copy of the
    template, so there it names a file beside the template -- a missing one, or
    someone else's (ledger 421). An ABSOLUTE name is not one of these: it is
    read in place from every copy.

    Parameters
    ----------
    sim_ws : str or pathlib.Path
        The simulation directory, holding ``mfsim.nam``.

    Returns
    -------
    list of tuple of (str, str)
        ``(name as written, what names it)`` -- the name file or package file, and
        the package it belongs to -- in the order MF6 meets them.
    """

    return _scan_simulation_inputs(sim_ws)[1]


def _outputs_outside_the_run(sim_ws: str | Path) -> list[tuple[str, str]]:
    """Each output name -- a ``FILEOUT``, or a model's ``LIST`` -- that writes outside ``sim_ws``.

    A RELATIVE name leading out of the directory (``../heads.hds``) and an ABSOLUTE
    name, wherever it points: MF6 writes either one outside the directory it runs
    in, so every PEST realization -- each in its own copy of the template, several
    at once under parallel workers -- writes the SAME file, and an absolute name
    into the run directory overwrites the run's own output (ledger 453). Read
    from every file :func:`_simulation_input_files` scans, and from each model name
    file's ``OPTIONS`` ``LIST <file>``, the listing file.

    Parameters
    ----------
    sim_ws : str or pathlib.Path
        The simulation directory, holding ``mfsim.nam``.

    Returns
    -------
    list of tuple of (str, str)
        ``(name as written, what names it)``, in the order MF6 meets them.
    """

    return _scan_simulation_inputs(sim_ws)[3]


def _output_directories(sim_ws: str | Path) -> list[Path]:
    """The directories inside ``sim_ws`` that a ``FILEOUT`` or a model's ``LIST`` writes into.

    MF6 does not create a directory: ``OC ... BUDGET FILEOUT out/m.cbc`` stops the run
    with ``Could not open "out/m.cbc"`` unless ``out/`` exists, and a name file's
    ``LIST out/m.lst`` stops it the same way. The PEST staging copy
    holds the run's INPUTS, so a directory holding only outputs was missing from the
    template and every realization failed (ledger 453); the build makes these.

    Parameters
    ----------
    sim_ws : str or pathlib.Path
        The simulation directory, holding ``mfsim.nam``.

    Returns
    -------
    list of pathlib.Path
        Each directory once, relative to ``sim_ws``, ``sim_ws`` itself excluded.
    """

    return list(dict.fromkeys(
        name.parent for name in _scan_simulation_inputs(sim_ws)[2] if name.parent != Path(".")
    ))


def _scan_simulation_inputs(
    sim_ws: str | Path,
) -> tuple[list[Path], list[tuple[str, str]], list[Path], list[tuple[str, str]]]:
    """The walk behind :func:`_simulation_input_files`, :func:`_inputs_outside_the_run`,
    :func:`_output_directories` and :func:`_outputs_outside_the_run`: the inputs, the
    inputs named outside, the outputs written inside (relative) and those outside."""

    root = Path(sim_ws)
    found: dict[Path, None] = {}
    outside: list[tuple[str, str]] = []
    written: dict[Path, None] = {}
    written_outside: list[tuple[str, str]] = []
    # (relative name, how the file is read, what it is): "simulation" / "model"
    # name files are read for their blocks, "input" files scanned for
    # FILEIN/OPEN/CLOSE, "data" files (OPEN/CLOSE) copied unread.
    pending: list[tuple[str, str, str]] = [("mfsim.nam", "simulation", "mfsim.nam")]

    def _record_output(name_written: str, owner: str) -> None:
        """File an output name as written inside ``sim_ws`` or outside it."""

        normal = Path(os.path.normpath(name_written))
        if normal.is_absolute() or normal.parts[:1] == ("..",):
            written_outside.append((name_written, owner))
        else:
            written[normal] = None

    while pending:
        name, kind, owner = pending.pop(0)
        relative = Path(os.path.normpath(name))
        if relative.is_absolute() or relative in found:
            continue
        if relative.parts[:1] == ("..",):
            outside.append((name, owner))
            continue
        path = root / relative
        if not path.is_file():
            continue
        found[relative] = None
        if kind == "data":
            continue
        text = _mf6_text(path)
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
                    ))
        named_by = "" if kind in ("simulation", "model") else f", named by {owner}"
        for keyword, token in _MF6_INPUT_REFERENCE.findall(text):
            pending.append((
                _mf6_field(token),
                "input" if keyword.upper() == "FILEIN" else "data",
                f"{keyword.upper()} in {relative}{named_by}",
            ))
        for token in _MF6_OUTPUT_REFERENCE.findall(text):
            _record_output(_mf6_field(token), f"FILEOUT in {relative}{named_by}")
    return list(found), outside, list(written), written_outside


def _record_control_file(archive: Path, noptmax: int, options: dict) -> None:
    """Add to ``archive``'s resume record what the run's control file said before the restart.

    ``noptmax`` is the run's TOTAL iterations and ``options`` maps each restart option
    to its value then (``None`` when unset); :meth:`PestProject._restore_resumed_control_file`
    writes them back (review of ledger 533).
    """

    path = archive / RESUME_RECORD
    record = json.loads(path.read_text(encoding="utf-8"))
    record["control_file"] = {
        "noptmax": int(noptmax),
        "options": {
            key: (None if value is None else value if isinstance(value, (bool, int, float)) else str(value))
            for key, value in options.items()
        },
    }
    path.write_text(json.dumps(record, indent=2), encoding="utf-8")


def build_forward_run_command(template_workspace: Path, interpreter: str) -> str:
    """Return the PEST++ model command that runs ``forward_run.py``.

    The forward-run interpreter is pinned to ``interpreter`` (this
    environment's Python) so PEST++ workers use the venv that has
    numpy/flopy/pyemu/myflopy, not a bare ``python`` that may resolve
    elsewhere. On Windows the quoted absolute path works directly. On POSIX,
    PEST++'s run manager neither shell-parses quotes nor accepts absolute
    command paths (it mangles the leading ``/`` and ``execv`` fails), so the
    command is a ``./run_forward.sh`` wrapper written into the template
    workspace that ``exec``s the absolute interpreter; PEST++ worker dirs are
    copies of the template, so the wrapper travels with them.
    """

    if os.name == "nt":
        return f'"{interpreter}" forward_run.py'
    wrapper = Path(template_workspace) / "run_forward.sh"
    wrapper.write_text(f'#!/bin/sh\nexec "{interpreter}" forward_run.py "$@"\n')
    wrapper.chmod(0o755)
    return "./run_forward.sh"


def _write_build_marker(directory: Path, role: str) -> None:
    """Mark ``directory`` as one a build wrote, in ``role`` (``"template"`` / ``"staging"``).

    The file names the directory's RESOLVED path, so a copy of it -- an IES master or
    worker directory, a backup -- is not taken for it (ledger 360).
    """

    directory = Path(directory).resolve()
    (directory / BUILD_MARKER).write_text(
        json.dumps(
            {
                "written_by": _MARKER_WRITER,
                "role": role,
                "directory": str(directory),
                "note": _MARKER_NOTE,
            },
            indent=2,
        ),
        encoding="utf-8",
    )


def _build_marker(directory: Path) -> dict | None:
    """The :data:`BUILD_MARKER` in ``directory``, or ``None`` when it holds no readable one."""

    try:
        marker = json.loads((Path(directory) / BUILD_MARKER).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        # Absent, a directory, unreadable, not UTF-8 or not JSON: no build wrote it.
        return None
    if not isinstance(marker, dict) or marker.get("written_by") != _MARKER_WRITER:
        return None
    return marker


def _written_by_a_build(directory: Path, role: str, *, naming: Path | None = None) -> bool:
    """Whether ``directory`` holds a :data:`BUILD_MARKER` in ``role`` naming ``naming``.

    ``naming`` defaults to ``directory`` itself (resolved by the caller): the directory a
    build wrote. Another path asks whether ``directory`` is a COPY of the one named --
    pyEMU's ``start_workers`` copies a template, marker and all, into the IES master and
    each worker (ledger 361).
    """

    marker = _build_marker(directory)
    return (
        marker is not None
        and marker.get("role") == role
        and marker.get("directory") == str(directory if naming is None else naming)
    )


def _relation(directory: Path, other: Path) -> str | None:
    """How resolved ``directory`` overlaps resolved ``other``: ``"IS"``, ``"CONTAINS"``, ``"sits INSIDE"`` or ``None``."""

    if directory == other:
        return "IS"
    if other.is_relative_to(directory):
        return "CONTAINS"
    if directory.is_relative_to(other):
        return "sits INSIDE"
    return None


@contextlib.contextmanager
def _working_directory(path: Path):
    """Run a block with ``path`` as the working directory, restoring the old one after.

    For ``PstFrom``'s constructor only: it opens ``PstFrom.log`` in the working
    directory and takes no argument to put it elsewhere (ledger 356). The working
    directory is the PROCESS's, so every path the block uses must be absolute.
    """

    previous = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def _import_pyemu():
    """Import pyEMU lazily with a workflow-oriented error message."""

    return require("pyemu", feature="PestProject workflows")


def _pyemu_warning_class(pyemu_module):
    """Return the installed ``PyemuWarning`` class if available."""

    try:
        return pyemu_module.pyemu_warnings.PyemuWarning
    except AttributeError:
        return Warning


class PestProject:
    """Build a pyEMU/PEST++ calibration workspace from a ``myflopy`` model.

    It compiles a few readable declarations straight to native
    ``pyemu.utils.PstFrom`` and lets pyEMU drive the forward run.

    **The front door is** ``model.pest("calib")``, which returns one of these
    with ``start_datetime`` filled in from the model's TDIS. Constructing
    ``PestProject(...)`` directly is the advanced path: everything else is the
    same -- ``workspace`` still defaults to ``<model workspace>.pest/<name>``, so
    the run is still discoverable through ``model.pest_runs`` -- but
    ``start_datetime`` has no TDIS fallback and you must supply it::

        cal = model.pest("calib", start_datetime="2020-01-01")
        cal.parameterize("k",        style="pilotpoints", pp_spacing=800.0, physical=(1e-3, 100))
        cal.parameterize("recharge", bounds=(0.5, 1.5))
        cal.parameterize("ghb.cond", bounds=(0.1, 10))
        cal.observe(head_targets)          # history-matching targets
        cal.forecast(prediction_targets)   # predictions of interest (zero weight)
        print(cal.settings())              # review the resolved configuration
        pst = cal.build("calib.pst")       # -> native .pst + forward_run.py
        cal.run_ies()                      # PESTPP-IES ensemble smoother

    Use :meth:`parameterize` (styles ``constant`` / ``zone`` / ``grid`` /
    ``pilotpoints`` across every target), :meth:`observe`, :meth:`forecast`,
    :meth:`build`, :meth:`run_ies`/:meth:`prior` (and :meth:`resume_ies` for an
    interrupted history match) and :meth:`settings`. The runs
    land beside the model and are reopened for review via ``model.pest_runs``
    (or ``run.pest_runs``) -> :meth:`~...runs.PestRunHandle.review`
    (:class:`~...ies.IesResults`).

    Parameters
    ----------
    model : SimulationBase
        A built ``myflopy`` model (``SimulationBase`` or compatible) exposing
        ``.sim`` (the FloPy simulation), ``.gwf``, ``.name``,
        ``.workspace`` and, for explicit ``pp_points=`` pilot points, any ``zones=`` and
        head or concentration targets given by cell id, ``.vor``. Its workspace
        -- the RUN directory, which must hold the simulation's files
        (``mfsim.nam``) -- is READ, never written: :meth:`build` copies the
        simulation's input files from it -- what ``mfsim.nam`` leads MF6 to read,
        none of its outputs, run metadata or old PEST folders (ledger 383) --
        to a staging directory beside the template, loads the simulation from that
        copy and stores every array and list there in external files
        (``set_all_data_external``), so the run's files and this object are
        left exactly as they were (ledger 356). The template is made from the
        run's FILES, not from this object: an edit made to it in memory after
        the run is not in the template. Observation targets read the run's
        output files at build time.
    name : str
        Project name (``str()`` is applied). It names the default template
        directory (each run of characters outside ``A-Za-z0-9_.-`` becomes one
        ``_``, leading and trailing ``_`` are stripped), the
        default ``.pst`` stem, the ``<name>.pst`` that :meth:`run_ies` and
        :meth:`prior` always write and run, their master directories
        (``<name>_ies_master`` / ``<name>_prior_master``, each with its workers in
        ``<master>.workers``), and the run as
        ``model.pest_runs`` lists it.
    workspace : str or pathlib.Path or None, optional
        The PEST template directory :meth:`build` writes. Default ``None`` uses
        ``<model workspace>.pest/<slug>`` (the slug as ``name`` describes), a
        sibling of the model directory where ``model.pest_runs`` finds it.
        **What a build deletes** (ledger 360): this directory and everything in
        it, and the staging directory it makes beside it, ``<workspace>.staging``
        (removed again when the build ends) -- each ONLY when it is absent, empty,
        or a directory a previous build wrote there in that role, which the build
        marks with a ``myflopy_pest_build.json`` naming the directory -- a visible
        file, so emptying the folder removes it too (a template built before that
        file existed counts when it IS the default template directory and its
        ``myflopy_pest_metadata.json`` names this ``name`` and this run directory;
        ledger 361). So re-running
        :meth:`build` on the same ``workspace`` rebuilds it, a serial
        :meth:`run_ies`/:meth:`prior`'s output in it included, while an existing
        folder a build did not write -- a Project root, a data folder, a copy of a
        template -- is refused (``FileExistsError``, naming it and what it holds)
        and left exactly as it was, and a template a serial run is still going in is
        refused (``RuntimeError``, ledger 532). Everything is checked before anything is
        deleted. A ``workspace`` (or its staging directory) that IS, CONTAINS or
        sits INSIDE the model's run directory is refused here and by
        :meth:`build` (``ValueError``) whatever it holds, because the build would
        delete the run or copy it into itself (ledger 107). Both refusals name the
        default and "a new or empty folder" as the spelling to use instead.
    start_datetime : str or pandas.Timestamp
        Simulation start date (e.g. ``"2020-01-01"``), passed to
        ``PstFrom(start_datetime=...)``, which parses it to a datetime to place
        time-varying parameters and observations on the time axis. Required
        here, keyword-only; ``model.pest(name)`` fills it from the model's
        TDIS ``start_date_time``.
    spatial_reference : object, optional
        Cell spatial reference passed to ``PstFrom(spatial_reference=...)``,
        which ``grid``-style parameters and ``correlation=`` geostatistics
        locate cells with. Default ``None`` uses, at :meth:`build`, the
        ``modelgrid`` of the flow model the build loads from its staging copy of
        the run -- ``model.gwf.modelgrid``'s mesh with the coordinates the run's
        files hold (FloPy writes them to 8 decimals), as the template's model
        reads them.
    zero_based : bool, optional
        Passed to ``PstFrom(zero_based=...)``: whether the index columns of the
        model's list files count from zero. Default ``False``, because MF6
        writes its ``layer cell`` columns one-based (pyEMU's own default is
        ``True``).
    longnames : bool, optional
        Passed to ``PstFrom(longnames=...)``: allow parameter and observation
        names longer than classic PEST's limits. Default ``True``.
    """

    def __init__(
        self,
        model,
        name: str,
        *,
        workspace: str | Path | None = None,
        start_datetime: str,
        spatial_reference=None,
        zero_based: bool = False,
        longnames: bool = True,
    ):
        """Create a calibration project around an existing model workspace.

        ``workspace`` is the PEST template directory. When omitted it defaults to
        ``<model workspace>.pest/<name>`` -- a SIBLING of the model directory, so
        the calibration lives beside the model it calibrates, is
        auto-discoverable via ``model.pest_runs`` / ``run.pest_runs``, and is
        not inside the tree the build copies (ledger 107). A ``workspace`` that
        is, contains or sits inside the model's run directory is refused here
        (``ValueError``, ledger 360). See the class docstring for the full
        reference, including what :meth:`build` deletes.
        """

        self.model = model
        self.name = str(name)
        self.original_workspace = Path(model.workspace)
        if workspace is None:
            # A SIBLING of the model directory, not a child: `PstFrom` copies the
            # model workspace wholesale via a bare `shutil.copytree` in a private
            # function, with no way to exclude a subtree, so a template inside it
            # gets copied into itself (ledger 107).
            workspace = default_pest_root(self.original_workspace) / _pest_run_slug(name)
        self.template_workspace = Path(workspace)
        # Refused here, the earliest it can be, and again by `build()` (ledger 360).
        self._refuse_a_workspace_over_the_run()
        self.start_datetime = start_datetime
        self.spatial_reference = spatial_reference
        self.zero_based = bool(zero_based)
        self.longnames = bool(longnames)
        self._observation_specs: list = []
        self._forecast_specs: list = []
        self._native_parameter_specs: list[NativeParameterSpec] = []
        self._native_parameter_frames: dict[str, object] = {}
        # One entry per pilot-point SET (call, layer file, zone), in build order
        # (ledger 525): its manifest record, parameter frame, template
        # and the spec that made it.
        self._pilot_point_sets: list[dict] = []
        self._capture_field_specs: dict[str, dict] = {}
        self._prepared_observations: list[dict] = []
        # The RIV files whose forward run carries rbot with a `riv.stage` offset, in
        # the build under way: emptied at each build, which makes a new template and
        # forward run -- kept across builds, a rebuild skipped every file and the
        # second template moved the stage alone (`_carry_river_bottom`).
        self._carried_river_files: set[str] = set()
        # The summed riv.stage lower bound already compiled onto each RIV file, in the
        # build under way (ledger 568).
        self._river_offset_lower_bounds: dict[str, float] = {}
        # The simulation LOADED from the build's staging copy of the run, while
        # `build()` runs; every recipe that reads a package's array reads it here,
        # never off `self.model` (ledger 356).
        self._staged_simulation = None
        self.pyemu = None
        self.pf = None
        self.pst = None

    # -- declarative facade -----------------------------------------------
    #
    # These methods compile directly to native ``pyemu.utils.PstFrom`` calls
    # (see native_parameters.py): ``parameterize``/``observe``/``forecast``
    # declare the calibration, ``build`` writes the ``.pst`` and forward run.

    def parameterize(
        self,
        target: str,
        *,
        style: str | None = None,
        bounds: tuple[float, float] | None = None,
        physical: tuple[float, float] | None = None,
        transform: str = "log",
        additive: bool | None = None,
        zones=None,
        layers=None,
        package: str | Sequence[str] | None = None,
        features: int | str | Sequence[int | str] | None = None,
        pp_spacing: float | None = None,
        pp_points=None,
        pp_nearest: int = DEFAULT_PP_NEAREST,
        correlation: float | None = None,
        anisotropy: float = 1.0,
        bearing: float = 0.0,
        nugget: float = 0.0,
        periods: str | Mapping[str, Sequence[int]] | None = None,
        name: str | None = None,
        capture: bool = False,
        **extra,
    ) -> NativeParameterSpec:
        """Declare a calibration parameter that compiles to ``pf.add_parameters``.

        Each call expands, at :meth:`build` time, into the correct
        ``pyemu.utils.PstFrom.add_parameters`` invocation against the external
        MODFLOW input file(s) for ``target`` -- except array-target pilot points,
        which myflopy builds itself (:mod:`~myflopy.modflow.mf6.pest.pilot_points`).
        Call it once per property you want to adjust; call it multiple times for
        the same target to stack styles (e.g. a broad ``"constant"`` plus a fine
        ``"grid"``, or a ``"constant"`` plus ``"pilotpoints"``), each call under
        its own ``name=`` (a repeated name is refused, ledger 523): every
        multiplier on a cell multiplies, pilot points included -- they are
        applied after pyEMU's multipliers, onto the file pyEMU wrote, so a
        constant x2 under pilot points at 1.0 writes 2x the base (ledger 525).

        Parameters
        ----------
        target : str
            Friendly name of the model property to calibrate (case and
            surrounding spaces ignored). Supported targets
            (with aliases): ``"k"`` (``"npf.k"``/``"hk"``/``"kh"``),
            ``"k33"`` (``"kv"``), ``"sto.ss"`` (``"ss"``/``"specific_storage"``),
            ``"sto.sy"`` (``"sy"``/``"specific_yield"``), ``"recharge"``
            (``"rch"``), ``"chd"``, ``"ghb.cond"``/``"ghb.bhead"`` (``"ghb"``),
            ``"drn.cond"``/``"drn.elev"`` (``"drn"``),
            ``"riv.cond"`` (``"riv"``)/``"riv.stage"``, ``"wel"`` (``"pumping"``),
            ``"uzf.vks"`` (``"vks"``) — UZF saturated vertical K, the one UZF
            quantity heads respond to —, ``"lak.bedleak"`` (``"bedleak"``) — the
            lake-bed leakance of each LAK connection —, ``"sfr.rhk"`` (``"rhk"``) —
            each reach's streambed K — and ``"porosity"``
            (``"mst.porosity"``/``"n"``), which lives on the GWT sibling of a
            coupled transport simulation and is resolved there automatically.
            ``k``, ``k33``, ``sto.ss``, ``sto.sy`` and ``porosity`` are ARRAY
            targets; ``recharge`` is an ARRAY target on an array RCH (READASARRAYS,
            what ``mf.rch.array`` writes) and a LIST target on a list RCH, its form
            read from the run's package; the rest scale a column of a LIST
            package's files. An unknown name raises ``KeyError`` here, listing the
            known ones.

            ``recharge`` on an array RCH scales the per-period ``recharge`` arrays
            -- one file per period BLOCK, which FloPy numbers by the block's
            ordinal, not its period, so each file's period is read from FloPy's
            records -- and never ``irch``; its period files share one parameter set
            (one per ``periods=`` group), and every style applies, pilot points
            included: one field, interpolated once per run, multiplies each period
            file (ledger 547). ``capture=True`` on a target of several
            period files raises ``NotImplementedError`` at :meth:`build`.

            ``sto.ss`` and ``sto.sy`` scale the STO griddata (under STO's
            ``STORAGECOEFFICIENT`` option ``ss`` holds a storage coefficient). Only
            transient data informs them: :meth:`build` refuses either on a run
            whose every stress period is steady state, and ``sto.sy`` where every
            cell of a parameter set -- the whole selection, a ``zones=`` label, or a
            layer for ``"grid"``/``"pilotpoints"`` -- has ``ICONVERT`` 0, since MF6
            reads SY only in a convertible cell (ledger 548).

            A target scales the package(s) of its TYPE -- ``riv.cond`` every RIV
            record, whatever the package is named -- found on the model the build
            stages from the run's files, and the external files FloPy's own records
            say those packages read (``<model>.<package file>_stress_period_data_<n>.txt``
            ...), never a file name guessed from the type (ledger 517).
            :meth:`build` raises ``ValueError`` when the model has no package of the
            type, or when the package is the array form (READASARRAYS) of a list
            target (``recharge`` on an array RCH). An ARRAY target whose package the
            run's files do not declare it on -- ``k33`` on an NPF written with no
            ``k33=``, where MF6 takes K33 from K -- is refused by :meth:`build` with
            a ``ValueError`` naming the target and the package.

            ``chd``, ``ghb.bhead``, ``drn.elev`` and ``riv.stage`` scale an
            ELEVATION (the package registry's ``FieldSpec.kind``) and are additive
            offsets in model length units (see ``additive``). ``riv.stage`` moves
            each record's river bottom by the same offset as its stage, so the
            declared depth survives every realization (ledgers 290/291's pair rule,
            ledger 519); :meth:`build` refuses a record whose stage is
            already below its bottom, which MF6 refuses at any icelltype, and a lower
            bound that would carry any record's bottom below its CELL's bottom --
            refused by MF6 at ``icelltype != 0``, over-draining at ``0`` (ledger
            290) -- naming how many, the worst record and the most negative lower
            bound that is safe; the lower bounds of several ``riv.stage``
            parameters on one file add (ledger 568).

            ``lak.bedleak`` and ``sfr.rhk`` are FEATURE targets (ledger
            558): their rows are lake connections and reaches, keyed
            by the lake and connection number (LAK) or the reach number (SFR) --
            never by the cell, which a lake connects to more than once and two
            reaches can share -- so ``features=`` and ``style="feature"`` act on
            them. ``lak.bedleak`` defaults to one multiplier per LAKE and takes
            ``"feature"`` or ``"constant"``; :meth:`build` refuses a ``NONE``
            bedleak on a lake it scales (MF6: no lakebed, nothing to multiply).
            ``sfr.rhk`` defaults to ``"constant"`` and also takes ``"grid"`` (one
            per reach), ``"feature"`` (one per named stream) and ``"zone"`` (one
            per ``zones=`` label of the reach's cell); :meth:`build` refuses a
            reach whose cellid is ``NONE``, whose row is one column short -- pyEMU,
            reading columns by position, multiplied its Manning's n. Under lake
            cells, K and the bed leakance trade off against each other: only
            data that see both -- a lake stage and heads under and beside it --
            separate them.

            Porosity is absent from the flow equation, so it is identifiable
            only from CONCENTRATION data — and since transport velocity is
            ``v = Ki/n``, K and porosity are near-collinear from concentration
            alone (measured cosine 0.98 on the canonical model). Estimate both
            only with heads in the mix: heads pin K, and concentration then pins
            porosity.
        style : str or None, optional
            Spatial parameterization style:

            - ``"constant"`` -- a single multiplier for the whole target;
            - ``"zone"`` -- one multiplier per ``zones`` value (pass ``zones=``);
            - ``"grid"`` -- one multiplier per list entry / cell;
            - ``"pilotpoints"`` -- pilot points interpolated to the cells (ARRAY
              targets only, ``recharge`` on an array RCH included);
            - ``"feature"`` -- one multiplier per lake (``lak.bedleak``) or per
              named stream, a reach BOUNDNAME (``sfr.rhk``), or per entry of
              ``features=`` -- FEATURE targets only.

            ``None`` uses the recipe default (``"feature"`` for ``lak.bedleak``,
            ``"constant"`` for every other target). A style the target does not
            take raises ``ValueError`` here, naming the ones it does --
            ``"feature"`` on a target with no features, ``"grid"``/``"zone"``/
            ``"pilotpoints"`` on ``lak.bedleak``, ``"pilotpoints"`` on ``sfr.rhk``
            (ledger 558).
            The value is lower-cased; anything but these five raises
            ``ValueError`` here -- pyEMU's own pilot-point spellings (``"pp"``,
            ``"pilot_points"``, ...), which reached pyEMU's structured-grid pilot
            points and died on a DISV grid, name ``"pilotpoints"`` (ledger
            520). It reaches pyEMU as ``par_type``, except
            ``"pilotpoints"`` on an array target, which myflopy builds itself;
            ``"pilotpoints"`` on a LIST target raises ``ValueError`` here (pyEMU's
            own pilot points reject a list file on every grid) -- use ``"grid"``,
            whose ``correlation=`` correlates the records; on ``recharge`` it raises
            when the RCH is the list form -- here when the model can list its
            packages, else at :meth:`build`.
            ``"grid"`` on
            *array* targets (K, K33) builds one multiplier per Voronoi/DISV
            cell, and an independent set per layer file, located through the
            model's spatial reference (the modelgrid is used automatically);
            pass ``correlation=`` to correlate them and ``layers=`` to restrict
            which layers. ``"grid"`` also works for list targets. A correlated
            ``"grid"`` block -- one layer file of an array target, or a list target
            -- of more than 5,000 parameters is refused by :meth:`build`: its
            prior is drawn from a dense covariance of every pair (3.3 GB for one
            20,323-cell layer); use ``"pilotpoints"`` (ledger 522).
            ``"pilotpoints"`` on an array target places pilot points
            (a ``pp_spacing=`` net or explicit ``pp_points=``; one set per layer
            file and per ``zones=`` label) and interpolates the multipliers to
            each cell from the ``pp_nearest`` nearest points of its set by inverse
            distance squared -- in log10 space for a log-transformed parameter --
            with weights computed once at build (ledger 526). It is the
            supported sparse-field style on Voronoi/DISV grids.
        bounds : tuple of (float, float) or None, optional
            ``(lower, upper)`` bounds on the adjustable value, pyEMU's
            ``lower_bound``/``upper_bound`` (each pilot point's ``parlbnd`` /
            ``parubnd``). For a multiplier these are factors (e.g. ``(0.2, 5)`` =
            5x down to 5x up); for an additive parameter they are offsets in model
            length units. ``None`` (default) is ``(0.5, 2.0)`` for a multiplier and
            ``(-2.0, 2.0)`` for an offset (ledger 518). Raises
            ``ValueError`` here when ``lower >= upper``, when either is not finite,
            when they EXCLUDE the value the parameter starts at -- ``1.0`` for a
            multiplier (every pilot point included), ``0.0`` for an offset, where
            PEST would begin outside its own bounds and the input value would never
            be tried -- or when a log-transformed multiplier's lower bound is not
            above 0.
        physical : tuple of (float, float) or None, optional
            ``(ult_lbound, ult_ubound)`` -- hard limits clamped onto the *final*
            model value after all multipliers are applied (pyEMU's
            ``ult_lbound``/``ult_ubound``; the pilot-point step, which multiplies
            after pyEMU has clamped its own product, clamps each cell it covers
            to the tightest limits of the pilot-point sets covering it AND of the
            pyEMU multipliers on its file, so a constant's limits hold under pilot
            points, and refuses at :meth:`build` covering limits that share no
            value -- ledger 528).
            Strongly recommended
            for multipliers so calibration cannot push K (etc.) to nonphysical
            values. ``None`` (default) leaves the input unclamped.
        transform : str, optional
            ``"log"`` (default, recommended for strictly positive properties
            like K and recharge) or ``"none"``; pyEMU also accepts
            ``"fixed"``, but not together with ``correlation=`` (the variogram
            conversion asserts ``"none"`` or ``"log"``). Forced to ``"none"``
            for additive parameters. Array-target pilot points interpolate in
            log10 space when it is ``"log"`` and linearly otherwise.
        additive : bool or None, optional
            Apply the parameter as an additive offset instead of a multiplier
            (pyEMU ``par_style="a"``). ``None`` (default) follows what the scaled
            field IS (the package registry's ``FieldSpec.kind``): an ELEVATION --
            ``chd``, ``ghb.bhead``, ``drn.elev``, ``riv.stage`` -- is an offset,
            everything else a multiplier. ``False`` on an elevation raises
            ``ValueError``: a multiplier on an elevation scales with the datum
            (ledger 518). ``True`` on array-target pilot points raises
            ``ValueError``: they MULTIPLY the array in their forward run.
        zones : ZoneSpec or array_like of int or None, optional
            Per-cell integer zone labels, normalized to the shape pyEMU wants
            for the target's family at :meth:`build` (see :class:`ZoneSpec`,
            which also reads a raster or polygon layer) and passed as pyEMU's
            ``zone_array``. With ``style="zone"`` each distinct label gets one
            parameter. With any style, an array target leaves cells labelled
            ``<= 0`` unparameterized; a list target makes zone ``0`` a real
            parameter and warns. A length that is not one label per cell raises
            ``ValueError`` at build; the labels are shaped on the grid of the model
            the build stages from the run's files (ledger 356). With
            ``style="constant"`` on a LIST target -- pyEMU reads ``zones=`` only for a
            ``"zone"`` or ``"grid"`` list parameter, so one number scaled every record
            -- it raises ``ValueError`` naming ``style="zone"``: here, or for
            ``recharge`` once its RCH's form is known (here when the model lists its
            packages, else at :meth:`build`; ledger 569).
            With array-target ``"pilotpoints"`` each label ``> 0`` gets its own
            pilot-point set on each selected layer -- its own net (or the
            explicit points inside it), interpolated onto its own cells only, so
            nothing blends across a zone edge, with its own prior covariance --
            while cells labelled ``<= 0`` keep their value (ledger
            527). An explicit point in a cell labelled ``<= 0``, a
            label ``> 0`` no explicit point lies in, and labels that are all
            ``<= 0`` raise ``ValueError`` at build. Different settings per zone
            are separate calls, each with a masked ``zones=`` and its own
            ``name=``. Default ``None`` passes no zones.
        layers : int or iterable of int or None, optional
            For multi-layer *array* targets, restrict the parameter to these
            zero-based model layers (e.g. ``layers=[0, 1]`` for the unconfined
            aquifers; a bare ``2`` is one layer, and numpy integers are taken as
            ``int``); other layers keep their input value. ``None`` (default)
            covers all layers. A value that is not a whole number at or above 0,
            or an empty list, raises ``ValueError`` here. Given for a LIST target it raises ``ValueError``
            here, naming ``package=``: a list target scales every record of its
            package's files, and ``zones=`` leaves none of them out (with
            ``style="zone"`` every label, ``0`` and negatives too, is a parameter of
            its own; ledger 559). ``recharge`` in either form, which MF6 applies to
            one cell per column, raises the same, naming ``zones=``. :meth:`build` raises ``ValueError`` when no requested
            layer has its own file, or the target is written as one whole-grid
            file.
        package : str or sequence of str or None, optional
            Which package(s) of the target's TYPE the parameter scales, by NAME
            (FloPy's ``package_name``, any case) -- ``package="crisp"`` on a model
            whose RIV packages are ``riv`` and ``crisp``. A sequence,
            ``package=("riv", "crisp")``, shares ONE parameter set across the
            packages it names. ``None`` (default) is the model's only package of
            the type; with several it raises ``SeveralPackagesError`` naming each
            and these spellings -- a package is never picked in silence (ledger
            352) -- here when the model's packages can be listed, else at
            :meth:`build`. A sequence whose packages' records differ in width -- one
            declares BOUNDNAMES or AUXILIARY variables the other does not -- raises
            ``ValueError`` at :meth:`build`, naming each package and the call per
            package that builds (pyEMU reads one set's files as one table; ledger
            559). A name that is no package of the type raises
            ``ValueError`` naming the ones there are; a name given twice, or an
            empty one, raises here. The default ``name`` gains ``_<package>``.
            A sequence with ``features=`` or ``style="feature"`` raises
            ``ValueError`` here: lakes and reaches are numbered and named within
            ONE package.
        features : int or str or sequence of int or str or None, optional
            The lakes (``lak.bedleak``) or reaches and named streams
            (``sfr.rhk``) the parameter covers -- each MF6's number for it,
            ZERO-based (``LakeStageTargets``' convention; text reading as a whole
            number is that number), or its boundname, any case; one may be given
            bare. With ``style="feature"`` each entry is ONE parameter, named
            ``<name>_<label>`` -- the entry's boundname (a lake given by number
            takes its boundname when it has one) or number, made a PEST name by
            the observation groups' rule (lower-cased, every run of characters but
            letters, digits, ``_``, ``.``, ``:``, ``-`` one ``_``); with
            ``"constant"`` the entries share one; with ``"grid"`` each reach among
            them is one. ``None`` (default): every lake / every reach, and for
            ``"feature"`` one parameter per lake (``lak.bedleak``) or per
            boundname (``sfr.rhk``; reaches without one are left alone, and an SFR
            with no BOUNDNAMES is refused, naming ``style="grid"``). Raises
            ``ValueError`` here on a target with no features, with
            ``style="zone"``, for an entry that is blank, a fraction, negative or
            a bool, and -- here when the model lists its packages, else at
            :meth:`build` -- for an entry naming no lake/reach/boundname of the
            package (listing those it has), a lake or reach named twice, or two
            labels making one PEST name (ledger 558).
        pp_spacing : float or None, optional
            Pilot-point net spacing for ``style="pilotpoints"`` on an array
            target, in MODEL LENGTH UNITS (feet on a model in feet), finite and
            above zero: the grid's cell centres are binned in squares this wide
            from the lowest centre's ``x`` and ``y``, and the cell centre nearest
            each occupied square's middle becomes a point -- per zone with
            ``zones=``. :meth:`build` records each set's point count on the
            returned spec (``pilot_point_sets``) and in :meth:`settings`.
            Default ``None``. Array-target pilot points take exactly one of it
            and ``pp_points`` -- neither, or both, raises ``ValueError`` here --
            and any other style raises ``ValueError`` here: it cannot act
            (ledger 520's rule). It replaced ``pp_space`` (multiples
            of the median cell width), which laid 4,501 points per layer on a
            refined 20,323-cell mesh (ledger 527).
        pp_points : geopandas.GeoDataFrame or array_like or None, optional
            Explicit pilot-point locations for ``style="pilotpoints"`` on an
            array target: a GeoDataFrame of Point geometries, or ``(x, y)``
            pairs (anything reshapable to ``(n, 2)``), in the model grid's own
            coordinates -- a GeoDataFrame's CRS is not read or reprojected. Each
            point takes the zone of the cell it lies in; a point in no cell, or
            given twice, raises ``ValueError`` at :meth:`build`. Default ``None``.
            Beside ``pp_spacing``, or with any other style, it raises
            ``ValueError`` here.
        pp_nearest : int, optional
            How many of its set's nearest pilot points each cell is interpolated
            from, by inverse distance squared (all of them when the set has
            fewer); a whole number, 1 or more. Default ``8`` (ledger
            526). A cell whose centre is a point takes that point's
            value. A value other than ``8`` raises ``ValueError`` here with any
            style but array-target ``"pilotpoints"``.
        correlation : float or None, optional
            Variogram range (model length units) of an exponential variogram
            (see :class:`~myflopy.modflow.mf6.pest.specs.ExpGeoStruct`). With
            ``style="grid"`` it becomes the parameters' pyEMU ``geostruct``, and
            :meth:`run_ies`/:meth:`prior` draw a spatially correlated prior
            ensemble from it (``ies_par_en``). With ``style="pilotpoints"`` on an
            array target it shapes only that prior draw -- the pilot values of
            each set (one layer, one zone) are drawn correlated with each other
            and not with another set's, while the interpolation to cells stays
            inverse-distance. ``None`` (default) builds no variogram. Given with
            ``"constant"``, ``"zone"`` or ``"feature"`` it raises ``ValueError``
            here: they have nothing to correlate (a feature parameter is keyed by a
            lake's or stream's number, not located at a cell). On a LIST target's
            ``"grid"`` pyEMU places each record at its cell's centre, and
            :meth:`build` refuses it (``ValueError``) where it cannot correlate them:
            two parameters at one point in one covariance group -- a column's
            records in two layers on DISV or DISU, where pyEMU keeps one group, or
            two reaches in one DIS cell and layer -- make the covariance singular
            (records sharing a whole cellid share ONE parameter, pyEMU's rule), and pyEMU's draw at :meth:`run_ies`/:meth:`prior` died with
            ``Matrix is not positive definite``; and where pyEMU, finding a shared
            point, splits the set into one group per record -- every reach of an
            ``sfr.rhk`` set on DISV -- the prior is drawn uncorrelated (both
            measured; ledger 558, ledger 570). On DIS a column's records
            in different layers are let through: pyEMU correlates each layer's records
            and not across layers. On ``sfr.rhk`` on DISU, where pyEMU cannot locate a
            reach, it is refused too (ledger 558).
        anisotropy : float, optional
            Anisotropy ratio of the variogram -- correlated this many
            times farther along the major axis than across it (``1.0`` =
            isotropic, the default). For an alluvial valley, e.g. ``5`` makes K
            correlate farther down-valley than across it. Acts only when
            ``correlation`` is given, and anything but ``1.0`` without it raises
            ``ValueError`` here; pyEMU asserts it is positive.
        bearing : float, optional
            Azimuth of the anisotropy major axis, in degrees east of north.
            Default ``0.0``. Anything else with ``anisotropy=1.0``, or without
            ``correlation``, raises ``ValueError`` here: it cannot act.
        nugget : float, optional
            Nugget (unresolved short-scale variance) of the variogram, relative
            to the variogram's contribution, which is fixed at ``1.0`` here: a
            nugget ``n`` is the share ``n / (1 + n)`` of each parameter's prior
            variance that is uncorrelated even with its nearest neighbour. The
            variance itself stays the one the bounds imply (log10 standard
            deviation ``log10(ub / lb) / 6`` for a log parameter); what moves is
            the correlation of two parameters ``h`` apart in the prior
            :meth:`run_ies`/:meth:`prior`/:meth:`draw_prior` draw,
            ``exp(-h / correlation) / (1 + n)`` -- so ``1.0`` halves every
            correlation (ledger 489; before, the draw discarded the nugget).
            The same rule for ``style="grid"`` and pilot points. Default ``0.0``.
            A negative or non-finite value raises ``ValueError`` here. Acts only
            when ``correlation`` is given; anything but ``0.0`` without it raises
            ``ValueError`` here.
        periods : str or mapping of str to sequence of int or None, optional
            Groups the target's PER-PERIOD files -- ``recharge`` in either form and
            the list targets -- into parameter sets, one ``add_parameters`` call per
            group, named ``<name>_<label>`` (its own PEST group), with the style's
            parameters within each (ledger 549):

            - ``None`` (default): one set over every period, the steady-state period
              included.
            - ``"month"``: one group per calendar month of the TRANSIENT periods,
              ``m01``..``m12`` (January = ``m01``), a period's month being its
              MIDPOINT's (ledger 529) on the run's calendar -- TDIS lengths and
              ``TIME_UNITS``, ``START_DATE_TIME`` (else this calibration's
              ``start_datetime``) -- so one month in every cycle shares one
              parameter. A STEADY-STATE period (the flow model's STO; a model with
              no STO is steady throughout) is in no group: only a parameter
              declared without ``periods`` reaches it -- stack a plain
              ``parameterize("recharge", name=...)`` for that (the user's rule, Q3).
            - ``{label: [zero-based periods]}``: explicit groups; each label is
              letters, digits and ``_`` (lower-cased); a period in no group is left
              alone.

            MF6 repeats a period's block until the next, so a file serves every
            period up to its package's next block; a file serving periods of two
            groups -- or of a group and none, a steady-state period under
            ``"month"`` included -- is refused at :meth:`build`, naming both
            periods and the fix (a block for every period, as
            ``mf.AverageYear(...).by_month(...)`` gives). Raises ``ValueError``
            here on a target whose files are not per period (``k``, ``sto.sy``,
            ``uzf.vks``, ``porosity``), for a string but ``"month"``, a label that
            is not letters/digits/``_`` or repeats, an empty group, a period that
            is not a whole number ``>= 0`` or is named twice; and at :meth:`build`
            for a period the run does not have, a group no file serves, ``"month"``
            on a transient period longer than 31 days, a TDIS with no
            ``TIME_UNITS``, or a TDIS ``START_DATE_TIME`` that differs from
            ``start_datetime``. A ``TypeError`` for anything but a string, a
            mapping or ``None``. The monthly multipliers are drawn independently in
            the prior; a temporal correlation between them is deferred. It
            replaced ``temporal=``, which nothing read.
        name : str or None, optional
            Parameter group / name base (pyEMU's ``par_name_base`` and
            ``pargp``). Default ``None`` uses the canonical target with dots
            removed (``"k"``, ``"ghbcond"``, ``"uzfvks"``), plus ``_<package>``
            when ``package=`` is given (``"rivcond_crisp"``,
            ``"rivcond_riv_crisp"`` for a sequence; ``"lakbedleak"``,
            ``"sfrrhk"``). ``style="feature"`` makes each lake's or stream's
            parameter and PEST group ``<name>_<label>`` (``lakbedleak_keevie``).
            A multi-layer array ``"grid"``
            parameter appends ``l<layer>`` per selected layer file, the
            ZERO-BASED layer (ledger 521); pilot points append ``l<layer>``
            (zero-based) and, with ``zones=``, ``z<zone>`` -- a set
            ``<name>l<layer>z<zone>``, its points ``...pp0000``, all lower-cased as
            PEST reads them; and a captured field is named ``<name>field``. A name
            another ``parameterize`` call of this project already holds (any case)
            raises ``ValueError`` here (ledger 523): stacking calls on one target,
            or pilot points declared layer by layer, needs a distinct name for
            each. With ``periods=`` each group's parameters and PEST group are
            ``<name>_<label>`` (``rch_m03``), and a pilot-point set
            ``<name>_<label>[l<layer>][z<zone>]``; a pilot-point set on period files
            carries no ``l<layer>``.
        capture : bool, optional
            Also record the *resolved* model-input field (after multipliers are
            applied) as zero-weight observations named from ``<name>field``
            (one prefix per layer file for a layered array), so each
            realization's per-cell property field is carried in the ensembles.
            Enables spatial parameter maps via ``IesResults.plot_field`` /
            ``IesResults.field``. Default ``False``. It costs one observation per
            CELL of every layer captured, in every ensemble every iteration writes:
            ``k33`` over Ten Trails' 5 x 20,323 cells added 101,615 and grew the
            template from 18 to 33 MB, and one 5-layer field at 200 realizations is
            a 352 MB observation-ensemble CSV per iteration. A capture covers the
            layers the parameter does (``layers=``); ``run_ies(ies_save_binary=True)`` writes the
            ensembles in binary, which the review reads too. ``True`` on
            ``lak.bedleak``/``sfr.rhk`` raises ``ValueError`` here: their rows are
            lake connections and reaches, not one value per cell (deferred, ledger
            558).
        **extra
            Advanced keyword arguments passed straight through to pyEMU's
            ``PstFrom.add_parameters`` for what the facade does not set
            (``mfile_fmt``, ``initial_value``, ``apply_order``, ...). A keyword
            the facade sets itself raises ``TypeError`` here, naming the argument
            that sets it: ``par_type`` (``style=``), ``par_name_base`` and ``pargp``
            (``name=``), ``lower_bound``/``upper_bound`` (``bounds=``),
            ``transform`` (``transform=``), ``ult_lbound``/``ult_ubound``
            (``physical=``), ``par_style`` (``additive=``), ``zone_array``
            (``zones=``), ``geostruct`` (``correlation=``), and ``index_cols``,
            ``use_cols``, ``use_rows``, which the target, ``features=`` and
            ``package=`` decide (ledger 571; they overrode the
            facade's own, and an ``index_cols`` or ``use_rows`` failed only after
            the staging). A keyword ``add_parameters`` does not take raises
            ``TypeError`` too. Given for array-target pilot points, which never call
            ``add_parameters``, they raise ``ValueError`` here. Stored on the
            returned spec as ``extra``.

        Returns
        -------
        NativeParameterSpec
            The recorded specification (also stored on the project). Inspect the
            full set at any time with :meth:`settings`.

        Raises
        ------
        KeyError
            For an unknown ``target``.
        ValueError
            For an argument that cannot act or bounds that cannot hold (above), a
            ``name`` already declared, or a ``package=`` naming no package of the
            type.
        TypeError
            For an ``**extra`` keyword the facade sets itself or
            ``PstFrom.add_parameters`` does not take (above).
        SeveralPackagesError
            When the model carries several packages of the target's type and no
            ``package=`` names one.

        Examples
        --------
        >>> cal.parameterize("k", style="constant", bounds=(0.2, 5), physical=(1e-3, 100))
        >>> cal.parameterize("recharge", bounds=(0.5, 1.5), physical=(0, 1e-2))
        >>> cal.parameterize("ghb.cond", bounds=(0.1, 10))
        >>> cal.parameterize("drn.elev", bounds=(-2, 2))   # additive offsets, in ft
        >>> cal.parameterize("riv.cond", package="crisp", bounds=(0.1, 10))
        >>> cal.parameterize("riv.stage", package=("riv", "crisp"))  # one shared offset
        >>> cal.parameterize("k", style="pilotpoints", name="kpp", pp_spacing=1500.0,
        ...                  zones=units, physical=(1e-3, 300))  # stacks on the constant
        >>> cal.parameterize("recharge", name="rch", bounds=(0.5, 2.0))  # every period
        >>> cal.parameterize("recharge", name="rchm", periods="month",   # x one per month
        ...                  bounds=(0.5, 2.0))
        >>> cal.parameterize("sto.sy", layers=[0], bounds=(0.5, 2.0), physical=(0.01, 0.35))
        >>> cal.parameterize("bedleak", bounds=(0.1, 10))     # one multiplier per lake
        >>> cal.parameterize("bedleak", features=["keevie"], name="kv")  # one lake
        >>> cal.parameterize("rhk", style="feature", bounds=(0.2, 5))   # per named stream
        """

        nugget = checked_nugget(nugget)
        spec = NativeParameterSpec(
            target=target,
            style=style,
            bounds=bounds,
            physical=physical,
            transform=transform,
            additive=additive,
            zones=zones,
            layers=layers,
            package=package,
            features=features,
            pp_spacing=pp_spacing,
            pp_points=pp_points,
            pp_nearest=pp_nearest,
            correlation=correlation,
            anisotropy=anisotropy,
            bearing=bearing,
            nugget=nugget,
            periods=periods,
            name=name,
            capture=capture,
            extra=extra,
        )
        self._refuse_a_repeated_parameter_name(spec)
        self._refuse_an_unnamed_package(spec)
        self._native_parameter_specs.append(spec)
        return spec

    def _refuse_a_repeated_parameter_name(self, spec: NativeParameterSpec) -> None:
        """Refuse a parameter name another :meth:`parameterize` call already holds (any case).

        The name keys pyEMU's parameter names and group, the project's frames and the
        build's metadata, so a second call under it overwrote the first's records:
        two pilot-point calls on one target left the first set's parameters out of
        the control file while its forward-run step kept reading a table frozen at
        1.0 (ledger 523).
        """

        wanted = str(spec.name).strip().lower()
        for earlier in self._native_parameter_specs:
            if str(earlier.name).strip().lower() == wanted:
                raise ValueError(
                    f"Parameter name {spec.name!r} is already declared, by "
                    f"parameterize({earlier.recipe.canonical!r}, style={earlier.style!r}). "
                    "Each call needs its own: pass name= (e.g. "
                    f"name={(spec.name + '_' + spec.style)!r})."
                )

    def _refuse_an_unnamed_package(self, spec: NativeParameterSpec) -> None:
        """Refuse, at the call, a ``package=`` the model's packages already show to be wrong.

        :meth:`build` decides on the staged model (:func:`select_packages`); this asks
        the model this project calibrates the same question, so several packages of the
        type and no ``package=`` -- or a name that is none of them -- is refused when it
        is written. A transport target, or a model that cannot list its packages, is
        left to :meth:`build`.
        """

        if spec.recipe.model != "flow":
            return
        try:
            held = model_packages_of_type(self.model, spec.recipe.package)
        except (AttributeError, TypeError, KeyError, OSError, FlopyException,
                MFDataException) as error:
            # A model that cannot list its packages (a test double, a run whose files
            # FloPy cannot read): the build's check on the staged model still runs.
            logger.debug("Package check of %r left to build(): %s", spec.name, error)
            return
        if held:
            # The package's FORM too: pilot points on a LIST RCH are refused here.
            packages = select_packages(self.model, spec.declared_recipe, spec.package)
            spec.resolve_form(packages)
            if spec.declared_recipe.feature is not None and len(packages) == 1:
                # A lake or stream features= names, checked against the package's
                # packagedata now; the build asks the staged one again (ledger
                # 558).
                resolve_features(spec, packages[0])

    def observe(self, targets, *, prefix: str | None = None):
        """Register history-matching observations from a ``myflopy`` target set.

        Observations are the measured data calibration tries to reproduce. At
        :meth:`build` time the targets are registered with pyEMU
        (``pf.add_observations``), their measured values and weights are written
        into the control file, and a post-processor is added to the forward run
        so simulated equivalents are regenerated on every model run.

        What the control file holds is exactly the targets (ledger 340). Each is
        matched to the simulated row whose index EQUALS its ``time`` read as a
        zero-based stress period (``1.0`` is period 1, ``1.5`` is refused), and a
        target no row holds -- a period the output did not save, a name with no
        simulated column, a head or concentration target on no active cell
        (ledger 341) -- makes :meth:`build` raise, naming it. A blank value
        is unmeasured. Every simulated row no target measures stays in the
        control file at weight 0, its ``obsval`` the base run's own value, so it
        is reported but never history-matched; it used to carry weight 1, which
        pulled the calibration toward the uncalibrated model at every period
        nobody measured. A set given no values at all registers every row at
        weight 0: its simulated values are captured, nothing is fitted.

        A FLOW is written to the control file in MF6's own sign, on both routes
        (ledger 350): an SFR reach's streamflow and a DRN zone's discharge are
        negative, water leaving the reach or the aquifer, so a gauged ``Q`` is
        entered as ``-Q`` -- the targets refuse a value above zero when they are
        built, which scored ``-2Q`` at a perfect fit.

        **A target's group is its PEST observation group** (ledger 358): each
        measured observation is in the group its target's ``group`` label names
        (a head, concentration or DRN-zone target's ``group_column``), made a PEST
        name -- lower-cased, every run of characters other than letters, digits,
        ``_``, ``.``, ``:`` and ``-`` one ``_`` -- and one with no label is in the
        set's prefix (``hds``, ``lak_stage``, ...), so PEST++ reports phi by your
        groups, and by data type out of the box. Every row no target measures is
        in ``<prefix>_unmeasured``, so it never joins a group of yours. A label
        that is not text, holds no letter or digit, starts as PESTPP-IES reads an
        inequality (``l_``, ``less``, ``g_``, ``greater``), takes a group the
        build gives a set, or becomes the same name as a different label makes
        :meth:`build` raise, naming it. pyEMU used to give each simulated column
        a group of its own, so a group was one well.

        **A set given ``inequality=`` is a one-sided BOUND** (ledger
        564): its measured observations are in PESTPP-IES's
        inequality group -- ``l_`` (``"at_most"``) or ``g_`` (``"at_least"``)
        before the group they would otherwise be in -- which scores 0 wherever the
        bound holds and the full weighted misfit where it does not. That is the
        ONE spelling of an inequality: a label or prefix starting ``l_``/``less``/
        ``g_``/``greater`` stays refused.

        Parameters
        ----------
        targets : target set or observation spec
            A high-level target set -- :class:`HeadTargets`,
            :class:`ConcTargets`, :class:`LakeStageTargets`,
            :class:`SfrStageTargets`, :class:`SfrFlowTargets`,
            :class:`DrnFlowTargets` or :class:`RivFlowTargets` -- which is
            wrapped in its matching ``*ObservationSpec``, or a pre-built
            observation spec (e.g. :class:`HeadTargetObservationSpec`),
            registered as it is. A set
            registered on a model and read back as ``model.targets.<name>`` (a
            bound helper) is unwrapped to its target set first. The
            matching simulated values are regenerated by a post-processor on
            every run. Anything else raises ``TypeError``.
        prefix : str or None, optional
            Observation-name prefix. ``None`` (default) uses the target kind's
            own prefix (``"hds"``, ``"conc"``, ``"lak_stage"``, ``"sfr_stage"``,
            ``"sfr_flow"``, ``"drn_flow"``, ``"riv_flow"``). Pass distinct prefixes to register
            more than one set of the same kind. It is also, made a PEST name,
            the observation group of the set's targets with no group label --
            so two sets under two prefixes are two groups -- and with
            ``_unmeasured`` added the group of its unmeasured rows; one starting
            ``l_``, ``less``, ``g_`` or ``greater`` (an inequality to
            PESTPP-IES) makes :meth:`build` raise (ledger 358). IGNORED when
            ``targets`` is already a spec: the spec's own ``prefix`` is used.

        Returns
        -------
        observation spec
            The registered spec -- the one passed in, or the
            ``*ObservationSpec`` (e.g. :class:`HeadTargetObservationSpec`)
            the target set was wrapped in.
        """

        spec = self._coerce_observation_spec(targets, prefix=prefix)
        self._observation_specs.append(spec)
        return spec

    def forecast(self, targets, *, prefix: str | None = None):
        """Register a prediction of interest as a zero-weight forecast.

        Forecasts use the *same* target objects as :meth:`observe`, but they are
        never history-matched: at :meth:`build` time their weight is set to zero
        and their observation names are recorded in
        ``pst.pestpp_options["forecasts"]``. The downstream uncertainty tools
        (IES ensembles, FOSM) report posterior uncertainty for exactly these
        quantities, so declare here whatever model output you ultimately care
        about predicting (a future head, a stream flux, a seepage rate).

        A forecast is the rows its targets NAME, matched as :meth:`observe`
        matches them (a target no row holds makes :meth:`build` raise); its value
        may be blank. The set's other rows -- the periods nobody forecast -- are
        zero-weighted too and are not forecasts; they used to stay at weight 1,
        history-matched to the base run's own values (ledger 340). A set given
        NO values forecasts every row of its table: its locations at every
        period.

        A forecast is in its set's prefix as its PEST observation group
        (``fore1``, ...), and the set's other rows in ``<prefix>_unmeasured``
        (ledger 358): a forecast never joins a group of measured observations,
        and its targets' ``group`` labels are not read.

        Parameters
        ----------
        targets : target set or observation spec
            Any target set :meth:`observe` accepts (:class:`HeadTargets`,
            :class:`ConcTargets`, :class:`LakeStageTargets`,
            :class:`SfrStageTargets`, :class:`SfrFlowTargets`,
            :class:`DrnFlowTargets`, :class:`RivFlowTargets`,
            ``model.targets.<name>`` included) or a
            pre-built observation spec,
            describing *where* and *when* the prediction is taken. Any values
            supplied are placeholders (e.g. an expected value), written as the
            observation's value; a blank keeps the base run's. The weight is
            forced to zero regardless. Anything else raises ``TypeError``, and a
            BOUND set -- one given ``inequality=`` -- raises ``ValueError``: a
            bound on a zero-weight quantity constrains nothing (ledger
            564).
        prefix : str or None, optional
            Observation-name prefix. ``None`` (default, and an empty string
            too) gives ``"fore1"``, ``"fore2"``, ... in call order, so forecast
            names never collide with calibration observations. IGNORED when
            ``targets`` is already a spec: the spec's own ``prefix`` is used,
            which for a default-built spec is the same as a calibration set's
            (``"hds"``, ...). It is the forecasts' observation group, refused as
            :meth:`observe`'s is.

        Returns
        -------
        observation spec
            The registered forecast spec -- the one passed in, or the
            ``*ObservationSpec`` the target set was wrapped in.
        """

        default = prefix or f"fore{len(self._forecast_specs) + 1}"
        spec = self._coerce_observation_spec(targets, prefix=default)
        inequality = getattr(spec.targets, "inequality", None)
        if inequality is not None:
            raise ValueError(
                f"forecast(): these {type(spec.targets).__name__} are a BOUND "
                f"(inequality={inequality!r}), and a forecast is a zero-weight quantity "
                "to predict -- a bound on it constrains nothing. Forecast the same "
                "locations from a set without inequality=, or observe() the bound to "
                "history-match it (ledger 564)."
            )
        self._forecast_specs.append(spec)
        return spec

    def _coerce_observation_spec(self, targets, *, prefix):
        """Turn a target set (or spec) into a registered observation spec.

        Accepts a pre-built observation spec, or any high-level target set
        (head / lake stage / SFR stage / SFR flow / DRN seepage), wrapping it in
        the matching spec -- a set registered on the model and read back as
        ``model.targets.<name>`` included, which is a bound helper around it and
        raised ``TypeError`` here. When ``prefix`` is ``None`` the spec's own default
        prefix is used (``"hds"``, ``"lak_stage"``, ``"sfr_flow"``, ...); pass distinct
        prefixes to register more than one set of the same kind.
        """

        targets = bare_targets(targets)
        if isinstance(targets, _OBSERVATION_SPEC_TYPES):
            return targets
        for target_type, spec_type in _TARGET_SPEC_TYPES:
            if isinstance(targets, target_type):
                if prefix is None:
                    return spec_type(targets=targets)
                return spec_type(targets=targets, prefix=prefix)
        raise TypeError(
            "observe()/forecast() accept a head/lake/SFR/DRN target set or a "
            f"pre-built observation spec; got {type(targets).__name__}."
        )

    def _default_template(self) -> Path:
        """The template ``model.pest(name)`` builds into: ``<model workspace>.pest/<slug>``."""

        return default_pest_root(self.original_workspace) / _pest_run_slug(self.name)

    def _where_to_build_instead(self) -> str:
        """The spellings a refusal names: the default (unless it IS the refused one) or a new folder."""

        default = self._default_template()
        lines = ["Build it somewhere else:"]
        if default.resolve() != self.template_workspace.resolve():
            lines.append(
                f"    model.pest({self.name!r}, ...)   # the default, beside the model: {default}"
            )
        lines.append(
            f"    model.pest({self.name!r}, workspace=<a new or empty folder>, ...)"
        )
        return "\n".join(lines)

    def _refuse_a_workspace_over_the_run(self) -> None:
        """Refuse a template or staging directory that is, contains or sits inside the run directory.

        A build deletes both directories and writes into them, so one that contains
        the run deletes it, and one inside it writes into the run, which a build only
        reads (ledger 356). While the build copied the WHOLE run directory (until
        ledger 383) one inside it also copied the run into itself --
        ``pest/<name>/pest/<name>/...``, 17 GB 40 levels deep (ledger 107). Refused
        whatever the directory holds, a build's marker included (ledger 360).
        Called when the project is made and again by :meth:`build`.

        Raises
        ------
        ValueError
            Naming the directory, the run directory and the spelling to use instead.
        """

        run = self.original_workspace.resolve()
        template = self.template_workspace.resolve()
        checks = [(template, "it")]
        # The staging directory only once the template passes: the file-system root has
        # no staging directory to name, and it contains every run.
        if _relation(template, run) is None:
            # Resolved itself too: `<template>.staging` may be a link to the run.
            staging = self._staging_workspace().resolve()
            checks.append((staging, f"its staging directory {staging}"))
        for directory, what in checks:
            relation = _relation(directory, run)
            if relation is None:
                continue
            raise ValueError(
                f"Refusing to build calibration {self.name!r} at workspace="
                f"{str(self.template_workspace)!r}: {what} {relation} the "
                f"model's run directory" + ("" if relation == "IS" else f" {run}")
                + ". A build deletes its workspace and the staging directory beside it, "
                "then writes into them, so the run would be deleted, or written into.\n\n"
                + self._where_to_build_instead()
            )

    def _clearance(self, directory: Path, role: str) -> str:
        """How a build may clear ``directory`` (resolved): ``"absent"``, ``"empty"`` or ``"built"``.

        ``"built"`` is a directory a previous build wrote in ``role`` -- its
        :data:`BUILD_MARKER` names it -- or, for the template, one built before the
        marker existed (:meth:`_a_template_built_before_the_marker`). Anything else
        is refused.

        Raises
        ------
        FileExistsError
            When ``directory`` is a file, or a non-empty directory no build wrote,
            naming it, what it holds and the spelling to use instead.
        """

        if not os.path.lexists(directory):
            return "absent"
        where = (
            f"Refusing to build calibration {self.name!r} at workspace="
            f"{str(self.template_workspace)!r}: "
            + ("it" if role == "template" else f"the staging directory beside it, {directory},")
        )
        if not directory.is_dir():
            raise FileExistsError(
                f"{where} is a FILE, and a build replaces it with a directory.\n\n"
                + self._where_to_build_instead()
            )
        entries = sorted(entry.name for entry in directory.iterdir())
        if not entries:
            return "empty"
        if _written_by_a_build(directory, role):
            return "built"
        if role == "template" and self._a_template_built_before_the_marker(directory):
            return "built"
        shown = ", ".join(entries[:5]) + (" ..." if len(entries) > 5 else "")
        raise FileExistsError(
            f"{where} already holds {len(entries)} entr{'y' if len(entries) == 1 else 'ies'} "
            f"({shown}), and no myflopy PEST build wrote it: it has no {BUILD_MARKER} "
            f"naming it. A build DELETES its {role} directory and everything in it, so "
            f"it deletes only one that is absent, empty, or a {role} directory a "
            "previous build wrote there.\n\n"
            + self._where_to_build_instead()
            + f"\nIf {directory} is an old {role} directory you no longer need, delete "
            "it yourself."
        )

    def _a_template_built_before_the_marker(self, directory: Path) -> bool:
        """A template a build wrote before :data:`BUILD_MARKER` existed: THIS calibration's, by its metadata.

        It has no marker (a marker, present, decides alone), it IS the directory
        ``model.pest(name)`` builds this calibration's template in, and its
        ``myflopy_pest_metadata.json`` records this project's name and this model's run
        directory -- so another calibration's template, and a COPY of this one anywhere
        else (an IES master, or a backup ``cp -r`` kept under its own name), are not
        taken for it (ledger 360). Until ledger 361 only the directory's NAME was
        compared, so a backup named as the template was deleted when given as
        ``workspace=``.
        """

        if os.path.lexists(directory / BUILD_MARKER) or directory != self._default_template().resolve():
            return False
        return self._metadata_names_this_calibration(directory)

    def _metadata_names_this_calibration(self, directory: Path) -> bool:
        """Whether ``directory``'s ``myflopy_pest_metadata.json`` names this project and this run directory.

        What a template a build wrote before ledger 360 -- and every copy pyEMU made of
        it -- carries in place of a :data:`BUILD_MARKER`; where it counts, the caller
        has already pinned the directory's location.
        """

        try:
            metadata = json.loads((directory / METADATA_FILENAME).read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return False
        return (
            isinstance(metadata, dict)
            and metadata.get("project_name") == self.name
            and metadata.get("original_workspace") == str(self.original_workspace.resolve())
        )

    def _clear_stale_template(self):
        """Delete the template and staging directories a previous build wrote -- and nothing else.

        A build writes its template at ``template_workspace`` and stages its copy of
        the run beside it (``<template>.staging``), so before it starts both must be
        gone. It deletes a directory only when it is absent, EMPTY, or one a
        previous build wrote in that role -- its :data:`BUILD_MARKER` names it, or,
        for a template built before the marker existed, its metadata names this
        calibration and this run (ledger 360). Everything is checked before anything
        is deleted, and ``PstFrom`` is made with ``remove_existing=False``, so it
        deletes nothing itself. Before ledger 360 ``workspace=`` was deleted
        WHATEVER it held -- a Project root holding the model's own run, the user's
        home -- by this method and again by ``PstFrom(remove_existing=True)``.

        A directory that is, contains or sits inside the run directory is refused
        whatever it holds (:meth:`_refuse_a_workspace_over_the_run`): one inside it
        is what copied the template into itself until the disk filled (ledger 107),
        one containing it would delete the run.

        A template a serial :meth:`run_ies`/:meth:`prior` is still running in is
        refused too, whatever it holds (ledger 532).

        Raises
        ------
        ValueError
            From :meth:`_refuse_a_workspace_over_the_run`.
        RuntimeError
            From :meth:`_refuse_a_live_run`: a run may still be going in the template.
        FileExistsError
            From :meth:`_clearance`: a directory no build wrote is in the way.
        """

        self._refuse_a_workspace_over_the_run()
        self._refuse_a_live_run(self.template_workspace, "template", method="build", instead="",
                                stacklevel=4)
        cleared = [
            (directory, self._clearance(directory, role))
            for directory, role in (
                (self.template_workspace.resolve(), "template"),
                (self._staging_workspace().resolve(), "staging"),
            )
        ]
        for directory, clearance in cleared:
            if clearance == "empty":
                directory.rmdir()   # refuses, rather than deletes, if a file appeared since
            elif clearance == "built":
                shutil.rmtree(directory)

    def _build_pstfrom(self, staging: Path):
        """Instantiate the underlying ``pyemu.utils.PstFrom`` object over the staged copy.

        ``PstFrom`` opens its log as ``PstFrom.log`` in the CURRENT WORKING
        DIRECTORY, unconditionally (``pyemu.Logger("PstFrom.log")`` in its
        ``__init__``), so every build left one in the caller's directory -- the
        repository root, after a test run. It is constructed with the staging
        directory as the working directory, and the log is then moved into the
        template, where the rest of the run's files are (ledger 356).
        """

        self.pyemu = _import_pyemu()
        from pyemu.utils.pst_from import PstFrom

        # Spatial reference enables geostatistical (grid / pilot-point)
        # parameters. Default to the model's (DISV/Voronoi) modelgrid -- pyEMU
        # reads cell centroids from it to build spatially-correlated parameters --
        # the STAGED model's, whose centres are the template's as written to its
        # files (ledger 356).
        spatial_reference = self.spatial_reference
        if spatial_reference is None:
            spatial_reference = pyemu_spatial_reference(
                getattr(staged_model(self, self.model.name), "modelgrid", None)
            )

        # Both resolved HERE, against the caller's working directory, never inside
        # the block: there a relative path resolves against the STAGING directory,
        # so a relative template -- the default whenever the model's workspace is
        # relative, as in every canonical notebook's `Path('../artifacts/...')` --
        # was made inside the staging copy and deleted with it (review of ledger
        # 356, ledger 357).
        staging = Path(staging).resolve()
        template = self.template_workspace.resolve()
        try:
            with _working_directory(staging):
                self.pf = PstFrom(
                    original_d=staging,
                    new_d=template,
                    # `_clear_stale_template` removed what may be removed; anything still
                    # there is refused by pyEMU rather than deleted (ledger 360).
                    remove_existing=False,
                    longnames=self.longnames,
                    spatial_reference=spatial_reference,
                    zero_based=self.zero_based,
                    start_datetime=self.start_datetime,
                    echo=False,
                    # Apply every multiplier file in the forward run's own process,
                    # never in a per-realization pool (ledger 511).
                    chunk_len=PYEMU_APPLY_CHUNK_LEN,
                )
        except BaseException:
            # `PstFrom` copies the staging directory into the template FIRST and can
            # refuse after -- a run holding `org/` or `mult/` (one made in a copy of a
            # template) does. That template is this build's: it carries the staging
            # copy's own marker, naming THIS staging directory, which nothing but
            # pyEMU's copy of it can. Marked as the template, the next build clears it
            # rather than refusing a folder this build wrote (ledger 361).
            if _written_by_a_build(template, "staging", naming=staging):
                _write_build_marker(template, "template")
            raise
        # First, so a build that fails from here on leaves a template the next one may clear.
        _write_build_marker(template, "template")
        self._adopt_pyemu_log(staging / PYEMU_LOG, template)

    def _adopt_pyemu_log(self, opened_at: Path, template: Path) -> None:
        """Move ``PstFrom``'s open log into ``template`` and keep it logging there.

        The handle is closed before the move -- Windows refuses to move or delete
        an open file -- and reopened for appending, so ``pf.draw`` and the rest of
        the build keep writing to the same log. pyEMU's ``Logger`` keeps the file
        as its plain ``f`` and ``filename`` attributes; were they to move, the log
        stays in the staging directory and goes with it, which still writes
        nothing into the caller's directory.
        """

        logger = getattr(self.pf, "logger", None)
        handle = getattr(logger, "f", None)
        if handle is None or not Path(opened_at).is_file():
            return
        handle.close()
        destination = Path(template) / PYEMU_LOG
        os.replace(opened_at, destination)
        # The logger owns this handle for the life of `cal.pf`, as pyEMU's own did.
        logger.f = open(destination, "a")  # noqa: SIM115
        logger.filename = str(destination)

    def _release_pyemu_log(self) -> None:
        """Close the previous build's ``PstFrom`` log, which lives in the template a rebuild deletes.

        Since ledger 356 the log is IN the template, and a second ``build()`` on one
        project removes the template first (``_clear_stale_template``; before ledger
        360 ``PstFrom(remove_existing=True)`` too); the previous ``cal.pf`` still held the file
        open, which Windows refuses to delete (ledger 357 -- reasoned from Windows'
        file semantics, not measured: there is no Windows machine here). That
        ``cal.pf`` is replaced by the rebuild; a rebuild that fails after this has
        already deleted the template the old ``cal.pf`` described.
        """

        handle = getattr(getattr(self.pf, "logger", None), "f", None)
        if handle is not None and not getattr(handle, "closed", True):
            handle.close()

    def _staging_workspace(self) -> Path:
        """Where the build stages its copy of the run: ``<template>.staging``, beside the template.

        Named from the RESOLVED template (ledger 361): ``workspace="."`` or ``""`` has no
        name of its own, and ``x/..`` named ``x/...staging``, INSIDE ``x``. A template
        with no name at all -- the file-system root -- is refused naming the spelling.
        """

        template = self.template_workspace.resolve()
        if not template.name:
            raise ValueError(
                f"Refusing to build calibration {self.name!r} at workspace="
                f"{str(self.template_workspace)!r}: it is the root of the file system.\n\n"
                + self._where_to_build_instead()
            )
        return template.with_name(f"{template.name}{STAGING_SUFFIX}")

    def _refuse_files_outside_the_run(self) -> None:
        """Refuse a file the simulation reads or writes outside the run directory.

        An INPUT named by a relative path out of it (ledger 421): the staging copy
        holds what the run directory holds; a name like ``../gwf/gwf.hds`` is not in
        it, and in every realization -- each its own copy of the template -- MF6
        resolves it against that copy's directory, so it reads a file beside the
        template (a missing one, or a stale or foreign one) where the run read the
        one beside the run. It was skipped in silence until then (ledger 383's "not
        copied"). An absolute input name is the same file from every copy, and is
        not refused.

        An OUTPUT -- a ``FILEOUT``, or a model name file's ``LIST`` -- named by a
        relative path out of it or by an absolute path (ledger 453): every realization writes that one file, several at once under
        parallel workers, and an absolute name into the run directory overwrites
        the run's own output.
        """

        run = Path(self.original_workspace).resolve()
        problems = []
        outside = _inputs_outside_the_run(self.original_workspace)
        if outside:
            listing = "; ".join(
                f"{name!r} ({owner}), which the run read as "
                f"{Path(os.path.normpath(run / name))}"
                for name, owner in outside
            )
            problems.append(
                f"It reads input from outside its run directory by a relative name: "
                f"{listing}. Every PEST realization runs the model in its own copy of "
                "the template, where that name reads a file beside the copy rather than "
                "this one. Name it by its ABSOLUTE path (read in place by every "
                "realization), or move the file into the run directory and name it "
                "relative to it."
            )
        written = _outputs_outside_the_run(self.original_workspace)
        if written:
            listing = "; ".join(
                f"{name!r} ({owner}), which the run wrote as "
                f"{Path(os.path.normpath(run / name))}"
                for name, owner in written
            )
            problems.append(
                f"It writes output outside its run directory: {listing}. Every PEST "
                "realization runs the model in its own copy of the template, several at "
                "once with parallel workers, and each would write that ONE file -- an "
                "absolute name into the run directory overwrites the run's own output. "
                "Name it relative to the run directory and inside it (e.g. by its file "
                "name alone)."
            )
        if not problems:
            return
        raise ValueError(
            f"The simulation in {run} cannot be copied for PEST as it stands. "
            + " ".join(problems)
            + " Then write and run the simulation again."
        )

    def _stage_model(self) -> Path:
        """Copy the run's INPUT files to a staging directory and make PEST's model of them there.

        Everything the build changes is changed in the COPY (ledger 356): the
        simulation is loaded from it, each declared array target is re-stored one
        array per layer (``relayer_array_target``), every array and list is moved
        into an external file, and the result is written back into the copy, which
        ``PstFrom`` then copies into the template. The run directory is only read,
        and ``self.model`` is never touched -- before, all of this was done to the
        caller's model object and written into the run directory.

        The copy is the simulation's own input files and nothing else of the run
        (``_simulation_input_files``, ledger 383): no output -- the build reads
        every base observation from the RUN directory's output files, and each
        realization's model writes its own -- no myflopy run metadata (``run.json``,
        ``.myflopy/``, a pickled model), and no PEST workspace a run directory holds.
        It copied the whole directory until then, so all of those went into every
        template and every IES worker. An input named by a relative path out of the
        run directory is not in the copy, and the build has refused it before this
        (``_refuse_files_outside_the_run``, ledger 421). The inputs are the
        run's FILES, which is what its output files -- the base observations PEST
        reads -- were computed from.
        FloPy's loader drops LAK ``IMPLICIT`` without a word (ledger 326), so the
        loaded simulation goes through ``keep_implicit``, as myflopy's own loaders
        do.

        Returns
        -------
        pathlib.Path
            The staging directory.

        Raises
        ------
        FileNotFoundError
            When the run directory holds no ``mfsim.nam``: the template is made
            from the run's files, and a model never written has none.
        """

        source = self.original_workspace
        if not (source / "mfsim.nam").is_file():
            raise FileNotFoundError(
                f"No MF6 simulation in {source} (no mfsim.nam). A PEST build makes its "
                "template from the run's own files -- the files its outputs, which the "
                "observations read, were computed from -- and never writes into the run "
                "directory. Write and run the model first (`run.execute()`, or "
                "`model.sim.write_simulation()` then `model.sim.run_simulation()`)."
            )
        staging = self._staging_workspace()
        started = time.perf_counter()
        # `_clear_stale_template` has cleared it; were anything there now, this
        # raises rather than deleting it (ledger 360). Marked BEFORE the copy, so a
        # build that dies copying leaves a directory the next build may clear.
        staging.mkdir(parents=True)
        _write_build_marker(staging, "staging")
        inputs = _simulation_input_files(source)
        for relative in inputs:
            destination = staging / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            # The file's content, never a link: FloPy writes the staged inputs back,
            # and through a link that would write into whatever it points at.
            shutil.copyfile(source / relative, destination)
        # MF6 opens an output in a directory, never makes one (ledger 453).
        for directory in _output_directories(source):
            (staging / directory).mkdir(parents=True, exist_ok=True)
        copied = time.perf_counter()
        simulation = flopy.mf6.MFSimulation.load(sim_ws=str(staging), verbosity_level=0)
        keep_implicit(simulation)
        self._staged_simulation = simulation
        # Before externalizing, not after: a griddata array supplied as a scalar
        # writes one whole-grid file that pyEMU cannot index against a DISV
        # spatial reference. See `relayer_array_target`.
        for spec in self._native_parameter_specs:
            # The FORM of the staged package decides the recipe (an array RCH's
            # per-period arrays, ledger 547) before anything reads it.
            spec.resolve_form(select_packages(
                flopy_model_for(self, spec.declared_recipe), spec.declared_recipe, spec.package))
            refuse_an_undeclared_array(self, spec.recipe, spec.package)
            relayer_array_target(self, spec.recipe, spec.package)
        if self._screened_head_sets():
            # A screened head set weights by each layer's K, read in every
            # realization from the per-layer file MF6 reads (ledger 556);
            # a K given as one whole-grid array has no such file until re-stored.
            relayer_array_target(self, resolve_target("k"))
        simulation.set_all_data_external()
        simulation.write_simulation(silent=True)
        logger.info(
            "Staged %s for PEST in %s: %d input files copied in %.2f s, loaded and "
            "externalized in %.2f s",
            source, staging, len(inputs), copied - started, time.perf_counter() - copied,
        )
        return staging

    @contextlib.contextmanager
    def _quiet_pyemu_context(self):
        """Suppress known harmless pyEMU chatter during workspace assembly."""

        pyemu_warning = _pyemu_warning_class(self.pyemu)
        stdout_buffer = io.StringIO()
        stderr_buffer = io.StringIO()
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=".*not being called directly.*",
                category=pyemu_warning,
            )
            warnings.filterwarnings(
                "ignore",
                message="no adjustable pars",
                category=pyemu_warning,
            )
            warnings.filterwarnings(
                "ignore",
                message="Setting an item of incompatible dtype is deprecated.*",
                category=FutureWarning,
            )
            with contextlib.redirect_stdout(stdout_buffer), contextlib.redirect_stderr(stderr_buffer):
                yield

    def _prepare_one_observation_spec(self, spec):
        """Create simulated files and register one observation spec with pyEMU."""

        if isinstance(spec, (HeadDifferenceObservationSpec, ConcDifferenceObservationSpec)):
            return prepare_difference_observations(self, spec)
        if isinstance(spec, ConcObservationSpec):
            return prepare_conc_observations(self, spec)
        if isinstance(spec, HeadTargetObservationSpec):
            return prepare_head_target_observations(self, spec)
        if isinstance(spec, LakeStageObservationSpec):
            return prepare_lake_stage_observations(self, spec)
        if isinstance(spec, SfrStageObservationSpec):
            return prepare_sfr_stage_observations(self, spec)
        if isinstance(spec, SfrFlowObservationSpec):
            return prepare_sfr_flow_observations(self, spec)
        if isinstance(spec, DrnFlowObservationSpec):
            return prepare_drn_flow_observations(self, spec)
        if isinstance(spec, RivFlowObservationSpec):
            return prepare_riv_flow_observations(self, spec)
        raise TypeError(f"Unsupported observation spec type: {type(spec).__name__}")

    def _refuse_shared_prefixes(self):
        """Refuse two observation or forecast sets under one prefix, before anything is written.

        The prefix names every file a set writes, so the second overwrote the
        first's table, instruction file and map, and the control file the build
        handed back could not run -- nor could its targets be told apart by the
        positional matching of ledger 340.
        """

        prefixes = [
            str(spec.prefix).strip().lower()
            for spec in [*self._observation_specs, *self._forecast_specs]
        ]
        repeated = sorted({prefix for prefix in prefixes if prefixes.count(prefix) > 1})
        if repeated:
            raise ValueError(
                f"Two observation or forecast sets share the prefix {repeated[0]!r}"
                + (f" (and {', '.join(map(repr, repeated[1:]))})" if len(repeated) > 1 else "")
                + ". The prefix names every file a set writes, so the second would "
                "overwrite the first. Pass observe(targets, prefix=...) a distinct "
                "prefix per set."
            )

    def _screened_head_sets(self) -> bool:
        """Whether any observed or forecast head set carries well screens (ledger 556)."""

        return any(
            (isinstance(spec, HeadTargetObservationSpec) and getattr(spec.targets, "screened", False))
            # A difference set of a screened source reads the same T-weighted heads
            # (ledger 565).
            or (isinstance(spec, HeadDifferenceObservationSpec) and spec.targets.source.screened)
            for spec in [*self._observation_specs, *self._forecast_specs]
        )

    def _prepare_observation_specs(self):
        """Create simulated observation files and register them with pyEMU."""

        self._refuse_shared_prefixes()
        prepared: list[dict] = []
        for spec in self._observation_specs:
            prepared.append(self._prepare_one_observation_spec(spec))
        for spec in self._forecast_specs:
            item = self._prepare_one_observation_spec(spec)
            item["is_forecast"] = True
            prepared.append(item)
        self._prepared_observations = prepared

    def _parameter_record(self, spec: NativeParameterSpec) -> dict:
        """What one :meth:`parameterize` call compiled to, for the build's metadata.

        Everything a later reader needs to put a realization's values back on the
        model without this object (ledger 524): the target, the packages
        and files it scales (with each file's package, layer and period, from FloPy's
        records), the style, transform, bounds, physical limits and whether it is an
        offset, the layers and zone labels, its parameter groups, and each of its
        pilot-point sets by the name the template's ``pilot_point_sets.json`` gives it
        (with its group, template, array file, layer and zone). ``parameter_sets`` was
        written before from a list nothing filled, so it was always empty.
        """

        frame = self._native_parameter_frames.get(spec.name)
        groups = (
            sorted({str(group) for group in frame["pargp"]})
            if frame is not None and "pargp" in getattr(frame, "columns", ()) else []
        )
        # One entry per pilot-point set of this call, as the template's
        # pilot_point_sets.json names it (its record there carries the points,
        # factors and limits a write-back rebuilds the field from; ledger 525).
        pilot_points = []
        for item in self._pilot_point_sets:
            if item["spec"] is not spec:
                continue
            record = item["record"]
            pilot_points.append({
                "name": record["name"],
                "pargp": record["name"],
                "template": record["template"],
                "array_files": list(record["array_files"]),
                "layer": record["layer"],
                "zone": record["zone"],
                "period_group": record["period_group"],
            })
            if record["name"] not in groups:
                groups.append(record["name"])
        zones = getattr(spec, "resolved_zone_labels", None)
        return {
            "name": spec.name,
            "target": spec.recipe.canonical,
            "family": spec.recipe.family,
            # Where the values land, for a write-back (ledger 560): the
            # model owning the files, the FloPy dataset scaled (the staged package's
            # form -- an array RCH's ``recharge``, a list one's ``stress_period_data``)
            # and, for a list target, the record fields the forward run writes.
            "model": model_name_for(self, spec.recipe),
            "variable": spec.recipe.variable,
            "fields": list(spec.recipe.written_fields),
            "package_type": spec.recipe.package,
            "package": list(spec.package) if isinstance(spec.package, tuple) else spec.package,
            "packages": sorted({target.package for target in spec.resolved_targets}),
            "style": spec.style,
            "transform": spec.resolved_transform,
            "bounds": [float(value) for value in spec.bounds],
            "physical": (
                [float(value) for value in spec.physical] if spec.physical is not None else None
            ),
            "additive": bool(spec.additive),
            "layers": list(spec.layers) if spec.layers is not None else None,
            "zones": zones,
            # A FEATURE target's (ledger 558): the lakes/reaches named as
            # given, and one record per parameter set they made -- its PEST group,
            # label, package, zero-based features and how many file rows it scales.
            "features": list(spec.features) if spec.features is not None else None,
            "feature_sets": [dict(record) for record in getattr(spec, "feature_sets", [])],
            "correlation": float(spec.correlation) if spec.correlation is not None else None,
            "files": [
                {"filename": target.filename, "package": target.package,
                 "layer": target.layer, "period": target.period}
                for target in spec.resolved_targets
            ],
            "periods": (
                {label: list(values) for label, values in spec.periods.items()}
                if isinstance(spec.periods, dict) else spec.periods
            ),
            # One per parameter set `periods=` made (ledger 549): its label
            # (None without `periods`), the PEST groups the build gave its parameters
            # (one per layer file of a grid array, one per pilot-point set), the
            # periods its files hold the values of, and the files.
            "period_groups": [
                {
                    "label": group.label,
                    "pargps": list(pargps),
                    "periods": [int(period) for period in group.periods],
                    "files": [target.filename for target in group.targets],
                }
                for group, pargps in zip(
                    spec.period_groups, spec.period_group_pargps, strict=True)
            ],
            "pargps": groups,
            "pilot_points": pilot_points,
            "capture": bool(spec.capture),
        }

    def _write_project_metadata(self, *, filename: str | Path | None = None):
        """Persist reopen-friendly run metadata inside the generated workspace.

        ``parameters`` holds one record per :meth:`parameterize` call
        (:meth:`_parameter_record`).
        """

        observation_sets = [
            item["metadata"]
            for item in self._prepared_observations
            if "metadata" in item and item["metadata"] is not None
        ]
        parameters = [self._parameter_record(spec) for spec in self._native_parameter_specs]

        metadata = {
            # 2: each parameter record names its model, dataset and written fields
            # (ledger 560).
            "version": 2,
            "project_name": self.name,
            "model_name": self.model.name,
            "original_workspace": str(self.original_workspace.resolve()),
            "pst_file": Path(filename).name if filename is not None else None,
            "observation_sets": observation_sets,
            "parameters": parameters,
            "capture_fields": list(self._capture_field_specs.values()),
        }
        metadata_path = self.template_workspace / METADATA_FILENAME
        metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
        return metadata_path

    def write(self, filename: str | Path | None = None):
        """Write the built control file to disk.

        Writes ``cal.pst`` -- including any edits made to it since
        :meth:`build` -- in the format pyEMU's ``Pst.write`` picks when given
        no version: version 2 when there are more than 10,000 parameters or
        observations or a version-2-only column is present, else version 1.
        :meth:`build`, :meth:`run_ies` and :meth:`prior` always write
        version 2.

        Parameters
        ----------
        filename : str or pathlib.Path or None, optional
            Control-file name. Only its final component is used: the file
            always lands in the template workspace. ``None`` (default) writes
            ``<name>.pst``, the file :meth:`run_ies`/:meth:`prior` run.

        Returns
        -------
        pathlib.Path
            The file written.

        Raises
        ------
        ValueError
            When :meth:`build` has not run yet.
        """

        if self.pst is None:
            raise ValueError("Call build() before write().")
        target = self.template_workspace / (Path(filename).name if filename else f"{self.name}.pst")
        self.pst.write(target)
        return target

    def draw_prior(self, num_reals: int = 1000, use_specsim: bool = False):
        """Draw a prior parameter ensemble using the underlying ``PstFrom``.

        Calls ``cal.pf.draw(num_reals=..., use_specsim=...)``, which draws every
        parameter in ``cal.pst`` from its bounds, correlated through any
        ``correlation=`` geostructs registered with ``PstFrom``. Array-target
        pilot points are added to ``cal.pst`` after ``PstFrom`` builds and have
        no geostruct there, so this draw includes them but draws each one
        independently; their ``correlation=`` is honoured only by the separate
        draw :meth:`run_ies`/:meth:`prior` make, which replaces these columns
        (ledger 385).

        Parameters
        ----------
        num_reals : int, optional
            Number of realizations to draw. Default ``1000``.
        use_specsim : bool, optional
            Draw the ``grid``-style parameters by pyEMU's spectral simulation
            instead of from the geostruct's covariance matrix. Default ``False``,
            pyEMU's own default and the draw :meth:`run_ies`/:meth:`prior` make,
            which works on every grid. ``True`` is possible only on a REGULAR
            structured grid -- the model's DIS grid (or ``spatial_reference=``)
            with one ``delr`` and one ``delc`` throughout -- and is refused
            (``ValueError``) on any other, a DISV/Voronoi grid included: pyEMU's
            spectral simulation needs row and column spacings, which a vertex
            grid does not have (ledger 384). The default was ``True`` until
            then, and every ``draw_prior()`` on a DISV model raised
            ``AttributeError`` from inside pyEMU.

        Returns
        -------
        pyemu.ParameterEnsemble or None
            The drawn ensemble, one row per realization; ``None`` (pyEMU's
            answer) when the control file has no adjustable parameter.

        Raises
        ------
        ValueError
            When :meth:`build` has not run yet; or when ``use_specsim`` is
            ``True`` and the spatial reference is not a regular structured grid.
        """

        if self.pf is None:
            raise ValueError("Call build() before draw_prior().")
        if use_specsim:
            self._refuse_specsim_off_a_regular_grid()
        return self.pf.draw(
            num_reals=num_reals, use_specsim=bool(use_specsim),
            sigma_range=PRIOR_SIGMA_RANGE, scale_offset=True,
        )

    def _refuse_specsim_off_a_regular_grid(self) -> None:
        """Refuse spectral simulation unless ``PstFrom``'s spatial reference is a regular structured grid.

        pyEMU's ``SpecSim2d`` reads the grid's ``delr``/``delc``. A vertex grid
        has neither: pyEMU then passed empty lists on and died with
        ``AttributeError: 'list' object has no attribute 'mean'``; an irregular
        structured grid gets pyEMU's own refusal. Both are refused here, naming
        the spelling that works (ledger 384).
        """

        import numpy as np

        pyemu = self.pyemu or _import_pyemu()
        reference = getattr(self.pf, "spatial_reference", None)
        delr = getattr(reference, "delr", None)
        delc = getattr(reference, "delc", None)
        if delr is not None and delc is not None:
            delr = np.atleast_1d(np.asarray(delr, dtype=float))
            delc = np.atleast_1d(np.asarray(delc, dtype=float))
            if delr.size and delc.size and pyemu.geostats.SpecSim2d.grid_is_regular(delr, delc):
                return
            what = "a structured grid whose rows or columns are not all one width"
        else:
            what = f"a {type(reference).__name__}, which has no rows or columns"
        raise ValueError(
            f"draw_prior(use_specsim=True): pyEMU's spectral simulation draws only on a "
            f"REGULAR structured grid (one delr, one delc), and this calibration's "
            f"spatial reference is {what}. Draw from the geostruct's covariance "
            f"instead: cal.draw_prior(num_reals, use_specsim=False) -- the default, and "
            f"the draw run_ies()/prior() make."
        )

    # -- native build path ------------------------------------------------

    def _drop_unused_multiplier_apply(self):
        """Remove pyEMU's multiplier-apply command when there are no multipliers.

        ``PstFrom.__init__`` inserts ``apply_list_and_array_pars`` into
        ``pre_py_cmds`` unconditionally (pst_from.py:296-302), but only WRITES
        the ``mult2model_info.csv`` manifest it reads when ``par_dfs`` is
        non-empty (pst_from.py:757-762). A calibration whose only parameter is
        pilot points never calls ``pf.add_parameters`` -- pilot points register
        through template files after ``build_pst`` -- so the manifest is absent
        and every forward run dies with ``FileNotFoundError:
        mult2model_info.csv`` from inside pyEMU, naming nothing the user wrote.

        An empty manifest is not an option: ``apply_list_and_array_pars`` asserts
        ``ddf.shape[0] > 0`` (helpers.py:1953) and has no zero-parameter path. So
        the command is dropped instead, keyed on the same ``par_dfs`` signal
        pyEMU itself branches on.
        """

        if getattr(self.pf, "par_dfs", None):
            return
        self.pf.pre_py_cmds = [
            command for command in self.pf.pre_py_cmds
            if "apply_list_and_array_pars" not in str(command)
        ]

    def _resolve_exe(self) -> str:
        """Return a forward-run-usable MF6 command (absolute path when found)."""

        exe = getattr(self.model.sim, "exe_name", "mf6") or "mf6"
        exe_path = Path(exe)
        if not exe_path.is_absolute():
            candidate = (Path(self.model.workspace) / exe).resolve()
            if candidate.exists():
                return str(candidate)
        if exe_path.exists():
            return str(exe_path.resolve())
        return exe

    def _geostruct_for(self, spec: NativeParameterSpec):
        """Build the pyEMU geostruct for a spatial parameter spec.

        Routes through the single :func:`~...geostats.build_geostruct` builder so
        there is one place that turns a variogram description into a pyEMU
        ``GeoStruct``. The spec's ``correlation`` is the range; ``anisotropy`` and
        ``bearing`` shape the ellipse (e.g. K correlated farther down an alluvial
        valley than across it); ``nugget`` is the unresolved short-scale variance.

        The structure is scaled to a sill of 1 (:func:`~...geostats.unit_sill`):
        pyEMU's draws replace the covariance's diagonal with the variance the bounds
        imply, so a nugget added on top of a sill-1 variogram was discarded there and
        the prior PESTPP-IES received was the same at every nugget (ledger 489).
        Scaled, the bounds give the total variance and the nugget its share.
        """

        return build_geostruct(
            unit_sill(ExpGeoStruct(
                range=float(spec.correlation),
                anisotropy=float(spec.anisotropy),
                bearing=float(spec.bearing),
                nugget=float(spec.nugget),
                transform=spec.resolved_transform,
            ))
        )

    def _attach_native_observation_postprocessors(self):
        """Add post-model functions that regenerate simulated observation CSVs.

        Head and concentration sets get one writer call each. Every named-series
        set -- lake stage, SFR stage/flow, DRN and RIV zones -- shares ONE call of
        ``write_named_series_targets``, which opens the run once for all of them:
        it was one call per set, each reloading the whole simulation, about 9 s a
        set per realization on an externalized Ten Trails template (ledger
        540).
        """

        helper_path = Path(__file__).with_name("forward_run.py")
        # Each period's length, so the forward run can check that the record it
        # reads ENDS its period, as the base run's build did (ledger 336). OC is
        # not a parameter, but ATS can move a STEPS save off a period's end in one
        # realization and not another.
        perlen = period_lengths(self.model)
        named_series: list[dict] = []
        for item in self._prepared_observations:
            head_config = item.get("forward_run_config")
            if head_config is not None:
                self.pf.add_py_function(
                    str(helper_path),
                    "_write_head_target_csv("
                    f"model_name='{self.model.name}', "
                    f"mapping_csv='{head_config['mapping_csv']}', "
                    f"output_csv='{head_config['output_csv']}', "
                    f"perlen={perlen!r}, "
                    # The file the OC declares (ledger 505).
                    f"head_file={head_config.get('head_file')!r}"
                    # A screened set's K files, read after the realization's
                    # parameters are applied (ledger 556).
                    + (f", k_sources={head_config['k_sources']!r}"
                       if head_config.get("k_sources") else "")
                    + ")",
                    is_pre_cmd=False,
                )
                continue
            conc_config = item.get("conc_forward_run_config")
            if conc_config is not None:
                self.pf.add_py_function(
                    str(helper_path),
                    "_write_conc_target_csv("
                    f"model_name='{conc_config['model_name']}', "
                    f"mapping_csv='{conc_config['mapping_csv']}', "
                    f"output_csv='{conc_config['output_csv']}', "
                    f"perlen={perlen!r}, "
                    f"conc_file={conc_config.get('conc_file')!r})",
                    is_pre_cmd=False,
                )
                continue
            difference_config = item.get("difference_forward_run_config")
            if difference_config is not None:
                # One writer for the base run and every realization (ledger
                # 565), on the one period-end reader.
                file_keyword = "head_file" if difference_config["field"] == "head" else "conc_file"
                self.pf.add_py_function(
                    str(helper_path),
                    "_write_target_difference_csv("
                    f"model_name={difference_config['model_name']!r}, "
                    f"mapping_csv={difference_config['mapping_csv']!r}, "
                    f"output_csv={difference_config['output_csv']!r}, "
                    f"periods_csv={difference_config['periods_csv']!r}, "
                    f"kind={difference_config['kind']!r}, "
                    f"field={difference_config['field']!r}, "
                    f"perlen={perlen!r}, "
                    f"{file_keyword}={difference_config.get(file_keyword)!r}"
                    + (f", k_sources={difference_config['k_sources']!r}"
                       if difference_config.get("k_sources") else "")
                    + ")",
                    is_pre_cmd=False,
                )
                continue
            series_config = item.get("named_series_forward_run_config")
            if series_config is not None:
                locations_file = series_config.get("locations_file")
                if locations_file is None:
                    raise ValueError(
                        f"Named-series observation set {item.get('prefix')!r} has no "
                        "saved location definition to regenerate simulated values from."
                    )
                package = series_config.get("package")
                named_series.append(
                    {
                        "kind": str(series_config["kind"]),
                        "locations_file": str(locations_file),
                        "output_csv": str(series_config["output_csv"]),
                        # The package the targets name (ledger 344), so every
                        # realization reads the observations the base run did --
                        # as a plain ``str``: the entry is written into the forward
                        # run by ``repr``, and a name read from a table is often
                        # ``np.str_``, whose repr is ``np.str_('crisp')``, a call
                        # the generated script would have to resolve.
                        "package": None if package is None else str(package),
                    }
                )
                continue
            raise NotImplementedError(
                f"Observation set {item.get('prefix')!r} has no forward-run "
                "post-processor; supported kinds are head-target and named-series "
                "(lake stage, SFR stage/flow, DRN seepage, RIV exchange) observations."
            )
        if named_series:
            # The calibrated model by NAME: without it the forward run reloads the
            # simulation's first GWF model, another model's lake or zone in a
            # simulation of several (ledger 340).
            self.pf.add_py_function(
                str(helper_path),
                f"write_named_series_targets(sets={named_series!r}, "
                f"model_name='{self.model.name}')",
                is_pre_cmd=False,
            )

    def _forecast_observation_names(self) -> list[str]:
        """The observations the forecast sets name, as ``finalize_observations`` read them.

        Each forecast target's observation is the one the instruction file reads at
        its row and column -- not a name rebuilt from its key's spelling, which
        dropped a forecast keyed ``1.0`` where pyEMU wrote ``per:1`` (ledger 340).
        """

        names: list[str] = []
        for item in self._prepared_observations:
            if item.get("is_forecast"):
                names.extend(item.get("forecast_names", []))
        return list(dict.fromkeys(names))

    def _apply_forecasts(self):
        """Zero forecast weights and register them as PEST++ forecasts."""

        names = [name for name in self._forecast_observation_names() if name in self.pst.observation_data.index]
        if not names:
            return
        self.pst.observation_data.loc[names, "weight"] = 0.0
        self.pst.pestpp_options["forecasts"] = ",".join(names)

    def _add_capture_field_observations(self, spec: NativeParameterSpec):
        """Record a parameter's resolved per-cell field as zero-weight observations."""

        recipe = spec.recipe
        base = f"{spec.name}field"
        # Each external array file is captured under its own prefix so that
        # multi-layer DISV K (one file per layer, each indexed 0..ncpl-1) does
        # not collide -- a single shared prefix would make layer 1 and layer 2
        # produce identical observation names.
        #
        # `layer_prefixes` stays EMPTY for a whole-grid array (MST porosity is
        # one file of nlay*ncpl values, with no LAYERED keyword). That is what
        # selects the flat `layer * ncpl + cell` reader in `IesResults._cell_map`
        # -- recording it as `{0: base}` instead would label every layer's
        # values "layer 0" and hand `plot_field` nlay*ncpl values for ncpl cells.
        layer_prefixes: dict[int, str] = {}
        for target in spec.resolved_targets:
            filename = target.filename
            if recipe.family == "array":
                # The layer FloPy's records name the file for (ledger 517),
                # never one parsed back out of the file name.
                prefix = f"{base}l{target.layer}" if target.layer is not None else base
                if target.layer is not None:
                    layer_prefixes[target.layer] = prefix
                self.pf.add_observations(filename, prefix=prefix)
            else:
                # The value column on THIS grid: one further right on DIS, whose
                # cellid is three columns (ledger 420).
                index_cols, use_col = list_columns(
                    recipe, flopy_model_for(self, recipe).modelgrid.grid_type
                )
                self.pf.add_observations(
                    filename, prefix=base, index_cols=index_cols, use_cols=[use_col]
                )
        self._capture_field_specs[spec.name] = {
            "prefix": base,
            "layer_prefixes": layer_prefixes,
            "target": recipe.canonical,
            "family": recipe.family,
            "name": spec.name,
        }

    def _finalize_capture_fields(self):
        """Zero the weights of captured-field observations after the pst is built."""

        if not self._capture_field_specs:
            return
        obs = self.pst.observation_data
        index = obs.index.to_series()
        for info in self._capture_field_specs.values():
            prefixes = list(info.get("layer_prefixes", {}).values()) or [info["prefix"]]
            for prefix in prefixes:
                mask = index.str.contains(f"oname:{prefix.lower()}_", regex=False)
                obs.loc[mask.to_numpy(), "weight"] = 0.0

    def build(self, filename: str | Path | None = None, *, noptmax: int = 0):
        """Compile the declarations into a runnable PEST(++) control file.

        In order it: deletes the template and staging directories a previous
        build wrote (and refuses, deleting nothing, when either is a folder no
        build wrote -- see ``workspace=`` on :class:`PestProject`, ledger 360);
        copies the simulation's input files -- what ``mfsim.nam`` leads MF6 to
        read, and nothing else of the run directory (ledger 383) -- to a
        staging directory beside the template, loads the simulation from that copy and writes it back with
        every array and list in an external file (``set_all_data_external``);
        creates a ``pyemu.utils.PstFrom`` over the staged copy; registers every
        :meth:`parameterize` call natively (``pf.add_parameters``) and every
        :meth:`observe`/:meth:`forecast` set; adds the MF6 run command and the
        observation post-processors to ``forward_run.py``; builds the ``.pst``;
        and removes the staging directory. The multiplier machinery is pyEMU's
        own ``apply_list_and_array_pars`` in the forward run, so parameter
        application is handled by pyEMU rather than hand-rolled code.

        **The run directory and the model object are left exactly as they were**
        (ledger 356): the build reads the run's files and output, and writes only
        into the template (and the staging directory it removes). It used to
        re-store the caller's model's arrays and write the simulation back into the
        run directory, inline arrays turned into external files. The template is
        made from the run's FILES -- the inputs its output files, which the base
        observations are read from, were computed from -- so an edit made to
        ``model`` in memory after the run does not reach it. pyEMU's
        ``PstFrom.log`` is written into the template, not the working directory.

        After this call the workspace contains a self-contained PEST setup:
        the ``.pst``, ``forward_run.py``, template/instruction files, and the
        ``mult/`` multiplier files. Run it with PEST++ (e.g. ``pestpp-ies`` /
        ``pestpp-glm``) or via :meth:`run_ies` / :meth:`prior`.

        Parameters
        ----------
        filename : str or pathlib.Path or None, optional
            Output ``.pst`` filename (default ``"<name>.pst"``). Only its final
            component is used: it is written into the template workspace and
            recorded in the run's metadata, whose stem is the case
            ``model.pest_runs[...].review()`` opens. :meth:`run_ies` and
            :meth:`prior` always write and run ``<name>.pst`` instead, so a
            stem other than ``name`` leaves ``review()`` reading the wrong
            case.
        noptmax : int, optional
            Initial ``NOPTMAX`` written to the control file. ``0`` (default)
            means "run the model once and compute residuals" -- the standard
            way to sanity-check a fresh setup before launching a real run.
            :meth:`run_ies`/:meth:`prior` overwrite it.

        Returns
        -------
        pyemu.Pst
            The built control object (also available as ``cal.pst``). The
            underlying ``PstFrom`` is ``cal.pf`` for advanced edits.

        Raises
        ------
        ValueError
            When an observation or forecast target matches no row of its
            simulated table (naming it and the rows the table has), or a head or
            concentration target's location lands on no active cell (ledger
            341); when a target's ``time`` is not a whole number, or two values
            name one target in one period; when two sets share a prefix; when
            a head or concentration ``simulated_values`` table's columns are not
            the forward run's, in its order (ledger 340); or when a target's
            group label, or a set's prefix, cannot be its PEST observation group
            (ledger 358; see :meth:`observe`). Also, before anything is deleted,
            when ``workspace=`` or the staging directory beside it is, contains or
            sits inside the model's run directory, whatever it holds (ledger 360);
            and, before anything is deleted, when an input the simulation reads is
            named by a RELATIVE path leading out of the run directory -- a GWT
            model's FMI reading ``../gwf/gwf.hds`` -- naming each and what names it
            (ledger 421). Every realization runs in its own copy of the
            template, where such a name reads a file beside the template; an
            absolute name is read in place and is not refused. Likewise when an
            output (``FILEOUT``) is named by a relative path leading out of the
            run directory or by an absolute path (ledger 453): every
            realization would write that one file -- and the same for a model name
            file's ``LIST <file>``, its listing (ledger 453).
        FileNotFoundError
            When the model's workspace holds no ``mfsim.nam``: the template is
            made from the run's files, and a model never written has none.
        FileExistsError
            When ``workspace=``, or the staging directory beside it, is a file or
            a non-empty folder no build wrote (ledger 360): the build deletes both
            directories first, so it refuses rather than delete what it did not
            write, before anything is deleted. Use the default, a new or empty
            folder, or delete it yourself.
        RuntimeError
            Before anything is deleted, when a serial :meth:`run_ies`/:meth:`prior`
            may still be running in the template (ledger 532).
        """

        target_name = Path(filename).name if filename else f"{self.name}.pst"
        self._refuse_shared_prefixes()
        self._refuse_files_outside_the_run()
        self._release_pyemu_log()
        self._clear_stale_template()
        staging = self._staging_workspace()
        try:
            self._stage_model()
            self._build_pstfrom(staging)
            self._prepare_observation_specs()
            self._carried_river_files = set()
            # The riv.stage lower bounds compiled onto each file so far (ledger
            # 568): offsets on one file add in the forward run.
            self._river_offset_lower_bounds = {}
            # Pilot-point sets are collected afresh by every build.
            self._pilot_point_sets = []
            with self._quiet_pyemu_context():
                for spec in self._native_parameter_specs:
                    if spec.recipe.family == "array" and spec.style == "pilotpoints":
                        add_pilot_point_parameter(self, spec)
                    else:
                        add_native_parameter(self, spec)
                    if spec.capture:
                        self._add_capture_field_observations(spec)
                # After every parameter: one pilot-point command per array file,
                # appended behind pyEMU's apply and multiplying what it wrote
                # (ledger 525).
                emit_pilot_point_commands(self)
            self.pf.mod_sys_cmds.append(self._resolve_exe())
            self._drop_unused_multiplier_apply()
            self._attach_native_observation_postprocessors()
            self._write_project_metadata(filename=target_name)
            with self._quiet_pyemu_context():
                self.pst = self.pf.build_pst(filename=target_name)
            if self._pilot_point_sets:
                register_pilot_point_parameters(self)
            self.pst.model_command = [
                build_forward_run_command(self.template_workspace, sys.executable)
            ]
            finalize_observations(self, self._prepared_observations)
            self._apply_forecasts()
            self._finalize_capture_fields()
            self.pst.control_data.noptmax = int(noptmax)
            with self._quiet_pyemu_context():
                self.pst.write(self.template_workspace / target_name, version=2)
        finally:
            # The template holds everything the staging copy was made for. Removed only
            # as the copy THIS build made and marked: a build that failed before making
            # it -- or on finding a folder already there -- leaves that folder alone.
            self._staged_simulation = None
            if _written_by_a_build(staging.resolve(), "staging"):
                shutil.rmtree(staging, ignore_errors=True)
        return self.pst

    def settings(self) -> PestSettings:
        """Return a printable snapshot of the resolved calibration configuration.

        The facade hides pyEMU boilerplate but never hides *state*: this lets
        you review exactly what was declared (every parameter with its style,
        bounds, physical limits, transform, ``periods=`` grouping and ``features=``;
        every observation and forecast)
        and -- once :meth:`build` has run -- the resulting control-file counts
        (number of parameters and groups, observations, nonzero-weight
        observations, forecasts, and ``NOPTMAX``).

        Returns
        -------
        PestSettings
            A dataclass with a readable ``str``/``repr`` (``print(cal.settings())``)
            plus ``.parameter_frame()`` / ``.observation_frame()`` accessors for
            the declared parameters and observations as ``pandas`` tables.

        Examples
        --------
        >>> print(cal.settings())          # before build: declared config
        >>> pst = cal.build("calib.pst")
        >>> print(cal.settings())          # after build: + control-file counts
        """

        parameters = [
            {
                "target": spec.recipe.canonical,
                "package": spec.package,
                "style": spec.style,
                "bounds": tuple(spec.bounds),
                "physical": tuple(spec.physical) if spec.physical is not None else None,
                "transform": spec.resolved_transform,
                "additive": bool(spec.additive),
                "name": spec.name,
                # `periods=` as declared: None, "month" or {label: periods}.
                "periods": spec.periods,
                # `features=` as declared, and after a build one row per parameter set
                # the features made (ledger 558).
                "features": list(spec.features) if spec.features is not None else None,
                "feature_sets": [dict(item) for item in spec.feature_sets],
                # After a build, array-target pilot points: one row per set with its
                # point count (ledger 527).
                "pilot_point_sets": [dict(item) for item in spec.pilot_point_sets],
            }
            for spec in self._native_parameter_specs
        ]

        def _obs_rows(specs, prepared_filter):
            """Summary rows for observation/forecast ``specs``, merged with prepared build info."""

            rows = []
            prepared = {
                item.get("prefix"): item
                for item in self._prepared_observations
                if bool(item.get("is_forecast")) == prepared_filter
            }
            for spec in specs:
                prefix = getattr(spec, "prefix", "?")
                item = prepared.get(prefix)
                target_frame = item.get("target_frame") if item else None
                rows.append(
                    {
                        "prefix": prefix,
                        "kind": type(spec).__name__.replace("ObservationSpec", "").lower(),
                        "n": int(len(target_frame)) if target_frame is not None else None,
                    }
                )
            return rows

        observations = _obs_rows(self._observation_specs, prepared_filter=False)
        forecasts = _obs_rows(self._forecast_specs, prepared_filter=True)

        built = self.pst is not None
        npar = npar_groups = nobs = nnz_obs = n_forecasts = noptmax = None
        if built:
            par_data = self.pst.parameter_data
            npar = int(self.pst.npar)
            npar_groups = int(par_data["pargp"].nunique())
            nobs = int(self.pst.nobs)
            nnz_obs = int(self.pst.nnz_obs)
            forecast_option = self.pst.pestpp_options.get("forecasts", "")
            n_forecasts = len([name for name in forecast_option.split(",") if name]) if forecast_option else 0
            noptmax = int(self.pst.control_data.noptmax)

        return PestSettings(
            name=self.name,
            model_name=self.model.name,
            original_workspace=str(self.original_workspace),
            template_workspace=str(self.template_workspace),
            start_datetime=str(self.start_datetime),
            parameters=parameters,
            observations=observations,
            forecasts=forecasts,
            built=built,
            npar=npar,
            npar_groups=npar_groups,
            nobs=nobs,
            nnz_obs=nnz_obs,
            n_forecasts=n_forecasts,
            noptmax=noptmax,
        )

    # -- running PEST++ ---------------------------------------------------

    def _resolve_pestpp(self, exe: str):
        """Return ``(exe_name, exe_dir)`` for a PEST++ tool next to the MF6 binary."""

        name = exe[:-4] if str(exe).lower().endswith(".exe") else str(exe)
        mf6 = Path(self._resolve_exe())
        if mf6.is_absolute() and mf6.parent.exists() and (mf6.parent / f"{name}.exe").exists():
            return name, str(mf6.parent)
        return name, None

    @contextlib.contextmanager
    def _pestpp_environment(self, exe_dir):
        """The environment PEST++ and its workers start in, for the duration of a run.

        ``PATH`` gains ``exe_dir`` so PEST++ tools resolve. ``PYTHONPATH`` is made
        absolute (:func:`myflopy._pythonpath.pythonpath_for_children`, the rule the
        test suite and the reproducers use): PEST++ starts the forward run in the
        template, where a relative entry -- a notebook started under
        ``PYTHONPATH=src`` -- names another directory, and the forward run imported
        another myflopy (ledger 454). Both are restored afterwards.
        """

        pythonpath = pythonpath_for_children()
        original = {name: os.environ.get(name) for name in ("PATH", "PYTHONPATH")}
        path = original["PATH"] or ""
        if exe_dir and exe_dir not in path:
            os.environ["PATH"] = exe_dir + os.pathsep + path
        if pythonpath is not None:
            os.environ["PYTHONPATH"] = pythonpath
        try:
            yield
        finally:
            for name, value in original.items():
                if value is None:
                    os.environ.pop(name, None)
                else:
                    os.environ[name] = value

    def run_ies(
        self,
        *,
        reals: int = 50,
        iterations: int = 3,
        workers: int | None = None,
        noise: bool = True,
        bad_phi_sigma: float | None = None,
        master_dir: str | Path | None = None,
        overwrite: bool = False,
        keep_workers: bool = False,
        exe: str = "pestpp-ies",
        overdue_giveup_fac: float | None = DEFAULT_OVERDUE_GIVEUP_FAC,
        **pestpp_options,
    ):
        """Run iterative ensemble smoother history matching with PESTPP-IES.

        This is the one-line entry point to ensemble calibration and "free"
        uncertainty analysis. It sets sensible PESTPP-IES options on the built
        control file, launches the run (serial, or parallel across ``workers``
        agents), and returns an :class:`~myflopy.modflow.mf6.pest.ies.IesResults`
        for assessing the outcome (phi convergence, ensembles vs observations,
        posterior forecast distributions).

        Call :meth:`build` first. Reasonable defaults make ``cal.run_ies()`` a
        valid first call; tune from there. A run interrupted part-way (a killed
        kernel, a reboot) continues from its last completed iteration with
        :meth:`resume_ies`.

        ``noise``, ``bad_phi_sigma`` and ``**pestpp_options`` are written onto
        ``cal.pst.pestpp_options`` and STAY there: ``noise=False``, a
        ``bad_phi_sigma`` or an extra option from an earlier
        :meth:`run_ies`/:meth:`prior` call is still set on the next one, whose
        defaults do not clear it. ``overdue_giveup_fac`` is set (or removed) by
        every call.

        The directory a run writes into holds that run's outputs only (ledger
        507). A serial run first deletes from the template every output
        an earlier run of this calibration left there (each iteration's ensembles
        and ``.base.par``, the phi files, the noise ensemble, a :meth:`resume_ies`
        archive -- except a file a PEST++ option names), so a rerun with fewer iterations, or a :meth:`prior`
        after a run, is not reviewed with the earlier run's last iteration as its
        posterior; it says so with a ``UserWarning`` (ledger 532). A
        parallel run's master is copied from the template with those outputs left
        out (ledger 532), and the template keeps them.

        **A run in progress is never touched** (ledger 532). Every launch
        writes a ``myflopy_pest_run.lock`` into the directories it runs in -- the
        template (serial), or the master and its worker directory (parallel) --
        naming the host, the process and when it started, and removes it when the
        run returns. A launch that would delete or reuse a directory whose lock names
        a process still running, a lock written on another host, or one where a PEST++
        process is working (read from ``/proc`` on Linux), is refused
        (``RuntimeError``) whatever ``overwrite`` says. A lock whose process is gone
        is STALE: it is said so (``UserWarning``) and ignored.

        Parameters
        ----------
        reals : int, optional
            Number of realizations in the ensemble (``ies_num_reals``). More is
            better for posterior coverage but costs more model runs; 50 (the
            default) is a reasonable starting point, production runs often use
            100-300. Also the size of the correlated prior drawn for
            ``correlation=`` parameters.
        iterations : int, optional
            Number of smoother iterations (``NOPTMAX``). Default ``3``, a
            common choice; PESTPP-IES may stop earlier if it converges.
        workers : int or None, optional
            Number of parallel agents. ``None`` (default) or ``1`` runs
            serially in the template workspace, where the results then are;
            ``>1`` starts that many workers with pyEMU's ``start_workers``
            and collects the results in ``master_dir``. A serial run deletes an
            earlier run's outputs from the template first (above).
        noise : bool, optional
            Whether to carry measurement noise into the posterior (recommended,
            the default ``True``, which sets nothing). ``False`` sets
            ``ies_no_noise``.
        bad_phi_sigma : float or None, optional
            Adaptive rejection threshold (``ies_bad_phi_sigma``);
            values of 1.5 (aggressive) to 2.5 (tolerant) help highly nonlinear
            problems by dropping lagging realizations. ``None`` (default) sets
            nothing.
        master_dir : str or pathlib.Path or None, optional
            Master directory for a parallel run, where its results are collected;
            ``None`` (default) uses ``<name>_ies_master`` beside the template
            workspace, which is where ``model.pest_runs`` looks. Its workers are made
            in ``<master_dir>.workers/`` as ``worker_0`` .. ``worker_<workers-1>``,
            so no two masters share one (ledger 532; until then every
            calibration of a model shared ``<model>.pest/worker_<n>``). Unused by a
            serial run. A launch deletes the master and each worker it uses, then
            copies the template into them, so each must be absent, empty, or a copy
            of THIS calibration's template; a master holding an earlier run's
            results is also refused unless ``overwrite=True``. Anything else is
            refused before anything is written (``FileExistsError``), a master or
            worker directory that is, contains or sits inside the run directory or
            the template is refused whatever it holds (``ValueError``), and one a
            run is still going in is refused (``RuntimeError``).
        overwrite : bool, optional
            Parallel run: ``True`` deletes a master holding an earlier run's results
            (finished, or interrupted) and runs again there. Default ``False``
            refuses it (``FileExistsError``), naming the folder, whether its run
            finished, and the spellings: this one, another ``master_dir=``, or
            :meth:`resume_ies` for an interrupted history match
            (ledger 532; ledger 361 deleted it unasked). It never reaches a master
            another calibration's template made, nor one a run is still going in.
            Serial run: the template's earlier results are deleted either way, and
            ``True`` only silences the ``UserWarning`` saying so.
        keep_workers : bool, optional
            Parallel run: keep the worker directories (``<master_dir>.workers/``)
            after the run. Default ``False`` deletes them when the run finishes
            cleanly -- each is a disposable copy of the template, one per worker, and
            :meth:`resume_ies` makes them again; a run that fails or is interrupted
            keeps them either way, for diagnosis (ledger 532). The master
            and its results are always kept. Unused by a serial run.
        exe : str, optional
            PEST++ IES executable name, default ``"pestpp-ies"`` (a trailing
            ``.exe`` is dropped). When the model's MF6 executable resolves to an
            absolute path and a file ``<exe>.exe`` sits beside it, that
            directory is prepended to ``PATH`` for the run (unless ``PATH``
            already contains it); otherwise the name is resolved on ``PATH`` as
            it is.
        overdue_giveup_fac : float or None, optional
            PESTPP-IES's ``overdue_giveup_fac``: in a PARALLEL run, a model run
            taking longer than this many times the average run is killed and its
            realization counted failed (dropped from the ensemble). Default
            :data:`DEFAULT_OVERDUE_GIVEUP_FAC`, ``5.0`` (ledger 510):
            PESTPP-IES resets an unset factor to ``2.0``, and on a loaded machine a
            NEWTON model's slow, hard-to-converge realizations are killed at twice
            the average -- selectively, so the ensemble loses them. A larger value
            tolerates slower runs; ``None`` removes the option, giving PESTPP-IES's
            own ``2.0``. Set on ``cal.pst.pestpp_options`` by every call. A serial
            run kills nothing.
        **pestpp_options
            Any additional ``pst.pestpp_options`` to set (e.g.
            ``ies_localizer="loc.mat"``), passed straight through and applied
            after the options above, so they override them. When a ``grid``
            parameter, or an array-target ``pilotpoints`` parameter, carries
            ``correlation=``, a call naming no ``ies_par_en`` draws a correlated
            prior ensemble into ``<name>.prior_par.csv`` and sets ``ies_par_en`` to
            it (replacing one an earlier call left on ``cal.pst``); an
            ``ies_par_en`` passed HERE is used as given and nothing is drawn (ledger
            509 -- it was replaced by the draw too). ``ies_restart_obs_en``
            -- the simulated values of an ``ies_par_en``'s realizations, which
            PESTPP-IES then does not rerun -- must come with the ``ies_par_en``
            that made them: set without one while a prior would be drawn, it is
            refused (``ValueError``).

        Returns
        -------
        IesResults
            Reader/visualizer for the completed run, opened on the directory
            holding its results with case ``name``.

        Raises
        ------
        ValueError
            When :meth:`build` has not run, or ``cal.pst`` has been put in
            ``regularization`` mode, which PESTPP-IES cannot parse; when
            ``ies_restart_obs_en`` is set while a correlated prior would be drawn
            into ``ies_par_en`` (ledger 509); or, for a parallel run, when the
            master or its worker directory is, contains or sits inside the run
            directory or the template (ledger 361).
        RuntimeError
            Before anything is written, when a run may still be going in the
            template, or (parallel) in the master or its worker directory: its
            lock names a process still running or written on another host, or a
            PEST++ process works there (ledger 532). Also, before PEST++
            starts, when a RELATIVE ``PYTHONPATH`` entry, read against the working
            directory, is not on ``sys.path`` -- the process has changed directory
            since it started. PEST++ is handed ``PYTHONPATH`` with every relative
            entry made absolute, so the forward run, started in the template,
            imports this process's myflopy
            (``myflopy._pythonpath.pythonpath_for_children``, ledger 454).
        FileExistsError
            For a parallel run, before anything is written, when the master or a
            worker directory is a file or a non-empty folder that is not a copy of
            this calibration's template (ledger 361), or the master holds an earlier
            run's results and ``overwrite`` is ``False`` (ledger 532).
        """

        results_dir = self._launch_pestpp_ies(
            reals=reals,
            noptmax=iterations,
            workers=workers,
            exe=exe,
            master_suffix="ies_master",
            master_dir=master_dir,
            noise=noise,
            bad_phi_sigma=bad_phi_sigma,
            overdue_giveup_fac=overdue_giveup_fac,
            overwrite=overwrite,
            keep_workers=keep_workers,
            pestpp_options=pestpp_options,
        )
        return self._review(results_dir)

    def prior(
        self,
        *,
        reals: int = 50,
        workers: int | None = None,
        noise: bool = True,
        master_dir: str | Path | None = None,
        overwrite: bool = False,
        keep_workers: bool = False,
        exe: str = "pestpp-ies",
        overdue_giveup_fac: float | None = DEFAULT_OVERDUE_GIVEUP_FAC,
        **pestpp_options,
    ):
        """Run the prior parameter ensemble once (prior Monte Carlo).

        Before history matching, this draws ``reals`` parameter sets from the
        prior and runs the model once for each -- "everything the model could
        plausibly do given only expert knowledge, without looking at the data."
        It is the cheap go/no-go check that the *prior brackets the observations*;
        if it does not, you have prior-data conflict (a model problem, not
        something history matching can fix), and history matching is premature.

        Mechanically this is PESTPP-IES with ``NOPTMAX=-1`` (evaluate the prior
        ensemble and stop). Call :meth:`build` first. An interrupted prior cannot be
        resumed -- its simulated values (``<name>.0.obs.csv``) are written only when
        the whole sweep is done -- so it is run again.

        Options persist on ``cal.pst`` between calls exactly as in
        :meth:`run_ies` (a ``bad_phi_sigma`` set there is still set here), and a
        serial prior deletes an earlier run's outputs from the template first, as
        :meth:`run_ies` describes -- a serial ``run_ies``'s results there included,
        with a ``UserWarning``. A run in progress is never touched, as
        :meth:`run_ies` describes (ledger 532).

        Parameters
        ----------
        reals : int, optional
            Number of prior realizations (``ies_num_reals``). Default ``50``,
            plenty for a visual bracketing check. Also the size of the
            correlated prior drawn for ``correlation=`` parameters.
        workers : int or None, optional
            Parallel agents. ``None`` (default) or ``1`` runs serially in the
            template workspace; ``>1`` starts that many workers and collects the
            results in ``master_dir``.
        noise : bool, optional
            Generate the measurement-noise ensemble too (default ``True``,
            which sets nothing). ``False`` sets ``ies_no_noise``.
        master_dir : str or pathlib.Path or None, optional
            Master directory for a parallel run. ``None`` (default) uses
            ``<name>_prior_master`` beside the template workspace, which is where
            ``model.pest_runs`` looks; the workers are made in
            ``<master_dir>.workers/``, so a prior and a history match of one
            calibration -- in two kernels, say -- never share one. Unused by a
            serial run. What a launch deletes there, and what is refused, is as for
            :meth:`run_ies`'s ``master_dir`` (ledger 361, ledger 532).
        overwrite : bool, optional
            ``True`` deletes a master holding an earlier prior's results and runs
            again there; default ``False`` refuses it (``FileExistsError``), exactly
            as :meth:`run_ies`'s ``overwrite``. For a serial run it only silences
            the ``UserWarning`` saying the template's earlier results were deleted.
        keep_workers : bool, optional
            Keep a parallel run's worker directories after it finishes cleanly;
            default ``False`` deletes them, exactly as :meth:`run_ies`'s
            ``keep_workers``. Unused by a serial run.
        exe : str, optional
            PEST++ IES executable, default ``"pestpp-ies"``, resolved as in
            :meth:`run_ies`.
        overdue_giveup_fac : float or None, optional
            How many times the average run a parallel model run may take before
            PESTPP-IES kills it and drops its realization; default
            :data:`DEFAULT_OVERDUE_GIVEUP_FAC`, ``5.0``, and ``None`` removes the
            option (PESTPP-IES's own ``2.0``), exactly as in :meth:`run_ies`
            (ledger 510).
        **pestpp_options
            Additional ``pst.pestpp_options``, passed straight through and
            applied after ``noise``; ``ies_par_en`` and ``ies_restart_obs_en``
            behave as in :meth:`run_ies` (a passed ``ies_par_en`` is used as given,
            and a restart ensemble without one is refused while a correlated prior
            would be drawn).

        Returns
        -------
        IesResults
            Use :meth:`~...ies.IesResults.plot_prior_vs_obs` (grey prior spaghetti
            vs measured) and :meth:`~...ies.IesResults.conflict` to inspect it.

        Raises
        ------
        ValueError
            When :meth:`build` has not run, or ``cal.pst`` has been put in
            ``regularization`` mode; for a restart ensemble a drawn prior would
            orphan, as in :meth:`run_ies` (ledger 509); or as
            :meth:`run_ies` raises it for a parallel run's master and worker
            directories (ledger 361).
        RuntimeError
            As :meth:`run_ies` raises it: a run may still be going in the template,
            the master or its worker directory (ledger 532), or a relative
            ``PYTHONPATH`` entry the process has changed directory away from
            (ledger 454).
        FileExistsError
            As :meth:`run_ies` raises it, for a parallel run's master and worker
            directories (ledger 361, ledger 532).
        """

        results_dir = self._launch_pestpp_ies(
            reals=reals,
            noptmax=-1,
            workers=workers,
            exe=exe,
            master_suffix="prior_master",
            master_dir=master_dir,
            noise=noise,
            bad_phi_sigma=None,
            overdue_giveup_fac=overdue_giveup_fac,
            overwrite=overwrite,
            keep_workers=keep_workers,
            pestpp_options=pestpp_options,
        )
        return self._review(results_dir)

    #: The PEST++ options :meth:`resume_ies` sets itself, to restart from the
    #: interrupted run's ensembles (ledger 533).
    _RESUME_OPTIONS = ("ies_par_en", "ies_restart_par_en", "ies_restart_obs_en", "ies_obs_en")

    def resume_ies(
        self,
        *,
        workers: int | None = None,
        iterations: int | None = None,
        master_dir: str | Path | None = None,
        keep_workers: bool = False,
        exe: str = "pestpp-ies",
        overdue_giveup_fac: float | None = DEFAULT_OVERDUE_GIVEUP_FAC,
        **pestpp_options,
    ):
        """Continue an interrupted :meth:`run_ies` from its last completed iteration.

        A multi-day history match stopped part-way -- a killed kernel, a reboot, a
        Ctrl-C -- keeps every iteration it completed. This restarts PESTPP-IES from
        the last one, in the same directory, with the same control file, and leaves
        one run behind: iterations ``0..N``, one phi file, reviewed exactly as an
        uninterrupted run (ledger 533). It also extends a FINISHED run by
        more iterations (``iterations=`` beyond its own).

        PESTPP-IES has no restart switch (``/r`` is refused: "/r option not supported
        by pestpp-ies", 5.2.16, measured); it restarts from ensembles. This archives
        the run's outputs into ``<name>.resumed_from_<n>/`` in the run directory,
        launches PESTPP-IES with ``ies_par_en`` the run's prior (iteration 0, so the
        regularization is the run's own), ``ies_restart_par_en`` /
        ``ies_restart_obs_en`` iteration ``n``'s parameters and simulated values
        (not rerun), ``ies_obs_en`` the run's noise ensemble, ``ies_initial_lambda``
        the Marquardt lambda the run's ``.rec`` shows in force for the next iteration,
        and ``NOPTMAX`` the iterations left; then -- when PESTPP-IES returns, finished
        or not -- moves the archived iterations ``0..n`` back, renumbers the resumed
        ones after them, joins the phi files (``total_runs`` continues) and puts the
        control file back as the run was launched: ``NOPTMAX`` its TOTAL iterations, the
        restart options as they were -- so a run interrupted again reads as interrupted
        and resumes to the same total. The archive keeps the interrupted run's ``.rec``
        and ``.log``. Measured on the canonical
        testing model, 6 realizations, killed after iteration 1 of 3
        (``docs/reproducers/ies_resume_interrupted_run.py``): the resumed iteration 2
        read phi 1.71477 against the uninterrupted run's 1.71448 -- 1.65128 without the
        lambda carried. Not carried: PESTPP-IES's other internal state (its
        best-phi history restarts at iteration ``n``), and run statistics in the
        ``.rec``.

        A resume killed too hard for that clean-up to run leaves the archive
        unspliced; the next :meth:`resume_ies` splices it first, and puts the control
        file back, before it reads what the run was launched for. A later serial
        :meth:`run_ies`/:meth:`prior` in the template deletes the archives with the
        earlier run's other outputs, so no resume splices them into a new run.

        Parameters
        ----------
        workers : int or None, optional
            Where the run to resume ran, and how this one runs: ``None`` (default)
            or ``1`` resumes a SERIAL run, in the template workspace; ``>1`` resumes
            a PARALLEL run in its master directory with that many workers (any
            number, not necessarily the interrupted run's), made afresh in
            ``<master_dir>.workers/`` from the template.
        iterations : int or None, optional
            The run's TOTAL number of iterations, counting those already completed.
            Default ``None`` uses the ``NOPTMAX`` the interrupted run was launched
            with. One no larger than the last completed iteration is refused
            (``ValueError``); a larger one extends a finished run.
        master_dir : str or pathlib.Path or None, optional
            The master directory of the parallel run to resume; ``None`` (default)
            is ``<name>_ies_master`` beside the template, :meth:`run_ies`'s default.
            Unused by a serial resume.
        keep_workers : bool, optional
            Keep the worker directories after the resumed run finishes cleanly;
            default ``False`` deletes them, as :meth:`run_ies`'s ``keep_workers``.
            Unused by a serial resume.
        exe : str, optional
            PEST++ IES executable, default ``"pestpp-ies"``, resolved as in
            :meth:`run_ies`.
        overdue_giveup_fac : float or None, optional
            As in :meth:`run_ies`; default :data:`DEFAULT_OVERDUE_GIVEUP_FAC`,
            ``5.0``, and ``None`` removes the option. Set on the run's own control
            file, not on ``cal.pst``.
        **pestpp_options
            Further ``pestpp_options`` for the resumed run, set on the interrupted
            run's own control file (whose options -- ``ies_num_reals``, noise,
            ``ies_bad_phi_sigma``, ... -- are kept), not on ``cal.pst``. The four
            options the restart sets, ``ies_par_en``, ``ies_restart_par_en``,
            ``ies_restart_obs_en`` and ``ies_obs_en``, are refused (``ValueError``).

        Returns
        -------
        IesResults
            The joined run, opened on its directory with case ``name``.

        Raises
        ------
        ValueError
            When :meth:`build` has not run; when an option the restart sets is
            passed; when the directory holds no ``<name>.pst`` (no run to resume --
            the message names where a run of the other kind is, if one is); when the
            run was a prior Monte Carlo (``NOPTMAX -1``); when the run's parameters or
            observations differ from ``cal.pst``'s (the template was rebuilt since:
            the workers, made from it, would run another calibration); when the run
            completed no iteration (nothing to restart from: run again with
            :meth:`run_ies`); or when ``iterations`` does not exceed the last
            completed one.
        RuntimeError
            When a run may still be going in the template, or (parallel) in the
            master or its worker directory (ledger 532); or for a relative
            ``PYTHONPATH`` entry, as in :meth:`run_ies` (ledger 454).
        FileExistsError
            As :meth:`run_ies` raises it, for a master or worker directory that is a
            file or not a copy of this calibration's template.
        """

        if self.pst is None:
            raise ValueError("Call build() before resuming a PESTPP-IES run.")
        refused = sorted(set(pestpp_options) & set(self._RESUME_OPTIONS))
        if refused:
            raise ValueError(
                f"resume_ies() sets {', '.join(refused)} itself, to restart from the "
                "interrupted run's own ensembles; pass them to run_ies() to start a run "
                "from ensembles of your choosing."
            )
        parallel = bool(workers and int(workers) > 1)
        case = self.name
        directory = (
            (Path(master_dir) if master_dir else self._default_master("run_ies")).resolve()
            if parallel else self.template_workspace.resolve()
        )
        self._refuse_a_live_run(self.template_workspace, "template", method="resume_ies",
                                instead="", stacklevel=3)
        if not (directory / f"{case}.pst").is_file():
            raise ValueError(self._nothing_to_resume(directory, parallel, workers))
        if parallel:
            self._check_parallel_directories(
                directory, int(workers), method="resume_ies", default=not master_dir,
                overwrite=True, resume=True, stacklevel=3,
            )
        pending = pending_resume(directory, case)
        if pending is not None:
            joined = splice_resumed_run(directory, case, pending)
            self._restore_resumed_control_file(directory, pending)
            message = (
                f"resume_ies(): joined the iterations a previous resume left unjoined in "
                f"{pending} back into {directory} (it was stopped before it could): "
                f"iterations {joined}."
            )
            logger.warning(message)
            warnings.warn(message, UserWarning, stacklevel=2)

        pyemu = self.pyemu or _import_pyemu()
        with self._quiet_pyemu_context():
            run_pst = pyemu.Pst(str(directory / f"{case}.pst"))
        planned = int(run_pst.control_data.noptmax)
        if planned == -1:
            raise ValueError(
                f"The run in {directory} is a prior Monte Carlo (NOPTMAX -1), which "
                "PESTPP-IES evaluates in one sweep and cannot resume. Run it again: "
                f"cal.prior(...{', workers=' + str(workers) if parallel else ''}"
                f"{', overwrite=True' if parallel else ''})."
            )
        if (
            list(run_pst.par_names) != list(self.pst.par_names)
            or list(run_pst.obs_names) != list(self.pst.obs_names)
        ):
            raise ValueError(
                f"The run in {directory} calibrated other parameters or observations than "
                f"this calibration's control file ({run_pst.npar} parameters and "
                f"{run_pst.nobs} observations there against {self.pst.npar} and "
                f"{self.pst.nobs} here): the template was rebuilt since. The resumed run's "
                "workers are made from the template, so they would run a different "
                "calibration. Start it again: "
                f"cal.run_ies(...{', workers=' + str(workers) + ', overwrite=True' if parallel else ''})."
            )
        binary = str(run_pst.pestpp_options.get("ies_save_binary", "")).strip().lower() in (
            "true", "1", "yes", "t")
        last = last_completed_iteration(directory, case, binary_first=binary)
        if last is None:
            raise ValueError(
                f"The run in {directory} completed no iteration -- not even its prior "
                "(iteration 0) -- so there is nothing to restart from. Start it again: "
                f"cal.run_ies(...{', workers=' + str(workers) + ', overwrite=True' if parallel else ''})."
            )
        total = planned if iterations is None else int(iterations)
        if total <= last:
            raise ValueError(
                f"The run in {directory} has completed iteration {last}, and "
                f"iterations={total} asks for no more. Pass a larger total, e.g. "
                f"cal.resume_ies(iterations={last + 1}"
                f"{', workers=' + str(workers) if parallel else ''}), to go further."
            )

        archive = archive_for_resume(directory, case, last)
        # What the run's control file said before the restart rewrites it, put back after
        # the splice (and by the next resume after a hard kill), so the directory reads as
        # ONE run launched for `total` iterations: left as the restart wrote it, NOPTMAX
        # was the iterations LEFT, and a run interrupted again read as finished and its
        # default resume was refused (review of ledger 533).
        _record_control_file(archive, total, {
            key: run_pst.pestpp_options.get(key)
            for key in (*self._RESUME_OPTIONS, "ies_initial_lambda")
        })
        options = {
            "ies_par_en": ensemble_file(archive, f"{case}.0.par", binary_first=binary),
            "ies_restart_par_en": ensemble_file(archive, f"{case}.{last}.par", binary_first=binary),
            "ies_restart_obs_en": ensemble_file(archive, f"{case}.{last}.obs", binary_first=binary),
            "ies_obs_en": ensemble_file(archive, f"{case}.obs+noise", binary_first=binary),
        }
        for key, path in options.items():
            if path is None:
                run_pst.pestpp_options.pop(key, None)
            else:
                run_pst.pestpp_options[key] = f"{archive.name}/{path.name}"
        carried = lambda_in_force(archive / f"{case}.rec", last)
        if carried is not None and carried > 0:
            run_pst.pestpp_options["ies_initial_lambda"] = carried
        if overdue_giveup_fac is None:
            run_pst.pestpp_options.pop("overdue_giveup_fac", None)
        else:
            run_pst.pestpp_options["overdue_giveup_fac"] = float(overdue_giveup_fac)
        for key, value in pestpp_options.items():
            run_pst.pestpp_options[key] = value
        run_pst.control_data.noptmax = total - last
        logger.info(
            "resume_ies(): restarting %r in %s from iteration %d for %d more iteration(s)",
            case, directory, last, total - last,
        )
        exe_name, exe_dir = self._resolve_pestpp(exe)
        try:
            with self._quiet_pyemu_context():
                run_pst.write(str(directory / f"{case}.pst"), version=2)
            if parallel:
                self._run_parallel(
                    pyemu, exe_name, exe_dir, f"{case}.pst", directory, int(workers),
                    method="resume_ies", keep_workers=keep_workers,
                )
            else:
                self._run_serial(pyemu, exe_name, exe_dir, f"{case}.pst", method="resume_ies")
        finally:
            splice_resumed_run(directory, case, archive)
            self._restore_resumed_control_file(directory, archive)
        return self._review(directory)

    def _restore_resumed_control_file(self, directory: Path, archive: Path) -> None:
        """Put back what :meth:`resume_ies` changed in the run's control file (review of ledger 533).

        ``NOPTMAX`` becomes the run's TOTAL iterations and the restart options
        (``ies_par_en``, ``ies_restart_par_en``, ``ies_restart_obs_en``, ``ies_obs_en``,
        ``ies_initial_lambda``) what they were before the restart, as
        :func:`_record_control_file` recorded them in ``archive`` -- so the joined
        directory reads as one run launched for that many iterations. Nothing is done
        for an archive that recorded none.
        """

        try:
            record = json.loads((archive / RESUME_RECORD).read_text(encoding="utf-8"))
        except (OSError, ValueError):
            # No record, or not JSON: nothing recorded to put back.
            return
        control = record.get("control_file") if isinstance(record, dict) else None
        if not isinstance(control, dict):
            return
        pyemu = self.pyemu or _import_pyemu()
        path = directory / f"{self.name}.pst"
        with self._quiet_pyemu_context():
            run_pst = pyemu.Pst(str(path))
            run_pst.control_data.noptmax = int(control["noptmax"])
            for key, value in control.get("options", {}).items():
                if value is None:
                    run_pst.pestpp_options.pop(key, None)
                else:
                    run_pst.pestpp_options[key] = value
            run_pst.write(str(path), version=2)

    def _review(self, directory: Path):
        """The :class:`~myflopy.modflow.mf6.pest.ies.IesResults` of the run in ``directory``, case ``name``."""

        from myflopy.modflow.mf6.pest.ies import IesResults

        return IesResults(directory, case_name=self.name, model=self.model)

    def _nothing_to_resume(self, directory: Path, parallel: bool, workers) -> str:
        """The refusal for :meth:`resume_ies` when ``directory`` holds no run, naming where one is."""

        case = self.name
        message = f"There is no run of {case!r} to resume in {directory}: it holds no {case}.pst."
        if parallel and (self.template_workspace / f"{case}.phi.actual.csv").is_file():
            message += (
                f" A SERIAL run is in the template {self.template_workspace}; resume it "
                "there: cal.resume_ies()."
            )
        elif not parallel and (self._default_master("run_ies") / f"{case}.pst").is_file():
            message += (
                f" A PARALLEL run's master is {self._default_master('run_ies')}; resume it "
                "with workers: cal.resume_ies(workers=<N>)."
            )
        return message

    def _refuse_regularization_mode(self) -> None:
        """Refuse a control file PESTPP-IES cannot parse, before launching it.

        Tikhonov regularization is a PESTPP-GLM concept. Applying pyEMU's
        ``zero_order_tikhonov``/``first_order_pearson_tikhonov`` to ``cal.pst``
        flips ``pestmode`` to ``regularization`` and adds prior-information
        equations, and PESTPP-IES does two things with that: it prints
        ``prior information equations not supported in ensemble methods,
        ignoring``, and -- because myflopy writes version-2 control files -- the
        ``* regularization`` keywords land in the control-data keyword block and
        the run dies with ``control file parsing error`` before a single model
        run. Measured 2026-07-30 on PEST++ 5.2.16.

        So the failure is not "regularization is ignored", it is "the run does
        not start", from an error naming a keyword the caller never typed. IES
        has its own knob for this -- ``run_ies(ies_reg_factor=...)`` -- which
        needs none of the above.
        """

        mode = str(getattr(self.pst.control_data, "pestmode", "")).strip().lower()
        if not mode.startswith("regul"):
            return
        raise ValueError(
            "This control file is in 'regularization' mode, which PESTPP-IES "
            "cannot parse in the version-2 format myflopy writes -- the run "
            "would fail with a control-file parsing error before any model run. "
            "Tikhonov regularization is a PESTPP-GLM concept; ensemble methods "
            "ignore prior-information equations entirely.\n\n"
            "For IES, the equivalent knob is a regularization FACTOR on the "
            "ensemble update:\n"
            "    cal.run_ies(reals=..., ies_reg_factor=0.1)\n"
            "and spatial structure belongs in the prior itself -- pass "
            "`correlation=` to a grid or pilot-point parameter."
        )

    def _inject_geostatistical_prior(self, reals: int) -> None:
        """Draw a geostatistically-correlated prior parameter ensemble.

        Grid / pilot-point parameters carry a variogram (``correlation=``). A
        plain bounds-based prior would treat each cell independently, so prior
        realizations come out as per-cell spatial noise and IES cannot recover a
        coherent field. Drawing from the prior covariance instead yields *smooth*
        prior realizations (handed to PESTPP-IES via ``ies_par_en``) -- this is
        what makes posterior property-pattern maps meaningful rather than
        laughable. No-op when there are no geostatistical parameters. The launch
        calls it only when the call names no ``ies_par_en`` of its own (ledger
        509).

        The ensemble holds EVERY adjustable parameter of ``cal.pst`` -- PESTPP-IES
        refuses one that lacks any -- so it starts from :meth:`draw_prior` (each
        ``grid`` parameter correlated through its geostruct, every other parameter
        drawn on its own from its bounds), and each correlated pilot-point set's
        columns are then REPLACED by its own correlated draw (ledger 385).
        Before, a correlated ``grid`` parameter beside correlated pilot points
        overlapped on every pilot column and pandas refused the join, and correlated
        pilot points beside any other parameter (with no correlated ``grid`` one)
        handed PESTPP-IES an ensemble of the pilot points alone, which it refused:
        either way ``run_ies()``/``prior()`` died before a model ran.
        """

        # Grid parameters live in pyEMU's PstFrom, so `pf.draw` correlates them.
        # Pilot points do NOT: they are PEST template parameters registered after
        # `build_pst`, so they need their own draw (see
        # `_draw_pilot_point_prior`).
        correlated = self._correlated_prior_specs()
        if self.pf is None or not correlated:
            return

        with self._quiet_pyemu_context():
            ensemble = self.draw_prior(int(reals), use_specsim=False)
        if ensemble is None:
            return
        ensemble = ensemble._df if hasattr(ensemble, "_df") else ensemble
        for spec in correlated:
            if spec.style != "pilotpoints":
                continue
            drawn = self._draw_pilot_point_prior(spec, int(reals))
            if drawn is None:
                continue
            # Both draws number their realizations 0..n-1; anything else would pair
            # one realization's pilot points with another's other parameters.
            if not ensemble.index.equals(drawn.index):
                raise RuntimeError(
                    f"The prior draw's realizations {list(ensemble.index[:5])} and the "
                    f"pilot-point draw's {list(drawn.index[:5])} of {spec.name!r} are not "
                    "the same realizations; refusing to pair them."
                )
            ensemble = ensemble.drop(columns=ensemble.columns.intersection(drawn.columns))
            ensemble = ensemble.join(drawn)
        if ensemble.empty:
            return

        ensemble_path = self.template_workspace / f"{self.name}.prior_par.csv"
        ensemble.to_csv(ensemble_path)
        self.pst.pestpp_options["ies_par_en"] = ensemble_path.name

    def _draw_pilot_point_prior(self, spec, reals: int):
        """A geostatistically correlated prior for one pilot-point spec.

        ``correlation=`` was documented for ``style="pilotpoints"`` since the
        module was written (``pilot_points.py:14``) and read nowhere in it: the
        pilot VALUES were drawn independently from their bounds by PESTPP-IES,
        and the code's stated reasoning was that IDW interpolation supplies the
        spatial structure. It supplies smoothness, but not the CORRELATION LENGTH
        the caller asked for -- the field's structure came from pilot spacing
        instead, so a 200 m and a 2000 m variogram produced the same prior.

        Drawing the pilot values from the variogram's covariance fixes that. It
        matters most where the data are thin: with a correlated prior, IES moves
        neighbouring points together and unobserved areas keep a plausible field,
        rather than each point wandering on its own.

        The covariance is SCALED by the variance the bounds imply
        (:data:`PRIOR_SIGMA_RANGE`), through pyEMU's own
        ``helpers.geostatistical_draws`` -- the rule ``PstFrom.draw`` applies to a
        correlated ``grid`` parameter. The unscaled variogram drew every pilot value
        with a log10 standard deviation of 1.0, so with bounds ``(0.5, 2)`` 77% of the
        draw lay outside them and PESTPP-IES, which holds the initial ensemble to the
        bounds, ran each at a bound (ledger 424).
        """

        sets = [item for item in self._pilot_point_sets if item["spec"] is spec]
        if not sets:
            return None
        import pandas as pd

        pyemu = self.pyemu or _import_pyemu()
        geostruct = self._geostruct_for(spec)
        # One covariance per pilot-point SET -- a (layer, zone) of this call: the sets
        # are separate parameter groups interpolated separately, and correlating
        # across them would assert a vertical structure, or one across a unit
        # boundary, the variogram never described (ledger 527). pyEMU skips
        # a set of one point, which then keeps the draw from its bounds that
        # `draw_prior` gave it. pyEMU's own geostatistical draw,
        # the one `pf.draw` makes for a correlated grid parameter: the variogram
        # (sill 1) SCALED by the variance the bounds imply at PRIOR_SIGMA_RANGE.
        # Drawn from the unscaled variogram, a pilot value's log10 standard
        # deviation was 1.0 whatever its bounds -- 0.10 for a grid parameter of the
        # same bounds (0.5, 2) -- so 77% of the draw lay outside them and
        # PESTPP-IES ran it AT a bound (ledger 424).
        layers = [
            pd.DataFrame({
                "parnme": [str(name) for name in item["frame"]["parnme"]],
                "x": item["frame"]["x"].to_numpy(dtype=float),
                "y": item["frame"]["y"].to_numpy(dtype=float),
            })
            for item in sets
        ]
        names = pd.Index([name for layer in layers for name in layer["parnme"]])
        with self._quiet_pyemu_context(), contextlib.redirect_stdout(io.StringIO()):
            drawn = pyemu.helpers.geostatistical_draws(
                self.pst, struct_dict={geostruct: layers}, num_reals=int(reals),
                sigma_range=PRIOR_SIGMA_RANGE, scale_offset=True, subset=names,
                verbose=False,
            )
        values = drawn._df if hasattr(drawn, "_df") else drawn
        kept = [name for name in names if name in values.columns]
        if not kept:
            return None
        return values.loc[:, kept]

    def _default_master(self, method: str) -> Path:
        """The master directory ``method`` uses when no ``master_dir=`` is given."""

        suffix = "prior_master" if method == "prior" else "ies_master"
        return self.template_workspace.parent / f"{self.name}_{suffix}"

    def _refuse_a_live_run(
        self, directory: Path, what: str, *, method: str, instead: str, stacklevel: int
    ) -> None:
        """Refuse when a run is going in ``directory``; warn when its lock is stale.

        ``what`` names the directory's role in the message. A run is going when the
        :data:`~myflopy.modflow.mf6.pest.run_state.RUN_LOCK` names a process still
        running, a lock cannot be checked (another host), or a PEST++ process works
        there (ledger 532). Refused whatever ``overwrite=`` says.
        ``stacklevel`` puts the stale-lock warning at the user's call.

        Raises
        ------
        RuntimeError
            Naming the directory, what holds it and the spellings to use.
        """

        activity = run_activity(directory)
        if activity.blocks:
            raise RuntimeError(
                f"Refusing to start {method}(): the {what} {activity.directory} belongs to a "
                f"PESTPP-IES run that may still be going -- {activity.describe()} Deleting "
                "or reusing it would destroy that run's results or kill it.\n\n"
                + (f"Wait for it to finish, or {instead}" if instead else "Wait for it to finish.")
                + f"\nIf that run is over, delete {activity.directory / RUN_LOCK}."
            )
        if activity.state == "stale":
            message = (
                f"{method}(): the {what} {activity.directory}: {activity.describe()} "
                "Treating it as free."
            )
            logger.warning(message)
            warnings.warn(message, UserWarning, stacklevel=stacklevel)

    def _check_parallel_directories(
        self,
        master: Path,
        workers: int,
        *,
        method: str,
        default: bool,
        overwrite: bool,
        resume: bool = False,
        stacklevel: int = 4,
    ) -> None:
        """Refuse a parallel run whose master or worker directories are not this run's to use.

        The master is ``master``; the workers are ``<master>.workers/worker_0`` ..
        ``worker_<workers-1>`` (:func:`~myflopy.modflow.mf6.pest.run_state.worker_root`),
        so no two masters -- two calibrations, or one calibration's prior and history
        match -- share a worker directory (ledger 532). Until then every
        calibration of a model shared ``<model>.pest/worker_<n>`` (ledger 361).

        A launch deletes the master and each worker it is about to use, so each must be
        absent, empty, or a COPY of this calibration's template -- the template's
        :data:`BUILD_MARKER` names it, and a copy of a template built before ledger 360,
        which has none, counts by its metadata in the default layout only. A master that
        holds an earlier run's RESULTS is refused unless ``overwrite``
        (ledger 532: until then a relaunch deleted them, ledger 361); ``resume`` uses
        it as it is. Any of them, the worker root included, that belongs to a run still
        going is refused whatever ``overwrite`` says. Everything is checked before
        anything is written or started.

        Raises
        ------
        ValueError
            When the master or the worker root is, contains or sits inside the run
            directory or the template.
        RuntimeError
            When a run is going in the template, the master or the worker root.
        FileExistsError
            When one is a file; a non-empty directory that is not a copy of this
            calibration's template; or a master holding an earlier run's results
            without ``overwrite``.
        """

        run = self.original_workspace.resolve()
        template = self.template_workspace.resolve()
        master = Path(master).resolve()
        root = worker_root(master)
        workers_at = [root / f"worker_{index}" for index in range(workers)]
        default_master = self._default_master(method)
        call = f"cal.{method}(workers={workers}"
        if resume:
            instead = (
                "Resume the run where its master is:\n"
                f"    {call}, master_dir=<its master directory>)"
            )
        else:
            lines = ["Run it somewhere else:"]
            if not default:
                lines.append(f"    {call}, ...)   # the default master: {default_master}")
            lines.append(
                f"    {call}, master_dir=<a new or empty folder>, ...)   # its workers go "
                "in <that folder>.workers"
            )
            instead = "\n".join(lines)
        where = f"Refusing to start {method}(workers={workers}): the "

        for what, directory in (("master directory", master), ("worker directory", root)):
            for other, whose in (
                (run, "the model's run directory"),
                (template, "this calibration's template"),
            ):
                relation = _relation(directory, other)
                if relation is not None:
                    raise ValueError(
                        f"{where}{what} {directory} {relation} {whose}"
                        + ("" if relation == "IS" else f" {other}")
                        + ". A parallel run deletes its master and worker directories "
                        f"and writes into them.\n\n{instead}"
                    )

        beside = (
            f"run beside it in another folder:\n    {call}, master_dir=<a new or empty "
            "folder>, ...)"
        )
        self._refuse_a_live_run(master, "master directory", method=method, instead=beside,
                                stacklevel=stacklevel + 1)
        self._refuse_a_live_run(root, "worker directory", method=method, instead=beside,
                                stacklevel=stacklevel + 1)

        for what, directory in (("master directory", master), ("worker directory", root),
                                *(("worker directory", worker) for worker in workers_at)):
            if os.path.lexists(directory) and not directory.is_dir():
                raise FileExistsError(f"{where}{what} {directory} is a FILE.\n\n{instead}")

        for what, directory in (("master directory", master),
                                *(("worker directory", worker) for worker in workers_at)):
            if not directory.is_dir():
                continue
            entries = sorted(entry.name for entry in directory.iterdir())
            if not entries:
                continue
            copy = _written_by_a_build(directory, "template", naming=template) or (
                default
                and not os.path.lexists(directory / BUILD_MARKER)
                and self._metadata_names_this_calibration(directory)
            )
            if not copy:
                shown = ", ".join(entries[:5]) + (" ..." if len(entries) > 5 else "")
                raise FileExistsError(
                    f"{where}{what} {directory} already holds {len(entries)} "
                    f"entr{'y' if len(entries) == 1 else 'ies'} ({shown}), and it is not a "
                    f"copy of this calibration's template: its {BUILD_MARKER} would name "
                    f"{template}. A parallel run deletes its master and worker directories "
                    "before it copies the template into them, so it uses only one that is "
                    "absent, empty, or a copy of THIS template.\n\n"
                    + instead
                    + f"\nIf {directory} is an old master or worker directory you no longer "
                    "need, delete it yourself."
                )
            if what == "master directory" and not resume and not overwrite:
                self._refuse_overwriting_results(master, workers, method=method, default=default)

    def _refuse_overwriting_results(self, master: Path, workers: int, *, method: str, default: bool) -> None:
        """Refuse a master holding an earlier run's results (ledger 532).

        Raises
        ------
        FileExistsError
            When ``master`` holds any file
            :func:`~myflopy.modflow.mf6.pest.runs.pestpp_ies_outputs` lists, naming the
            folder, whether its run finished, and the spellings: ``overwrite=True``,
            another ``master_dir=``, and ``resume_ies`` for an interrupted IES run.
        """

        outputs = pestpp_ies_outputs(master, self.name)
        if not outputs:
            return
        done = recorded_iterations(master, self.name)
        planned = None
        try:
            for line in (master / f"{self.name}.pst").read_text(encoding="utf-8").splitlines():
                words = line.split()
                if len(words) >= 2 and words[0].lower() == "noptmax":
                    planned = int(float(words[1]))
                    break
        except (OSError, ValueError):
            planned = None
        master_arg = "" if default else f", master_dir={str(master)!r}"
        lines = [
            f"    cal.{method}(workers={workers}{master_arg}, overwrite=True, ...)   "
            "# delete those results and run again",
            f"    cal.{method}(workers={workers}, master_dir=<a new or empty folder>, ...)   "
            "# keep them",
        ]
        if planned is not None and planned > 0 and (not done or max(done) < planned):
            state = (
                f"an INTERRUPTED history match: iterations {done} of the {planned} it was "
                "launched for" if done else
                f"an INTERRUPTED history match that recorded no iteration of the {planned} "
                "it was launched for"
            )
            if done:
                lines.insert(0, f"    cal.resume_ies(workers={workers}{master_arg})   "
                                "# continue it from its last completed iteration")
        elif planned == -1:
            state = "a prior Monte Carlo's results"
        else:
            state = (
                f"a finished history match (iterations {done})" if done
                else "an earlier run's results"
            )
        raise FileExistsError(
            f"Refusing to start {method}(workers={workers}): the master directory "
            f"{master} holds {state} -- {len(outputs)} output file(s) of calibration "
            f"{self.name!r}, e.g. {outputs[0].name}. A run deletes its master directory "
            "before it starts, and a multi-day run's results are not deleted unasked.\n\n"
            + "\n".join(lines)
        )

    def _launch_pestpp_ies(
        self,
        *,
        reals: int,
        noptmax: int,
        workers: int | None,
        exe: str,
        master_suffix: str,
        master_dir: str | Path | None = None,
        noise: bool = True,
        bad_phi_sigma: float | None = None,
        overdue_giveup_fac: float | None = DEFAULT_OVERDUE_GIVEUP_FAC,
        overwrite: bool = False,
        keep_workers: bool = False,
        pestpp_options: dict | None = None,
    ) -> Path:
        """Configure options, write the pst, run PESTPP-IES, and return the results dir.

        The results directory holds THIS run's outputs only (ledger 507): a
        serial run first deletes the earlier run's outputs from the template
        (:func:`~myflopy.modflow.mf6.pest.runs.pestpp_ies_outputs`, bar a file a
        PEST++ option names), saying so (ledger 532); a parallel run makes its
        master itself, a copy of the template with those outputs removed, and hands it
        to pyEMU's ``start_workers(reuse_master=True)`` (ledger 532).
        """

        if self.pst is None:
            raise ValueError("Call build() before running PESTPP-IES.")
        method = "run_ies" if master_suffix == "ies_master" else "prior"
        parallel = bool(workers and int(workers) > 1)
        self._refuse_a_live_run(self.template_workspace, "template", method=method, instead="",
                                stacklevel=4)
        if parallel:
            master = (Path(master_dir) if master_dir else self._default_master(method)).resolve()
            # Before anything is written or started (ledger 361, ledger 532).
            self._check_parallel_directories(
                master,
                int(workers),
                method=method,
                default=not master_dir,
                overwrite=overwrite,
            )
        caller_options = dict(pestpp_options or {})
        case = f"{self.name}.pst"
        self.pst.pestpp_options["ies_num_reals"] = int(reals)
        if not noise:
            self.pst.pestpp_options["ies_no_noise"] = True
        if bad_phi_sigma is not None:
            self.pst.pestpp_options["ies_bad_phi_sigma"] = float(bad_phi_sigma)
        if overdue_giveup_fac is None:
            self.pst.pestpp_options.pop("overdue_giveup_fac", None)
        else:
            self.pst.pestpp_options["overdue_giveup_fac"] = float(overdue_giveup_fac)
        for key, value in caller_options.items():
            self.pst.pestpp_options[key] = value
        self._refuse_regularization_mode()
        # A caller's own `ies_par_en` is used as given; the correlated prior is drawn
        # only when the call names none (ledger 509).
        if "ies_par_en" not in caller_options:
            self._refuse_a_restart_from_other_parameters(method)
            self._inject_geostatistical_prior(int(reals))
        self.pst.control_data.noptmax = int(noptmax)

        if not parallel:
            # The run writes into the template, where an earlier run of this case left
            # its outputs; one it does not overwrite -- an iteration past this run's
            # last, a noise ensemble a `noise=False` run never writes -- would be
            # reviewed as this run's (ledger 507).
            self._clear_earlier_run_outputs(self.template_workspace, announce=not overwrite, method=method)
        with self._quiet_pyemu_context():
            self.pst.write(str(self.template_workspace / case), version=2)

        exe_name, exe_dir = self._resolve_pestpp(exe)
        pyemu = self.pyemu or _import_pyemu()
        if parallel:
            # The master is made HERE, not by pyEMU: a copy of the template without an
            # earlier serial run's outputs, which pyEMU's own copy carried into it
            # (ledger 532; ledger 507 removed them after the run instead).
            if master.exists():
                shutil.rmtree(master)
            shutil.copytree(
                self.template_workspace, master, symlinks=True,
                ignore=shutil.ignore_patterns(RUN_LOCK),
            )
            self._clear_earlier_run_outputs(master, announce=False, method=method)
            self._run_parallel(
                pyemu, exe_name, exe_dir, case, master, int(workers),
                method=method, keep_workers=keep_workers,
            )
            return master

        self._run_serial(pyemu, exe_name, exe_dir, case, method=method)
        return self.template_workspace

    def _run_serial(self, pyemu, exe_name: str, exe_dir, case: str, *, method: str) -> None:
        """Run PESTPP-IES in the template, holding its run lock while it goes."""

        write_run_lock(self.template_workspace, method=method, calibration=self.name)
        try:
            with self._pestpp_environment(exe_dir):
                pyemu.os_utils.run(f"{exe_name} {case}", cwd=str(self.template_workspace))
        finally:
            release_run_lock(self.template_workspace)

    def _run_parallel(
        self, pyemu, exe_name: str, exe_dir, case: str, master: Path, workers: int,
        *, method: str, keep_workers: bool,
    ) -> None:
        """Run PESTPP-IES in ``master`` with ``workers`` workers in ``<master>.workers``.

        Both directories carry a run lock while it goes. The workers are deleted when the
        run finishes cleanly, unless ``keep_workers``; a run that fails or is interrupted
        keeps them, for diagnosis (ledger 532). pyEMU's own clean-up
        (``cleanup=True``) deleted them whenever the master exited, failed or not.
        """

        root = worker_root(master)
        root.mkdir(exist_ok=True)
        write_run_lock(master, method=method, calibration=self.name)
        write_run_lock(root, method=method, calibration=self.name)
        finished = False
        try:
            with self._pestpp_environment(exe_dir):
                pyemu.os_utils.start_workers(
                    str(self.template_workspace),
                    exe_name,
                    case,
                    num_workers=workers,
                    worker_root=str(root),
                    master_dir=str(master),
                    reuse_master=True,
                    cleanup=False,
                )
            finished = True
        finally:
            release_run_lock(master)
            release_run_lock(root)
            if finished and not keep_workers:
                self._remove_workers(root)
            elif any(root.glob("worker_*")):
                logger.info(
                    "kept the worker directories of %s() in %s%s", method, root,
                    "" if finished else " (the run did not finish cleanly)",
                )

    def _remove_workers(self, root: Path) -> None:
        """Delete each ``worker_<n>`` in ``root`` that is a copy of this template, then ``root`` if empty."""

        template = self.template_workspace.resolve()
        removed = 0
        for worker in sorted(root.glob("worker_*")):
            if worker.is_dir() and _written_by_a_build(worker, "template", naming=template):
                shutil.rmtree(worker)
                removed += 1
        with contextlib.suppress(OSError):
            root.rmdir()   # refuses, rather than deletes, when anything else is in it
        logger.info("removed %d worker director%s of a finished run from %s",
                    removed, "y" if removed == 1 else "ies", root)

    def _clear_earlier_run_outputs(self, directory: Path, *, announce: bool = False,
                                   method: str = "run_ies") -> None:
        """Delete an earlier PESTPP-IES run's outputs of this case from ``directory``.

        Every file :func:`~myflopy.modflow.mf6.pest.runs.pestpp_ies_outputs` lists,
        and every ``<name>.resumed_from_<n>/`` archive a :meth:`resume_ies` made,
        except one a PEST++ option of this run names (``ies_par_en=<name>.0.par.csv``
        reuses an earlier prior from the template itself), which PESTPP-IES reads
        before it writes (ledger 507). With ``announce`` the deletion is said through
        ``warnings.warn`` as well as the logger (ledger 532; ledger 288: a log
        line alone reaches nobody) -- a serial run deleting the results an earlier serial
        run left in the template.
        """

        named = {
            Path(str(value)).name
            for value in self.pst.pestpp_options.values()
            if isinstance(value, (str, os.PathLike))
        }
        removed = []
        for path in pestpp_ies_outputs(directory, self.name):
            if path.name in named:
                continue
            path.unlink()
            removed.append(path.name)
        # And what a resume_ies of an earlier run left: its archive of that run's
        # .rec/.log -- or, when that resume was killed before it could join it back,
        # the earlier run's iterations, which the next resume_ies would splice into
        # THIS run (review of ledger 533). A parallel master is copied from
        # the template, so it would carry them too.
        for archive in sorted(directory.glob(f"{self.name}.resumed_from_*")):
            if archive.is_dir() and (archive / RESUME_RECORD).is_file():
                shutil.rmtree(archive)
                removed.append(archive.name + "/")
        if not removed:
            return
        shown = ", ".join(removed[:8]) + (" ..." if len(removed) > 8 else "")
        if announce:
            message = (
                f"{method}() deleted {len(removed)} output file(s) of an earlier PESTPP-IES "
                f"run of {self.name!r} from the template {directory} before this serial run "
                f"writes there ({shown}). A serial run's results last until the next serial "
                f"run; run in parallel (workers=N) to keep each run in its own master "
                f"directory, or pass overwrite=True to say this deletion is meant."
            )
            logger.warning(message)
            warnings.warn(message, UserWarning, stacklevel=4)
        else:
            logger.info(
                "removed %d output file(s) of an earlier PESTPP-IES run of %r from %s "
                "before this run writes there: %s", len(removed), self.name, directory, shown,
            )

    def _refuse_a_restart_from_other_parameters(self, method: str) -> None:
        """Refuse a restart observation ensemble the correlated prior draw would orphan.

        ``ies_restart_obs_en`` hands PESTPP-IES the simulated values of the parameter
        ensemble ``ies_par_en`` names, realization by realization; it skips running
        them. When a ``grid`` or array ``pilotpoints`` parameter carries
        ``correlation=``, a launch that names no ``ies_par_en`` DRAWS a new one, so the
        restart ensemble would be paired with parameters it was never run with -- the
        realization names, 0..n-1 in both, match, so nothing failed (ledger
        509).

        Raises
        ------
        ValueError
            When ``ies_restart_obs_en`` is set on ``cal.pst`` (by this call or an
            earlier one) and a correlated prior would be drawn.
        """

        restart = self.pst.pestpp_options.get("ies_restart_obs_en")
        if not restart or self.pf is None or not self._correlated_prior_specs():
            return
        names = ", ".join(repr(spec.name) for spec in self._correlated_prior_specs())
        raise ValueError(
            f"Refusing to start {method}(): ies_restart_obs_en={restart!r} is set, so "
            "PESTPP-IES would take those simulated values as the runs of the parameter "
            "ensemble ies_par_en names -- but this call names no ies_par_en, and the "
            f"correlated parameter(s) {names} would have a NEW prior ensemble drawn "
            "into ies_par_en. The observations would belong to other parameters, "
            "realization by realization.\n\n"
            "Name the parameter ensemble those observations were run with:\n"
            f"    cal.{method}(ies_par_en=<its .par.csv>, ies_restart_obs_en={restart!r}, ...)\n"
            "or, to draw a new prior and run it, drop the restart ensemble (it persists "
            "on cal.pst from the call that set it):\n"
            "    del cal.pst.pestpp_options['ies_restart_obs_en']"
        )

    def _correlated_prior_specs(self) -> list:
        """The parameter specs whose prior :meth:`_inject_geostatistical_prior` draws correlated."""

        return [
            spec
            for spec in self._native_parameter_specs
            if spec.style in ("grid", "pilotpoints") and spec.correlation is not None
        ]
