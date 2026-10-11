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
from pathlib import Path

import flopy

from myflopy._logging import get_logger
from myflopy._optional import require
from myflopy.modflow.mf6.lak_implicit import keep_implicit
from myflopy.modflow.mf6.observations import (
    ConcTargets,
    DrnFlowTargets,
    HeadTargets,
    LakeStageTargets,
    SfrFlowTargets,
    SfrStageTargets,
)
from myflopy.modflow.mf6.period_end import period_lengths
from myflopy.modflow.mf6.pest.geostats import build_geostruct
from myflopy.modflow.mf6.pest.native_parameters import (
    NativeParameterSpec,
    add_native_parameter,
    relayer_array_target,
    staged_model,
)
from myflopy.modflow.mf6.pest.observations import (
    finalize_observations,
    prepare_conc_observations,
    prepare_drn_flow_observations,
    prepare_head_target_observations,
    prepare_lake_stage_observations,
    prepare_sfr_flow_observations,
    prepare_sfr_stage_observations,
)
from myflopy.modflow.mf6.pest.pilot_points import (
    add_pilot_point_parameter,
    register_pilot_point_parameters,
)
from myflopy.modflow.mf6.pest.runs import default_pest_root
from myflopy.modflow.mf6.pest.specs import (
    ConcObservationSpec,
    DrnFlowObservationSpec,
    ExpGeoStruct,
    HeadTargetObservationSpec,
    LakeStageObservationSpec,
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

#: The file a build writes into its template, and into its staging copy, naming the
#: directory it was written for and that directory's role: the one sign by which a
#: later build may DELETE that directory (ledger 360). A build deletes a directory
#: only when it is absent, empty, or carries this file naming it -- a folder the build
#: did not write is refused, whatever else it holds. VISIBLE, not a dot-file (ledger
#: 361): a hidden marker survived ``rm <template>/*`` and Finder's select-all delete,
#: so a folder emptied that way and refilled with the user's files still carried it,
#: and the next build deleted them.
BUILD_MARKER = "myflopy_pest_build.json"

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
    LakeStageObservationSpec,
    SfrStageObservationSpec,
    SfrFlowObservationSpec,
    DrnFlowObservationSpec,
)

# High-level target sets and the observation spec each is wrapped in, so
# observe()/forecast() accept head, lake-stage, SFR stage/flow and DRN seepage
# targets uniformly (each spec carries its own default prefix).
_TARGET_SPEC_TYPES = (
    # ConcTargets subclasses HeadTargets, so it MUST come first -- an isinstance
    # ladder would otherwise wrap every concentration target as a head spec.
    (ConcTargets, ConcObservationSpec),
    (HeadTargets, HeadTargetObservationSpec),
    (LakeStageTargets, LakeStageObservationSpec),
    (SfrStageTargets, SfrStageObservationSpec),
    (SfrFlowTargets, SfrFlowObservationSpec),
    (DrnFlowTargets, DrnFlowObservationSpec),
)


def _pest_run_slug(name: str) -> str:
    """Filesystem-safe directory name for a PEST run's default workspace."""

    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(name)).strip("_")
    return slug or "pest_run"


#: Name-file blocks whose entries name an INPUT file in their second field: the
#: simulation name file's ``timing``/``models``/``exchanges``/``solutiongroup`` and a
#: model name file's ``packages``. Never ``options``, whose ``LIST`` names an OUTPUT.
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
    (read in place from anywhere), outside ``sim_ws``, or not a file there.

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

    root = Path(sim_ws)
    found: dict[Path, None] = {}
    # (relative name, how the file is read): "simulation" / "model" name files
    # are read for their blocks, "input" files scanned for FILEIN/OPEN/CLOSE,
    # "data" files (OPEN/CLOSE) copied unread.
    pending: list[tuple[str, str]] = [("mfsim.nam", "simulation")]
    while pending:
        name, kind = pending.pop(0)
        relative = Path(os.path.normpath(name))
        if relative.is_absolute() or relative.parts[:1] == ("..",) or relative in found:
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
                elif (
                    block in _NAME_FILE_BLOCKS
                    and len(fields) > 1
                    and _MF6_FILE_TYPE.match(fields[0])
                ):
                    pending.append((fields[1], "model" if block == "models" else "input"))
        for keyword, token in _MF6_INPUT_REFERENCE.findall(text):
            pending.append((_mf6_field(token), "input" if keyword.upper() == "FILEIN" else "data"))
    return list(found)


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
        cal.parameterize("k",        style="pilotpoints", pp_space=8, physical=(1e-3, 100))
        cal.parameterize("recharge", bounds=(0.5, 1.5))
        cal.parameterize("ghb.cond", bounds=(0.1, 10))
        cal.observe(head_targets)          # history-matching targets
        cal.forecast(prediction_targets)   # predictions of interest (zero weight)
        print(cal.settings())              # review the resolved configuration
        pst = cal.build("calib.pst")       # -> native .pst + forward_run.py
        cal.run_ies()                      # PESTPP-IES ensemble smoother

    Use :meth:`parameterize` (styles ``constant`` / ``zone`` / ``grid`` /
    ``pilotpoints`` across every target), :meth:`observe`, :meth:`forecast`,
    :meth:`build`, :meth:`run_ies`/:meth:`prior` and :meth:`settings`. The runs
    land beside the model and are reopened for review via ``model.pest_runs``
    (or ``run.pest_runs``) -> :meth:`~...runs.PestRunHandle.review`
    (:class:`~...ies.IesResults`).

    Parameters
    ----------
    model : SimulationBase
        A built ``myflopy`` model (``SimulationBase`` or compatible) exposing
        ``.sim`` (the FloPy simulation), ``.gwf``, ``.name``,
        ``.workspace`` and, for ``pp_space=`` pilot points, any ``zones=`` and
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
        (``<name>_ies_master`` / ``<name>_prior_master``), and the run as
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
        and left exactly as it was. Everything is checked before anything is
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
        self._parameter_specs: list = []
        self._observation_specs: list = []
        self._forecast_specs: list = []
        self._native_parameter_specs: list[NativeParameterSpec] = []
        self._native_parameter_frames: dict[str, object] = {}
        self._pilot_point_frames: dict[str, list] = {}
        self._capture_field_specs: dict[str, dict] = {}
        self._prepared_observations: list[dict] = []
        self._prepared_parameters: dict[str, dict] = {}
        self._parameter_frames: dict[str, object] = {}
        self._registered_geostructs: dict[str, object] = {}
        self._forward_run_config: dict = {}
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
        bounds: tuple[float, float] = (0.5, 2.0),
        physical: tuple[float, float] | None = None,
        transform: str = "log",
        additive: bool | None = None,
        zones=None,
        layers=None,
        pp_space: int | None = None,
        pp_points=None,
        correlation: float | None = None,
        anisotropy: float = 1.0,
        bearing: float = 0.0,
        nugget: float = 0.0,
        temporal: float | None = None,
        name: str | None = None,
        capture: bool = False,
        **extra,
    ) -> NativeParameterSpec:
        """Declare a calibration parameter that compiles to ``pf.add_parameters``.

        Each call expands, at :meth:`build` time, into the correct
        ``pyemu.utils.PstFrom.add_parameters`` invocation against the external
        MODFLOW input file(s) for ``target``. Call it once per property you want
        to adjust; call it multiple times for the same target to stack styles
        (e.g. a broad ``"constant"`` plus a fine ``"grid"``).

        Parameters
        ----------
        target : str
            Friendly name of the model property to calibrate (case and
            surrounding spaces ignored). Supported targets
            (with aliases): ``"k"`` (``"npf.k"``/``"hk"``/``"kh"``),
            ``"k33"`` (``"kv"``), ``"recharge"`` (``"rch"``),
            ``"chd"``, ``"ghb.cond"``/``"ghb.bhead"`` (``"ghb"``),
            ``"drn.cond"``/``"drn.elev"`` (``"drn"``), ``"wel"`` (``"pumping"``),
            ``"uzf.vks"`` (``"vks"``) — UZF saturated vertical K, the one UZF
            quantity heads respond to — and ``"porosity"``
            (``"mst.porosity"``/``"n"``), which lives on the GWT sibling of a
            coupled transport simulation and is resolved there automatically.
            ``k``, ``k33`` and ``porosity`` are ARRAY targets; the rest scale a
            column of a LIST package's files. An unknown name raises
            ``KeyError`` here, listing the known ones. Each target scales the
            external files FloPy writes for it (``recharge``, for instance,
            reads ``<model>.rch_stress_period_data_*.txt``), and :meth:`build`
            raises ``FileNotFoundError`` when the model has none.

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
            - ``"pilotpoints"`` -- geostatistical pilot points.

            ``None`` uses the recipe default (``"constant"`` for every target).
            The value is lower-cased and not checked against these four: it
            reaches pyEMU as ``par_type``, except ``"pilotpoints"`` on an array
            target, which myflopy builds itself. ``"grid"`` on
            *array* targets (K, K33) builds one multiplier per Voronoi/DISV
            cell, and an independent set per layer file, located through the
            model's spatial reference (the modelgrid is used automatically);
            pass ``correlation=`` to correlate them and ``layers=`` to restrict
            which layers. ``"grid"`` also works for list targets.
            ``"pilotpoints"`` on an array target places a pilot-point net
            (``pp_space=`` spacing or explicit ``pp_points=``) and interpolates
            the multipliers to cells by inverse-distance weighting in the
            forward run -- the supported sparse-field style on Voronoi/DISV
            grids.
        bounds : tuple of (float, float), optional
            ``(lower, upper)`` bounds on the adjustable value, pyEMU's
            ``lower_bound``/``upper_bound`` (each pilot point's ``parlbnd`` /
            ``parubnd``). Default ``(0.5, 2.0)``. For a multiplier
            these are factors (e.g. ``(0.2, 5)`` = 5x down to 5x up); for an
            additive parameter they are offsets in model units.
        physical : tuple of (float, float) or None, optional
            ``(ult_lbound, ult_ubound)`` -- hard limits clamped onto the *final*
            model value after all multipliers are applied (pyEMU's
            ``ult_lbound``/``ult_ubound``; the pilot-point interpolation clamps
            to them itself). Strongly recommended
            for multipliers so calibration cannot push K (etc.) to nonphysical
            values. ``None`` (default) leaves the input unclamped.
        transform : str, optional
            ``"log"`` (default, recommended for strictly positive properties
            like K and recharge) or ``"none"``; pyEMU also accepts
            ``"fixed"``, but not together with ``correlation=`` (the variogram
            conversion asserts ``"none"`` or ``"log"``). Forced to ``"none"``
            for additive parameters.
        additive : bool or None, optional
            Apply the parameter as an additive offset instead of a multiplier
            (pyEMU ``par_style="a"``). ``None`` (default) uses the recipe
            default (drain elevation is additive; everything else is
            multiplicative). Array-target pilot points still MULTIPLY in their
            forward run; for them it only forces the transform to ``"none"``.
        zones : ZoneSpec or array_like of int or None, optional
            Per-cell integer zone labels, normalized to the shape pyEMU wants
            for the target's family at :meth:`build` (see :class:`ZoneSpec`,
            which also reads a raster or polygon layer) and passed as pyEMU's
            ``zone_array``. With ``style="zone"`` each distinct label gets one
            parameter. With any style, an array target leaves cells labelled
            ``<= 0`` unparameterized; a list target makes zone ``0`` a real
            parameter and warns. A length that is not one label per cell raises
            ``ValueError`` at build. Default ``None`` passes no zones.
            Array-target pilot points ignore it.
        layers : iterable of int or None, optional
            For multi-layer *array* targets, restrict the parameter to these
            zero-based model layers (e.g. ``layers=[0, 1]`` for the unconfined
            aquifers); other layers keep their input value. ``None`` (default)
            covers all layers. Read only for array targets -- a list target
            ignores it. :meth:`build` raises ``ValueError`` when no requested
            layer has its own file, or the target is written as one whole-grid
            file.
        pp_space : int or None, optional
            Pilot-point spacing for ``style="pilotpoints"`` on an array target,
            in multiples of the median cell width (the square root of the median
            ``model.vor`` cell area): the grid is binned at that spacing and the
            cell centre nearest each occupied bin's middle becomes a point.
            Truncated with ``int``; ``None`` (default) or ``0`` means ``5``.
            Ignored when ``pp_points`` is given.
        pp_points : geopandas.GeoDataFrame or array_like or None, optional
            Explicit pilot-point locations for ``style="pilotpoints"`` on an
            array target: a GeoDataFrame of Point geometries, or ``(x, y)``
            pairs (anything reshapable to ``(n, 2)``), in the model grid's own
            coordinates -- a GeoDataFrame's CRS is not read or reprojected.
            Takes precedence over ``pp_space``. Default ``None``.
        correlation : float or None, optional
            Variogram range (model length units) of an exponential variogram
            (see :class:`~myflopy.modflow.mf6.pest.specs.ExpGeoStruct`). With
            ``style="grid"`` it becomes the parameters' pyEMU ``geostruct``, and
            :meth:`run_ies`/:meth:`prior` draw a spatially correlated prior
            ensemble from it (``ies_par_en``). With ``style="pilotpoints"`` on an
            array target it shapes only that prior draw -- the pilot values are
            drawn correlated, while the interpolation to cells stays
            inverse-distance. ``None`` (default) builds no variogram; ignored
            for ``"constant"``/``"zone"``.
        anisotropy : float, optional
            Anisotropy ratio of the variogram -- correlated this many
            times farther along the major axis than across it (``1.0`` =
            isotropic, the default). For an alluvial valley, e.g. ``5`` makes K
            correlate farther down-valley than across it. Acts only when
            ``correlation`` is given; pyEMU asserts it is positive.
        bearing : float, optional
            Azimuth of the anisotropy major axis, in degrees east of north.
            Default ``0.0``. Ignored when ``anisotropy`` is ``1.0``, and
            without ``correlation``.
        nugget : float, optional
            Nugget (unresolved short-scale variance) of the variogram, added to
            the covariance diagonal on top of the variogram's contribution,
            which is fixed at ``1.0`` here (so the total sill is
            ``1.0 + nugget``). Default ``0.0``. Acts only when ``correlation``
            is given.
        temporal : float or None, optional
            Temporal correlation range (days) for time-varying list packages.
            Recorded on the returned spec and read by nothing yet, so it has no
            effect on the build. Default ``None``.
        name : str or None, optional
            Parameter group / name base (pyEMU's ``par_name_base`` and
            ``pargp``). Default ``None`` uses the canonical target with dots
            removed (``"k"``, ``"ghbcond"``, ``"uzfvks"``). A multi-layer
            array ``"grid"`` parameter appends ``l1``, ``l2``, ... per selected
            layer file, pilot points append ``l<layer>`` (zero-based), and a
            captured field is named ``<name>field``. Give each call a distinct
            name when stacking calls on one target: the project's
            per-parameter records are keyed by it, so a second call under the
            same name overwrites the first's.
        capture : bool, optional
            Also record the *resolved* model-input field (after multipliers are
            applied) as zero-weight observations named from ``<name>field``
            (one prefix per layer file for a layered array), so each
            realization's per-cell property field is carried in the ensembles.
            Enables spatial parameter maps via ``IesResults.plot_field`` /
            ``IesResults.field``. Default ``False``.
        **extra
            Advanced keyword arguments passed straight through to pyEMU's
            ``PstFrom.add_parameters`` for cases the facade does not cover.
            They are applied LAST, so a key here overrides what the facade
            derived (``par_type``, the bounds, ``transform``, ``zone_array``,
            ``geostruct``, ...). Not forwarded for array-target pilot points,
            which do not go through ``add_parameters``. Stored on the returned
            spec as ``extra``.

        Returns
        -------
        NativeParameterSpec
            The recorded specification (also stored on the project). Inspect the
            full set at any time with :meth:`settings`.

        Examples
        --------
        >>> cal.parameterize("k", style="constant", bounds=(0.2, 5), physical=(1e-3, 100))
        >>> cal.parameterize("recharge", bounds=(0.5, 1.5), physical=(0, 1e-2))
        >>> cal.parameterize("ghb.cond", bounds=(0.1, 10))
        >>> cal.parameterize("drn.elev", bounds=(-2, 2))   # additive offsets
        """

        spec = NativeParameterSpec(
            target=target,
            style=style,
            bounds=bounds,
            physical=physical,
            transform=transform,
            additive=additive,
            zones=zones,
            layers=tuple(layers) if layers is not None else None,
            pp_space=pp_space,
            pp_points=pp_points,
            correlation=correlation,
            anisotropy=anisotropy,
            bearing=bearing,
            nugget=nugget,
            temporal=temporal,
            name=name,
            capture=capture,
            extra=extra,
        )
        self._native_parameter_specs.append(spec)
        return spec

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

        Parameters
        ----------
        targets : target set or observation spec
            A high-level target set -- :class:`HeadTargets`,
            :class:`ConcTargets`, :class:`LakeStageTargets`,
            :class:`SfrStageTargets`, :class:`SfrFlowTargets` or
            :class:`DrnFlowTargets` -- which is wrapped in its matching
            ``*ObservationSpec``, or a pre-built observation spec (e.g.
            :class:`HeadTargetObservationSpec`), registered as it is. The
            matching simulated values are regenerated by a post-processor on
            every run. Anything else raises ``TypeError``.
        prefix : str or None, optional
            Observation-name prefix. ``None`` (default) uses the target kind's
            own prefix (``"hds"``, ``"conc"``, ``"lak_stage"``, ``"sfr_stage"``,
            ``"sfr_flow"``, ``"drn_flow"``). Pass distinct prefixes to register
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
            :class:`DrnFlowTargets`) or a pre-built observation spec,
            describing *where* and *when* the prediction is taken. Any values
            supplied are placeholders (e.g. an expected value), written as the
            observation's value; a blank keeps the base run's. The weight is
            forced to zero regardless. Anything else raises ``TypeError``.
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
        self._forecast_specs.append(spec)
        return spec

    def _coerce_observation_spec(self, targets, *, prefix):
        """Turn a target set (or spec) into a registered observation spec.

        Accepts a pre-built observation spec, or any high-level target set
        (head / lake stage / SFR stage / SFR flow / DRN seepage), wrapping it in
        the matching spec. When ``prefix`` is ``None`` the spec's own default
        prefix is used (``"hds"``, ``"lak_stage"``, ``"sfr_flow"``, ...); pass distinct
        prefixes to register more than one set of the same kind.
        """

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

        Raises
        ------
        ValueError
            From :meth:`_refuse_a_workspace_over_the_run`.
        FileExistsError
            From :meth:`_clearance`: a directory no build wrote is in the way.
        """

        self._refuse_a_workspace_over_the_run()
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
            spatial_reference = getattr(staged_model(self, self.model.name), "modelgrid", None)

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
        template and every IES worker. The inputs are the run's FILES, which is what
        its output files -- the base observations PEST reads -- were computed from.
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
        copied = time.perf_counter()
        simulation = flopy.mf6.MFSimulation.load(sim_ws=str(staging), verbosity_level=0)
        keep_implicit(simulation)
        self._staged_simulation = simulation
        # Before externalizing, not after: a griddata array supplied as a scalar
        # writes one whole-grid file that pyEMU cannot index against a DISV
        # spatial reference. See `relayer_array_target`.
        for spec in self._native_parameter_specs:
            relayer_array_target(self, spec.recipe)
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

    def _write_project_metadata(self, *, filename: str | Path | None = None):
        """Persist reopen-friendly run metadata inside the generated workspace."""

        observation_sets = [
            item["metadata"]
            for item in self._prepared_observations
            if "metadata" in item and item["metadata"] is not None
        ]
        parameter_sets = []
        for spec in self._parameter_specs:
            prepared = self._prepared_parameters.get(spec.name, {})
            config = prepared.get("config", {})
            parameter_sets.append(
                {
                    "name": spec.name,
                    "type": type(spec).__name__,
                    "kind": config.get("kind"),
                    "transform": getattr(spec, "transform", None),
                    "bounds_mode": getattr(spec, "bounds_mode", None),
                    "bounds": list(spec.bounds) if getattr(spec, "bounds", None) is not None else None,
                    "parameter_space": getattr(spec, "parameter_space", None),
                    "parameter_style": getattr(spec, "parameter_style", None),
                    "files": {key: value for key, value in config.items() if key.endswith("_csv")},
                }
            )

        metadata = {
            "version": 1,
            "project_name": self.name,
            "model_name": self.model.name,
            "original_workspace": str(self.original_workspace.resolve()),
            "pst_file": Path(filename).name if filename is not None else None,
            "observation_sets": observation_sets,
            "parameter_sets": parameter_sets,
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
        return self.pf.draw(num_reals=num_reals, use_specsim=bool(use_specsim))

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
        """

        return build_geostruct(
            ExpGeoStruct(
                range=float(spec.correlation),
                anisotropy=float(spec.anisotropy),
                bearing=float(spec.bearing),
                nugget=float(spec.nugget),
                transform=spec.resolved_transform,
            )
        )

    def _attach_native_observation_postprocessors(self):
        """Add post-model functions that regenerate simulated observation CSVs."""

        helper_path = Path(__file__).with_name("forward_run.py")
        # Each period's length, so the forward run can check that the record it
        # reads ENDS its period, as the base run's build did (ledger 336). OC is
        # not a parameter, but ATS can move a STEPS save off a period's end in one
        # realization and not another.
        perlen = period_lengths(self.model)
        for item in self._prepared_observations:
            head_config = item.get("forward_run_config")
            if head_config is not None:
                self.pf.add_py_function(
                    str(helper_path),
                    "_write_head_target_csv("
                    f"model_name='{self.model.name}', "
                    f"mapping_csv='{head_config['mapping_csv']}', "
                    f"output_csv='{head_config['output_csv']}', "
                    f"perlen={perlen!r})",
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
                    f"perlen={perlen!r})",
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
                # The calibrated model by NAME: without it the forward run reloads
                # the simulation's first GWF model, another model's lake or zone
                # in a simulation of several (ledger 340).
                self.pf.add_py_function(
                    str(helper_path),
                    "write_named_series_targets("
                    f"kind='{series_config['kind']}', "
                    f"locations_file='{locations_file}', "
                    f"output_csv='{series_config['output_csv']}', "
                    f"model_name='{self.model.name}', "
                    f"package={series_config.get('package')!r})",
                    is_pre_cmd=False,
                )
                continue
            raise NotImplementedError(
                f"Observation set {item.get('prefix')!r} has no forward-run "
                "post-processor; supported kinds are head-target and named-series "
                "(lake stage, SFR stage/flow, DRN seepage) observations."
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
        for filename in spec.resolved_files:
            if recipe.family == "array":
                match = re.search(r"_layer(\d+)\.txt$", str(filename))
                prefix = f"{base}l{int(match.group(1)) - 1}" if match else base
                if match:
                    layer_prefixes[int(match.group(1)) - 1] = prefix
                self.pf.add_observations(filename, prefix=prefix)
            else:
                self.pf.add_observations(
                    filename, prefix=base, index_cols=[0, 1], use_cols=[recipe.use_col]
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
            sits inside the model's run directory, whatever it holds (ledger 360).
        FileNotFoundError
            When the model's workspace holds no ``mfsim.nam``: the template is
            made from the run's files, and a model never written has none.
        FileExistsError
            When ``workspace=``, or the staging directory beside it, is a file or
            a non-empty folder no build wrote (ledger 360): the build deletes both
            directories first, so it refuses rather than delete what it did not
            write, before anything is deleted. Use the default, a new or empty
            folder, or delete it yourself.
        """

        target_name = Path(filename).name if filename else f"{self.name}.pst"
        self._refuse_shared_prefixes()
        self._release_pyemu_log()
        self._clear_stale_template()
        staging = self._staging_workspace()
        try:
            self._stage_model()
            self._build_pstfrom(staging)
            self._prepare_observation_specs()
            with self._quiet_pyemu_context():
                for spec in self._native_parameter_specs:
                    if spec.recipe.family == "array" and spec.style == "pilotpoints":
                        add_pilot_point_parameter(self, spec)
                    else:
                        add_native_parameter(self, spec)
                    if spec.capture:
                        self._add_capture_field_observations(spec)
            self.pf.mod_sys_cmds.append(self._resolve_exe())
            self._drop_unused_multiplier_apply()
            self._attach_native_observation_postprocessors()
            self._write_project_metadata(filename=target_name)
            with self._quiet_pyemu_context():
                self.pst = self.pf.build_pst(filename=target_name)
            if self._pilot_point_frames:
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
        bounds, physical limits, and transform; every observation and forecast)
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
                "style": spec.style,
                "bounds": tuple(spec.bounds),
                "physical": tuple(spec.physical) if spec.physical is not None else None,
                "transform": spec.resolved_transform,
                "additive": bool(spec.additive),
                "name": spec.name,
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
    def _augmented_path(self, exe_dir):
        """Temporarily prepend ``exe_dir`` to ``PATH`` so PEST++ tools resolve."""

        import os

        original = os.environ.get("PATH", "")
        if exe_dir and exe_dir not in original:
            os.environ["PATH"] = exe_dir + os.pathsep + original
        try:
            yield
        finally:
            os.environ["PATH"] = original

    def run_ies(
        self,
        *,
        reals: int = 50,
        iterations: int = 3,
        workers: int | None = None,
        noise: bool = True,
        bad_phi_sigma: float | None = None,
        master_dir: str | Path | None = None,
        exe: str = "pestpp-ies",
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
        valid first call; tune from there.

        ``noise``, ``bad_phi_sigma`` and ``**pestpp_options`` are written onto
        ``cal.pst.pestpp_options`` and STAY there: ``noise=False``, a
        ``bad_phi_sigma`` or an extra option from an earlier
        :meth:`run_ies`/:meth:`prior` call is still set on the next one, whose
        defaults do not clear it.

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
            and collects the results in ``master_dir``.
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
            Master directory for a parallel run; ``None`` (default) uses
            ``<name>_ies_master`` beside the template workspace, which is where
            ``model.pest_runs`` looks. The workers are created in its parent
            directory, as ``worker_0`` .. ``worker_<workers-1>``. Unused by a
            serial run. **pyEMU deletes the master and worker directories before
            it copies the template into them** (ledger 361), so each must be
            absent, empty, or a copy of THIS calibration's template -- a master or
            worker a previous run of it left, deleted with its results; anything
            else is refused before pyEMU starts (``FileExistsError``), and a
            master or worker that is, contains or sits inside the run directory
            or the template is refused whatever it holds (``ValueError``).
        exe : str, optional
            PEST++ IES executable name, default ``"pestpp-ies"`` (a trailing
            ``.exe`` is dropped). When the model's MF6 executable resolves to an
            absolute path and a file ``<exe>.exe`` sits beside it, that
            directory is prepended to ``PATH`` for the run (unless ``PATH``
            already contains it); otherwise the name is resolved on ``PATH`` as
            it is.
        **pestpp_options
            Any additional ``pst.pestpp_options`` to set (e.g.
            ``ies_localizer="loc.mat"``), passed straight through and applied
            after the options above, so they override them. ``ies_par_en`` is
            then replaced by the drawn prior when a ``grid`` parameter, or an
            array-target ``pilotpoints`` parameter, carries ``correlation=``.

        Returns
        -------
        IesResults
            Reader/visualizer for the completed run, opened on the directory
            holding its results with case ``name``.

        Raises
        ------
        ValueError
            When :meth:`build` has not run, or ``cal.pst`` has been put in
            ``regularization`` mode, which PESTPP-IES cannot parse; or, for a
            parallel run, when the master or a worker directory is, contains or
            sits inside the run directory or the template, or the master is named
            as one of its own workers (ledger 361).
        FileExistsError
            For a parallel run, before pyEMU starts, when the master or a worker
            directory is a file or a non-empty folder that is not a copy of this
            calibration's template (ledger 361): pyEMU would delete it.
        """

        from myflopy.modflow.mf6.pest.ies import IesResults

        results_dir = self._launch_pestpp_ies(
            reals=reals,
            noptmax=iterations,
            workers=workers,
            exe=exe,
            master_suffix="ies_master",
            master_dir=master_dir,
            noise=noise,
            bad_phi_sigma=bad_phi_sigma,
            pestpp_options=pestpp_options,
        )
        return IesResults(results_dir, case_name=self.name, model=self.model)

    def prior(
        self,
        *,
        reals: int = 50,
        workers: int | None = None,
        noise: bool = True,
        master_dir: str | Path | None = None,
        exe: str = "pestpp-ies",
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
        ensemble and stop). Call :meth:`build` first.

        Options persist on ``cal.pst`` between calls exactly as in
        :meth:`run_ies` (a ``bad_phi_sigma`` set there is still set here).

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
            ``model.pest_runs`` looks; the workers are created in its parent
            directory. Unused by a serial run. What pyEMU deletes there, and what
            is refused, is as for :meth:`run_ies`'s ``master_dir`` (ledger 361).
        exe : str, optional
            PEST++ IES executable, default ``"pestpp-ies"``, resolved as in
            :meth:`run_ies`.
        **pestpp_options
            Additional ``pst.pestpp_options``, passed straight through and
            applied after ``noise``; ``ies_par_en`` is replaced as in
            :meth:`run_ies`.

        Returns
        -------
        IesResults
            Use :meth:`~...ies.IesResults.plot_prior_vs_obs` (grey prior spaghetti
            vs measured) and :meth:`~...ies.IesResults.conflict` to inspect it.

        Raises
        ------
        ValueError
            When :meth:`build` has not run, or ``cal.pst`` has been put in
            ``regularization`` mode; or as :meth:`run_ies` raises it for a
            parallel run's master and worker directories (ledger 361).
        FileExistsError
            As :meth:`run_ies` raises it, for a parallel run's master and worker
            directories (ledger 361).
        """

        from myflopy.modflow.mf6.pest.ies import IesResults

        results_dir = self._launch_pestpp_ies(
            reals=reals,
            noptmax=-1,
            workers=workers,
            exe=exe,
            master_suffix="prior_master",
            master_dir=master_dir,
            noise=noise,
            bad_phi_sigma=None,
            pestpp_options=pestpp_options,
        )
        return IesResults(results_dir, case_name=self.name, model=self.model)

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
        laughable. No-op when there are no geostatistical parameters.

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
        correlated = [
            spec
            for spec in self._native_parameter_specs
            if spec.style in ("grid", "pilotpoints") and spec.correlation is not None
        ]
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
        """

        frames = self._pilot_point_frames.get(spec.name)
        if not frames:
            return None
        import pandas as pd

        pyemu = self.pyemu or _import_pyemu()
        geostruct = self._geostruct_for(spec)
        columns = []
        for frame, _ in frames:
            # One covariance per LAYER's point set: the layers are separate
            # parameter groups here, and correlating across them would assert a
            # vertical structure the variogram never described.
            names = [str(name) for name in frame["parnme"]]
            with self._quiet_pyemu_context():
                cov = geostruct.covariance_matrix(
                    frame["x"].to_numpy(dtype=float),
                    frame["y"].to_numpy(dtype=float),
                    names=names,
                )
                drawn = pyemu.ParameterEnsemble.from_gaussian_draw(
                    pst=self.pst, cov=cov, num_reals=int(reals), by_groups=False,
                    fill=False,
                )
            values = drawn._df if hasattr(drawn, "_df") else drawn
            columns.append(values.loc[:, [name for name in names if name in values.columns]])
        if not columns:
            return None
        return pd.concat(columns, axis=1)

    def _refuse_what_start_workers_would_delete(
        self, master: Path, workers: int, *, method: str, default: bool
    ) -> None:
        """Refuse a parallel run whose master or worker directories hold anything but a copy of this template.

        pyEMU's ``start_workers`` DELETES ``master_dir`` and ``worker_0`` ..
        ``worker_<workers-1>`` in its parent (the ``worker_root`` this passes) whenever
        they exist, whatever they hold, then copies the template into each. Before
        ledger 361 nothing checked them: a ``master_dir=`` naming a folder of the
        user's, or the run directory, was deleted, and in the DEFAULT layout every
        calibration of one model shares ``<model>.pest/`` as its worker root, so a
        second calibration's parallel run deleted the first's workers while it ran.

        A run deletes only one that is absent, empty, or a COPY of this calibration's
        template -- pyEMU copies the template's :data:`BUILD_MARKER` with it, naming
        the template, so a master or worker a previous run of THIS calibration left
        counts and another calibration's does not. A copy of a template built before
        ledger 360 has no marker; it counts in the default layout (``master_dir``
        omitted) when its metadata names this calibration and this run. Everything is
        checked before pyEMU is started.

        Raises
        ------
        ValueError
            When the master or a worker directory is, contains or sits inside the run
            directory or the template, or the master is one of its own workers.
        FileExistsError
            When one is a file, or a non-empty directory that is not a copy of this
            calibration's template, naming it, what it holds and the spelling to use.
        """

        run = self.original_workspace.resolve()
        template = self.template_workspace.resolve()
        master = Path(master)
        directories = [("master directory", master.resolve())] + [
            ("worker directory", (master.parent / f"worker_{index}").resolve())
            for index in range(workers)
        ]
        default_master = (
            self.template_workspace.parent
            / f"{self.name}_{'ies_master' if method == 'run_ies' else 'prior_master'}"
        )
        lines = ["Run it somewhere else:"]
        if not default:
            lines.append(
                f"    cal.{method}(workers={workers}, ...)   # the default: {default_master}, "
                "its workers beside it"
            )
        lines.append(
            f"    cal.{method}(workers={workers}, master_dir=<a new or empty folder>, ...)"
            f"   # its workers go beside it, worker_0 .. worker_{workers - 1}"
        )
        instead = "\n".join(lines)
        where = f"Refusing to start {method}(workers={workers}): the "
        deletes = (
            "pyEMU's start_workers DELETES the master directory and worker_0 .. "
            f"worker_{workers - 1} beside it whenever they exist, whatever they hold, then "
            "copies the template into each."
        )

        for what, directory in directories:
            for other, whose in (
                (run, "the model's run directory"),
                (template, "this calibration's template"),
            ):
                relation = _relation(directory, other)
                if relation is not None:
                    raise ValueError(
                        f"{where}{what} {directory} {relation} {whose}"
                        + ("" if relation == "IS" else f" {other}")
                        + f". {deletes}\n\n{instead}"
                    )
        if directories[0][1] in {directory for _, directory in directories[1:]}:
            raise ValueError(
                f"{where}master directory {directories[0][1]} is also one of its worker "
                f"directories, worker_0 .. worker_{workers - 1}. {deletes}\n\n{instead}"
            )

        for what, directory in directories:
            if not os.path.lexists(directory):
                continue
            if not directory.is_dir():
                raise FileExistsError(
                    f"{where}{what} {directory} is a FILE. {deletes}\n\n{instead}"
                )
            entries = sorted(entry.name for entry in directory.iterdir())
            if not entries or _written_by_a_build(directory, "template", naming=template):
                continue
            if (
                default
                and not os.path.lexists(directory / BUILD_MARKER)
                and self._metadata_names_this_calibration(directory)
            ):
                continue
            shown = ", ".join(entries[:5]) + (" ..." if len(entries) > 5 else "")
            raise FileExistsError(
                f"{where}{what} {directory} already holds {len(entries)} "
                f"entr{'y' if len(entries) == 1 else 'ies'} ({shown}), and it is not a "
                f"copy of this calibration's template: its {BUILD_MARKER} would name "
                f"{template}. {deletes} So a run deletes only one that is absent, empty, "
                "or a copy pyEMU made of THIS template; a copy of another calibration's "
                "template may be that calibration's run, still going.\n\n"
                + instead
                + f"\nIf {directory} is an old master or worker directory you no longer "
                "need, delete it yourself."
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
        pestpp_options: dict | None = None,
    ) -> Path:
        """Configure options, write the pst, run PESTPP-IES, and return the results dir."""

        if self.pst is None:
            raise ValueError("Call build() before running PESTPP-IES.")
        parallel = bool(workers and int(workers) > 1)
        if parallel:
            master = Path(master_dir) if master_dir else self.template_workspace.parent / f"{self.name}_{master_suffix}"
            # Before anything is written or started: pyEMU deletes these (ledger 361).
            self._refuse_what_start_workers_would_delete(
                master,
                int(workers),
                method="run_ies" if master_suffix == "ies_master" else "prior",
                default=not master_dir,
            )
        case = f"{self.name}.pst"
        self.pst.pestpp_options["ies_num_reals"] = int(reals)
        if not noise:
            self.pst.pestpp_options["ies_no_noise"] = True
        if bad_phi_sigma is not None:
            self.pst.pestpp_options["ies_bad_phi_sigma"] = float(bad_phi_sigma)
        for key, value in (pestpp_options or {}).items():
            self.pst.pestpp_options[key] = value
        self._refuse_regularization_mode()
        self._inject_geostatistical_prior(int(reals))
        self.pst.control_data.noptmax = int(noptmax)
        with self._quiet_pyemu_context():
            self.pst.write(str(self.template_workspace / case), version=2)

        exe_name, exe_dir = self._resolve_pestpp(exe)
        pyemu = self.pyemu or _import_pyemu()
        if parallel:
            with self._augmented_path(exe_dir):
                pyemu.os_utils.start_workers(
                    str(self.template_workspace),
                    exe_name,
                    case,
                    num_workers=int(workers),
                    worker_root=str(master.parent),
                    master_dir=str(master),
                )
            return master

        with self._augmented_path(exe_dir):
            pyemu.os_utils.run(f"{exe_name} {case}", cwd=str(self.template_workspace))
        return self.template_workspace
