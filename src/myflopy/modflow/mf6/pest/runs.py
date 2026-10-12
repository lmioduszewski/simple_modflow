"""Discover and open the PEST runs that have been done on a model.

A *PEST run* is one calibration setup (a ``PstFrom`` build, identified by its
``myflopy_pest_metadata.json``) together with whatever PESTPP-IES executions have
been launched from it (a prior Monte Carlo and/or an IES history match, each in
its own master directory). :func:`find_pest_runs` discovers them under a root;
:class:`PestRunHandle` opens one for review via the native :class:`IesResults`.

This is the workflow glue behind ``model.pest_runs`` and ``run.pest_runs``: by
default :class:`~...project.PestProject` writes under ``<model workspace>.pest/``
-- a SIBLING of the model directory -- so the runs sit next to the model and are
found automatically. Calibrations built before 2026-07-29 live in
``<model workspace>/pest/`` and are still discovered; see
:func:`pest_run_roots`.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

METADATA_FILENAME = "myflopy_pest_metadata.json"

#: The marker a build writes into its template, naming the template's resolved path
#: (ledger 360). pyEMU copies it into every IES master and worker, so a directory
#: whose marker names ANOTHER directory is a copy (ledger 361).
BUILD_MARKER = "myflopy_pest_build.json"

#: The execution kinds :meth:`PestProject.run_ies` and :meth:`PestProject.prior` launch,
#: each with the suffix of its default parallel master directory, ``<name>_<suffix>``.
EXECUTION_MASTERS = {"ies": "ies_master", "prior": "prior_master"}

#: Every file a PESTPP-IES run of case ``<case>`` writes that a review reads, after the
#: ``<case>.`` prefix: one per iteration (``<case>.<n>.obs.csv``, ``.par.jcb``,
#: ``.base.par``, ...), the phi files, the noise ensemble and the base run's final
#: parameters (ledger 507).
_RUN_OUTPUT = re.compile(r"(\d+\..+|phi\..+|obs\+noise\..+|base\.(par|rei))")

#: Appended to the model directory's NAME (not joined as a child) to give the
#: default calibration root. A sibling, because ``PstFrom`` copies the model
#: workspace wholesale and cannot be told to skip a subdirectory -- see
#: compromise ledger 107.
PEST_ROOT_SUFFIX = ".pest"

#: Where calibrations used to live: a child of the model workspace, and
#: therefore inside the tree ``PstFrom`` copies. Still searched so runs already
#: on disk keep reviewing.
LEGACY_PEST_DIRNAME = "pest"


def default_pest_root(model_workspace) -> Path:
    """The directory ``model.pest(name)`` builds into by default."""

    workspace = Path(model_workspace)
    return workspace.parent / f"{workspace.name}{PEST_ROOT_SUFFIX}"


def pest_run_roots(model_workspace) -> tuple[Path, ...]:
    """Every directory a calibration of this model may live in, newest first."""

    workspace = Path(model_workspace)
    return (default_pest_root(workspace), workspace / LEGACY_PEST_DIRNAME)


def pestpp_ies_outputs(directory, case: str) -> list[Path]:
    """The files in ``directory`` a PESTPP-IES run of ``case`` wrote that a review reads.

    One per iteration (``<case>.<n>.obs.csv``/``.par.csv`` and their ``.jcb`` forms,
    ``<case>.<n>.base.par``, ...), the phi files (``<case>.phi.*``), the noise ensemble
    (``<case>.obs+noise.*``) and ``<case>.base.par``/``.base.rei``. Not the control file,
    its external tables, the ``.rec`` record or ``<case>.prior_par.csv``. What a launch
    clears before a serial run, and what a parallel run's master inherited from the
    template (ledger 507).

    Parameters
    ----------
    directory : str or pathlib.Path
        The directory to list (not recursed into). One that does not exist gives an
        empty list.
    case : str
        The control-file stem the run used (``<case>.pst``); matched exactly, so the
        outputs of a case ``valley`` are not those of ``valley_ss``.

    Returns
    -------
    list of pathlib.Path
        The matching FILES, sorted by name.
    """

    directory = Path(directory)
    if not directory.is_dir():
        return []
    prefix = f"{case}."
    return sorted(
        path for path in directory.iterdir()
        if path.name.startswith(prefix)
        and _RUN_OUTPUT.fullmatch(path.name[len(prefix):])
        and path.is_file()
    )


def execution_kind(directory, case: str) -> str | None:
    """Which execution a directory's ``<case>.pst`` launched: ``"ies"``, ``"prior"`` or ``None``.

    Read from the control file the launch wrote (``NOPTMAX -1`` is
    :meth:`PestProject.prior`, anything else :meth:`PestProject.run_ies`), and only when
    ``<case>.phi.actual.csv`` sits beside it -- a run that wrote no phi is not one a
    review can open. ``None`` when either file is absent or the control file names no
    ``noptmax``.

    Parameters
    ----------
    directory : str or pathlib.Path
        The template (a serial run) or a master directory (a parallel one).
    case : str
        The control-file stem the run used.

    Returns
    -------
    str or None
        ``"prior"``, ``"ies"`` or ``None``.
    """

    directory = Path(directory)
    if not (directory / f"{case}.phi.actual.csv").is_file():
        return None
    try:
        with open(directory / f"{case}.pst", encoding="utf-8") as handle:
            for line in handle:
                words = line.split()
                if len(words) >= 2 and words[0].lower() == "noptmax":
                    return "prior" if int(float(words[1])) == -1 else "ies"
    except (OSError, ValueError):
        return None
    return None


def _marker_names(directory: Path) -> str | None:
    """The directory a :data:`BUILD_MARKER` in ``directory`` names, or ``None`` when it holds none."""

    try:
        marker = json.loads((directory / BUILD_MARKER).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    named = marker.get("directory") if isinstance(marker, dict) else None
    return str(named) if named else None


def _resolved(directory) -> str:
    """``directory``'s resolved path, as a :data:`BUILD_MARKER` records it."""

    return str(Path(directory).resolve())


def _a_master_of(master: Path, template_dir: Path) -> bool:
    """Was ``master`` copied from the template ``template_dir``?

    Yes when its build marker names ``template_dir`` -- or names what the TEMPLATE's
    own marker names, which is the same template after the calibration root was moved
    or copied elsewhere, since a marker keeps the path the build wrote (ledger
    516) -- or when it holds no marker (a template built before ledger 360).
    No: a master copied from ANOTHER template, whose marker names that one.
    """

    named = _marker_names(master)
    return named is None or named in (_resolved(template_dir), _marker_names(template_dir))


def _results_written_at(directory: Path, case: str) -> float:
    """When ``directory``'s run of ``case`` last wrote its phi (else its control file); 0 when neither exists."""

    for name in (f"{case}.phi.actual.csv", f"{case}.pst"):
        try:
            return (directory / name).stat().st_mtime
        except OSError:
            continue
    return 0.0


def find_model_pest_runs(model_workspace, **kwargs) -> list[PestRunHandle]:
    """Discover a model's PEST runs across the current and legacy roots."""

    handles: list[PestRunHandle] = []
    seen: set[Path] = set()
    for root in pest_run_roots(model_workspace):
        for handle in find_pest_runs(root, **kwargs):
            # A template and its legacy copy would otherwise both be listed.
            key = handle.template_dir.resolve()
            if key in seen:
                continue
            seen.add(key)
            handles.append(handle)
    return sorted(handles, key=lambda handle: handle.name)


@dataclass
class PestRunHandle:
    """A discovered PEST run: its template build plus any IES/prior executions.

    Handed out by :func:`find_pest_runs`, ``model.pest_runs`` and
    ``run.pest_runs``, which fill every field from what is on disk; you do not
    normally construct one.

    Attributes
    ----------
    name : str
        The PEST run (project) name: the metadata's ``project_name``, else the
        template directory's name.
    model_name : str or None
        Name of the model this run calibrates, from the build metadata's
        ``model_name``; ``None`` when the metadata records none.
    template_dir : pathlib.Path
        The ``PstFrom`` template directory (holds the ``.pst`` and
        ``myflopy_pest_metadata.json``). A SERIAL ``run_ies``/``prior`` writes
        its results here too.
    executions : dict of str to pathlib.Path, optional
        The directory holding each execution's results, keyed by kind, ``"ies"``
        (:meth:`PestProject.run_ies`) or ``"prior"`` (:meth:`PestProject.prior`). A
        PARALLEL run's is its default master, found by its EXACT name
        ``<name>_ies_master``/``<name>_prior_master`` beside the template (holding a
        ``.pst``, and not a copy of another calibration's template: its build marker
        names this template, or what this template's own marker names after the
        root was moved, or it holds none; ledger 516); a SERIAL run's is
        the template, its kind read from the ``<name>.pst`` the launch wrote there
        (``NOPTMAX -1`` is a prior), counted only once it holds a phi file. When both
        hold a run of one kind, the one whose phi was written last. A master passed as
        ``master_dir=`` elsewhere is not found. Default: empty (built, not run).
        Before ledger 508 this was ``masters``: every sibling whose name
        started ``<name>_`` and contained ``master``, classified by a substring, and
        serial runs were never listed.
    case : str, optional
        The control-file stem ``review`` opens (``<case>.pst``): ``name``, the case
        PESTPP-IES is launched under (``<name>.pst``) whatever file
        ``cal.build(filename)`` wrote. Default ``""``, which lets ``review`` discover
        the case from the ``.pst`` files it finds. It was the build's file stem until
        ledger 508, which opened no results whenever the two differed.
    """

    name: str
    model_name: str | None
    template_dir: Path
    executions: dict[str, Path] = field(default_factory=dict)
    case: str = ""
    _model: object = field(default=None, repr=False)
    # Deferred alternative to ``_model``: a zero-argument callable returning the
    # model, or None. Discovering runs must stay side-effect free -- listing what
    # was calibrated should not build a model view, cache it, or raise because a
    # run holds several models -- so the model is resolved only when
    # :meth:`review` actually needs the grid.
    _model_factory: object = field(default=None, repr=False)

    @property
    def kinds(self) -> list[str]:
        """The execution kinds available to :meth:`review` (e.g. ``["ies", "prior"]``)."""

        return sorted(self.executions)

    def review(self, kind: str = "ies", *, model=None):
        """Open one execution of this PEST run as :class:`IesResults`.

        Parameters
        ----------
        kind : str, optional
            Which execution to open, a key of :attr:`executions`: ``"ies"``
            (history match, default) or ``"prior"`` (prior Monte Carlo), serial or
            parallel. A kind the run does not hold RAISES ``ValueError`` naming the
            kinds it does (ledger 508) -- it fell back, silently, to
            another execution, so a run holding only a prior reviewed its prior Monte
            Carlo as the history match. A run with no execution at all (built, not
            run) opens :attr:`template_dir`.
        model : SimulationBase or None, optional
            A myflopy model carrying the grid, enabling spatial parameter maps
            (:meth:`IesResults.plot_field`). Default ``None`` uses the model the
            handle was discovered with -- ``model.pest_runs`` attaches it
            directly, ``run.pest_runs`` resolves it here on first use. A run
            holding several models cannot pick one, so pass ``model=`` yourself
            there.

        Returns
        -------
        IesResults
            Opened on the chosen directory with ``case_name=`` :attr:`case`
            (discovered from the directory's ``.pst`` files when empty).

        Raises
        ------
        ValueError
            When the run holds executions but none of ``kind``.
        """

        from myflopy.modflow.mf6.pest.ies import open_ies_run

        if not self.executions:
            target = self.template_dir
        elif kind in self.executions:
            target = self.executions[kind]
        else:
            held = ", ".join(f"review(kind={name!r})" for name in self.kinds)
            raise ValueError(
                f"The PEST run {self.name!r} holds no {kind!r} execution; it holds "
                f"{', '.join(map(repr, self.kinds))}. Open one of them: {held}."
            )
        if model is None:
            model = self._model
        if model is None and self._model_factory is not None:
            model = self._model_factory()
        return open_ies_run(target, case_name=self.case or None, model=model)

    def __repr__(self) -> str:
        """Show the run name, model, and the execution kinds discovered on disk."""

        kinds = ", ".join(self.kinds) or "built (not run)"
        model = f" on {self.model_name}" if self.model_name else ""
        return f"<PestRun {self.name!r}{model}: {kinds}>"


def _is_inside_another_template(meta_path: Path, root: Path) -> bool:
    """Is this build's metadata nested inside another build's directory?"""

    for ancestor in meta_path.parent.parents:
        if (ancestor / METADATA_FILENAME).exists():
            return True
        if ancestor == root or ancestor == ancestor.parent:
            return False
    return False


#: An IES worker directory pyEMU's ``start_workers`` makes, ``worker_0`` .. ``worker_<n>``.
_WORKER_DIRECTORY = re.compile(r"worker_\d+")


def _executions(template_dir: Path, name: str) -> dict[str, Path]:
    """Each execution kind of the calibration ``name`` and the directory holding its results."""

    found: dict[str, list[Path]] = {}
    serial = execution_kind(template_dir, name)
    if serial is not None:
        found.setdefault(serial, []).append(template_dir)
    for kind, suffix in EXECUTION_MASTERS.items():
        master = template_dir.parent / f"{name}_{suffix}"
        if (
            master.is_dir()
            and any(master.glob("*.pst"))
            and _a_master_of(master, template_dir)
        ):
            found.setdefault(kind, []).append(master)
    return {
        kind: max(places, key=lambda place: _results_written_at(place, name))
        for kind, places in found.items()
    }


def find_pest_runs(root, *, model_name: str | None = None, model=None,
                   model_factory=None) -> list[PestRunHandle]:
    """Discover the PEST runs under ``root``.

    Scans for ``myflopy_pest_metadata.json`` (one per ``PstFrom`` build) and pairs
    each build with its executions: a serial run in the template, and the parallel
    masters named exactly ``<name>_ies_master`` / ``<name>_prior_master`` beside it
    (:attr:`PestRunHandle.executions`).

    Parameters
    ----------
    root : str or pathlib.Path
        Directory to scan (recursively). A root that does not exist gives an
        empty list. A metadata file that is a COPY of a template is skipped: one in a
        ``*_master`` directory, one in a directory whose ``myflopy_pest_build.json``
        names ANOTHER build directory under ``root`` (an IES worker or master pyEMU
        copied, a backup), one in a ``worker_<n>`` directory whose marker does not
        name it or that holds none (a template built before ledger 360 left them
        unmarked; ledger 508), and one nested inside another build's
        directory. A template whose marker names a path no build under ``root``
        occupies is NOT a copy: it is the template itself after its calibration
        root was moved or copied elsewhere, the marker keeping the path its build
        wrote (ledger 516). An unreadable one is skipped silently. ``model.pest_runs`` calls
        :func:`find_model_pest_runs` instead, which scans every root in
        :func:`pest_run_roots`.
    model_name : str or None, optional
        When given, keep only runs whose build metadata's ``model_name``
        equals it. Default ``None`` keeps every run.
    model : SimulationBase or None, optional
        Model attached to each handle so ``handle.review()`` enables spatial
        maps without re-passing it. Default ``None``.
    model_factory : callable or None, optional
        Zero-argument callable resolving that model on demand, for callers whose
        model cannot be resolved side-effect free while merely listing runs (see
        ``Run.pest_runs``). Called by ``review()`` only when neither
        ``model`` here nor ``review(model=)`` supplied one. Default ``None``.

    Returns
    -------
    list[PestRunHandle]
        One handle per discovered build, in the sorted order of their metadata
        file paths (:func:`find_model_pest_runs` sorts by name).
    """

    root = Path(root)
    if not root.exists():
        return []

    found = sorted(root.rglob(METADATA_FILENAME))
    # Every build directory this scan sees, so a copy is told from a moved original:
    # a marker keeps the path its build wrote, so after the calibration root is moved
    # or copied elsewhere EVERY template's marker names a directory that is not its
    # own, and only a copy's names one that is still here (ledger 516).
    builds = {_resolved(meta_path.parent) for meta_path in found}
    handles: list[PestRunHandle] = []
    for meta_path in found:
        # PESTPP-IES execution dirs (``<run>_ies_master`` / ``<run>_prior_master``,
        # created for parallel runs) are clones of the template and carry a copy
        # of its metadata. The run is represented by its template dir, so skip the
        # master copies to avoid listing the same run several times.
        if meta_path.parent.name.endswith("_master"):
            continue
        # So are the worker directories (`<master>.workers/worker_<n>`,
        # ledger 532), which a failed or interrupted parallel run, or one launched
        # with `keep_workers=True`, leaves behind, and a backup copy:
        # each one's marker names the template it was copied from (ledger
        # 508). A copy is one whose marker names ANOTHER build of this scan
        # -- never merely another path, which every template of a moved or copied
        # root names (ledger 516). A `worker_<n>` is pyEMU's name for a
        # worker, so one whose marker does not name itself is a copy whatever it
        # names (a template built before the marker existed left them unmarked).
        directory = _resolved(meta_path.parent)
        named = _marker_names(meta_path.parent)
        if named is not None and named != directory and named in builds:
            continue
        if _WORKER_DIRECTORY.fullmatch(meta_path.parent.name) and named != directory:
            continue
        # Same reason, one level up: a template built before ledger 383 is a
        # COPY of the whole model workspace, so a calibration that lived inside
        # that workspace (the pre-2026-07-29 default) was copied into every
        # template built then. A run nested inside another run's template is that
        # copy, never a run of its own. (A template now holds the run's MF6
        # inputs only, so no new one carries a nested calibration.)
        if _is_inside_another_template(meta_path, root):
            continue
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
        except (ValueError, OSError):
            continue
        name = meta.get("project_name") or meta_path.parent.name
        run_model = meta.get("model_name")
        if model_name is not None and run_model != model_name:
            continue
        template_dir = meta_path.parent

        handles.append(
            PestRunHandle(
                name=name,
                model_name=run_model,
                template_dir=template_dir,
                executions=_executions(template_dir, name),
                # PESTPP-IES runs `<name>.pst` whatever file `cal.build(filename)`
                # wrote (`PestProject._launch_pestpp_ies`), so its outputs are
                # `<name>.*` (ledger 508).
                case=name,
                _model=model,
                _model_factory=model_factory,
            )
        )
    return handles
