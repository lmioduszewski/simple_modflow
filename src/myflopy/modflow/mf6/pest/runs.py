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
from dataclasses import dataclass, field
from pathlib import Path

METADATA_FILENAME = "myflopy_pest_metadata.json"

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
    masters : dict of str to pathlib.Path, optional
        Parallel-execution master directories found beside the template (a
        sibling directory named ``<name>_*master*`` holding a ``.pst``), keyed
        by kind: ``"prior"`` when the directory name contains ``prior``,
        ``"ies"`` when it contains ``ies``, else ``"run"``. Two directories of
        one kind keep the later name in sort order. Default: empty.
    case : str, optional
        The control-file stem ``review`` opens (``<case>.pst``): the stem of
        the build metadata's ``pst_file`` -- the file ``cal.build(filename)``
        wrote -- else ``name``. Default ``""``, which lets ``review`` discover
        the case from the ``.pst`` files it finds.
    """

    name: str
    model_name: str | None
    template_dir: Path
    masters: dict[str, Path] = field(default_factory=dict)
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

        return sorted(self.masters)

    def review(self, kind: str = "ies", *, model=None):
        """Open one execution of this PEST run as :class:`IesResults`.

        Parameters
        ----------
        kind : str, optional
            Which execution to open, a key of :attr:`masters`: ``"ies"``
            (history match, default), ``"prior"`` (prior Monte Carlo) or
            ``"run"``. When there is no master of that kind it falls back,
            SILENTLY, to the first master in ``masters`` order, and to
            :attr:`template_dir` when there is none. Only a PARALLEL run gets a
            master, so ``"ies"`` after a serial ``run_ies`` plus a parallel
            ``prior`` opens the prior's directory, not the template holding the
            IES results.
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
        """

        from myflopy.modflow.mf6.pest.ies import open_ies_run

        target = self.masters.get(kind)
        if target is None:
            target = next(iter(self.masters.values()), self.template_dir)
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


def _master_kind(dirname: str) -> str:
    """Classify a PEST master directory name as ``"prior"``, ``"ies"``, or ``"run"``."""

    lower = dirname.lower()
    if "prior" in lower:
        return "prior"
    if "ies" in lower:
        return "ies"
    return "run"


def _is_inside_another_template(meta_path: Path, root: Path) -> bool:
    """Is this build's metadata nested inside another build's directory?"""

    for ancestor in meta_path.parent.parents:
        if (ancestor / METADATA_FILENAME).exists():
            return True
        if ancestor == root or ancestor == ancestor.parent:
            return False
    return False


def find_pest_runs(root, *, model_name: str | None = None, model=None,
                   model_factory=None) -> list[PestRunHandle]:
    """Discover the PEST runs under ``root``.

    Scans for ``myflopy_pest_metadata.json`` (one per ``PstFrom`` build) and pairs
    each build with its sibling execution directories (``<name>_ies_master`` /
    ``<name>_prior_master``).

    Parameters
    ----------
    root : str or pathlib.Path
        Directory to scan (recursively). A root that does not exist gives an
        empty list. A metadata file inside a ``*_master`` directory, or nested
        inside another build's directory, is a copy and is skipped; an
        unreadable one is skipped silently. ``model.pest_runs`` calls
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

    handles: list[PestRunHandle] = []
    for meta_path in sorted(root.rglob(METADATA_FILENAME)):
        # PESTPP-IES execution dirs (``<run>_ies_master`` / ``<run>_prior_master``,
        # created for parallel runs) are clones of the template and carry a copy
        # of its metadata. The run is represented by its template dir, so skip the
        # master copies to avoid listing the same run several times.
        if meta_path.parent.name.endswith("_master"):
            continue
        # Same reason, one level up: a template is a COPY of the model workspace,
        # so a calibration that lived inside that workspace (the pre-2026-07-29
        # default) is copied into every template built afterwards. A run nested
        # inside another run's template is that copy, never a run of its own.
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
        case = Path(meta.get("pst_file") or "").stem or name

        masters: dict[str, Path] = {}
        for candidate in sorted(template_dir.parent.glob(f"{name}_*master*")):
            if candidate.is_dir() and any(candidate.glob("*.pst")):
                masters[_master_kind(candidate.name)] = candidate

        handles.append(
            PestRunHandle(
                name=name,
                model_name=run_model,
                template_dir=template_dir,
                masters=masters,
                case=case,
                _model=model,
                _model_factory=model_factory,
            )
        )
    return handles
