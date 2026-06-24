"""Discover and open the PEST runs that have been done on a model.

A *PEST run* is one calibration setup (a ``PstFrom`` build, identified by its
``myflopy_pest_metadata.json``) together with whatever PESTPP-IES executions have
been launched from it (a prior Monte Carlo and/or an IES history match, each in
its own master directory). :func:`find_pest_runs` discovers them under a root;
:class:`PestRunHandle` opens one for review via the native :class:`IesResults`.

This is the workflow glue behind ``model.pest_runs`` and ``run.pest_runs``: by
default :class:`~...project.PestProject` writes under ``<model workspace>/pest/``,
so the runs sit beside the model and are found automatically.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

METADATA_FILENAME = "myflopy_pest_metadata.json"


@dataclass
class PestRunHandle:
    """A discovered PEST run: its template build plus any IES/prior executions.

    Attributes
    ----------
    name
        The PEST run (project) name.
    model_name
        Name of the model this run calibrates (from the build metadata).
    template_dir
        The ``PstFrom`` template directory (holds the ``.pst`` and metadata).
    masters
        Discovered execution directories keyed by kind -- e.g. ``{"ies": ...,
        "prior": ...}``.
    case
        The control-file stem (``<case>.pst``).
    """

    name: str
    model_name: str | None
    template_dir: Path
    masters: dict[str, Path] = field(default_factory=dict)
    case: str = ""
    _model: object = field(default=None, repr=False)

    @property
    def kinds(self) -> list[str]:
        """The execution kinds available to :meth:`review` (e.g. ``["ies", "prior"]``)."""

        return sorted(self.masters)

    def review(self, kind: str = "ies", *, model=None):
        """Open one execution of this PEST run as :class:`IesResults`.

        Parameters
        ----------
        kind
            Which execution to open: ``"ies"`` (history match, default) or
            ``"prior"`` (prior Monte Carlo). Falls back to whatever is available.
        model
            A myflopy model carrying the grid, enabling spatial parameter maps
            (:meth:`IesResults.plot_field`). Defaults to the model the handle was
            discovered from (``model.pest_runs``).
        """

        from myflopy.modflow.mf6.pest.ies import open_ies_run

        target = self.masters.get(kind)
        if target is None:
            target = next(iter(self.masters.values()), self.template_dir)
        return open_ies_run(target, case_name=self.case or None, model=model or self._model)

    def __repr__(self) -> str:
        kinds = ", ".join(self.kinds) or "built (not run)"
        model = f" on {self.model_name}" if self.model_name else ""
        return f"<PestRun {self.name!r}{model}: {kinds}>"


def _master_kind(dirname: str) -> str:
    lower = dirname.lower()
    if "prior" in lower:
        return "prior"
    if "ies" in lower:
        return "ies"
    return "run"


def find_pest_runs(root, *, model_name: str | None = None, model=None) -> list[PestRunHandle]:
    """Discover the PEST runs under ``root``.

    Scans for ``myflopy_pest_metadata.json`` (one per ``PstFrom`` build) and pairs
    each build with its sibling execution directories (``<name>_ies_master`` /
    ``<name>_prior_master``).

    Parameters
    ----------
    root
        Directory to scan (recursively). For ``model.pest_runs`` this is
        ``<model workspace>/pest``.
    model_name
        When given, keep only runs whose build metadata calibrates this model.
    model
        Optional model attached to each handle so ``handle.review()`` enables
        spatial maps without re-passing it.

    Returns
    -------
    list[PestRunHandle]
        One handle per discovered build, sorted by name.
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
            )
        )
    return handles
