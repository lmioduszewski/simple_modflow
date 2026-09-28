"""`ModelGroup` — compare multiple models with one lazy API."""

from __future__ import annotations

import pandas as pd

from myflopy._defaults import DEFAULT_CRS
from myflopy._deprecation import deprecated_instance_getattr
from myflopy.project.group._shared import _coerce_models, _grid_signature_for_model
from myflopy.project.group.budget import GroupBudget
from myflopy.project.group.conc import GroupConc
from myflopy.project.group.inputs import GroupPackageInputs
from myflopy.project.group.packages import GroupOutputs, GroupPackages
from myflopy.project.group.spatial import GroupHeads
from myflopy.project.group.temp import GroupTemp
from myflopy.project.group.uzf import GroupUzfInputs


class ModelGroup:
    """Collection of models with group-wise get/compare helpers.

    Parameters
    ----------
    models : Mapping or Sequence
        Mapping of names to models/workspace directories, or a sequence of
        model objects/workspace directories. An entry that is a ``str`` or
        :class:`~pathlib.Path` is opened with :func:`load_mf6_run` as a lazy,
        file-backed :class:`LoadedMf6Run`; any other entry is used as given (a
        built or loaded model). A mapping's keys become the group keys through
        ``str()``, in the mapping's order. For a sequence, each loaded model's
        ``name`` attribute is the group key, so two models of the same ``name``
        keep only the LAST. A bare ``str``/``bytes`` or any other type raises
        ``TypeError``; an empty mapping or sequence raises ``ValueError``.
    reference : str or None, optional
        Group key of the reference model: the baseline for the ``compare()``
        methods and for :meth:`diff`, the ``is_reference`` column of
        :meth:`summary`, and the anchor for ``shared_grid``. ``None`` (the
        default) takes the first model. A key not in the group raises
        ``KeyError``.
    crs : str, optional
        CRS (any string pyproj accepts, e.g. ``"EPSG:2927"``) given to
        :func:`load_mf6_run` for the workspace-path entries; model objects are
        not touched. Default :data:`~myflopy._defaults.DEFAULT_CRS`
        (``"EPSG:2927"``). It is always passed explicitly, so the CRS a myflopy
        run recorded in its ``run.json`` is NOT read here, unlike
        ``load_mf6_run(ws)`` with its own ``crs=None`` default.
    verbosity_level : int, optional
        Given to :func:`load_mf6_run` for the workspace-path entries. ``0`` (the
        default) is silent; above 0 prints the loader's progress messages and is
        passed to FloPy's ``MFSimulation.load`` when the files are read.
    shared_grid : bool, optional
        If ``True``, confirm every model has the same discretization and reuse
        the reference model's lazily built Voronoi/grid view (``.vor``) for the
        whole group, which avoids rebuilding identical grids. A model's grid is
        compared by its ``grid_signature()`` (the SHA-256 of a
        :class:`LoadedMf6Run`'s ``.dis``/``.disv``/``.disu`` file), else by a
        hash of the ``<name>.<grid_type>`` file in its workspace with blank and
        ``#`` comment lines dropped. A model with neither, or one whose signature
        differs from the reference's, raises ``ValueError``. Default ``False``.
    """

    def __init__(
        self,
        models,
        *,
        reference: str | None = None,
        crs: str = DEFAULT_CRS,
        verbosity_level: int = 0,
        shared_grid: bool = False,
    ):
        """Load ``models`` into the group, pick the reference, and wire the accessors.

        Coerces the ``models`` mapping/sequence to loaded models, validates the
        chosen ``reference``, optionally reconciles a shared grid, and constructs
        the heads/outputs/package accessors. See the class docstring for parameters.
        """

        self.models = _coerce_models(models, crs=crs, verbosity_level=verbosity_level)
        if not self.models:
            raise ValueError("ModelGroup requires at least one model.")

        self.reference = reference or next(iter(self.models))
        if self.reference not in self.models:
            raise KeyError(f"Reference model {self.reference!r} is not in the group.")
        self.shared_grid = bool(shared_grid)

        if self.shared_grid:
            self._configure_shared_grid()

        self.hds = GroupHeads(self)
        # Constructed eagerly like every other accessor, and deliberately NOT
        # kind-gated here: the gate lives on the member models' readers, so a
        # group of flow models still builds and only errors if someone asks it
        # for concentration.
        self.conc = GroupConc(self)
        self.temp = GroupTemp(self)
        self.outputs = GroupOutputs(self)
        self._rch = GroupPackageInputs(self, "rch")
        self._chd = GroupPackageInputs(self, "chd")
        self._drn = GroupPackageInputs(self, "drn")
        self._ghb = GroupPackageInputs(self, "ghb")
        self._riv = GroupPackageInputs(self, "riv")
        self._wel = GroupPackageInputs(self, "wel")
        self._evt = GroupPackageInputs(self, "evt")
        self._uzf = GroupUzfInputs(self)
        self.packages = GroupPackages(self)

    # -- deprecated flat input shortcuts --------------------------------------
    # These duplicate ``group.packages.<pkg>.inputs`` (which mirrors the
    # single-model ``model.<pkg>.inputs``) and have no single-model
    # equivalent, so they are deprecated in favor of the mirrored path.
    # Resolved only via __getattr__ so they stay out of dir()/completion (D12).
    __getattr__ = deprecated_instance_getattr(
        {
            "rch": ("_rch", "group.packages.rch.inputs", "0.1.0"),
            "chd": ("_chd", "group.packages.chd.inputs", "0.1.0"),
            "drn": ("_drn", "group.packages.drn.inputs", "0.1.0"),
            "ghb": ("_ghb", "group.packages.ghb.inputs", "0.1.0"),
            "wel": ("_wel", "group.packages.wel.inputs", "0.1.0"),
            "uzf": ("_uzf", "group.packages.uzf.inputs", "0.1.0"),
        },
        "myflopy.project.model_group.ModelGroup",
    )

    @property
    def run_ids(self) -> list[str]:
        """Return the ordered model names in the group."""

        return list(self.models)

    def summary(self) -> pd.DataFrame:
        """Return per-model summaries stacked into one table."""

        frames = []
        for model_name, model in self.models.items():
            frame = model.summary().copy()
            frame["group_model"] = model_name
            frame["is_reference"] = model_name == self.reference
            frames.append(frame)
        return pd.concat(frames, ignore_index=True)

    def bud(self, package: str | None = None) -> GroupBudget:
        """Return a grouped budget accessor.

        Parameters
        ----------
        package : str or None, optional
            Default package for the accessor's ``get()``/``compare()``, as each
            member model's ``bud(package)`` takes it (e.g. ``"rch"`` or
            ``"drn"``); lower-cased. ``package=`` on those calls overrides it.
            ``None`` (the default) sets none, and ``get()``/``compare()`` then
            raise ``ValueError`` unless they are given ``package=``.
        """

        return GroupBudget(self, package=package)

    def diff(self):
        """Return a :class:`~myflopy.project.model_diff.ModelDiff` for this group.

        Reference-star: every other model is compared against the group's
        reference. See :class:`~myflopy.project.model_diff.ModelDiff` for the
        structural + value difference tiers, ``summary()``, and ``report()``.

        Returns
        -------
        ModelDiff
            The single ``diff`` verb over this group.

        Examples
        --------
        >>> group = mf.ModelGroup({"base": run_a, "variant": run_b}, reference="base")
        >>> group.diff().report()                       # Markdown setup diff
        >>> group.diff().hds.map("variant", per=8)       # Δhead map vs reference
        """

        from myflopy.project.model_diff import ModelDiff

        return ModelDiff(self)

    def _configure_shared_grid(self):
        """Enable shared lazy grid reuse across models with identical grids."""

        anchor_name = self.reference
        anchor_model = self.models[anchor_name]
        anchor_signature = _grid_signature_for_model(anchor_model)
        if anchor_signature is None:
            raise ValueError(
                "shared_grid=True requires discretization files or grid signatures "
                f"for reference model {anchor_name!r}."
            )

        for model_name, model in self.models.items():
            signature = _grid_signature_for_model(model)
            if signature is None:
                raise ValueError(
                    "shared_grid=True requires discretization files or grid signatures "
                    f"for model {model_name!r}."
                )
            if signature != anchor_signature:
                raise ValueError(
                    "shared_grid=True requires all models to use the same discretization. "
                    f"Model {model_name!r} does not match reference model {anchor_name!r}."
                )

        for model_name, model in self.models.items():
            if model_name == anchor_name:
                continue
            model._shared_vor_source = anchor_model
