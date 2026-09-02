"""Unified multi-model difference engine -- the ``diff`` verb.

:class:`ModelDiff` is the single entry point for "what is different between these
models?" It is built from a :class:`~myflopy.project.model_group.ModelGroup`
(reference-star topology: the group's reference vs each other model) and reports
tiers of difference. Phase 1 ships two:

* **structural** -- which packages, and which cells within a package, exist in one
  model but not another. This is a true set difference: unlike a value-only
  compare (which inner-joins on shared cells) it surfaces cells that appear in
  only one model instead of silently dropping them.
* **value** -- numeric field differences on the cells shared by both models
  (this reuses the group's existing aligned-value ``compare``).

Two doors reach the same engine, and both return a :class:`ModelDiff`:

* ``model.diff(other)`` / ``model.diff([m2, m3])`` -- quick and ad hoc; ``self``
  is the reference. Accepts model objects or workspace paths.
* ``ModelGroup(...).diff()`` -- durable and named, with explicit reference control.

All stress-period and layer arguments are zero-based, matching the rest of the
myflopy API.
"""

from __future__ import annotations

from collections import Counter
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from myflopy._logging import get_logger
from myflopy.modflow.mf6.package_explorer import (
    DiffSpatialView,
    LeafFieldSugar,
    get_default_package_value_column,
    get_package_input_field_names,
)
from myflopy.modflow.mf6.package_registry import _PACKAGE_EXPLORER_SPECS
from myflopy.modflow.mf6.package_tables import (
    PACKAGE_TABLE_UNAVAILABLE,
    build_cell_package_input_table,
    build_lak_connection_table,
    build_sfr_reach_table,
)
from myflopy.project.model_group import GroupPackageInputs
from myflopy.project.model_results_diff import (
    BudgetResultDiff,
    CellResultsDiffNamespace,
    ConcResultDiff,
    HeadsResultDiff,
    LakResultsDiffNamespace,
    MvrResultDiff,
    SfrResultsDiffNamespace,
    TempResultDiff,
    UzfResultsDiffNamespace,
)

if TYPE_CHECKING:  # pragma: no cover - typing only
    from myflopy.project.model_group import ModelGroup

# Cell-based stress-period BC packages the Phase-1 diff understands, derived
# from the registry (plan 4.7.3). The ModelGroup's GroupPackageInputs accessors
# (project/group/core.py) are still written out by hand and must stay in step;
# ``test_diff_tier_covers_exactly_the_group_cell_bc_accessors`` pins that, and
# is now the guard on the GROUP side rather than on both.
_DIFF_PACKAGES: tuple[str, ...] = tuple(
    name for name, spec in _PACKAGE_EXPLORER_SPECS.items() if spec.tiers.diffable
)

# Advanced packages diffed by connection/reach geometry (Phase 3).
_CONNECTION_PACKAGES: tuple[str, ...] = tuple(
    name
    for name, spec in _PACKAGE_EXPLORER_SPECS.items()
    if spec.tiers.connection_diffable
)
_LAK_IDENTITY = (
    "lake", "layer", "cell", "claktype", "belev", "telev", "connlen", "connwidth",
)
_SFR_IDENTITY = ("reach", "layer", "cell", "rlen")

logger = get_logger(__name__)

#: What "read a finished model's outputs" can raise when the run never finished.
#:
#: FileNotFoundError is the obvious one, but the case this fallback is really
#: written for -- a run MF6 opened and then aborted -- leaves a ZERO-BYTE or
#: truncated .hds/.cbc, and flopy raises ValueError for those ("file is empty",
#: "min() iterable argument is empty"), not a file error.
_RESULTS_UNAVAILABLE = (
    OSError, AttributeError, ValueError, IndexError, KeyError,
)

_SUMMARY_COLUMNS = [
    "package",
    "model",
    "present_in_reference",
    "present_in_model",
    "cells_only_in_reference",
    "cells_only_in_model",
    "cells_shared",
    "value_cells_changed",
    "identical",
]


def _package_present(model, package_name: str) -> bool:
    """Return whether ``model`` has a package of the given type attached."""

    try:
        names = [str(name).lower() for name in model.package_names]
    except (AttributeError, TypeError):  # pragma: no cover - defensive
        # The group accepts duck-typed model objects, so `package_names` may be
        # absent (AttributeError) or not iterable (TypeError).
        logger.debug("cannot list packages on %r", model)
        return False
    pkg = package_name.lower()
    return any(name == pkg or name.startswith(pkg) for name in names)


def _package_keys(model, package_name: str, *, per=None, layer=None):
    """Return ``(present, keys)`` for a BC package on ``model``.

    ``keys`` is a set of zero-based ``(per, layer, cell)`` tuples identifying
    every cell the package applies to (optionally filtered by ``per``/``layer``).
    A package that is absent, empty, or unreadable yields an empty set; the
    ``present`` flag distinguishes "attached but no cells" from "not attached".
    """

    try:
        table = build_cell_package_input_table(model, package_name, per=per, layer=layer)
    except PACKAGE_TABLE_UNAVAILABLE:
        # NOTE the asymmetry this leaves, which narrowing made visible:
        # `_package_present` matches by PREFIX, so an array-form package
        # (RCHA/EVTA on an externally loaded model) reports as present while its
        # cell keys come back empty -- two such models compare as identical. The
        # cell tier cannot read array-form packages at all; telling those apart
        # needs a third state ("present, unreadable by this tier") rather than a
        # wider except. Recorded in the ledger rather than fixed here.
        logger.debug(
            "no readable %s cells on %s", package_name,
            getattr(model, "name", model), exc_info=True,
        )
        return _package_present(model, package_name), set()
    if table is None or table.empty:
        return _package_present(model, package_name), set()
    n = len(table)
    per_col = table["per"].to_numpy(dtype=int) if "per" in table.columns else np.zeros(n, dtype=int)
    layer_col = table["layer"].to_numpy(dtype=int) if "layer" in table.columns else np.zeros(n, dtype=int)
    cell_col = table["cell"].to_numpy(dtype=int)
    keys = set(zip(per_col.tolist(), layer_col.tolist(), cell_col.tolist(), strict=False))
    return True, keys


class PackageDiff(LeafFieldSugar, DiffSpatialView):
    """Difference view for one cell-based BC package across the group.

    Reached via ``diff.packages.<package>.inputs`` (e.g. ``.ghb.inputs``).
    Exposes the structural tier (:meth:`cells`), the value tier
    (:meth:`values`), the unified delta grammar (``map``/``plot``/``mosaic``/
    ``animate``, with registry fields as first-class nodes -- ``inputs.cond``
    or ``field="cond"``), and per-model counts (:meth:`summary`), always
    relative to the group's reference model.
    """

    def __init__(self, diff: ModelDiff, package_name: str, field_name: str | None = None):
        """Bind one BC package's diff view; ``field_name`` pins it to a single field."""

        self._diff = diff
        self.group = diff.group
        self.package_name = str(package_name).lower()
        self.field_name = None if field_name is None else str(field_name).lower()
        self._inputs = GroupPackageInputs(self.group, self.package_name)

    # -- field nodes (LeafFieldSugar hooks) ------------------------------------
    def _field_names(self) -> list[str]:
        """The registry-declared input field names for this package."""

        return get_package_input_field_names(self.package_name)

    def _field_node(self, name: str) -> PackageDiff:
        """Return a copy of this diff view pinned to field ``name``."""

        return PackageDiff(self._diff, self.package_name, field_name=name)

    # -- spatial-view hooks (delta maps vs the reference) ---------------------
    def _spatial_map(self, *, per=0, layer=0, model=None, **kwargs):
        """Delta choropleth (compared model - reference) for one model."""

        if self.field_name is not None:
            kwargs.setdefault("value_column", self.field_name)
        return self._inputs.compare_map(model_name=model, per=per, layer=layer, **kwargs)

    def _spatial_models(self):
        """The non-reference model names -- one delta panel each."""

        return [name for name in self.group.models if name != self.group.reference]

    def _spatial_reference_model(self):
        """The reference model whose grid/layers frame the delta maps."""

        return self.group.models[self.group.reference]

    def _spatial_periods(self):
        """Stress periods present in this package's data (delegates to the group inputs)."""

        return self._inputs._spatial_periods()

    def _spatial_layers(self):
        """Layers present in this package's data (delegates to the group inputs)."""

        return self._inputs._spatial_layers()

    def _spatial_value_label(self):
        """Label for the mapped quantity: the pinned field, else the package name."""

        return self.field_name or self.package_name

    # -- series hooks: plot() draws the field's Δ by period per model ----------
    def _series_table(self) -> pd.DataFrame:
        """The aligned value-diff table backing ``plot()`` (the group's compare)."""

        return self._inputs.compare()

    def _series_value_column(self, frame) -> str:
        """The ``<field>_diff`` column ``plot()`` draws (validated against ``frame``)."""

        field = self.field_name or get_default_package_value_column(self.package_name)
        column = f"{field}_diff"
        if column not in getattr(frame, "columns", []):
            raise KeyError(
                f"Aligned diff column {column!r} was not found; available: "
                f"{list(getattr(frame, 'columns', []))}."
            )
        return column

    def _targets(self, model_name=None) -> list[str]:
        """Return the non-reference model names to diff (or just ``model_name``)."""

        if model_name is not None:
            name = str(model_name)
            if name not in self.group.models:
                raise KeyError(f"Model {name!r} is not in the group.")
            if name == self.group.reference:
                raise ValueError("The reference model cannot be diffed against itself.")
            return [name]
        return [name for name in self.group.models if name != self.group.reference]

    # -- value tier (reuses the group's aligned-value compare) ----------------
    def values(self, model_name=None, **kwargs) -> pd.DataFrame:
        """Aligned value table: each numeric field as ``x`` / ``reference_x`` /
        ``x_diff`` on the cells shared with the reference model.

        ``model_name`` (the first positional argument) restricts the table to one
        compared model; ``per`` / ``layer`` / ``cells`` filter as on the group
        input accessor. With no ``model_name`` every non-reference model is
        returned.
        """

        return self._inputs.compare(model_name=model_name, **kwargs)

    # -- structural tier (the new set-difference engine) ----------------------
    def cells(
        self,
        *,
        model_name=None,
        per=None,
        layer=None,
        include_shared: bool = False,
    ) -> pd.DataFrame:
        """Structural set difference of applied cells vs the reference model.

        Returns a tidy frame with columns ``model, per, layer, cell,
        membership`` where ``membership`` is ``only_in_reference`` (the cell is
        in the reference but not this model), ``only_in_model`` (in this model
        but not the reference), or -- when ``include_shared=True`` --
        ``shared``. This is the difference a value-only compare cannot show,
        because its inner join drops cells present in only one model.
        """

        reference = self.group.reference
        _, ref_keys = _package_keys(
            self.group.models[reference], self.package_name, per=per, layer=layer
        )
        rows: list[tuple] = []
        for name in self._targets(model_name):
            _, model_keys = _package_keys(
                self.group.models[name], self.package_name, per=per, layer=layer
            )
            for key in sorted(ref_keys - model_keys):
                rows.append((name, *key, "only_in_reference"))
            for key in sorted(model_keys - ref_keys):
                rows.append((name, *key, "only_in_model"))
            if include_shared:
                for key in sorted(ref_keys & model_keys):
                    rows.append((name, *key, "shared"))
        return pd.DataFrame(rows, columns=["model", "per", "layer", "cell", "membership"])

    def summary(self, *, model_name=None) -> pd.DataFrame:
        """Per-model structural + value difference counts for this package."""

        reference = self.group.reference
        ref_present, ref_keys = _package_keys(self.group.models[reference], self.package_name)
        rows = []
        for name in self._targets(model_name):
            model_present, model_keys = _package_keys(self.group.models[name], self.package_name)
            only_ref = len(ref_keys - model_keys)
            only_model = len(model_keys - ref_keys)
            shared = len(ref_keys & model_keys)
            changed = self._value_cells_changed(name) if shared else 0
            identical = (
                ref_present == model_present
                and only_ref == 0
                and only_model == 0
                and changed == 0
            )
            rows.append(
                {
                    "package": self.package_name,
                    "model": name,
                    "present_in_reference": ref_present,
                    "present_in_model": model_present,
                    "cells_only_in_reference": only_ref,
                    "cells_only_in_model": only_model,
                    "cells_shared": shared,
                    "value_cells_changed": changed,
                    "identical": identical,
                }
            )
        return pd.DataFrame(rows, columns=_SUMMARY_COLUMNS)

    def _value_cells_changed(self, model_name: str) -> int:
        """Count shared cells whose numeric fields differ from the reference."""

        try:
            comparison = self._inputs.compare(model_name=model_name)
        except PACKAGE_TABLE_UNAVAILABLE:
            # Returning 0 here reads as "no values changed", so this handler
            # firing at all is a degradation worth seeing in a debug log: it is
            # what made a group with one package-less member report every other
            # member as identical to the reference (now fixed in
            # `GroupPackageInputs.get`, which skips such members).
            logger.debug(
                "cannot compare %s values for %s; reporting no value changes",
                self.package_name, model_name, exc_info=True,
            )
            return 0
        if comparison.empty:
            return 0
        diff_columns = [column for column in comparison.columns if column.endswith("_diff")]
        if not diff_columns:
            return 0
        nonzero = (comparison[diff_columns].fillna(0.0) != 0.0).any(axis=1)
        return int(nonzero.sum())


class ConnectionDiff:
    """Connection/reach geometry set difference for an advanced package.

    LAK connections have no stable per-connection key (a cell carries one
    vertical plus N horizontal connections, and their order is not guaranteed),
    so the comparison is a **multiset** difference of full geometry tuples: a
    connection whose geometry changed appears as one ``only_in_reference`` row
    (the old geometry) and one ``only_in_model`` row (the new). SFR reaches are
    handled the same way, keyed by reach number plus geometry. Float geometry is
    rounded (``round_to``) so build noise does not read as a difference.
    """

    _package_name: str = ""
    _identity_columns: tuple = ()
    _row_label: str = "connections"

    def __init__(self, diff: ModelDiff):
        """Bind a connection/reach geometry diff to a group via its ``ModelDiff``."""

        self._diff = diff
        self.group = diff.group

    @staticmethod
    def _build_table(model):  # pragma: no cover - overridden
        """Return one model's connection/reach geometry table (subclass override)."""

        raise NotImplementedError

    def _targets(self, model_name=None) -> list[str]:
        """The non-reference model names to diff (or just ``model_name``, validated)."""

        if model_name is not None:
            name = str(model_name)
            if name not in self.group.models:
                raise KeyError(f"Model {name!r} is not in the group.")
            if name == self.group.reference:
                raise ValueError("The reference model cannot be diffed against itself.")
            return [name]
        return [name for name in self.group.models if name != self.group.reference]

    def _counter(self, model_name, *, round_to):
        """Return ``(present, Counter-of-geometry-tuples, columns)`` for a model."""

        model = self.group.models[model_name]
        try:
            table = type(self)._build_table(model)
        except PACKAGE_TABLE_UNAVAILABLE:
            logger.debug(
                "no readable %s geometry on %s; comparing it as having none",
                self._package_name, model_name, exc_info=True,
            )
            return _package_present(model, self._package_name), Counter(), ()
        if table is None or table.empty:
            return _package_present(model, self._package_name), Counter(), ()
        columns = tuple(col for col in self._identity_columns if col in table.columns)
        frame = table.loc[:, list(columns)].copy()
        for column in columns:
            if pd.api.types.is_float_dtype(frame[column]):
                frame[column] = frame[column].round(round_to)
        tuples = [tuple(row) for row in frame.itertuples(index=False, name=None)]
        return True, Counter(tuples), columns

    def rows(self, *, model_name=None, round_to: int = 6, include_shared: bool = False):
        """Return the connections/reaches that differ from the reference network.

        Columns: ``model``, the geometry columns, and ``membership``
        (``only_in_reference`` / ``only_in_model`` / ``shared``).
        """

        _, ref_counter, ref_columns = self._counter(self.group.reference, round_to=round_to)
        collected: list[tuple] = []
        columns = ref_columns
        for name in self._targets(model_name):
            _, counter, model_columns = self._counter(name, round_to=round_to)
            if not columns:
                columns = model_columns
            for tup, count in sorted((ref_counter - counter).items()):
                collected.extend([(name, *tup, "only_in_reference")] * count)
            for tup, count in sorted((counter - ref_counter).items()):
                collected.extend([(name, *tup, "only_in_model")] * count)
            if include_shared:
                for tup, count in sorted((ref_counter & counter).items()):
                    collected.extend([(name, *tup, "shared")] * count)
        resolved = list(columns) if columns else list(self._identity_columns)
        return pd.DataFrame(collected, columns=["model", *resolved, "membership"])

    def summary(self, *, model_name=None, round_to: int = 6) -> pd.DataFrame:
        """Per-model connection/reach difference counts vs the reference."""

        ref_present, ref_counter, _ = self._counter(self.group.reference, round_to=round_to)
        rows = []
        for name in self._targets(model_name):
            present, counter, _ = self._counter(name, round_to=round_to)
            only_ref = int(sum((ref_counter - counter).values()))
            only_model = int(sum((counter - ref_counter).values()))
            shared = int(sum((ref_counter & counter).values()))
            identical = (ref_present == present) and only_ref == 0 and only_model == 0
            rows.append(
                {
                    "package": self._package_name,
                    "model": name,
                    "present_in_reference": ref_present,
                    "present_in_model": present,
                    "only_in_reference": only_ref,
                    "only_in_model": only_model,
                    "shared": shared,
                    "identical": identical,
                }
            )
        return pd.DataFrame(
            rows,
            columns=[
                "package", "model", "present_in_reference", "present_in_model",
                "only_in_reference", "only_in_model", "shared", "identical",
            ],
        )


class LakConnectionDiff(ConnectionDiff):
    """LAK lake-connection geometry difference (``diff.packages.lak``)."""

    _package_name = "lak"
    _identity_columns = _LAK_IDENTITY
    _row_label = "connections"

    @staticmethod
    def _build_table(model):
        """Return the model's LAK connection-geometry table."""

        return build_lak_connection_table(model)

    def connections(self, **kwargs) -> pd.DataFrame:
        """Per-connection geometry set difference vs the reference lake network."""

        return self.rows(**kwargs)


class SfrReachDiff(ConnectionDiff):
    """SFR reach geometry difference (``diff.packages.sfr``)."""

    _package_name = "sfr"
    _identity_columns = _SFR_IDENTITY
    _row_label = "reaches"

    @staticmethod
    def _build_table(model):
        """Return the model's SFR reach-geometry table."""

        return build_sfr_reach_table(model)

    def reaches(self, **kwargs) -> pd.DataFrame:
        """Per-reach geometry set difference vs the reference stream network."""

        return self.rows(**kwargs)

    # A reach is a stream's connection to a cell; expose both names.
    connections = reaches


def _resolve_package_diff(diff: ModelDiff, name: str):
    """Return the *inputs-tier* diff accessor for a package name.

    Internal: this backs ``summary()``/``report()``, which aggregate the
    setup-tier differences. The public tree exposes the same accessors as
    ``diff.packages.<pkg>.inputs``.
    """

    pkg = str(name).lower()
    if pkg in _DIFF_PACKAGES:
        return PackageDiff(diff, pkg)
    if pkg == "lak":
        return LakConnectionDiff(diff)
    if pkg == "sfr":
        return SfrReachDiff(diff)
    raise AttributeError(
        f"ModelDiff does not diff package {name!r}; supported: "
        f"{', '.join((*_DIFF_PACKAGES, *_CONNECTION_PACKAGES))}."
    )


class _BcPackageDiffNode:
    """``diff.packages.<pkg>`` for a BC package: ``.inputs`` / ``.results``.

    Mirrors the single-model/group tree shape -- the declared-data difference
    lives under ``inputs`` and the computed-output difference under ``results``.
    """

    def __init__(self, diff: ModelDiff, package_name: str):
        """Bind a BC package's ``.inputs``/``.results`` diff node to ``package_name``."""

        self._diff = diff
        self.package_name = str(package_name).lower()

    @property
    def inputs(self) -> PackageDiff:
        """Declared stress-period data difference (cells + values)."""

        return PackageDiff(self._diff, self.package_name)

    @property
    def results(self) -> CellResultsDiffNamespace:
        """Computed cell-budget difference (field ``q``); requires runs."""

        return CellResultsDiffNamespace(self._diff, self.package_name)

    def __repr__(self) -> str:
        """Show the package name and the two sub-nodes available."""

        return f"<diff.packages.{self.package_name}: .inputs / .results>"


class _LakDiffNode:
    """``diff.packages.lak``: connection-geometry inputs + q/stage results."""

    def __init__(self, diff: ModelDiff):
        """Bind the LAK ``.inputs``/``.results`` diff node to a group's ``ModelDiff``."""

        self._diff = diff
        self.package_name = "lak"

    @property
    def inputs(self) -> LakConnectionDiff:
        """Lake connection-geometry difference (the LAK input tier)."""

        return LakConnectionDiff(self._diff)

    @property
    def results(self) -> LakResultsDiffNamespace:
        """Computed lake results difference -- fields ``q`` and ``stage``."""

        return LakResultsDiffNamespace(self._diff)

    def __repr__(self) -> str:
        """Show the LAK diff node's available sub-nodes."""

        return "<diff.packages.lak: .inputs / .results>"


class _SfrDiffNode:
    """``diff.packages.sfr``: reach-geometry inputs + q/stage results."""

    def __init__(self, diff: ModelDiff):
        """Bind the SFR ``.inputs``/``.results`` diff node to a group's ``ModelDiff``."""

        self._diff = diff
        self.package_name = "sfr"

    @property
    def inputs(self) -> SfrReachDiff:
        """Stream reach-geometry difference (the SFR input tier)."""

        return SfrReachDiff(self._diff)

    @property
    def results(self) -> SfrResultsDiffNamespace:
        """Computed stream results difference -- fields ``q`` and ``stage``."""

        return SfrResultsDiffNamespace(self._diff)

    def __repr__(self) -> str:
        """Show the SFR diff node's available sub-nodes."""

        return "<diff.packages.sfr: .inputs / .results>"


class _UzfDiffNode:
    """``diff.packages.uzf``: results only (UZF input diffing is not built)."""

    def __init__(self, diff: ModelDiff):
        """Bind the UZF ``.results``-only diff node to a group's ``ModelDiff``."""

        self._diff = diff
        self.package_name = "uzf"

    @property
    def results(self) -> UzfResultsDiffNamespace:
        """Computed UZF results difference -- fields ``gwrch`` and ``sat``."""

        return UzfResultsDiffNamespace(self._diff)

    def __repr__(self) -> str:
        """Show the UZF diff node's available sub-node."""

        return "<diff.packages.uzf: .results>"


class _MvrDiffNode:
    """``diff.packages.mvr``: results only (mover options diff in ``config``)."""

    def __init__(self, diff: ModelDiff):
        """Bind the MVR ``.results``-only diff node to a group's ``ModelDiff``."""

        self._diff = diff
        self.package_name = "mvr"

    @property
    def results(self) -> MvrResultDiff:
        """Mover-flow difference per moved package and direction."""

        return MvrResultDiff(self._diff)

    def __repr__(self) -> str:
        """Show the MVR diff node's available sub-node."""

        return "<diff.packages.mvr: .results>"


class _PackageDiffNamespace:
    """Attribute access ``diff.packages.<package>`` -> a package diff node.

    Every node mirrors the single-model/group tree: ``.inputs`` for declared
    data, ``.results`` for computed outputs. The supported packages are
    declared as explicit properties (below) so IDEs / static analyzers can
    discover them; ``__getattr__`` remains as a fallback that raises a helpful
    error for unsupported package names.
    """

    def __init__(self, diff: ModelDiff):
        """Bind the ``diff.packages`` namespace to a group's ``ModelDiff``."""

        self._diff = diff

    # -- explicit, IDE-discoverable package accessors -------------------------
    @property
    def rch(self) -> _BcPackageDiffNode:
        """Recharge (RCH) difference node."""
        return _BcPackageDiffNode(self._diff, "rch")

    @property
    def chd(self) -> _BcPackageDiffNode:
        """Constant-head (CHD) difference node."""
        return _BcPackageDiffNode(self._diff, "chd")

    @property
    def drn(self) -> _BcPackageDiffNode:
        """Drain (DRN) difference node."""
        return _BcPackageDiffNode(self._diff, "drn")

    @property
    def ghb(self) -> _BcPackageDiffNode:
        """General-head-boundary (GHB) difference node."""
        return _BcPackageDiffNode(self._diff, "ghb")

    @property
    def riv(self) -> _BcPackageDiffNode:
        """River (RIV) difference node."""
        return _BcPackageDiffNode(self._diff, "riv")

    @property
    def wel(self) -> _BcPackageDiffNode:
        """Well (WEL) difference node."""
        return _BcPackageDiffNode(self._diff, "wel")

    @property
    def evt(self) -> _BcPackageDiffNode:
        """Evapotranspiration (EVT) difference node."""
        return _BcPackageDiffNode(self._diff, "evt")

    @property
    def lak(self) -> _LakDiffNode:
        """Lake (LAK) difference node."""
        return _LakDiffNode(self._diff)

    @property
    def sfr(self) -> _SfrDiffNode:
        """Stream (SFR) difference node."""
        return _SfrDiffNode(self._diff)

    @property
    def uzf(self) -> _UzfDiffNode:
        """Unsaturated-zone (UZF) difference node."""
        return _UzfDiffNode(self._diff)

    @property
    def mvr(self) -> _MvrDiffNode:
        """Mover (MVR) difference node."""
        return _MvrDiffNode(self._diff)

    def __getattr__(self, name: str):
        """Fallback package lookup -> a BC diff node, or a helpful ``AttributeError``."""

        # Fallback for any package name not declared above -> helpful AttributeError.
        pkg = str(name).lower()
        if pkg in _DIFF_PACKAGES:
            return _BcPackageDiffNode(self._diff, pkg)
        raise AttributeError(
            f"ModelDiff does not diff package {name!r}; supported: "
            f"{', '.join((*_DIFF_PACKAGES, *_CONNECTION_PACKAGES, 'uzf', 'mvr'))}."
        )

    def __dir__(self):
        """Advertise the diffable package names for tab-completion."""

        known = (
            set(self._diff.package_names)
            | set(self._diff.connection_package_names)
            | {"uzf", "mvr"}
        )
        return sorted(set(super().__dir__()) | known)


class _FocusedModelDiff:
    """A :class:`ModelDiff` narrowed to a single non-reference model."""

    def __init__(self, diff: ModelDiff, model_name: str):
        """Narrow a ``ModelDiff`` to one non-reference model (validated, not the reference)."""

        self._diff = diff
        self.model_name = str(model_name)
        if self.model_name not in diff.group.models:
            raise KeyError(f"Model {self.model_name!r} is not in the group.")
        if self.model_name == diff.group.reference:
            raise ValueError("The reference model cannot be diffed against itself.")

    def summary(self) -> pd.DataFrame:
        """Per-package difference counts for this one model vs the reference."""

        return self._diff.summary(model_name=self.model_name)

    def report(self) -> str:
        """The Markdown faithful-copy report narrowed to this one model."""

        return self._diff._render_report(model_name=self.model_name)

    def cells(self, package: str, **kwargs) -> pd.DataFrame:
        """Structural cell set-difference for one ``package`` on this model."""

        return _resolve_package_diff(self._diff, package).cells(
            model_name=self.model_name, **kwargs
        )

    def values(self, package: str, **kwargs) -> pd.DataFrame:
        """Aligned value difference for one ``package`` on this model."""

        return _resolve_package_diff(self._diff, package).values(
            model_name=self.model_name, **kwargs
        )


_CONFIG_ABSENT = "<absent>"


def _config_values_equal(left, right) -> bool:
    """Compare two normalized config values robustly."""

    if left is right:
        return True
    try:
        return bool(left == right)
    except (ValueError, TypeError):  # pragma: no cover - defensive
        # A numpy array compares elementwise, so `bool(...)` on the result
        # raises ValueError ("truth value of an array is ambiguous"); TypeError
        # covers types whose __eq__ refuses the comparison outright. Comparing
        # the reprs is the documented fallback.
        return repr(left) == repr(right)


class ConfigDiff:
    """Configuration-tier difference: tdis / ims / oc / package OPTIONS.

    Compares each model's normalized settings (``model.config``) against the
    reference model's, surfacing every setting whose value differs or is present
    in only one model. This answers "are these runs configured the same?" --
    solver block, timing (including per-period ``nstp``/``tsmult``), and package
    options -- the part a value/structural cell diff cannot see.
    """

    def __init__(self, diff: ModelDiff):
        """Bind the configuration-tier diff to a group via its ``ModelDiff``."""

        self._diff = diff
        self.group = diff.group

    def _targets(self, model_name=None) -> list[str]:
        """The non-reference model names to diff (or just ``model_name``, validated)."""

        if model_name is not None:
            name = str(model_name)
            if name not in self.group.models:
                raise KeyError(f"Model {name!r} is not in the group.")
            if name == self.group.reference:
                raise ValueError("The reference model cannot be diffed against itself.")
            return [name]
        return [name for name in self.group.models if name != self.group.reference]

    def _settings_map(self, model_name: str) -> dict:
        """A ``{(section, setting): value}`` map of one model's normalized config."""

        frame = self.group.models[model_name].config.settings()
        return {
            (row.section, row.setting): row.value
            for row in frame.itertuples(index=False)
        }

    def settings(self, *, model_name=None, section=None) -> pd.DataFrame:
        """Return the settings that differ from the reference model.

        Columns: ``model, section, setting, reference_value, model_value``. A
        setting present in only one model shows ``'<absent>'`` on the other side.
        Optionally restrict to one ``section`` (e.g. ``"ims"`` or ``"tdis"``).
        """

        wanted_section = None if section is None else str(section).lower()
        reference_map = self._settings_map(self.group.reference)
        rows = []
        for name in self._targets(model_name):
            model_map = self._settings_map(name)
            for key in sorted(set(reference_map) | set(model_map)):
                section_name, setting = key
                if wanted_section is not None and section_name != wanted_section:
                    continue
                ref_value = reference_map.get(key, _CONFIG_ABSENT)
                model_value = model_map.get(key, _CONFIG_ABSENT)
                if _config_values_equal(ref_value, model_value):
                    continue
                rows.append(
                    {
                        "model": name,
                        "section": section_name,
                        "setting": setting,
                        "reference_value": ref_value,
                        "model_value": model_value,
                    }
                )
        return pd.DataFrame(
            rows,
            columns=["model", "section", "setting", "reference_value", "model_value"],
        )

    def summary(self, *, model_name=None) -> pd.DataFrame:
        """Per-model count of differing configuration settings."""

        differences = self.settings(model_name=model_name)
        rows = []
        for name in self._targets(model_name):
            count = 0 if differences.empty else int((differences["model"] == name).sum())
            rows.append(
                {"model": name, "settings_differing": count, "identical": count == 0}
            )
        return pd.DataFrame(rows, columns=["model", "settings_differing", "identical"])


class ModelDiff:
    """Reference-star difference across a :class:`ModelGroup` (the ``diff`` verb)."""

    def __init__(self, group: ModelGroup):
        """Wrap a :class:`ModelGroup` as the ``diff`` engine (its reference is the baseline)."""

        self.group = group

    @property
    def reference(self) -> str:
        """The baseline model every other model is compared against."""

        return self.group.reference

    @property
    def model_names(self) -> list[str]:
        """The non-reference model names being compared to the reference."""

        return [name for name in self.group.models if name != self.group.reference]

    @property
    def package_names(self) -> list[str]:
        """BC packages present in at least one model in the group."""

        return [
            pkg
            for pkg in _DIFF_PACKAGES
            if any(_package_present(model, pkg) for model in self.group.models.values())
        ]

    @property
    def connection_package_names(self) -> list[str]:
        """Connection-based packages (lak/sfr) present in at least one model."""

        return [
            pkg
            for pkg in _CONNECTION_PACKAGES
            if any(_package_present(model, pkg) for model in self.group.models.values())
        ]

    @property
    def packages(self) -> _PackageDiffNamespace:
        """Namespace for per-package diffs: ``diff.packages.ghb.cells()`` etc."""

        return _PackageDiffNamespace(self)

    @property
    def config(self) -> ConfigDiff:
        """Configuration-tier diff (tdis / ims / oc / package options)."""

        return ConfigDiff(self)

    @property
    def hds(self) -> HeadsResultDiff:
        """Head-difference leaf (Δhead maps/plots/sections vs the reference).

        Mirrors ``model.hds`` / ``group.hds``; requires completed runs.
        """

        return HeadsResultDiff(self)

    @property
    def conc(self) -> ConcResultDiff:
        """Concentration-difference leaf (Δconc maps/plots vs the reference).

        Mirrors ``model.conc`` / ``group.conc``; every member must be a GWT model.
        """

        return ConcResultDiff(self)

    @property
    def temp(self) -> TempResultDiff:
        """Temperature-difference leaf (Δtemp maps/plots vs the reference).

        Mirrors ``model.temp`` / ``group.temp``; every member must be a GWE model.
        """

        return TempResultDiff(self)

    @property
    def bud(self) -> BudgetResultDiff:
        """Volumetric (listing) budget difference per term vs the reference."""

        return BudgetResultDiff(self)

    def package(self, name: str):
        """Return the diff node for one package -- same as ``diff.packages.<name>``.

        Every node has ``.inputs`` (declared-data difference) and, where the
        package produces cell output, ``.results`` (computed difference).
        """

        return getattr(self.packages, str(name).lower())

    def model(self, name: str) -> _FocusedModelDiff:
        """Focus the diff on a single non-reference model."""

        return _FocusedModelDiff(self, name)

    def summary(self, *, model_name=None) -> pd.DataFrame:
        """Package x model matrix of structural + value difference counts.

        One row per (package, non-reference model). ``identical`` is ``True``
        when a package matches the reference in presence, cell membership, and
        all field values.
        """

        frames = [
            _resolve_package_diff(self, pkg).summary(model_name=model_name)
            for pkg in self.package_names
        ]
        frames = [frame for frame in frames if not frame.empty]
        if not frames:
            return pd.DataFrame(columns=_SUMMARY_COLUMNS)
        return pd.concat(frames, ignore_index=True)

    def report(self, *, results: bool = False) -> str:
        """Return a readable Markdown 'faithful-copy' report.

        By default this covers how the models are *set up* (packages, config,
        connection geometry). Pass ``results=True`` to also compare computed
        outputs (heads, budget) -- this reads output files and requires runs.
        """

        return self._render_report(results=results)

    def _render_report(self, *, model_name=None, results: bool = False) -> str:
        """Assemble the Markdown diff report (package/config/connection, optional results).

        Backs both :meth:`report` (all non-reference models) and the focused
        report (one ``model_name``); each model is flagged identical or its
        differing tiers are tabulated.
        """

        reference = self.group.reference
        summary = self.summary(model_name=model_name)
        config = self.config
        targets = [model_name] if model_name is not None else self.model_names
        lines = [f"# Model diff -- reference: `{reference}`", ""]
        for name in targets:
            block = summary[summary["model"] == name]
            config_diffs = config.settings(model_name=name)
            connection_rows = []
            for pkg in self.connection_package_names:
                connection_summary = _resolve_package_diff(self, pkg).summary(model_name=name)
                if not connection_summary.empty and not bool(
                    connection_summary.iloc[0]["identical"]
                ):
                    connection_rows.append(connection_summary.iloc[0])
            results_block = self._results_report_block(name) if results else None
            package_identical = block.empty or bool(block["identical"].all())
            config_identical = config_diffs.empty
            connections_identical = len(connection_rows) == 0
            results_identical = results_block is None or results_block["identical"]
            if (
                package_identical
                and config_identical
                and connections_identical
                and results_identical
            ):
                lines.append(f"## `{name}` -- identical to reference")
                if results_block is not None and results_block["lines"]:
                    lines.append("")
                    lines.extend(results_block["lines"])
                lines.append("")
                continue
            lines.append(f"## `{name}` -- differs from reference")
            lines.append("")
            if not package_identical:
                lines.append("Package differences:")
                lines.append("")
                lines.append(
                    "| package | present (ref/model) | only-ref cells | only-model cells "
                    "| shared | value cells changed |"
                )
                lines.append("|---|---|---|---|---|---|")
                for _, row in block.iterrows():
                    lines.append(
                        f"| {row['package']} "
                        f"| {row['present_in_reference']}/{row['present_in_model']} "
                        f"| {row['cells_only_in_reference']} "
                        f"| {row['cells_only_in_model']} "
                        f"| {row['cells_shared']} "
                        f"| {row['value_cells_changed']} |"
                    )
                lines.append("")
            if not config_identical:
                lines.append("Configuration differences:")
                lines.append("")
                lines.append("| section | setting | reference | model |")
                lines.append("|---|---|---|---|")
                for _, row in config_diffs.iterrows():
                    lines.append(
                        f"| {row['section']} | {row['setting']} "
                        f"| {row['reference_value']} | {row['model_value']} |"
                    )
                lines.append("")
            if connection_rows:
                lines.append("Connection differences:")
                lines.append("")
                lines.append(
                    "| package | present (ref/model) | only-ref | only-model | shared |"
                )
                lines.append("|---|---|---|---|---|")
                for row in connection_rows:
                    lines.append(
                        f"| {row['package']} "
                        f"| {row['present_in_reference']}/{row['present_in_model']} "
                        f"| {row['only_in_reference']} "
                        f"| {row['only_in_model']} "
                        f"| {row['shared']} |"
                    )
                lines.append("")
            if results_block is not None and results_block["lines"]:
                lines.extend(results_block["lines"])
                lines.append("")
        return "\n".join(lines)

    def _results_report_block(self, model_name: str) -> dict:
        """Build the optional 'Results differences' report lines for one model.

        Returns ``{"identical": bool, "lines": [...]}``. Results require completed
        runs; if outputs are missing the block notes that and does not affect the
        identical decision.
        """

        try:
            heads = self.hds.summary(model_name=model_name)
            budget = self.bud.summary(model_name=model_name)
        except _RESULTS_UNAVAILABLE:
            logger.debug(
                "no comparable outputs for %s", model_name, exc_info=True,
            )
            return {"identical": True, "lines": ["_Results differences: outputs unavailable (models not run)._"]}

        heads_ok = heads.empty or bool(heads["within_tolerance"].all())
        budget_ok = budget.empty or bool(budget["within_tolerance"].all())
        lines = ["Results differences:", ""]
        lines.append("| output | within tolerance | detail |")
        lines.append("|---|---|---|")
        if not heads.empty:
            row = heads.iloc[0]
            lines.append(
                f"| heads | {bool(row['within_tolerance'])} "
                f"| max|diff|={row['max_abs_diff']:.4g} at cell {row['argmax_cell']} "
                f"layer {row['argmax_layer']} kstpkper {row['argmax_kstpkper']} |"
            )
        breached = budget[~budget["within_tolerance"]] if not budget.empty else budget
        if not budget.empty:
            if breached.empty:
                lines.append("| budget | True | all terms within tolerance |")
            else:
                worst = breached.reindex(
                    breached["max_abs_diff"].abs().sort_values(ascending=False).index
                ).iloc[0]
                lines.append(
                    f"| budget | False "
                    f"| {worst['term']} diff_total={worst['diff_total']:.4g} "
                    f"({worst['pct_change']:.1f}%) |"
                )
        return {"identical": heads_ok and budget_ok, "lines": lines}

    def __repr__(self) -> str:
        """Show the reference and the compared model names."""

        return (
            f"ModelDiff(reference={self.group.reference!r}, "
            f"models={self.model_names!r})"
        )
