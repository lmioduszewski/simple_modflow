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
from myflopy.modflow.mf6.package_explorer_utils import period_number
from myflopy.modflow.mf6.package_registry import (
    _PACKAGE_EXPLORER_SPECS,
    input_field_agg,
    model_packages_of_type,
)
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
# Every package type ``diff.packages`` answers -- what ``ModelDiff.package(name)``
# accepts, so a namespace attribute name (``"_diff"``) is not taken for one.
_DIFFED_PACKAGE_TYPES: tuple[str, ...] = (*_DIFF_PACKAGES, *_CONNECTION_PACKAGES, "uzf", "mvr")
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
    """Return whether ``model`` carries a package of TYPE ``package_name``.

    Found by TYPE through ``package_registry.model_packages_of_type`` (ledger 352's
    one rule), never by name: a DRN named ``seeps`` is a DRN, and a GHB named
    ``lake_edge`` is not a LAK. It matched a NAME prefix before, so such a DRN was
    never compared and ``report()`` read "identical to reference".
    """

    return bool(model_packages_of_type(model, package_name))


def _package_keys(model, package_name: str, *, per=None, layer=None):
    """Return ``(present, keys)`` for a BC package on ``model``.

    ``keys`` is a set of zero-based ``(per, layer, cell)`` tuples identifying
    every cell the package applies to (optionally filtered by ``per``/``layer``).
    A package that is absent, empty, or unreadable yields an empty set; the
    ``present`` flag distinguishes "attached but no cells" from "not attached".
    ``per`` reads the records MF6 runs in that period; a mistyped one is refused
    HERE, before the handler below would read its ``TypeError`` as "no package"
    (ledger 449).
    """

    per = None if per is None else period_number(per)
    try:
        table = build_cell_package_input_table(model, package_name, per=per, layer=layer)
    except PACKAGE_TABLE_UNAVAILABLE:
        # An array-form RCHA/EVTA used to be unreadable here while `_package_present`
        # (a PREFIX match) reported it present, so two such models compared as
        # identical. The table reads every package of the type now, array form
        # included (ledger 229), so its cells are keys like any other.
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

    def _series_default_agg(self) -> str:
        """How ``plot()`` combines a period's differences: by the field's kind (ledger 472).

        The SUM of a conductance's differences is the change in the total; an
        elevation's are averaged. Every field was summed before.
        """

        field = self.field_name or get_default_package_value_column(self.package_name)
        return input_field_agg(self.package_name, field)

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

        The group's input compare (``group.packages.<pkg>.inputs.compare``): every
        member's stress-period rows of this package -- every package of the type,
        each row tagged ``pname`` -- inner-joined to the reference's on ``per``,
        ``layer``, ``cell`` and every non-numeric column (``package``, ``pname``,
        a boundname), so a cell in only one model is absent here (:meth:`cells`
        shows it) and so is a row whose package NAME differs. A cell holding
        several records of one package pairs them in file order, one row per
        record (it paired each with each until ledger 472: two records on
        each side gave four rows). ``x_diff`` is model minus reference, in the
        field's own units.

        Parameters
        ----------
        model_name : str or None, optional
            Keep one compared model, by its name in the group (the first
            positional argument). Default ``None`` returns every non-reference
            model. It is a FILTER, not checked: a name not in the group, or the
            reference's, gives an empty frame where :meth:`cells` and
            :meth:`summary` raise.
        **kwargs
            Forwarded to the group compare: ``per`` (int, a zero-based stress
            period: the records MF6 runs in it, the last block at or before it;
            ledger 449),
            ``layer`` (int or sequence of int, zero-based) and ``cells`` (int or
            sequence of int, zero-based in-layer), combined with AND; default
            ``None`` keeps everything. Anything else raises ``TypeError``. The
            pinned field (``inputs.cond``) does not narrow the table: every numeric
            field is diffed.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``reference_model``, the join columns, then ``x``,
            ``reference_x`` and ``x_diff`` for each numeric field (``cond``,
            ``bhead``, ...). Empty (no columns) when no member holds rows the
            filters keep.
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

        A cell is its ``(per, layer, cell)`` key, pooled over every package of the
        type (two CHD packages on one cell are one key), so a cell moved from one
        package of the type to another is not a difference here. A package a
        member does not carry, or cannot be read, contributes no keys.

        Parameters
        ----------
        model_name : str or None, optional
            Compare one model, by its name in the group. Default ``None`` compares
            every non-reference model. A name not in the group raises ``KeyError``;
            the reference's raises ``ValueError``.
        per : int or None, optional
            Keep one zero-based stress period: the cells MF6 runs in it in each
            model -- the last block at or before ``per``, none under an empty
            block (ledger 449). Default ``None`` uses every block as written.
        layer : int or sequence of int or None, optional
            Keep these zero-based layers. Default ``None`` keeps every layer.
        include_shared : bool, optional
            Default ``False`` lists only the cells in one model and not the other;
            ``True`` adds the cells in both, as ``membership="shared"``.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``per``, ``layer``, ``cell`` (zero-based ints) and
            ``membership``, sorted by key within each membership.
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
        """Per-model structural + value difference counts for this package.

        Over every period and layer: the cell keys of :meth:`cells`, and the
        shared rows of :meth:`values` whose numeric fields differ at all (no
        tolerance; a ``NaN`` on either side counts as no change).

        Parameters
        ----------
        model_name : str or None, optional
            Summarize one model, by its name in the group. Default ``None``
            summarizes every non-reference model. A name not in the group raises
            ``KeyError``; the reference's raises ``ValueError``.

        Returns
        -------
        pandas.DataFrame
            One row per compared model: ``package``, ``model``,
            ``present_in_reference``, ``present_in_model``,
            ``cells_only_in_reference``, ``cells_only_in_model``, ``cells_shared``,
            ``value_cells_changed`` (the shared ``(per, layer, cell)`` keys, as
            ``cells_shared`` counts them, with any record whose value differs --
            pooled over the type's packages and its records, so one cell counts
            once per period; it counted ROWS of :meth:`values` until ledger
            472, so a cell of two records a side, both changed, read 4) and
            ``identical`` (presence, cells and values all agree).
        """

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
        # CELLS, as `cells_shared` counts them -- a (per, layer, cell) key -- not
        # rows: a cell of two changed records counted 2, and 4 when records
        # cross-paired (ledger 472).
        keys = [column for column in ("per", "layer", "cell") if column in comparison.columns]
        return int(len(comparison.loc[nonzero, keys].drop_duplicates()))


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

        A multiset difference of whole geometry rows -- LAK ``lake``, ``layer``,
        ``cell``, ``claktype``, ``belev``, ``telev``, ``connlen``, ``connwidth``;
        SFR ``reach``, ``layer``, ``cell``, ``rlen`` (zero-based ids, model length
        units) -- so a connection whose geometry changed is one
        ``only_in_reference`` row and one ``only_in_model`` row. A member without
        the package, or whose table cannot be read, counts as holding no rows, so
        every row of the other side is listed as one-sided.

        Parameters
        ----------
        model_name : str or None, optional
            Compare one model, by its name in the group. Default ``None`` compares
            every non-reference model. A name not in the group raises ``KeyError``;
            the reference's raises ``ValueError``.
        round_to : int, optional
            Decimal places every FLOAT geometry column is rounded to on both sides
            before rows are compared (``pandas.Series.round``; a negative value
            rounds to tens, hundreds, ...), so build noise below it is not a
            difference. Integer and text columns are compared exactly. Default
            ``6``. The rows returned carry the ROUNDED values.
        include_shared : bool, optional
            Default ``False`` lists only the rows in one network and not the
            other; ``True`` adds the rows in both, as ``membership="shared"``.

        Returns
        -------
        pandas.DataFrame
            ``model``, the identity columns the tables carry, and
            ``membership``; a row repeated in a network is repeated here as often
            as its count differs.
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
        """Per-model connection/reach difference counts vs the reference.

        Parameters
        ----------
        model_name : str or None, optional
            Summarize one model, by its name in the group. Default ``None``
            summarizes every non-reference model. A name not in the group raises
            ``KeyError``; the reference's raises ``ValueError``.
        round_to : int, optional
            Decimal places the float geometry columns are rounded to before rows
            are matched, as in :meth:`rows`. Default ``6``.

        Returns
        -------
        pandas.DataFrame
            One row per compared model: ``package``, ``model``,
            ``present_in_reference``, ``present_in_model``, ``only_in_reference``,
            ``only_in_model``, ``shared`` (row counts) and ``identical``
            (presence agrees and no row is on one side only).
        """

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
        """Per-connection geometry set difference vs the reference lake network.

        Parameters
        ----------
        **kwargs
            Forwarded to :meth:`rows`: ``model_name``, ``round_to`` and
            ``include_shared``, keyword-only. Anything else raises ``TypeError``.

        Returns
        -------
        pandas.DataFrame
            :meth:`rows`' table: ``model``, the LAK connection-geometry columns and
            ``membership``.
        """

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
        """Per-reach geometry set difference vs the reference stream network.

        Also reachable as ``connections`` (a reach is a stream's connection to a
        cell).

        Parameters
        ----------
        **kwargs
            Forwarded to :meth:`rows`: ``model_name``, ``round_to`` and
            ``include_shared``, keyword-only. Anything else raises ``TypeError``.

        Returns
        -------
        pandas.DataFrame
            :meth:`rows`' table: ``model``, ``reach``, ``layer``, ``cell``,
            ``rlen`` and ``membership``.
        """

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
            f"{', '.join(_DIFFED_PACKAGE_TYPES)}."
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
        """Structural cell set-difference for one ``package`` on this model.

        :meth:`PackageDiff.cells` with ``model_name`` set to this model.

        Parameters
        ----------
        package : str
            A cell-based BC package type, any case: ``rch``, ``chd``, ``drn``,
            ``ghb``, ``riv``, ``wel`` or ``evt``. ``lak`` and ``sfr`` raise
            ``AttributeError`` naming no alternative (their diff has no ``cells``;
            ``diff.packages.lak.inputs.rows(model_name=...)`` is the connection
            difference), and any other name raises ``AttributeError`` listing the
            supported ones.
        **kwargs
            Forwarded to :meth:`PackageDiff.cells`: ``per``, ``layer`` and
            ``include_shared``. ``model_name`` is this model's and may not be
            passed again (``TypeError``); anything else raises ``TypeError``.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``per``, ``layer``, ``cell`` and ``membership``.
        """

        return _resolve_package_diff(self._diff, package).cells(
            model_name=self.model_name, **kwargs
        )

    def values(self, package: str, **kwargs) -> pd.DataFrame:
        """Aligned value difference for one ``package`` on this model.

        :meth:`PackageDiff.values` with ``model_name`` set to this model.

        Parameters
        ----------
        package : str
            A cell-based BC package type, any case, as in :meth:`cells`; ``lak``
            and ``sfr`` raise ``AttributeError`` (their diff has no ``values``),
            any other name ``AttributeError`` listing the supported ones.
        **kwargs
            Forwarded to :meth:`PackageDiff.values` and on to the group compare:
            ``per``, ``layer`` and ``cells``. ``model_name`` is this model's and may
            not be passed again (``TypeError``); anything else raises
            ``TypeError``.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``reference_model``, the join columns, and ``x`` /
            ``reference_x`` / ``x_diff`` for each numeric field.
        """

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

        Each member's ``model.config.settings()`` is compared key by key on
        ``(section, setting)``; values that compare equal (``==``, or the same
        ``repr`` where ``==`` cannot answer) are left out.

        Parameters
        ----------
        model_name : str or None, optional
            Compare one model, by its name in the group. Default ``None`` compares
            every non-reference model. A name not in the group raises ``KeyError``;
            the reference's raises ``ValueError``.
        section : str or None, optional
            Keep one section, any case: the lower-case names ``model.config``
            writes -- ``"tdis"``, ``"ims"``, ``"oc"``, ``"model"`` (the name-file
            options such as ``NEWTON``), and one per package type present
            (``"npf"``, ``"disv"``, ``"chd"``, ...). Default ``None`` keeps every
            section. A name no member has is NOT refused: it gives an empty frame,
            the same answer as "no difference".

        Returns
        -------
        pandas.DataFrame
            One row per differing setting per compared model, in ``(section,
            setting)`` order.
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
        """Per-model count of differing configuration settings.

        Parameters
        ----------
        model_name : str or None, optional
            Summarize one model, by its name in the group. Default ``None``
            summarizes every non-reference model. A name not in the group raises
            ``KeyError``; the reference's raises ``ValueError``.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``settings_differing`` (the rows :meth:`settings` gives it,
            every section) and ``identical`` (none differ).
        """

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
        """Wrap a :class:`ModelGroup` as the ``diff`` engine (its reference is the baseline).

        Parameters
        ----------
        group : ModelGroup
            The models to compare: its ``reference`` is the baseline and every
            other member is diffed against it. A group of ONE model is refused
            with ``ValueError`` naming the spelling with a second model: it has
            nothing to compare, and used to return an empty diff that read as "no
            differences" (ledger 373).
        """

        if len(group.models) < 2:
            only = next(iter(group.models), None)
            raise ValueError(
                f"A diff compares models against the reference {group.reference!r}, and "
                f"this group holds only {only!r}: give it a second model -- "
                f"ModelGroup({{{only!r}: m1, 'other': m2}}).diff(), or "
                "model.diff(other)."
            )
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
        """BC package types some model in the group carries, found by TYPE, whatever named."""

        return [
            pkg
            for pkg in _DIFF_PACKAGES
            if any(_package_present(model, pkg) for model in self.group.models.values())
        ]

    @property
    def connection_package_names(self) -> list[str]:
        """Connection-based package types (lak/sfr) some model carries, found by TYPE."""

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

        Parameters
        ----------
        name : str
            A package type, any case: ``rch``, ``chd``, ``drn``, ``ghb``, ``riv``,
            ``wel``, ``evt``, ``lak``, ``sfr``, ``uzf`` or ``mvr``. Any other name
            raises ``AttributeError`` listing those -- the namespace's own
            attribute names included (``"_diff"`` returned this ``ModelDiff``
            until 2026-10-11). Whether a member carries the package is not
            checked here.

        Returns
        -------
        object
            The node ``diff.packages.<name>`` returns: ``.inputs`` and
            ``.results`` for the BC packages, LAK and SFR, ``.results`` alone for
            UZF and MVR.
        """

        package = str(name).lower()
        if package not in _DIFFED_PACKAGE_TYPES:
            raise AttributeError(
                f"ModelDiff does not diff package {name!r}; supported: "
                f"{', '.join(_DIFFED_PACKAGE_TYPES)}."
            )
        return getattr(self.packages, package)

    def model(self, name: str) -> _FocusedModelDiff:
        """Focus the diff on a single non-reference model.

        Parameters
        ----------
        name : str
            The model's name in the group (converted with ``str``). A name not in
            the group raises ``KeyError``; the reference's raises ``ValueError``.

        Returns
        -------
        _FocusedModelDiff
            ``summary()``, ``report()``, ``cells(package)`` and
            ``values(package)``, each narrowed to this model.
        """

        return _FocusedModelDiff(self, name)

    def summary(self, *, model_name=None) -> pd.DataFrame:
        """Package x model matrix of structural + value difference counts.

        One row per (package, non-reference model). ``identical`` is ``True``
        when a package matches the reference in presence, cell membership, and
        all field values.

        The cell-based BC packages only (:attr:`package_names`), each through its
        :meth:`PackageDiff.summary`; LAK/SFR geometry, configuration and results
        are :meth:`report`'s. A type is listed when some member carries a package
        of that TYPE, whatever it is named (``package_registry.
        model_packages_of_type``, ledger 352): a DRN named ``seeps`` is a DRN. It
        was found by a NAME prefix until 2026-10-11, so that DRN was left out --
        measured, a conductance changed 10 to 30 there and :meth:`report` read
        "identical to reference" -- and a GHB named ``lake_edge`` counted as a
        LAK (``docs/reproducers/diff_results_tier.py``).

        Parameters
        ----------
        model_name : str or None, optional
            Summarize one model, by its name in the group. Default ``None``
            summarizes every non-reference model. A name not in the group raises
            ``KeyError`` and the reference's ``ValueError`` -- but only when
            :attr:`package_names` lists a type; with none, the empty table comes
            back whatever it says.

        Returns
        -------
        pandas.DataFrame
            :meth:`PackageDiff.summary`'s columns, one row per listed package and
            compared model; empty when :attr:`package_names` is.
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

        Parameters
        ----------
        results : bool, optional
            Default ``False`` reports packages (:meth:`summary`), configuration
            (``config.settings()``) and LAK/SFR connection geometry -- the package
            and geometry tiers finding a package by its TYPE, whatever it is named
            (:meth:`summary`). ``True`` adds,
            per model, ``hds.summary()`` and ``bud.summary()`` at their default
            tolerances (``atol`` = ``rtol`` = 1e-3) as a "Results differences"
            table, and a breach makes the model "differs" -- a cell dry or
            inactive in one model only is one (``n_one_sided``, ledger 484),
            and the heads row says how many there are; their caveats included
            (FloPy's ``tslen`` and
            ``time_step`` columns counted as budget terms, so a change of ``nstp``
            alone breaches the budget). When either read fails
            with a missing or unreadable output (``OSError``, ``ValueError``,
            ``KeyError``, ...), the report says "outputs unavailable" for that
            model and the results do not count against "identical".

        Returns
        -------
        str
            Markdown: a heading per compared model, "identical to reference" or
            "differs from reference" with the differing tables.
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
                f"layer {row['argmax_layer']} kstpkper {row['argmax_kstpkper']}"
                + (
                    f"; {int(row['n_one_sided'])} cell(s) with no value in one model only"
                    if int(row.get("n_one_sided", 0) or 0)
                    else ""
                )
                + " |"
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
