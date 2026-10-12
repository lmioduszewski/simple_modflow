"""Grouped package-input exploration (`GroupPackageInputs` + field views)."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

import pandas as pd

from myflopy._logging import get_logger
from myflopy.modflow.mf6.package_explorer import (
    LeafFieldSugar,
    build_cell_input_map_payload,
    build_cell_package_input_table,
    build_group_input_compare_map_payload,
    get_default_group_compare_colorscale,
    get_default_package_colorscale,
    get_default_package_value_column,
    get_package_input_field_names,
)
from myflopy.modflow.mf6.package_explorer_utils import period_number
from myflopy.modflow.mf6.package_registry import input_field_agg
from myflopy.modflow.mf6.package_tables import PACKAGE_TABLE_UNAVAILABLE
from myflopy.modflow.utils.datatypes.hover import cell_input_hover, compare_hover
from myflopy.project.group._shared import (
    _ensure_group_map_compatible,
    _filter_group_input_table,
    _resolve_group_compare_target,
)
from myflopy.project.group.spatial import _GroupSpatialView

if TYPE_CHECKING:
    from myflopy.project.group.core import ModelGroup

logger = get_logger(__name__)


class GroupPackageInputField(_GroupSpatialView):
    """One input field of a grouped package, pinned for the unified grammar.

    ``group.packages.ghb.inputs.cond`` -- same verbs as the parent accessor
    but every panel/series draws this field.
    """

    def __init__(self, parent: GroupPackageInputs, field_name: str):
        """Pin a grouped package-inputs accessor to one field for the unified grammar."""

        self._parent = parent
        self.group = parent.group
        self.package_name = parent.package_name
        self.field_name = str(field_name).lower()

    def get(self, **kwargs) -> pd.DataFrame:
        """Return the aligned input rows narrowed to this field.

        The parent's :meth:`GroupPackageInputs.get` rows, keeping only ``model``,
        ``package``, ``pname``, ``per``, ``layer``, ``cell`` and this field's
        column (those of them the table has).

        Parameters
        ----------
        **kwargs
            Forwarded unchanged to :meth:`GroupPackageInputs.get`, which takes
            ``model_name``, ``per``, ``layer`` and ``cells`` (see there); any other
            name raises ``TypeError``.
        """

        frame = self._parent.get(**kwargs)
        keep = [
            column
            for column in ("model", "package", "pname", "per", "layer", "cell", self.field_name)
            if column in getattr(frame, "columns", [])
        ]
        return frame[keep] if keep else frame

    def summary(self, **kwargs) -> pd.DataFrame:
        """Return the parent's per-model coverage summary.

        The package's summary, not this field's: :meth:`GroupPackageInputs.summary`
        counts the cells and periods holding a record, whichever fields it sets.

        Parameters
        ----------
        **kwargs
            Forwarded unchanged to :meth:`GroupPackageInputs.summary`, which takes
            ``model_name``, ``per`` and ``layer`` (see there); any other name
            raises ``TypeError``.
        """

        return self._parent.summary(**kwargs)

    def _spatial_map(self, **kwargs):
        """Draw one model's map for this pinned field (delegates to the parent accessor)."""

        kwargs.setdefault("value_column", self.field_name)
        return self._parent._spatial_map(**kwargs)

    def _spatial_value_label(self) -> str:
        """The mapped quantity's label -- this pinned field name."""

        return self.field_name

    def _series_value_column(self, frame) -> str:
        """The column ``plot()`` draws -- this pinned field name."""

        return self.field_name

    def _series_default_agg(self) -> str:
        """How ``plot()`` combines a period's records: by the field's kind (ledger 472)."""

        return input_field_agg(self.package_name, self.field_name)


class GroupPackageInputs(LeafFieldSugar, _GroupSpatialView):
    """Grouped input accessor for simple cell-based MF6 stress-period packages.

    Inherits the unified grammar (``map``/``plot``/``mosaic``/``animate``;
    ``mosaic()`` defaults to one panel per model). Registry-backed fields are
    first-class nodes (``inputs.cond``) and every verb takes ``field=`` as
    sugar over them; without ``field=`` the package default field is drawn.
    """

    def __init__(self, group: ModelGroup, package_name: str):
        """Bind a grouped BC-package inputs accessor to ``group`` for one package."""

        self.group = group
        self.package_name = str(package_name).lower()

    def _field_names(self) -> list[str]:
        """The registry-declared input field names for this package."""

        return get_package_input_field_names(self.package_name)

    def _field_node(self, name: str) -> GroupPackageInputField:
        """Return this accessor pinned to field ``name``."""

        return GroupPackageInputField(self, name)

    def get(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """Return aligned stress-period input data for all models.

        Each member's ``build_cell_package_input_table`` for this package type --
        every package of the type, array form included, each row tagged with its
        package's ``pname`` (ledger 229) -- plus a ``model`` column. A member that
        does not carry the package is left out (logged at DEBUG), not an error.

        Parameters
        ----------
        model_name : str or None, optional
            Keep one member's rows, by its name in the group. Default ``None``:
            every member. A name that is not in the group matches nothing and
            returns an empty frame.
        per : int or None, optional
            One stress period, zero-based: each member's records MF6 RUNS in it --
            the last block at or before ``per``, none under an empty block -- as
            the single-model ``model.<pkg>.inputs.get(per=)`` reads them (ledger
            449). A non-integer raises ``TypeError`` (refused before a
            member is read, so it is never mistaken for a missing package); a
            period past a member's ``nper`` raises ``IndexError``. Default
            ``None``: every block as written, each row's ``per`` its block.
        layer : int or list of int or None, optional
            Keep the rows in these zero-based layers (an int, list or tuple; a
            NumPy array raises ``TypeError``). Default ``None``: every layer.
        cells : int or list of int or None, optional
            Keep the rows on these zero-based cell numbers within a layer (the
            ``cell`` column, not the model-wide node). Default ``None``: every
            cell.

        Returns
        -------
        pandas.DataFrame
            The package's own record fields (``cond``, ``bhead``, ... in model
            units) with ``package``, ``pname``, ``per``, ``layer``, ``cell`` and
            ``model``. Empty when no member carries the package.
        """

        # Refused HERE, before the per-member handler below reads a TypeError as
        # "this member carries no package" (ledger 449).
        per = None if per is None else period_number(per)
        frames = []
        for current_model_name, model in self.group.models.items():
            # A group member that does not carry this package is skipped, NOT an
            # error. It used to raise, and because `compare()` (and so
            # `ModelDiff._value_cells_changed`) builds the whole group before
            # filtering to one model, ONE member without the package discarded
            # the comparison for every other member -- reporting models with
            # genuinely different values as identical to the reference.
            try:
                frame = build_cell_package_input_table(model, self.package_name, per=per)
            except PACKAGE_TABLE_UNAVAILABLE:
                logger.debug(
                    "%s carries no readable %s; leaving it out of the group table",
                    current_model_name, self.package_name, exc_info=True,
                )
                continue
            frame["model"] = current_model_name
            frames.append(frame)

        if not frames:
            return pd.DataFrame()
        combined = pd.concat(frames, ignore_index=True)
        return _filter_group_input_table(
            combined,
            model_name=model_name,
            per=per,
            layer=layer,
            cells=cells,
        )

    def summary(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """Return a per-model summary of this package's input coverage.

        One row per member holding a record: how many distinct cells and stress
        periods its rows cover.

        Parameters
        ----------
        model_name : str or None, optional
            Summarize one member, by its name in the group. Default ``None``:
            every member. A name that is not in the group returns an empty frame.
        per : int or None, optional
            Count only the records in force in one stress period, zero-based (see
            :meth:`get`). Default ``None``: every block as written.
        layer : int or list of int or None, optional
            Count only the rows in these zero-based layers. Default ``None``:
            every layer.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``package``, ``cells`` (distinct values of the ``cell``
            column -- a cell number used in two layers counts once) and
            ``periods`` (distinct written periods). Empty, with those columns,
            when :meth:`get` is.
        """

        data = self.get(model_name=model_name, per=per, layer=layer)
        columns = ["model", "package", "cells", "periods"]
        if data.empty:
            return pd.DataFrame(columns=columns)
        rows = []
        for current_model_name, sub in data.groupby("model"):
            rows.append(
                {
                    "model": current_model_name,
                    "package": self.package_name,
                    "cells": int(sub["cell"].nunique()) if "cell" in sub.columns else 0,
                    "periods": int(sub["per"].nunique()) if "per" in sub.columns else 0,
                }
            )
        return pd.DataFrame(rows, columns=columns)

    def compare(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """Compare input values for all models against the reference model.

        Every NUMERIC column of :meth:`get` but ``per``/``layer``/``cell`` is
        differenced; each non-reference member's rows are paired with the
        reference's on every remaining column -- ``package``, ``pname``, ``per``,
        ``layer``, ``cell`` and any text column such as ``boundname`` -- by an
        inner join, so a row with no identical partner is dropped (a member whose
        package carries another name pairs nothing). A cell holding SEVERAL
        records with the same keys pairs them in file order, the member's first
        with the reference's first and so on -- one row per record, a record with
        no counterpart dropped. It paired every one with every one until ledger
        472: two records each gave four rows, and a change of 1.0 in each
        drew 4.0 on the map.

        Parameters
        ----------
        model_name : str or None, optional
            Keep one compared member's rows, by its name, after the comparison.
            Default ``None``: every non-reference member. The reference's own
            name, or a name not in the group, returns an empty frame.
        per : int or None, optional
            Compare the records in force in one stress period, zero-based (see
            :meth:`get`). Default ``None``: every block as written.
        layer : int or list of int or None, optional
            Compare only the rows in these zero-based layers. Default ``None``:
            every layer.
        cells : int or list of int or None, optional
            Compare only these zero-based cell numbers within a layer. Default
            ``None``: every cell.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``reference_model``, the join keys, then for each value
            column ``<col>``, ``reference_<col>`` and ``<col>_diff`` (member minus
            reference, model units). Empty when :meth:`get` is.
        """

        data = self.get(per=per, layer=layer, cells=cells)
        if data.empty:
            return pd.DataFrame()

        reference = self.group.reference
        value_columns = [
            column
            for column in data.columns
            if column not in {"model", "per", "layer", "cell"}
            and pd.api.types.is_numeric_dtype(data[column])
        ]
        key_columns = [column for column in data.columns if column not in {"model", *value_columns}]
        # Records PAIR, in file order within their keys: an inner join on the keys
        # alone gave a cell of k records k*k rows (ledger 472).
        data = data.assign(
            _record=data.groupby(["model", *key_columns], dropna=False).cumcount()
        )
        ref = data[data["model"] == reference].drop(columns=["model"]).copy()
        rename_map = {column: f"reference_{column}" for column in value_columns}
        ref = ref.rename(columns=rename_map)
        comp = data[data["model"] != reference].merge(
            ref, on=[*key_columns, "_record"], how="inner"
        ).drop(columns=["_record"])
        comp["reference_model"] = reference
        for column in value_columns:
            comp[f"{column}_diff"] = comp[column].astype(float) - comp[f"reference_{column}"].astype(float)
        ordered = ["model", "reference_model", *key_columns]
        for column in value_columns:
            ordered.extend([column, f"reference_{column}", f"{column}_diff"])
        result = comp[ordered]
        return _filter_group_input_table(
            result,
            model_name=model_name,
            per=per,
            layer=layer,
            cells=cells,
        )

    def _spatial_map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        model=None,
        value_column: str | None = None,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str | None = None,
        colorscale: str | None = None,
        **kwargs,
    ):
        """Render one model's package-input field as a raw-value ``Choro``.

        ``agg=None`` combines a cell's records by the drawn field's kind (ledger
        472): it was ``"sum"``, so two ``bhead`` records of 5 and 7 drew 12.
        """

        target_name = self._group_target(model)
        target_model = self.group.models[target_name]
        selected = self.get(model_name=target_name, per=per, layer=layer)
        fallback = get_default_package_value_column(self.package_name)
        source = selected if not selected.empty else self.get(model_name=target_name)
        numeric_columns = [
            column
            for column in source.columns
            if column not in {"model", "package", "per", "layer", "cell"}
            and pd.api.types.is_numeric_dtype(source[column])
        ]
        chosen_value_column = (
            value_column
            or (fallback if fallback in source.columns else None)
            or (numeric_columns[0] if numeric_columns else None)
        )
        if chosen_value_column is None:
            raise ValueError(f"Could not infer a numeric value column for package {self.package_name!r}.")
        values, hover = build_cell_input_map_payload(
            selected,
            ncpl=target_model.vor.ncpl,
            value_column=chosen_value_column,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=input_field_agg(self.package_name, chosen_value_column, agg),
        )
        kwargs.setdefault("hover_spec", cell_input_hover(chosen_value_column))
        return target_model.plot.map(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or get_default_package_colorscale(self.package_name),
            **kwargs,
        )

    def compare_map(
        self,
        *,
        model_name: str | None = None,
        per: int = 0,
        layer: int = 0,
        value_column: str | None = None,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str | None = None,
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build a diff choropleth against the group's reference model.

        Draws ``<value_column>_diff`` from :meth:`compare` (member minus
        reference) on the reference model's grid, one value per cell.

        Parameters
        ----------
        model_name : str or None, optional
            The non-reference member to compare. Default ``None`` takes the
            group's one non-reference member and raises ``ValueError`` when there
            are several; the reference itself raises ``ValueError`` and a name not
            in the group ``KeyError``. The member must have the reference's
            discretization FILE content, read from each workspace (a changed
            ``botm`` or ``idomain`` refuses too), grid type and ``ncpl``, or
            ``ValueError``.
        per : int, default 0
            The stress period whose records in force are compared, zero-based
            (see :meth:`get`). A period with no rows draws every cell at
            ``fill_value``.
        layer : int, default 0
            The zero-based layer whose rows are drawn.
        value_column : str or None, optional
            The record field whose difference is drawn (``"cond"``, ``"bhead"``,
            ...). Default ``None``: the package's registry default field when the
            comparison has it, else its first numeric column; ``ValueError`` when
            there is none, and ``KeyError`` when the field has no ``_diff``
            column.
        multiplier : float, default 1.0
            Scales every difference BEFORE ``agg`` -- a unit conversion or a
            sign flip -- and never the ``fill_value`` cells (it acted after
            ``agg`` and scaled them too until ledger 473).
        fill_value : float, default 0.0
            The value of a cell with no comparison row. ``float("nan")`` leaves
            it blank. A cell whose rows hold no difference at all is ``NaN``, not
            this (ledger 329).
        agg : str or None, optional
            How a cell's paired differences combine -- any pandas aggregation
            name. Default ``None``: by the field's kind
            (``package_registry.input_field_agg``, ledger 297) -- the SUM of a
            conductance's or a rate's differences (the change in the cell's
            total), the MEAN of an elevation's. A field the registry does not
            declare raises ``ValueError`` asking for ``agg``. It was ``"sum"``
            for every field, over cross-paired rows, until ledger 472: a
            cell of two ``bhead`` records each moved 1.0 drew 4.0.
        colorscale : str or list or None, optional
            Any Plotly colorscale name or explicit stops. Default ``None``: the
            group diff scale -- red at the negative end, white at zero, blue at
            the positive end.
        **kwargs
            Forwarded to the reference model's ``model.plot.map(...)``: its named
            drawing options (``contours``, ``select``, ``features``, ``extent``,
            ``backend``, ...) and its Plotly trace tail. ``zmin``/``zmax``
            default to plus/minus the largest absolute value drawn (a filled
            cell's included) and ``hover_spec`` to the comparison hover; passing
            them overrides those. The names this method fills itself beyond its
            own parameters (``type``, ``custom_zs``, ``custom_hover``,
            ``hover_heads``, ``hover_ks``) raise ``TypeError``; ``values=``
            REPLACES the drawn differences while the hover keeps them.

        Returns
        -------
        Choro or matplotlib.figure.Figure
            The map picture, or a bare Matplotlib figure with ``backend="mpl"``.
        """

        target_name = _resolve_group_compare_target(self.group, model_name)
        _ensure_group_map_compatible(self.group, target_name)
        reference_model = self.group.models[self.group.reference]
        selected = self.compare(model_name=target_name, per=per, layer=layer)
        fallback = get_default_package_value_column(self.package_name)
        source = selected if not selected.empty else self.compare(model_name=target_name)
        numeric_columns = [
            column
            for column in source.columns
            if column not in {"per", "layer", "cell", "model", "reference_model", "package"}
            and pd.api.types.is_numeric_dtype(source[column])
            and not column.startswith("reference_")
            and not column.endswith("_diff")
        ]
        chosen_value_column = (
            value_column
            or (fallback if fallback in source.columns else None)
            or (numeric_columns[0] if numeric_columns else None)
        )
        if chosen_value_column is None:
            raise ValueError(f"Could not infer a numeric value column for package {self.package_name!r}.")
        diff_column = f"{chosen_value_column}_diff"
        if diff_column not in source.columns:
            raise KeyError(f"Difference column {diff_column!r} was not found.")
        values, hover, absmax = build_group_input_compare_map_payload(
            selected,
            ncpl=reference_model.vor.ncpl,
            value_column=chosen_value_column,
            diff_column=diff_column,
            model_name=target_name,
            reference_model=self.group.reference,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=input_field_agg(self.package_name, chosen_value_column, agg),
        )
        kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
        kwargs.setdefault("zmax", absmax if absmax > 0 else None)
        kwargs.setdefault("hover_spec", compare_hover(chosen_value_column, diff_column))
        return reference_model.plot.map(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or get_default_group_compare_colorscale(),
            **kwargs,
        )


