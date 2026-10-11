"""Grouped per-cell package results (base class + namespace)."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

import pandas as pd

from myflopy.modflow.mf6.package_explorer import (
    FieldMappable,
    _symmetric_color_limit,
    budget_value_units,
    build_budget_result_table,
    build_cell_input_map_payload,
    build_group_input_compare_map_payload,
    get_default_group_compare_colorscale,
    get_default_package_colorscale,
)
from myflopy.modflow.utils.datatypes.hover import compare_hover, result_hover
from myflopy.project.group._shared import (
    _ensure_group_map_compatible,
    _filter_group_input_table,
    _reduce_to_period_end,
    _resolve_group_compare_target,
    _stable_compare_keys,
)
from myflopy.project.group.spatial import _GroupSpatialView

if TYPE_CHECKING:
    from myflopy.project.group.core import ModelGroup


class GroupCellPackageResults(_GroupSpatialView):
    """Grouped accessor for cell-based package result tables.

    Inherits the unified ``map`` / ``mosaic`` / ``animate`` grammar (``mosaic()``
    defaults to one panel per model).
    """

    def __init__(self, group: ModelGroup, package_name: str, *, budget_text: str, value_name: str):
        """Bind a grouped cell-budget result accessor: one ``budget_text`` term + value column."""

        self.group = group
        self.package_name = str(package_name).lower()
        self.budget_text = str(budget_text)
        self.value_name = str(value_name)

    def get(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """Return aligned package-result rows for all models in the group.

        Each member's ``build_budget_result_table`` for this term: one row per
        record entry of the cell-budget file, EVERY package of the type read and
        tagged with its ``pname`` (ledger 352), every saved step kept, plus a
        ``model`` column.

        Parameters
        ----------
        model_name : str or None, optional
            Keep one member's rows, by its name in the group. Default ``None``:
            every member. A name that is not in the group matches nothing and
            returns an empty frame.
        per : int or None, optional
            Keep one stress period, zero-based: a FILTER, so every step the
            period saved is kept (:meth:`compare` and the maps read its last).
            Default ``None``: every period. A period nothing saved returns an
            empty frame.
        layer : int or list of int or None, optional
            Keep the rows in these zero-based layers (an int, list or tuple; a
            NumPy array raises ``TypeError``). Default ``None``: every layer.
        cells : int or list of int or None, optional
            Keep the rows on these zero-based cell numbers within a layer (the
            ``cell`` column, not the model-wide ``node``). Default ``None``: every
            cell.

        Returns
        -------
        pandas.DataFrame
            ``node``/``node2`` (zero-based), the term's value column (named by
            this accessor's ``value_name`` -- a flow such as ``q_gwf`` in model
            length^3/time with MF6's sign), ``layer``, ``cell``, ``package``,
            ``pname``, ``kstpkper`` (zero-based, on that member's budget-file
            clock), ``per`` and ``model``. The UZF ``sat`` accessor reads NPF's
            ``DATA-SAT`` record, whose table carries TWO ``sat`` columns -- the
            record's flow (0.0) renamed, then the dimensionless saturation -- so
            its :meth:`compare`, ``map`` and :meth:`compare_map` raise (ledger
            439).
        """

        frames = []
        for current_model_name, model in self.group.models.items():
            frame = build_budget_result_table(
                model,
                budget_text=self.budget_text,
                package_name=self.package_name,
                value_name=self.value_name,
            )
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

    def compare(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """Compare grouped package results against the reference model.

        Each member is first reduced to its own LAST saved step of every period
        (``_reduce_to_period_end``, ledger 337), so members of different ``nstp``
        align by stress period. Rows are then paired with the reference's on the
        integer and text identity columns -- ``node``, ``node2``, ``layer``,
        ``cell``, ``package``, ``pname``, ``per`` -- by an inner join, so a row
        with no partner (a cell only one member has a record on) is dropped;
        ``kstpkper`` and every float column but the value are left out of the
        key AND of the result.

        Parameters
        ----------
        model_name : str or None, optional
            Keep one compared member's rows, by its name, after the comparison.
            Default ``None``: every non-reference member. The reference's own
            name, or a name not in the group, returns an empty frame.
        per : int or None, optional
            Compare one stress period, zero-based, each member at its own last
            saved step of it. Default ``None``: every period.
        layer : int or list of int or None, optional
            Compare only the rows in these zero-based layers. Default ``None``:
            every layer.
        cells : int or list of int or None, optional
            Compare only these zero-based cell numbers within a layer. Default
            ``None``: every cell.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``reference_model``, the join keys, the value column,
            ``reference_<value>`` and ``<value>_diff`` (member minus reference, in
            the value's units). Empty when :meth:`get` is.
        """

        data = self.get(per=per, layer=layer, cells=cells)
        if data.empty:
            return pd.DataFrame()
        data = _reduce_to_period_end(data)

        reference = self.group.reference
        value_columns = [self.value_name]
        key_columns = _stable_compare_keys(data, value_columns)
        ref = data[data["model"] == reference][key_columns + value_columns].copy()
        rename_map = {column: f"reference_{column}" for column in value_columns}
        ref = ref.rename(columns=rename_map)
        comp = data[data["model"] != reference].merge(ref, on=key_columns, how="inner")
        comp["reference_model"] = reference
        comp[f"{self.value_name}_diff"] = (
            comp[self.value_name].astype(float) - comp[f"reference_{self.value_name}"].astype(float)
        )
        ordered = [
            "model",
            "reference_model",
            *key_columns,
            self.value_name,
            f"reference_{self.value_name}",
            f"{self.value_name}_diff",
        ]
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
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Render one model's package result field as a ``Choro``."""

        target_name = self._group_target(model)
        target_model = self.group.models[target_name]
        selected = self.get(model_name=target_name, per=per, layer=layer)
        values, hover = build_cell_input_map_payload(
            selected,
            ncpl=target_model.vor.ncpl,
            value_column=self.value_name,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        # The signed exchange columns name their reference frame (q_gwf / q_lake);
        # "q" is only the legacy raw name. All three take the diverging RdBu scale
        # with a symmetric, zero-centred range.
        is_exchange = self.value_name in {"q", "q_gwf", "q_lake"}
        if is_exchange:
            absmax = _symmetric_color_limit(values)
            kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
            kwargs.setdefault("zmax", absmax if absmax > 0 else None)
            kwargs.setdefault("zmid", 0.0)
        kwargs.setdefault(
            "hover_spec",
            result_hover(
                self.value_name,
                title=f"{self.package_name.upper()} {self.value_name}",
                units=(
                    {self.value_name: budget_value_units(target_model)}
                    if is_exchange
                    else None
                ),
            ),
        )
        return target_model.plot.map(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=(
                colorscale
                or ("RdBu" if is_exchange else None)
                or get_default_package_colorscale(self.package_name)
                or "earth"
            ),
            **kwargs,
        )

    def compare_map(
        self,
        *,
        model_name: str | None = None,
        per: int = 0,
        layer: int = 0,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build a diff choropleth for one grouped package result field.

        Draws ``<value>_diff`` from :meth:`compare` (member minus reference, each
        at its own last saved step of ``per``) on the reference model's grid. The
        public route is the one ``diff()`` verb --
        ``group.diff().packages.<pkg>.results.map(...)``.

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
            The stress period drawn, zero-based -- each member's own last saved
            step of it. A period with no rows draws every cell at ``fill_value``.
        layer : int, default 0
            The zero-based layer whose rows are drawn.
        multiplier : float, default 1.0
            Scales every drawn difference AFTER ``agg`` -- a unit conversion or a
            sign flip -- and the filled cells with them, except when the
            selection is empty (every cell is then ``fill_value`` as given).
        fill_value : float, default 0.0
            The value of a cell with no comparison row. ``float("nan")`` leaves
            it blank. A cell whose rows hold no difference at all is ``NaN``, not
            this (ledger 329).
        agg : str, default "sum"
            How several rows on one cell combine -- two drains in one cell, or
            two packages of the type -- as any pandas aggregation name. ``"sum"``
            is right for a flow, which MF6 superposes.
        colorscale : str or list or None, optional
            Any Plotly colorscale name or explicit stops. Default ``None``: the
            group diff scale -- red at the negative end, white at zero, blue at
            the positive end.
        **kwargs
            Forwarded to the reference model's ``model.plot.map(...)``: its named
            drawing options (``contours``, ``select``, ``features``, ``extent``,
            ``backend``, ...) and its Plotly trace tail. ``zmin``/``zmax`` default
            to plus/minus the largest absolute value drawn (a filled cell's
            included), ``zmid`` to 0 and ``hover_spec`` to the comparison hover;
            passing them overrides those. The names this method fills itself
            beyond its own parameters (``type``, ``custom_zs``, ``custom_hover``,
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
        diff_column = f"{self.value_name}_diff"
        values, hover, absmax = build_group_input_compare_map_payload(
            selected,
            ncpl=reference_model.vor.ncpl,
            value_column=self.value_name,
            diff_column=diff_column,
            model_name=target_name,
            reference_model=self.group.reference,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
        kwargs.setdefault("zmax", absmax if absmax > 0 else None)
        kwargs.setdefault("zmid", 0.0)
        kwargs.setdefault(
            "hover_spec",
            compare_hover(
                self.value_name,
                diff_column,
                title=f"Δ {self.package_name.upper()} {self.value_name} vs reference",
            ),
        )
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

    # NOTE: the series view is the unified grammar's ``plot()`` (SpatialView) --
    # one line per model (and per cell with ``cells=[...]``).


class GroupCellPackageResultsNamespace(FieldMappable):
    """Namespace for grouped cell-based package result accessors.

    ``results.map()`` maps the exchange field ``q`` for one model;
    ``results.mosaic(by="model")`` and ``results.animate(...)`` inherit the
    unified grammar. (SFR/LAK add a ``stage`` field via their subclasses.)
    """

    _default_field = "q"

    def _field_names(self):
        """The single mappable grouped result field: exchange ``q``."""

        return ["q"]

    def __init__(self, result_accessor: GroupCellPackageResults):
        """Wrap a grouped cell-budget result accessor as a ``field=``-aware namespace."""

        self._result_accessor = result_accessor

    @property
    def q(self) -> GroupCellPackageResults:
        """Return the primary grouped package-exchange result accessor."""

        return self._result_accessor


