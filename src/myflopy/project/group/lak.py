"""Grouped LAK budget/stage/connection results + outputs + accessor."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from myflopy.modflow.mf6.hydraulics import unit_labels
from myflopy.modflow.mf6.package_explorer import (
    LakStageResultsExplorer,
    _exchange_colorscale,
    _normalize_connection_type_filter,
    _symmetric_color_limit,
    build_cell_input_map_payload,
    build_group_input_compare_map_payload,
    build_lak_budget_result_table,
    build_lak_connection_table,
    build_lak_q_map_payload,
    build_lak_stage_result_table,
    get_default_group_compare_colorscale,
)
from myflopy.modflow.mf6.package_plotting import (
    _apply_backend,
    refuse_noun_parameters,
    resolve_noun_hover,
)
from myflopy.modflow.mf6.package_surface_water import join_lak_stage
from myflopy.modflow.utils.datatypes.hover import (
    cell_input_hover,
    compare_hover,
    dimension_units,
    lak_hover,
)
from myflopy.project.group._shared import (
    _ensure_group_map_compatible,
    _filter_group_input_table,
    _reduce_to_period_end,
    _resolve_group_compare_target,
)
from myflopy.project.group.results import GroupCellPackageResults, GroupCellPackageResultsNamespace
from myflopy.project.group.spatial import _GroupSpatialView

if TYPE_CHECKING:
    from myflopy.project.group.core import ModelGroup


def _one_stage_per_lake(table: pd.DataFrame) -> pd.DataFrame:
    """A LAK stage table reduced to one row per saved step per lake.

    ``build_lak_stage_result_table`` repeats a lake's one stage on every cell it
    connects to; these are copies (the file holds one stage per lake per record), so
    the first is the lake's stage -- a stage is not added (ledger 297).

    Parameters
    ----------
    table : pandas.DataFrame
        ``build_lak_stage_result_table(model)``.

    Returns
    -------
    pandas.DataFrame
        ``kstpkper``, ``per``, ``lake``, ``stage``, in the table's order.
    """

    frame = table.drop_duplicates(subset=["kstpkper", "lake"])
    return frame.loc[:, ["kstpkper", "per", "lake", "stage"]].reset_index(drop=True)


class GroupLakBudgetResults(GroupCellPackageResults):
    """Grouped LAK exchange accessor with area-normalized maps."""

    def __init__(self, group: ModelGroup, *, budget_text: str = "GWF", value_name: str = "q_lake"):
        """Bind a grouped LAK exchange accessor (defaults to the ``GWF`` term's ``q_lake``).

        LAK's exchange is feature-referenced, so the emitted column is
        ``q_lake`` (accessor stays ``results.q``).
        """

        super().__init__(group, "lak", budget_text=budget_text, value_name=value_name)

    def get(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
        connection_type: str | Sequence[str] | None = None,
    ) -> pd.DataFrame:
        """Return aligned LAK exchange rows, including ``q_per_area``.

        Each member's ``build_lak_budget_result_table``: one row per lake
        connection per saved step of the LAK package's own budget file (its
        ``GWF`` records), every saved step kept, joined to the connection
        geometry, plus a ``model`` column. Every member must carry exactly one
        LAK package: a member with none raises ``AttributeError`` and one with
        several ``SeveralPackagesError`` (each numbers its own lakes, ledger 352),
        either of which ends the whole call -- the group view takes no
        ``package=`` (ledger 439).

        Parameters
        ----------
        model_name : str or None, optional
            Keep one member's rows, by its name in the group. Default ``None``:
            every member. A name that is not in the group matches nothing and
            returns an empty frame.
        per : int or None, optional
            Keep one stress period, zero-based: a FILTER, so every step the
            period saved is kept. Default ``None``: every period.
        layer : int or list of int or None, optional
            Keep the connections in these zero-based layers (an int, list or
            tuple; a NumPy array raises ``TypeError``). Default ``None``: every
            layer.
        cells : int or list of int or None, optional
            Keep the connections on these zero-based cell numbers within a layer.
            Default ``None``: every cell.
        connection_type : str or list of str or None, optional
            Keep only connections of these MF6 types (``"VERTICAL"``,
            ``"HORIZONTAL"``, ``"EMBEDDEDH"``, ``"EMBEDDEDV"``), compared upper-cased
            with each connection's ``claktype``. A spelling that is no MF6 type
            matches nothing and returns an empty frame, without an error.
            Default ``None``: every connection.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``package``, ``kstpkper``, ``per``, ``lake`` and ``iconn``
            (zero-based), ``layer``, ``cell``, the connection's ``claktype``,
            ``belev``, ``telev``, ``connlen``, ``connwidth`` and ``bedleak``,
            ``node``/``node2``, ``q_lake`` -- the lake's side of the exchange, so
            POSITIVE where the lake gains from the aquifer (model
            length^3/time) -- ``FLOW-AREA``/``flow_area``, the WETTED area MF6
            wrote for that step (length^2; it moves with the stage and is 0 for a
            dry connection), and ``q_per_area`` = ``q_lake / flow_area``
            (length/time, ``NaN`` where the area is 0).
        """

        frames = []
        for current_model_name, model in self.group.models.items():
            frame = build_lak_budget_result_table(
                model,
                budget_text=self.budget_text,
                value_name=self.value_name,
            )
            frame["model"] = current_model_name
            frames.append(frame)

        if not frames:
            return pd.DataFrame()
        combined = pd.concat(frames, ignore_index=True)
        combined = _filter_group_input_table(
            combined,
            model_name=model_name,
            per=per,
            layer=layer,
            cells=cells,
        )
        connection_types = _normalize_connection_type_filter(connection_type)
        if connection_types is not None and "claktype" in combined.columns:
            combined = combined.loc[combined["claktype"].astype("string").str.upper().isin(connection_types)].copy()
        return combined.reset_index(drop=True)

    def compare(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
        connection_type: str | Sequence[str] | None = None,
    ) -> pd.DataFrame:
        """Compare grouped LAK exchange rows against the reference model.

        Each non-reference member's rows are paired with the reference's on
        EVERY column of :meth:`get` but ``model``, ``q_lake`` and ``q_per_area``
        by an inner join -- ``kstpkper`` and the float columns ``flow_area``/
        ``FLOW-AREA`` and the geometry included -- and NOT reduced to period ends
        first. So ``kstpkper`` pairs the reference's step ``(k, p)`` with the
        other member's step ``(k, p)``: with different ``nstp`` that is not the
        same time in the period (with one step against two, the reference's
        period end against the other's FIRST step). And because MF6's wetted
        ``FLOW-AREA`` moves with the stage, a connection whose wetted area
        differs between the members has no partner and is dropped (ledger
        439).

        Parameters
        ----------
        model_name : str or None, optional
            Keep one compared member's rows, by its name, after the comparison.
            Default ``None``: every non-reference member. The reference's own
            name, or a name not in the group, returns an empty frame.
        per : int or None, optional
            Compare one stress period, zero-based -- every step it saved, paired
            step by step as above. Default ``None``: every period.
        layer : int or list of int or None, optional
            Compare only the connections in these zero-based layers. Default
            ``None``: every layer.
        cells : int or list of int or None, optional
            Compare only these zero-based cell numbers within a layer. Default
            ``None``: every cell.
        connection_type : str or list of str or None, optional
            Keep only connections of these MF6 types (``"VERTICAL"``,
            ``"HORIZONTAL"``, ``"EMBEDDEDH"``, ``"EMBEDDEDV"``), compared upper-cased
            with each connection's ``claktype``. A spelling that is no MF6 type
            matches nothing and compares nothing, without an error.
            Default ``None``: every connection.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``reference_model``, the join keys, ``q_lake``,
            ``reference_q_lake``, ``q_lake_diff``, and the same three for
            ``q_per_area`` -- each difference member minus reference, the lake's
            sign. Empty when :meth:`get` is.
        """

        data = self.get(per=per, layer=layer, cells=cells, connection_type=connection_type)
        if data.empty:
            return pd.DataFrame()

        reference = self.group.reference
        q = self.value_name  # feature-referenced exchange column (q_lake)
        value_columns = [q, "q_per_area"] if "q_per_area" in data.columns else [q]
        key_columns = [column for column in data.columns if column not in {"model", *value_columns}]
        ref = data[data["model"] == reference].drop(columns=["model"]).copy()
        rename_map = {column: f"reference_{column}" for column in value_columns}
        ref = ref.rename(columns=rename_map)
        comp = data[data["model"] != reference].merge(ref, on=key_columns, how="inner")
        comp["reference_model"] = reference
        comp[f"{q}_diff"] = comp[q].astype(float) - comp[f"reference_{q}"].astype(float)
        if "q_per_area" in value_columns:
            comp["q_per_area_diff"] = comp["q_per_area"].astype(float) - comp["reference_q_per_area"].astype(float)
        ordered = ["model", "reference_model", *key_columns, q, f"reference_{q}", f"{q}_diff"]
        if "q_per_area" in value_columns:
            ordered.extend(["q_per_area", "reference_q_per_area", "q_per_area_diff"])
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
        connection_type: str | Sequence[str] | None = None,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Render one model's LAK exchange (area-normalized) as a ``Choro``."""

        del agg
        target_name = self._group_target(model)
        target_model = self.group.models[target_name]
        selected = join_lak_stage(
            target_model,
            self.get(model_name=target_name, per=per, layer=layer, connection_type=connection_type),
            per=per,
        )
        values, hover = build_lak_q_map_payload(
            selected,
            ncpl=target_model.vor.ncpl,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            value_column=self.value_name,
        )
        absmax = _symmetric_color_limit(values)
        kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
        kwargs.setdefault("zmax", absmax if absmax > 0 else None)
        kwargs.setdefault("zmid", 0.0)
        length, time = unit_labels(target_model)
        kwargs.setdefault("hover_spec", lak_hover(length=length, time=time))
        return target_model.plot.map(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            # match the SFR convention: gaining (negative q) blue, losing red
            # LAK exchange is feature-referenced (gaining is POSITIVE), so blue
            # sits at the positive end -- matching the single-model LAK map. Using
            # the raw blue-at-negative scale here inverted the group map (gaining
            # drew red) while the single-model map drew it blue.
            colorscale=colorscale or _exchange_colorscale("feature"),
            **kwargs,
        )

    def compare_map(
        self,
        *,
        model_name: str | None = None,
        per: int = 0,
        layer: int = 0,
        connection_type: str | Sequence[str] | None = None,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build a grouped LAK diff map using normalized exchange per area.

        Each member's connections of period ``per`` (its own LAST saved step,
        ledger 348) are reduced to one value per cell, ``sum(q_lake) /
        sum(flow_area)`` over the connections in the cell -- the flow per unit
        wetted area, positive where the lake gains -- and the map draws member
        minus reference on the reference's grid, for the cells both members
        connect a lake to. (Unlike :meth:`compare`, this does align each member
        at its own period end.)

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
            The zero-based layer whose connections are drawn.
        connection_type : str or list of str or None, optional
            Keep only connections of these MF6 types (``"VERTICAL"``,
            ``"HORIZONTAL"``, ``"EMBEDDEDH"``, ``"EMBEDDEDV"``), compared upper-cased
            with each connection's ``claktype``. A spelling that is no MF6 type
            matches nothing and draws every cell at ``fill_value``, without an error.
            Default ``None``: every connection.
        multiplier : float, default 1.0
            Scales every drawn difference -- a unit conversion or a sign flip --
            and never the ``fill_value`` cells (it scaled them too until ledger
            473).
        fill_value : float, default 0.0
            The value of a cell that is not a lake cell of both members -- one
            alone connects a lake to it, or neither does. ``float("nan")`` leaves
            it blank. A cell whose connections are all dry (wetted area 0) in
            either member is ``NaN``.
        agg : str, default "sum"
            Accepted and IGNORED: the per-cell reduction is always the
            area-weighted ``sum(q_lake) / sum(flow_area)`` above, whatever is
            passed (ledger 439).
        colorscale : str or list or None, optional
            Any Plotly colorscale name or explicit stops. Default ``None``: the
            group diff scale -- red at the negative end, white at zero, blue at
            the positive end (a lake gaining MORE than the reference's is
            positive, so blue).
        **kwargs
            Forwarded to the reference model's ``model.plot.map(...)``: its named
            drawing options (``contours``, ``select``, ``features``, ``extent``,
            ``backend``, ...) and its Plotly trace tail. ``zmin``/``zmax`` default
            to plus/minus the largest absolute value drawn (a filled cell's
            included), ``zmid`` to 0 and ``hover_spec`` to the comparison hover,
            whose units are length/time in the reference model's own units
            (ledger 434; a bare number when the model declares none); passing
            them overrides those. The names this method fills itself beyond its
            own parameters (``type``, ``custom_zs``, ``custom_hover``,
            ``hover_heads``, ``hover_ks``) raise ``TypeError``; ``values=``
            REPLACES the drawn differences while the hover keeps them.

        Returns
        -------
        Choro or matplotlib.figure.Figure
            The map picture, or a bare Matplotlib figure with ``backend="mpl"``.
        """

        del agg
        target_name = _resolve_group_compare_target(self.group, model_name)
        _ensure_group_map_compatible(self.group, target_name)
        reference_model = self.group.models[self.group.reference]
        # Each model's last saved step of the period (ledger 348): the per-cell sums
        # below added every step the period saved.
        data = _reduce_to_period_end(
            self.get(per=per, layer=layer, connection_type=connection_type)
        )

        def _normalized_by_cell(frame: pd.DataFrame, current_model_name: str) -> pd.DataFrame:
            """Per-cell exchange per unit connection area (``sum(q)/sum(flow_area)``) for one model."""

            selected = frame[frame["model"] == current_model_name].copy()
            if selected.empty:
                return pd.DataFrame(columns=["cell", "q_per_area"])
            q = self.value_name
            selected[q] = pd.to_numeric(selected[q], errors="coerce")
            selected["flow_area"] = pd.to_numeric(selected["flow_area"], errors="coerce")
            grouped = selected.groupby("cell", as_index=False).agg({q: "sum", "flow_area": "sum"})
            grouped["q_per_area"] = np.where(grouped["flow_area"] > 0.0, grouped[q] / grouped["flow_area"], np.nan)
            return grouped[["cell", "q_per_area"]]

        reference = _normalized_by_cell(data, self.group.reference).rename(
            columns={"q_per_area": "reference_q_per_area"}
        )
        target = _normalized_by_cell(data, target_name)
        comparison = target.merge(reference, on="cell", how="inner")
        comparison["model"] = target_name
        comparison["reference_model"] = self.group.reference
        comparison["per"] = int(per)
        comparison["layer"] = int(layer)
        comparison["q_per_area_diff"] = (
            comparison["q_per_area"].astype(float) - comparison["reference_q_per_area"].astype(float)
        )
        values, hover, absmax = build_group_input_compare_map_payload(
            comparison,
            ncpl=reference_model.vor.ncpl,
            value_column="q_per_area",
            diff_column="q_per_area_diff",
            model_name=target_name,
            reference_model=self.group.reference,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg="sum",
        )
        kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
        kwargs.setdefault("zmax", absmax if absmax > 0 else None)
        kwargs.setdefault("zmid", 0.0)
        kwargs.setdefault(
            "hover_spec",
            compare_hover(
                "q_per_area",
                "q_per_area_diff",
                title="Δ lake exchange vs reference",
                # q / area is L/T, in the drawn model's own units (ledger 434).
                units=dimension_units(*unit_labels(reference_model), {
                    "q_per_area": (1, -1),
                    "reference_q_per_area": (1, -1),
                    "q_per_area_diff": (1, -1),
                }),
                labels={"q_per_area_diff": "Δ q / area"},
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


class GroupLakOutputs:
    """Lake-output accessor for :class:`ModelGroup`."""

    def __init__(self, group: ModelGroup):
        """Bind the grouped lake-output accessor to ``group``."""

        self.group = group

    def stage(self) -> pd.DataFrame:
        """Return every saved lake stage of every member as one long table.

        The group twin of ``model.outputs.lak.stage``: EVERY record the stage file
        saved, one row per lake, labelled by the ``kstpkper`` in that file's own
        record headers (``build_lak_stage_result_table``, ledger 348). It paired
        the records with ``model.kstpkper`` -- the head file's clock, which MF6
        gates on the same ``SAVE HEAD`` flag, so the labels were right on a file
        MF6 wrote -- and silently cut the records to that list's length; a stage
        file whose record count differs from its headers now raises instead. For
        one value per period use ``group.packages.lak.results.stage.get()``.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``kstpkper`` (zero-based tuple), ``lake`` (zero-based) and
            ``stage`` (model length units, ``NaN`` where a lake has none, ledger
            329), members in group order, each in its file's order.

        Raises
        ------
        ValueError
            When a member's stage file does not line up with its record headers.
        """

        rows: list[pd.DataFrame] = []
        for model_name, model in self.group.models.items():
            frame = _one_stage_per_lake(build_lak_stage_result_table(model))
            frame["model"] = model_name
            rows.append(frame)

        if not rows:
            return pd.DataFrame(columns=["model", "kstpkper", "lake", "stage"])
        return pd.concat(rows, ignore_index=True)[["model", "kstpkper", "lake", "stage"]]


class GroupLakStageResults(_GroupSpatialView):
    """Grouped accessor for lake stages and stage comparisons.

    Inherits the :class:`SpatialView` grammar (``map``/``mosaic``/``animate``)
    with the model axis being the group's members; each panel delegates to the
    single-model LAK stage explorer so grouped stage maps match the single-model
    ``lak.results.stage.map`` exactly.
    """

    def __init__(self, group: ModelGroup):
        """Bind the grouped lake stage-results accessor to ``group``."""

        self.group = group

    def _spatial_value_label(self) -> str:
        """The mapped quantity's label -- lake stage."""

        return "stage"

    def _spatial_map(self, *, per: int = 0, layer: int = 0, model=None, **kwargs):
        """Render one model's lake stage broadcast to connected cells."""

        target_model = self.group.models[self._group_target(model)]
        return LakStageResultsExplorer(target_model).map(
            per=int(per), layer=int(layer), backend="plotly", **kwargs
        )

    def get(
        self,
        *,
        model_name: str | None = None,
        lake: int | None = None,
        per: int | None = None,
    ) -> pd.DataFrame:
        """Return each member's lake stages, one per lake per stress period.

        Each period is its LAST saved step on each member's own stage-file clock
        (``_reduce_to_period_end``, ledger 337's rule for a group table), so two
        members with different ``nstp`` align on the period. Until ledger 348
        ``per`` was the index of the saved RECORD (ledger 301's stage clock): at
        ``nstp=3`` ``per=1`` read step 2 of period 0, and a member saving more
        steps than the reference was compared step-for-record against it.

        Parameters
        ----------
        model_name : str or None, optional
            Keep one member, by its name in the group. Default ``None``: every
            member.
        lake : int or None, optional
            Keep one lake, zero-based. Default ``None``: every lake.
        per : int or None, optional
            Keep one stress period, zero-based -- its last saved step. Default
            ``None``: every period.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``kstpkper`` (the step read, zero-based), ``per``, ``lake``
            and ``stage`` (model length units, ``NaN`` where a lake has none).

        Raises
        ------
        ValueError
            When a member's stage file does not line up with its record headers.
        """

        rows: list[pd.DataFrame] = []
        for current_model_name, model in self.group.models.items():
            frame = _one_stage_per_lake(build_lak_stage_result_table(model))
            frame["model"] = current_model_name
            rows.append(frame)

        if not rows:
            return pd.DataFrame(columns=["model", "kstpkper", "per", "lake", "stage"])
        combined = _reduce_to_period_end(pd.concat(rows, ignore_index=True))
        combined = combined[["model", "kstpkper", "per", "lake", "stage"]]
        if model_name is not None:
            combined = combined.loc[combined["model"] == str(model_name)].copy()
        if lake is not None:
            combined = combined.loc[combined["lake"] == int(lake)].copy()
        if per is not None:
            combined = combined.loc[combined["per"] == int(per)].copy()
        return combined.reset_index(drop=True)

    def compare(
        self,
        *,
        model_name: str | None = None,
        lake: int | None = None,
        per: int | None = None,
    ) -> pd.DataFrame:
        """Compare lake stages against the reference model, period by period.

        Rows are aligned on the stress PERIOD and the lake, each model read at its
        own last saved step of the period (:meth:`get`) -- not on the saved-record
        index they were aligned on until ledger 348.

        Parameters
        ----------
        model_name : str or None, optional
            Keep one compared (non-reference) member. Default ``None``: every one.
        lake : int or None, optional
            Keep one lake, zero-based. Default ``None``: every lake.
        per : int or None, optional
            Keep one stress period, zero-based. Default ``None``: every period.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``reference_model``, ``per``, ``lake``, ``stage``,
            ``reference_stage`` and ``stage_diff`` (model minus reference).
        """

        data = self.get(lake=lake, per=per)
        if data.empty:
            return pd.DataFrame(columns=["model", "reference_model", "per", "lake", "stage", "reference_stage", "stage_diff"])

        reference = self.group.reference
        ref = (
            data.loc[data["model"] == reference, ["per", "lake", "stage"]]
            .rename(columns={"stage": "reference_stage"})
            .copy()
        )
        comp = data.loc[data["model"] != reference].merge(ref, on=["per", "lake"], how="inner")
        comp["reference_model"] = reference
        comp["stage_diff"] = comp["stage"].astype(float) - comp["reference_stage"].astype(float)
        if model_name is not None:
            comp = comp.loc[comp["model"] == str(model_name)].copy()
        return comp[["model", "reference_model", "per", "lake", "stage", "reference_stage", "stage_diff"]]

    # NOTE: the series view is the unified grammar's ``plot()`` (SpatialView) --
    # one line per model and lake.


class GroupLakConnections:
    """Grouped accessor for lake-connection geometry."""

    def __init__(self, group: ModelGroup):
        """Bind the grouped lake-connection-geometry accessor to ``group``."""

        self.group = group

    def get(
        self,
        *,
        model_name: str | None = None,
        lake: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """Return aligned LAK connection rows for all models.

        Each member's ``build_lak_connection_table``: one row per
        ``connectiondata`` record of its LAK package, plus a ``model`` column.
        Every member must carry exactly one LAK package: a member with none
        raises ``AttributeError`` and one with several ``SeveralPackagesError``
        (ledger 352), either of which ends the whole call.

        Parameters
        ----------
        model_name : str or None, optional
            Keep one member's rows, by its name in the group. Default ``None``:
            every member. A name that is not in the group matches nothing and
            returns an empty frame.
        lake : int or None, optional
            Keep one lake's connections, by its zero-based lake number
            (``ifno``). Default ``None``: every lake.
        layer : int or list of int or None, optional
            Keep the connections in these zero-based layers (an int, list or
            tuple; a NumPy array raises ``TypeError``). Default ``None``: every
            layer.
        cells : int or list of int or None, optional
            Keep the connections on these zero-based cell numbers within a layer.
            Default ``None``: every cell.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``package``, ``lake``, ``iconn``, ``layer``, ``cell``,
            ``claktype``, ``belev``/``telev``/``connlen``/``connwidth`` (model
            length) as declared, ``bedleak`` and ``connection_area`` (length^2:
            the cell's plan area for a vertical connection, ``connwidth * (telev -
            belev)`` for a horizontal one -- the full face, not MF6's wetted
            ``FLOW-AREA`` -- and ``NaN`` for an embedded one; see
            ``model.lak.connections.get``).
        """

        frames = []
        for current_model_name, model in self.group.models.items():
            frame = build_lak_connection_table(model)
            frame["model"] = current_model_name
            frames.append(frame)

        if not frames:
            return pd.DataFrame()
        combined = pd.concat(frames, ignore_index=True)
        combined = _filter_group_input_table(
            combined,
            model_name=model_name,
            per=None,
            layer=layer,
            cells=cells,
        )
        if lake is not None and "lake" in combined.columns:
            combined = combined.loc[combined["lake"] == int(lake)].copy()
        return combined.reset_index(drop=True)

    def map(
        self,
        *,
        # -- which model, which connections -----------------------------------
        model_name: str | None = None,
        lake: int | None = None,
        layer: int = 0,
        # -- table to cells ---------------------------------------------------
        value_column: str = "connection_area",
        agg: str = "sum",
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        # -- colour -----------------------------------------------------------
        zmin: float | None = None,
        zmax: float | None = None,
        colorscale: str | list | tuple | None = None,
        logscale: bool = False,
        # -- contours ---------------------------------------------------------
        contours: bool | str = False,
        contour_values=None,
        contour_levels: int | float | list = 10,
        contour_color: str = "black",
        contour_width: float = 1.5,
        contour_name: str | None = None,
        contour_clip: bool = True,
        contour_resolution: int = 150,
        contour_method: str = "linear",
        # -- highlighting -----------------------------------------------------
        select=None,
        select_style: str = "outline",
        select_color: str | None = None,
        # -- overlays and framing ---------------------------------------------
        features=None,
        hillshade_path=None,
        fit_bounds: bool = True,
        bounds_padding: float = 0.05,
        extent=None,
        # -- hover ------------------------------------------------------------
        hover=None,
        show_layer_elevs: bool | None = None,
        hover_surfaces: bool | None = None,
        # -- renderer ---------------------------------------------------------
        backend: str = "plotly",
        **trace_kwargs,
    ):
        """LAK connection geometry for ONE member model, as a choropleth.

        A GROUP noun drawing a single picture. The group's axis is the model
        list, but a choropleth has one grid, so ``model_name`` picks the member
        and everything else describes that member's map. The rows come from
        :meth:`get` -- one per LAK connection -- and are reduced to one value
        per cell by ``agg``, which is why this signature carries
        ``value_column``/``agg``/``multiplier``/``fill_value`` on top of the
        shared drawing parameters. Use ``mosaic`` on a results noun when you
        want every member side by side; connection geometry is static input, so
        the useful comparison is one member at a time against another.

        Every parameter is named rather than swept into ``**kwargs``: PyCharm
        and Pylance read the ``def`` line and never run the module, so a
        parameter that arrives through a tail is one no editor can ever offer
        (plan 8.8). ``backend`` in particular used to reach the renderer only
        through that tail.

        The free verb's time and head arguments -- ``kstpkper``,
        ``per_timestep``, ``bgs``, ``hover_layers``, ``show_mounding`` -- are
        not named here: connection geometry has no time axis, so each was
        measured to leave the figure byte-identical, and ``show_mounding`` is
        worse than inert: it injects a head-derived row into a geometry tooltip
        on every cell. Passed anyway, each raises ``TypeError`` (ledger 465).
        ``show_layer_elevs`` and ``hover_surfaces`` ARE named: the layer frame is
        the grid's, so they act on any map (ledger 306).

        Parameters
        ----------
        model_name : str, optional
            Which member of the group to draw. Defaults to the group's
            reference model (``group.reference``). Pass a member name to draw
            that model's connections instead -- the usual reason is to look at
            the same lake on a refined grid. A name that is not in the group
            raises ``KeyError`` rather than silently falling back.
        lake : int, optional
            Draw only one lake's connections, by zero-based lake number
            (``ifno``). Defaults to ``None``, meaning every lake in the
            package, which on a multi-lake model paints them all into one
            colour range. Set it when one lake's geometry is the question, or
            when a large lake's areas are flattening a small one.
        layer : int, default 0
            Zero-based layer whose connections are drawn. A lake connects
            downward through several layers; only the rows in this layer are
            mapped, so pass the layer the connections you care about actually
            sit in -- ``get()``'s ``layer`` column lists them.
        value_column : str, default 'connection_area'
            Which connection column to colour by. ``connection_area`` is the
            physical interface area (plan-view cell area for a vertical
            connection, ``connwidth * (telev - belev)`` for a horizontal one).
            Any other numeric column of :meth:`get` works -- ``belev``,
            ``telev``, ``connlen``, ``connwidth`` -- and a name that is not in
            the table raises ``KeyError`` naming it.
        agg : {'sum', 'mean', 'min', 'max', 'first', 'last'}, default 'sum'
            How several connections landing in ONE cell are reduced to one
            number. ``sum`` is right for ``connection_area``, which is
            extensive -- a cell touched by two lakes offers both areas. Switch
            to ``mean`` for an intensive column such as ``belev`` or
            ``connlen``, where adding two elevations produces a meaningless
            number.
        multiplier : float, default 1.0
            Scale every drawn value -- unit conversion (``0.3048`` for feet to
            metres on an elevation column), or a sign flip. Applied to each
            record BEFORE ``agg``, so ``agg="max"`` with ``multiplier=-1`` draws
            the largest negated value, and never to the ``fill_value`` cells (it
            acted after ``agg`` and scaled the fill until ledger 473).
        fill_value : float, default 0.0
            Value given to cells with no LAK connection in this layer. The
            default colours them zero, which reads as "no lake here" on an
            area map; pass ``float("nan")`` to leave them uncoloured instead,
            which is what you want on an elevation column where zero is a real
            elevation and would distort the colour range.
        **trace_kwargs
            Every other keyword, forwarded to the member's ``model.plot.map``
            after three checks: ``values``, ``type``, ``custom_hover``,
            ``hover_heads`` and ``hover_ks`` raise ``TypeError`` (this noun
            supplies its own field and hover); ``per`` raises ``TypeError`` naming
            ``group.packages.lak.results`` (connection geometry has no time
            axis); and ``hover_spec`` is read as the old spelling of ``hover``,
            ``hover`` winning when both are given. ``custom_zs``, which this noun
            computes, raises ``TypeError`` as a second value. A free-verb
            parameter this noun does not name -- ``zoom``, ``hover_fields``, ...
            -- reaches the verb by that name, while a per-layer field's ``bgs``,
            ``kstpkper``, ``per_timestep``, ``hover_layers`` and
            ``show_mounding`` raise ``TypeError`` (ledger 465). The rest
            are properties of the Plotly
            ``Choroplethmap`` trace -- ``zmid``, ``colorbar``, ``reversescale``,
            ``showscale`` -- checked only when the figure is built, so a
            misspelled one raises ``ValueError`` at ``.fig`` or ``.show()``, and
            IGNORED under ``backend="mpl"``, which never builds that trace.

            .. generated from myflopy.plot.map -- edit there, then run scripts/derive_docstrings.py
        zmin, zmax : float, optional
            Fixed color-scale limits. Set both to hold the scale steady across
            frames or panels; ``zmin >= zmax`` raises rather than rendering one flat
            color. With neither, the range comes from the data.
        colorscale : str or list of (float, str), optional
            Any Plotly colorscale name -- ``'earth'`` (the default), ``'Viridis'``,
            ``'Turbo'``, ``'RdBu'``, ...; append ``_r`` to reverse one -- or explicit
            ``[[position, color], ...]`` stops. The static backend builds its
            colormap from the same stops, so a name renders the same way on both.
            An unknown name warns and falls back to ``'earth'``.
        logscale : bool, default False
            Color on a log10 scale, for a field that spans orders of magnitude
            (hydraulic conductivity, recharge). The colorbar is labelled in real
            units at whole decades on both backends and the hover keeps the real
            values. ``zmin``/``zmax`` (and a ``zmid``) are then in DATA units --
            ``zmin=1e-4, zmax=10`` -- and a non-positive limit raises. Non-positive
            values draw as gaps.
        contours : bool or str, default False
            Overlay contour lines. ``True`` contours the mapped values; a string
            names a different field to contour instead.
        contour_values : sequence of float, optional
            Contour a supplied array rather than the mapped one.
        contour_levels : int or float or list of float, default 10
            A count of levels, a fixed interval, or explicit level values.
        contour_color : str, default 'black'
            Line colour for the contour trace.
        contour_width : float, default 1.5
            Line width for the contour trace.
        contour_name : str, optional
            Legend name for the contour trace.
        contour_clip : bool, default True
            Clip contours to the active domain instead of the full grid extent.
        contour_resolution : int, default 150
            Grid size used to build the contours, per axis. Higher is smoother and
            slower. **Only acts under ``contour_method="cubic"``**, which is the one
            that interpolates onto a square grid before contouring; the default
            linear method triangulates the cell centres directly, so there is no
            grid for this to size and passing it changes nothing (measured: 338
            contour points at 40, 150 and 400 alike).
        contour_method : {'linear', 'cubic'}, default 'linear'
            How the scattered cell values become contours. ``'linear'``
            triangulates the cell centres and contours the triangulation --
            fast, exact at the centres, and faceted. ``'cubic'`` interpolates onto
            a ``contour_resolution``-square grid first (Clough-Tocher) and contours
            that -- smoother, slower, and able to overshoot between cells.
            ``'tri'``/``'tricontour'`` and ``'clough'``/``'clough_tocher'``/
            ``'cloughtocher'`` are accepted as aliases. Anything else raises naming
            both; ``'nearest'`` in particular was documented here for a while and
            has never been implemented.
        select : sequence of int, ndarray, str, Path or geometry, optional
            Cells to highlight. Cell indices, a boolean mask of length ``ncpl``, the
            name of a registered model region (``"all_streams"``), a path to a
            vector file, or a shapely/GeoPandas geometry, by coverage. Highlights
            nothing (with a warning) when the selection is empty.
        select_style : {'outline', 'dim', 'both'}, default 'outline'
            How the highlight is drawn. ``'outline'`` traces the dissolved boundary
            of the selection and leaves every cell at full opacity. ``'dim'`` fades
            the *unselected* cells to 20% instead, which suits a bare grid where the
            selection is the subject, but on a field it costs a measured 6.9x of
            readable contrast everywhere you did not select -- and one box/lasso
            gesture in the browser overwrites it. ``'both'`` draws each.
        select_color : str, optional
            Highlight colour, defaulting to ``viz.PALETTE.highlight``. Worth setting
            when the default red collides with a red-blue diverging colorscale.
        features : path, GeoDataFrame, geometry, Feature, or list/dict of those, optional
            Custom geometry drawn over the cells -- wells, pond footprints, a
            boundary. A path to any vector file GeoPandas reads (name a GeoPackage
            layer with ``Feature(path, layer=...)``), a ``GeoDataFrame`` or
            ``GeoSeries``, a shapely geometry, or a :class:`myflopy.plot.Feature`
            carrying its own ``color``/``width``/``size``/``fill``/``label``/
            ``hover``. A list draws several layers; a dict names them. Points
            become markers, lines polylines, polygons outlines (or fills). Each
            layer gets a legend entry, a stable house color and a hover listing its
            attribute columns. A source without a CRS is taken to be in the grid's.
            Drawn by both backends and kept by ``mosaic`` and ``animate``.
        hillshade_path : Path, optional
            A hillshade GeoTIFF to draw beneath the cells for topographic context.
        fit_bounds : bool, default True
            Fit the initial view to the grid extent rather than using ``zoom``. The
            fit is for the figure's own size, or plotly's default 700x450 when it
            names none -- a larger display only adds margin around the grid -- and
            ``.save()`` to an image refits for the export's size.
        bounds_padding : float, default 0.05
            Fractional padding around the fitted bounds.
        extent : tuple, points, geometry, GeoDataFrame, area or str, optional
            The window the map OPENS on, on both backends -- the one way to frame a
            static (matplotlib) map, which cannot be zoomed. Four numbers
            ``(xmin, ymin, xmax, ymax)`` in the grid's CRS (shapely's order, NOT
            FloPy's ``(xmin, xmax, ymin, ymax)``); a run of ``(x, y)`` points (their
            bounding box); any shapely geometry -- a polygon, a transect line -- by its
            bounds; a GeoDataFrame, reprojected from its own CRS; or the name of an
            area of interest (``model.aoi``). It sets the VIEW only: no cell is
            cropped and the colour scale still spans the whole field, so pass
            ``zmin``/``zmax`` to tighten it for a close-up. ``bounds_padding`` pads it
            as it pads the grid. A box off the grid raises, printing both boxes and
            both CRSs. In a ``mosaic``, one extent frames every panel and a LIST of
            them frames each panel on its own.
        hover : HoverSpec, optional
            Replace the sectioned hover outright. See
            :mod:`myflopy.modflow.utils.datatypes.hover`.
        hover_surfaces : bool, optional
            Add the model-top / layer-bottom table to the sectioned hover.
        show_layer_elevs : bool, optional
            Add model-top and per-layer-bottom rows to the hover. Defaults to
            whether the grid actually carries layer elevations (``vor.gdf_topbtm``),
            because forcing it on a grid without them raises.
        backend : {'plotly', 'mpl'}, default 'plotly'
            Which renderer draws the map. ``'plotly'`` returns the interactive
            ``Choro`` picture -- pan, zoom, hover, a basemap. ``'mpl'`` returns a
            static :class:`matplotlib.figure.Figure` instead, for a report, a
            multi-panel figure of your own, or anywhere a live figure is not wanted.
            Accepts ``'interactive'`` and ``'matplotlib'``/``'static'`` as aliases;
            anything else raises rather than being ignored.

            The switch changes the RENDERER, never the subject: both backends draw
            this same map, and a ``colorscale`` -- named or as stops -- renders the
            same ramp on each (the static colormap is built from Plotly's own
            stops). One difference is worth knowing before you rely on it: the
            Matplotlib branch draws in **model coordinates** with no basemap, so
            ``bgs`` and ``zoom`` have nothing to act on there, while ``features``,
            the ``select`` outline and a log colorbar all carry across. Pass explicit
            stops for a diverging scale when the direction carries meaning: the
            stops state which end is red, which a bare name cannot (ledger 253).

            ``backend='mpl'`` and ``.plot_mpl()`` on the returned picture are the
            same renderer reached two ways. Prefer the parameter: it is the spelling
            the whole grammar shares, so it also works on the nouns
            (``model.hds.map(backend='mpl')``) and on the composers
            (``mosaic``/``animate``), where there is no intermediate picture to call
            a method on.

            .. end generated

        Returns
        -------
        Choro or matplotlib.figure.Figure
            A :class:`~myflopy.viz.Picture` under the default backend --
            ``.show()``, ``.save(path)``, ``.html(path)`` -- or a bare
            Matplotlib figure with ``backend="mpl"``.

        See Also
        --------
        get : the connection rows behind the picture, as a DataFrame, aligned
            across every model in the group.

        Examples
        --------
        >>> group.packages.lak.connections.map()
        >>> group.packages.lak.connections.map(model_name="refined")
        >>> group.packages.lak.connections.map(lake=0, layer=1)
        >>> group.packages.lak.connections.map(
        ...     value_column="belev", agg="mean", fill_value=float("nan")
        ... )
        >>> group.packages.lak.connections.map(logscale=True, colorscale="Viridis")
        >>> group.packages.lak.connections.map(backend="mpl").savefig("lak_conn.png")
        """

        refuse_noun_parameters(
            "group.packages.lak.connections",
            value_column,
            trace_kwargs,
        )
        hover = resolve_noun_hover(hover, trace_kwargs)
        if "per" in trace_kwargs:
            # Every SIBLING noun takes `per`, so reaching for it here is the
            # natural mistake -- but connection geometry is static and the
            # forward hardcodes `per=0`, so the name collides in the tail and
            # dies as "ModelPlots.map() got multiple values for keyword
            # argument 'per'", naming a class the caller never typed.
            raise TypeError(
                "group.packages.lak.connections.map() has no per=: LAK "
                "connection geometry does not change with time, so there is "
                "one picture for every stress period. Use "
                "`group.packages.lak.results.<noun>` for a quantity that does."
            )
        target_name = self.group.reference if model_name is None else str(model_name)
        if target_name not in self.group.models:
            raise KeyError(f"Model {target_name!r} is not in the group.")
        target_model = self.group.models[target_name]
        selected = self.get(model_name=target_name, lake=lake, layer=layer)
        if value_column not in selected.columns:
            raise KeyError(f"LAK connection column {value_column!r} was not found.")
        values, hover_table = build_cell_input_map_payload(
            selected,
            ncpl=target_model.vor.ncpl,
            value_column=value_column,
            per=None,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        choro = target_model.plot.map(
            per=0,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover_table,
            hover_heads=False,
            hover_ks=False,
            # The noun's OWN spec is the picture's BASE hover, so it goes in
            # `hover_spec=`; `hover=` is the call-site override slot and stays the
            # caller's. Putting the default in the override slot renders the same
            # but leaves `choro.hover_spec` None, which disagrees with `model.hds`
            # and with what this noun did before the 8.8 conversion.
            hover=hover,
            show_layer_elevs=show_layer_elevs,
            hover_surfaces=hover_surfaces,
            hover_spec=cell_input_hover(value_column),
            zmin=zmin,
            zmax=zmax,
            colorscale=colorscale or "earth",
            logscale=logscale,
            contours=contours,
            contour_values=contour_values,
            contour_levels=contour_levels,
            contour_color=contour_color,
            contour_width=contour_width,
            contour_name=contour_name,
            contour_clip=contour_clip,
            contour_resolution=contour_resolution,
            contour_method=contour_method,
            select=select,
            select_style=select_style,
            select_color=select_color,
            features=features,
            hillshade_path=hillshade_path,
            fit_bounds=fit_bounds,
            bounds_padding=bounds_padding,
            extent=extent,
            **trace_kwargs,
        )
        return _apply_backend(choro, backend)


class GroupLakResultsNamespace(GroupCellPackageResultsNamespace):
    """Namespace for grouped LAK result accessors."""

    def _field_names(self):
        """The mappable grouped LAK result fields: exchange ``q`` and ``stage``."""

        return ["q", "stage"]

    @property
    def stage(self) -> GroupLakStageResults:
        """Return grouped LAK stage helpers."""

        return GroupLakStageResults(self._result_accessor.group)

    @property
    def q(self) -> GroupLakBudgetResults:
        """Return grouped LAK exchange helpers with area-normalized map behavior."""

        return self._result_accessor


class GroupLakPackageAccessor:
    """Namespace for grouped LAK geometry and result helpers."""

    def __init__(self, group: ModelGroup, results_namespace: GroupLakResultsNamespace):
        """Bind the grouped LAK package accessor (connections + results) to ``group``."""

        self.group = group
        self._results_namespace = results_namespace

    @property
    def connections(self) -> GroupLakConnections:
        """Return grouped LAK connection-geometry helpers."""

        return GroupLakConnections(self.group)

    @property
    def results(self) -> GroupLakResultsNamespace:
        """Return grouped LAK result helpers."""

        return self._results_namespace




from myflopy.modflow.mf6.package_plotting import NOUN_MAP_PARAMS  # noqa: E402
from myflopy.plot import inherit_map_docs  # noqa: E402


# --- the shared map reference, spliced onto this module's noun ----------------
#
# `myflopy.plot` (layer 7) splices every noun BELOW it at import. This module is
# layer 9, so reaching down for the reference is this module's job: doing it from
# `plot` pulls in `myflopy.project` (layer 13) while `simulation.base` is still
# half-built, and the import fails outright. So the reference is fetched HERE
# instead, at the bottom of the module -- module-level (this file is layer 9 and
# both targets are below it, so the ratchet stays where it was) but placed after
# the class, which must exist before it can be spliced.
def _inherit_group_lak_docs() -> None:
    """Give ``GroupLakConnections.map`` the free verb's entries for what it takes."""

    # `keep` is exactly the shared drawing set. Deliberately NOT `| {"per",
    # "layer"}`:
    #   `per`   -- this noun has no `per` parameter (connection geometry is
    #              static, so the forward hardcodes `per=0`). Inheriting its
    #              entry documented a knob that does not exist, and
    #              `map(per=3)` reached the tail and died as
    #              `ModelPlots.map() got multiple values for keyword argument
    #              'per'` -- a TypeError naming a class the caller never typed,
    #              the exact defect `refuse_noun_parameters` was written against.
    #   `layer` -- documented locally above, with what it means for a LAK
    #              connection. Inheriting it too put two `layer` entries in one
    #              Parameters block.
    inherit_map_docs(GroupLakConnections.map, keep=set(NOUN_MAP_PARAMS))


_inherit_group_lak_docs()
