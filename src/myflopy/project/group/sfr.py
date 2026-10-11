"""Grouped SFR budget/stage results + namespace."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from myflopy.modflow.mf6.hydraulics import unit_labels
from myflopy.modflow.mf6.package_explorer import (
    SfrStageResultsExplorer,
    _exchange_colorscale,
    _symmetric_color_limit,
    build_group_input_compare_map_payload,
    build_sfr_budget_result_table,
    build_sfr_q_map_payload,
    build_sfr_stage_result_table,
    get_default_group_compare_colorscale,
)
from myflopy.modflow.mf6.package_surface_water import join_sfr_stage
from myflopy.modflow.utils.datatypes.hover import compare_hover, dimension_units, sfr_hover
from myflopy.project.group._shared import (
    _ensure_group_map_compatible,
    _filter_group_input_table,
    _reduce_to_period_end,
    _resolve_group_compare_target,
    _stable_compare_keys,
)
from myflopy.project.group.results import GroupCellPackageResults, GroupCellPackageResultsNamespace
from myflopy.project.group.spatial import _GroupSpatialView

if TYPE_CHECKING:
    from myflopy.project.group.core import ModelGroup


class GroupSfrBudgetResults(GroupCellPackageResults):
    """Grouped SFR exchange accessor with reach-length-normalized maps."""

    def __init__(self, group: ModelGroup, *, budget_text: str = "SFR", value_name: str = "q_gwf"):
        """Bind a grouped SFR exchange accessor (defaults to the ``SFR`` term's ``q_gwf``).

        The SFR exchange is aquifer-referenced, so the emitted column is
        ``q_gwf`` (accessor stays ``results.q``).
        """

        super().__init__(group, "sfr", budget_text=budget_text, value_name=value_name)

    def get(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """Return aligned SFR exchange rows, including ``q_per_length``.

        Each member's ``build_sfr_budget_result_table``: one row per reach per
        saved step of the cell-budget file's ``SFR`` record, every saved step
        kept, plus a ``model`` column. A member with no SFR package adds no rows;
        one with SEVERAL raises ``SeveralPackagesError`` -- each numbers its own
        reaches (ledger 352), and the group view takes no ``package=``.

        Parameters
        ----------
        model_name : str or None, optional
            Keep one member's rows, by its name in the group. Default ``None``:
            every member. A name that is not in the group matches nothing and
            returns an empty frame.
        per : int or None, optional
            Keep one stress period, zero-based: a FILTER, so every step the
            period saved is kept (:meth:`compare` and the maps read its last).
            Default ``None``: every period.
        layer : int or list of int or None, optional
            Keep the reaches in these zero-based layers (an int, list or tuple; a
            NumPy array raises ``TypeError``). Default ``None``: every layer.
        cells : int or list of int or None, optional
            Keep the reaches on these zero-based cell numbers within a layer.
            Default ``None``: every cell.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``package``, ``pname``, ``kstpkper``, ``per``, ``reach``,
            ``layer``, ``cell``, ``node``/``node2``, the reach geometry (``rlen``
            and ``distance_start/mid/end``, model length), ``q_gwf`` -- MF6's raw
            sign, flow FROM the reach TO the aquifer, so NEGATIVE where the reach
            gains (model length^3/time) -- and ``q_per_length`` = ``q_gwf / rlen``
            (length^2/time).
        """

        frames = []
        for current_model_name, model in self.group.models.items():
            frame = build_sfr_budget_result_table(
                model,
                budget_text=self.budget_text,
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
        """Compare grouped SFR exchange rows against the reference model.

        Each member is first reduced to its own LAST saved step of every period
        (``_reduce_to_period_end``, ledger 337), so members of different ``nstp``
        align by stress period. Rows are paired with the reference's on the
        integer and text identity columns -- ``package``, ``pname``, ``per``,
        ``reach``, ``layer``, ``cell``, ``node``, ``node2`` -- by an inner join;
        the float reach geometry is left out of the key and of the result, so a
        reach whose length differs by rounding still pairs.

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
            Compare only the reaches in these zero-based layers. Default
            ``None``: every layer.
        cells : int or list of int or None, optional
            Compare only these zero-based cell numbers within a layer. Default
            ``None``: every cell.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``reference_model``, the join keys, ``q_gwf``,
            ``reference_q_gwf``, ``q_gwf_diff``, and the same three for
            ``q_per_length`` -- each difference member minus reference, in MF6's
            sign. Empty when :meth:`get` is.
        """

        data = self.get(per=per, layer=layer, cells=cells)
        if data.empty:
            return pd.DataFrame()
        data = _reduce_to_period_end(data)

        reference = self.group.reference
        q = self.value_name  # aquifer-referenced exchange column (q_gwf)
        value_columns = [q, "q_per_length"] if "q_per_length" in data.columns else [q]
        key_columns = _stable_compare_keys(data, value_columns)
        ref = data[data["model"] == reference][key_columns + value_columns].copy()
        rename_map = {column: f"reference_{column}" for column in value_columns}
        ref = ref.rename(columns=rename_map)
        comp = data[data["model"] != reference].merge(ref, on=key_columns, how="inner")
        comp["reference_model"] = reference
        comp[f"{q}_diff"] = comp[q].astype(float) - comp[f"reference_{q}"].astype(float)
        if "q_per_length" in value_columns:
            comp["q_per_length_diff"] = (
                comp["q_per_length"].astype(float) - comp["reference_q_per_length"].astype(float)
            )
        ordered = [
            "model",
            "reference_model",
            *key_columns,
            q,
            f"reference_{q}",
            f"{q}_diff",
        ]
        if "q_per_length" in value_columns:
            ordered.extend(["q_per_length", "reference_q_per_length", "q_per_length_diff"])
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
        """Render one model's SFR exchange (reach-length-normalized) as a ``Choro``."""

        del agg
        target_name = self._group_target(model)
        target_model = self.group.models[target_name]
        selected = join_sfr_stage(
            target_model, self.get(model_name=target_name, per=per, layer=layer), per=per
        )
        values, hover = build_sfr_q_map_payload(
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
        kwargs.setdefault("hover_spec", sfr_hover(length=length, time=time))
        return target_model.plot.map(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            # SFR exchange is aquifer-referenced (gaining is NEGATIVE), so blue
            # sits at the negative end -- matching the single-model SFR map.
            colorscale=colorscale or _exchange_colorscale("gwf"),
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
        """Build a grouped SFR diff map using normalized exchange per unit length.

        Each member's rows of period ``per`` (its own last saved step) are
        reduced to one value per cell, ``sum(q_gwf) / sum(rlen)`` over the
        reaches in the cell -- the flow per unit reach length, in MF6's sign --
        and the map draws member minus reference on the reference's grid, for
        the cells both members have a reach in.

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
            The zero-based layer whose reaches are drawn.
        multiplier : float, default 1.0
            Scales every drawn difference -- a unit conversion or a sign flip --
            and the filled cells with them, except when no cell has a reach in
            both members (every cell is then ``fill_value`` as given).
        fill_value : float, default 0.0
            The value of a cell that lacks a reach in either member -- in one
            alone or in neither. ``float("nan")`` leaves it blank. A cell whose
            reaches have no length is ``NaN``.
        agg : str, default "sum"
            Accepted and IGNORED: the per-cell reduction is always the
            length-weighted ``sum(q_gwf) / sum(rlen)`` above, whatever is passed
            (ledger 439).
        colorscale : str or list or None, optional
            Any Plotly colorscale name or explicit stops. Default ``None``: the
            group diff scale -- red at the negative end, white at zero, blue at
            the positive end (in MF6's sign a reach gaining MORE than the
            reference's is negative, so red).
        **kwargs
            Forwarded to the reference model's ``model.plot.map(...)``: its named
            drawing options (``contours``, ``select``, ``features``, ``extent``,
            ``backend``, ...) and its Plotly trace tail. ``zmin``/``zmax`` default
            to plus/minus the largest absolute value drawn (a filled cell's
            included), ``zmid`` to 0 and ``hover_spec`` to the comparison hover,
            whose units are length²/time in the reference model's own units
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
        data = _reduce_to_period_end(self.get(per=per, layer=layer))

        def _normalized_by_cell(frame: pd.DataFrame, current_model_name: str) -> pd.DataFrame:
            """Per-cell exchange per unit reach length (``sum(q)/sum(rlen)``) for one model."""

            selected = frame[frame["model"] == current_model_name].copy()
            if selected.empty:
                return pd.DataFrame(columns=["cell", "q_per_length"])
            q = self.value_name
            selected[q] = pd.to_numeric(selected[q], errors="coerce")
            selected["rlen"] = pd.to_numeric(selected["rlen"], errors="coerce")
            grouped = selected.groupby("cell", as_index=False).agg({q: "sum", "rlen": "sum"})
            grouped["q_per_length"] = np.where(grouped["rlen"] > 0.0, grouped[q] / grouped["rlen"], np.nan)
            return grouped[["cell", "q_per_length"]]

        reference = _normalized_by_cell(data, self.group.reference).rename(
            columns={"q_per_length": "reference_q_per_length"}
        )
        target = _normalized_by_cell(data, target_name)
        comparison = target.merge(reference, on="cell", how="inner")
        comparison["model"] = target_name
        comparison["reference_model"] = self.group.reference
        comparison["per"] = int(per)
        comparison["layer"] = int(layer)
        comparison["q_per_length_diff"] = (
            comparison["q_per_length"].astype(float) - comparison["reference_q_per_length"].astype(float)
        )
        values, hover, absmax = build_group_input_compare_map_payload(
            comparison,
            ncpl=reference_model.vor.ncpl,
            value_column="q_per_length",
            diff_column="q_per_length_diff",
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
                "q_per_length",
                "q_per_length_diff",
                title="Δ stream exchange vs reference",
                # q / length is L²/T, in the drawn model's own units (ledger 434).
                units=dimension_units(*unit_labels(reference_model), {
                    "q_per_length": (2, -1),
                    "reference_q_per_length": (2, -1),
                    "q_per_length_diff": (2, -1),
                }),
                labels={"q_per_length_diff": "Δ q / length"},
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

class GroupSfrStageResults(_GroupSpatialView):
    """Grouped accessor for SFR reach stages and stage comparisons.

    Mirrors :class:`GroupLakStageResults` for streams (keyed by ``reach``),
    reading ``model.outputs.sfr.stage``. Inherits the :class:`SpatialView`
    grammar; each panel delegates to the single-model SFR stage explorer.
    """

    def __init__(self, group: ModelGroup):
        """Bind the grouped SFR reach stage-results accessor to ``group``."""

        self.group = group

    def _spatial_value_label(self) -> str:
        """The mapped quantity's label -- reach stage."""

        return "stage"

    def _spatial_map(self, *, per: int = 0, layer: int = 0, model=None, **kwargs):
        """Render one model's reach stage broadcast to reach cells."""

        target_model = self.group.models[self._group_target(model)]
        return SfrStageResultsExplorer(target_model).map(
            per=int(per), layer=int(layer), backend="plotly", **kwargs
        )

    def get(
        self,
        *,
        model_name: str | None = None,
        reach: int | None = None,
        per: int | None = None,
    ) -> pd.DataFrame:
        """Return each member's reach stages, one per reach per stress period.

        Read from ``build_sfr_stage_result_table`` -- the table the single-model
        noun draws, labelled by the stage file's own ``kstpkper`` -- with each
        period at its LAST saved step on each member's own clock
        (``_reduce_to_period_end``, ledger 337's rule for a group table). Until
        ledger 348 it decoded ``model.outputs.sfr.stage.get()`` -- a REACH-by-time
        frame -- as LAK's time-by-lake array, so ``per`` held the reach number and
        ``reach`` the saved record (ledger 301's "transposed decode"): measured on
        a four-reach stream at ``nstp=3``, ``get(per=1, reach=0)`` read reach 1 at
        the first saved step.

        Parameters
        ----------
        model_name : str or None, optional
            Keep one member, by its name in the group. Default ``None``: every
            member.
        reach : int or None, optional
            Keep one reach, zero-based. Default ``None``: every reach.
        per : int or None, optional
            Keep one stress period, zero-based -- its last saved step. Default
            ``None``: every period.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``kstpkper`` (the step read, zero-based), ``per``, ``reach``
            and ``stage`` (model length units, ``NaN`` where a reach has none --
            dry or inactive, ledger 353).
        """

        rows: list[pd.DataFrame] = []
        for current_model_name, model in self.group.models.items():
            table = build_sfr_stage_result_table(model)
            frame = table.drop_duplicates(subset=["kstpkper", "reach"]).loc[
                :, ["kstpkper", "per", "reach", "stage"]
            ]
            frame["model"] = current_model_name
            rows.append(frame)

        if not rows:
            return pd.DataFrame(columns=["model", "kstpkper", "per", "reach", "stage"])
        combined = _reduce_to_period_end(pd.concat(rows, ignore_index=True))
        combined = combined[["model", "kstpkper", "per", "reach", "stage"]]
        if model_name is not None:
            combined = combined.loc[combined["model"] == str(model_name)].copy()
        if reach is not None:
            combined = combined.loc[combined["reach"] == int(reach)].copy()
        if per is not None:
            combined = combined.loc[combined["per"] == int(per)].copy()
        return combined.reset_index(drop=True)

    def compare(
        self,
        *,
        model_name: str | None = None,
        reach: int | None = None,
        per: int | None = None,
    ) -> pd.DataFrame:
        """Compare SFR reach stages against the reference model, period by period.

        Rows are aligned on the stress PERIOD and the reach, each model read at its
        own last saved step of the period (:meth:`get`); until ledger 348 the
        "period" was a reach number and the "reach" a saved record.

        Parameters
        ----------
        model_name : str or None, optional
            Keep one compared (non-reference) member. Default ``None``: every one.
        reach : int or None, optional
            Keep one reach, zero-based. Default ``None``: every reach.
        per : int or None, optional
            Keep one stress period, zero-based. Default ``None``: every period.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``reference_model``, ``per``, ``reach``, ``stage``,
            ``reference_stage`` and ``stage_diff`` (model minus reference).
        """

        columns = [
            "model", "reference_model", "per", "reach",
            "stage", "reference_stage", "stage_diff",
        ]
        data = self.get(reach=reach, per=per)
        if data.empty:
            return pd.DataFrame(columns=columns)

        reference = self.group.reference
        ref = (
            data.loc[data["model"] == reference, ["per", "reach", "stage"]]
            .rename(columns={"stage": "reference_stage"})
            .copy()
        )
        comp = data.loc[data["model"] != reference].merge(ref, on=["per", "reach"], how="inner")
        comp["reference_model"] = reference
        comp["stage_diff"] = comp["stage"].astype(float) - comp["reference_stage"].astype(float)
        if model_name is not None:
            comp = comp.loc[comp["model"] == str(model_name)].copy()
        return comp[columns]


class GroupSfrResultsNamespace(GroupCellPackageResultsNamespace):
    """Namespace for grouped SFR result accessors."""

    def _field_names(self):
        """The mappable grouped SFR result fields: exchange ``q`` and ``stage``."""

        return ["q", "stage"]

    @property
    def q(self) -> GroupSfrBudgetResults:
        """Return grouped SFR exchange helpers with normalized map behavior."""

        return self._result_accessor

    @property
    def stage(self) -> GroupSfrStageResults:
        """Return grouped SFR reach-stage helpers."""

        return GroupSfrStageResults(self._result_accessor.group)


