"""Results-tier difference for :class:`~myflopy.project.model_diff.ModelDiff`.

Where the other tiers compare how models are *set up* (inputs, config, connection
geometry), the results tier compares what the runs *computed* -- the outputs.
This requires both models to have completed runs (heads and the listing budget on
disk).

Two questions at once ("both" framing):

* **regression** -- with the same inputs the outputs should match; a
  ``within_tolerance`` flag (absolute + relative tolerance) says whether they do.
* **scenario** -- with different inputs, *how much* do the outputs differ, and
  *where / when* (max |diff| and its location).

These accessors surface in the unified tree as ``diff.hds`` (heads), ``diff.bud``
(volumetric budget), and ``diff.packages.<pkg>.results`` (per-package cell
budgets, UZF gwrch/sat, lake/SFR q + stage, MVR). All stress-period/time
references are zero-based.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

from myflopy._logging import get_logger
from myflopy.modflow.mf6.package_budget import (
    build_lak_budget_term_table,
    build_sfr_budget_term_table,
)
from myflopy.modflow.mf6.package_explorer import (
    DiffSpatialView,
    FieldMappable,
    LineProfileView,
)
from myflopy.modflow.mf6.package_registry import (
    _PACKAGE_EXPLORER_SPECS,
    model_packages_of_type,
)
from myflopy.modflow.utils.datatypes.xsections import XSection
from myflopy.project.group._shared import _reduce_to_period_end, _stable_compare_keys
from myflopy.project.group.lak import GroupLakStageResults
from myflopy.project.group.results import GroupCellPackageResults
from myflopy.project.group.sfr import GroupSfrStageResults

logger = get_logger(__name__)

# Default tolerances for the ``within_tolerance`` test: |diff| <= atol + rtol*|ref|.
_DEFAULT_ATOL = 1e-3
_DEFAULT_RTOL = 1e-3


def _within_tolerance(diff, reference, atol: float, rtol: float) -> np.ndarray:
    """Elementwise ``|diff| <= atol + rtol*|reference|``."""

    diff = np.asarray(diff, dtype=float)
    reference = np.asarray(reference, dtype=float)
    return np.abs(diff) <= (atol + rtol * np.abs(reference))


class _ResultDiffBase:
    """Shared target resolution for results accessors."""

    def __init__(self, diff):
        """Bind a results-tier accessor to a group via its ``ModelDiff``."""

        self._diff = diff
        self.group = diff.group

    def _targets(self, model_name=None) -> list[str]:
        """The non-reference model names to compare (or just ``model_name``, validated)."""

        if model_name is not None:
            name = str(model_name)
            if name not in self.group.models:
                raise KeyError(f"Model {name!r} is not in the group.")
            if name == self.group.reference:
                raise ValueError("The reference model cannot be diffed against itself.")
            return [name]
        return [name for name in self.group.models if name != self.group.reference]


class HeadsResultDiff(_ResultDiffBase, LineProfileView, DiffSpatialView):
    """Head-difference results vs the reference model (per cell / layer / time).

    ``map`` / ``mosaic`` / ``animate`` render Δhead (model - reference) --
    ``diff.hds.map("F9b", per=8, layer=1)``, ``diff.hds.mosaic(by="model")``.

    Configured by the three class attributes below, so the transport twins
    (:class:`ConcResultDiff`, :class:`TempResultDiff`) are the same class with a
    different field rather than a copy. Mirrors ``_GroupFieldView`` on the group
    side, which they all delegate to.
    """

    #: the ``ModelGroup`` accessor this leaf diffs (``group.hds``/``conc``/``temp``)
    group_attribute = "hds"
    #: the value column in that accessor's compare table
    value_column = "elev"
    #: short label for the differenced quantity
    value_label = "head"

    def _noun_path(self) -> str:
        """``diff.hds`` / ``diff.conc`` / ``diff.temp``, for messages."""

        return f"diff.{self.group_attribute}"

    @property
    def _field(self):
        """The grouped field accessor backing this diff leaf."""

        return getattr(self.group, self.group_attribute)

    @property
    def _reference_column(self) -> str:
        """The compare table's reference column for this field."""

        return f"reference_{self.value_column}"

    def get(self, *, model_name=None, per=None, layer=None, cells=None) -> pd.DataFrame:
        """Return aligned differences: value / reference value / ``diff``.

        The group's compare (``group.hds``/``conc``/``temp``): each member's saved
        values reduced to its OWN last saved step of each stress period, then
        inner-joined to the reference's on ``(per, layer, cell)`` -- so members
        with different ``nstp`` line up by period, and a cell or period one side
        lacks is absent. ``diff`` is model minus reference, in the field's unit
        (head: the model's length unit; ``conc``: the GWT model's concentration
        unit; ``temp``: the GWE model's temperature unit).

        A cell MF6 wrote as dry (-1e30) or inactive (+1e30) has no value, as on
        ``model.hds.get()`` (``headsplus._is_sentinel``, ledger 315): its value
        is ``NaN``, the row kept, and so is its ``diff``. Until ledger 484
        the sentinel entered ``diff`` as a number -- one dry cell against a wet
        reference read a ``diff`` of -1e30 (``docs/reproducers/
        diff_results_tier.py``); :meth:`summary` counts such a row in
        ``n_one_sided``.

        Parameters
        ----------
        model_name : str or None, optional
            Keep one compared model, by its name in the group. Default ``None``
            keeps every non-reference model. A name not in the group raises
            ``KeyError``; the reference's raises ``ValueError``. Every member is
            read either way.
        per : int or None, optional
            Keep one zero-based stress period, each member read at its own last
            saved step of it (:func:`~myflopy.modflow.mf6.period_end.last_saved_step`,
            ledger 337). A period a member's file saved nothing in raises
            ``ValueError`` naming the periods it has. Default ``None`` keeps every
            period.
        layer : int or list of int or None, optional
            Keep these zero-based layers. A ``list`` is several layers; any other
            value goes through ``int()``, so a TUPLE raises ``TypeError``. Default
            ``None`` keeps every layer.
        cells : list of int or None, optional
            Keep these zero-based in-layer cells (any iterable of ints). A bare
            ``int`` raises ``TypeError`` -- pass ``[5]``. Default ``None`` keeps
            every cell.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``reference_model``, ``kstpkper`` (the compared member's own
            step), ``per``, ``layer``, ``cell``, the value column (``elev``,
            ``conc`` or ``temp``), ``reference_<value>`` and ``diff`` -- ``NaN``
            where a side has no value.

        Raises
        ------
        KeyError
            If ``model_name`` names no model in the group.
        ValueError
            If ``model_name`` is the reference, or ``per`` names a period a
            member saved nothing in.
        """

        frame = self._field.compare(per=per, layer=layer, cells=cells)
        if model_name is not None:
            targets = self._targets(model_name)
            frame = frame[frame["model"].isin(targets)]
        return frame.reset_index(drop=True)

    # -- spatial-view hooks (delta head maps) ---------------------------------
    def _spatial_map(self, *, per=0, layer=0, model=None, **kwargs):
        """Δhead choropleth (compared model - reference) for one model."""

        return self._field.compare_map(model_name=model, per=per, layer=layer, **kwargs)

    def _spatial_models(self):
        """The non-reference model names -- one Δhead panel each."""

        return self._targets()

    def _spatial_reference_model(self):
        """The reference model whose grid/layers frame the Δhead maps."""

        return self.group.models[self.group.reference]

    def _spatial_value_label(self):
        """The mapped quantity's label."""

        return self.value_label

    # -- series hooks: plot() draws mean Δhead by period, one line per model --
    def _series_value_column(self, frame) -> str:
        """The column ``plot()`` draws: the aligned head difference (``diff``)."""

        return "diff"

    def _series_default_agg(self) -> str:
        """Collapse cells within a line by mean (averaging Δhead, not summing)."""

        return "mean"

    def _sections(
        self,
        model=None,
        *,
        line=None,
        cells: int | list[int] | None = None,
        per: int | None = None,
        layer: int | list[int] = 0,
        **kwargs,
    ):
        """Return reference + compared-model :class:`XSection` objects.

        The grammar's ``xs`` verbs overlay the reference model's profile with
        each target model's along the same ``line`` (or ``cells`` path), so
        head differences read directly off the profiles; ``model`` narrows to
        one compared model.
        """

        names = [self.group.reference, *self._targets(model)]
        return {
            name: XSection(
                model=self.group.models[name],
                section_name=name,
                line=line,
                cells=cells,
                per=per,
                layer=layer,
                **kwargs,
            )
            for name in names
        }

    def summary(
        self,
        *,
        model_name=None,
        atol: float = _DEFAULT_ATOL,
        rtol: float = _DEFAULT_RTOL,
        per=None,
        layer=None,
    ) -> pd.DataFrame:
        """One row per model: max/mean |Δhead|, RMSE, where/when the max is, and
        whether every cell is within tolerance.

        Over the rows :meth:`get` returns, as a stage diff counts a dry lake
        (ledger 329). A cell dry or inactive in a model has no value, so the
        statistics run over the rows where BOTH models have one, and ``n`` counts
        those; ``n_one_sided`` counts the rows where exactly one does -- a cell
        dry in one model and wet in the other, a difference no number measures --
        and any such row makes ``within_tolerance`` False. A row with no value on
        either side agrees and counts in neither. Until ledger 484 the
        sentinel's +-1e30 set ``max_abs_diff``, ``mean_abs_diff`` and ``rmse``.

        Parameters
        ----------
        model_name : str or None, optional
            Summarize one compared model, by its name in the group. Default
            ``None`` summarizes every non-reference model. A name not in the group
            raises ``KeyError``; the reference's raises ``ValueError``.
        atol : float, optional
            Absolute tolerance, in the field's unit (head: the model's length unit;
            ``diff.conc``/``diff.temp``: concentration/temperature). A row is
            within tolerance when ``|diff| <= atol + rtol * |reference value|``,
            and the model when every row is. Default ``1e-3``.
        rtol : float, optional
            Relative tolerance, dimensionless, multiplying ``|reference value|`` in
            the same test. Default ``1e-3``.
        per : int or None, optional
            One zero-based stress period, each member at its own last saved step
            of it, as in :meth:`get` (a period a member did not save raises
            ``ValueError``). Default ``None`` uses every period.
        layer : int or list of int or None, optional
            These zero-based layers, as in :meth:`get` (a tuple raises
            ``TypeError``). Default ``None`` uses every layer.

        Returns
        -------
        pandas.DataFrame
            One row per compared model: ``model``, ``n`` (rows where both have a
            value), ``n_one_sided``, ``max_abs_diff``, ``mean_abs_diff``,
            ``rmse``, ``argmax_cell``, ``argmax_layer``, ``argmax_kstpkper``
            (where the largest ``|diff|`` is; the first such row on a tie) and
            ``within_tolerance``. With no row where both have a value, the three
            statistics are ``NaN``, the two argmax cell columns ``-1`` and
            ``argmax_kstpkper`` ``None``. Empty when there is nothing to compare.
        """

        data = self.get(model_name=model_name, per=per, layer=layer)
        columns = [
            "model", "n", "n_one_sided", "max_abs_diff", "mean_abs_diff", "rmse",
            "argmax_cell", "argmax_layer", "argmax_kstpkper", "within_tolerance",
        ]
        if data.empty:
            return pd.DataFrame(columns=columns)
        rows = []
        for name, sub in data.groupby("model"):
            has_value = np.isfinite(sub[self.value_column].to_numpy(dtype=float))
            has_reference = np.isfinite(sub[self._reference_column].to_numpy(dtype=float))
            one_sided = int((has_value ^ has_reference).sum())
            compared = sub.loc[has_value & has_reference]
            diff = compared["diff"].to_numpy(dtype=float)
            reference = compared[self._reference_column].to_numpy(dtype=float)
            abs_diff = np.abs(diff)
            imax = int(np.argmax(abs_diff)) if abs_diff.size else 0
            rows.append(
                {
                    "model": name,
                    "n": int(abs_diff.size),
                    "n_one_sided": one_sided,
                    "max_abs_diff": float(abs_diff.max()) if abs_diff.size else float("nan"),
                    "mean_abs_diff": float(abs_diff.mean()) if abs_diff.size else float("nan"),
                    "rmse": float(np.sqrt(np.mean(diff**2))) if diff.size else float("nan"),
                    "argmax_cell": int(compared["cell"].iloc[imax]) if abs_diff.size else -1,
                    "argmax_layer": int(compared["layer"].iloc[imax]) if abs_diff.size else -1,
                    "argmax_kstpkper": compared["kstpkper"].iloc[imax] if abs_diff.size else None,
                    "within_tolerance": bool(
                        one_sided == 0
                        and np.all(_within_tolerance(diff, reference, atol, rtol))
                    ),
                }
            )
        return pd.DataFrame(rows, columns=columns)


class ConcResultDiff(HeadsResultDiff):
    """Δconcentration results vs the reference model (``diff.conc``).

    Same class as :class:`HeadsResultDiff`, pointed at the group's GWT
    concentration accessor. ``diff.conc.map("variant", per=8)`` is the public
    route to a Δconc choropleth.
    """

    group_attribute = "conc"
    value_column = "conc"
    value_label = "conc"


class TempResultDiff(HeadsResultDiff):
    """Δtemperature results vs the reference model (``diff.temp``)."""

    group_attribute = "temp"
    value_column = "temp"
    value_label = "temp"


class BudgetResultDiff(_ResultDiffBase):
    """Overall (volumetric listing) budget difference vs the reference model.

    Compares the incremental per-timestep budget term by term, aligned on
    simulation time (``totim``) so models with the same physics but different
    time discretization still line up where their timesteps coincide.
    """

    #: FloPy's bookkeeping columns in ``get_incremental()``: the time step, its
    #: stress period and its length -- what a row IS, not a budget term. They were
    #: melted and differenced like terms until ledger 485, so two runs differing
    #: only in ``nstp`` read ``tslen`` at -50% and ``time_step`` out of tolerance.
    _BOOKKEEPING = ("time_step", "stress_period", "tslen")

    def _model_budget(self, model_name: str) -> pd.DataFrame:
        """One model's incremental listing budget, melted to ``totim`` / ``term`` / ``value``.

        Carries the step's ``per`` and ``kstpkper`` (zero-based) when FloPy's table
        has them, and never melts a bookkeeping column into a term.
        """

        # ``budget_incremental`` is a property on real models (returns a
        # DataFrame); tolerate a callable too so fakes/other readers work.
        raw = self.group.models[model_name].budget_incremental
        raw = raw() if callable(raw) else raw
        frame = pd.DataFrame(raw).copy()
        if "totim" not in frame.columns:
            frame = frame.reset_index().rename(columns={"index": "totim"})
        if {"time_step", "stress_period"} <= set(frame.columns):
            frame["per"] = frame["stress_period"].astype(int)
            frame["kstpkper"] = [
                (int(kstp), int(kper))
                for kstp, kper in zip(frame["time_step"], frame["stress_period"], strict=True)
            ]
        else:
            frame["per"] = None
            frame["kstpkper"] = None
        keys = ("totim", "per", "kstpkper", *self._BOOKKEEPING)
        terms = [
            column
            for column in frame.columns
            if column not in keys and pd.api.types.is_numeric_dtype(frame[column])
        ]
        return frame.melt(
            id_vars=["totim", "per", "kstpkper"], value_vars=terms,
            var_name="term", value_name="value",
        )

    def get(self, *, model_name=None) -> pd.DataFrame:
        """Return aligned per-term budget differences per timestep.

        Each member's incremental listing budget (``model.budget_incremental``,
        FloPy's ``get_incremental()``: the rate for each time step, in the model's
        volume per time unit) is melted to one row per ``(totim, term)`` and
        inner-joined to the reference's on an equal ``totim``, so a time step only
        one model's listing file holds is absent. A "term" is every numeric
        column of that table bar ``totim`` and FloPy's bookkeeping columns -- the
        budget's ``<TERM>_IN``/``<TERM>_OUT``, ``TOTAL_IN``/``TOTAL_OUT``,
        ``IN-OUT`` and ``PERCENT_DISCREPANCY``. ``time_step``, ``stress_period``
        and ``tslen`` say what a row IS; they label it (``per``, ``kstpkper``)
        and are not differenced -- they were, as terms, until ledger 485
        (``docs/reproducers/diff_results_tier.py``).

        Parameters
        ----------
        model_name : str or None, optional
            Compare one model, by its name in the group. Default ``None`` compares
            every non-reference model. A name not in the group raises ``KeyError``;
            the reference's raises ``ValueError``.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``totim``, ``per`` and ``kstpkper`` (the compared model's
            zero-based stress period and ``(kstp, kper)`` of that time; ``None``
            where its table carries no such columns), ``term``,
            ``reference_value``, ``model_value`` and ``diff`` (model minus
            reference).
        """

        reference = self._model_budget(self.group.reference).drop(columns=["per", "kstpkper"])
        columns = [
            "model", "totim", "per", "kstpkper", "term",
            "reference_value", "model_value", "diff",
        ]
        frames = []
        for name in self._targets(model_name):
            model = self._model_budget(name)
            merged = reference.merge(
                model, on=["totim", "term"], suffixes=("_reference", "_model"), how="inner"
            )
            merged["model"] = name
            merged["diff"] = merged["value_model"].astype(float) - merged[
                "value_reference"
            ].astype(float)
            frames.append(
                merged.rename(
                    columns={
                        "value_reference": "reference_value",
                        "value_model": "model_value",
                    }
                )[columns]
            )
        if not frames:
            return pd.DataFrame(columns=columns)
        return pd.concat(frames, ignore_index=True)

    def summary(
        self,
        *,
        model_name=None,
        atol: float = _DEFAULT_ATOL,
        rtol: float = _DEFAULT_RTOL,
    ) -> pd.DataFrame:
        """One row per (model, budget term): totals, percent change, worst
        per-timestep difference, and whether every timestep is within tolerance.

        Over :meth:`get`'s rows: budget terms only. FloPy's ``time_step``/
        ``stress_period``/``tslen`` had rows too until ledger 485, so two
        members differing only in ``nstp`` read ``tslen`` at -50% and
        ``time_step`` out of tolerance.

        Parameters
        ----------
        model_name : str or None, optional
            Summarize one model, by its name in the group. Default ``None``
            summarizes every non-reference model. A name not in the group raises
            ``KeyError``; the reference's raises ``ValueError``.
        atol : float, optional
            Absolute tolerance, in each term's own unit (a flow term's volume per
            time; ``PERCENT_DISCREPANCY``'s percent). A time step is within
            tolerance when ``|diff| <= atol + rtol * |reference value|``, and the
            term when every aligned step is. Default ``1e-3``.
        rtol : float, optional
            Relative tolerance, dimensionless, multiplying ``|reference value|`` in
            the same test. Default ``1e-3``.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``term``, ``reference_total`` and ``model_total`` (the SUM of
            the aligned steps' rates, not a volume), ``diff_total``,
            ``pct_change`` (``NaN`` when the reference total is 0),
            ``max_abs_diff`` (worst step) and ``within_tolerance``; sorted by model
            and term. Empty when there is nothing to compare.
        """

        data = self.get(model_name=model_name)
        columns = [
            "model", "term", "reference_total", "model_total", "diff_total",
            "pct_change", "max_abs_diff", "within_tolerance",
        ]
        if data.empty:
            return pd.DataFrame(columns=columns)
        rows = []
        for (name, term), sub in data.groupby(["model", "term"]):
            reference_total = float(sub["reference_value"].sum())
            model_total = float(sub["model_value"].sum())
            diff_total = model_total - reference_total
            pct = (diff_total / reference_total * 100.0) if reference_total != 0 else float("nan")
            rows.append(
                {
                    "model": name,
                    "term": term,
                    "reference_total": reference_total,
                    "model_total": model_total,
                    "diff_total": diff_total,
                    "pct_change": pct,
                    "max_abs_diff": float(sub["diff"].abs().max()),
                    "within_tolerance": bool(
                        np.all(
                            _within_tolerance(
                                sub["diff"].to_numpy(dtype=float),
                                sub["reference_value"].to_numpy(dtype=float),
                                atol,
                                rtol,
                            )
                        )
                    ),
                }
            )
        return pd.DataFrame(rows, columns=columns).sort_values(
            ["model", "term"]
        ).reset_index(drop=True)


class CellBudgetResultDiff(_ResultDiffBase, DiffSpatialView):
    """Per-package cell-budget difference vs the reference model.

    Wraps a grouped cell-budget accessor (GHB/DRN/… leakage, SFR/LAK exchange,
    UZF recharge/saturation) and adds Δ summary stats + a within-tolerance flag on
    its value column (``q``, ``gwrch``, ``sat``, …), plus the unified delta
    ``map`` / ``mosaic`` / ``animate`` grammar.
    """

    def __init__(self, diff, accessor):
        """Wrap a grouped cell-budget ``accessor``, taking its package + value names."""

        super().__init__(diff)
        self._accessor = accessor
        self.package_name = accessor.package_name
        self.value_name = accessor.value_name

    def get(self, *, model_name=None, per=None, layer=None, cells=None) -> pd.DataFrame:
        """Return aligned per-cell budget differences (value / reference / diff).

        The wrapped group accessor's ``compare``. The value column, and how rows
        are aligned, depend on which accessor this is:

        - a BC package's ``results.q`` (``rch``/``chd``/``drn``/``ghb``/``riv``/
          ``wel``/``evt``), SFR's ``results.q`` and UZF's ``gwrch``/``sat``: each
          member reduced to its OWN last saved step of each stress period, then
          joined to the reference on every non-float column but ``kstpkper``
          (``per``, ``layer``, ``cell``, ``node``, ``pname``, a reach ...), so
          members with different ``nstp`` line up by period. ``q_gwf`` is MF6's
          sign, positive INTO the aquifer, in model volume per time (SFR adds
          ``q_per_length``); ``gwrch`` is UZF's recharge to the water table, same
          sign and unit. ``sat`` reads MF6's ``DATA-SAT`` (NPF's, under
          ``SAVE_SATURATION``, for every cell): the saturated fraction, 0-1,
          from the record's ``SAT`` auxiliary (its placeholder flow is dropped;
          it was renamed ``sat`` beside it until ledger 471, and this
          raised).
        - LAK's ``results.q``: reduced the same way, and joined on the
          connection's identity -- ``package``, ``per``, ``lake``, ``iconn``,
          ``layer``, ``cell``, ``claktype`` -- never on a float. ``q_lake`` is the
          lake's frame, positive INTO the lake; ``q_per_area`` and MF6's wetted
          ``flow_area`` are differenced beside it. Until ledger 485 it was
          joined on the exact ``kstpkper`` and on every float column: a member at
          ``nstp=2`` was read at the FIRST step of each period against a
          reference at ``nstp=1``, and a connection whose wetted area moved had
          no partner (2 of 72 at period 3, ``docs/reproducers/
          diff_results_tier.py``).
        - MVR (:class:`MvrResultDiff`): the reduced form, column ``q``.

        Parameters
        ----------
        model_name : str or None, optional
            Keep one compared model, by its name in the group. Default ``None``
            keeps every non-reference model. A name not in the group raises
            ``KeyError`` and the reference's ``ValueError`` -- only when there are
            rows to filter; an empty comparison comes back unchecked.
        per : int or None, optional
            Keep one zero-based stress period, each member at its last saved step
            of it. A period no member saved gives an empty frame, not an error.
            Default ``None`` keeps every period.
        layer : int or sequence of int or None, optional
            Keep these zero-based layers. Default ``None`` keeps every layer.
        cells : int or sequence of int or None, optional
            Keep these zero-based in-layer cells. Default ``None`` keeps every
            cell.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``reference_model``, the join columns, ``<value>``,
            ``reference_<value>`` and ``<value>_diff`` (model minus reference).
            Empty, with no columns, when no row is selected -- no member has a
            record of the term, or the filters keep none; empty WITH the columns
            when rows were read but none pairs with a reference row.
        """

        frame = self._accessor.compare(per=per, layer=layer, cells=cells)
        if frame.empty:
            return frame
        if model_name is not None:
            frame = frame[frame["model"].isin(self._targets(model_name))]
        return frame.reset_index(drop=True)

    # -- series hook: plot() draws the Δ column, one line per model -----------
    def _series_value_column(self, frame) -> str:
        """The column ``plot()`` draws: the aligned ``<value>_diff``."""

        return f"{self.value_name}_diff"

    def summary(
        self,
        *,
        model_name=None,
        atol: float = _DEFAULT_ATOL,
        rtol: float = _DEFAULT_RTOL,
        per=None,
        layer=None,
        cells=None,
    ) -> pd.DataFrame:
        """One row per model: max/mean/RMSE of the per-cell flow Δ, where the max
        is, and whether every cell is within tolerance.

        Over the rows :meth:`get` returns (its alignment rules apply).

        Parameters
        ----------
        model_name : str or None, optional
            Summarize one compared model, by its name in the group. Default
            ``None`` summarizes every non-reference model. Checked as in
            :meth:`get`: ``KeyError``/``ValueError`` only when there are rows.
        atol : float, optional
            Absolute tolerance, in the value's unit (a flow's model volume per
            time; ``sat`` raises before it is used, :meth:`get`). A row is within
            tolerance when
            ``|<value>_diff| <= atol + rtol * |reference value|``, and the model
            when every row is. Default ``1e-3``.
        rtol : float, optional
            Relative tolerance, dimensionless, multiplying ``|reference value|`` in
            the same test. Default ``1e-3``.
        per : int or None, optional
            One zero-based stress period, as in :meth:`get`. Default ``None`` uses
            every period.
        layer : int or sequence of int or None, optional
            These zero-based layers, as in :meth:`get`. Default ``None`` uses every
            layer.
        cells : int or sequence of int or None, optional
            These zero-based in-layer cells, as in :meth:`get`. Default ``None``
            uses every cell.

        Returns
        -------
        pandas.DataFrame
            One row per compared model: ``model``, ``package``, ``value`` (the
            value column's name), ``n``, ``max_abs_diff``, ``mean_abs_diff``,
            ``rmse``, ``argmax_cell`` (the first row with the largest ``|diff|``;
            ``-1`` without a ``cell`` column) and ``within_tolerance``. Empty when
            there is nothing to compare.
        """

        data = self.get(model_name=model_name, per=per, layer=layer, cells=cells)
        diff_column = f"{self.value_name}_diff"
        reference_column = f"reference_{self.value_name}"
        columns = [
            "model", "package", "value", "n", "max_abs_diff", "mean_abs_diff",
            "rmse", "argmax_cell", "within_tolerance",
        ]
        if data.empty or diff_column not in data.columns:
            return pd.DataFrame(columns=columns)
        rows = []
        for name, sub in data.groupby("model"):
            diff = sub[diff_column].to_numpy(dtype=float)
            reference = (
                sub[reference_column].to_numpy(dtype=float)
                if reference_column in sub.columns
                else np.zeros_like(diff)
            )
            abs_diff = np.abs(diff)
            imax = int(np.argmax(abs_diff)) if abs_diff.size else 0
            rows.append(
                {
                    "model": name,
                    "package": self.package_name,
                    "value": self.value_name,
                    "n": int(abs_diff.size),
                    "max_abs_diff": float(abs_diff.max()) if abs_diff.size else 0.0,
                    "mean_abs_diff": float(abs_diff.mean()) if abs_diff.size else 0.0,
                    "rmse": float(np.sqrt(np.mean(diff**2))) if diff.size else 0.0,
                    "argmax_cell": (
                        int(sub["cell"].iloc[imax])
                        if "cell" in sub.columns and abs_diff.size
                        else -1
                    ),
                    "within_tolerance": bool(
                        np.all(_within_tolerance(diff, reference, atol, rtol))
                    ),
                }
            )
        return pd.DataFrame(rows, columns=columns)

    # -- spatial-view hooks (delta per-cell flux maps) ------------------------
    def _spatial_map(self, *, per=0, layer=0, model=None, **kwargs):
        """Δ leakage/exchange choropleth for one compared model."""

        return self._accessor.compare_map(model_name=model, per=per, layer=layer, **kwargs)

    def _spatial_models(self):
        """The non-reference model names -- one Δflux panel each."""

        return self._targets()

    def _spatial_reference_model(self):
        """The reference model whose grid/layers frame the Δflux maps."""

        return self.group.models[self.group.reference]

    def _spatial_periods(self):
        """Stress periods present in this term's data (delegates to the accessor)."""

        return self._accessor._spatial_periods()

    def _spatial_layers(self):
        """Layers present in this term's data (delegates to the accessor)."""

        return self._accessor._spatial_layers()

    def _spatial_value_label(self):
        """The mapped quantity's label -- this term's value name (``q``, ``gwrch``, ...)."""

        return self.value_name


# Packages with a per-cell budget diff via group.packages.<pkg>.results.q.
# Derived from the registry (plan 4.7.3) -- adding a package to the descriptor
# with ``results_diffable=True`` is now the only edit needed.
_CELL_BUDGET_PACKAGES = tuple(
    name
    for name, spec in _PACKAGE_EXPLORER_SPECS.items()
    if spec.tiers.results_diffable
)


class CellResultsDiffNamespace(FieldMappable):
    """``diff.packages.<pkg>.results`` for cell-budget BC packages.

    One field, ``q`` (the package-groundwater exchange Δ), matching the
    single-model and group results namespaces -- so ``results.map()`` /
    ``results.q.map("F9b")`` follow the unified grammar with Δ content.
    """

    _default_field = "q"

    def __init__(self, diff, package_name: str):
        """Bind a cell-budget results-diff namespace to one package."""

        self._diff = diff
        self.group = diff.group
        self.package_name = str(package_name).lower()

    def _field_names(self):
        """The single mappable results-diff field: exchange ``q``."""

        return ["q"]

    @property
    def q(self) -> CellBudgetResultDiff:
        """Per-cell exchange difference vs the reference model."""

        try:
            accessor = getattr(self.group.packages, self.package_name)
        except AttributeError as exc:
            raise AttributeError(
                f"No grouped results accessor for package {self.package_name!r}."
            ) from exc
        cell = getattr(getattr(accessor, "results", None), "q", None)
        if cell is None:
            raise AttributeError(
                f"Package {self.package_name!r} has no cell-budget results."
            )
        return CellBudgetResultDiff(self._diff, cell)


class UzfResultsDiffNamespace(FieldMappable):
    """``diff.packages.uzf.results`` -- fields ``gwrch`` (default) and ``sat``."""

    _default_field = "gwrch"

    def __init__(self, diff):
        """Bind the UZF results-diff namespace to a group via its ``ModelDiff``."""

        self._diff = diff
        self.group = diff.group

    def _field_names(self):
        """The mappable UZF results-diff fields: ``gwrch`` and ``sat``."""

        return ["gwrch", "sat"]

    def _uzf_results(self):
        """The group's UZF results accessor (source of the diffed fields)."""

        return self.group.packages.uzf.results

    @property
    def gwrch(self) -> CellBudgetResultDiff:
        """Groundwater-recharge difference vs the reference model."""

        return CellBudgetResultDiff(self._diff, self._uzf_results().gwrch)

    @property
    def sat(self) -> CellBudgetResultDiff:
        """Cell saturation difference vs the reference model.

        MF6's ``DATA-SAT`` record, which NPF writes under ``SAVE_SATURATION`` for
        EVERY groundwater cell (the saturated fraction of the cell, 0-1) -- not
        the unsaturated zone's saturation, though it sits under ``uzf``. Its
        ``get()``, ``summary()`` and ``map()`` difference that fraction (they
        raised until ledger 471, which left one ``sat`` column).
        """

        return CellBudgetResultDiff(self._diff, self._uzf_results().sat)


class StageResultDiff(_ResultDiffBase):
    """Stage-difference results (LAK lakes or SFR reaches) vs the reference.

    Wraps a grouped stage accessor keyed by a feature id (``lake`` or ``reach``)
    and adds per-model Δstage stats, the worst feature/period, and a
    within-tolerance flag.
    """

    def __init__(self, diff, accessor, *, entity: str):
        """Wrap a grouped stage ``accessor`` keyed by ``entity`` (``"lake"`` or ``"reach"``)."""

        super().__init__(diff)
        self._accessor = accessor
        self._entity = entity  # "lake" | "reach"

    def get(self, *, model_name=None, per=None, entity=None) -> pd.DataFrame:
        """Return aligned stage differences (per / feature / stage_diff).

        The group's stage compare: each member's stage file, masked (a dry or
        inactive lake or reach is ``NaN``, ledgers 329 and 353), reduced to its
        own last saved step of each stress period and inner-joined to the
        reference's on ``(per, <entity>)``. Also reachable as ``timeseries``. A
        member carrying several LAK (SFR) packages raises
        ``mf.SeveralPackagesError``: this diff cannot name one.

        Parameters
        ----------
        model_name : str or None, optional
            Keep one compared model, by its name in the group. Default ``None``
            keeps every non-reference model. A name not in the group raises
            ``KeyError`` and the reference's ``ValueError`` -- only when there are
            stage rows; an empty comparison comes back unchecked.
        per : int or None, optional
            Keep one zero-based stress period, each member at its own last saved
            step of it. A period no member saved gives an empty frame. Default
            ``None`` keeps every period.
        entity : int or None, optional
            Keep one feature, zero-based: the lake number on
            ``diff.packages.lak.results.stage``, the reach number on
            ``diff.packages.sfr.results.stage``. Default ``None`` keeps every one.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``reference_model``, ``per``, ``lake`` or ``reach``,
            ``stage``, ``reference_stage`` (model length units, ``NaN`` where a
            side has none) and ``stage_diff`` (model minus reference).
        """

        kwargs = {"per": per}
        if entity is not None:
            kwargs[self._entity] = entity
        frame = self._accessor.compare(**kwargs)
        if model_name is not None and not frame.empty:
            frame = frame[frame["model"].isin(self._targets(model_name))]
        return frame.reset_index(drop=True)

    # A stage series over stress periods; alias reads naturally for timeseries use.
    timeseries = get

    def summary(
        self,
        *,
        model_name=None,
        atol: float = _DEFAULT_ATOL,
        rtol: float = _DEFAULT_RTOL,
        per=None,
    ) -> pd.DataFrame:
        """One row per model: max/mean/RMSE Δstage, the worst feature/period, within tolerance.

        A lake or reach with no stage -- dry or inactive, ``NaN`` in the stage
        table (ledgers 329, 353) -- has no Δstage to average. The statistics run
        over the rows where BOTH models have a stage, and ``n`` counts those.
        ``n_one_sided`` counts the rows where exactly one of them does: one model's
        lake or reach is dry where the other's holds water, a difference no number
        measures, so any such row makes ``within_tolerance`` False. Rows where
        neither has a stage agree and count in neither column.

        Parameters
        ----------
        model_name : str, optional
            One compared model, by its name in the group; ``None`` (default)
            summarizes every non-reference member.
        atol : float, default 1e-3
            Absolute tolerance, in the model's length unit: a row is within
            tolerance when ``|Δstage| <= atol + rtol * |reference stage|``. Only
            rows where both models have a stage are tested.
        rtol : float, default 1e-3
            Relative tolerance, dimensionless, multiplying ``|reference stage|``
            in the same test as ``atol``.
        per : int, optional
            Restrict to one zero-based stress period, each model read at its own
            last saved step of it (the grouped stage table's ``per``; until ledger
            348 the index of the saved stage record for LAK, and a REACH number for
            SFR). ``None`` (default) uses every period.

        Returns
        -------
        pandas.DataFrame
            One row per compared model: ``model``, ``n``, ``n_one_sided``,
            ``max_abs_diff``, ``mean_abs_diff``, ``rmse``, ``argmax_<entity>``
            (``lake`` or ``reach``), ``argmax_per`` and ``within_tolerance``. With
            no row where both have a stage, the three statistics are ``NaN`` and
            the two argmax columns ``-1``. Empty when there is nothing to compare.

        Raises
        ------
        KeyError
            If ``model_name`` names no model in the group. Only when the group
            has stage rows to compare; with none, the result is the empty table
            whatever ``model_name`` says.
        ValueError
            If ``model_name`` is the group's reference model, which cannot be
            diffed against itself. Only when there are stage rows, as above.

        Examples
        --------
        >>> group = mf.ModelGroup({"base": base, "pumped": pumped}, reference="base")
        >>> stage = group.diff().packages.lak.results.stage
        >>> stage.summary()                           # one row: "pumped" vs "base"
        >>> stage.summary(model_name="pumped", atol=0.01, per=3)
        """

        data = self.get(model_name=model_name, per=per)
        argmax_column = f"argmax_{self._entity}"
        columns = [
            "model", "n", "n_one_sided", "max_abs_diff", "mean_abs_diff", "rmse",
            argmax_column, "argmax_per", "within_tolerance",
        ]
        if data.empty or "stage_diff" not in data.columns:
            return pd.DataFrame(columns=columns)
        rows = []
        for name, sub in data.groupby("model"):
            has_stage = np.isfinite(sub["stage"].to_numpy(dtype=float))
            has_reference = np.isfinite(sub["reference_stage"].to_numpy(dtype=float))
            both = has_stage & has_reference
            one_sided = int((has_stage ^ has_reference).sum())
            compared = sub.loc[both]
            diff = compared["stage_diff"].to_numpy(dtype=float)
            reference = compared["reference_stage"].to_numpy(dtype=float)
            abs_diff = np.abs(diff)
            imax = int(np.argmax(abs_diff)) if abs_diff.size else 0
            rows.append(
                {
                    "model": name,
                    "n": int(abs_diff.size),
                    "n_one_sided": one_sided,
                    "max_abs_diff": float(abs_diff.max()) if abs_diff.size else float("nan"),
                    "mean_abs_diff": float(abs_diff.mean()) if abs_diff.size else float("nan"),
                    "rmse": float(np.sqrt(np.mean(diff**2))) if diff.size else float("nan"),
                    argmax_column: int(compared[self._entity].iloc[imax]) if abs_diff.size else -1,
                    "argmax_per": int(compared["per"].iloc[imax]) if abs_diff.size else -1,
                    "within_tolerance": bool(
                        one_sided == 0
                        and np.all(_within_tolerance(diff, reference, atol, rtol))
                    ),
                }
            )
        return pd.DataFrame(rows, columns=columns)


def _declares_mover(model, package_type: str) -> bool:
    """Whether ``model`` carries a package of TYPE ``package_type`` declaring ``MOVER``."""

    for package in model_packages_of_type(model, package_type):
        option = getattr(package, "mover", None)
        value = option.get_data() if hasattr(option, "get_data") else option
        if value:
            return True
    return False


class MvrResultDiff(_ResultDiffBase):
    """Mover (MVR) difference vs the reference model, read where MF6 writes it.

    One package type and direction at a time: the water a package hands the mover
    (``"to_mvr"``) or takes from it (``"from_mvr"``), each member at its own last
    saved step of every period, differenced against the reference's. MF6 writes the
    two in different files by package (ledger 486):

    - a DRN, GHB, RIV or WEL only PROVIDES: ``<TYPE>-TO-MVR`` in the groundwater
      budget file, one row per boundary cell;
    - SFR and LAK provide and receive: ``TO-MVR`` and ``FROM-MVR`` in the package's
      OWN budget file, one row per reach or lake;
    - UZF writes ``UZF-GWD TO-MVR`` in the groundwater budget and ``FROM-MVR`` /
      ``REJ-INF-TO-MVR`` in its own, which is not read here: ``package="uzf"`` is
      refused, and :meth:`summary` warns when a member's UZF declares ``MOVER``.

    Until ledger 486 it read records named ``FROM-MVR``/``TO-MVR`` from the
    groundwater budget file, where MF6 writes neither, so every table was empty
    although water moved.
    """

    # Packages MF6 accepts as mover providers/receivers, i.e. exactly those
    # whose FloPy constructor exposes a ``mover`` option (CHD, EVT and RCH have
    # none). Derived from the registry as of 4.7.3; ``rch`` was wrongly listed
    # here for years (ledger entry 14) and deriving it makes that class of
    # error unrepresentable.
    _MOVER_PACKAGES = tuple(
        name
        for name, spec in _PACKAGE_EXPLORER_SPECS.items()
        if spec.capabilities.mover
    )
    #: the packages that write their mover terms in their own budget file, per feature
    _FEATURE_TABLES = {"sfr": ("reach", "sfr"), "lak": ("lake", "lak")}
    _DIRECTIONS = (("from_mvr", "FROM-MVR"), ("to_mvr", "TO-MVR"))
    _COLUMNS = [
        "model", "package", "direction", "value", "n",
        "max_abs_diff", "mean_abs_diff", "rmse", "argmax_cell", "argmax_feature",
        "within_tolerance",
    ]

    def _checked(self, package: str, direction: str) -> tuple[str, str]:
        """``(package type, direction)`` lower-cased, refused by name when MF6 writes no such term."""

        name = str(package).lower()
        if name not in self._MOVER_PACKAGES:
            raise ValueError(
                f"{package!r} is no mover package type; MF6 moves water for "
                f"{', '.join(self._MOVER_PACKAGES)}."
            )
        if name == "uzf":
            raise ValueError(
                "UZF's mover terms are not read by this diff: MF6 writes 'UZF-GWD TO-MVR' "
                "in the groundwater budget file and 'FROM-MVR' / 'REJ-INF-TO-MVR' in "
                "UZF's own budget file, which myflopy has no reader for."
            )
        directions = dict(self._DIRECTIONS)
        if direction not in directions:
            raise ValueError(
                f"direction={direction!r}; spell it 'to_mvr' (water the package hands "
                "the mover) or 'from_mvr' (water it takes from the mover)."
            )
        if direction == "from_mvr" and name not in self._FEATURE_TABLES:
            raise ValueError(
                f"A {name.upper()} only provides water to the mover: MF6 writes "
                f"'{name.upper()}-TO-MVR' for it and no FROM-MVR. Use direction='to_mvr'; "
                "only SFR and LAK (and UZF) receive."
            )
        return name, direction

    def _feature_compare(self, package: str, term: str, *, per=None) -> pd.DataFrame:
        """One SFR/LAK mover term from each member's package budget file, paired by period and feature."""

        builder = {"sfr": build_sfr_budget_term_table, "lak": build_lak_budget_term_table}[package]
        frames = []
        for name, model in self.group.models.items():
            if not model_packages_of_type(model, package):
                continue
            frame = builder(model, term=term)
            if frame.empty:
                continue
            frame = frame.copy()
            frame["model"] = name
            frames.append(frame)
        if not frames:
            return pd.DataFrame()
        data = _reduce_to_period_end(pd.concat(frames, ignore_index=True))
        if per is not None:
            data = data[data["per"] == int(per)]
        key_columns = _stable_compare_keys(data, ["q"])
        reference = self.group.reference
        ref = data[data["model"] == reference][[*key_columns, "q"]].rename(
            columns={"q": "reference_q"}
        )
        comp = data[data["model"] != reference].merge(ref, on=key_columns, how="inner")
        comp["reference_model"] = reference
        comp["q_diff"] = comp["q"].astype(float) - comp["reference_q"].astype(float)
        ordered = ["model", "reference_model", *key_columns, "q", "reference_q", "q_diff"]
        return comp[ordered].reset_index(drop=True)

    def get(
        self,
        *,
        package: str,
        direction: str = "to_mvr",
        model_name=None,
        per=None,
        layer=None,
        cells=None,
    ) -> pd.DataFrame:
        """Return per-row mover-flow differences for one package type and direction.

        Each member's mover term, reduced to its own last saved step of every
        stress period (``period_end.at_last_saved_steps``, ledger 348) and paired
        with the reference's on its integer and text identity (the period, and the
        cell and package NAME for a boundary; the reach or lake for SFR and LAK)
        by an inner join. ``q`` is MF6's own number in model volume per time:
        NEGATIVE for water a package hands the mover, POSITIVE for water it takes.
        A member that does not carry the package type contributes no rows.

        Parameters
        ----------
        package : str
            The package TYPE, any case: ``"drn"``, ``"ghb"``, ``"riv"``, ``"wel"``,
            ``"sfr"`` or ``"lak"`` -- every package of the type is read, each row
            naming its own (``pname``) for a boundary; SFR and LAK read the
            model's ONE package of the type and raise ``SeveralPackagesError``
            with several. ``"uzf"`` raises ``ValueError`` (its own budget file is
            not read), as does any other name.
        direction : str, default "to_mvr"
            ``"to_mvr"``, the water the package hands the mover (DRN/GHB/RIV/WEL:
            ``<TYPE>-TO-MVR`` in the groundwater budget file; SFR/LAK: ``TO-MVR``
            in the package's budget file), or ``"from_mvr"``, the water it takes
            (SFR/LAK ``FROM-MVR``). A boundary package with ``"from_mvr"`` raises
            ``ValueError``: MF6 writes no such term, since only SFR, LAK, UZF
            and MAW receive. Any other spelling raises ``ValueError``. The default
            was ``"from_mvr"`` until ledger 486.
        model_name : str or None, optional
            Keep one compared model, by its name in the group. Default ``None``
            keeps every non-reference model. A name not in the group raises
            ``KeyError`` and the reference's ``ValueError`` -- only when there are
            rows to filter.
        per : int or None, optional
            Keep one zero-based stress period, each member at its own last saved
            step of it. Default ``None`` keeps every period.
        layer : int or sequence of int or None, optional
            Keep these zero-based layers of a boundary package's cells. Default
            ``None`` keeps every layer. Given for SFR or LAK it raises
            ``ValueError``: their rows are reaches and lakes, not cells.
        cells : int or sequence of int or None, optional
            Keep these zero-based in-layer cells of a boundary package. Default
            ``None`` keeps every cell; for SFR or LAK it raises ``ValueError`` as
            ``layer`` does.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``reference_model``, the identity columns (a boundary's
            ``node``, ``node2``, ``layer``, ``cell``, ``package``, ``pname``,
            ``per``; SFR's ``package``, ``term``, ``per``, ``reach``, ``layer``,
            ``cell``, ``node``, ``node2``; LAK's the same with ``lake`` and no
            cell), ``q``, ``reference_q`` and ``q_diff`` (model minus reference).
            Empty, with no columns, when no member wrote the term.
        """

        package, direction = self._checked(package, direction)
        term = dict(self._DIRECTIONS)[direction]
        if package in self._FEATURE_TABLES:
            if layer is not None or cells is not None:
                raise ValueError(
                    f"A {package.upper()} mover row is a "
                    f"{self._FEATURE_TABLES[package][0]}, not a cell: layer= and cells= "
                    "select a boundary package's cells only."
                )
            frame = self._feature_compare(package, term, per=per)
        else:
            accessor = GroupCellPackageResults(
                self.group, package, budget_text=f"{package.upper()}-TO-MVR", value_name="q"
            )
            frame = CellBudgetResultDiff(self._diff, accessor).get(
                per=per, layer=layer, cells=cells
            )
        if model_name is not None and not frame.empty:
            frame = frame[frame["model"].isin(self._targets(model_name))]
        return frame.reset_index(drop=True)

    def summary(
        self, *, model_name=None, atol: float = _DEFAULT_ATOL, rtol: float = _DEFAULT_RTOL
    ) -> pd.DataFrame:
        """Per (model, package, direction) mover-flow difference statistics and tolerance.

        Reads :meth:`get` for every package type some member declares ``MOVER``
        on -- DRN, GHB, RIV and WEL in the direction they have (``"to_mvr"``),
        SFR and LAK in both -- and leaves out a pair no member wrote a row for. A
        member's UZF declaring ``MOVER`` is NOT read (:meth:`get`): it is warned
        about, by ``warnings.warn`` and the log, and left out.

        Parameters
        ----------
        model_name : str or None, optional
            Summarize one compared model, by its name in the group. Default
            ``None`` summarizes every non-reference model. A name not in the
            group raises ``KeyError`` and the reference's ``ValueError`` -- when
            some term has rows.
        atol : float, optional
            Absolute tolerance, in model volume per time, for
            ``|q_diff| <= atol + rtol * |reference q|``, row by row. Default
            ``1e-3``.
        rtol : float, optional
            Relative tolerance, dimensionless, in the same test. Default ``1e-3``.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``package``, ``direction``, ``value`` (``"q"``), ``n``,
            ``max_abs_diff``, ``mean_abs_diff``, ``rmse``, ``argmax_cell`` (the
            in-layer cell of the largest ``|q_diff|`` for a boundary package, else
            ``-1``), ``argmax_feature`` (its reach or lake for SFR/LAK, else
            ``-1``) and ``within_tolerance``; one row per model, package and
            direction with rows. Empty when no member moves water.
        """

        models = list(self.group.models.values())
        if any(_declares_mover(model, "uzf") for model in models):
            message = (
                "diff.packages.mvr.results.summary() leaves out UZF's mover terms: MF6 "
                "writes them in UZF's own budget file, which this diff does not read."
            )
            warnings.warn(message, UserWarning, stacklevel=2)
            logger.warning(message)
        rows = []
        for package in self._MOVER_PACKAGES:
            if package == "uzf" or not any(_declares_mover(model, package) for model in models):
                continue
            directions = (
                self._DIRECTIONS if package in self._FEATURE_TABLES else (("to_mvr", "TO-MVR"),)
            )
            feature = self._FEATURE_TABLES.get(package, (None,))[0]
            for direction, _term in directions:
                data = self.get(package=package, direction=direction, model_name=model_name)
                if data.empty:
                    continue
                for name, sub in data.groupby("model"):
                    diff = sub["q_diff"].to_numpy(dtype=float)
                    reference = sub["reference_q"].to_numpy(dtype=float)
                    abs_diff = np.abs(diff)
                    imax = int(np.argmax(abs_diff))
                    rows.append(
                        {
                            "model": name,
                            "package": package,
                            "direction": direction,
                            "value": "q",
                            "n": int(abs_diff.size),
                            "max_abs_diff": float(abs_diff.max()),
                            "mean_abs_diff": float(abs_diff.mean()),
                            "rmse": float(np.sqrt(np.mean(diff**2))),
                            "argmax_cell": (
                                int(sub["cell"].iloc[imax]) if feature is None else -1
                            ),
                            "argmax_feature": (
                                -1 if feature is None else int(sub[feature].iloc[imax])
                            ),
                            "within_tolerance": bool(
                                np.all(_within_tolerance(diff, reference, atol, rtol))
                            ),
                        }
                    )
        return pd.DataFrame(rows, columns=self._COLUMNS)


class LakResultsDiffNamespace(CellResultsDiffNamespace):
    """``diff.packages.lak.results`` -- fields ``q`` (default) and ``stage``."""

    def __init__(self, diff):
        """Bind the LAK results-diff namespace (fields ``q`` + ``stage``) to the group."""

        super().__init__(diff, "lak")

    def _field_names(self):
        """The mappable LAK results-diff fields: exchange ``q`` and ``stage``."""

        return ["q", "stage"]

    @property
    def stage(self) -> StageResultDiff:
        """Lake-stage difference vs the reference model (per lake / period)."""

        return StageResultDiff(self._diff, GroupLakStageResults(self.group), entity="lake")


class SfrResultsDiffNamespace(CellResultsDiffNamespace):
    """``diff.packages.sfr.results`` -- fields ``q`` (default) and ``stage``."""

    def __init__(self, diff):
        """Bind the SFR results-diff namespace (fields ``q`` + ``stage``) to the group."""

        super().__init__(diff, "sfr")

    def _field_names(self):
        """The mappable SFR results-diff fields: exchange ``q`` and ``stage``."""

        return ["q", "stage"]

    @property
    def stage(self) -> StageResultDiff:
        """Reach-stage difference vs the reference model (per reach / period)."""

        return StageResultDiff(self._diff, GroupSfrStageResults(self.group), entity="reach")
