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

import numpy as np
import pandas as pd

from myflopy.modflow.mf6.package_explorer import DiffSpatialView, FieldMappable
from myflopy.modflow.utils.datatypes.xsections import XSection
from myflopy.project.group.lak import GroupLakStageResults
from myflopy.project.group.results import GroupCellPackageResults
from myflopy.project.group.sfr import GroupSfrStageResults

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


class HeadsResultDiff(_ResultDiffBase, DiffSpatialView):
    """Head-difference results vs the reference model (per cell / layer / time).

    ``map`` / ``mosaic`` / ``animate`` render Δhead (model - reference) --
    ``diff.hds.map("F9b", per=8, layer=1)``, ``diff.hds.mosaic(by="model")``.
    """

    def get(self, *, model_name=None, per=None, layer=None, cells=None) -> pd.DataFrame:
        """Return aligned head differences: ``elev`` / ``reference_elev`` / ``diff``."""

        frame = self.group.hds.compare(per=per, layer=layer, cells=cells)
        if model_name is not None:
            targets = self._targets(model_name)
            frame = frame[frame["model"].isin(targets)]
        return frame.reset_index(drop=True)

    # -- spatial-view hooks (delta head maps) ---------------------------------
    def _spatial_map(self, *, per=0, layer=0, model=None, **kwargs):
        """Δhead choropleth (compared model - reference) for one model."""

        return self.group.hds.compare_map(model_name=model, per=per, layer=layer, **kwargs)

    def _spatial_models(self):
        """The non-reference model names -- one Δhead panel each."""

        return self._targets()

    def _spatial_reference_model(self):
        """The reference model whose grid/layers frame the Δhead maps."""

        return self.group.models[self.group.reference]

    def _spatial_value_label(self):
        """The mapped quantity's label -- head difference."""

        return "head"

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
        whether every cell is within tolerance."""

        data = self.get(model_name=model_name, per=per, layer=layer)
        columns = [
            "model", "n", "max_abs_diff", "mean_abs_diff", "rmse",
            "argmax_cell", "argmax_layer", "argmax_kstpkper", "within_tolerance",
        ]
        if data.empty:
            return pd.DataFrame(columns=columns)
        rows = []
        for name, sub in data.groupby("model"):
            diff = sub["diff"].to_numpy(dtype=float)
            reference = sub["reference_elev"].to_numpy(dtype=float)
            abs_diff = np.abs(diff)
            imax = int(np.argmax(abs_diff)) if abs_diff.size else 0
            rows.append(
                {
                    "model": name,
                    "n": int(abs_diff.size),
                    "max_abs_diff": float(abs_diff.max()) if abs_diff.size else 0.0,
                    "mean_abs_diff": float(abs_diff.mean()) if abs_diff.size else 0.0,
                    "rmse": float(np.sqrt(np.mean(diff**2))) if diff.size else 0.0,
                    "argmax_cell": int(sub["cell"].iloc[imax]) if abs_diff.size else -1,
                    "argmax_layer": int(sub["layer"].iloc[imax]) if abs_diff.size else -1,
                    "argmax_kstpkper": sub["kstpkper"].iloc[imax] if abs_diff.size else None,
                    "within_tolerance": bool(
                        np.all(_within_tolerance(diff, reference, atol, rtol))
                    ),
                }
            )
        return pd.DataFrame(rows, columns=columns)


class BudgetResultDiff(_ResultDiffBase):
    """Overall (volumetric listing) budget difference vs the reference model.

    Compares the incremental per-timestep budget term by term, aligned on
    simulation time (``totim``) so models with the same physics but different
    time discretization still line up where their timesteps coincide.
    """

    def _model_budget(self, model_name: str) -> pd.DataFrame:
        """One model's incremental listing budget, melted to ``totim`` / ``term`` / ``value``."""

        # ``budget_incremental`` is a property on real models (returns a
        # DataFrame); tolerate a callable too so fakes/other readers work.
        raw = self.group.models[model_name].budget_incremental
        raw = raw() if callable(raw) else raw
        frame = pd.DataFrame(raw).copy()
        if "totim" not in frame.columns:
            frame = frame.reset_index().rename(columns={"index": "totim"})
        terms = [
            column
            for column in frame.columns
            if column != "totim" and pd.api.types.is_numeric_dtype(frame[column])
        ]
        return frame.melt(
            id_vars=["totim"], value_vars=terms, var_name="term", value_name="value"
        )

    def get(self, *, model_name=None) -> pd.DataFrame:
        """Return aligned per-term budget differences per timestep."""

        reference = self._model_budget(self.group.reference)
        columns = ["model", "totim", "term", "reference_value", "model_value", "diff"]
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
        per-timestep difference, and whether every timestep is within tolerance."""

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
        """Return aligned per-cell budget differences (value / reference / diff)."""

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
        is, and whether every cell is within tolerance."""

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
_CELL_BUDGET_PACKAGES = ("ghb", "drn", "chd", "riv", "wel", "rch", "evt", "sfr", "lak")


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
        """Unsaturated-zone saturation difference vs the reference model."""

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
        """Return aligned stage differences (per / feature / stage_diff)."""

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
        """One row per model: max/mean/RMSE Δstage, the worst feature/period, within tolerance."""

        data = self.get(model_name=model_name, per=per)
        argmax_column = f"argmax_{self._entity}"
        columns = [
            "model", "n", "max_abs_diff", "mean_abs_diff", "rmse",
            argmax_column, "argmax_per", "within_tolerance",
        ]
        if data.empty or "stage_diff" not in data.columns:
            return pd.DataFrame(columns=columns)
        rows = []
        for name, sub in data.groupby("model"):
            diff = sub["stage_diff"].to_numpy(dtype=float)
            reference = sub["reference_stage"].to_numpy(dtype=float)
            abs_diff = np.abs(diff)
            imax = int(np.argmax(abs_diff)) if abs_diff.size else 0
            rows.append(
                {
                    "model": name,
                    "n": int(abs_diff.size),
                    "max_abs_diff": float(abs_diff.max()) if abs_diff.size else 0.0,
                    "mean_abs_diff": float(abs_diff.mean()) if abs_diff.size else 0.0,
                    "rmse": float(np.sqrt(np.mean(diff**2))) if diff.size else 0.0,
                    argmax_column: int(sub[self._entity].iloc[imax]) if abs_diff.size else -1,
                    "argmax_per": int(sub["per"].iloc[imax]) if abs_diff.size else -1,
                    "within_tolerance": bool(
                        np.all(_within_tolerance(diff, reference, atol, rtol))
                    ),
                }
            )
        return pd.DataFrame(rows, columns=columns)


class MvrResultDiff(_ResultDiffBase):
    """Mover (MVR) difference vs the reference model.

    MF6 realizes the mover as ``FROM-MVR`` / ``TO-MVR`` budget terms in the
    packages it moves water between (no standalone MVR output), so this compares
    those terms per package and direction across models.
    """

    # Packages MF6 accepts as mover providers/receivers (FloPy exposes a
    # ``mover`` option on each). EVT is absent because it has none; "rch" is a
    # pre-existing entry that also has none -- harmless, the term lookup below
    # is try-wrapped and simply finds nothing.
    _MOVER_PACKAGES = ("lak", "sfr", "uzf", "drn", "ghb", "riv", "wel", "rch")
    _DIRECTIONS = (("from_mvr", "FROM-MVR"), ("to_mvr", "TO-MVR"))
    _COLUMNS = [
        "model", "package", "direction", "value", "n",
        "max_abs_diff", "mean_abs_diff", "rmse", "argmax_cell", "within_tolerance",
    ]

    def _cell_diff(self, package: str, term: str) -> CellBudgetResultDiff:
        """Build a cell-budget diff for one package's mover ``term`` (FROM/TO-MVR)."""

        accessor = GroupCellPackageResults(
            self.group, package, budget_text=term, value_name="q"
        )
        return CellBudgetResultDiff(self._diff, accessor)

    def get(self, *, package: str, direction: str = "from_mvr", model_name=None, **kwargs):
        """Return per-cell mover flow differences for one package + direction."""

        term = dict(self._DIRECTIONS)[direction]
        return self._cell_diff(package.lower(), term).get(model_name=model_name, **kwargs)

    def summary(
        self, *, model_name=None, atol: float = _DEFAULT_ATOL, rtol: float = _DEFAULT_RTOL
    ) -> pd.DataFrame:
        """Per (model, package, direction) mover-flow Δ stats + within tolerance.

        Packages without a mover term (or without results) are skipped.
        """

        frames = []
        for package in self._MOVER_PACKAGES:
            for direction, term in self._DIRECTIONS:
                try:
                    part = self._cell_diff(package, term).summary(
                        model_name=model_name, atol=atol, rtol=rtol
                    )
                except Exception:
                    continue
                if part.empty:
                    continue
                part = part.copy()
                part["direction"] = direction
                frames.append(part)
        if not frames:
            return pd.DataFrame(columns=self._COLUMNS)
        return pd.concat(frames, ignore_index=True)[self._COLUMNS]


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
