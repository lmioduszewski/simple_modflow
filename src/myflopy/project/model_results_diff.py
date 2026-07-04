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

Phase 5a covers **heads** and the overall (volumetric) **budget**; later
sub-phases add per-package cell budgets, UZF, lake/SFR stage & flow, and MVR.
All stress-period/time references are zero-based.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

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
        self._diff = diff
        self.group = diff.group

    def _targets(self, model_name=None) -> list[str]:
        if model_name is not None:
            name = str(model_name)
            if name not in self.group.models:
                raise KeyError(f"Model {name!r} is not in the group.")
            if name == self.group.reference:
                raise ValueError("The reference model cannot be diffed against itself.")
            return [name]
        return [name for name in self.group.models if name != self.group.reference]


class HeadsResultDiff(_ResultDiffBase):
    """Head-difference results vs the reference model (per cell / layer / time)."""

    def get(self, *, model_name=None, per=None, layer=None, cells=None) -> pd.DataFrame:
        """Return aligned head differences: ``elev`` / ``reference_elev`` / ``diff``."""

        frame = self.group.hds.compare(per=per, layer=layer, cells=cells)
        if model_name is not None:
            targets = self._targets(model_name)
            frame = frame[frame["model"].isin(targets)]
        return frame.reset_index(drop=True)

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


class CellBudgetResultDiff(_ResultDiffBase):
    """Per-package cell-budget difference vs the reference model.

    Wraps a grouped cell-budget accessor (GHB/DRN/… leakage, SFR/LAK exchange,
    UZF recharge/saturation) and adds Δ summary stats + a within-tolerance flag on
    its value column (``q``, ``gwrch``, ``sat``, …).
    """

    def __init__(self, diff, accessor):
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


# Packages with a per-cell budget diff via group.packages.<pkg>.results.q.
_CELL_BUDGET_PACKAGES = ("ghb", "drn", "chd", "wel", "rch", "sfr", "lak")


class _ResultsPackageNamespace:
    """``diff.results.packages.<pkg>`` -> :class:`CellBudgetResultDiff`."""

    def __init__(self, diff):
        self._diff = diff
        self.group = diff.group

    def __getattr__(self, name: str) -> CellBudgetResultDiff:
        pkg = str(name).lower()
        try:
            accessor = getattr(self.group.packages, pkg)
        except AttributeError as exc:
            raise AttributeError(
                f"No grouped results accessor for package {name!r}."
            ) from exc
        cell = getattr(getattr(accessor, "results", None), "q", None)
        if cell is None:
            raise AttributeError(f"Package {name!r} has no cell-budget results.")
        return CellBudgetResultDiff(self._diff, cell)

    def __dir__(self):
        return sorted(set(super().__dir__()) | set(_CELL_BUDGET_PACKAGES))


class _ResultsUzfNamespace:
    """``diff.results.uzf.gwrch`` / ``.sat`` -> :class:`CellBudgetResultDiff`."""

    def __init__(self, diff):
        self._diff = diff
        self.group = diff.group

    def _uzf_results(self):
        return self.group.packages.uzf.results

    @property
    def gwrch(self) -> CellBudgetResultDiff:
        return CellBudgetResultDiff(self._diff, self._uzf_results().gwrch)

    @property
    def sat(self) -> CellBudgetResultDiff:
        return CellBudgetResultDiff(self._diff, self._uzf_results().sat)


class ResultsDiff:
    """Results-tier facade for :class:`~myflopy.project.model_diff.ModelDiff`.

    Compares computed outputs between models -- ``heads``, the overall ``budget``,
    per-package cell budgets (``packages.<pkg>``), and ``uzf`` recharge/saturation.
    Requires completed runs; accessors raise if outputs are missing.
    """

    def __init__(self, diff):
        self._diff = diff
        self.group = diff.group

    @property
    def heads(self) -> HeadsResultDiff:
        return HeadsResultDiff(self._diff)

    @property
    def budget(self) -> BudgetResultDiff:
        return BudgetResultDiff(self._diff)

    @property
    def packages(self) -> _ResultsPackageNamespace:
        return _ResultsPackageNamespace(self._diff)

    @property
    def uzf(self) -> _ResultsUzfNamespace:
        return _ResultsUzfNamespace(self._diff)
