"""Reopen and evaluate completed ``simple_modflow`` PEST workspaces."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.geometry import Polygon

from simple_modflow.modflow.mf6.observations import HeadTargets
from simple_modflow.modflow.mf6.pest.forward_run import _apply_drain_specs, _apply_k_specs
from simple_modflow.project.run_model import LoadedMf6Run, load_mf6_run

METADATA_FILENAME = "simple_modflow_pest_metadata.json"
MATERIALIZATION_FILENAME = "simple_modflow_calibrated_materialization.json"
MATERIALIZATION_VERSION = 2


@dataclass
class PestRunReview:
    """Convenience bundle for a reopened calibration review."""

    run: "PestRunResults"
    targets: HeadTargets
    baseline_model: LoadedMf6Run
    calibrated_model: LoadedMf6Run
    residual_compare: pd.DataFrame
    stats: pd.DataFrame
    k_geodata: gpd.GeoDataFrame


def _is_model_workspace(path: Path) -> bool:
    """Return ``True`` when a directory looks like an MF6 model workspace."""

    return path.is_dir() and any(item.name.lower() != "mfsim.nam" for item in path.glob("*.nam"))


def _is_pest_workspace(path: Path) -> bool:
    """Return ``True`` when a directory looks like a generated PEST workspace."""

    return path.is_dir() and (path / "pest_forward_config.json").exists() and any(path.glob("*.pst"))


def _discover_pest_workspace(workspace: Path) -> Path:
    """Discover the actual PEST workspace from a root or direct path."""

    workspace = Path(workspace)
    if _is_pest_workspace(workspace):
        return workspace

    # Prefer completed master runs when both the template `pest/` workspace and
    # a finished `pest_master/` workspace are present under one artifact root.
    candidate_names = ("pest_master", "pest")
    checked: list[Path] = [workspace]
    for name in candidate_names:
        nested = workspace / name
        checked.append(nested)
        if _is_pest_workspace(nested):
            return nested

    raise FileNotFoundError(
        "Could not locate a PEST workspace. Provide either a generated "
        "`pest/` or `pest_master/` directory, or its artifact root. "
        f"Checked: {', '.join(str(path) for path in checked)}"
    )


def _discover_baseline_workspace(root: Path, pest_workspace: Path) -> Path:
    """Discover the paired baseline MF6 workspace for a PEST run."""

    root = Path(root)
    candidate_roots = []
    model_dir = root / "model"
    if model_dir.is_dir():
        candidate_roots.append(model_dir)
    sibling_model = pest_workspace.parent / "model"
    if sibling_model.is_dir() and sibling_model not in candidate_roots:
        candidate_roots.append(sibling_model)
    if _is_model_workspace(root):
        candidate_roots.append(root)

    for candidate in candidate_roots:
        if _is_model_workspace(candidate):
            return candidate
        child_dirs = [child for child in candidate.iterdir() if child.is_dir()]
        if len(child_dirs) == 1:
            return child_dirs[0]
        child_workspaces = [child for child in child_dirs if _is_model_workspace(child)]
        if len(child_workspaces) == 1:
            return child_workspaces[0]
        if len(child_workspaces) > 1 and candidate.name.lower() == "model":
            return child_workspaces[0]

    recursive = []
    for candidate in candidate_roots:
        if not candidate.exists():
            continue
        for nam_path in candidate.rglob("*.nam"):
            if nam_path.name.lower() == "mfsim.nam":
                continue
            recursive.append(nam_path.parent)
    recursive = list(dict.fromkeys(recursive))
    if recursive:
        return recursive[0]

    raise FileNotFoundError(
        "Could not locate a baseline MF6 workspace for the provided PEST run. "
        f"Checked around: {root}"
    )


def _discover_final_parameter_file(pest_workspace: Path) -> Path | None:
    """Return the newest final-parameter file, if one exists."""

    candidates = sorted(pest_workspace.glob("*.par"), key=lambda path: path.stat().st_mtime, reverse=True)
    return candidates[0] if candidates else None


def _discover_metadata_file(pest_workspace: Path) -> Path | None:
    """Return the saved simple_modflow PEST metadata file when present."""

    candidate = pest_workspace / METADATA_FILENAME
    return candidate if candidate.exists() else None


def _discover_materialization_file(pest_workspace: Path) -> Path | None:
    """Return the calibrated-materialization marker when present."""

    candidate = pest_workspace / MATERIALIZATION_FILENAME
    return candidate if candidate.exists() else None


def _file_sha256(path: Path) -> str:
    """Return the SHA-256 hash of one file."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _clear_runtime_caches(model):
    """Clear cached result readers after rewriting or rerunning a workspace."""

    for attr in ("_hds", "_bud", "_times", "_obs", "_kstpkper"):
        if hasattr(model, attr):
            setattr(model, attr, None)


def _read_parameter_file(par_path: Path) -> pd.Series:
    """Read a PEST-style parameter file into a lower-case name/value series."""

    rows: list[tuple[str, float]] = []
    with open(par_path, "r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.lower() == "single point":
                continue
            parts = stripped.split()
            if len(parts) < 2:
                continue
            rows.append((parts[0].lower(), float(parts[1])))
    if not rows:
        raise ValueError(f"No parameter rows found in {par_path}")
    return pd.Series({name: value for name, value in rows}, name="value")


def _write_final_parameter_csvs(
    pest_workspace: Path,
    config: dict,
    final_parameters: pd.Series,
) -> dict:
    """Create temporary parameter CSVs populated with final PEST values."""

    config = json.loads(json.dumps(config))
    for spec_group in (config.get("k_specs", []), config.get("drain_specs", [])):
        for spec in spec_group:
            for key in ("points_meta_csv", "cells_meta_csv", "support_csv"):
                if key in spec:
                    spec[key] = str((pest_workspace / spec[key]).resolve())
            if "parameter_csv" not in spec:
                continue
            parameter_csv = pest_workspace / spec["parameter_csv"]
            frame = pd.read_csv(parameter_csv)
            frame["parnme"] = frame["parnme"].astype(str).str.lower()
            frame["value"] = frame["parnme"].map(final_parameters).fillna(frame["value"])
            final_csv = parameter_csv.with_name(f"{parameter_csv.stem}.final.csv")
            frame.to_csv(final_csv, index=False)
            spec["parameter_csv"] = str(final_csv.resolve())
    return config


@dataclass
class PestRunResults:
    """Reopen a PEST workspace and evaluate baseline vs calibrated models.

    Parameters
    ----------
    root
        Either the artifact root containing ``model/`` and ``pest/`` or
        ``pest_master/`` folders, or one of those generated workspaces directly.
    crs
        CRS passed through when reopening file-backed models.
    verbosity_level
        FloPy/simple_modflow load verbosity.
    """

    root: str | Path
    crs: str = "EPSG:2927"
    verbosity_level: int = 0

    def __post_init__(self):
        self.root = Path(self.root)
        self.pest_workspace = _discover_pest_workspace(self.root)
        self.baseline_workspace = _discover_baseline_workspace(self.root, self.pest_workspace)
        self.pst_path = sorted(self.pest_workspace.glob("*.pst"), key=lambda path: path.stat().st_mtime, reverse=True)[0]
        self.config_path = self.pest_workspace / "pest_forward_config.json"
        self.metadata_path = _discover_metadata_file(self.pest_workspace)
        self.materialization_path = _discover_materialization_file(self.pest_workspace)
        self.par_path = _discover_final_parameter_file(self.pest_workspace)
        self.rec_path = next(iter(sorted(self.pest_workspace.glob("*.rec"), key=lambda path: path.stat().st_mtime, reverse=True)), None)
        self._baseline_model: LoadedMf6Run | None = None
        self._calibrated_model: LoadedMf6Run | None = None
        self._final_parameters: pd.Series | None = None
        self._config_cache: dict | None = None
        self._metadata_cache: dict | None = None
        self._materialization_cache: dict | None = None

    @classmethod
    def from_workspace(cls, workspace: str | Path, *, crs: str = "EPSG:2927", verbosity_level: int = 0):
        """Open a PEST artifact root or direct ``pest`` workspace."""

        return cls(root=workspace, crs=crs, verbosity_level=verbosity_level)

    @property
    def has_final_parameters(self) -> bool:
        """Whether a final ``.par`` file is available for calibrated reopening."""

        return self.par_path is not None and self.par_path.exists()

    @property
    def final_parameters(self) -> pd.Series:
        """Return the final PEST parameter values."""

        if not self.has_final_parameters:
            raise FileNotFoundError(
                "No final .par file was found in the PEST workspace. "
                "This run likely did not finish successfully."
            )
        if self._final_parameters is None:
            self._final_parameters = _read_parameter_file(self.par_path)
        return self._final_parameters

    @property
    def forward_config(self) -> dict:
        """Return the saved forward-run configuration."""

        if self._config_cache is None:
            if not self.config_path.exists():
                raise FileNotFoundError(f"Missing pest_forward_config.json in {self.pest_workspace}")
            self._config_cache = json.loads(self.config_path.read_text(encoding="utf-8"))
        return self._config_cache

    @property
    def metadata(self) -> dict:
        """Return saved simple_modflow PEST metadata when available."""

        if self._metadata_cache is None:
            if self.metadata_path is None:
                self._metadata_cache = {}
            else:
                self._metadata_cache = json.loads(self.metadata_path.read_text(encoding="utf-8"))
        return self._metadata_cache

    @property
    def materialization_info(self) -> dict:
        """Return calibrated-output materialization metadata when present."""

        if self._materialization_cache is None:
            if self.materialization_path is None:
                self._materialization_cache = {}
            else:
                self._materialization_cache = json.loads(self.materialization_path.read_text(encoding="utf-8"))
        return self._materialization_cache

    def summary(self) -> pd.DataFrame:
        """Return a one-row summary of the discovered PEST result workspace."""

        return pd.DataFrame(
            [
                {
                    "root": str(self.root),
                    "baseline_workspace": str(self.baseline_workspace),
                    "pest_workspace": str(self.pest_workspace),
                    "pst": self.pst_path.name,
                    "par": self.par_path.name if self.par_path else None,
                    "rec": self.rec_path.name if self.rec_path else None,
                    "metadata": self.metadata_path.name if self.metadata_path else None,
                    "materialization": self.materialization_path.name if self.materialization_path else None,
                    "n_saved_head_target_sets": len(self.saved_observation_sets(kind="head_targets")),
                    "has_final_parameters": self.has_final_parameters,
                }
            ]
        )

    def _current_par_hash(self) -> str:
        """Return the hash of the final parameter file used for calibrated review."""

        if not self.has_final_parameters:
            raise FileNotFoundError("No final .par file is available for hashing.")
        return _file_sha256(self.par_path)

    def _materialization_is_current(self) -> bool:
        """Whether calibrated outputs already reflect the current final ``.par`` file."""

        if not self.has_final_parameters:
            return False
        info = self.materialization_info
        if not info:
            return False
        return (
            info.get("version") == MATERIALIZATION_VERSION
            and info.get("par_hash") == self._current_par_hash()
        )

    def saved_observation_sets(self, *, kind: str | None = None) -> list[dict]:
        """Return saved observation metadata entries from the workspace manifest."""

        observation_sets = list(self.metadata.get("observation_sets", []))
        if kind is None:
            return observation_sets
        kind_key = str(kind).strip().lower()
        return [
            item
            for item in observation_sets
            if str(item.get("kind", "")).strip().lower() == kind_key
        ]

    def load_head_targets(self, *, prefix: str | None = None) -> HeadTargets:
        """Load a saved head-target dataset from the workspace metadata."""

        entries = self.saved_observation_sets(kind="head_targets")
        if prefix is not None:
            prefix_key = str(prefix).strip().lower()
            entries = [
                item
                for item in entries
                if str(item.get("prefix", "")).strip().lower() == prefix_key
            ]
        if not entries:
            raise FileNotFoundError(
                "No saved head-target metadata was found in the PEST workspace."
            )
        if len(entries) > 1:
            raise ValueError(
                "Multiple saved head-target sets are available. Pass prefix=... to select one."
            )
        entry = entries[0]
        return HeadTargets(
            locations=self.pest_workspace / entry["locations_file"],
            values=self.pest_workspace / entry["values_file"],
            name_column=entry.get("name_column", "name"),
            layer_column=entry.get("layer_column", "layer"),
            group_column=entry.get("group_column", "group"),
            weight_column=entry.get("weight_column", "weight"),
            time_column=entry.get("time_column", "time"),
            value_column=entry.get("value_column", "head"),
        )

    def load_baseline_model(self) -> LoadedMf6Run:
        """Open the paired baseline model workspace as a file-backed model."""

        if self._baseline_model is None:
            self._baseline_model = load_mf6_run(
                self.baseline_workspace,
                crs=self.crs,
                verbosity_level=self.verbosity_level,
            )
        return self._baseline_model

    def load_calibrated_model(self) -> LoadedMf6Run:
        """Open the calibrated model state with final PEST parameters applied.

        The calibrated workspace outputs are materialized from the final
        ``.par`` file the first time this is called for a given parameter hash.
        This keeps reopened residual comparisons tied to the final parameter
        values rather than stale head files left behind by prior worker runs.
        """

        if self._calibrated_model is not None:
            return self._calibrated_model
        config = _write_final_parameter_csvs(
            self.pest_workspace,
            self.forward_config,
            self.final_parameters,
        )
        calibrated = load_mf6_run(
            self.pest_workspace,
            crs=self.crs,
            verbosity_level=self.verbosity_level,
        )
        calibrated.load_all()

        k_base = np.asarray(calibrated.gwf.npf.k.array, dtype=float).copy()
        for spec in config.get("k_specs", []):
            cells_meta = pd.read_csv(spec["cells_meta_csv"])
            for layer, group in cells_meta.groupby("layer", dropna=False):
                layer = int(layer)
                cell_idx = group["cell"].astype(int).to_numpy()
                base_vals = pd.to_numeric(group["base_k"], errors="coerce").to_numpy(dtype=float)
                k_base[layer, cell_idx] = base_vals
        calibrated.gwf.npf.k.set_data(k_base)

        _apply_k_specs(calibrated.gwf, config.get("k_specs", []))
        _apply_drain_specs(calibrated.gwf, config.get("drain_specs", []))

        if not self._materialization_is_current():
            _clear_runtime_caches(calibrated)
            calibrated.sim.write_simulation()
            success, buff = calibrated.sim.run_simulation(silent=False, report=True)
            if not success:
                tail = ""
                if isinstance(buff, (list, tuple)):
                    tail = "\n".join(str(item) for item in buff[-20:])
                else:
                    tail = str(buff)
                raise RuntimeError(
                    "Could not materialize calibrated MF6 outputs from the final .par file.\n"
                    f"Workspace: {self.pest_workspace}\n"
                    f"Run output tail:\n{tail}"
                )
            materialization = {
                "version": MATERIALIZATION_VERSION,
                "par_file": self.par_path.name if self.par_path else None,
                "par_hash": self._current_par_hash(),
            }
            materialization_path = self.pest_workspace / MATERIALIZATION_FILENAME
            materialization_path.write_text(json.dumps(materialization, indent=2), encoding="utf-8")
            self.materialization_path = materialization_path
            self._materialization_cache = materialization

        calibrated._pest_run_info = {
            "root": self.root,
            "baseline_workspace": self.baseline_workspace,
            "pest_workspace": self.pest_workspace,
            "pst": self.pst_path,
            "par": self.par_path,
            "rec": self.rec_path,
        }
        calibrated._pest_k_initial = np.asarray(k_base, dtype=float)
        _clear_runtime_caches(calibrated)
        self._calibrated_model = calibrated
        return calibrated

    def k_geodata(self) -> gpd.GeoDataFrame:
        """Return a GeoDataFrame of initial/final K values for the calibrated run."""

        calibrated = self.load_calibrated_model()
        k_initial = np.asarray(calibrated._pest_k_initial, dtype=float)[0, :]
        k_final = np.asarray(calibrated.gwf.npf.k.array, dtype=float)[0, :]
        k_ratio = k_final / k_initial

        ncpl = int(calibrated.gwf.modelgrid.ncpl)
        polygons = [Polygon(calibrated.gwf.modelgrid.get_cell_vertices(i)) for i in range(ncpl)]
        x = np.asarray(calibrated.gwf.modelgrid.xcellcenters, dtype=float).reshape(-1)
        y = np.asarray(calibrated.gwf.modelgrid.ycellcenters, dtype=float).reshape(-1)
        return gpd.GeoDataFrame(
            {
                "cell": np.arange(ncpl, dtype=int),
                "x": x,
                "y": y,
                "k_initial": k_initial,
                "k_final": k_final,
                "k_ratio": k_ratio,
            },
            geometry=polygons,
            crs=getattr(calibrated.gwf.modelgrid, "crs", None),
        )

    def compare_head_targets(self, targets: HeadTargets | None = None, *, prefix: str | None = None) -> pd.DataFrame:
        """Compare one set of head targets against baseline and calibrated models."""

        if targets is None:
            targets = self.load_head_targets(prefix=prefix)
        baseline_compare = targets.compare(self.load_baseline_model()).rename(
            columns={
                "sim_head": "sim_head_baseline",
                "residual": "residual_baseline",
                "abs_residual": "abs_residual_baseline",
            }
        )
        calibrated_compare = targets.compare(self.load_calibrated_model()).rename(
            columns={
                "sim_head": "sim_head_calibrated",
                "residual": "residual_calibrated",
                "abs_residual": "abs_residual_calibrated",
            }
        )
        merged = baseline_compare.merge(
            calibrated_compare[
                [
                    "name",
                    "time",
                    "sim_head_calibrated",
                    "residual_calibrated",
                    "abs_residual_calibrated",
                ]
            ],
            on=["name", "time"],
            how="inner",
        )
        merged["delta_residual"] = merged["residual_calibrated"] - merged["residual_baseline"]
        merged["abs_residual_improvement"] = (
            merged["abs_residual_baseline"] - merged["abs_residual_calibrated"]
        )
        return merged

    def compare_head_target_stats(
        self,
        targets: HeadTargets | None = None,
        *,
        prefix: str | None = None,
    ) -> pd.DataFrame:
        """Return summary residual stats before and after calibration."""

        if targets is None:
            targets = self.load_head_targets(prefix=prefix)
        baseline_stats = targets.stats(self.load_baseline_model()).rename(
            columns={
                "mean_error": "mean_error_baseline",
                "mae": "mae_baseline",
                "rmse": "rmse_baseline",
            }
        )
        calibrated_stats = targets.stats(self.load_calibrated_model()).rename(
            columns={
                "mean_error": "mean_error_calibrated",
                "mae": "mae_calibrated",
                "rmse": "rmse_calibrated",
            }
        )
        stats = baseline_stats.merge(calibrated_stats, on="n", how="inner")
        stats["mae_improvement"] = stats["mae_baseline"] - stats["mae_calibrated"]
        stats["rmse_improvement"] = stats["rmse_baseline"] - stats["rmse_calibrated"]
        return stats

    def review(
        self,
        targets: HeadTargets | None = None,
        *,
        prefix: str | None = None,
    ) -> PestRunReview:
        """Assemble the common review objects for a completed calibration run."""

        if targets is None:
            targets = self.load_head_targets(prefix=prefix)
        baseline_model = self.load_baseline_model()
        calibrated_model = self.load_calibrated_model()
        residual_compare = self.compare_head_targets(targets)
        stats = self.compare_head_target_stats(targets)
        k_geodata = self.k_geodata()
        return PestRunReview(
            run=self,
            targets=targets,
            baseline_model=baseline_model,
            calibrated_model=calibrated_model,
            residual_compare=residual_compare,
            stats=stats,
            k_geodata=k_geodata,
        )

    def plot_k(self, *, ax=None, legend: bool = True, **kwargs):
        """Plot final calibrated ``K`` values and return the matplotlib axis."""

        k_gdf = self.k_geodata()
        if ax is None:
            _, ax = plt.subplots()
        plot_kwargs = {"column": "k_final", "legend": legend, **kwargs}
        k_gdf.plot(ax=ax, **plot_kwargs)
        ax.set_title("Final K")
        ax.set_axis_off()
        return ax

    def plot_k_ratio(self, *, ax=None, legend: bool = True, **kwargs):
        """Plot ``K_final / K_initial`` and return the matplotlib axis."""

        k_gdf = self.k_geodata()
        if ax is None:
            _, ax = plt.subplots()
        plot_kwargs = {"column": "k_ratio", "legend": legend, **kwargs}
        k_gdf.plot(ax=ax, **plot_kwargs)
        ax.set_title("K final / K initial")
        ax.set_axis_off()
        return ax

    def plot_residuals_by_period(
        self,
        targets: HeadTargets | None = None,
        *,
        prefix: str | None = None,
        ax=None,
    ):
        """Plot baseline vs calibrated MAE by period/time."""

        residual_compare = self.compare_head_targets(targets, prefix=prefix)
        period_stats = (
            residual_compare.groupby("time")
            .agg(
                mae_baseline=("abs_residual_baseline", "mean"),
                mae_calibrated=("abs_residual_calibrated", "mean"),
            )
            .reset_index()
            .sort_values("time")
        )
        if ax is None:
            _, ax = plt.subplots()
        period_stats.plot(x="time", y="mae_baseline", marker="o", ax=ax, label="Baseline MAE")
        period_stats.plot(x="time", y="mae_calibrated", marker="o", ax=ax, label="Calibrated MAE")
        ax.set_title("Residual MAE by period")
        ax.set_xlabel("Time / period")
        ax.set_ylabel("Mean absolute error")
        ax.grid(True, alpha=0.3)
        return ax

    def plot_obs_vs_sim(
        self,
        targets: HeadTargets | None = None,
        *,
        prefix: str | None = None,
        ax=None,
    ):
        """Plot observed heads against baseline and calibrated simulated heads."""

        residual_compare = self.compare_head_targets(targets, prefix=prefix)
        if ax is None:
            _, ax = plt.subplots()
        ax.scatter(
            residual_compare["head_target"],
            residual_compare["sim_head_baseline"],
            label="Baseline",
            alpha=0.8,
        )
        ax.scatter(
            residual_compare["head_target"],
            residual_compare["sim_head_calibrated"],
            label="Calibrated",
            alpha=0.8,
        )
        values = pd.concat(
            [
                residual_compare["head_target"],
                residual_compare["sim_head_baseline"],
                residual_compare["sim_head_calibrated"],
            ],
            axis=0,
        ).dropna()
        if not values.empty:
            lower = float(values.min())
            upper = float(values.max())
            ax.plot([lower, upper], [lower, upper], linestyle="--", color="black", linewidth=1)
        ax.set_title("Observed vs simulated heads")
        ax.set_xlabel("Observed head")
        ax.set_ylabel("Simulated head")
        ax.grid(True, alpha=0.3)
        ax.legend()
        return ax

    def plot_well_timeseries(
        self,
        name: str,
        targets: HeadTargets | None = None,
        *,
        prefix: str | None = None,
        ax=None,
    ):
        """Plot one observation location through time for baseline and calibrated runs."""

        residual_compare = self.compare_head_targets(targets, prefix=prefix)
        well_frame = residual_compare.loc[
            residual_compare["name"].astype(str).str.lower() == str(name).strip().lower()
        ].sort_values("time")
        if well_frame.empty:
            raise ValueError(f"No residual rows found for observation name {name!r}.")
        if ax is None:
            _, ax = plt.subplots()
        well_frame.plot(x="time", y="head_target", marker="o", ax=ax, label="Target")
        well_frame.plot(x="time", y="sim_head_baseline", marker="o", ax=ax, label="Baseline")
        well_frame.plot(x="time", y="sim_head_calibrated", marker="o", ax=ax, label="Calibrated")
        ax.set_title(str(well_frame['name'].iloc[0]))
        ax.set_xlabel("Time / period")
        ax.set_ylabel("Head")
        ax.grid(True, alpha=0.3)
        return ax


def open_pest_run(workspace: str | Path, *, crs: str = "EPSG:2927", verbosity_level: int = 0) -> PestRunResults:
    """Open a PEST artifact root or direct ``pest`` workspace for evaluation."""

    return PestRunResults.from_workspace(workspace, crs=crs, verbosity_level=verbosity_level)
