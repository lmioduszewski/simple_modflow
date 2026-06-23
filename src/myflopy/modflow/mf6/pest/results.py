"""Reopen and evaluate completed ``myflopy`` PEST workspaces."""

from __future__ import annotations
from myflopy.viz import mpl_axes

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.io as pio
from shapely.geometry import Polygon

from myflopy.modflow.calcs.calibration import CalibrationPlot
from myflopy.modflow.mf6.observations import (
    DrnFlowTargets,
    HeadTargets,
    LakeStageTargets,
    SfrFlowTargets,
    SfrStageTargets,
)
from myflopy.modflow.mf6.pest.forward_run import _apply_drain_specs, _apply_k_specs
from myflopy.project.run_model import LoadedMf6Run, load_mf6_run

METADATA_FILENAME = "myflopy_pest_metadata.json"
MATERIALIZATION_FILENAME = "myflopy_calibrated_materialization.json"
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
    """Return the saved myflopy PEST metadata file when present."""

    candidate = pest_workspace / METADATA_FILENAME
    return candidate if candidate.exists() else None


def _read_saved_locations(path: Path) -> pd.DataFrame | gpd.GeoDataFrame:
    """Read one saved target-definition table from CSV or GPKG."""

    path = Path(path)
    if path.suffix.lower() == ".gpkg":
        return gpd.read_file(path)
    frame = pd.read_csv(path)
    if "cells" in frame.columns:
        frame["cells"] = frame["cells"].fillna("").apply(
            lambda text: [int(value) for value in str(text).split(",") if str(value).strip() != ""]
        )
    return frame


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


def _sanitize_filename_token(value: str) -> str:
    """Return a filesystem-friendly token for figure and export filenames."""

    token = "".join(character if str(character).isalnum() else "_" for character in str(value).strip())
    token = token.strip("_")
    return token or "item"


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
        FloPy/myflopy load verbosity.
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
        """Return saved myflopy PEST metadata when available."""

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

    def _saved_observation_entry(self, *, kind: str, prefix: str | None = None) -> dict:
        """Return one saved observation metadata entry by kind/prefix."""

        display_kind = {
            "head_targets": "head-target",
            "lake_stage": "lake-stage",
            "sfr_stage": "SFR-stage",
            "sfr_flow": "SFR-flow",
            "drn_flow": "DRN-flow",
        }.get(str(kind), str(kind).replace("_", "-"))
        entries = self.saved_observation_sets(kind=kind)
        if prefix is not None:
            prefix_key = str(prefix).strip().lower()
            entries = [
                item
                for item in entries
                if str(item.get("prefix", "")).strip().lower() == prefix_key
            ]
        if not entries:
            raise FileNotFoundError(f"No saved {display_kind} metadata was found in the PEST workspace.")
        if len(entries) > 1:
            raise ValueError(f"Multiple saved {display_kind} sets are available. Pass prefix=... to select one.")
        return entries[0]

    def load_head_targets(self, *, prefix: str | None = None) -> HeadTargets:
        """Load a saved head-target dataset from the workspace metadata."""

        entry = self._saved_observation_entry(kind="head_targets", prefix=prefix)
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

    def load_lake_stage_targets(self, *, prefix: str | None = None) -> LakeStageTargets:
        """Load a saved lake-stage target dataset from the workspace metadata."""

        entry = self._saved_observation_entry(kind="lake_stage", prefix=prefix)
        locations = _read_saved_locations(self.pest_workspace / entry["locations_file"])
        values = self.pest_workspace / entry["values_file"]
        return LakeStageTargets(
            locations=locations,
            values=values,
            time_column=entry.get("time_column", "time"),
            value_column=entry.get("value_column", "stage"),
        )

    def load_sfr_stage_targets(self, *, prefix: str | None = None) -> SfrStageTargets:
        """Load a saved SFR-stage target dataset from the workspace metadata."""

        entry = self._saved_observation_entry(kind="sfr_stage", prefix=prefix)
        locations = _read_saved_locations(self.pest_workspace / entry["locations_file"])
        values = self.pest_workspace / entry["values_file"]
        return SfrStageTargets(
            locations=locations,
            values=values,
            time_column=entry.get("time_column", "time"),
            value_column=entry.get("value_column", "stage_target"),
        )

    def load_sfr_flow_targets(self, *, prefix: str | None = None) -> SfrFlowTargets:
        """Load a saved SFR-flow target dataset from the workspace metadata."""

        entry = self._saved_observation_entry(kind="sfr_flow", prefix=prefix)
        locations = _read_saved_locations(self.pest_workspace / entry["locations_file"])
        values = self.pest_workspace / entry["values_file"]
        return SfrFlowTargets(
            locations=locations,
            values=values,
            time_column=entry.get("time_column", "time"),
            value_column=entry.get("value_column", "flow_target"),
        )

    def load_drn_flow_targets(self, *, prefix: str | None = None) -> DrnFlowTargets:
        """Load a saved DRN seepage-zone target dataset from the workspace metadata."""

        entry = self._saved_observation_entry(kind="drn_flow", prefix=prefix)
        locations = _read_saved_locations(self.pest_workspace / entry["locations_file"])
        values = self.pest_workspace / entry["values_file"]
        return DrnFlowTargets(
            locations=locations,
            values=values,
            time_column=entry.get("time_column", "time"),
            value_column=entry.get("value_column", "flow_target"),
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
        gdf = gpd.GeoDataFrame(
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
        if gdf.crs is None and getattr(calibrated, "_crs", None) is not None:
            gdf = gdf.set_crs(calibrated._crs, allow_override=True)
        return gdf

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
            _, ax = mpl_axes()
        plot_kwargs = {"column": "k_final", "legend": legend, **kwargs}
        k_gdf.plot(ax=ax, **plot_kwargs)
        ax.set_title("Final K")
        ax.set_axis_off()
        return ax

    def plot_k_ratio(self, *, ax=None, legend: bool = True, **kwargs):
        """Plot ``K_final / K_initial`` and return the matplotlib axis."""

        k_gdf = self.k_geodata()
        if ax is None:
            _, ax = mpl_axes()
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
        backend: str = "plotly",
    ):
        """Return a by-period residual summary plot for baseline vs calibrated fits.

        ``backend`` selects ``"plotly"`` (interactive, default) or
        ``"matplotlib"`` (static matplotlib/seaborn).
        """

        residual_compare = self.compare_head_targets(targets, prefix=prefix)
        current = residual_compare.rename(
            columns={
                "sim_head_calibrated": "sim_head",
                "residual_calibrated": "residual",
                "abs_residual_calibrated": "abs_residual",
            }
        )
        baseline = residual_compare.rename(
            columns={
                "sim_head_baseline": "sim_head",
                "residual_baseline": "residual",
                "abs_residual_baseline": "abs_residual",
            }
        )
        return CalibrationPlot.from_residuals_by_period(
            current,
            baseline_compare=baseline,
            backend=backend,
        )

    def plot_obs_vs_sim(
        self,
        targets: HeadTargets | None = None,
        *,
        prefix: str | None = None,
        backend: str = "plotly",
    ):
        """Return an observed-vs-simulated plot for baseline and calibrated fits.

        ``backend`` selects ``"plotly"`` (interactive, default) or
        ``"matplotlib"`` (static matplotlib/seaborn).
        """

        residual_compare = self.compare_head_targets(targets, prefix=prefix)
        current = residual_compare.rename(
            columns={
                "sim_head_calibrated": "sim_head",
                "residual_calibrated": "residual",
                "abs_residual_calibrated": "abs_residual",
            }
        )
        baseline = residual_compare.rename(
            columns={
                "sim_head_baseline": "sim_head",
                "residual_baseline": "residual",
                "abs_residual_baseline": "abs_residual",
            }
        )
        return CalibrationPlot.from_obs_vs_sim(
            current,
            baseline_compare=baseline,
            backend=backend,
        )

    def plot_well_timeseries(
        self,
        name: str,
        targets: HeadTargets | None = None,
        *,
        prefix: str | None = None,
        backend: str = "plotly",
    ):
        """Return one observation location through time for baseline and calibrated runs.

        ``backend`` selects ``"plotly"`` (interactive, default) or
        ``"matplotlib"`` (static matplotlib/seaborn).
        """

        residual_compare = self.compare_head_targets(targets, prefix=prefix)
        current = residual_compare.rename(
            columns={
                "sim_head_calibrated": "sim_head",
                "residual_calibrated": "residual",
                "abs_residual_calibrated": "abs_residual",
            }
        )
        baseline = residual_compare.rename(
            columns={
                "sim_head_baseline": "sim_head",
                "residual_baseline": "residual",
                "abs_residual_baseline": "abs_residual",
            }
        )
        return CalibrationPlot.from_timeseries(
            current,
            name=name,
            baseline_compare=baseline,
            backend=backend,
        )

    def export_review(
        self,
        folder: str | Path,
        targets: HeadTargets | None = None,
        *,
        prefix: str | None = None,
        well_names: list[str] | tuple[str, ...] | None = None,
        timeseries_names: list[str] | tuple[str, ...] | None = None,
        max_wells: int = 5,
    ) -> dict:
        """Export a compact review bundle of tables, geodata, and figures."""

        export_root = Path(folder)
        export_root.mkdir(parents=True, exist_ok=True)
        figures_dir = export_root / "figures"
        figures_dir.mkdir(parents=True, exist_ok=True)

        review = self.review(targets=targets, prefix=prefix)
        residual_compare = review.residual_compare
        stats = review.stats
        k_gdf = review.k_geodata.copy()
        if k_gdf.crs is None and self.crs is not None:
            k_gdf = k_gdf.set_crs(self.crs, allow_override=True)

        residual_compare_path = export_root / "residual_compare.csv"
        stats_path = export_root / "residual_stats.csv"
        k_gdf_path = export_root / "k_review.gpkg"

        residual_compare.to_csv(residual_compare_path, index=False)
        stats.to_csv(stats_path, index=False)
        k_gdf.to_file(k_gdf_path, driver="GPKG")

        fig_obs_vs_sim = self.plot_obs_vs_sim(review.targets)
        fig_period = self.plot_residuals_by_period(review.targets)
        obs_vs_sim_path = figures_dir / "obs_vs_sim.html"
        residuals_by_period_path = figures_dir / "residuals_by_period.html"
        pio.write_html(fig_obs_vs_sim, file=obs_vs_sim_path, include_plotlyjs="cdn")
        pio.write_html(fig_period, file=residuals_by_period_path, include_plotlyjs="cdn")

        ax_k = self.plot_k()
        k_png_path = figures_dir / "k_final.png"
        ax_k.figure.savefig(k_png_path, dpi=200, bbox_inches="tight")
        plt.close(ax_k.figure)

        ax_k_ratio = self.plot_k_ratio()
        k_ratio_png_path = figures_dir / "k_ratio.png"
        ax_k_ratio.figure.savefig(k_ratio_png_path, dpi=200, bbox_inches="tight")
        plt.close(ax_k_ratio.figure)

        unique_names = residual_compare["name"].astype(str).dropna().drop_duplicates().tolist()
        if well_names is None and timeseries_names is not None:
            well_names = [str(name) for name in timeseries_names]

        if well_names is None:
            selected_wells = unique_names[: max(0, int(max_wells))]
        else:
            selected_wells = [str(name) for name in well_names]

        exported_well_paths: dict[str, str] = {}
        for name in selected_wells:
            fig = self.plot_well_timeseries(name, review.targets)
            filename = f"timeseries_{_sanitize_filename_token(name)}.html"
            path = figures_dir / filename
            pio.write_html(fig, file=path, include_plotlyjs="cdn")
            exported_well_paths[name] = str(path)

        manifest = {
            "export_root": str(export_root),
            "figures_dir": str(figures_dir),
            "residual_compare_csv": str(residual_compare_path),
            "residual_stats_csv": str(stats_path),
            "k_review_gpkg": str(k_gdf_path),
            "obs_vs_sim_html": str(obs_vs_sim_path),
            "residuals_by_period_html": str(residuals_by_period_path),
            "k_final_png": str(k_png_path),
            "k_ratio_png": str(k_ratio_png_path),
            "well_timeseries_html": exported_well_paths,
            "n_exported_wells": len(exported_well_paths),
        }
        manifest_path = export_root / "review_manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        manifest["manifest_json"] = str(manifest_path)
        return manifest


def open_pest_run(workspace: str | Path, *, crs: str = "EPSG:2927", verbosity_level: int = 0) -> PestRunResults:
    """Open a PEST artifact root or direct ``pest`` workspace for evaluation."""

    return PestRunResults.from_workspace(workspace, crs=crs, verbosity_level=verbosity_level)
