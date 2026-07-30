from __future__ import annotations

import os
import subprocess
import sys
import warnings
from pathlib import Path
from types import SimpleNamespace

os.environ.setdefault("MPLBACKEND", "Agg")

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import Point, Polygon

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

PYEMU_DEPS = Path(r"C:\Users\lukem\Python\GMDSI_notebooks-main\dependencies")
PYEMU_PACKAGE = PYEMU_DEPS / "pyemu"
try:
    import pyemu as _pyemu  # noqa: F401
except Exception:
    for dep_path in (PYEMU_PACKAGE, PYEMU_DEPS):
        if dep_path.exists() and str(dep_path) not in sys.path:
            sys.path.insert(0, str(dep_path))

from myflopy import viz  # noqa: E402
from myflopy.modflow.calcs.calibration import CalibrationPlot  # noqa: E402
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig  # noqa: E402
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus  # noqa: E402
from myflopy.modflow.mf6.observations import (  # noqa: E402
    DrnFlowTargets,
    HeadTargets,
    LakeStageTargets,
    SfrFlowTargets,
    SfrStageTargets,
    TargetRegistry,
)
from myflopy.modflow.mf6.pest.gis import derive_bounds  # noqa: E402
from myflopy.modflow.mf6.pest.observations import (  # noqa: E402
    prepare_drn_flow_observations,
    prepare_lake_stage_observations,
    prepare_sfr_flow_observations,
    prepare_sfr_stage_observations,
)
from myflopy.modflow.mf6.pest.project import PestProject  # noqa: E402
from myflopy.modflow.mf6.pest.specs import (  # noqa: E402
    DrnFlowObservationSpec,
    LakeStageObservationSpec,
    SfrFlowObservationSpec,
    SfrStageObservationSpec,
)
from myflopy.modflow.mf6.simulation.base import SimulationBase  # noqa: E402
from myflopy.modflow.mf6.simulation.discretization import (  # noqa: E402
    DisvGrid,
    TemporalDiscretization,
)
from myflopy.modflow.mf6.simulation.packages import (  # noqa: E402
    CHD,
    InitialConditions,
    KFlow,
    OutputControl,
    Recharge,
    Storage,
)


def _project_temp_dir(name: str) -> Path:
    root = ROOT / ".pytest-work" / name
    if root.exists():
        import shutil

        shutil.rmtree(root, ignore_errors=True)
    root.mkdir(parents=True, exist_ok=True)
    return root


def _write_gpkg(path: Path, gdf: gpd.GeoDataFrame) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()
    gdf.to_file(path, driver="GPKG")
    return path


def _write_synthetic_par(path: Path, parameter_frame: pd.DataFrame, value_by_prefix: dict[str, float]) -> Path:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("single point\n")
        for row in parameter_frame.itertuples(index=False):
            parnme = str(row.parnme)
            lower = parnme.lower()
            value = float(row.parval1)
            for prefix, replacement in value_by_prefix.items():
                if lower.startswith(prefix.lower()):
                    value = float(replacement)
                    break
            handle.write(f"    {lower:<30} {value:<20.12g} {float(row.parval1):<20.12g} 0.0\n")
    return path


def _two_cell_vor_clockwise():
    verts = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [2.0, 0.0],
            [2.0, 1.0],
        ],
        dtype=float,
    )
    iverts = [[0, 3, 2, 1], [1, 2, 5, 4]]
    xcyc = np.array([[0.5, 0.5], [1.5, 0.5]], dtype=float)
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)


def _two_cell_model(name: str, workspace: Path, *, nper: int = 2):
    vor = _two_cell_vor_clockwise()
    vor.gdf_topbtm = gpd.GeoDataFrame(
        {
            "geometry": vor.gdf_vorPolys.geometry,
            0: [10.0, 12.0],
            1: [2.0, 4.0],
        },
        geometry="geometry",
        crs=vor.crs,
    )
    model = SimulationBase(name=name, mf_folder_path=workspace, vor=vor, nper=nper)
    DisvGrid(vor=vor, model=model, top=[10.0, 12.0], bottom=[[2.0, 4.0]], nlay=1, idomain=[[1, 1]])
    TemporalDiscretization(model=model, per_len=1, num_steps=1, multiplier=1.0)
    KFlow(model=model, k=[1.0, 1.0], k33_vert=[0.25, 0.4])
    return model, vor


def _build_two_cell_pest_forward_model(name: str, workspace: Path):
    vor = _two_cell_vor_clockwise()
    vor.gdf_topbtm = gpd.GeoDataFrame(
        {
            "geometry": vor.gdf_vorPolys.geometry,
            0: [10.0, 10.0],
            1: [0.0, 0.0],
        },
        geometry="geometry",
        crs=vor.crs,
    )
    model = SimulationBase(name=name, mf_folder_path=workspace, vor=vor, nper=1)
    DisvGrid(vor=vor, model=model, top=[10.0, 10.0], bottom=[[0.0, 0.0]], nlay=1, idomain=[[1, 1]])
    TemporalDiscretization(model=model, per_len=1, num_steps=1, multiplier=1.0)
    InitialConditions(model=model, vor=vor, nlay=1, strt=[10.0, 9.0])
    KFlow(model=model, k=[5.0, 5.0], k33_vert=[1.0, 1.0], save_specific_discharge=False)
    Storage(model=model, sto_steady={0: True}, sto_transient={})
    OutputControl(model=model)
    CHD(model=model, stress_period_data={0: [[(0, 0), 10.0]]})
    return model, vor


def _fake_model_with_heads():
    vor = _two_cell_vor_clockwise()
    index = pd.MultiIndex.from_tuples(
        [((0, 0), 0, 0), ((0, 0), 0, 1), ((0, 1), 0, 0), ((0, 1), 0, 1)],
        names=["kstpkper", "layer", "cell"],
    )
    heads = pd.DataFrame({"elev": [10.0, 9.5, 10.25, 9.25]}, index=index)
    return SimpleNamespace(vor=vor, all_heads=heads)


def _fake_surface_water_model():
    vor = _two_cell_vor_clockwise()
    lake_stage = pd.DataFrame(
        {
            "per": [0, 1, 0, 1],
            "lake": [0, 0, 1, 1],
            "stage": [100.0, 101.0, 97.5, 97.0],
        }
    )
    sfr_stage = pd.DataFrame(
        {
            "per": [0, 1],
            "reach": [0, 0],
            "stage": [12.0, 11.5],
        }
    )
    sfr_flow = pd.DataFrame(
        {
            "per": [0, 1],
            "reach": [0, 0],
            "q": [1.25, 1.5],
        }
    )
    drn_df = pd.DataFrame(
        {
            "kstpkper": [(0, 0), (0, 0), (0, 1), (0, 1)],
            "node": [0, 1, 0, 1],
            "q": [0.2, 0.8, 0.25, 0.85],
        }
    ).set_index(["kstpkper", "node"])

    model = SimpleNamespace(
        vor=vor,
        gwf=SimpleNamespace(lak="lak_owner", sfr="sfr_owner", drn="drn_owner"),
        packages=SimpleNamespace(
            lak=SimpleNamespace(results=SimpleNamespace(stage=SimpleNamespace(get=lambda: lake_stage.copy()))),
            sfr=SimpleNamespace(
                results=SimpleNamespace(
                    stage=SimpleNamespace(get=lambda: sfr_stage.copy()),
                    q=SimpleNamespace(get=lambda: sfr_flow.copy()),
                )
            ),
        ),
    )
    model.bud = lambda package: SimpleNamespace(df=drn_df.copy())
    model.targets = TargetRegistry(model)
    return model


def test_head_targets_normalize_match_and_compare():
    locations = gpd.GeoDataFrame(
        {
            "name": ["OBS_A", "OBS_B"],
            "layer": [0, 0],
            "group": ["wells", "wells"],
            "weight": [1.0, 2.0],
        },
        geometry=[Point(0.5, 0.5), Point(1.5, 0.5)],
        crs="EPSG:2927",
    )
    values = pd.DataFrame(
        {
            "per": [0, 1],
            "OBS_A": [9.75, 10.0],
            "OBS_B": [9.0, 9.5],
        }
    )
    targets = HeadTargets(locations=locations, values=values, time_column="per")

    long = targets.to_long()
    assert {"name", "time", "head_target", "weight", "group", "layer"}.issubset(long.columns)
    assert len(long) == 4

    fake_model = _fake_model_with_heads()
    matches = targets.match_to_model(fake_model)
    assert matches.sort_values("name")["cell"].tolist() == [0, 1]

    comparison = targets.compare(fake_model)
    assert pytest.approx(comparison.loc[comparison["name"] == "OBS_A", "sim_head"].iloc[0]) == 10.0
    assert "residual" in comparison.columns
    stats = targets.stats(fake_model)
    assert stats.loc[0, "n"] == 4


def test_model_bound_targets_registry_and_calibration_plot():
    locations = gpd.GeoDataFrame(
        {
            "name": ["OBS_A", "OBS_B"],
            "layer": [0, 0],
            "group": ["wells", "wells"],
            "weight": [1.0, 2.0],
        },
        geometry=[Point(0.5, 0.5), Point(1.5, 0.5)],
        crs="EPSG:2927",
    )
    values = pd.DataFrame({"per": [0, 1], "OBS_A": [9.75, 10.0], "OBS_B": [9.0, 9.5]})
    targets = HeadTargets(locations=locations, values=values, time_column="per")
    fake_model = _fake_model_with_heads()
    fake_model.targets = TargetRegistry(fake_model)
    fake_model.targets.heads = targets

    assert "heads" in fake_model.targets
    assert fake_model.targets.keys() == ["heads"]
    bound = fake_model.targets.heads
    assert bound.targets is targets

    standalone = targets.compare(fake_model)
    bound_compare = bound.compare()
    pd.testing.assert_frame_equal(bound_compare, standalone)

    fig = bound.calibration_plot()
    assert isinstance(fig, CalibrationPlot)
    assert len(fig.data) >= 2

    ax_locations = bound.plot.locations()
    assert ax_locations.get_title() == "Head Target Locations"

    obs_fig = bound.plot.obs_vs_sim()
    assert isinstance(obs_fig, CalibrationPlot)
    assert obs_fig.layout.title.text == "Observed vs simulated heads"
    assert len(obs_fig.data) == 2

    ts_fig = bound.plot.timeseries("OBS_A")
    assert isinstance(ts_fig, CalibrationPlot)
    assert ts_fig.layout.title.text == "OBS_A"
    assert len(ts_fig.data) == 2

    ax_heads = bound.plot.calibration(type="heads")
    assert isinstance(ax_heads, CalibrationPlot)

    period_fig = bound.plot.residuals_by_period()
    assert isinstance(period_fig, CalibrationPlot)
    assert period_fig.layout.title.text == "Residual MAE by period"
    assert len(period_fig.data) == 1

    # The same plots are available as static matplotlib/seaborn figures.
    from matplotlib.figure import Figure as MplFigure

    assert isinstance(bound.plot.obs_vs_sim(backend="matplotlib"), MplFigure)
    assert isinstance(bound.plot.timeseries("OBS_A", backend="matplotlib"), MplFigure)
    assert isinstance(bound.plot.residuals_by_period(backend="matplotlib"), MplFigure)

    baseline_model = _fake_model_with_heads()
    baseline_model.all_heads = baseline_model.all_heads.copy()
    baseline_model.all_heads.loc[:, "elev"] = [9.8, 9.3, 10.0, 9.0]
    obs_fig_with_baseline = bound.plot.obs_vs_sim(baseline=baseline_model)
    assert len(obs_fig_with_baseline.data) == 3
    ts_fig_with_baseline = bound.plot.timeseries("OBS_A", baseline=baseline_model)
    assert len(ts_fig_with_baseline.data) == 3
    period_fig_with_baseline = bound.plot.residuals_by_period(baseline=baseline_model)
    assert len(period_fig_with_baseline.data) == 2


def test_calibration_plot_target_driven_constructors():
    compare = pd.DataFrame(
        {
            "name": ["OBS_A", "OBS_A", "OBS_B", "OBS_B"],
            "time": [0, 1, 0, 1],
            "per": [0, 1, 0, 1],
            "layer": [0, 0, 0, 0],
            "head_target": [10.0, 10.2, 9.5, 9.7],
            "sim_head": [9.8, 10.1, 9.4, 9.9],
            "residual": [-0.2, -0.1, -0.1, 0.2],
            "abs_residual": [0.2, 0.1, 0.1, 0.2],
        }
    )
    baseline = compare.copy()
    baseline["sim_head"] = [9.4, 9.5, 9.0, 9.2]
    baseline["residual"] = baseline["sim_head"] - baseline["head_target"]
    baseline["abs_residual"] = baseline["residual"].abs()

    obs_vs_sim = CalibrationPlot.from_obs_vs_sim(compare, baseline_compare=baseline)
    assert obs_vs_sim.layout.title.text == "Observed vs simulated heads"
    assert len(obs_vs_sim.data) == 3

    ts_fig = CalibrationPlot.from_timeseries(compare, name="OBS_A", baseline_compare=baseline)
    assert ts_fig.layout.title.text == "OBS_A"
    assert len(ts_fig.data) == 3

    period_fig = CalibrationPlot.from_residuals_by_period(compare, baseline_compare=baseline)
    assert period_fig.layout.title.text == "Residual MAE by period"
    assert len(period_fig.data) == 2

    heads_fig = CalibrationPlot.from_compare(compare, type="heads")
    assert isinstance(heads_fig, CalibrationPlot)
    assert len(heads_fig.data) == 4


def test_head_targets_can_build_and_attach_flopy_obs(monkeypatch):
    workspace = _project_temp_dir("head_targets_obs_attach")
    model, _ = _two_cell_model("head_targets_obs_attach", workspace)
    locations = gpd.GeoDataFrame(
        {"name": ["OBS_A", "OBS_B"], "layer": [0, 0]},
        geometry=[Point(0.5, 0.5), Point(1.5, 0.5)],
        crs="EPSG:2927",
    )
    values = pd.DataFrame({"per": [0], "OBS_A": [10.0], "OBS_B": [9.0]})
    targets = HeadTargets(locations=locations, values=values, time_column="per")

    obs_dict = targets.to_flopy_obs(model, csv_name="cumb_obs.csv")
    assert list(obs_dict) == ["cumb_obs.csv"]
    assert obs_dict["cumb_obs.csv"] == [
        ("OBS_A", "HEAD", (0, 0)),
        ("OBS_B", "HEAD", (0, 1)),
    ]

    calls = {}

    def _fake_obs(owner, pname, continuous, filename):
        calls["owner"] = owner
        calls["pname"] = pname
        calls["continuous"] = continuous
        calls["filename"] = filename
        return SimpleNamespace(owner=owner, pname=pname, continuous=continuous, filename=filename)

    monkeypatch.setattr("flopy.mf6.modflow.mfutlobs.ModflowUtlobs", _fake_obs)
    package = targets.attach_flopy_obs(
        model,
        pname="gwf_obs",
        filename="cumb_heads.obs",
        csv_name="cumb_obs.csv",
    )
    assert package.filename == "cumb_heads.obs"
    assert calls["owner"] is model.gwf
    assert calls["continuous"] == obs_dict


def test_head_targets_from_cells_supports_dict_and_compare():
    fake_model = _fake_model_with_heads()
    targets = HeadTargets.from_cells(
        cells=[0, 1],
        names=["OBS_A", "OBS_B"],
        layers=0,
        times=[0, 1],
        values={
            "OBS_A": [9.75, 10.0],
            "OBS_B": [9.0, 9.5],
        },
        time_column="per",
    )

    long = targets.to_long()
    assert {"cell", "name", "time", "head_target"}.issubset(long.columns)
    assert sorted(long["cell"].unique().tolist()) == [0, 1]

    matches = targets.match_to_model(fake_model)
    assert matches.sort_values("name")["cell"].tolist() == [0, 1]

    comparison = targets.compare(fake_model)
    assert pytest.approx(comparison.loc[comparison["name"] == "OBS_B", "sim_head"].iloc[0]) == 9.5

    obs_dict = targets.to_flopy_obs(fake_model, csv_name="heads.csv")
    assert obs_dict["heads.csv"] == [
        ("OBS_A", "HEAD", (0, 0)),
        ("OBS_B", "HEAD", (0, 1)),
    ]


def test_head_targets_from_cells_supports_array_and_series():
    fake_model = _fake_model_with_heads()
    single = HeadTargets.from_cells(
        cells=[0],
        names=["OBS_A"],
        values=pd.Series([9.75, 10.0], index=[0, 1], name="head"),
        time_column="per",
    )
    assert len(single.to_long()) == 2
    assert single.compare(fake_model)["sim_head"].notna().all()

    multi = HeadTargets.from_cells(
        cells=[0, 1],
        names=["OBS_A", "OBS_B"],
        values=[[9.75, 9.0], [10.0, 9.5]],
        times=[0, 1],
        time_column="per",
    )
    assert len(multi.to_long()) == 4
    assert sorted(multi.to_wide().columns.tolist()) == ["OBS_A", "OBS_B", "time"]


def test_head_targets_from_records_supports_cell_records():
    fake_model = _fake_model_with_heads()
    targets = HeadTargets.from_records(
        [
            {"name": "OBS_A", "cell": 0, "layer": 0, "time": 0, "head": 9.75},
            {"name": "OBS_A", "cell": 0, "layer": 0, "time": 1, "head": 10.0},
            {"name": "OBS_B", "cell": 1, "layer": 0, "time": 0, "head": 9.0},
            {"name": "OBS_B", "cell": 1, "layer": 0, "time": 1, "head": 9.5},
        ],
        value_column="head",
    )
    comparison = targets.compare(fake_model)
    assert len(comparison) == 4
    assert comparison["cell"].tolist()[:2] == [0, 0]


def test_head_targets_plain_constructor_accepts_cell_dict_and_series_inputs():
    fake_model = _fake_model_with_heads()
    locations = {"OBS_A": 0, "OBS_B": 1}
    targets = HeadTargets(
        locations=locations,
        values={"per": [0, 1], "OBS_A": [9.75, 10.0], "OBS_B": [9.0, 9.5]},
        time_column="per",
    )
    assert sorted(targets.to_long()["cell"].unique().tolist()) == [0, 1]
    assert len(targets.compare(fake_model)) == 4

    series_targets = HeadTargets(
        locations=pd.Series([0], index=["OBS_A"]),
        values=pd.Series([9.75, 10.0], index=[0, 1], name="head"),
        time_column="per",
    )
    comparison = series_targets.compare(fake_model)
    assert comparison["sim_head"].notna().all()


def test_lake_stage_targets_can_build_and_attach_flopy_obs(monkeypatch):
    targets = LakeStageTargets({"deep_lake": 0, "shallow lake": 1})
    obs_dict = targets.to_flopy_obs(csv_name="lak_obs.csv")
    assert list(obs_dict) == ["lak_obs.csv"]
    assert obs_dict["lak_obs.csv"] == [
        ("deep_lake", "STAGE", 1),
        ("shallow_lake", "STAGE", 2),
    ]

    calls = {}
    fake_model = SimpleNamespace(gwf=SimpleNamespace(lak=object()))

    def _fake_obs(owner, pname, continuous, filename):
        calls["owner"] = owner
        calls["pname"] = pname
        calls["continuous"] = continuous
        calls["filename"] = filename
        return SimpleNamespace(owner=owner, pname=pname, continuous=continuous, filename=filename)

    monkeypatch.setattr("flopy.mf6.modflow.mfutlobs.ModflowUtlobs", _fake_obs)
    package = targets.attach_flopy_obs(fake_model, filename="cumb_lak.obs", csv_name="lak_obs.csv")
    assert package.filename == "cumb_lak.obs"
    assert calls["owner"] is fake_model.gwf.lak
    assert calls["continuous"] == obs_dict


def test_lake_stage_targets_from_series_and_records():
    series_targets = LakeStageTargets.from_series(
        lake="deep_lake",
        lake_id=0,
        values=pd.Series([766.1, 766.3], index=[0, 1], name="stage"),
        time_column="per",
    )
    long = series_targets.get()
    assert {"name", "lake", "time", "stage"}.issubset(long.columns)
    assert series_targets.summary().loc[0, "n_rows"] == 2

    record_targets = LakeStageTargets.from_records(
        [
            {"name": "deep_lake", "lake": 0, "time": 0, "stage": 766.1},
            {"name": "deep_lake", "lake": 0, "time": 1, "stage": 766.3},
            {"name": "shallow_lake", "lake": 1, "time": 0, "stage": 760.0},
        ]
    )
    assert record_targets.summary().loc[0, "n_lakes"] == 2
    obs_dict = record_targets.to_flopy_obs(csv_name="lak_obs.csv")
    assert len(obs_dict["lak_obs.csv"]) == 2


def test_lake_stage_targets_plain_constructor_accepts_simple_inputs():
    targets = LakeStageTargets(
        locations={"deep_lake": 0},
        values=[766.1, 766.3],
        times=[0, 1],
        time_column="per",
    )
    long = targets.get()
    assert long["name"].tolist() == ["deep_lake", "deep_lake"]
    assert long["stage"].tolist() == [766.1, 766.3]

    series_targets = LakeStageTargets(
        locations=pd.Series([0], index=["deep_lake"]),
        values=pd.Series([766.1, 766.3], index=[0, 1], name="stage"),
        time_column="per",
    )
    assert series_targets.summary().loc[0, "n_rows"] == 2


def test_named_surface_water_and_drn_targets_are_model_bound_and_plot_ready(monkeypatch):
    model = _fake_surface_water_model()

    lake_targets = LakeStageTargets(
        locations={"deep_lake": 0, "shallow_lake": 1},
        values={"per": [0, 1], "deep_lake": [99.5, 100.5], "shallow_lake": [97.75, 96.9]},
        time_column="per",
    )
    sfr_stage_targets = SfrStageTargets(
        locations={"reach_001": 0},
        values={"per": [0, 1], "reach_001": [11.75, 11.0]},
        time_column="per",
    )
    sfr_flow_targets = SfrFlowTargets(
        locations={"reach_001": 0},
        values={"per": [0, 1], "reach_001": [1.0, 1.8]},
        time_column="per",
        value_column="flow",
    )
    drn_targets = DrnFlowTargets(
        locations={"zone_west": [0], "zone_east": [1]},
        values={"per": [0, 1], "zone_west": [0.15, 0.3], "zone_east": [0.75, 0.9]},
        time_column="per",
        value_column="flow",
    )

    model.targets.lake_stage = lake_targets
    model.targets.sfr_stage = sfr_stage_targets
    model.targets.sfr_flow = sfr_flow_targets
    model.targets.drn_flow = drn_targets

    lake_compare = model.targets.lake_stage.compare()
    assert pytest.approx(lake_compare.loc[lake_compare["name"] == "deep_lake", "sim_stage"].iloc[0]) == 100.0
    assert model.targets.lake_stage.stats().loc[0, "n"] == 4
    lake_plot = model.targets.lake_stage.plot.obs_vs_sim()
    assert isinstance(lake_plot, CalibrationPlot)
    assert lake_plot.layout.title.text == "Observed vs simulated lake stage"

    sfr_stage_compare = model.targets.sfr_stage.compare()
    assert pytest.approx(sfr_stage_compare["sim_stage"].iloc[0]) == 12.0
    sfr_stage_plot = model.targets.sfr_stage.plot.timeseries("reach_001")
    assert isinstance(sfr_stage_plot, CalibrationPlot)
    assert len(sfr_stage_plot.data) == 2

    sfr_flow_compare = model.targets.sfr_flow.compare()
    assert pytest.approx(sfr_flow_compare["sim_flow"].iloc[1]) == 1.5
    sfr_flow_plot = model.targets.sfr_flow.plot.residuals_by_period()
    assert isinstance(sfr_flow_plot, CalibrationPlot)
    assert sfr_flow_plot.layout.title.text == "Residual MAE by period"

    drn_compare = model.targets.drn_flow.compare()
    assert pytest.approx(drn_compare.loc[drn_compare["name"] == "zone_east", "sim_flow"].iloc[0]) == 0.8
    assert model.targets.drn_flow.summary().loc[0, "n_cells"] == 2
    drn_plot = model.targets.drn_flow.plot.obs_vs_sim()
    assert isinstance(drn_plot, CalibrationPlot)
    assert drn_plot.layout.title.text == "Observed vs simulated DRN seepage"
    for compare in (lake_compare, sfr_stage_compare, sfr_flow_compare, drn_compare):
        assert "time" in compare.columns
        assert "time_x" not in compare.columns
        assert "time_y" not in compare.columns

    calls = {}

    def _fake_obs(owner, pname, continuous, filename):
        calls[filename] = {"owner": owner, "pname": pname, "continuous": continuous}
        return SimpleNamespace(owner=owner, pname=pname, continuous=continuous, filename=filename)

    monkeypatch.setattr("flopy.mf6.modflow.mfutlobs.ModflowUtlobs", _fake_obs)
    model.targets.sfr_stage.attach_flopy_obs(filename="sfr_stage.obs")
    model.targets.sfr_flow.attach_flopy_obs(filename="sfr_flow.obs")
    model.targets.drn_flow.attach_flopy_obs(filename="drn_flow.obs")

    assert calls["sfr_stage.obs"]["owner"] == "sfr_owner"
    assert calls["sfr_flow.obs"]["owner"] == "sfr_owner"
    assert calls["drn_flow.obs"]["owner"] == "drn_owner"
    assert calls["sfr_flow.obs"]["continuous"]["sfr_flow.csv"] == [
        ("reach_001", "DOWNSTREAM-FLOW", 1)
    ]
    drn_continuous = calls["drn_flow.obs"]["continuous"]
    assert list(drn_continuous) == ["drn_flow.csv"]
    assert len(drn_continuous["drn_flow.csv"]) == 2


def _fake_flow_ja_face_model():
    """A model whose SFR budget exposes FLOW-JA-FACE routing flow (+ TO-MVR).

    Exercises the *primary* path of ``SfrFlowTargets._package_flow_table`` — the one
    a real, properly-run model hits — which the mock in ``_fake_surface_water_model``
    (no ``.outputs``) never reaches. Node/rno numbering is 1-based, matching MF6.
    """

    flow_ja_face = pd.DataFrame(
        {
            "kstpkper": [(0, 0), (0, 0), (0, 0), (0, 0)],
            "node": [1, 1, 2, 3],       # 1-based; -> reaches 0, 0, 1, 2
            "node2": [2, 3, 1, 2],
            "q": [-4.0, -2.0, 4.0, -1.5],  # outbound rows are negative
        }
    )
    to_mvr = pd.DataFrame({"kstpkper": [(0, 0)], "node": [3], "q": [0.5]})

    def _get(text):
        if text == "FLOW-JA-FACE":
            return flow_ja_face.copy()
        if text == "TO-MVR":
            return to_mvr.copy()
        return None

    return SimpleNamespace(
        outputs=SimpleNamespace(sfr=SimpleNamespace(bud=SimpleNamespace(get=_get)))
    )


def test_sfr_flow_table_uses_flow_ja_face_routing_and_aggregates():
    # Covers the primary FLOW-JA-FACE branch: outbound rows summed per reach, plus
    # TO-MVR folded in. (Diverging reach 0 has two outbound rows; reach 2 also moves
    # water via TO-MVR; the inbound-only reach 1 is excluded.)
    targets = SfrFlowTargets(
        locations={"reach_000": 0, "reach_002": 2},
        values={"per": [0], "reach_000": [0.0], "reach_002": [0.0]},
        time_column="per",
        value_column="flow",
    )
    frame = targets._package_flow_table(_fake_flow_ja_face_model())
    by_reach = frame.set_index("reach")["sim_flow"]
    assert pytest.approx(by_reach.loc[0]) == 6.0   # 4.0 + 2.0 outbound rows summed
    assert pytest.approx(by_reach.loc[2]) == 2.0   # 1.5 outbound + 0.5 TO-MVR
    assert 1 not in by_reach.index                  # inbound-only reach dropped


def test_sfr_flow_table_refuses_to_substitute_exchange_for_missing_routing():
    # Ledger 49a: when FLOW-JA-FACE routing flow is unavailable on a real model
    # (e.g. SFRBudget.get returned the raw record list on a cadence mismatch), the
    # table must NOT relabel the stream-aquifer exchange (sfr.results.q) as sim_flow.
    exchange = pd.DataFrame({"per": [0], "reach": [0], "q_gwf": [123.0]})

    def _get(text):
        return None  # routing flow unavailable -> not a DataFrame

    model = SimpleNamespace(
        outputs=SimpleNamespace(sfr=SimpleNamespace(bud=SimpleNamespace(get=_get))),
        packages=SimpleNamespace(
            sfr=SimpleNamespace(
                results=SimpleNamespace(q=SimpleNamespace(get=lambda: exchange.copy()))
            )
        ),
    )
    targets = SfrFlowTargets(
        locations={"reach_000": 0},
        values={"per": [0], "reach_000": [0.0]},
        time_column="per",
        value_column="flow",
    )
    with pytest.warns(UserWarning, match="routing flow"):
        frame = targets._package_flow_table(model)
    assert frame.empty                              # no simulated values, not wrong ones
    assert list(frame.columns) == ["per", "reach", "sim_flow"]
    assert 123.0 not in frame["sim_flow"].values    # exchange never leaks in as flow


def test_surface_water_targets_prefer_direct_mf6_observation_csvs():
    model = _fake_surface_water_model()
    workspace = _project_temp_dir("surface_water_obs_csv_preference")
    model.workspace = workspace

    pd.DataFrame(
        {
            "time": [0, 1],
            "deep_lake": [101.25, 101.75],
            "shallow_lake": [96.5, 96.0],
        }
    ).to_csv(workspace / "gold_lakes.csv", index=False)
    pd.DataFrame(
        {
            "time": [0, 1],
            "reach_001": [12.4, 11.9],
        }
    ).to_csv(workspace / "gold_sfr_stage.csv", index=False)

    lake_targets = LakeStageTargets(
        locations={"deep_lake": 0, "shallow_lake": 1},
        values={"per": [0, 1], "deep_lake": [99.5, 100.5], "shallow_lake": [97.75, 96.9]},
        time_column="per",
    )
    sfr_stage_targets = SfrStageTargets(
        locations={"reach_001": 0},
        values={"per": [0, 1], "reach_001": [11.75, 11.0]},
        time_column="per",
    )

    lake_compare = lake_targets.compare(model)
    assert pytest.approx(
        lake_compare.loc[(lake_compare["name"] == "deep_lake") & (lake_compare["per"] == 0), "sim_stage"].iloc[0]
    ) == 101.25
    assert pytest.approx(
        lake_compare.loc[(lake_compare["name"] == "shallow_lake") & (lake_compare["per"] == 1), "sim_stage"].iloc[0]
    ) == 96.0

    sfr_stage_compare = sfr_stage_targets.compare(model)
    assert pytest.approx(
        sfr_stage_compare.loc[sfr_stage_compare["per"] == 0, "sim_stage"].iloc[0]
    ) == 12.4
    assert pytest.approx(
        sfr_stage_compare.loc[sfr_stage_compare["per"] == 1, "sim_stage"].iloc[0]
    ) == 11.9


def test_drn_flow_targets_accept_polygon_zones():
    model = _fake_surface_water_model()
    zone_targets = DrnFlowTargets(
        locations=gpd.GeoDataFrame(
            {"name": ["zone_west"]},
            geometry=[Polygon([(0.0, 0.0), (0.95, 0.0), (0.95, 0.95), (0.0, 0.95)])],
            crs=model.vor.crs,
        ),
        values={"per": [0, 1], "zone_west": [0.15, 0.3]},
        time_column="per",
        value_column="flow",
    )
    compare = zone_targets.compare(model)
    assert pytest.approx(compare["sim_flow"].iloc[0]) == 0.2
    zones = zone_targets.zone_definitions(model)
    assert zones["cells"].iloc[0] == [0]


def test_named_surface_water_and_drn_pest_observation_builders():
    model = _fake_surface_water_model()
    tmp_path = _project_temp_dir("named_surface_water_obs_builders")
    project = SimpleNamespace(
        model=model,
        template_workspace=tmp_path,
        pf=SimpleNamespace(calls=[]),
    )

    def _add_observations(*args, **kwargs):
        project.pf.calls.append((args, kwargs))

    project.pf.add_observations = _add_observations

    lake_targets = LakeStageTargets(
        locations={"deep_lake": 0},
        values={"per": [0, 1], "deep_lake": [99.5, 100.5]},
        time_column="per",
    )
    sfr_stage_targets = SfrStageTargets(
        locations={"reach_001": 0},
        values={"per": [0, 1], "reach_001": [11.75, 11.0]},
        time_column="per",
    )
    sfr_flow_targets = SfrFlowTargets(
        locations={"reach_001": 0},
        values={"per": [0, 1], "reach_001": [1.0, 1.8]},
        time_column="per",
        value_column="flow",
    )
    drn_targets = DrnFlowTargets(
        locations={"zone_west": [0], "zone_east": [1]},
        values={"per": [0, 1], "zone_west": [0.15, 0.3], "zone_east": [0.75, 0.9]},
        time_column="per",
        value_column="flow",
    )

    lake_prepared = prepare_lake_stage_observations(
        project,
        LakeStageObservationSpec(targets=lake_targets, prefix="lak_stage"),
    )
    sfr_stage_prepared = prepare_sfr_stage_observations(
        project,
        SfrStageObservationSpec(targets=sfr_stage_targets, prefix="sfr_stage"),
    )
    sfr_flow_prepared = prepare_sfr_flow_observations(
        project,
        SfrFlowObservationSpec(targets=sfr_flow_targets, prefix="sfr_flow"),
    )
    drn_prepared = prepare_drn_flow_observations(
        project,
        DrnFlowObservationSpec(targets=drn_targets, prefix="drn_flow"),
    )

    assert len(project.pf.calls) == 4
    assert lake_prepared["metadata"]["kind"] == "lake_stage"
    assert sfr_stage_prepared["metadata"]["kind"] == "sfr_stage"
    assert sfr_flow_prepared["metadata"]["kind"] == "sfr_flow"
    assert drn_prepared["metadata"]["kind"] == "drn_flow"
    assert lake_prepared["named_series_forward_run_config"]["output_csv"] == "lak_stage_simulated_lake_stage.csv"
    assert sfr_stage_prepared["named_series_forward_run_config"]["output_csv"] == "sfr_stage_simulated_sfr_stage.csv"
    assert sfr_flow_prepared["named_series_forward_run_config"]["output_csv"] == "sfr_flow_simulated_sfr_flow.csv"
    assert drn_prepared["named_series_forward_run_config"]["output_csv"] == "drn_flow_simulated_drn_flow.csv"
    assert (tmp_path / "lak_stage_target_values.csv").exists()
    assert (tmp_path / "sfr_stage_target_values.csv").exists()
    assert (tmp_path / "sfr_flow_target_values.csv").exists()
    assert (tmp_path / "drn_flow_target_values.csv").exists()


def test_drn_flow_targets_accept_csv_style_cell_lists():
    targets = DrnFlowTargets(
        locations=pd.DataFrame(
            {
                "name": ["zone_west", "zone_east"],
                "cells": ["0,1,2", "5,6"],
                "weight": [1.0, 0.5],
            }
        ),
        values={"per": [0], "zone_west": [0.2], "zone_east": [0.7]},
        time_column="per",
        value_column="flow",
    )

    definitions = targets.get()
    assert definitions.loc[definitions["name"] == "zone_west", "cells"].iloc[0] == [0, 1, 2]
    assert definitions.loc[definitions["name"] == "zone_east", "cells"].iloc[0] == [5, 6]


def test_derive_bounds_supports_all_first_slice_modes():
    base = pd.Series([10.0, 20.0])
    lower, upper = derive_bounds(base, (1.0, 5.0), "absolute")
    assert lower.tolist() == [1.0, 1.0]
    assert upper.tolist() == [5.0, 5.0]

    lower, upper = derive_bounds(base, (0.5, 2.0), "multiplier")
    assert lower.tolist() == [5.0, 10.0]
    assert upper.tolist() == [20.0, 40.0]

    lower, upper = derive_bounds(
        base,
        None,
        "from_columns",
        lower_bound_column=pd.Series([7.0, 8.0]),
        upper_bound_column=pd.Series([11.0, 22.0]),
    )
    assert lower.tolist() == [7.0, 8.0]
    assert upper.tolist() == [11.0, 22.0]

    lower, upper = derive_bounds(
        base,
        None,
        "multiplier_from_columns",
        lower_bound_column=pd.Series([0.1, 0.25]),
        upper_bound_column=pd.Series([10.0, 2.0]),
    )
    assert lower.tolist() == [1.0, 5.0]
    assert upper.tolist() == [100.0, 40.0]


def test_flatten_array_file_handles_mf6_wrapped_ragged_arrays(tmp_path):
    from myflopy.modflow.mf6.pest.native_parameters import _flatten_array_file

    # MF6 wraps arrays at a fixed width (e.g. 20 values/line), leaving a short
    # final line. numpy.loadtxt (used by pyEMU) rejects that ragged layout.
    values = np.arange(24, dtype=float)
    path = tmp_path / "npf_k.txt"
    path.write_text(
        " ".join(str(v) for v in values[:20]) + "\n"
        + " ".join(str(v) for v in values[20:]) + "\n"
    )
    with pytest.raises(ValueError):
        np.loadtxt(path)

    _flatten_array_file(path)

    loaded = np.loadtxt(path)
    assert loaded.shape == (24,)
    assert np.allclose(np.sort(loaded), values)


def test_model_pest_factory_binds_model_and_default_workspace(tmp_path):
    # The model is the single front door: model.pest(...) is the write-side
    # companion to model.pest_runs (the read side). It returns a PestProject
    # already bound to the model with the auto-discoverable default workspace.
    model, _ = _build_two_cell_pest_forward_model("pest_factory", tmp_path / "ws")

    cal = model.pest("calib", start_datetime="2021-03-01")
    assert isinstance(cal, PestProject)
    assert cal.model is model
    assert cal.start_datetime == "2021-03-01"
    # Default workspace lands BESIDE the model -- a sibling directory, not a
    # child -- so PstFrom's copy of the model workspace cannot swallow it
    # (ledger 107), while model.pest_runs still finds it.
    assert cal.template_workspace == model.workspace.parent / f"{model.workspace.name}.pest" / "calib"
    assert model.workspace not in cal.template_workspace.parents

    # Extra kwargs thread through to PestProject (e.g. an explicit workspace).
    custom = tmp_path / "elsewhere"
    cal2 = model.pest("calib2", start_datetime="2021-03-01", workspace=custom)
    assert cal2.template_workspace == custom


def test_model_pest_factory_infers_and_requires_start_datetime(tmp_path):
    model, _ = _build_two_cell_pest_forward_model("pest_factory_dt", tmp_path / "ws")

    # No start_datetime and the TDIS carries none -> a clear, actionable error.
    with pytest.raises(ValueError, match="start_datetime"):
        model.pest("calib")

    # Once the model's TDIS has a start date, the factory infers it.
    model.sim.tdis.start_date_time = "2019-06-01"
    cal = model.pest("calib")
    assert str(cal.start_datetime).startswith("2019-06-01")


def test_observe_accepts_named_series_targets_with_default_prefix(tmp_path):
    # observe() should wrap any high-level target set (not just HeadTargets) in
    # its matching spec, using the spec's own default prefix.
    model, _ = _build_two_cell_pest_forward_model("pest_obs_coerce", tmp_path / "ws")
    cal = model.pest("calib", start_datetime="2020-01-01")

    lake = LakeStageTargets(
        locations={"deep_lake": 0},
        values={"per": [0, 1], "deep_lake": [99.5, 100.5]},
        time_column="per",
    )
    spec = cal.observe(lake)  # bare target set: no spec, no prefix
    assert isinstance(spec, LakeStageObservationSpec)
    assert spec.targets is lake
    assert spec.prefix == "stage"  # the spec's own default, not "hds"
    assert cal._observation_specs == [spec]

    # Unsupported inputs are rejected without pointing at the deleted legacy API.
    with pytest.raises(TypeError) as excinfo:
        cal.observe(object())
    assert "legacy" not in str(excinfo.value).lower()
    assert "build_pst" not in str(excinfo.value)


def test_named_series_observation_postprocessor_is_wired(tmp_path):
    # The native build wires forward-run post-processors for head targets AND
    # named-series (lake/SFR/DRN) observations -- the latter was the gap.
    model, _ = _build_two_cell_pest_forward_model("pest_obs_wire", tmp_path / "ws")
    cal = model.pest("calib", start_datetime="2020-01-01")

    recorded: list[str] = []
    cal.pf = SimpleNamespace(
        add_py_function=lambda path, call, is_pre_cmd: recorded.append(call)
    )
    cal._prepared_observations = [
        {"prefix": "hds", "forward_run_config": {"mapping_csv": "m.csv", "output_csv": "h.csv"}},
        {
            "prefix": "lak_stage",
            "named_series_forward_run_config": {
                "kind": "lake_stage",
                "locations_file": "lak_stage_target_locations.gpkg",
                "output_csv": "lak_stage_simulated_lake_stage.csv",
            },
        },
    ]
    cal._attach_native_observation_postprocessors()

    assert any("_write_head_target_csv(" in call for call in recorded)
    series_call = next(call for call in recorded if "write_named_series_targets(" in call)
    assert "kind='lake_stage'" in series_call
    assert "locations_file='lak_stage_target_locations.gpkg'" in series_call
    assert "output_csv='lak_stage_simulated_lake_stage.csv'" in series_call

    # An unknown kind raises a clear error, not one pointing at the deleted build.
    cal._prepared_observations = [{"prefix": "mystery"}]
    with pytest.raises(NotImplementedError) as excinfo:
        cal._attach_native_observation_postprocessors()
    assert "build_pst" not in str(excinfo.value)


def test_parameterize_records_grid_anisotropy(tmp_path):
    # parameterize() carries the variogram anisotropy/bearing/nugget onto the spec.
    model, _ = _build_two_cell_pest_forward_model("pest_param_aniso", tmp_path / "ws")
    cal = model.pest("calib", start_datetime="2020-01-01")
    spec = cal.parameterize("k", style="grid", correlation=300.0, anisotropy=4.0, bearing=15.0)
    assert spec.anisotropy == 4.0
    assert spec.bearing == 15.0
    assert spec.nugget == 0.0  # isotropic-friendly default


def test_geostruct_for_routes_through_build_geostruct_with_anisotropy(tmp_path):
    # _geostruct_for delegates to the single build_geostruct builder and honours
    # anisotropy/bearing/nugget -- no longer the isotropic-only inline version.
    pytest.importorskip("pyemu")
    from myflopy.modflow.mf6.pest.native_parameters import NativeParameterSpec

    model, _ = _build_two_cell_pest_forward_model("pest_geostruct", tmp_path / "ws")
    cal = model.pest("calib", start_datetime="2020-01-01")

    spec = NativeParameterSpec(
        target="k",
        style="grid",
        correlation=500.0,
        anisotropy=5.0,
        bearing=30.0,
        nugget=0.1,
        transform="log",
    )
    geostruct = cal._geostruct_for(spec)
    vario = geostruct.variograms[0]
    assert vario.a == 500.0
    assert vario.anisotropy == 5.0
    assert vario.bearing == 30.0
    assert geostruct.nugget == 0.1
    assert geostruct.transform == "log"


def test_find_pest_runs_dedupes_parallel_master_copies(tmp_path):
    # Parallel PESTPP-IES clones the template into <run>_ies_master / _prior_master,
    # each carrying a copy of the metadata. find_pest_runs must list the run once
    # (with the masters as its execution kinds), not once per metadata copy.
    import json

    from myflopy.modflow.mf6.pest.runs import find_pest_runs

    pest = tmp_path / "pest"
    meta = {"project_name": "calib", "model_name": "m", "pst_file": "calib.pst"}
    for sub in ("calib", "calib_ies_master", "calib_prior_master"):
        d = pest / sub
        d.mkdir(parents=True)
        (d / "myflopy_pest_metadata.json").write_text(json.dumps(meta), encoding="utf-8")
        (d / "calib.pst").write_text("", encoding="utf-8")

    runs = find_pest_runs(pest)
    assert [r.name for r in runs] == ["calib"]  # one run, not three
    assert runs[0].kinds == ["ies", "prior"]    # masters become its execution kinds


# --- where calibrations live (ledger 107) ------------------------------------


def _write_stub_run(directory, name, model_name="m"):
    """The two files find_pest_runs keys on, for one build."""

    import json

    directory.mkdir(parents=True, exist_ok=True)
    (directory / "myflopy_pest_metadata.json").write_text(
        json.dumps({"project_name": name, "model_name": model_name,
                    "pst_file": f"{name}.pst"}),
        encoding="utf-8",
    )
    (directory / f"{name}.pst").write_text("", encoding="utf-8")
    return directory


def test_discovery_finds_runs_in_both_the_current_and_legacy_roots(tmp_path):
    """The default moved from `<ws>/pest` (inside the tree PstFrom copies) to
    `<ws>.pest` (a sibling). Calibrations already on disk are in the old place
    and must keep reviewing -- silently dropping them would look like the runs
    were never done."""

    from myflopy.modflow.mf6.pest.runs import find_model_pest_runs, pest_run_roots

    workspace = tmp_path / "model"
    workspace.mkdir()
    current, legacy = pest_run_roots(workspace)
    assert current == tmp_path / "model.pest"
    assert legacy == workspace / "pest"

    _write_stub_run(current / "new_run", "new_run")
    _write_stub_run(legacy / "old_run", "old_run")

    assert [run.name for run in find_model_pest_runs(workspace)] == [
        "new_run", "old_run"
    ]
    assert [run.name for run in find_model_pest_runs(workspace, model_name="other")] == []


def test_a_legacy_run_copied_into_a_template_is_not_listed_twice(tmp_path):
    """A template is a COPY of the model workspace, so a calibration that lived
    inside it (the pre-2026-07-29 default) rides along into every template built
    afterwards. Listing that copy as a run of its own shows the same
    calibration twice, and ``review()`` on the copy opens a directory nothing
    ever ran in."""

    import shutil

    from myflopy.modflow.mf6.pest.runs import find_model_pest_runs

    workspace = tmp_path / "model"
    workspace.mkdir()
    _write_stub_run(workspace / "pest" / "old_run", "old_run")

    template = tmp_path / "model.pest" / "new_run"
    shutil.copytree(workspace, template)          # what PstFrom does
    _write_stub_run(template, "new_run")
    assert (template / "pest" / "old_run").exists(), "fixture must carry the copy"

    assert [run.name for run in find_model_pest_runs(workspace)] == [
        "new_run", "old_run"
    ]


def test_a_template_inside_the_copied_workspace_still_refuses_to_nest(tmp_path):
    """`_clear_stale_template` had no test at all, and it is what stands between
    a rebuild and the 17 GB recursion of ledger 107.

    The default no longer lands inside the copied tree, so this drives the case
    that remains: an explicit `workspace=` under the model directory, with
    another calibration already there. Deleting that one is not an option -- its
    IES master holds finished results -- so the build must refuse.
    """

    from myflopy.modflow.mf6.pest.project import PestProject

    model_workspace = tmp_path / "model"
    model_workspace.mkdir()
    pest_root = model_workspace / "pest"

    project = object.__new__(PestProject)
    project.name = "second"
    project.original_workspace = model_workspace
    project.template_workspace = pest_root / "second"

    # Alone in the directory: the stale template AND the now-empty parent go, so
    # the copy sees no `pest/` at all -- an EMPTY one is enough to recurse.
    project.template_workspace.mkdir(parents=True)
    (project.template_workspace / "leftover.txt").write_text("x", encoding="utf-8")
    project._clear_stale_template()
    assert not pest_root.exists()

    # With somebody else's finished run there, refuse rather than delete it.
    _write_stub_run(pest_root / "first", "first")
    project.template_workspace.mkdir(parents=True)
    with pytest.raises(RuntimeError, match="still holds 1 other PEST run"):
        project._clear_stale_template()
    assert (pest_root / "first").exists(), "another run's results were deleted"

    # The default location is outside the copied tree, so it never gets here.
    project.template_workspace = tmp_path / "model.pest" / "second"
    project.template_workspace.mkdir(parents=True)
    project._clear_stale_template()
    assert (pest_root / "first").exists()


def test_two_named_calibrations_on_one_model_coexist_at_the_default_location():
    """The capability ledger 107 gave up to stop the disk filling: with the
    template inside the copied tree, a second differently-named run at the
    default location could not be built at all, so `model.pest_runs` with
    several runs required explicit workspaces.

    Also pins the constraint that makes discovery work: the template must stay
    OUTSIDE the model workspace (or the copy recurses), and any IES master must
    stay a SIBLING of the template (or `find_pest_runs` reports a finished run
    as "built (not run)").
    """

    pytest.importorskip("pyemu")
    pytest.importorskip("flopy")

    workspace = _project_temp_dir("two_calibrations")
    model, vor = _build_two_cell_pest_forward_model("two_cal", workspace / "model")
    Recharge(model=model, vor=vor, rch_dict={0: [[(0, 0), 0.001], [(0, 1), 0.001]]})
    assert model.run_simulation()[0] is True

    obs_path = _write_gpkg(
        workspace / "two_cal_obs.gpkg",
        gpd.GeoDataFrame(
            {"name": ["OBS_A"], "layer": [0], "weight": [1.0]},
            geometry=[Point(0.5, 0.5)],
            crs=model.vor.crs,
        ),
    )
    targets = HeadTargets(
        locations=obs_path,
        values=pd.DataFrame({"per": [0], "OBS_A": [10.0]}),
        time_column="per",
    )

    for name in ("first_calib", "second_calib"):
        cal = model.pest(name, start_datetime="2024-01-01")
        cal.parameterize("k", style="constant", bounds=(0.2, 5.0), physical=(1e-3, 100.0))
        cal.observe(targets)
        cal.build(f"{name}.pst", noptmax=0)
        assert model.workspace not in cal.template_workspace.parents, (
            "the template is inside the tree PstFrom copies; it will recurse"
        )
        # Where run_ies puts its master, spelled the way project.py spells it.
        assert (cal.template_workspace.parent / f"{name}_ies_master").parent == (
            cal.template_workspace.parent
        )

    # Nothing was copied into the model workspace, at any depth.
    assert not list(Path(model.workspace).rglob("myflopy_pest_metadata.json"))

    discovered = model.pest_runs
    assert [run.name for run in discovered] == ["first_calib", "second_calib"]
    assert all(run.model_name == model.name for run in discovered)


def test_native_pstfrom_parameterize_build_and_forward_run_end_to_end():
    pytest.importorskip("pyemu")
    pytest.importorskip("flopy")
    from pyemu.pst.pst_utils import write_to_template

    workspace = _project_temp_dir("native_pstfrom")
    model, vor = _build_two_cell_pest_forward_model("native_pstfrom", workspace / "model")
    Recharge(model=model, vor=vor, rch_dict={0: [[(0, 0), 0.001], [(0, 1), 0.001]]})
    success, _ = model.run_simulation()
    assert success is True

    obs_path = _write_gpkg(
        workspace / "native_obs.gpkg",
        gpd.GeoDataFrame(
            {"name": ["OBS_A", "OBS_B"], "layer": [0, 0], "weight": [1.0, 1.0]},
            geometry=[Point(0.5, 0.5), Point(1.5, 0.5)],
            crs=model.vor.crs,
        ),
    )
    targets = HeadTargets(
        locations=obs_path,
        values=pd.DataFrame({"per": [0], "OBS_A": [10.0], "OBS_B": [8.75]}),
        time_column="per",
    )
    fore_path = _write_gpkg(
        workspace / "native_fore.gpkg",
        gpd.GeoDataFrame(
            {"name": ["PRED"], "layer": [0], "weight": [1.0]},
            geometry=[Point(1.5, 0.5)],
            crs=model.vor.crs,
        ),
    )
    forecast = HeadTargets(
        locations=fore_path,
        values=pd.DataFrame({"per": [0], "PRED": [9.0]}),
        time_column="per",
    )

    cal = model.pest(
        "native_demo",
        workspace=workspace / "template",
        start_datetime="2024-01-01",
    )
    cal.parameterize("k", style="constant", bounds=(0.2, 5.0), physical=(1e-3, 100.0))
    cal.parameterize("recharge", style="constant", bounds=(0.5, 1.5), physical=(0.0, 0.01))
    cal.observe(targets)
    cal.forecast(forecast)

    settings_before = cal.settings()
    assert settings_before.built is False
    assert len(settings_before.parameters) == 2

    pst = cal.build("native_demo.pst", noptmax=0)

    template = workspace / "template"
    forward_run_text = (template / "forward_run.py").read_text(encoding="utf-8")
    # The whole point of the native path: pyEMU's own apply is kept, not stripped.
    assert "apply_list_and_array_pars" in forward_run_text
    call_lines = [
        line for line in forward_run_text.splitlines()
        if "_write_head_target_csv(" in line and "def " not in line
    ]
    assert any("hds_simulated_heads.csv" in line for line in call_lines)
    assert any("fore1_simulated_heads.csv" in line for line in call_lines)
    assert pst.pestpp_options.get("forecasts")

    settings_after = cal.settings()
    assert settings_after.built is True
    assert settings_after.npar == 2
    assert settings_after.nnz_obs == 2
    assert settings_after.n_forecasts == 1

    # A unit-multiplier forward run must reproduce the baseline model inputs.
    subprocess.run(
        [sys.executable, "forward_run.py"],
        cwd=template, check=True, capture_output=True, text=True,
    )
    assert (template / "mult2model_info.csv").exists()
    k_file = template / f"{model.name}.npf_k.txt"
    assert np.allclose(np.loadtxt(k_file), [5.0, 5.0])
    assert (template / "hds_simulated_heads.csv").exists()
    assert (template / "fore1_simulated_heads.csv").exists()

    # Now prove the multipliers actually scale the model inputs.
    par = pst.parameter_data
    par.loc[par.index[par["pargp"] == "k"], "parval1"] = 2.0
    par.loc[par.index[par["pargp"] == "recharge"], "parval1"] = 0.5
    for tpl, inp in zip(pst.template_files, pst.input_files, strict=False):
        write_to_template(par["parval1"], str(template / tpl), str(template / inp))
    subprocess.run(
        [sys.executable, "forward_run.py"],
        cwd=template, check=True, capture_output=True, text=True,
    )
    assert np.allclose(np.loadtxt(k_file), [10.0, 10.0])
    rch = pd.read_csv(
        template / f"{model.name}.rch_stress_period_data_1.txt", sep=r"\s+", header=None
    )
    assert np.allclose(rch.iloc[:, 2].to_numpy(), [0.0005, 0.0005])


def test_named_series_drn_observation_forward_run_end_to_end():
    # End-to-end proof that named-series (lake/SFR/DRN) observations regenerate
    # their simulated values in the native forward run. The post-processor
    # reloads the model (mf.load_mf6_run) and rebuilds the target series -- a path
    # the head-target forward run never exercises. DRN is the lightest runnable
    # named-series fixture.
    pytest.importorskip("pyemu")
    pytest.importorskip("flopy")
    import flopy

    workspace = _project_temp_dir("named_series_drn")
    model, vor = _build_two_cell_pest_forward_model("ns_drn", workspace / "model")
    # A drain on the downstream cell seeps groundwater -> a DRN-flow target.
    flopy.mf6.ModflowGwfdrn(
        model.gwf,
        stress_period_data={0: [[(0, 1), 8.0, 1.0]]},
        save_flows=True,
        pname="drn",
    )
    success, _ = model.run_simulation()
    assert success is True
    assert not model.bud("drn").df.empty  # the drain is seeping

    drn_targets = DrnFlowTargets(
        locations={"seep": [1]},
        values=pd.DataFrame({"per": [0], "seep": [0.5]}),
        time_column="per",
        value_column="flow",
    )

    cal = model.pest("ns_drn_demo", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.2, 5.0), physical=(1e-3, 100.0))
    cal.observe(drn_targets)  # a bare named-series target set
    template = cal.template_workspace
    cal.build("ns_drn_demo.pst", noptmax=0)

    # The named-series post-processor is wired into the forward run (not the
    # deleted legacy build_pst path).
    forward_run_text = (template / "forward_run.py").read_text(encoding="utf-8")
    assert "def write_named_series_targets(" in forward_run_text
    assert "kind='drn_flow'" in forward_run_text

    # Run it: MF6 + the post-processor that reloads the model and regenerates the
    # DRN seepage series pyEMU reads back as observations.
    subprocess.run(
        [sys.executable, "forward_run.py"],
        cwd=template, check=True, capture_output=True, text=True,
    )
    sim_csv = template / "drn_flow_simulated_drn_flow.csv"
    assert sim_csv.exists()
    simulated = pd.read_csv(sim_csv)
    assert "seep" in simulated.columns
    assert len(simulated) >= 1


def test_prior_monte_carlo_and_conflict_end_to_end():
    pytest.importorskip("pyemu")
    pytest.importorskip("flopy")
    import plotly.graph_objects as go
    from matplotlib.figure import Figure as MplFigure

    workspace = _project_temp_dir("prior_mc")
    model, vor = _build_two_cell_pest_forward_model("prior_mc", workspace / "model")
    Recharge(model=model, vor=vor, rch_dict={0: [[(0, 0), 0.001], [(0, 1), 0.001]]})
    success, _ = model.run_simulation()
    assert success is True

    obs_path = _write_gpkg(
        workspace / "prior_obs.gpkg",
        gpd.GeoDataFrame(
            {"name": ["OBS_A", "OBS_B"], "layer": [0, 0], "weight": [1.0, 1.0]},
            geometry=[Point(0.5, 0.5), Point(1.5, 0.5)],
            crs=model.vor.crs,
        ),
    )
    targets = HeadTargets(
        locations=obs_path,
        values=pd.DataFrame({"per": [0], "OBS_A": [9.5], "OBS_B": [8.75]}),
        time_column="per",
    )

    cal = model.pest(
        "prior_demo",
        workspace=workspace / "template",
        start_datetime="2024-01-01",
    )
    cal.parameterize("k", style="constant", bounds=(0.2, 5.0), physical=(1e-3, 100.0))
    cal.observe(targets)
    cal.build("prior_demo.pst", noptmax=0)

    # Prior Monte Carlo: run the prior ensemble once (NOPTMAX=-1), no iterations.
    prior = cal.prior(reals=6)

    assert prior.iterations == [0]
    assert prior.prior._df.shape[0] >= 6

    conflict = prior.conflict()
    assert {"measured", "prior_lo", "prior_hi", "in_conflict"}.issubset(conflict.columns)
    assert len(conflict) == 2

    assert isinstance(prior.plot_prior_vs_obs(), go.Figure)
    assert isinstance(prior.plot_prior_vs_obs(backend="matplotlib"), MplFigure)
    assert isinstance(prior.plot_conflict(), go.Figure)
    assert isinstance(prior.plot_conflict(backend="matplotlib"), MplFigure)


def test_run_ies_end_to_end_and_assess_with_ies_results():
    pytest.importorskip("pyemu")
    pytest.importorskip("flopy")
    import plotly.graph_objects as go

    from myflopy.modflow.mf6.pest.ies import IesForecast, IesResults

    workspace = _project_temp_dir("run_ies")
    model, vor = _build_two_cell_pest_forward_model("run_ies", workspace / "model")
    Recharge(model=model, vor=vor, rch_dict={0: [[(0, 0), 0.001], [(0, 1), 0.001]]})
    success, _ = model.run_simulation()
    assert success is True

    obs_path = _write_gpkg(
        workspace / "ies_obs.gpkg",
        gpd.GeoDataFrame(
            {"name": ["OBS_A", "OBS_B"], "layer": [0, 0], "weight": [1.0, 1.0]},
            geometry=[Point(0.5, 0.5), Point(1.5, 0.5)],
            crs=model.vor.crs,
        ),
    )
    targets = HeadTargets(
        locations=obs_path,
        values=pd.DataFrame({"per": [0], "OBS_A": [9.5], "OBS_B": [8.75]}),
        time_column="per",
    )
    fore_path = _write_gpkg(
        workspace / "ies_fore.gpkg",
        gpd.GeoDataFrame(
            {"name": ["PRED"], "layer": [0], "weight": [1.0]},
            geometry=[Point(1.5, 0.5)],
            crs=model.vor.crs,
        ),
    )
    forecast = HeadTargets(
        locations=fore_path,
        values=pd.DataFrame({"per": [0], "PRED": [9.0]}),
        time_column="per",
    )

    # No explicit workspace: it defaults to <model workspace>.pest/<name>, so the
    # run is auto-discoverable via model.pest_runs (the workflow integration).
    cal = model.pest("run_ies_demo", start_datetime="2024-01-01")
    assert cal.template_workspace == (
        model.workspace.parent / f"{model.workspace.name}.pest" / "run_ies_demo"
    )
    cal.parameterize("k", style="constant", bounds=(0.2, 5.0), physical=(1e-3, 100.0))
    cal.parameterize("recharge", style="grid", bounds=(0.5, 1.5), physical=(0.0, 1e-2))
    cal.observe(targets)
    cal.forecast(forecast)
    cal.build("run_ies_demo.pst", noptmax=0)

    # iterations=1 + lambda_scale_fac=1.0 cut pestpp forward runs 82 -> 20
    # while still exercising prior + posterior; workers=4 also covers the
    # parallel PESTPP-IES agent path (previously untested).
    # ies_lambda_mults=1.0 + lambda_scale_fac=1.0 tests a single upgrade
    # candidate per iteration: 82 -> 12 pestpp forward runs. Deterministic
    # (fixed demo seed + pestpp's fixed default ies seed).
    ies = cal.run_ies(
        reals=6, iterations=1, noise=True, workers=4,
        ies_lambda_mults=1.0, lambda_scale_fac=1.0,
    )

    assert isinstance(ies, IesResults)
    assert ies.iterations[0] == 0
    assert ies.posterior_iteration >= 1
    # prior and posterior observation ensembles share columns; >= reals rows
    assert ies.prior._df.shape[0] >= 6
    assert ies.posterior._df.shape[1] == ies.prior._df.shape[1]
    assert ies.noise is not None

    # forecast access + uncertainty summary
    assert len(ies.forecast_names) == 1
    fc = ies.forecast("pred")
    assert isinstance(fc, IesForecast)
    summary = fc.summary()
    assert {"prior_std", "posterior_std", "truth"}.issubset(summary.index)
    assert not ies.forecasts().empty

    # headline plots are Plotly figures by default, matplotlib on request
    from matplotlib.figure import Figure as MplFigure

    assert isinstance(ies.plot_phi(), go.Figure)
    assert isinstance(ies.plot_vs_obs(), go.Figure)
    assert isinstance(fc.plot(), go.Figure)
    assert isinstance(ies.plot_phi(backend="matplotlib"), MplFigure)
    assert isinstance(ies.plot_vs_obs(backend="matplotlib"), MplFigure)
    assert isinstance(fc.plot(backend="matplotlib"), MplFigure)

    # phi & weight diagnostics (both backends)
    assert isinstance(ies.plot_phi_distribution(), go.Figure)
    assert isinstance(ies.plot_phi_distribution(backend="matplotlib"), MplFigure)
    bounds_table = ies.parameters_at_bounds()
    assert "pct_at_bound" in bounds_table.columns and not bounds_table.empty
    assert isinstance(ies.plot_parameters_at_bounds(), go.Figure)
    assert isinstance(ies.plot_parameters_at_bounds(backend="matplotlib"), MplFigure)
    contributions = ies.phi_contributions()
    assert not contributions.empty
    assert isinstance(ies.plot_phi_contributions(), go.Figure)
    assert isinstance(ies.plot_phi_contributions(kind="pie"), go.Figure)
    assert isinstance(ies.plot_phi_contributions(kind="pie", backend="matplotlib"), MplFigure)

    # base realization is the recommended single parameter set
    assert ies.best() == "base"

    # one-shot HTML report
    report = ies.report(workspace / "ies_report.html")
    assert report.exists()

    # settings snapshot reflects the run configuration
    assert ies.settings.num_reals == 6
    assert ies.settings.n_forecasts == 1

    # Workflow integration: the completed run is discoverable from the model and
    # reopens to the same IesResults review.
    discovered = model.pest_runs
    assert [run.name for run in discovered] == ["run_ies_demo"]
    handle = discovered[0]
    assert handle.model_name == model.name
    reopened = handle.review(model=model)
    assert isinstance(reopened, IesResults)
    assert reopened.iterations == ies.iterations
    assert isinstance(reopened.plot_phi(), go.Figure)


def test_ies_capture_field_and_spatial_maps_end_to_end():
    pytest.importorskip("pyemu")
    pytest.importorskip("flopy")
    import plotly.graph_objects as go

    workspace = _project_temp_dir("ies_capture")
    model, vor = _build_two_cell_pest_forward_model("ies_capture", workspace / "model")
    Recharge(model=model, vor=vor, rch_dict={0: [[(0, 0), 0.001], [(0, 1), 0.001]]})
    success, _ = model.run_simulation()
    assert success is True

    obs_path = _write_gpkg(
        workspace / "cap_obs.gpkg",
        gpd.GeoDataFrame(
            {"name": ["OBS_A", "OBS_B"], "layer": [0, 0], "weight": [1.0, 1.0]},
            geometry=[Point(0.5, 0.5), Point(1.5, 0.5)],
            crs=model.vor.crs,
        ),
    )
    targets = HeadTargets(
        locations=obs_path,
        values=pd.DataFrame({"per": [0], "OBS_A": [9.5], "OBS_B": [8.75]}),
        time_column="per",
    )

    cal = model.pest(
        "cap_demo",
        workspace=workspace / "template",
        start_datetime="2024-01-01",
    )
    cal.parameterize("k", style="constant", bounds=(0.2, 5.0), physical=(1e-3, 100.0), capture=True)
    cal.observe(targets)
    pst = cal.build("cap_demo.pst", noptmax=0)

    # The resolved K field is captured as zero-weight observations.
    captured = pst.observation_data[pst.observation_data.index.str.contains("oname:kfield", regex=False)]
    assert len(captured) == int(model.vor.ncpl)
    assert bool((captured["weight"] == 0).all())

    # workers=4 keeps this off the parallel suite's critical path; the
    # serial pestpp runner path stays covered by the prior-MC test. The
    # trimmed lambda sweep cuts pestpp forward runs 44 -> 20.
    ies = cal.run_ies(
        reals=6, iterations=1, workers=4, ies_lambda_mults=1.0, lambda_scale_fac=1.0
    )

    assert any(info["target"] == "k" for info in ies.capture_fields)
    field = ies.field("k")
    assert len(field) == int(model.vor.ncpl)
    assert {"prior_mean", "posterior_mean", "posterior_std", "change"}.issubset(field.columns)

    assert isinstance(ies.plot_field("k", stat="mean", which="posterior"), go.Figure)
    assert isinstance(ies.plot_field("k", stat="change"), go.Figure)
    assert isinstance(ies.plot_field("k", stat="std"), go.Figure)

    from matplotlib.figure import Figure as MplFigure

    assert isinstance(ies.plot_field("k", stat="mean", backend="matplotlib"), MplFigure)
    assert isinstance(ies.plot_field("k", stat="change", backend="matplotlib"), MplFigure)
    assert isinstance(ies.plot_field("k", stat="std", backend="matplotlib"), MplFigure)
    assert isinstance(ies.plot_field("k", stat="reduction", backend="matplotlib"), MplFigure)

    # -- the color policy and hover survive the real pyEMU/pestpp round trip ----
    # (the arithmetic is pinned fast in test_colorscale_policy.py; what is only
    # provable here is that a REAL captured field reaches the trace intact.)
    change = ies.plot_field("k", stat="change").data[0]
    # a diverging ratio ships as explicit stops, never the name 'RdBu', which
    # renders mirrored between the two backends
    assert not isinstance(change.colorscale, str)
    assert change.zmin == -change.zmax
    # the RAW ratio is in the hover even though log10 is what gets plotted
    assert change.customdata is not None
    assert "posterior / prior" in change.hovertemplate
    # run context rides the title, because the hover footer cannot carry it
    assert "realizations" in ies.plot_field("k", stat="change").layout.title.text

    mean = ies.plot_field("k", stat="mean").data[0]
    assert tuple(mean.colorbar.ticktext)  # decades relabeled out of log10 units

    # an explicit kwarg overrides the policy rather than raising "multiple values"
    override = ies.plot_field("k", stat="change", colorscale="earth").data[0]
    assert override.colorscale[0][1] != change.colorscale[0][1]

    # -- 6.4B: uncertainty reduction, mosaics, residuals on the same real run --
    assert "reduction" in field.columns
    reduction = ies.plot_field("k", stat="reduction").data[0]
    # 1 - post_sd/prior_sd is derived from the SAME two columns the frame shows
    assert reduction.customdata is not None
    assert "realizations" in ies.plot_field("k", stat="reduction").layout.title.text

    mosaic = ies.plot_field_mosaic("k", stat="mean")
    # both panels pooled onto one color axis -- that shared scale is the whole
    # point of a mosaic, and its ticks must still read in real units
    assert len([trace for trace in mosaic.data if getattr(trace, "z", None) is not None]) == 2
    assert all(trace.coloraxis == "coloraxis" for trace in mosaic.data
               if getattr(trace, "z", None) is not None)
    assert tuple(mosaic.layout.coloraxis.colorbar.ticktext)

    # the head targets this run was built with are recoverable from disk alone
    residuals = ies.obs_residuals()
    assert not residuals.empty
    assert {"OBS_A".lower(), "OBS_B".lower()} == set(residuals["location"])
    assert residuals["x"].notna().all()  # joined back to the saved gpkg
    residual_map = ies.plot_obs_residuals()
    markers = [trace for trace in residual_map.data if trace.type == "scattermap"]
    assert len(markers) == 1 and len(markers[0].lon) == len(residuals)
    assert isinstance(ies.plot_obs_residuals(backend="matplotlib"), MplFigure)


def test_canonical_calibration_demo_builds_native_pst_with_multilayer_k(tmp_path):
    """The canonical calibration demo wires the canonical valley model into a
    native PstFrom build, including multi-layer DISV K resolved per layer."""

    from myflopy.modflow.mf6.canonical_calibration import (
        build_canonical_calibration_demo,
    )

    demo = build_canonical_calibration_demo(
        tmp_path / "model", config=CanonicalModelConfig.testing(), n_head_wells=8
    )
    assert demo.model.gwf.modelgrid.nlay == 4
    assert demo.head_targets.locations_gdf.shape[0] == 8
    assert demo.forecast_targets.locations_gdf.shape[0] == 1
    assert demo.start_k_factor == 3.0

    cal = demo.model.pest(
        "cc",
        workspace=tmp_path / "template",
        start_datetime="2024-01-01",
    )
    kspec = cal.parameterize("k", style="constant", bounds=(0.1, 10.0), physical=(0.1, 300.0))
    cal.parameterize("recharge", style="constant", bounds=(0.3, 3.0), physical=(0.0, 1e-2))
    cal.observe(demo.head_targets)
    cal.forecast(demo.forecast_targets)
    pst = cal.build("cc.pst", noptmax=0)

    # Two adjustable parameters: one constant K multiplier + one constant recharge.
    assert pst.npar_adj == 2
    # Eight wells x every stress period of nonzero-weight head observations.
    assert pst.nnz_obs == 8 * demo.model.nper
    # The multi-layer DISV K array resolves to one external file per layer
    # (npf_k_layer1..4), and never grabs the k33 files.
    assert len(kspec.resolved_files) == demo.model.gwf.modelgrid.nlay
    assert all("npf_k_layer" in name and "k33" not in name for name in kspec.resolved_files)


def test_grid_k_parameterization_on_voronoi_with_capture(tmp_path):
    """Native grid K builds one geostatistical parameter per Voronoi cell, per
    selected layer, and capture records a non-colliding per-layer K field."""

    from myflopy.modflow.mf6.canonical_calibration import (
        build_canonical_calibration_demo,
    )

    demo = build_canonical_calibration_demo(
        tmp_path / "model",
        config=CanonicalModelConfig.testing(),
        n_head_wells=6,
        start_k_constant=10.0,
        start_k_layers=(0, 1),
    )
    cal = demo.model.pest(
        "gk", workspace=tmp_path / "tmpl",
        start_datetime="2024-01-01",
    )
    kspec = cal.parameterize("k", style="grid", layers=[0, 1], correlation=600.0,
                             bounds=(0.02, 50.0), physical=(0.001, 300.0), capture=True)
    cal.parameterize("recharge", style="constant", bounds=(0.3, 3.0), physical=(0.0, 1e-2))
    cal.observe(demo.head_targets)
    pst = cal.build("gk.pst", noptmax=0)

    ncpl = int(demo.model.vor.ncpl)
    # Two unconfined layer files -> one grid parameter per cell on each, + recharge.
    assert len(kspec.resolved_files) == 2
    assert all("npf_k_layer" in name for name in kspec.resolved_files)
    assert pst.npar_adj == 2 * ncpl + 1
    assert pst.parameter_data["pargp"].nunique() == 3  # kl1, kl2, recharge

    # Capture obs: per-cell K field per layer, distinct (non-colliding) names,
    # all zero-weight.
    capture = pst.observation_data.index.str.contains("kfield")
    assert int(capture.sum()) == 2 * ncpl
    assert (pst.observation_data.loc[capture, "weight"] == 0.0).all()


def test_pilot_point_k_parameterization_on_voronoi(tmp_path):
    """Native pilot-point K places points, registers one parameter each, and the
    IDW forward run applies them and runs MF6 -- pyEMU's pilot points do not work
    on unstructured grids, so the facade uses inverse-distance weighting."""

    import subprocess

    from myflopy.modflow.mf6.canonical_calibration import (
        build_canonical_calibration_demo,
    )

    demo = build_canonical_calibration_demo(
        tmp_path / "model",
        config=CanonicalModelConfig.testing(),
        n_head_wells=6,
        start_k_constant=10.0,
        start_k_layers=(0, 1),
    )
    cal = demo.model.pest(
        "ppk", workspace=tmp_path / "tmpl",
        start_datetime="2024-01-01",
    )
    # pp_space=4 on the 2,100 m testing-profile domain -> ~5x5 net (>= 10 pp);
    # the old pp_space=8 was sized for the 5,000 m validation domain.
    cal.parameterize("k", style="pilotpoints", pp_space=4, layers=[0, 1],
                     bounds=(0.05, 20.0), physical=(0.001, 300.0), capture=True)
    cal.parameterize("recharge", style="constant", bounds=(0.3, 3.0), physical=(0.0, 1e-2))
    cal.observe(demo.head_targets)
    pst = cal.build("ppk.pst", noptmax=0)

    # One parameter per pilot point on each of the two layers, plus recharge.
    assert pst.parameter_data["pargp"].nunique() == 3  # kl0, kl1, recharge
    n_pp = int((pst.parameter_data["pargp"] == "kl0").sum())
    assert n_pp >= 10
    assert int((pst.parameter_data["pargp"] == "kl1").sum()) == n_pp
    assert pst.npar_adj == 2 * n_pp + 1
    # Per-cell captured K field on both layers.
    assert int(pst.observation_data.index.str.contains("kfield").sum()) == 2 * int(demo.model.vor.ncpl)

    # The IDW forward run applies pilot points and runs MF6.
    result = subprocess.run(
        [sys.executable, "forward_run.py"], cwd=cal.template_workspace,
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout[-2000:] + "\n" + result.stderr[-2000:]
    assert (cal.template_workspace / "hds_simulated_heads.csv").exists()


# --- pilot points scale the array they NAME (2026-07-30) ---------------------


@pytest.mark.slow
def test_pilot_points_scale_their_own_target_not_the_k_field(tmp_path):
    """`add_pilot_point_parameter` read `gwf.npf.k` as the interpolation base for
    EVERY target, then wrote the result to the target's own file. So
    `parameterize("k33", style="pilotpoints")` replaced K33 with horizontal K --
    measured 30x on the canonical model (2.86 -> 85.9), destroying vertical
    anisotropy while the forward run exited 0 and MODFLOW reported normal
    termination.

    The earlier guard (ledger 110) only refused `recipe.model != "flow"`, so it
    caught transport targets and missed k33 entirely -- and every future flow
    array target would have inherited the same bug.
    """

    pytest.importorskip("pyemu")
    from myflopy.modflow.mf6.canonical_calibration import (
        build_canonical_calibration_demo,
    )

    demo = build_canonical_calibration_demo(
        tmp_path / "model", config=CanonicalModelConfig.testing(), n_head_wells=4
    )
    model = demo.model
    k = np.asarray(model.gwf.npf.k.get_data(), dtype=float)
    k33 = np.asarray(model.gwf.npf.k33.get_data(), dtype=float)
    assert not np.allclose(k[0], k33[0]), (
        "fixture is useless unless K and K33 differ"
    )

    cal = model.pest("k33pp", start_datetime="2024-01-01")
    cal.parameterize("k33", style="pilotpoints", pp_space=6,
                     bounds=(0.5, 2.0), physical=(1e-6, 1e3))
    cal.parameterize("recharge", style="constant", bounds=(0.5, 1.5),
                     physical=(0.0, 1e-2))
    cal.observe(demo.head_targets)
    cal.build("k33pp.pst", noptmax=0)

    result = subprocess.run(
        [sys.executable, "forward_run.py"], cwd=cal.template_workspace,
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout[-1500:] + result.stderr[-1500:]

    written = np.loadtxt(cal.template_workspace / f"{model.name}.npf_k33_layer1.txt")
    assert np.allclose(written, k33[0], rtol=1e-6), (
        "unit multipliers must reproduce the target's OWN array"
    )
    assert not np.allclose(written, k[0], rtol=1e-3), (
        "K33 was overwritten with the horizontal K field"
    )


def test_pilot_points_refuse_an_array_target_with_no_declared_base_array():
    """The replacement guard: pilot points multiply a base array read off the
    model, so an ARRAY recipe that does not name one must refuse rather than
    fall back to whatever package the code happened to hardcode.

    Driven through a synthetic recipe because all three shipped array targets
    declare package/variable -- which is the point: the guard now describes the
    real precondition instead of naming one model kind.
    """

    from myflopy.modflow.mf6.pest import native_parameters as np_mod

    incomplete = np_mod._Recipe("mystery", "array", "{model}.mystery.txt")
    original = dict(np_mod._RECIPES)
    np_mod._RECIPES["mystery"] = incomplete
    try:
        with pytest.raises(NotImplementedError, match="which model array"):
            np_mod.NativeParameterSpec(target="mystery", style="pilotpoints")
        # ...but the same recipe is fine on the styles that need no base array.
        assert np_mod.NativeParameterSpec(target="mystery", style="constant").style == "constant"
    finally:
        np_mod._RECIPES.clear()
        np_mod._RECIPES.update(original)

    # The array targets that DO declare one are accepted -- including the
    # transport target the previous guard refused for a reason that the
    # recipe-resolved base array has now removed.
    for target in ("k", "k33", "porosity"):
        assert np_mod.NativeParameterSpec(target=target, style="pilotpoints").style == "pilotpoints"


@pytest.mark.slow
def test_a_pilot_points_only_calibration_runs(tmp_path):
    """pyEMU inserts `apply_list_and_array_pars` into every forward run but only
    writes the `mult2model_info.csv` it reads when `pf.add_parameters` was
    called. Pilot points register through template files instead, so a
    calibration whose ONLY parameter is pilot points died with
    `FileNotFoundError: mult2model_info.csv` from inside pyEMU.

    It stayed hidden because every example pairs pilot points with a second
    parameter -- but "calibrate K with pilot points" is a perfectly ordinary
    thing to ask for.
    """

    pytest.importorskip("pyemu")
    from myflopy.modflow.mf6.canonical_calibration import (
        build_canonical_calibration_demo,
    )

    demo = build_canonical_calibration_demo(
        tmp_path / "model", config=CanonicalModelConfig.testing(), n_head_wells=4
    )
    model = demo.model
    k = np.asarray(model.gwf.npf.k.get_data(), dtype=float)

    cal = model.pest("pponly", start_datetime="2024-01-01")
    cal.parameterize("k", style="pilotpoints", pp_space=6,
                     bounds=(0.5, 2.0), physical=(1e-3, 300.0))
    cal.observe(demo.head_targets)
    pst = cal.build("pponly.pst", noptmax=0)

    assert pst.npar_adj > 0
    result = subprocess.run(
        [sys.executable, "forward_run.py"], cwd=cal.template_workspace,
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout[-1500:] + result.stderr[-1500:]
    written = np.loadtxt(cal.template_workspace / f"{model.name}.npf_k_layer1.txt")
    assert np.allclose(written, k[0], rtol=1e-6)


def test_build_forward_run_command_is_pestpp_compatible(tmp_path):
    """The PST model command must survive PEST++'s run manager on both OSes.

    POSIX PEST++ neither shell-parses quotes nor accepts absolute command
    paths (it mangles the leading '/' and execv fails), so the command must
    be a ./wrapper script; Windows uses the quoted interpreter directly.
    """

    from myflopy.modflow.mf6.pest.project import build_forward_run_command

    command = build_forward_run_command(tmp_path, "/opt/env/bin/python")
    if os.name == "nt":
        assert command == '"/opt/env/bin/python" forward_run.py'
    else:
        assert command == "./run_forward.sh"
        wrapper = tmp_path / "run_forward.sh"
        content = wrapper.read_text()
        assert content.startswith("#!/bin/sh")
        assert 'exec "/opt/env/bin/python" forward_run.py' in content
        assert wrapper.stat().st_mode & 0o111, "wrapper must be executable"
        assert not command.startswith("/"), "absolute commands break POSIX PEST++"


# --- observation residual maps (6.4B) ----------------------------------------
#
# These drive IesResults straight off a hand-built run directory: the residual
# path reads the control file, ONE posterior ensemble, and the location
# snapshots, so it does not need pestpp-ies to have run. That keeps the join
# logic -- which is where the real risk is -- on the fast suite, with the slow
# end-to-end test proving the same code against a genuine run.


def _stub_ies_results(tmp_path, *, heads=True, zones=False, forecast=False,
                      conc=False, unmeasured_periods=0, zone_cells="0,1",
                      zone_residual=-6.0, head_residuals=(2.0, -1.0),
                      conc_residual=0.02):
    """A minimal completed-run stand-in for the residual-map surface.

    Writes the same snapshot files a real build leaves beside the control file,
    so the join logic under test is the real one. ``unmeasured_periods`` adds
    the rows pyEMU creates for output times that were never measured — weight
    1.0, ``obsval`` equal to the model's own simulated value — which is what a
    real run looks like when targets cover only some stress periods.
    """

    from myflopy.modflow.mf6.pest.ies import IesResults

    workspace = tmp_path / "master"
    workspace.mkdir(parents=True, exist_ok=True)
    vor = _two_cell_vor_clockwise()
    observation_sets, names, values = [], [], {}

    if heads:
        frame = gpd.GeoDataFrame(
            {"name": ["OBS_A", "OBS_B"], "layer": [0, 1],
             "group": ["g", "g"], "weight": [1.0, 1.0]},
            geometry=[Point(0.5, 0.5), Point(1.5, 0.5)],
            crs="EPSG:2927",
        )
        frame.to_file(workspace / "hds_target_locations.gpkg", driver="GPKG")
        pd.DataFrame({"name": ["OBS_A", "OBS_B"], "layer": [0, 1], "cell": [0, 1]}).to_csv(
            workspace / "hds_head_target_map.csv", index=False
        )
        observation_sets.append({
            "kind": "head_targets", "prefix": "hds",
            "locations_file": "hds_target_locations.gpkg",
            "values_file": "hds_target_values.csv",
        })
        measured_rows = []
        # Two MEASURED periods each, so the mean-over-time aggregation runs.
        for location, measured, offset in (
            ("obs_a", 10.0, head_residuals[0]), ("obs_b", 20.0, head_residuals[1])
        ):
            for period in (0, 1):
                name = f"oname:hds_otype:lst_usecol:{location}_per:{period}"
                names.append((name, measured, 1.0))
                values[name] = measured + offset
                measured_rows.append({"time": period, "name": location, "head": measured})
            for period in range(2, 2 + unmeasured_periods):
                # pyEMU's phantom row: obsval IS the simulated value, weight 1.0,
                # and it is absent from the target-values snapshot.
                name = f"oname:hds_otype:lst_usecol:{location}_per:{period}"
                names.append((name, 999.0, 1.0))
                values[name] = 999.0
        pd.DataFrame(measured_rows).to_csv(
            workspace / "hds_target_values.csv", index=False
        )

    if forecast:
        # cal.forecast() registers an ordinary observation set at weight 0.
        frame = gpd.GeoDataFrame(
            {"name": ["FORE_1"], "layer": [0], "group": ["f"], "weight": [0.0]},
            geometry=[Point(1.5, 0.5)], crs="EPSG:2927",
        )
        frame.to_file(workspace / "fore1_target_locations.gpkg", driver="GPKG")
        observation_sets.append({
            "kind": "head_targets", "prefix": "fore1",
            "locations_file": "fore1_target_locations.gpkg",
            "values_file": "fore1_target_values.csv",
        })
        pd.DataFrame([{"time": 0, "name": "fore_1", "head": 5.0}]).to_csv(
            workspace / "fore1_target_values.csv", index=False
        )
        name = "oname:fore1_otype:lst_usecol:fore_1_per:0"
        names.append((name, 5.0, 0.0))
        values[name] = 105.0  # a huge "residual" that must never reach the map

    if conc:
        # Written in the CURRENT metadata shape -- `geometry`, `mapping_file`
        # and `value_column` all declared -- while the head block above stays
        # in the older shape, so one fixture covers both the declared path and
        # the kind-keyed fallback that keeps on-disk runs reviewable.
        frame = gpd.GeoDataFrame(
            {"name": ["MW_1"], "layer": [0], "group": ["c"], "weight": [1.0]},
            geometry=[Point(1.5, 0.5)],
            crs="EPSG:2927",
        )
        frame.to_file(workspace / "conc_target_locations.gpkg", driver="GPKG")
        # Deliberately NOT the `_head_target_map.csv` name the old code
        # hardcoded: finding this file proves the mapping name is read from
        # metadata rather than assumed.
        pd.DataFrame({"name": ["MW_1"], "layer": [0], "cell": [1]}).to_csv(
            workspace / "conc_conc_target_map.csv", index=False
        )
        observation_sets.append({
            "kind": "conc_targets", "prefix": "conc", "geometry": "points",
            "locations_file": "conc_target_locations.gpkg",
            "mapping_file": "conc_conc_target_map.csv",
            "values_file": "conc_target_values.csv",
            "value_column": "conc",
        })
        pd.DataFrame([{"time": 0, "name": "mw_1", "conc": 1.0}]).to_csv(
            workspace / "conc_target_values.csv", index=False
        )
        name = "oname:conc_otype:lst_usecol:mw_1_per:0"
        names.append((name, 1.0, 1.0))
        values[name] = 1.0 + conc_residual

    if zones:
        frame = gpd.GeoDataFrame(
            {"name": ["ZONE_1"], "group": ["z"], "weight": [1.0], "cells": [zone_cells]},
            geometry=[Polygon([(0, 0), (2, 0), (2, 1), (0, 1)])],
            crs="EPSG:2927",
        )
        frame.to_file(workspace / "drn_flow_target_locations.gpkg", driver="GPKG")
        observation_sets.append({
            "kind": "drn_flow", "prefix": "drn_flow",
            "locations_file": "drn_flow_target_locations.gpkg",
            "values_file": "drn_flow_target_values.csv",
        })
        pd.DataFrame([{"time": 1.0, "name": "zone_1", "flow_target": 100.0}]).to_csv(
            workspace / "drn_flow_target_values.csv", index=False
        )
        name = "oname:drn_flow_otype:lst_usecol:zone_1_time:1.0"
        names.append((name, 100.0, 1.0))
        values[name] = 100.0 + zone_residual

    observation_data = pd.DataFrame(
        {"obsnme": [n for n, _, _ in names],
         "obsval": [v for _, v, _ in names],
         "weight": [w for _, _, w in names],
         "obgnme": [n.rsplit("_", 1)[0] for n, _, _ in names]}
    ).set_index("obsnme")
    posterior = pd.DataFrame(
        [values, {k: v + 0.5 for k, v in values.items()}], index=["base", "1"]
    )

    class _Stub(IesResults):
        def __init__(self):
            self.workspace = workspace
            self.model = SimpleNamespace(vor=vor)
            self.pst = SimpleNamespace(observation_data=observation_data)
            self.__dict__["_metadata"] = {"observation_sets": observation_sets}

        @property
        def posterior(self):
            return SimpleNamespace(_df=posterior)

        @property
        def posterior_iteration(self):
            return 3

    return _Stub()


def test_obs_residuals_joins_pest_names_back_to_their_grid_locations(tmp_path):
    """The whole point: a bare obs name recovers where on the grid it lives.

    ``pst.try_parse_name_metadata()``'s own ``usecol`` column truncates at the
    first underscore -- ``obs_a`` arrives as ``obs`` -- so a location name with
    an underscore in it (i.e. nearly all of them) needs the name parsed whole.
    """

    frame = _stub_ies_results(tmp_path).obs_residuals()

    assert list(frame["location"]) == ["obs_a", "obs_b"]
    assert list(frame["cell"]) == [0, 1]
    assert frame.loc[frame["location"] == "obs_a", "x"].item() == pytest.approx(0.5)
    # simulated - measured, averaged over the two periods each was measured in.
    assert list(frame["residual"]) == pytest.approx([2.0, -1.0])
    assert list(frame["n"]) == [2, 2]


def test_a_residual_is_simulated_minus_measured_not_the_other_way(tmp_path):
    """The sign convention is named, never inferred.

    PEST's own ``.res`` file reports ``measured - modelled``; myflopy reports
    ``simulated - measured`` (what ``phi_contributions`` already uses). A model
    simulating 12 where 10 was measured is over-simulating, so the residual is
    POSITIVE -- flip this and every color on the map inverts.
    """

    frame = _stub_ies_results(tmp_path).obs_residuals()
    over = frame.loc[frame["location"] == "obs_a"].iloc[0]

    assert over["simulated"] == pytest.approx(12.0)
    assert over["measured"] == pytest.approx(10.0)
    assert over["residual"] == pytest.approx(2.0)


def test_residuals_fall_back_to_the_ensemble_mean_when_base_is_absent(tmp_path):
    """Mirrors phi_contributions: an unknown realization label is not an error."""

    results = _stub_ies_results(tmp_path)
    base = results.obs_residuals(realization="base")
    mean = results.obs_residuals(realization="nope")

    # The stub's second realization is +0.5 everywhere, so the mean sits halfway.
    assert mean.loc[0, "residual"] == pytest.approx(base.loc[0, "residual"] + 0.25)


def test_drn_zones_color_the_cells_they_cover(tmp_path):
    """A zone target has no point; it paints its snapshotted cell list."""

    frame = _stub_ies_results(tmp_path, heads=False, zones=True).obs_residuals()

    assert list(frame["location"]) == ["zone_1"]
    assert list(frame.loc[0, "cells"]) == [0, 1]
    assert frame.loc[0, "residual"] == pytest.approx(-6.0)


def test_points_and_zones_read_on_one_scale(tmp_path):
    """A point and the cell under it must mean the same thing at the same color.

    Both halves of the figure are built from one symmetric limit over ALL
    residuals; separate autoscaling would let a +2 point and a +2 cell render
    differently on the same map.
    """

    results = _stub_ies_results(tmp_path, heads=True, zones=True)
    figure = results.plot_obs_residuals()
    cells = next(trace for trace in figure.data if getattr(trace, "z", None) is not None)
    markers = next(trace for trace in figure.data if trace.type == "scattermap")

    assert cells.zmin == pytest.approx(markers.marker.cmin)
    assert cells.zmax == pytest.approx(markers.marker.cmax)
    assert cells.zmin == pytest.approx(-cells.zmax)
    assert cells.colorscale == markers.marker.colorscale
    # Exactly one scale bar: the cells carry it whenever there are zones.
    assert markers.marker.showscale is False


def test_a_heads_only_residual_map_still_shows_a_scale(tmp_path):
    """With no zones every cell is NaN, so the markers must carry the legend.

    Otherwise the figure has an empty colorbar over blank cells and no key at
    all for the only data actually drawn.
    """

    figure = _stub_ies_results(tmp_path, heads=True, zones=False).plot_obs_residuals()
    cells = next(trace for trace in figure.data if getattr(trace, "z", None) is not None)
    markers = next(trace for trace in figure.data if trace.type == "scattermap")

    assert cells.showscale is False
    assert markers.marker.showscale is True
    assert np.isnan(np.asarray(cells.z, dtype=float)).all()


def test_the_static_residual_map_draws_its_points_in_model_coordinates(tmp_path):
    """plot_mpl ignores Choro overlays entirely, so the points are drawn onto
    its axes directly -- in MODEL coordinates, not the lon/lat the Plotly path
    needs. Getting that frame wrong puts every observation off the map."""

    figure = _stub_ies_results(tmp_path).plot_obs_residuals(backend="matplotlib")
    scatter = figure.axes[0].collections[-1]
    offsets = np.asarray(scatter.get_offsets())

    assert offsets.tolist() == [[0.5, 0.5], [1.5, 0.5]]
    assert scatter.get_clim() == pytest.approx((-2.0, 2.0))


def test_a_run_with_no_locatable_observations_says_so(tmp_path):
    """Lake and SFR targets record only a lake or reach number. Refusing with a
    named reason beats drawing an empty grid that looks like zero residuals."""

    results = _stub_ies_results(tmp_path, heads=False, zones=False)

    assert results.obs_residuals().empty
    with pytest.raises(ValueError, match="lake or reach number"):
        results.plot_obs_residuals()


# --- concentration on the residual map (ledger 105) --------------------------


def test_concentration_targets_reach_the_residual_map(tmp_path):
    """Concentration is sampled at a point exactly as head is, so it belongs on
    the misfit map. Before ledger 105 the coordinates were computed by
    ``match_to_model`` and then dropped, and the map silently omitted every
    transport observation -- a calibration whose only data was concentration
    drew a blank grid."""

    results = _stub_ies_results(tmp_path, heads=False, conc=True)
    frame = results.obs_residuals()

    assert list(frame["location"]) == ["mw_1"]
    assert list(frame["kind"]) == ["conc_targets"]
    assert frame["x"].item() == pytest.approx(1.5)
    # Read from `mapping_file`, not the head suffix the old code hardcoded.
    assert frame["cell"].item() == 1
    assert frame["residual"].item() == pytest.approx(0.02)

    figure = results.plot_obs_residuals()
    markers = next(trace for trace in figure.data if trace.type == "scattermap")
    assert len(markers.lon) == 1


def test_a_family_that_declares_no_geometry_is_still_skipped(tmp_path):
    """Dispatch moved from kind to ``geometry``, and it must stay just as
    closed: lake and SFR sets snapshot a lake or reach NUMBER, and placing one
    on the grid would put a residual somewhere it was never measured."""

    results = _stub_ies_results(tmp_path, heads=True)
    results.__dict__["_metadata"]["observation_sets"].append({
        "kind": "lake_stage", "prefix": "lak",
        "locations_file": "hds_target_locations.gpkg",
        "values_file": "hds_target_values.csv",
    })

    assert set(results.obs_residuals()["prefix"]) == {"hds"}


def test_mixing_observation_kinds_on_one_color_scale_warns(tmp_path):
    """Heads are a length and concentration a mass per volume. On one shared
    diverging scale the bigger family sets the limit and the other renders
    uniformly white -- which reads as a perfect fit, not as an unreadable
    figure. Warn and name the way out."""

    results = _stub_ies_results(tmp_path, heads=True, conc=True)

    with pytest.warns(UserWarning, match="prefix="):
        results.plot_obs_residuals()

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        figure = results.plot_obs_residuals(prefix="conc")

    markers = next(trace for trace in figure.data if trace.type == "scattermap")
    assert len(markers.lon) == 1, "prefix= must draw only the selected family"
    assert markers.marker.cmax == pytest.approx(0.02), (
        "the color limit must come from the selected family alone"
    )


def test_selecting_a_prefix_that_was_never_recorded_names_the_options(tmp_path):
    """A silent empty map would look like a converged model with no misfit."""

    results = _stub_ies_results(tmp_path, heads=True, conc=True)

    with pytest.raises(ValueError, match=r"\['conc', 'hds'\]"):
        results.obs_residuals(prefix="nope")


def test_a_mosaic_refuses_the_stats_that_have_no_prior_and_posterior_form(tmp_path):
    """``change``, ``reduction`` and ``base`` already summarize BOTH ensembles.

    Composing one of them "prior vs posterior" would draw the same map twice.
    The message names plot_field so the reader is not left guessing.
    """

    results = _stub_ies_results(tmp_path)
    for stat in ("change", "reduction", "base"):
        with pytest.raises(ValueError, match="plot_field"):
            results.plot_field_mosaic("k", stat=stat)


def test_a_mosaic_refuses_matplotlib_rather_than_silently_dropping_panels(tmp_path):
    """viz.mosaic composes Plotly subplots only; the message names the way out."""

    with pytest.raises(ValueError, match="backend='matplotlib'"):
        _stub_ies_results(tmp_path).plot_field_mosaic("k", backend="matplotlib")


def test_a_mosaic_of_one_panel_is_not_a_mosaic(tmp_path):
    with pytest.raises(ValueError, match="at least two panels"):
        _stub_ies_results(tmp_path).plot_field_mosaic("k", which=("posterior",))


_UNSET = object()


def _four_cell_vor_clockwise():
    """A 1x4 strip, for asserting on a field where only SOME cells were captured."""

    verts = np.array(
        [[float(x), y] for x in range(5) for y in (0.0, 1.0)], dtype=float
    )
    # vertex ids: cell i uses (2i, 2i+1, 2i+3, 2i+2) -> clockwise from bottom-left
    iverts = [[2 * i, 2 * i + 1, 2 * i + 3, 2 * i + 2] for i in range(4)]
    xcyc = np.array([[i + 0.5, 0.5] for i in range(4)], dtype=float)
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)


def _stub_ies_with_capture(tmp_path, prior_sds, posterior_sds, means=(1.0, 1.0), *,
                           prior_means=None, posterior_means=None, cells=(0, 1),
                           vor=None, realizations=("base", "1"), family="array",
                           model=_UNSET):
    """A stub carrying a CAPTURED FIELD, for the `field()` statistics and maps.

    Two realizations per ensemble, with the per-cell spread chosen by the caller
    so the derived columns have known values.

    `prior_means`/`posterior_means` both default to `means`; pass them
    separately to make the two ensembles DIFFER, which is what `which=`,
    `change` and a two-panel mosaic are about. `cells` shorter than the grid
    leaves the rest uncaptured, exercising the NaN padding onto the full grid.
    """

    from myflopy.modflow.mf6.pest.ies import IesResults

    workspace = tmp_path / "capture"
    workspace.mkdir(parents=True, exist_ok=True)
    names = [f"oname:kfieldl0_otype:arr_i:{cell}_j:0" for cell in cells]

    def _ensemble(cell_means, sds):
        # Two realizations symmetric about the mean -> ddof=1 std is exactly sd.
        rows = {
            name: [mean - sd / 2 ** 0.5, mean + sd / 2 ** 0.5]
            for name, mean, sd in zip(names, cell_means, sds, strict=True)
        }
        return pd.DataFrame(rows, index=list(realizations))

    prior = _ensemble(prior_means or means, prior_sds)
    posterior = _ensemble(posterior_means or means, posterior_sds)
    observation_data = pd.DataFrame(
        {"obsnme": names, "obsval": [0.0] * len(names), "weight": [0.0] * len(names)}
    ).set_index("obsnme")

    class _Stub(IesResults):
        def __init__(self):
            self.workspace = workspace
            self.model = (
                SimpleNamespace(vor=vor or _two_cell_vor_clockwise())
                if model is _UNSET else model
            )
            self.pst = SimpleNamespace(observation_data=observation_data)
            self.__dict__["_metadata"] = {"capture_fields": [
                {"prefix": "kfield", "layer_prefixes": {"0": "kfieldl0"},
                 "target": "k", "family": family, "name": "k"}
            ]}

        @property
        def prior(self):
            return SimpleNamespace(_df=prior)

        @property
        def posterior(self):
            return SimpleNamespace(_df=posterior)

        @property
        def prior_iteration(self):
            return 0

        @property
        def posterior_iteration(self):
            return 3

    return _Stub()


def test_uncertainty_reduction_is_the_share_of_the_prior_spread_removed(tmp_path):
    """`1 - posterior_sd / prior_sd`, and nothing else.

    Cell 0's spread was halved (reduction 0.5); cell 1's was untouched
    (reduction 0). Inverting the ratio, or dropping the `1 -`, changes both
    numbers -- which the color-policy tests cannot see, because they are handed
    value lists rather than computing them.
    """

    frame = _stub_ies_with_capture(tmp_path, prior_sds=(2.0, 1.0),
                                   posterior_sds=(1.0, 1.0)).field("k")

    assert list(frame["prior_std"]) == pytest.approx([2.0, 1.0])
    assert list(frame["posterior_std"]) == pytest.approx([1.0, 1.0])
    assert list(frame["reduction"]) == pytest.approx([0.5, 0.0])


def test_a_posterior_spread_that_grew_reads_as_a_negative_reduction(tmp_path):
    """The case `_field_map_policy` has a dedicated fallback branch for, which
    nothing else ever produces from real ensembles: more uncertain after
    calibration than before."""

    frame = _stub_ies_with_capture(
        tmp_path, prior_sds=(1.0, 1.0), posterior_sds=(2.0, 0.5)
    ).field("k")

    assert frame.loc[0, "reduction"] == pytest.approx(-1.0)  # spread doubled
    assert frame.loc[1, "reduction"] == pytest.approx(0.5)   # spread halved


def test_phantom_observations_never_reach_a_residual(tmp_path):
    """pyEMU makes one observation per simulated ROW, not per measured row.

    Unmeasured times keep weight 1.0 and an obsval equal to the model's own
    output, so averaging them in drags every residual toward zero -- and with
    enough of them it flips the sign, painting the map the opposite color. They
    are excluded by joining back to the target-values snapshot.
    """

    honest = _stub_ies_results(tmp_path / "a").obs_residuals()
    padded = _stub_ies_results(tmp_path / "b", unmeasured_periods=4).obs_residuals()

    assert list(padded["residual"]) == pytest.approx(list(honest["residual"]))
    assert list(padded["n"]) == [2, 2]  # not 6: the four phantoms are not data


def test_forecasts_are_not_drawn_as_calibration_misfit(tmp_path):
    """`cal.forecast(...)` registers an ordinary observation set at weight 0.

    Nothing in the persisted metadata marks it as a forecast, so it arrives
    looking exactly like a head target. Drawing it on a figure captioned "where
    is the model biased" would be wrong twice over: it was never fitted, and its
    huge apparent residual would set the whole map's color scale.
    """

    results = _stub_ies_results(tmp_path, forecast=True)
    frame = results.obs_residuals()

    assert "fore_1" not in set(frame["location"])
    assert set(frame["prefix"]) == {"hds"}
    figure = results.plot_obs_residuals()
    cells = next(t for t in figure.data if getattr(t, "z", None) is not None)
    assert cells.zmax == pytest.approx(2.0)  # not 100, the forecast's "residual"


def test_a_duplicated_location_snapshot_does_not_multiply_the_markers(tmp_path):
    """Two rows for one name (duplicated source row, two screens) used to
    cross-join into n**2 rows against a single pyEMU observation."""

    results = _stub_ies_results(tmp_path)
    mapping = pd.read_csv(results.workspace / "hds_head_target_map.csv")
    pd.concat([mapping, mapping]).to_csv(
        results.workspace / "hds_head_target_map.csv", index=False
    )

    assert len(results.obs_residuals()) == 2


def test_one_degenerate_zone_does_not_take_down_the_whole_frame(tmp_path):
    """A DRN zone that intersects no cells is snapshotted as an empty cell list,
    which round-trips through CSV as NaN. Parsing that unguarded raised, killing
    the residuals for every head target in the run too."""

    frame = _stub_ies_results(tmp_path, zones=True, zone_cells="").obs_residuals()

    assert list(frame.loc[frame["location"] == "zone_1", "cells"].item()) == []
    assert "obs_a" in set(frame["location"])  # the heads survived


def test_the_shared_limit_is_set_by_whichever_half_is_larger(tmp_path):
    """Pins the UNION, not merely that the two traces agree with each other.

    Here the largest residual is a head point, so a limit computed from the
    zones alone would be too small and the point would clip.
    """

    results = _stub_ies_results(tmp_path, zones=True, head_residuals=(12.0, -1.0),
                                zone_residual=-6.0)
    figure = results.plot_obs_residuals()
    cells = next(t for t in figure.data if getattr(t, "z", None) is not None)
    markers = next(t for t in figure.data if t.type == "scattermap")

    assert cells.zmax == pytest.approx(12.0)
    assert markers.marker.cmax == pytest.approx(12.0)


def test_the_static_map_draws_points_on_the_map_axes_not_the_colorbar(tmp_path):
    """With zones present the figure has two axes, and only then can an
    axes[-1] slip put every observation onto the colorbar instead of the map."""

    figure = _stub_ies_results(tmp_path, zones=True).plot_obs_residuals(
        backend="matplotlib"
    )
    map_axes = figure.axes[0]

    assert len(figure.axes) == 2  # map + colorbar, so the two indices differ
    # The cells' patches plus our points, both on the MAP axes. Scattering onto
    # axes[-1] leaves this at 1 and the offsets below become the patches'.
    assert len(map_axes.collections) == 2
    offsets = np.asarray(map_axes.collections[-1].get_offsets())
    assert offsets.tolist() == [[0.5, 0.5], [1.5, 0.5]]


def test_a_heads_only_static_map_still_has_a_color_key(tmp_path):
    """plot_mpl ignores Choro overlays, so the scattered points cannot carry a
    legend the way the Plotly markers do -- the cells' bar is the only key this
    backend has, and it must be drawn even over an all-NaN cell column."""

    figure = _stub_ies_results(tmp_path).plot_obs_residuals(backend="matplotlib")

    assert len(figure.axes) == 2
    assert figure.axes[-1].get_ylim() == pytest.approx((-2.0, 2.0))


def test_a_mixed_case_observation_name_still_joins(tmp_path):
    """The snapshot side is force-lowercased, so if a control file ever carries
    mixed-case `usecol` values the inner merge would drop EVERY row -- returning
    an empty frame and blaming lake/SFR targets for it. Both sides lower."""

    def _shout(name):
        return name.replace("usecol:obs_a", "usecol:OBS_A")

    results = _stub_ies_results(tmp_path)
    # Both sides of the control file, exactly as a real run would carry them --
    # renaming only the pst index would drop the rows before the join under test.
    results.pst = SimpleNamespace(
        observation_data=results.pst.observation_data.rename(index=_shout)
    )
    posterior = results.posterior._df.rename(columns=_shout)
    type(results).posterior = property(lambda self: SimpleNamespace(_df=posterior))

    assert "obs_a" in set(results.obs_residuals()["location"])


def test_a_linear_field_mosaic_is_not_labelled_in_decades(tmp_path):
    """The colorbar callback keys on the EFFECTIVE logscale, not on the stat.

    `mean` is log by default, but an explicit `logscale=False` makes the pooled
    data linear -- where decade ticks are not merely wrong, they can overflow
    (`10 ** 314`) and take the whole figure down.
    """

    results = _stub_ies_with_capture(tmp_path, prior_sds=(2.0, 1.0),
                                     posterior_sds=(1.0, 1.0), means=(1e2, 3e2))
    figure = results.plot_field_mosaic("k", stat="mean", logscale=False)

    assert figure.layout.coloraxis.colorbar.ticktext is None
    assert figure.layout.coloraxis.cmax == pytest.approx(300.0)  # real units, not log


def test_a_linear_single_field_map_is_not_labelled_in_decades(tmp_path):
    """Same override on the single map: the policy's decade colorbar carries
    tickvals in LOG space, which over linear data crush every label into the
    bottom of the bar."""

    results = _stub_ies_with_capture(tmp_path, prior_sds=(2.0, 1.0),
                                     posterior_sds=(1.0, 1.0), means=(1.5, 4.5))
    trace = results.plot_field("k", stat="mean", logscale=False).data[0]

    assert trace.colorbar.ticktext is None
    assert trace.colorbar.tickvals is None


# -- 6.4C: fast coverage for the field maps ---------------------------------
# The capture stub above drives these without pyEMU, pestpp or MODFLOW; the slow
# end-to-end test proves the same path against a real run.


def test_a_prior_field_map_maps_the_prior_and_says_so(tmp_path):
    """`which=` picks the ensemble, and nothing pinned it before.

    Until the stub could give prior and posterior DIFFERENT means, both panels
    carried identical values -- so swapping `which` (or ignoring it) passed the
    whole suite. The title has to move too, or a reader cannot tell the two
    figures apart.
    """

    results = _stub_ies_with_capture(
        tmp_path, prior_sds=(4.0, 4.0), posterior_sds=(1.0, 1.0),
        prior_means=(0.001, 100.0), posterior_means=(0.01, 10.0),
    )

    prior = results.plot_field("k", stat="mean", which="prior")
    posterior = results.plot_field("k", stat="mean", which="posterior")

    # mean is log-scaled by policy, so these are log10 of the real means
    assert list(prior.data[0].z) == pytest.approx([-3.0, 2.0])
    assert list(posterior.data[0].z) == pytest.approx([-2.0, 1.0])
    # the label has to move with the data, iteration included -- a map showing
    # the prior under the posterior's iteration number is worse than no label
    assert "(prior)" in prior.layout.title.text
    assert "iteration 0" in prior.layout.title.text
    assert "(posterior)" in posterior.layout.title.text
    assert "iteration 3" in posterior.layout.title.text

    # `std` selects a column the same way, and is linear -- so a `which` that is
    # honored only on the log-scaled `mean` path would slip through above
    assert list(results.plot_field("k", stat="std", which="prior").data[0].z) \
        == pytest.approx([4.0, 4.0])
    assert list(results.plot_field("k", stat="std", which="posterior").data[0].z) \
        == pytest.approx([1.0, 1.0])


def test_uncaptured_cells_stay_blank_on_a_partly_captured_field(tmp_path):
    """`field()` has one row per CAPTURED cell, but every trace column must be
    exactly ncpl long -- so the values are scattered onto the full grid with NaN
    left behind. Both the stub and the end-to-end test used to capture every
    cell, so this padding ran in no test at all."""

    results = _stub_ies_with_capture(
        tmp_path, prior_sds=(1.5, 2.5), posterior_sds=(1.5, 2.5),
        means=(2.0, 20.0), cells=(0, 3), vor=_four_cell_vor_clockwise(),
    )

    frame = results.field("k")
    assert list(frame["cell"]) == [0, 3]
    assert len(frame) == 2 < int(results.model.vor.ncpl)

    trace = results.plot_field("k", stat="mean").data[0]
    assert len(trace.z) == 4 and len(trace.customdata) == 4
    assert trace.z[0] == pytest.approx(np.log10(2.0))
    assert trace.z[3] == pytest.approx(np.log10(20.0))

    # Asserted on a LINEAR stat: `mean` is log-scaled, and log10(0) is blanked
    # anyway, so a zero-filled pad would read as NaN there and this would pass
    # against a padding that quietly invents "0" for cells nobody captured.
    spread = results.plot_field("k", stat="std").data[0]
    assert list(spread.z)[0] == pytest.approx(1.5)
    assert list(spread.z)[3] == pytest.approx(2.5)
    assert np.isnan(spread.z[1]) and np.isnan(spread.z[2])


def test_a_field_frame_carries_cell_as_a_column_not_an_index(tmp_path):
    """The docstring promises this explicitly -- callers merge the frame on
    `cell`, which silently produces nothing if it is the index instead."""

    frame = _stub_ies_with_capture(
        tmp_path, prior_sds=(1.0, 1.0), posterior_sds=(1.0, 1.0),
        cells=(3, 0), vor=_four_cell_vor_clockwise(),
    ).field("k")

    assert "cell" in frame.columns
    assert list(frame.index) == list(range(len(frame)))
    assert list(frame["cell"]) == [0, 3]  # sorted, whatever order they arrived in


def test_a_change_map_plots_the_log_ratio_but_hovers_the_raw_ratio(tmp_path):
    """What is plotted is log10(post/prior) so the scale is symmetric about no
    change; what the reader is shown on hover is the ratio itself. Conflating
    the two turns "10x more conductive" into "1"."""

    results = _stub_ies_with_capture(
        tmp_path, prior_sds=(1.0, 1.0), posterior_sds=(1.0, 1.0),
        prior_means=(1.0, 10.0), posterior_means=(10.0, 1.0),
    )

    frame = results.field("k")
    assert list(frame["change"]) == pytest.approx([10.0, 0.1])

    trace = results.plot_field("k", stat="change").data[0]
    assert list(trace.z) == pytest.approx([1.0, -1.0])  # log10(10), log10(0.1)
    # ...while the hover carries the ratio itself (customdata[1], per the template;
    # [0] is the cell id). Pre-formatted strings, hence the float().
    assert "posterior / prior" in trace.hovertemplate
    assert [float(row[1]) for row in trace.customdata] == pytest.approx([10.0, 0.1])


def test_a_mosaic_panel_shows_its_own_ensemble_on_the_pooled_scale(tmp_path):
    """The point of a mosaic is that the two panels differ but are comparable:
    each draws its own ensemble, both on one shared color axis spanning both."""

    results = _stub_ies_with_capture(
        tmp_path, prior_sds=(1.0, 1.0), posterior_sds=(1.0, 1.0),
        prior_means=(0.001, 100.0), posterior_means=(0.01, 10.0),
    )

    figure = results.plot_field_mosaic("k", stat="mean")
    panels = [trace for trace in figure.data if getattr(trace, "z", None) is not None]

    assert len(panels) == 2
    assert list(panels[0].z) == pytest.approx([-3.0, 2.0])   # prior
    assert list(panels[1].z) == pytest.approx([-2.0, 1.0])   # posterior
    # pooled across BOTH panels, not just the first
    assert figure.layout.coloraxis.cmin == pytest.approx(-3.0)
    assert figure.layout.coloraxis.cmax == pytest.approx(2.0)


def test_a_static_field_map_reads_in_real_units_and_a_linear_one_does_not(tmp_path):
    """`plot_field`'s own call of `_relabel_log_colorbar` -- the helper is tested
    directly in test_colorscale_policy.py, but the wiring that decides WHEN to
    apply it is only asserted here. A log map whose bar reads -2..1 is wrong."""

    results = _stub_ies_with_capture(tmp_path, prior_sds=(0.5, 2.0),
                                     posterior_sds=(0.5, 2.0), means=(0.01, 10.0))

    logged = results.plot_field("k", stat="mean", backend="matplotlib")
    labels = [t.get_text() for t in logged.axes[-1].get_yticklabels() if t.get_text()]
    # Which decades matplotlib places inside the limits is its business; that
    # they are DECADES and not the log10 exponents (-2..1) is the contract.
    assert labels and all(float(text) in (0.001, 0.01, 0.1, 1, 10, 100) for text in labels)
    assert {"0.1", "10"} <= set(labels)

    # std is linear by policy (a spread is legitimately 0), so it is not relabelled
    linear = results.plot_field("k", stat="std", backend="matplotlib")
    spread = [t.get_text() for t in linear.axes[-1].get_yticklabels() if t.get_text()]
    assert any(float(text) not in (0.001, 0.01, 0.1, 1, 10, 100) for text in spread)


def test_a_field_map_without_a_grid_names_the_way_to_attach_one(tmp_path):
    """A run reopened without its model can still show tables; the map has to
    say how to get the grid rather than failing on an attribute."""

    results = _stub_ies_with_capture(tmp_path, prior_sds=(1.0, 1.0),
                                     posterior_sds=(1.0, 1.0), model=None)

    with pytest.raises(ValueError, match="Spatial maps need the model grid"):
        results.field("k")


def test_an_unavailable_field_lists_what_is_available(tmp_path):
    """Both lookups name the alternatives -- an empty message here means opening
    a notebook and grepping the metadata by hand."""

    results = _stub_ies_with_capture(tmp_path, prior_sds=(1.0, 1.0),
                                     posterior_sds=(1.0, 1.0))

    with pytest.raises(ValueError, match=r"Captured layers: \[0\]"):
        results.field("k", layer=1)
    with pytest.raises(KeyError, match="Available"):
        results.field("zzz")


def test_a_list_field_is_refused_by_name_rather_than_mapped_wrong(tmp_path):
    """Only array families (K, K33) have one value per cell. A list field would
    otherwise be scattered onto cell ids it has no relationship to."""

    results = _stub_ies_with_capture(tmp_path, prior_sds=(1.0, 1.0),
                                     posterior_sds=(1.0, 1.0), family="list")

    with pytest.raises(NotImplementedError, match="list field"):
        results.field("k")


def test_a_run_without_a_base_realization_has_no_base_stat(tmp_path):
    """`base` is the minimum-error-variance realization and is simply absent
    from some runs -- the column is skipped and the stat refused by name."""

    results = _stub_ies_with_capture(tmp_path, prior_sds=(1.0, 1.0),
                                     posterior_sds=(1.0, 1.0),
                                     realizations=("0", "1"))

    assert "base" not in results.field("k").columns
    with pytest.raises(ValueError, match="or unavailable"):
        results.plot_field("k", stat="base")


# -- 6.4C: a run with no registered forecasts (ledger 73) --------------------
# `cal.forecast(...)` is optional, so `pestpp_options` carries no "forecasts"
# key at all on an ordinary calibration -- and pyEMU answers None there.


def _forecastless_ies(tmp_path, *, forecasts=None):
    """An IesResults over a REAL pyemu Pst that registers no forecasts.

    Deliberately not `_stub_ies_results`: its SimpleNamespace `pst` has no
    `pestpp_options` at all, so it would raise AttributeError and never reach
    the defect this pins.
    """

    pytest.importorskip("pyemu")
    import pyemu

    from myflopy.modflow.mf6.pest.ies import IesResults

    pst = pyemu.Pst.from_par_obs_names(par_names=["p1", "p2"], obs_names=["o1", "o2"])
    if forecasts is not None:
        pst.pestpp_options["forecasts"] = forecasts

    class _Stub(IesResults):
        def __init__(self):
            self.workspace = tmp_path
            self.model = None
            self.case = "noforecast"
            self.pst = pst
            self.__dict__["_metadata"] = {}

        @property
        def iterations(self):
            return [0, 1]

        @property
        def noise(self):
            return None

    return _Stub()


def test_forecast_names_is_empty_not_none_when_no_forecasts_registered(tmp_path):
    """The regression anchor for ledger 73. pyEMU returns None (not []) when the
    `forecasts` option was never written, and `list(None)` is a TypeError -- so
    `settings`, `report()`, `forecasts()` and `forecast()` all died on the most
    ordinary kind of run there is."""

    assert _forecastless_ies(tmp_path).forecast_names == []


def test_settings_renders_on_a_run_with_no_forecasts(tmp_path):
    settings = _forecastless_ies(tmp_path).settings

    assert settings.n_forecasts == 0
    assert "forecasts: 0" in str(settings)


def test_forecasts_table_is_empty_not_an_error_when_none_registered(tmp_path):
    """`forecasts()` already carried an `if not rows` guard -- it was simply
    unreachable, because the property raised before the loop could run."""

    assert _forecastless_ies(tmp_path).forecasts().empty


def test_forecast_lookup_reports_the_available_names_when_there_are_none(tmp_path):
    """The authored KeyError names what you could have asked for; the TypeError
    it used to raise instead named nothing."""

    with pytest.raises(KeyError, match=r"Available: \[\]"):
        _forecastless_ies(tmp_path).forecast("nope")


def test_an_empty_forecasts_option_is_not_read_as_one_nameless_forecast(tmp_path):
    """pyEMU splits an empty option string into `[""]`. Taken at face value that
    reports one forecast and sends `report()` looking up an ensemble column
    named "" -- so the falsy names are filtered, not merely None-guarded."""

    assert _forecastless_ies(tmp_path, forecasts="").forecast_names == []


def test_report_writes_html_on_a_run_with_no_forecasts(tmp_path):
    """Covers BOTH of report()'s forecast sites: the loop, and the settings
    block at the very end -- which is built after every figure, so fixing only
    the loop would still blow up on the last line."""

    results = _forecastless_ies(tmp_path)
    # Stand in for the phi/obs figures, which need ensemble CSVs on disk; the
    # forecast handling under test is independent of them.
    blank = viz.Fig()
    for name in ("plot_phi", "plot_phi_distribution", "plot_vs_obs",
                 "plot_phi_contributions", "plot_parameters_at_bounds"):
        setattr(results, name, lambda *a, **k: blank)

    path = results.report(tmp_path / "review.html")

    assert path.exists()
    assert "forecasts: 0" in path.read_text(encoding="utf-8")


def _run_with_a_pest_build(tmp_path, monkeypatch, model_result):
    """A `Run` over a workspace holding one discovered PEST build.

    `model_result` is called in place of `Run.model()` — return a stand-in model
    or raise, to drive the two branches of `Run._pest_review_model`. `Run` is a
    slots dataclass, so the patch goes on the CLASS.
    """

    from myflopy.modflow.mf6.pest import ies as ies_module
    from myflopy.workspace import Run

    template = tmp_path / "pest" / "demo"
    template.mkdir(parents=True)
    (template / "myflopy_pest_metadata.json").write_text(
        '{"project_name": "demo", "model_name": "m", "pst_file": "demo.pst"}',
        encoding="utf-8",
    )

    opened = {}
    monkeypatch.setattr(ies_module, "open_ies_run",
                        lambda target, **kwargs: opened.update(kwargs))
    monkeypatch.setattr(Run, "model", lambda self, name=None: model_result())

    return Run(name="r", workspace=tmp_path), opened


def test_a_loaded_run_hands_its_model_to_review_without_resolving_it_to_list(
    tmp_path, monkeypatch
):
    """`run.pest_runs` used to attach no model, so `review()` returned an
    IesResults that refused every spatial map -- while both PEST notebooks said
    it was interchangeable with `model.pest_runs`.

    Goes through `Run.pest_runs` itself, not `find_pest_runs`: the wiring IS the
    fix, and a test that calls `find_pest_runs(model_factory=...)` by hand passes
    just as happily with the fix reverted.
    """

    resolved = []
    sentinel = object()
    run, opened = _run_with_a_pest_build(
        tmp_path, monkeypatch, lambda: (resolved.append(1), sentinel)[1]
    )

    handles = run.pest_runs
    assert len(handles) == 1
    assert resolved == []  # listing must not resolve a model

    handles[0].review()
    assert resolved == [1] and opened["model"] is sentinel


def test_a_run_that_cannot_name_one_model_reviews_without_a_grid(tmp_path, monkeypatch):
    """A run holding several models cannot pick one, and guessing would be worse
    than declining -- `review()` still opens, and the spatial maps raise their
    own message naming the fix."""

    def _ambiguous():
        raise ValueError("Model name is required because this run has multiple models")

    run, opened = _run_with_a_pest_build(tmp_path, monkeypatch, _ambiguous)
    run.pest_runs[0].review()

    assert opened["model"] is None


def test_an_explicit_review_model_beats_the_run_default(tmp_path, monkeypatch):
    """...and passing one explicitly must not resolve the run's model at all."""

    resolved = []
    chosen = object()
    run, opened = _run_with_a_pest_build(
        tmp_path, monkeypatch, lambda: resolved.append(1)
    )

    run.pest_runs[0].review(model=chosen)

    assert opened["model"] is chosen
    assert resolved == []
