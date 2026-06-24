from __future__ import annotations

import os
import subprocess
import sys
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

from myflopy.modflow.calcs.calibration import CalibrationPlot  # noqa: E402
from myflopy.modflow.mf6.observations import (  # noqa: E402
    DrnFlowTargets,
    HeadTargets,
    LakeStageTargets,
    SfrFlowTargets,
    SfrStageTargets,
    TargetRegistry,
)
from myflopy.modflow.mf6.pest.observations import (  # noqa: E402
    prepare_drn_flow_observations,
    prepare_lake_stage_observations,
    prepare_sfr_flow_observations,
    prepare_sfr_stage_observations,
)
from myflopy.modflow.mf6.pest.gis import derive_bounds  # noqa: E402
from myflopy.modflow.mf6.pest.project import PestProject  # noqa: E402
from myflopy.modflow.mf6.pest.specs import (  # noqa: E402
    DrnFlowObservationSpec,
    HeadTargetObservationSpec,
    ExpGeoStruct,
    LakeStageObservationSpec,
    SfrFlowObservationSpec,
    SfrStageObservationSpec,
    VectorParameterSource,
)
from myflopy.modflow.mf6.simulation.base import SimulationBase  # noqa: E402
from myflopy.modflow.mf6.simulation.discretization import DisvGrid, TemporalDiscretization  # noqa: E402
from myflopy.modflow.mf6.simulation.packages import (  # noqa: E402
    CHD,
    Drains,
    InitialConditions,
    KFlow,
    OutputControl,
    Recharge,
    Storage,
)
from myflopy.modflow.mf6.drn import DRNFromVector  # noqa: E402
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus  # noqa: E402
from myflopy.project.run_model import LoadedMf6Run  # noqa: E402


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
    # Default workspace lands beside the model so model.pest_runs finds it.
    assert cal.template_workspace == model.workspace / "pest" / "calib"

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

    cal = PestProject(
        model=model,
        name="native_demo",
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
    for tpl, inp in zip(pst.template_files, pst.input_files):
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

    cal = PestProject(
        model=model,
        name="prior_demo",
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

    # No explicit workspace: it defaults to <model workspace>/pest/<name>, so the
    # run is auto-discoverable via model.pest_runs (the workflow integration).
    cal = PestProject(
        model=model,
        name="run_ies_demo",
        start_datetime="2024-01-01",
    )
    assert cal.template_workspace == model.workspace / "pest" / "run_ies_demo"
    cal.parameterize("k", style="constant", bounds=(0.2, 5.0), physical=(1e-3, 100.0))
    cal.parameterize("recharge", style="grid", bounds=(0.5, 1.5), physical=(0.0, 1e-2))
    cal.observe(targets)
    cal.forecast(forecast)
    cal.build("run_ies_demo.pst", noptmax=0)

    ies = cal.run_ies(reals=6, iterations=2, noise=True)

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
    assert isinstance(ies.plot_parameters_at_bounds(backend="matplotlib"), MplFigure)
    contributions = ies.phi_contributions()
    assert not contributions.empty
    assert isinstance(ies.plot_phi_contributions(), go.Figure)
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

    cal = PestProject(
        model=model,
        name="cap_demo",
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

    ies = cal.run_ies(reals=6, iterations=1)

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


def test_canonical_calibration_demo_builds_native_pst_with_multilayer_k(tmp_path):
    """The canonical calibration demo wires the canonical valley model into a
    native PstFrom build, including multi-layer DISV K resolved per layer."""

    from myflopy.modflow.mf6.canonical_calibration import (
        build_canonical_calibration_demo,
    )
    from myflopy.modflow.mf6.pest import PestProject

    demo = build_canonical_calibration_demo(tmp_path / "model", n_head_wells=8)
    assert demo.model.gwf.modelgrid.nlay == 4
    assert demo.head_targets.locations_gdf.shape[0] == 8
    assert demo.forecast_targets.locations_gdf.shape[0] == 1
    assert demo.start_k_factor == 3.0

    cal = PestProject(
        model=demo.model,
        name="cc",
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
    from myflopy.modflow.mf6.pest import PestProject

    demo = build_canonical_calibration_demo(
        tmp_path / "model", n_head_wells=6, start_k_constant=10.0, start_k_layers=(0, 1)
    )
    cal = PestProject(
        model=demo.model, name="gk", workspace=tmp_path / "tmpl",
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
    from myflopy.modflow.mf6.pest import PestProject

    demo = build_canonical_calibration_demo(
        tmp_path / "model", n_head_wells=6, start_k_constant=10.0, start_k_layers=(0, 1)
    )
    cal = PestProject(
        model=demo.model, name="ppk", workspace=tmp_path / "tmpl",
        start_datetime="2024-01-01",
    )
    cal.parameterize("k", style="pilotpoints", pp_space=8, layers=[0, 1],
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
