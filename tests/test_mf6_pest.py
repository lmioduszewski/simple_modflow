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

from simple_modflow.modflow.calcs.calibration import CalibrationPlot  # noqa: E402
from simple_modflow.modflow.mf6.observations import HeadTargets, LakeStageTargets, TargetRegistry  # noqa: E402
from simple_modflow.modflow.mf6.pest.gis import derive_bounds  # noqa: E402
from simple_modflow.modflow.mf6.pest.parameters import build_k_pilotpoint_frame  # noqa: E402
from simple_modflow.modflow.mf6.pest.project import PestProject  # noqa: E402
from simple_modflow.modflow.mf6.pest.results import PestRunResults, PestRunReview, open_pest_run  # noqa: E402
from simple_modflow.modflow.mf6.pest.specs import (  # noqa: E402
    HeadTargetObservationSpec,
    DrainConductanceParameter,
    DrainElevationParameter,
    ExpGeoStruct,
    KPilotPointParameter,
    VectorParameterSource,
)
from simple_modflow.modflow.mf6.simulation.base import SimulationBase  # noqa: E402
from simple_modflow.modflow.mf6.simulation.discretization import DisvGrid, TemporalDiscretization  # noqa: E402
from simple_modflow.modflow.mf6.simulation.packages import (  # noqa: E402
    CHD,
    Drains,
    InitialConditions,
    KFlow,
    OutputControl,
    Storage,
)
from simple_modflow.modflow.mf6.drn import DRNFromVector  # noqa: E402
from simple_modflow.modflow.mf6.grid.voronoi import VoronoiGridPlus  # noqa: E402
from simple_modflow.project.run_model import LoadedMf6Run  # noqa: E402


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

    ax_obs = bound.plot.obs_vs_sim()
    assert ax_obs.get_title() == "Observed vs simulated heads"

    ax_ts = bound.plot.timeseries("OBS_A")
    assert ax_ts.get_title() == "OBS_A"
    assert len(ax_ts.lines) == 2

    ax_heads = bound.plot.calibration(type="heads")
    assert isinstance(ax_heads, CalibrationPlot)

    baseline_model = _fake_model_with_heads()
    baseline_model.all_heads = baseline_model.all_heads.copy()
    baseline_model.all_heads.loc[:, "elev"] = [9.8, 9.3, 10.0, 9.0]
    ax_obs_with_baseline = bound.plot.obs_vs_sim(baseline=baseline_model)
    assert len(ax_obs_with_baseline.collections) == 2
    ax_ts_with_baseline = bound.plot.timeseries("OBS_A", baseline=baseline_model)
    assert len(ax_ts_with_baseline.lines) == 3


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
        ("deep_lake", "STAGE", (0,)),
        ("shallow_lake", "STAGE", (1,)),
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


def test_build_k_pilotpoint_frame_from_vector_source():
    workspace = _project_temp_dir("pest_k_frame")
    model, _ = _two_cell_model("pest_k_frame", workspace)
    k_path = _write_gpkg(
        workspace / "hk.gpkg",
        gpd.GeoDataFrame(
            {
                "name": ["zone_a", "zone_b"],
                "unit": ["a", "b"],
                "k": [10.0, 20.0],
            },
            geometry=[
                Polygon([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]),
                Polygon([(1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0)]),
            ],
            crs=model.vor.crs,
        ),
    )

    project = PestProject(model=model, name="pest_k_frame", workspace=workspace / "template", start_datetime="2024-01-01")
    spec = KPilotPointParameter(
        name="hk",
        source=VectorParameterSource(path=k_path, value_column="k", zone_column="unit"),
        bounds=(0.25, 4.0),
        bounds_mode="multiplier",
        transform="log",
        pp_spacing=0.75,
        geostruct=ExpGeoStruct(range=100.0, transform="log"),
    )
    frame = build_k_pilotpoint_frame(project, spec)
    assert not frame.empty
    assert {"parnme", "pargp", "parval1", "parlbnd", "parubnd", "base_value", "zone"}.issubset(frame.columns)
    assert set(frame["zone"]) == {"a", "b"}
    assert set(frame["parval1"]) == {1.0}
    assert set(frame["parlbnd"]) == {0.25}
    assert set(frame["parubnd"]) == {4.0}


def test_pest_project_builds_first_slice_control_file_smoke():
    pyemu = pytest.importorskip("pyemu")
    workspace = _project_temp_dir("pest_project_smoke")
    model, _ = _two_cell_model("pest_project_smoke", workspace / "model")

    points_path = _write_gpkg(
        workspace / "targets.gpkg",
        gpd.GeoDataFrame(
            {
                "name": ["OBS_A", "OBS_B"],
                "layer": [0, 0],
                "weight": [1.0, 2.0],
            },
            geometry=[Point(0.5, 0.5), Point(1.5, 0.5)],
            crs=model.vor.crs,
        ),
    )
    target_values = pd.DataFrame({"per": [0, 1], "OBS_A": [10.0, 10.5], "OBS_B": [9.5, 9.0]})
    simulated_values = pd.DataFrame({"per": [0, 1], "OBS_A": [10.1, 10.4], "OBS_B": [9.4, 9.2]})
    targets = HeadTargets(locations=points_path, values=target_values, time_column="per")

    hk_path = _write_gpkg(
        workspace / "hk.gpkg",
        gpd.GeoDataFrame(
            {"name": ["hk_zone"], "unit": ["all"], "k": [15.0]},
            geometry=[Polygon([(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)])],
            crs=model.vor.crs,
        ),
    )
    drn_path = _write_gpkg(
        workspace / "drn.gpkg",
        gpd.GeoDataFrame(
            {
                "name": ["west", "east"],
                "group": ["main", "main"],
                "layer": [1, 1],
                "height": [9.0, 8.5],
                "elev": [9.0, 8.5],
                "cond": [25.0, 30.0],
                "min_elev": [0.0, 0.0],
            },
            geometry=[Point(0.5, 0.5).buffer(0.1), Point(1.5, 0.5).buffer(0.1)],
            crs=model.vor.crs,
        ),
    )
    model.vor.gdf_topbtm = gpd.GeoDataFrame(
        {
            "geometry": model.vor.gdf_vorPolys.geometry,
            0: [10.0, 10.0],
            1: [0.0, 0.0],
        },
        geometry="geometry",
        crs=model.vor.crs,
    )
    drn_builder = DRNFromVector(model=model, vor=model.vor, shp_gpkg=drn_path, uid="name")
    drn_spd = drn_builder.from_vector(
        fields={
            "name": "name",
            "height_over_btm": "height",
            "conductance": "cond",
            "layer": "layer",
            "min_elev": "min_elev",
        }
    )
    Drains(model=model, stress_period_data=drn_spd)

    pest = PestProject(
        model=model,
        name="cumberland_slice",
        workspace=workspace / "template",
        start_datetime="2024-01-01",
    )
    pest.add_parameter(
        KPilotPointParameter(
            name="hk",
            source=VectorParameterSource(path=hk_path, value_column="k", zone_column="unit"),
            bounds=(0.25, 4.0),
            bounds_mode="multiplier",
            transform="log",
            pp_spacing=1.0,
            geostruct=ExpGeoStruct(range=10.0, transform="log"),
        )
    )
    pest.add_parameter(
        DrainElevationParameter(
            name="drn_elev",
            source=VectorParameterSource(
                path=drn_path,
                value_column="elev",
                feature_id_column="name",
                group_column="group",
            ),
            bounds=(-2.0, 2.0),
            bounds_mode="absolute",
        )
    )
    pest.add_parameter(
        DrainConductanceParameter(
            name="drn_cond",
            source=VectorParameterSource(
                path=drn_path,
                value_column="cond",
                feature_id_column="name",
                group_column="group",
            ),
            bounds=(0.5, 2.0),
            bounds_mode="multiplier",
            transform="log",
        )
    )
    pest.add_observation(HeadTargetObservationSpec(targets=targets, simulated_values=simulated_values))

    pst = pest.build_pst("slice.pst")
    assert pst is not None
    assert "hk_pp_0000" in pst.parameter_data.index
    assert any(name.startswith("drn_elev_") for name in pst.parameter_data.index)
    assert any(name.startswith("drn_cond_") for name in pst.parameter_data.index)
    obsname = "oname:hds_otype:lst_usecol:obs_a_per:0"
    assert obsname in pst.observation_data.index
    assert pytest.approx(float(pst.observation_data.loc[obsname, "obsval"])) == 10.0
    assert pytest.approx(float(pst.observation_data.loc[obsname, "weight"])) == 1.0


def test_pest_forward_run_applies_k_and_drn_parameters_end_to_end():
    pytest.importorskip("pyemu")
    flopy = pytest.importorskip("flopy")
    workspace = _project_temp_dir("pest_forward_run")
    model, vor = _build_two_cell_pest_forward_model("pest_forward_run", workspace / "model")

    hk_path = _write_gpkg(
        workspace / "hk.gpkg",
        gpd.GeoDataFrame(
            {"name": ["hk_zone"], "unit": ["all"], "k": [5.0]},
            geometry=[Polygon([(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)])],
            crs=model.vor.crs,
        ),
    )
    drn_path = _write_gpkg(
        workspace / "drn.gpkg",
        gpd.GeoDataFrame(
            {
                "name": ["east_drn"],
                "group": ["main"],
                "layer": [1],
                "height": [8.5],
                "elev": [8.5],
                "cond": [20.0],
                "min_elev": [0.0],
            },
            geometry=[Polygon([(1.1, 0.1), (1.9, 0.1), (1.9, 0.9), (1.1, 0.9)])],
            crs=model.vor.crs,
        ),
    )
    drn_builder = DRNFromVector(model=model, vor=vor, shp_gpkg=drn_path, uid="name")
    drn_spd = drn_builder.from_vector(
        fields={
            "name": "name",
            "height_over_btm": "height",
            "conductance": "cond",
            "layer": "layer",
            "min_elev": "min_elev",
        }
    )
    Drains(model=model, stress_period_data=drn_spd)
    success, _ = model.run_simulation()
    assert success is True

    points_path = _write_gpkg(
        workspace / "targets.gpkg",
        gpd.GeoDataFrame(
            {"name": ["OBS_A", "OBS_B"], "layer": [0, 0], "weight": [1.0, 1.0]},
            geometry=[Point(0.5, 0.5), Point(1.5, 0.5)],
            crs=model.vor.crs,
        ),
    )
    targets = HeadTargets(
        locations=points_path,
        values=pd.DataFrame({"per": [0], "OBS_A": [10.0], "OBS_B": [8.75]}),
        time_column="per",
    )

    pest = PestProject(
        model=model,
        name="cumberland_style_forward",
        workspace=workspace / "template",
        start_datetime="2024-01-01",
    )
    pest.add_parameter(
        KPilotPointParameter(
            name="hk",
            source=VectorParameterSource(path=hk_path, value_column="k", zone_column="unit"),
            bounds=(0.25, 4.0),
            bounds_mode="multiplier",
            transform="log",
            pp_spacing=5.0,
            geostruct=ExpGeoStruct(range=10.0, transform="log"),
        )
    )
    pest.add_parameter(
        DrainElevationParameter(
            name="drn_elev",
            source=VectorParameterSource(
                path=drn_path,
                value_column="elev",
                feature_id_column="name",
                group_column="group",
                layer_column="layer",
            ),
            bounds=(-2.0, 2.0),
            bounds_mode="absolute",
        )
    )
    pest.add_parameter(
        DrainConductanceParameter(
            name="drn_cond",
            source=VectorParameterSource(
                path=drn_path,
                value_column="cond",
                feature_id_column="name",
                group_column="group",
                layer_column="layer",
            ),
            bounds=(0.25, 4.0),
            bounds_mode="multiplier",
            transform="log",
        )
    )
    pest.add_observation(HeadTargetObservationSpec(targets=targets))
    pest.build_pst("forward_slice.pst")

    template = workspace / "template"
    forward_run_text = (template / "forward_run.py").read_text(encoding="utf-8")
    assert "def _write_simulation_with_retry" in forward_run_text
    assert "_write_simulation_with_retry(sim)" in forward_run_text
    assert "error removing tmp file" not in forward_run_text
    initial_heads = pd.read_csv(template / "hds_simulated_heads.csv")

    hk_frame = pd.read_csv(template / "hk_pilot_points.csv")
    hk_frame["value"] = 2.0
    hk_frame.to_csv(template / "hk_pilot_points.csv", index=False)

    drn_elev_frame = pd.read_csv(template / "drn_elev_drain_elevation.csv")
    drn_elev_frame["value"] = 0.75
    drn_elev_frame.to_csv(template / "drn_elev_drain_elevation.csv", index=False)

    drn_cond_frame = pd.read_csv(template / "drn_cond_drain_conductance.csv")
    drn_cond_frame["value"] = 0.5
    drn_cond_frame.to_csv(template / "drn_cond_drain_conductance.csv", index=False)

    subprocess.run(
        [sys.executable, "forward_run.py"],
        cwd=template,
        check=True,
        capture_output=True,
        text=True,
    )

    rerun_heads = pd.read_csv(template / "hds_simulated_heads.csv")
    assert not rerun_heads.equals(initial_heads)

    sim = flopy.mf6.MFSimulation.load(sim_ws=str(template), verbosity_level=0)
    gwf = sim.get_model(model.name)
    k_array = np.asarray(gwf.npf.k.array, dtype=float)
    assert np.allclose(k_array[0, :], [10.0, 10.0])

    drn_data = gwf.drn.stress_period_data.data[0]
    assert len(drn_data) == 1
    assert pytest.approx(float(drn_data[0]["elev"])) == 9.25
    assert pytest.approx(float(drn_data[0]["cond"])) == 10.0


def test_pest_run_results_reopen_completed_artifact_and_compare_heads():
    pytest.importorskip("pyemu")
    workspace = _project_temp_dir("pest_run_results")
    artifact_root = workspace / "artifact"
    model, vor = _build_two_cell_pest_forward_model("pest_run_results", artifact_root / "model")

    hk_path = _write_gpkg(
        workspace / "hk_results.gpkg",
        gpd.GeoDataFrame(
            {"name": ["hk_zone"], "unit": ["all"], "k": [5.0]},
            geometry=[Polygon([(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)])],
            crs=model.vor.crs,
        ),
    )
    drn_path = _write_gpkg(
        workspace / "drn_results.gpkg",
        gpd.GeoDataFrame(
            {
                "name": ["east_drn"],
                "group": ["main"],
                "layer": [1],
                "height": [8.5],
                "elev": [8.5],
                "cond": [20.0],
                "min_elev": [0.0],
            },
            geometry=[Polygon([(1.1, 0.1), (1.9, 0.1), (1.9, 0.9), (1.1, 0.9)])],
            crs=model.vor.crs,
        ),
    )
    drn_builder = DRNFromVector(model=model, vor=vor, shp_gpkg=drn_path, uid="name")
    drn_spd = drn_builder.from_vector(
        fields={
            "name": "name",
            "height_over_btm": "height",
            "conductance": "cond",
            "layer": "layer",
            "min_elev": "min_elev",
        }
    )
    Drains(model=model, stress_period_data=drn_spd)
    success, _ = model.run_simulation()
    assert success is True

    points_path = _write_gpkg(
        workspace / "targets_results.gpkg",
        gpd.GeoDataFrame(
            {"name": ["OBS_A", "OBS_B"], "layer": [0, 0], "weight": [1.0, 1.0]},
            geometry=[Point(0.5, 0.5), Point(1.5, 0.5)],
            crs=model.vor.crs,
        ),
    )
    targets = HeadTargets(
        locations=points_path,
        values=pd.DataFrame({"per": [0], "OBS_A": [10.0], "OBS_B": [8.75]}),
        time_column="per",
    )

    pest = PestProject(
        model=model,
        name="pest_run_results",
        workspace=artifact_root / "pest",
        start_datetime="2024-01-01",
    )
    pest.add_parameter(
        KPilotPointParameter(
            name="hk",
            source=VectorParameterSource(path=hk_path, value_column="k", zone_column="unit"),
            bounds=(0.25, 4.0),
            bounds_mode="multiplier",
            transform="log",
            pp_spacing=5.0,
            geostruct=ExpGeoStruct(range=10.0, transform="log"),
        )
    )
    pest.add_parameter(
        DrainConductanceParameter(
            name="drn_cond",
            source=VectorParameterSource(
                path=drn_path,
                value_column="cond",
                feature_id_column="name",
                group_column="group",
                layer_column="layer",
            ),
            bounds=(0.25, 4.0),
            bounds_mode="multiplier",
            transform="log",
        )
    )
    pest.add_observation(HeadTargetObservationSpec(targets=targets))
    pst = pest.build_pst("results_slice.pst")
    metadata_path = artifact_root / "pest" / "simple_modflow_pest_metadata.json"
    assert metadata_path.exists()

    _write_synthetic_par(
        artifact_root / "pest" / "results_slice.par",
        pst.parameter_data.reset_index()[["parnme", "parval1"]],
        {"hk_pp": 2.0, "drn_cond": 0.5},
    )

    results = open_pest_run(artifact_root)
    summary = results.summary()
    assert summary.loc[0, "has_final_parameters"] is True or bool(summary.loc[0, "has_final_parameters"]) is True
    assert summary.loc[0, "n_saved_head_target_sets"] == 1

    baseline_model = results.load_baseline_model()
    calibrated_model = results.load_calibrated_model()
    assert isinstance(results, PestRunResults)
    assert isinstance(baseline_model, LoadedMf6Run)
    assert isinstance(calibrated_model, LoadedMf6Run)
    materialization_path = artifact_root / "pest" / "simple_modflow_calibrated_materialization.json"
    assert materialization_path.exists()
    saved_targets = results.load_head_targets()
    assert isinstance(saved_targets, HeadTargets)
    assert saved_targets.summary().loc[0, "n_rows"] == targets.summary().loc[0, "n_rows"]

    k_frame = results.k_geodata()
    assert not k_frame.empty
    assert np.allclose(k_frame["k_initial"].to_numpy(dtype=float), [5.0, 5.0])
    assert np.allclose(k_frame["k_final"].to_numpy(dtype=float), [10.0, 10.0])
    assert np.allclose(k_frame["k_ratio"].to_numpy(dtype=float), [2.0, 2.0])

    drn_data = calibrated_model.gwf.drn.stress_period_data.data[0]
    assert pytest.approx(float(drn_data[0]["cond"])) == 10.0

    residual_compare = results.compare_head_targets(targets)
    assert "abs_residual_improvement" in residual_compare.columns
    assert not np.allclose(
        residual_compare["sim_head_baseline"].to_numpy(dtype=float),
        residual_compare["sim_head_calibrated"].to_numpy(dtype=float),
    )
    stats_compare = results.compare_head_target_stats(targets)
    assert {"mae_baseline", "mae_calibrated", "rmse_baseline", "rmse_calibrated"}.issubset(stats_compare.columns)
    residual_compare_auto = results.compare_head_targets()
    assert len(residual_compare_auto) == len(residual_compare)
    stats_compare_auto = results.compare_head_target_stats()
    assert stats_compare_auto.loc[0, "n"] == stats_compare.loc[0, "n"]

    review = results.review()
    assert isinstance(review, PestRunReview)
    assert len(review.residual_compare) == len(residual_compare)
    assert not review.k_geodata.empty

    ax_k = results.plot_k()
    assert ax_k.get_title() == "Final K"
    ax_ratio = results.plot_k_ratio()
    assert ax_ratio.get_title() == "K final / K initial"
    ax_period = results.plot_residuals_by_period()
    assert ax_period.get_title() == "Residual MAE by period"
    ax_scatter = results.plot_obs_vs_sim()
    assert ax_scatter.get_title() == "Observed vs simulated heads"
    ax_well = results.plot_well_timeseries("OBS_A")
    assert ax_well.get_title() == "OBS_A"


def test_pest_run_results_raise_clean_error_when_final_par_is_missing():
    pytest.importorskip("pyemu")
    workspace = _project_temp_dir("pest_run_results_missing_par")
    artifact_root = workspace / "artifact"
    model, _ = _build_two_cell_pest_forward_model("pest_missing_par", artifact_root / "model")

    points_path = _write_gpkg(
        workspace / "targets_missing_par.gpkg",
        gpd.GeoDataFrame(
            {"name": ["OBS_A"], "layer": [0], "weight": [1.0]},
            geometry=[Point(0.5, 0.5)],
            crs=model.vor.crs,
        ),
    )
    targets = HeadTargets(
        locations=points_path,
        values=pd.DataFrame({"per": [0], "OBS_A": [10.0]}),
        time_column="per",
    )
    simulated_values = pd.DataFrame({"per": [0], "OBS_A": [10.0]})
    pest = PestProject(
        model=model,
        name="pest_missing_par",
        workspace=artifact_root / "pest",
        start_datetime="2024-01-01",
    )
    pest.add_observation(HeadTargetObservationSpec(targets=targets, simulated_values=simulated_values))
    pest.build_pst("missing_par_slice.pst")

    results = open_pest_run(artifact_root)
    assert results.has_final_parameters is False
    with pytest.raises(FileNotFoundError, match="No final \\.par file"):
        results.load_calibrated_model()


def test_pest_run_results_raise_clean_error_when_saved_targets_are_missing():
    pytest.importorskip("pyemu")
    workspace = _project_temp_dir("pest_run_results_missing_metadata")
    artifact_root = workspace / "artifact"
    model, _ = _build_two_cell_pest_forward_model("pest_missing_metadata", artifact_root / "model")

    points_path = _write_gpkg(
        workspace / "targets_missing_metadata.gpkg",
        gpd.GeoDataFrame(
            {"name": ["OBS_A"], "layer": [0], "weight": [1.0]},
            geometry=[Point(0.5, 0.5)],
            crs=model.vor.crs,
        ),
    )
    targets = HeadTargets(
        locations=points_path,
        values=pd.DataFrame({"per": [0], "OBS_A": [10.0]}),
        time_column="per",
    )
    simulated_values = pd.DataFrame({"per": [0], "OBS_A": [10.0]})
    pest = PestProject(
        model=model,
        name="pest_missing_metadata",
        workspace=artifact_root / "pest",
        start_datetime="2024-01-01",
    )
    pest.add_observation(HeadTargetObservationSpec(targets=targets, simulated_values=simulated_values))
    pst = pest.build_pst("missing_metadata_slice.pst")

    _write_synthetic_par(
        artifact_root / "pest" / "missing_metadata_slice.par",
        pst.parameter_data.reset_index()[["parnme", "parval1"]],
        {},
    )
    metadata_path = artifact_root / "pest" / "simple_modflow_pest_metadata.json"
    assert metadata_path.exists()
    metadata_path.unlink()

    results = open_pest_run(artifact_root)
    with pytest.raises(FileNotFoundError, match="No saved head-target metadata"):
        results.load_head_targets()
    with pytest.raises(FileNotFoundError, match="No saved head-target metadata"):
        results.compare_head_targets()
