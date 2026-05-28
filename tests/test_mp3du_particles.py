from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import Point, Polygon

from simple_modflow.modflow.mp3du import prepare_particle_tracking, run_particle_tracking
from simple_modflow.modflow.mp3du.particles import ParticleTrackingInput, ParticleTrackingResult


class _DummyVor:
    def __init__(self, mapped_cells):
        self._mapped_cells = pd.Series(mapped_cells)
        flat_cells: list[int] = []
        for value in self._mapped_cells.tolist():
            if isinstance(value, (list, tuple, np.ndarray, pd.Series)):
                flat_cells.extend(int(cell) for cell in value)
            else:
                flat_cells.append(int(value))
        self.ncpl = (max(flat_cells) + 1) if flat_cells else 0
        self.gdf_vorPolys = gpd.GeoDataFrame(
            geometry=[Point(float(i), 0.0).buffer(0.25) for i in range(self.ncpl)],
            crs="EPSG:2926",
        )

    def get_vor_cells_as_series(self, frame):
        return self._mapped_cells.iloc[: len(frame)]


class _DummyPackage:
    def __init__(self, package_type: str, **attrs):
        self.package_type = package_type
        for key, value in attrs.items():
            setattr(self, key, value)


class _DummyModel:
    def __init__(self, workspace: Path, particle_cells, package_dict=None, idomain=None):
        if idomain is None:
            idomain = np.array([[1, 1, 1]], dtype=int)
        gwf = SimpleNamespace(
            modelgrid=SimpleNamespace(nlay=1, idomain=idomain, ncpl=idomain.shape[-1]),
            package_dict=package_dict or {},
        )
        self.name = "demo_model"
        self.model_output_folder_path = workspace
        self.gwf = gwf
        self.vor = _DummyVor(particle_cells)


def _write_particle_file(path: Path, **columns) -> Path:
    count = len(next(iter(columns.values()))) if columns else 1
    gdf = gpd.GeoDataFrame(
        columns,
        geometry=[Point(float(i), 0.0) for i in range(count)],
        crs="EPSG:2926",
    )
    gdf.to_file(path, driver="GPKG")
    return path


def test_mp3du_particle_field_inference_and_json(tmp_path):
    particle_path = _write_particle_file(
        tmp_path / "particles.gpkg",
        Node=[10, 11],
        TimeRel=[3000, 3000],
        ZLoc=[1.0, 1.0],
        LocName=[101, 102],
    )
    model = _DummyModel(tmp_path, particle_cells=[10, 11])
    pti = ParticleTrackingInput(
        model=model,
        output_path=tmp_path,
        porosities_by_layer=[0.25],
        particle_shp=particle_path,
        direction="BACKWARD",
        simulation_end_time=3650.0,
    )

    assert pti.particle_field_map == {
        "CELLID_ATTR": "Node",
        "TIME_ATTR": "TimeRel",
        "ZLOC_ATTR": "ZLoc",
        "ADDTL_ATTR": ["LocName"],
    }

    json_path = pti.create_json_file()
    config = json.loads(json_path.read_text(encoding="utf-8"))
    pathline = config["SIMULATIONS"][0]["PATHLINE"]
    shapefile = pathline["PARTICLE_START_LOCATIONS"]["SHAPEFILE"]

    assert pathline["DIRECTION"] == "BACKWARD"
    assert pathline["SIMULATION_END_TIME"] == 3650.0
    assert shapefile["CELLID_ATTR"] == "Node"
    assert shapefile["TIME_ATTR"] == "TimeRel"
    assert shapefile["ZLOC_ATTR"] == "ZLoc"
    assert shapefile["ADDTL_ATTR"] == ["LocName"]


def test_mp3du_particle_field_generation_accepts_attrless_vector_input(tmp_path):
    particle_path = _write_particle_file(tmp_path / "particles_missing.gpkg")
    model = _DummyModel(tmp_path, particle_cells=[0])
    pti = ParticleTrackingInput(
        model=model,
        output_path=tmp_path,
        porosities_by_layer=[0.25],
        particle_shp=particle_path,
    )

    assert pti.particle_field_map == {
        "CELLID_ATTR": "P3D_CellID",
        "TIME_ATTR": "TimeRel",
        "ZLOC_ATTR": "ZLoc",
        "ADDTL_ATTR": ["LocName", "Cell0"],
    }
    assert pti.particle_input_path.name == "mp3du_particles_from_vector.shp"


def test_mp3du_start_cell_diagnostics_report_boundary_overlap(tmp_path):
    particle_path = _write_particle_file(
        tmp_path / "particles_diag.gpkg",
        Node=[0, 1, 2],
        TimeRel=[0, 0, 0],
        ZLoc=[0.5, 0.5, 0.5],
    )
    drn_dtype = np.dtype([("cellid", object), ("elev", float), ("cond", float)])
    uzf_dtype = np.dtype([("cellid", object), ("flag", int)])
    drn_data = np.array([((0, 0), 100.0, 1.0)], dtype=drn_dtype)
    uzf_data = np.array([((0, 1), 1)], dtype=uzf_dtype)
    package_dict = {
        "drn": _DummyPackage("drn", stress_period_data=SimpleNamespace(data={0: drn_data})),
        "uzf": _DummyPackage("uzf", packagedata=SimpleNamespace(array=uzf_data)),
    }
    model = _DummyModel(
        tmp_path,
        particle_cells=[0, 1, 2],
        package_dict=package_dict,
        idomain=np.array([[1, 0, 1]], dtype=int),
    )
    pti = ParticleTrackingInput(
        model=model,
        output_path=tmp_path,
        porosities_by_layer=[0.25],
        particle_shp=particle_path,
    )

    summary = pti.get_start_cell_diagnostics()

    assert summary["mapped_particles"] == 3
    assert summary["inactive_particles"] == 1
    assert summary["boundary_packages"]["DRN"]["count"] == 1
    assert summary["boundary_packages"]["UZF"]["count"] == 1
    assert summary["packages_without_iface_override"] == ["UZF"]


def test_mp3du_start_cell_diagnostics_respect_one_based_cellids(tmp_path):
    particle_path = _write_particle_file(
        tmp_path / "particles_one_based.gpkg",
        P3D_CellID=[1, 2, 3],
        TimeRel=[0, 0, 0],
        ZLoc=[0.5, 0.5, 0.5],
    )
    drn_dtype = np.dtype([("cellid", object), ("elev", float), ("cond", float)])
    drn_data = np.array([((0, 0), 100.0, 1.0)], dtype=drn_dtype)
    package_dict = {
        "drn": _DummyPackage("drn", stress_period_data=SimpleNamespace(data={0: drn_data})),
    }
    model = _DummyModel(tmp_path, particle_cells=[0, 1, 2], package_dict=package_dict)
    pti = ParticleTrackingInput(
        model=model,
        output_path=tmp_path,
        porosities_by_layer=[0.25],
        particle_shp=particle_path,
        particle_field_map={"CELLID_ATTR": "P3D_CellID", "TIME_ATTR": "TimeRel", "ZLOC_ATTR": "ZLoc"},
        cellid_index_base=1,
    )

    summary = pti.get_start_cell_diagnostics()

    assert summary["cellid_index_base"] == 1
    assert summary["declared_vs_geometry_mismatches"] == 0
    assert summary["boundary_packages"]["DRN"]["count"] == 1


def test_mp3du_generated_particle_file_from_cell_ids(tmp_path):
    model = _DummyModel(tmp_path, particle_cells=[0, 2, 4], idomain=np.array([[1, 1, 1, 1, 1]], dtype=int))
    pti = ParticleTrackingInput(
        model=model,
        output_path=tmp_path,
        porosities_by_layer=[0.25],
        particle_cells=[0, 2, 4],
    )

    particle_path = pti.particle_input_path
    summary = pti.get_start_cell_diagnostics()

    assert particle_path.exists()
    assert particle_path.suffix == ".shp"
    assert pti.particle_field_map == {
        "CELLID_ATTR": "P3D_CellID",
        "TIME_ATTR": "TimeRel",
        "ZLOC_ATTR": "ZLoc",
        "ADDTL_ATTR": ["LocName", "Cell0"],
    }
    assert summary["particle_file"].endswith("mp3du_particle_cells.shp")
    assert summary["mapped_particles"] == 3
    assert summary["declared_vs_geometry_mismatches"] == 0


def test_mp3du_generated_particle_file_from_polygon_vector(tmp_path):
    polygon_path = tmp_path / "selection.gpkg"
    polygon_gdf = gpd.GeoDataFrame(
        geometry=[Polygon([(0.0, -1.0), (2.5, -1.0), (2.5, 1.0), (0.0, 1.0)])],
        crs="EPSG:2926",
    )
    polygon_gdf.to_file(polygon_path, driver="GPKG")
    model = _DummyModel(tmp_path, particle_cells=[[0, 1, 2]], idomain=np.array([[1, 1, 1]], dtype=int))
    pti = ParticleTrackingInput(
        model=model,
        output_path=tmp_path,
        porosities_by_layer=[0.25],
        particle_shp=polygon_path,
    )

    particle_path = pti.particle_input_path
    summary = pti.get_start_cell_diagnostics()

    assert particle_path.exists()
    assert particle_path.name == "mp3du_particles_from_vector.shp"
    assert summary["cellid_index_base"] == 1
    assert summary["mapped_particles"] == 3
    assert summary["declared_vs_geometry_mismatches"] == 0


def test_prepare_particle_tracking_accepts_cell_ids(monkeypatch, tmp_path):
    model = _DummyModel(tmp_path, particle_cells=[0, 2, 4], idomain=np.array([[1, 1, 1, 1, 1]], dtype=int))
    monkeypatch.setattr(
        ParticleTrackingInput,
        "run",
        lambda self, **kwargs: ParticleTrackingResult(
            json_file=tmp_path / "mp3du_input.json",
            path_file=tmp_path / "mp3du.p3d",
            start_cell_diagnostics={"mapped_particles": len(self.particle_cells or [])},
            diagnostics_file=tmp_path / "mp3du_diagnostics.json",
        ),
    )
    tracker, result = prepare_particle_tracking(
        model=model,
        particles=[0, 2, 4],
        porosity=0.25,
        output_path=tmp_path,
        execute=False,
        convert_output=False,
    )

    assert isinstance(tracker, ParticleTrackingInput)
    assert tracker.particle_cells == [0, 2, 4]
    assert tracker.porosities_by_layer == [0.25]
    assert result["start_cell_diagnostics"]["mapped_particles"] == 3


def test_run_particle_tracking_returns_result_for_vector_input(monkeypatch, tmp_path):
    particle_path = _write_particle_file(
        tmp_path / "particles_run.gpkg",
        Node=[10, 11],
        TimeRel=[3000, 3000],
        ZLoc=[1.0, 1.0],
    )
    model = _DummyModel(tmp_path, particle_cells=[10, 11])
    diagnostics_path = tmp_path / "mp3du_diagnostics.json"
    diagnostics_path.write_text("{}", encoding="utf-8")
    monkeypatch.setattr(
        ParticleTrackingInput,
        "run",
        lambda self, **kwargs: ParticleTrackingResult(
            json_file=tmp_path / "mp3du_input.json",
            path_file=tmp_path / "mp3du.p3d",
            start_cell_diagnostics={"mapped_particles": 2},
            diagnostics_file=diagnostics_path,
        ),
    )
    result = run_particle_tracking(
        model=model,
        particles=particle_path,
        porosity=0.25,
        output_path=tmp_path,
        convert_output=False,
        write_diagnostics=True,
        execute=False,
    )

    assert result["start_cell_diagnostics"]["mapped_particles"] == 2
    assert result["diagnostics_file"].exists()
    assert isinstance(result, ParticleTrackingResult)


def test_mp3du_endpoint_summary_counts_termination_reasons(tmp_path):
    endpoint_path = tmp_path / "endpoint.gpkg"
    endpoint = gpd.GeoDataFrame(
        {
            "PTERM": [
                "Left Model Domain",
                "Internal sink/source: {DRN}",
                "Left Model Domain",
            ]
        },
        geometry=[Point(0.0, 0.0), Point(1.0, 0.0), Point(2.0, 0.0)],
        crs="EPSG:2926",
    )
    endpoint.to_file(endpoint_path, driver="GPKG")

    particle_path = _write_particle_file(
        tmp_path / "particles_endpoint.gpkg",
        Node=[0],
        TimeRel=[0],
        ZLoc=[0.5],
    )
    model = _DummyModel(tmp_path, particle_cells=[0])
    pti = ParticleTrackingInput(
        model=model,
        output_path=tmp_path,
        porosities_by_layer=[0.25],
        particle_shp=particle_path,
    )

    assert pti.summarize_endpoint_output(endpoint_path) == {
        "Internal sink/source: {DRN}": 1,
        "Left Model Domain": 2,
    }


def test_mp3du_result_supports_attribute_and_mapping_access(tmp_path):
    result = ParticleTrackingResult(
        json_file=tmp_path / "mp3du_input.json",
        path_file=tmp_path / "mp3du.p3d",
        diagnostics_file=tmp_path / "mp3du_diagnostics.json",
        start_cell_diagnostics={"mapped_particles": 2},
        endpoint_summary={"Left Model Domain": 1},
    )

    assert result.json_file == tmp_path / "mp3du_input.json"
    assert result["endpoint_summary"] == {"Left Model Domain": 1}
    assert result.to_dict()["diagnostics_file"] == tmp_path / "mp3du_diagnostics.json"


def test_mp3du_run_defaults_execute_and_convert_output(monkeypatch, tmp_path):
    particle_path = _write_particle_file(
        tmp_path / "particles_defaults.gpkg",
        Node=[0],
        TimeRel=[0],
        ZLoc=[0.5],
    )
    model = _DummyModel(tmp_path, particle_cells=[0])
    pti = ParticleTrackingInput(
        model=model,
        output_path=tmp_path,
        porosities_by_layer=[0.25],
        particle_shp=particle_path,
    )

    calls = []

    monkeypatch.setattr(pti, "validate_inputs", lambda **kwargs: None)
    monkeypatch.setattr(pti, "create_modflow_input_files", lambda: None)
    monkeypatch.setattr(pti, "get_start_cell_diagnostics", lambda: {"mapped_particles": 1})
    monkeypatch.setattr(pti, "create_gsf_file", lambda: tmp_path / "grb_to_gsf.json")
    monkeypatch.setattr(pti, "create_path_file", lambda: tmp_path / "mp3du.p3d")
    monkeypatch.setattr(pti, "create_json_file", lambda: tmp_path / "mp3du_input.json")
    monkeypatch.setattr(pti, "create_output_json_file", lambda: tmp_path / "P3DOutput_json.json")
    monkeypatch.setattr(pti, "summarize_endpoint_output", lambda path: {"Left Model Domain": 1})
    monkeypatch.setattr(pti, "write_diagnostics_file", lambda **kwargs: tmp_path / "mp3du_diagnostics.json")
    monkeypatch.setattr(pti, "run_mp3du", lambda path: calls.append(("run_mp3du", path)))
    monkeypatch.setattr(pti, "run_output_conversion", lambda path: calls.append(("run_output_conversion", path)))
    endpoint_path = tmp_path / pti.output_names["ENDPOINT"]
    endpoint_path.write_text("placeholder", encoding="utf-8")

    result = pti.run()

    assert calls == [
        ("run_mp3du", tmp_path / "mp3du_input.json"),
        ("run_output_conversion", tmp_path / "P3DOutput_json.json"),
    ]
    assert isinstance(result, ParticleTrackingResult)
    assert result.endpoint_summary == {"Left Model Domain": 1}


def test_mp3du_subprocess_wrappers_raise_clear_errors(monkeypatch, tmp_path):
    particle_path = _write_particle_file(
        tmp_path / "particles_errors.gpkg",
        Node=[0],
        TimeRel=[0],
        ZLoc=[0.5],
    )
    model = _DummyModel(tmp_path, particle_cells=[0])
    pti = ParticleTrackingInput(
        model=model,
        output_path=tmp_path,
        porosities_by_layer=[0.25],
        particle_shp=particle_path,
    )

    class _Result:
        def __init__(self, returncode, stderr):
            self.returncode = returncode
            self.stderr = stderr
            self.stdout = ""

    monkeypatch.setattr(
        "simple_modflow.modflow.mp3du.particles.subprocess.run",
        lambda *args, **kwargs: _Result(1, "boom"),
    )

    with pytest.raises(RuntimeError, match="writeP3DGSF.exe"):
        pti.create_gsf_file()
    with pytest.raises(RuntimeError, match="mp3du.exe"):
        pti.run_mp3du(tmp_path / "mp3du_input.json")
    with pytest.raises(RuntimeError, match="writeP3DOutput.exe"):
        pti.run_output_conversion(tmp_path / "P3DOutput_json.json")


def test_legacy_prt_names_lazy_load_with_deprecation_warning():
    with pytest.deprecated_call():
        from simple_modflow.modflow.mp3du.particles import PRT

    assert PRT.__module__.endswith("legacy_prt")
