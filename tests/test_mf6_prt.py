from __future__ import annotations

import os
import shutil
import sys
import time
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")

import numpy as np
import pandas as pd
import pytest
import flopy

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from simple_modflow.modflow.mf6.grid.voronoi import VoronoiGridPlus  # noqa: E402
from simple_modflow.modflow.mf6.prt import (  # noqa: E402
    PRTProject,
    PRTReleasePoints,
    PRTRunResults,
    open_prt_run,
)
from simple_modflow.modflow.mf6.canonical_example import representative_cells  # noqa: E402
from simple_modflow.modflow.mf6.simulation.base import SimulationBase  # noqa: E402
from simple_modflow.modflow.mf6.simulation.discretization import (  # noqa: E402
    DisvGrid,
    TemporalDiscretization,
)
from simple_modflow.modflow.mf6.simulation.packages import (  # noqa: E402
    CHD,
    InitialConditions,
    KFlow,
    OutputControl,
    Storage,
)


def _workspace(name: str) -> Path:
    path = ROOT / ".pytest-work" / f"{name}_{time.time_ns()}"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _flow_model(name: str, workspace: Path):
    verts = np.array(
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [2.0, 0.0], [2.0, 1.0]],
        dtype=float,
    )
    vor = VoronoiGridPlus(
        verts=verts,
        iverts=[[0, 3, 2, 1], [1, 2, 5, 4]],
        xcyc=np.array([[0.5, 0.5], [1.5, 0.5]], dtype=float),
    )
    model = SimulationBase(name=name, mf_folder_path=workspace, vor=vor, nper=1)
    DisvGrid(vor=vor, model=model, top=[10.0, 10.0], bottom=[[0.0, 0.0]], nlay=1, idomain=[[1, 1]])
    TemporalDiscretization(model=model, per_len=1.0, num_steps=1, multiplier=1.0)
    InitialConditions(model=model, vor=vor, nlay=1, strt=[10.0, 9.0])
    KFlow(model=model, k=[10.0, 10.0], k33_vert=[1.0, 1.0], save_specific_discharge=True)
    Storage(model=model, sto_steady={0: True})
    CHD(model=model, stress_period_data={0: [((0, 0), 10.0), ((0, 1), 9.0)]})
    OutputControl(model=model)
    return model


def test_prt_release_points_from_cells_and_points():
    workspace = _workspace("prt_release_points")
    try:
        model = _flow_model("prt_release_points", workspace)
        from_cells = PRTReleasePoints.from_cells(model, [0, 1], local_z=0.25)
        from_points = PRTReleasePoints.from_points(model, [(0.5, 0.5), (1.5, 0.5)], local_z=0.75)

        assert from_cells.packagedata[0] == (0, (0, 0), 0.5, 0.5, 0.25)
        assert from_points.packagedata[1] == (1, (0, 1), 1.5, 0.5, 0.75)
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_prt_release_points_support_structured_grid_cells():
    workspace = _workspace("prt_structured_release_points")
    try:
        model = SimulationBase(name="structured_flow", mf_folder_path=workspace, nper=1)
        flopy.mf6.ModflowGwfdis(
            model.gwf,
            nlay=1,
            nrow=2,
            ncol=2,
            delr=10.0,
            delc=10.0,
            top=10.0,
            botm=0.0,
        )
        from_cells = PRTReleasePoints.from_cells(model, [(0, 1)], local_z=0.5)
        from_points = PRTReleasePoints.from_points(model, [(15.0, 15.0)], local_z=0.5)

        assert from_cells.packagedata[0][1] == (0, 0, 1)
        assert len(from_points.packagedata[0][1]) == 3
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_prt_project_copies_disv_timing_and_writes_inputs():
    workspace = _workspace("prt_build")
    try:
        model = _flow_model("prt_build_flow", workspace)
        project = PRTProject(
            model,
            workspace=workspace / "prt",
            release_points=PRTReleasePoints.from_cells(model, [0]),
            porosity=0.25,
        )
        project.write()

        assert len(project.name) <= 16
        assert model.particle_tracking is model.particle_tracking
        bound_project = model.particle_tracking.prt(
            workspace=workspace / "bound_prt",
            release_points=PRTReleasePoints.from_cells(model, [1]),
        )
        assert isinstance(bound_project, PRTProject)
        assert project.prt.disv is not None
        assert project.prt.disv.ncpl.get_data() == 2
        assert project.prt.prp.nreleasepts.get_data() == 1
        fmi_data = project.prt.fmi.packagedata.get_data()
        assert fmi_data.shape[0] == 3
        assert set(fmi_data["flowtype"]) == {"GWFHEAD", "GWFBUDGET", "GWFGRID"}
        assert (workspace / "prt" / "mfsim.nam").exists()
        assert (workspace / "prt" / f"{project.name}.prp").exists()
        assert (workspace / "prt" / f"{project.name}.fmi").exists()
        assert not project.prt.prp.extend_tracking.get_data()

        extended = PRTProject(
            model,
            workspace=workspace / "extended_prt",
            release_points=PRTReleasePoints.from_cells(model, [0]),
            extend_tracking=True,
        )
        assert extended.prt.prp.extend_tracking.get_data()
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_prt_project_rejects_explicit_name_over_mf6_limit():
    workspace = _workspace("prt_name_limit")
    try:
        model = _flow_model("flow", workspace)
        with pytest.raises(ValueError, match="16"):
            PRTProject(
                model,
                workspace=workspace / "prt",
                name="this_prt_name_is_too_long",
                release_points=PRTReleasePoints.from_cells(model, [0]),
            )
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_open_prt_run_reads_track_csv_and_terminal_points():
    workspace = _workspace("prt_open")
    try:
        model = _flow_model("prt_open_flow", workspace)
        prt_workspace = workspace / "prt"
        prt_workspace.mkdir()
        track_path = prt_workspace / "tracking.trk.csv"
        pd.DataFrame(
            {
                "kper": [0, 0, 0],
                "kstp": [0, 0, 0],
                "imdl": [1, 1, 1],
                "iprp": [1, 1, 1],
                "irpt": [0, 0, 1],
                "ilay": [1, 1, 1],
                "icell": [1, 2, 2],
                "izone": [0, 0, 0],
                "istatus": [0, 0, 0],
                "x": [0.5, 1.0, 1.5],
                "y": [0.5, 0.5, 0.5],
                "z": [5.0, 5.0, 5.0],
                "t": [0.0, 1.0, 1.0],
                "trelease": [0.0, 0.0, 0.0],
                "ireason": [0, 3, 3],
                "name": ["p0", "p0", "p1"],
            }
        ).to_csv(track_path, index=False)

        result = open_prt_run(model, prt_workspace)
        assert isinstance(result, PRTRunResults)
        assert result.engine == "mf6-prt"
        assert len(result.pathlines) == 3
        first_load = result.pathlines
        assert len(result.terminal_points) == 2
        fig, ax = result.plot_map()
        assert fig is not None
        assert ax.get_title() == "Particle pathlines"

        pd.DataFrame({"ireason": [3], "x": [0.0], "y": [0.0]}).to_csv(track_path, index=False)
        assert result.pathlines is first_load
        assert len(result.refresh().pathlines) == 1
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_prt_run_requires_completed_flow_outputs():
    workspace = _workspace("prt_missing_flow_outputs")
    try:
        model = _flow_model("prt_missing_flow_outputs", workspace)
        project = PRTProject(
            model,
            workspace=workspace / "prt",
            release_points=PRTReleasePoints.from_cells(model, [0]),
        )

        with pytest.raises(FileNotFoundError, match="completed GWF head and budget outputs"):
            project.run()
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_real_gwf_to_prt_run_and_shared_scene():
    workspace = _workspace("prt_real_run")
    try:
        model = _flow_model("prt_real_flow", workspace)
        success = model.run_simulation()
        assert success

        project = PRTProject(
            model,
            workspace=workspace / "prt",
            release_points=PRTReleasePoints.from_cells(model, [0], local_z=0.5),
            porosity=0.25,
            stop_at_weak_sink=False,
        )
        result = project.run(silent=True)

        assert result.success
        assert result.track_csv_path.exists()
        assert not result.pathlines.empty
        scene = result.scene(off_screen=True)
        try:
            assert len(scene.meshes) >= 2
        finally:
            scene.plotter.close()
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


@pytest.mark.canonical
def test_canonical_prt_run_stops_at_flow_end_and_writes_pathlines(
    canonical_run,
    canonical_config,
    tmp_path,
):
    project = PRTProject(
        canonical_run,
        workspace=tmp_path / "canonical_prt",
        release_points=PRTReleasePoints.from_cells(
            canonical_run,
            representative_cells(canonical_config)["releases"],
        ),
        porosity=0.25,
    )

    result = project.run(silent=True)

    assert result.success
    assert result.track_csv_path.stat().st_size > 0
    assert not result.pathlines.empty
