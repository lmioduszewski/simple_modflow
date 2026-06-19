from __future__ import annotations

import json
import pickle
import shutil
from pathlib import Path

import numpy as np
import pytest

from examples.mf6.visualization_prt_master_support import (
    MasterExampleConfig,
    build_transient_model,
)
from myflopy.modflow.mf6.interactive_plotting import ModelMapStyle


ROOT = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.canonical


def test_master_example_defaults_to_large_complex_model():
    config = MasterExampleConfig()

    assert config.ncpl >= 10_000
    assert config.nlay >= 4
    assert config.nper >= 6
    assert config.ncpl * config.nlay >= 40_000


def test_canonical_fixture_enforces_complete_model_contract(canonical_model):
    canonical_model.canonical_contract.validate(canonical_model)
    assert set(canonical_model.gwf.package_names) >= canonical_model.canonical_contract.required_packages
    icelltype = np.asarray(canonical_model.gwf.npf.icelltype.get_data())
    assert tuple(int(icelltype[layer].reshape(-1)[0]) for layer in range(4)) == (1, 1, 0, 0)
    conductivity = np.asarray(canonical_model.gwf.npf.k.get_data())
    assert np.median(conductivity[2]) < np.median(conductivity[1]) * 0.01
    assert canonical_model.vor.gdf_vorPolys.area.std() > 0
    areas = canonical_model.vor.gdf_vorPolys.area
    assert areas.loc[canonical_model.get_region_cells("all_lakes")].median() < areas.median() * 0.90
    assert areas.loc[canonical_model.get_region_cells("all_streams")].median() < areas.median() * 0.90


def test_master_example_full_profile_boundary_heads_stay_above_cell_bottoms():
    workspace = ROOT / ".pytest-work" / "master_example_full_boundary_check"
    shutil.rmtree(workspace, ignore_errors=True)
    try:
        model = build_transient_model(workspace, config=MasterExampleConfig(), name="master_bounds")
        bottoms = model.gwf.disv.botm.get_data()
        for records in model.gwf.chd.stress_period_data.get_data().values():
            for record in records:
                layer, cell = record["cellid"]
                assert float(record["head"]) > float(bottoms[layer, cell])
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_master_example_validation_profile_runs_complex_package_topology():
    workspace = ROOT / ".pytest-work" / "master_example_validation"
    shutil.rmtree(workspace, ignore_errors=True)
    try:
        config = MasterExampleConfig.validation()
        model = build_transient_model(workspace, config=config, name="master_valid")

        required = {
            "disv", "npf", "sto", "chd", "ghb", "drn", "rch", "wel",
            "lak", "sfr", "mvr", "uzf", "oc",
        }
        assert required <= set(model.gwf.package_names)
        assert model.gwf.modelgrid.nlay == config.nlay
        assert model.gwf.modelgrid.ncpl == config.ncpl
        assert model.get_region_cells("all_lakes")
        assert model.get_region_cells("all_streams")
        assert model.get_region_cells("uzf_active")

        success, report = model.run_simulation()
        assert success, "\n".join(report[-20:])
        with (model.model_output_folder_path / f"{model.name}.model").open("rb") as file:
            reopened = pickle.load(file)
        assert set(reopened.targets.keys()) == {"drn_flow", "heads", "lake_stage", "sfr_flow", "sfr_stage"}
        heads = model.gwf.output.head()
        assert len(heads.get_kstpkper()) == config.nper * config.steps_per_period
        assert heads.get_data().shape == (config.nlay, 1, config.ncpl)
        assert (model.model_output_folder_path / "master_heads.csv").exists()
        assert (model.model_output_folder_path / "master_sfr_stage.csv").exists()
        assert (model.model_output_folder_path / "master_sfr_flow.csv").exists()
        assert (model.model_output_folder_path / "master_drn.csv").exists()

        simulated = {
            "heads": model.targets.heads.targets.simulated_heads(model),
            "lake_stage": model.targets.lake_stage.targets.simulated_series(model),
            "sfr_stage": model.targets.sfr_stage.targets.simulated_series(model),
            "sfr_flow": model.targets.sfr_flow.targets.simulated_series(model),
            "drn_flow": model.targets.drn_flow.targets.simulated_series(model),
        }
        assert all(not frame.empty for frame in simulated.values())
        assert not model.packages.lak.results.stage.get().empty
        assert not model.packages.sfr.results.q.get().empty
        assert not model.outputs.uzf.ifno_to_cellid.empty
        from myflopy.modflow.mf6.canonical import canonical_head_signals

        head_signals = canonical_head_signals(model)
        assert min(head_signals["spatial_range_by_layer"]) > 10.0
        assert max(head_signals["temporal_range_by_layer"]) > 3.0
        assert head_signals["maximum_drawdown_by_layer"][3] > 1.0
        from myflopy.modflow.mf6.canonical import canonical_feature_signals

        feature_signals = canonical_feature_signals(model)
        assert feature_signals["head_change"]["pond_mound"] > 0.5
        assert feature_signals["head_change"]["unconfined_pumping"] > 1.0
        assert feature_signals["head_change"]["confined_pumping"] > 1.0
        assert feature_signals["seepage_peak"]["unconfined_seepage"] > 0.0
        assert feature_signals["seepage_peak"]["confined_seepage"] > 0.0
        from myflopy.modflow.mf6.canonical import canonical_sfr_signals

        sfr_signals = canonical_sfr_signals(model)
        assert sfr_signals["reach_count"] >= 40
        assert sfr_signals["wet_reach_count"] == sfr_signals["reach_count"]
        assert sfr_signals["minimum_depth"] > 0.05
        assert sfr_signals["gaining_reach_count"] >= 3
        assert sfr_signals["losing_reach_count"] >= 3
        assert sfr_signals["minimum_routed_flow"] > 50_000.0

        slider = model.visualize.head_map_slider_html(
            workspace / "validation_head_slider.html",
            kstpkpers=heads.get_kstpkper()[:2],
            layer=0,
            style=ModelMapStyle(show_contours=False, dpi=40, figsize=(4, 3)),
            embed_frames=False,
        )
        assert slider.frame_count == 2
        assert len(list(slider.frame_directory.glob("*.png"))) == 2
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_master_notebook_documents_full_integration_surface():
    path = ROOT / "examples" / "mf6" / "notebooks" / "master_large_model_visualization_prt_workflow.ipynb"
    notebook = json.loads(path.read_text(encoding="utf-8"))
    source = "\n".join(
        "".join(cell.get("source", []))
        for cell in notebook["cells"]
    )

    for required in (
        "10,000-cell",
        "head_map_slider_html",
        "head_layer_mosaic_slider_html",
        "cross_section_slider_html",
        "plotly_head_map_animation",
        "plotly_cross_section_animation",
        "PRTReleasePoints",
        "open_prt_run",
        "export_3d_html",
    ):
        assert required in source


def test_canonical_master_notebook_set_uses_one_builder():
    notebook_root = ROOT / "examples" / "mf6" / "notebooks"
    # Only the numbered canonical run set (canonical_00..canonical_04) uses the
    # shared mf.build_canonical_model() builder. canonical_model_template.ipynb is
    # an editable scaffold that demonstrates the builder API directly, so it is
    # intentionally excluded here.
    paths = sorted(
        path
        for path in notebook_root.glob("canonical_[0-9]*.ipynb")
        if not path.name.endswith(".executed.ipynb")
    )
    assert len(paths) == 5
    combined = ""
    for path in paths:
        notebook = json.loads(path.read_text(encoding="utf-8"))
        source = "\n".join("".join(cell.get("source", [])) for cell in notebook["cells"])
        assert "build_canonical_model" in source
        assert "notebook_header" in source
        assert sum(cell["cell_type"] == "markdown" for cell in notebook["cells"]) >= 2
        assert sum(cell["cell_type"] == "code" for cell in notebook["cells"]) >= 2
        combined += source
    for required in (
        "infiltration",
        "seepage",
        "observations",
        "visual",
        "PRT",
        "parallel",
        "PEST",
    ):
        assert required.lower() in combined.lower()


def test_canonical_pest_notebook_covers_full_build_run_and_review_workflow():
    path = ROOT / "examples" / "mf6" / "notebooks" / "canonical_04_pest_and_results.ipynb"
    notebook = json.loads(path.read_text(encoding="utf-8"))
    source = "\n".join("".join(cell.get("source", [])) for cell in notebook["cells"])

    assert sum(cell["cell_type"] == "markdown" for cell in notebook["cells"]) >= 9
    assert sum(cell["cell_type"] == "code" for cell in notebook["cells"]) >= 9
    assert not any(
        output.get("output_type") == "error"
        for cell in notebook["cells"]
        for output in cell.get("outputs", [])
    )
    for required in (
        "build_canonical_model",
        "build_and_optionally_run_gold_standard_demo",
        "RUN_PESTPP",
        "NOPTMAX",
        "N_WORKERS",
        "forward_run.py",
        "parameter_data",
        "observation_data",
        "load_head_targets",
        "load_lake_stage_targets",
        "load_sfr_stage_targets",
        "load_sfr_flow_targets",
        "load_drn_flow_targets",
        "review()",
        "plot_obs_vs_sim",
        "plot_residuals_by_period",
        "plot_well_timeseries",
        "export_review",
    ):
        assert required in source
