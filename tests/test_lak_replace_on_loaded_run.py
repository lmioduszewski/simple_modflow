"""End-to-end: replace the LAK package in a model opened from a run folder.

Proves the documented workflow for swapping an advanced package on a
disk-loaded (flopy-backed) model: build + run -> ``mf.load_run`` -> drop the
old LAK -> build a fresh ``mf.lak()`` against the reconstructed grid -> re-run,
and confirms the new package (its starting stage) survives the re-run and lands
in the written input files.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import box

import myflopy as mf
from myflopy.layers import Array
from myflopy.modflow.mf6.canonical_example import (
    CanonicalModelConfig,
    rectangular_voronoi,
)


def _lake_gdf(vor, name: str = "testlake") -> gpd.GeoDataFrame:
    """A single lake polygon over the middle of the grid."""
    xmin, ymin, xmax, ymax = vor.gdf_vorPolys.total_bounds
    cx, cy = (xmin + xmax) / 2, (ymin + ymax) / 2
    w, h = (xmax - xmin) * 0.18, (ymax - ymin) * 0.18
    return gpd.GeoDataFrame(
        {"name": [name], "geometry": [box(cx - w, cy - h, cx + w, cy + h)]},
        crs=vor.crs,
    )


def _build_lake_model(root, *, starting_stage: float, bed_leakance: float):
    """Build + run a tiny 1-layer steady model that contains a LAK package."""
    vor = rectangular_voronoi(CanonicalModelConfig(nrow=10, ncol=10, nlay=1, nper=1))
    gp = vor.get_disv_gridprops()
    ncpl = int(vor.ncpl)
    layers = (
        mf.LayerStack(vor, top=Array(np.full(ncpl, 100.0)), length_units="feet")
        .add("aq", thickness=40.0)
        .build(attach=True)
    )
    ctx = mf.ModelContext(grid=vor, domain=layers.idomain, surfaces=vor.gdf_topbtm)
    lak = mf.lak(
        context=ctx, nper=1, lakes=_lake_gdf(vor), lake_id_field="name",
        starting_stage={"testlake": starting_stage}, lake_bottom={"testlake": 80.0},
        bed_leakance=bed_leakance, connection_modes="automatic",
        status={"testlake": ["ACTIVE"]},
    )
    flow = mf.gwf("lkm", context=ctx, newtonoptions="UNDER_RELAXATION", packages=[
        mf.disv(nlay=1, ncpl=gp["ncpl"], nvert=len(gp["vertices"]),
                vertices=gp["vertices"], cell2d=gp["cell2d"], top=layers.top,
                botm=layers.botm, idomain=layers.idomain, length_units="FEET"),
        mf.ic(strt=95.0), mf.npf(k=10.0), mf.sto(steady_state={0: True}),
        mf.chd(stress_period_data={0: [[(0, 0), 96.0], [(0, ncpl - 1), 90.0]]}),
        lak,
        mf.oc(head_filerecord="lkm.hds", budget_filerecord="lkm.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    sim = mf.SimulationSpec("base", models=(flow,), packages=(
        mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
        mf.ims(models=("lkm",), complexity="MODERATE",
               outer_maximum=100, inner_maximum=100),
    ))
    project = mf.Project(root)
    project.add_simulation(sim)
    return project.run("base", "base")


def _lak_strt(model) -> float:
    """The first lake's starting stage as currently set on the flopy model."""
    return float(model.gwf.lak.packagedata.array["strt"][0])


@pytest.mark.slow
def test_replace_lak_on_loaded_run(tmp_path):
    # 1) build + run a model containing a LAK package.
    run1 = _build_lake_model(tmp_path / "proj", starting_stage=92.0, bed_leakance=0.2)
    assert run1.success

    # 2) reopen it from the run folder; load packages up front so the later
    #    run_simulation() load_all() cannot reload the old LAK over our swap.
    model = mf.load_run(run1.workspace).model()
    model.load_all()
    assert "lak" in [p.lower() for p in model.package_names]
    assert _lak_strt(model) == pytest.approx(92.0)        # the original package

    # 3) build a NEW lak against the reconstructed grid and swap it in.
    vor = model.vor
    ctx = mf.ModelContext(grid=vor, surfaces=vor.gdf_topbtm)
    new_lak = mf.lak(
        context=ctx, nper=model.gwf.modeltime.nper, lakes=_lake_gdf(vor),
        lake_id_field="name", starting_stage={"testlake": 88.0},
        lake_bottom={"testlake": 80.0}, bed_leakance=0.05,
        connection_modes="automatic", status={"testlake": ["ACTIVE"]},
    )
    model.gwf.remove_package("lak")
    new_lak.build(model.gwf)
    assert _lak_strt(model) == pytest.approx(88.0)        # the new package is in place

    # 4) re-run from the workspace.
    success, _ = model.run_simulation()
    assert success
    assert _lak_strt(model) == pytest.approx(88.0)        # swap survived run_simulation()

    # 5) the swap actually persisted to the written input files.
    reread = mf.load_run(run1.workspace).model()
    reread.load_all()
    assert _lak_strt(reread) == pytest.approx(88.0)

    # heads were regenerated on the same grid.
    import flopy.utils as fu
    heads = fu.HeadFile(str(run1.workspace / "lkm.hds")).get_data()
    assert heads.shape[-1] == int(vor.ncpl)
