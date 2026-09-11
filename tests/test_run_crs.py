"""A reopened run is drawn in the CRS it ran in.

``run.json`` recorded no CRS, and a reopened run's models came back through
``load_mf6_run(crs="EPSG:2927")``: a run built on a UTM grid reopened labelled State
Plane, and every map of it was drawn in the wrong place. A plain MF6 workspace --
no manifest -- still defaults to EPSG:2927, on purpose: location-agnostic models
sit at the origin there.
"""

from __future__ import annotations

import json

import numpy as np
import pyproj
import pytest

import myflopy as mf
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus


def _epsg(crs) -> int:
    return pyproj.CRS(crs).to_epsg()


@pytest.fixture(scope="module")
def utm_run(tmp_path_factory):
    n, d = 3, 100.0
    verts = np.array([[i * d, j * d] for j in range(n + 1) for i in range(n + 1)], dtype=float)
    iverts = [
        [j * (n + 1) + i, (j + 1) * (n + 1) + i, (j + 1) * (n + 1) + i + 1, j * (n + 1) + i + 1]
        for j in range(n)
        for i in range(n)
    ]
    xcyc = np.array([[(i + 0.5) * d, (j + 0.5) * d] for j in range(n) for i in range(n)])
    vor = VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc, crs="EPSG:32610")
    model = mf.gwf(
        "flow",
        context=mf.ModelContext(grid=vor),
        packages=[
            mf.disv(nlay=1, top=10.0, botm=[0.0], **vor.get_disv_gridprops()),
            mf.ic(strt=5.0),
            mf.npf(k=1.0),
            mf.chd(stress_period_data={0: [[(0, 0), 6.0], [(0, 8), 4.0]]}),
            mf.oc(),
        ],
    )
    project = mf.Project(tmp_path_factory.mktemp("crs"), name="crs")
    run = project.prepare_run("utm", mf.simulation(model, name="utm"))
    success, report = run.execute(silent=True)
    assert success, report
    return run


def test_the_manifest_records_each_models_crs(utm_run):
    manifest = json.loads((utm_run.workspace / "run.json").read_text())
    assert {name: _epsg(crs) for name, crs in manifest["crs"].items()} == {"flow": 32610}


def test_a_reopened_run_is_drawn_in_the_crs_it_ran_in(utm_run):
    assert _epsg(mf.load_run(utm_run.workspace).model("flow").vor.crs) == 32610


def test_load_mf6_run_reads_it_too_and_an_explicit_crs_still_wins(utm_run):
    assert _epsg(mf.load_mf6_run(utm_run.workspace).vor.crs) == 32610
    assert _epsg(mf.load_mf6_run(utm_run.workspace, crs="EPSG:2927").vor.crs) == 2927


def test_a_reopened_run_keeps_the_block_when_its_manifest_is_rewritten(utm_run):
    reopened = mf.load_run(utm_run.workspace)
    reopened._record_status(reopened.status)
    manifest = json.loads((utm_run.workspace / "run.json").read_text())
    assert {name: _epsg(crs) for name, crs in manifest["crs"].items()} == {"flow": 32610}


def test_a_manifest_without_the_block_still_opens_in_2927(utm_run, tmp_path):
    import shutil

    copy = tmp_path / "old"
    shutil.copytree(utm_run.workspace, copy)
    manifest = json.loads((copy / "run.json").read_text())
    manifest.pop("crs", None)
    (copy / "run.json").write_text(json.dumps(manifest))
    assert _epsg(mf.load_mf6_run(copy).vor.crs) == 2927
