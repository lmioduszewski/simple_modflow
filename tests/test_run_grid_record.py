"""A run records the grid each model was built on, and ``runs.summary()`` shows it (ledger 316).

The same ``GridSpec`` resolves to a different grid before and after ledger 316's repair,
and nothing else a run records covers a grid -- no digest, and ``built_from`` hashes
package sources only -- so a run whose packages are all lazy rebuilt on the repaired grid
without a word (measured on the L: CHD inflow +3.81%).
"""

from __future__ import annotations

import json

import pytest
from shapely.geometry import Polygon

import myflopy as mf

L_SHAPE = Polygon([(0, 0), (1000, 0), (1000, 400), (400, 400), (400, 1000), (0, 1000)])


@pytest.fixture(scope="module")
def boundary(tmp_path_factory):
    import geopandas as gpd

    path = tmp_path_factory.mktemp("l") / "boundary.gpkg"
    gpd.GeoDataFrame({"name": ["domain"]}, geometry=[L_SHAPE], crs=2927).to_file(path)
    return path


def _stack():
    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    stack.add("aquifer", thickness=30.0)
    return stack


def _project(root, boundary, **grid_options):
    project = mf.Project(root, name="p")
    model = mf.gwf("m", packages=[mf.disv(_stack()), mf.npf(k=1.0), mf.ic(strt=90.0)]).with_grid(
        mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec(path=boundary), max_area=5000.0, **grid_options)
    )
    project.add_simulation(mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])]))
    project.prepare_run("r", "s")
    return project


def test_a_run_records_its_grid_and_what_the_repair_did(tmp_path, boundary):
    from myflopy.grid_identity import grid_pin

    project = _project(tmp_path / "p", boundary)
    run = project.runs["r"]
    record = run.grids["m"]
    assert record["tessellation"] == "dual"
    assert record["repaired"] == 1 and record["rebuilt"] is False
    resolved = mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec(path=boundary), max_area=5000.0).resolve(
        workspace=tmp_path / "check")
    assert record["pin"] == grid_pin(resolved) != resolved.tessellation_report.legacy_pin
    manifest = json.loads((project.runs_dir / "r" / "run.json").read_text())
    assert manifest["metadata"]["grids"]["m"]["pin"] == record["pin"]
    assert project.runs.summary().set_index("key").loc["r", "grid"] == "m: 1 cell repaired"


def test_the_escape_hatch_says_so(tmp_path, boundary):
    project = _project(tmp_path / "p", boundary, tessellation="flopy")
    assert project.runs["r"].grids["m"]["tessellation"] == "flopy"
    assert project.runs.summary().set_index("key").loc["r", "grid"] == "m: FloPy's cells"


def test_a_run_from_before_the_record_reads_unrecorded(tmp_path, boundary):
    """Built before runs recorded a grid: it may sit on FloPy's cells, holes and all."""

    project = _project(tmp_path / "p", boundary)
    manifest_path = project.runs_dir / "r" / "run.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["metadata"].pop("grids")
    manifest_path.write_text(json.dumps(manifest))
    assert project.refresh_runs().summary().set_index("key").loc["r", "grid"] == "unrecorded"


def test_a_model_without_a_grid_has_nothing_to_say(tmp_path):
    project = mf.Project(tmp_path / "p", name="p")
    model = mf.gwf("flow", packages=[mf.dis(nlay=1, nrow=2, ncol=2, delr=1.0, delc=1.0, top=10.0, botm=0.0),
                                     mf.ic(strt=1.0), mf.npf(k=1.0)])
    project.add_simulation(mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])]))
    project.prepare_run("r", "s")
    manifest_path = project.runs_dir / "r" / "run.json"
    manifest = json.loads(manifest_path.read_text())
    assert manifest["metadata"]["grids"] == {}
    assert project.runs.summary().set_index("key").loc["r", "grid"] == "—"
    manifest["metadata"].pop("grids")
    manifest_path.write_text(json.dumps(manifest))
    assert project.refresh_runs().summary().set_index("key").loc["r", "grid"] == "—"
