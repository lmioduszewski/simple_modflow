"""`GridSpec.resolve()` writes only where it was told to (ledger 599).

With no ``workspace=`` it meshed into ``<cwd>/_triangle`` -- a side effect in whatever
directory the process started in, shared by every grid resolved there, and reached
silently by ``with_grid(remap=True)`` and ``coverage_on`` given a ``GridSpec``. The user's
rule decided it: raise naming the spelling rather than guess. The library's own callers
name a directory instead -- a temporary one where only the polygons are wanted, the
simulation's ``_grid/<model>`` in a build.

Every test here runs in an EMPTY working directory so that anything left there shows.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import box

import myflopy as mf
from myflopy.modflow.mf6.grid.triangle import TriangleGrid
from myflopy.specs import ModelContext


@pytest.fixture
def boundary(tmp_path):
    path = tmp_path / "inputs" / "dom.gpkg"
    path.parent.mkdir()
    gpd.GeoDataFrame({"n": ["d"]}, geometry=[box(0, 0, 600, 600)], crs=2927).to_file(path)
    return path


@pytest.fixture
def cwd(tmp_path, monkeypatch):
    """An empty working directory: whatever resolving leaves in it is visible."""

    here = tmp_path / "cwd"
    here.mkdir()
    monkeypatch.chdir(here)
    return here


def _spec(boundary, max_area=40_000.0, **options):
    return mf.GridSpec.voronoi(
        boundary=mf.GeoPackageSourceSpec(path=boundary), max_area=max_area, **options
    )


def _stack():
    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    stack.add("upper", thickness=20.0)
    stack.add("lower", thickness=30.0)
    return stack


def _scratch_left():
    return sorted(
        path.name for path in Path(tempfile.gettempdir()).iterdir()
        if path.name.startswith("myflopy-grid-")
    )


# --- the public call --------------------------------------------------------- #


def test_meshing_with_no_workspace_is_refused_naming_it(boundary, cwd):
    with pytest.raises(ValueError, match=r"spec\.resolve\(workspace=\.\.\.\)"):
        _spec(boundary).resolve()

    assert list(cwd.iterdir()) == [], "refused BEFORE anything was written"


@pytest.mark.parametrize("where", ["argument", "options"])
def test_a_NAMED_workspace_is_where_the_files_go(boundary, cwd, tmp_path, where):
    """Both spellings that name one are accepted: the refusal is for naming NONE."""

    ws = tmp_path / "grid"
    if where == "argument":
        grid = _spec(boundary).resolve(workspace=ws)
    else:
        grid = _spec(boundary, model_ws=str(ws)).resolve()

    assert grid.ncpl > 0
    assert (ws / "_triangle.0.poly").exists()
    assert list(cwd.iterdir()) == []


def test_build_False_writes_nothing_and_needs_no_workspace(boundary, cwd):
    """The documented inspect-before-meshing call, kept working."""

    tri = _spec(boundary).resolve(build=False)

    assert isinstance(tri, TriangleGrid)
    assert list(cwd.iterdir()) == []


def test_a_python_spec_with_no_directory_is_refused(cwd):
    """Its builder is handed a workspace that is CREATED: `<cwd>/_grid/<name>`."""

    (cwd / "build_grid.py").write_text(
        "def build_grid(*, project_root, workspace, spec):\n"
        "    class Grid: pass\n"
        "    grid = Grid(); grid.workspace = workspace\n"
        "    return grid\n"
    )
    spec = mf.GridSpec.python("build_grid.py", function="build_grid", name="g")

    with pytest.raises(ValueError, match=r"spec\.resolve\(workspace=\.\.\.\)"):
        spec.resolve()
    assert sorted(path.name for path in cwd.iterdir()) == ["build_grid.py"]

    grid = spec.resolve(project_root=cwd)
    assert grid.workspace == cwd / "_grid" / "g", "a named root is a named directory"


# --- the library's own callers name a directory ------------------------------ #


def test_remap_onto_a_RECIPE_leaves_nothing_behind(boundary, cwd, tmp_path):
    """It meshed the recipe into `<cwd>/_triangle` to read its polygons, and kept it."""

    old = _spec(boundary).resolve(workspace=tmp_path / "old")
    model = mf.gwf("m", grid=old, packages=[mf.disv(_stack()), mf.npf(k=1.0),
                                             mf.ic(strt=np.full((2, old.ncpl), 90.0))])
    before = _scratch_left()

    moved = model.with_grid(_spec(boundary, 12_000.0), remap=True)

    assert next(p for p in moved.packages if p.name == "ic").options["strt"].shape[1] > old.ncpl
    assert list(cwd.iterdir()) == []
    assert _scratch_left() == before, "and the temporary mesh is gone"


def test_coverage_of_a_RECIPE_leaves_nothing_behind(boundary, cwd, tmp_path):
    old = _spec(boundary).resolve(workspace=tmp_path / "old")
    before = _scratch_left()

    report = mf.gwf("m", context=ModelContext(grid=old)).coverage_on(_spec(boundary, 12_000.0))

    assert report.uncovered_area == pytest.approx(0.0, abs=1e-6)
    assert list(cwd.iterdir()) == []
    assert _scratch_left() == before


def test_a_build_naming_no_workspace_meshes_beside_the_simulation(boundary, cwd):
    """`<sim_path>/_grid/<model>`, the documented layout -- not `<cwd>/_triangle`."""

    model = mf.gwf("m", packages=[mf.disv(_stack()), mf.npf(k=1.0), mf.ic(strt=90.0),
                                  mf.oc(head_filerecord="m.hds")]).with_grid(_spec(boundary))
    sim = mf.simulation(model, name="s",
                        packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])

    sim.build_flopy()

    assert (cwd / "_grid" / "m" / "_triangle.0.poly").exists()
    assert not (cwd / "_triangle").exists()


@pytest.mark.parametrize(
    "call, verb",
    [
        (lambda model, ref: model.with_grid(ref, remap=True), "with_grid(remap=True)"),
        (lambda model, ref: model.coverage_on(ref), "coverage_on"),
    ],
    ids=["remap", "coverage"],
)
def test_the_refusal_names_the_call_that_reached_it(boundary, tmp_path, call, verb):
    """It said "coverage needs ..." when the remap reached it."""

    old = _spec(boundary).resolve(workspace=tmp_path / "old")
    model = mf.gwf("m", grid=old, packages=[mf.ic(strt=np.full((2, old.ncpl), 90.0))])

    with pytest.raises(TypeError) as caught:
        call(model, mf.grid_ref("refined"))

    assert str(caught.value).startswith(f"{verb} needs a built grid or a GridSpec")
