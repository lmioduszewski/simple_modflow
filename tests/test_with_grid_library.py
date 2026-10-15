"""A move classifies what an `mf.ref` NAMES, which needs the library: `with_grid(library=)`.

`with_grid`/`move_report` classified the `PackageRef` itself. A ref carries no options,
inputs or metadata, so `package_move` answered `unchanged` -- "names no cells" -- for
every one, and `with_grid` never stamped it (its gate was `isinstance(package,
PackageSpec)`). Measured 2026-10-15: a library `mf.ic(strt=<(2, 28) array>)` used as
`mf.ref("ic_regrid")`, moved between two 28-cell meshes, read `unchanged` and BUILT with
the old mesh's heads -- the route the docs taught for saveable remapped heads until
`mf.Heads` (ledger 596), and still teach for a concentration or temperature array.
Ledger 281(a) had fixed exactly this for `Run.resolution` only; the spec side has no
library, so it is handed one.

The same measurement found three more silent paths, all closed here:

* a model on a DEFERRED grid (`mf.grid_ref`, a `GridSpec` recipe) had no mesh in hand,
  so an INLINE per-cell array moved as `unchanged` too -- the limit ledgers 271 and 295
  recorded. The declared grid is resolved at the call now;
* `with_grid(remap=True)` skipped a ref without a word (`_remap_onto` had the same
  gate): the ref read `unchanged` where the inline twin read `remapped`;
* reading `grid=` as the declared mesh let a second move's remap reach a package the
  FIRST move had stamped: the report read `remapped` while the surviving stamp refused
  the build.

User decisions, 2026-10-15: `library=project` on `with_grid` and `move_report`, refused
naming it when a ref or a `grid_ref` cannot be classified without it; resolve the
deferred declared grid at the call; a move that must CHANGE what a ref names puts the
changed copy INLINE in place of the ref.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import box

import myflopy as mf
from myflopy.grid_identity import INLINED_FROM_REF, MOVED_FROM, grid_pin
from myflopy.modflow.mf6.canonical_example import (
    CanonicalModelConfig,
    rectangular_voronoi,
)
from myflopy.specs import ModelContext, PackageRef, PackageSpec
from myflopy.workspace import Project


def _stack():
    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    stack.add("upper", thickness=20.0)
    stack.add("lower", thickness=30.0)
    return stack


@pytest.fixture(scope="module")
def meshes():
    """Two DIFFERENT meshes with the SAME cell count -- the shape no range check sees."""

    a = rectangular_voronoi(CanonicalModelConfig(nrow=4, ncol=7, nlay=2, nper=1))
    b = rectangular_voronoi(CanonicalModelConfig(nrow=7, ncol=4, nlay=2, nper=1))
    assert a.ncpl == b.ncpl == 28
    return a, b


@pytest.fixture
def strt(meshes):
    a, _b = meshes
    return np.linspace(50.0, 90.0, 2 * a.ncpl).reshape(2, a.ncpl)


def _sim(*packages, **gwf):
    model = mf.gwf(
        "m",
        packages=[mf.disv(_stack()), *packages, mf.oc(head_filerecord="m.hds")],
        **gwf,
    )
    return mf.simulation(
        model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])]
    )


def _on_mesh(grid):
    return {"context": ModelContext(grid=grid, domain=np.ones((2, grid.ncpl), dtype=int))}


def _project(tmp_path, name="p", **packages):
    project = Project(tmp_path / name, name=name)
    for key, package in packages.items():
        project.add_package(key, package)
    return project


def _rows(sim, library=None):
    return {row["package"]: row for row in sim.move_report(library=library)}


def _built_strt(run):
    return run.built.built_model("m").model.get_package("ic").strt.array


# --- the headline: a ref is classified as what it NAMES --------------------------- #


_PAYLOADS = {
    "per-cell array": lambda strt: mf.ic(strt=strt),
    "typed cellids": lambda strt: mf.wel(stress_period_data={0: [[(0, 3), -50.0]]}),
}


@pytest.mark.parametrize("payload", sorted(_PAYLOADS))
def test_a_ref_reads_what_its_inline_twin_reads(tmp_path, meshes, strt, payload):
    """THE SHAPE OF THE DECLARATION MUST NOT CHANGE THE VERDICT (ledger 281's invariant).

    The same package, once inline and once through `mf.ref`, is one model -- so the
    report before the move, the report after it and the build must all say the same
    thing about both. Before the fix the ref read `unchanged` at every step and built.
    """

    a, b = meshes
    package = _PAYLOADS[payload](strt)
    extra = [mf.ic(strt=90.0)] if payload == "typed cellids" else []
    project = _project(tmp_path, lib=package)
    inline = _sim(mf.npf(k=1.0), package, *extra, **_on_mesh(a))
    by_ref = _sim(mf.npf(k=1.0), mf.ref("lib"), *extra, **_on_mesh(a))

    assert _rows(inline)[package.name]["how"] == "unverifiable"
    assert _rows(by_ref, project)["lib"]["how"] == "unverifiable"

    moved = by_ref.with_grid(b, library=project)
    entry = next(e for e in moved.model("m").packages if getattr(e, "key", None) == "lib")
    assert isinstance(entry, PackageRef), "the ref stays a ref -- the library is untouched"
    assert entry.moved_from == grid_pin(a)
    assert MOVED_FROM not in (project.package("lib").metadata or {})
    assert _rows(moved)["lib"]["how"] == "unverifiable", "a stamp reads without the library"

    with pytest.raises(ValueError, match=r"declared as mf\.ref\('lib'\), was moved"):
        project.prepare_run("r", moved)
    with pytest.raises(ValueError, match=r"was moved to this grid by `with_grid`"):
        project.prepare_run("r_inline", inline.with_grid(b), overwrite=True)


def test_the_documented_heads_route_no_longer_builds_the_old_heads(tmp_path, meshes, strt):
    """The case that was measured: `ic_regrid` built on the OLD mesh's values, clean.

    The documented route for a saveable remapped array (heads until `mf.Heads`, ledger
    596; a concentration or temperature array still).
    """

    a, b = meshes
    project = _project(tmp_path, ic_regrid=mf.ic(strt=strt))
    sim = _sim(mf.npf(k=1.0), mf.ref("ic_regrid"), **_on_mesh(a))

    moved = sim.with_grid(b, library=project)

    with pytest.raises(ValueError, match=r"name 'ic_regrid' in `grid_free=`"):
        project.prepare_run("r", moved)


def test_a_ref_to_a_scalar_package_is_not_refused(tmp_path, meshes):
    """The other direction: a ref that names no cells must still move freely.

    Stamping every ref "to be safe" would refuse this -- the conservative design that was
    rejected -- and it is the common case: a library `npf(k=1.0)` means the same on any
    mesh.
    """

    a, b = meshes
    project = _project(tmp_path, npf_scalar=mf.npf(k=1.0), ic_flat=mf.ic(strt=90.0))
    sim = _sim(mf.ref("npf_scalar"), mf.ref("ic_flat"), **_on_mesh(a))

    assert {_rows(sim, project)[key]["how"] for key in ("npf_scalar", "ic_flat")} == {
        "unchanged"
    }
    moved = sim.with_grid(b, library=project)
    assert all(
        entry.moved_from is None
        for entry in moved.model("m").packages
        if isinstance(entry, PackageRef)
    )
    run = project.prepare_run("r", moved)
    assert run.built.built_model("m").model.get_package("disv").ncpl.get_data() == 28


def test_grid_free_takes_the_ref_KEY(tmp_path, meshes, strt):
    """A ref is named by the key the caller typed, not by the package it resolves to."""

    a, b = meshes
    project = _project(tmp_path, ic_regrid=mf.ic(strt=strt))
    sim = _sim(mf.npf(k=1.0), mf.ref("ic_regrid"), **_on_mesh(a))

    assert _rows(sim, project)["ic_regrid"]["how"] == "unverifiable"
    rows = sim.move_report(library=project, grid_free=("ic_regrid",))
    assert {r["package"]: r["how"] for r in rows}["ic_regrid"] == "unchanged"
    moved = sim.with_grid(b, library=project, grid_free=("ic_regrid",))
    assert next(e for e in moved.model("m").packages if isinstance(e, PackageRef)).moved_from is None
    project.prepare_run("r", moved)  # asserted mesh-independent, so it builds


# --- without the library: refused, naming it -------------------------------------- #


def test_a_ref_without_the_library_is_refused_naming_it(meshes, strt):
    """Every unclassifiable ref at once -- naming the first sent the caller back per ref."""

    a, b = meshes
    sim = _sim(mf.ref("npf_scalar"), mf.ref("ic_regrid"), **_on_mesh(a))

    for call in (lambda: sim.with_grid(b), lambda: sim.move_report()):
        with pytest.raises(
            ValueError,
            match=r"mf\.ref\('npf_scalar'\), mf\.ref\('ic_regrid'\).*library=project",
        ):
            call()


def test_the_FIRST_grid_needs_no_library(meshes):
    """`mf.gwf(...).with_grid(grid)` gives a model its first grid: nothing to stamp."""

    _a, b = meshes
    model = mf.gwf("m", packages=[mf.ref("npf_scalar")]).with_grid(b)
    assert model.packages == (mf.ref("npf_scalar"),)


def test_library_must_be_something_that_holds_one(meshes):
    a, b = meshes
    with pytest.raises(TypeError, match=r"library=\.\.\.\) takes the Project"):
        _sim(mf.npf(k=1.0), **_on_mesh(a)).with_grid(b, library={"npf": mf.npf(k=1.0)})


# --- a DEFERRED declared grid is resolved at the call ----------------------------- #


@pytest.mark.parametrize("form", ["inline", "ref"])
def test_a_grid_ref_model_is_measured_against_the_mesh_it_names(tmp_path, meshes, strt, form):
    """Ledger 271/295's recorded limit, and it was NOT about refs: an inline array too.

    Measured before the fix: `mf.gwf(..., grid=mf.grid_ref("mesh"))` holding
    `mf.ic(strt=<(2, 28) array>)` -- inline or by ref -- read `unchanged` across a
    same-size move and BUILT on the old mesh's values.
    """

    a, b = meshes
    project = _project(tmp_path, ic_regrid=mf.ic(strt=strt))
    project.add_grid("mesh", mf.GridSpec.from_object(a))
    ic = mf.ref("ic_regrid") if form == "ref" else mf.ic(strt=strt)
    label = "ic_regrid" if form == "ref" else "ic"
    sim = _sim(mf.npf(k=1.0), ic, grid=mf.grid_ref("mesh"))

    with pytest.raises(ValueError, match=r"mf\.grid_ref\('mesh'\).*library=project"):
        sim.with_grid(b)
    assert _rows(sim, project)[label]["how"] == "unverifiable"

    moved = sim.with_grid(b, library=project)
    with pytest.raises(ValueError, match=r"was moved to this grid by `with_grid`"):
        project.prepare_run("r", moved)


def test_a_grid_RECIPE_model_is_measured_too_and_leaves_no_mesh_behind(
    tmp_path, monkeypatch
):
    """A `GridSpec.voronoi` recipe is meshed once, in a directory removed on return.

    No library is needed -- the recipe carries its own sources -- and a move must not
    write a mesh into the working directory (`GridSpec.resolve` with no `workspace=`
    falls back to `<cwd>/_triangle`). The verdict and the stamp are the point: this
    move changes the cell count, so FloPy refuses the wrong-length array at build on its
    own (ledger 282); a same-size move is where only the stamp stands in the way, and
    the `grid_ref` test above builds that one.
    """

    boundary = tmp_path / "dom.gpkg"
    gpd.GeoDataFrame({"n": ["d"]}, geometry=[box(0, 0, 600, 600)], crs=2927).to_file(boundary)

    def recipe(max_area):
        return mf.GridSpec.voronoi(
            boundary=mf.GeoPackageSourceSpec(path=boundary), max_area=max_area
        )

    ncpl = recipe(40_000.0).resolve(workspace=tmp_path / "probe").ncpl
    sim = _sim(mf.npf(k=np.full((2, ncpl), 3.0)), mf.ic(strt=90.0), grid=recipe(40_000.0))
    cwd = tmp_path / "cwd"
    cwd.mkdir()
    monkeypatch.chdir(cwd)

    assert _rows(sim)["npf"]["how"] == "unverifiable"
    moved = sim.with_grid(recipe(12_000.0))

    npf = next(p for p in moved.model("m").packages if p.name == "npf")
    assert MOVED_FROM in npf.metadata
    assert list(cwd.iterdir()) == [], "the declared mesh was resolved somewhere temporary"


# --- a move that CHANGES what a ref names puts the copy inline --------------------- #


def test_remap_carries_a_ref_array_as_an_INLINE_copy(tmp_path, meshes, strt):
    """`remap=True` skipped a ref in silence: `unchanged` where the inline twin read `remapped`."""

    a, b = meshes
    project = _project(tmp_path, ic_regrid=mf.ic(strt=strt))
    sim = _sim(mf.npf(k=1.0), mf.ref("ic_regrid"), **_on_mesh(a))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        moved = sim.with_grid(b, remap=True, library=project)

    ic = next(e for e in moved.model("m").packages if e.name == "ic")
    assert isinstance(ic, PackageSpec), "the remapped copy is not the library's package"
    assert ic.metadata[INLINED_FROM_REF] == "ic_regrid"
    row = _rows(moved)["ic"]
    assert row["how"] == "remapped"
    assert "library package 'ic_regrid' is unchanged" in row["why"]
    np.testing.assert_array_equal(project.package("ic_regrid").options["strt"], strt)

    built = _built_strt(project.prepare_run("r", moved))
    assert built.shape == strt.shape and not np.array_equal(built, strt)


def test_an_eager_gpkg_ref_is_re_declared_as_an_inline_lazy_copy(tmp_path, meshes):
    """Ledger 296's exact move, reached through a ref: inline, `rebuilt`, and it builds.

    Before, the library's eager `.gpkg` reached the move as a bare ref, read `unchanged`,
    and its own pin refused the build -- a report and a build that disagreed.
    """

    a, _b = meshes
    b = rectangular_voronoi(CanonicalModelConfig(nrow=9, ncol=9, nlay=2, nper=1))
    xmin, ymin, xmax, ymax = a.gdf_vorPolys.total_bounds
    path = tmp_path / "bc.gpkg"
    gpd.GeoDataFrame(
        {"name": ["west"], "layer": [1], "head": [95.0], "cond": [2.0]},
        geometry=[box(xmin, ymin, xmin + (xmax - xmin) * 0.3, ymax)], crs=a.crs,
    ).to_file(path, layer="ghb", driver="GPKG")
    ghb = mf.ghb.gpkg(path, layer="ghb", nper=1, head="head", conductance="cond",
                      context=_on_mesh(a)["context"])
    project = _project(tmp_path, **{"ghb/west": ghb})
    sim = _sim(mf.npf(k=1.0), mf.ic(strt=95.0), mf.ref("ghb/west"), **_on_mesh(a))

    assert _rows(sim, project)["ghb/west"]["how"] == "stale"
    moved = sim.with_grid(mf.GridSpec.from_object(b), library=project)

    copy = next(e for e in moved.model("m").packages if e.name == "ghb")
    assert isinstance(copy, PackageSpec) and copy.inputs
    assert _rows(moved)["ghb"]["how"] == "rebuilt"
    records = moved.build_flopy(tmp_path / "ws").simulation.get_model("m").get_package(
        "ghb"
    ).stress_period_data.get_data(0)
    assert len(records) > 0


# --- the stamp is durable and is not overwritten ---------------------------------- #


def test_a_stamped_ref_survives_a_save_and_still_refuses(tmp_path, meshes, strt):
    a, b = meshes
    project = _project(tmp_path, ic_regrid=mf.ic(strt=strt))
    project.add_grid("old", mf.GridSpec.from_object(a))
    project.add_grid("new", mf.GridSpec.from_object(b))
    sim = _sim(mf.npf(k=1.0), mf.ref("ic_regrid"), grid=mf.grid_ref("old"))
    moved = sim.with_grid(mf.grid_ref("new"), library=project)

    assert mf.ref("ic_regrid").to_dict() == {"kind": "PackageRef", "key": "ic_regrid"}
    again = mf.SimulationSpec.from_dict(moved.to_dict())
    refs = [e for e in again.model("m").packages if isinstance(e, PackageRef)]
    assert refs == [PackageRef("ic_regrid", moved_from=grid_pin(a))]

    project.add_simulation(moved)
    report = project.save()
    assert report.skipped == ()
    reloaded = Project.load(tmp_path / "p")
    with pytest.raises(ValueError, match=r"declared as mf\.ref\('ic_regrid'\)"):
        reloaded.prepare_run("r", reloaded.simulation("s"))


def test_run_resolution_reads_the_stamp_on_the_ref(tmp_path, meshes, strt):
    """`Run.resolution` resolves refs from the run's library copy, which carries no stamp."""

    a, b = meshes
    project = _project(tmp_path, ic_regrid=mf.ic(strt=strt))
    moved = _sim(mf.npf(k=1.0), mf.ref("ic_regrid"), **_on_mesh(a)).with_grid(
        b, library=project
    )

    run = project.prepare_run("r", moved, build=False)
    rows = {row["package"]: row["how"] for row in run.resolution}
    assert rows["ic_regrid"] == "unverifiable"


@pytest.mark.parametrize("form", ["inline", "ref"])
def test_a_second_move_keeps_the_FIRST_stamp(tmp_path, meshes, strt, form):
    """The values still mean the mesh they were declared on, and the refusal names it."""

    a, b = meshes
    c = rectangular_voronoi(CanonicalModelConfig(nrow=2, ncol=14, nlay=2, nper=1))
    project = _project(tmp_path, ic_regrid=mf.ic(strt=strt))
    ic = mf.ref("ic_regrid") if form == "ref" else mf.ic(strt=strt)
    once = _sim(mf.npf(k=1.0), ic, **_on_mesh(a)).with_grid(b, library=project)
    twice = once.with_grid(c, library=project)

    entry = twice.model("m").packages[2]
    stamp = entry.moved_from if form == "ref" else entry.metadata[MOVED_FROM]
    assert stamp == grid_pin(a) != grid_pin(b)


def test_a_second_moves_remap_does_not_carry_what_the_first_move_stamped(
    tmp_path, meshes, strt
):
    """Newly reachable, so closed in the same pass.

    A move with no `remap=` stamps a per-cell `ic`; its values still mean mesh A. Once
    the declared grid is read from `grid=` (B, a built grid), a SECOND move's remap had
    a mesh in hand and weighted A's values by B's polygons. Measured without the guard:
    the report read `remapped` and the build REFUSED (the first stamp survives the
    remap) -- not a silent build, but the report and the build disagreeing about the
    same package, which is what `move_report` exists to prevent. It stays stamped.
    """

    a, b = meshes
    c = rectangular_voronoi(CanonicalModelConfig(nrow=2, ncol=14, nlay=2, nper=1))
    once = _sim(mf.npf(k=1.0), mf.ic(strt=strt), **_on_mesh(a)).with_grid(b)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        twice = once.with_grid(c, remap=True)

    assert _rows(twice)["ic"]["how"] == "unverifiable"
    with pytest.raises(ValueError, match=r"'ic' was moved to this grid by `with_grid`"):
        twice.build_flopy(tmp_path / "ws")


def test_nothing_in_the_library_changes(tmp_path, meshes, strt):
    """The library package is shared with every other model and run that names it."""

    a, b = meshes
    project = _project(tmp_path, ic_regrid=mf.ic(strt=strt), npf_scalar=mf.npf(k=1.0))
    library = {key: project.package(key) for key in ("ic_regrid", "npf_scalar")}
    sim = _sim(mf.ref("npf_scalar"), mf.ref("ic_regrid"), **_on_mesh(a))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sim.with_grid(b, library=project)
        sim.with_grid(b, library=project, remap=True)
    sim.move_report(library=project)

    for key, package in library.items():
        assert project.package(key) is package
        assert MOVED_FROM not in (package.metadata or {})
    np.testing.assert_array_equal(library["ic_regrid"].options["strt"], strt)
    assert not any(Path(project.root).glob("_grid*")), "no mesh written into the project"
