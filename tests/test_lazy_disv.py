"""The discretization is a recipe, not a copy of the mesh.

`mf.disv(stack)` on a gridless `LayerStack` used to raise; it now emits a LAZY spec
that resolves the stack against whatever grid the build supplies. That is the first
real remap -- one declaration builds correctly on two different meshes -- and it is
also what finally lets a LAYERED MODEL BE SAVED: `mf.disv(result)` holds resolved
`top`/`botm`/`idomain` arrays, `to_dict` refuses an ndarray (ledger 193), and so
`Project.save()` persisted zero simulations for every layered model.

Two traps this module pins, both of which were live in the first version:

- `ModelSpec.build` passes a LOCAL `context` to the guard and to `BuiltModel`, set
  before the package loop. Without re-reading `model.myflopy_context` afterwards, the
  guard condemns a CORRECT lazy build -- it raises "context domain is shaped (2, 28)
  against ncpl=92" whether or not the DISV just published a correct (2, 92) one.
- The run must freeze what the lazy package BUILT, not its recipe. Freezing the
  declared spec stores a `LayerStack` whose top is still `mf.Raster("gis/ground.tif")`,
  so a reopened run goes back to the files the snapshot exists to make unnecessary
  (ledger 251, 270).
"""

from __future__ import annotations

import shutil

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import box

import myflopy as mf
from myflopy.project.freeze import thaw_simulation
from myflopy.specs import ModelContext
from myflopy.workspace import Project


@pytest.fixture
def boundary(tmp_path):
    """A 600 ft square the deferred grids are cut from."""

    path = tmp_path / "dom.gpkg"
    gpd.GeoDataFrame(
        {"n": ["d"]}, geometry=[box(0, 0, 600, 600)], crs=2927
    ).to_file(path)
    return path


def _stack(top=None):
    """A gridless two-unit declaration."""

    stack = mf.LayerStack(top=top if top is not None else mf.Surface.flat(100.0))
    stack.add("upper", thickness=20.0)
    stack.add("lower", thickness=30.0)
    return stack


def _deferred_model(boundary, max_area, *, stack=None, context=None):
    model = mf.gwf(
        "m",
        context=context,
        packages=[
            mf.disv(stack if stack is not None else _stack()),
            mf.npf(k=1.0), mf.ic(strt=90.0), mf.oc(head_filerecord="m.hds"),
        ],
    )
    return model.with_grid(
        mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec(path=boundary),
                            max_area=max_area)
    )


def _simulation(boundary, max_area, name="s", **kwargs):
    return mf.simulation(
        _deferred_model(boundary, max_area, **kwargs), name=name,
        packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])],
    )


# --- the remap -------------------------------------------------------------- #


@pytest.mark.parametrize("max_area, expected_ncpl", [(40_000.0, 28), (12_000.0, 92)])
def test_one_declaration_builds_on_two_different_grids(
    tmp_path, boundary, max_area, expected_ncpl
):
    """The point of the whole commit, stated as plainly as it can be."""

    built = _simulation(boundary, max_area).build_flopy(tmp_path / f"ws{expected_ncpl}")

    disv = built.built_model("m").model.get_package("disv")
    assert int(disv.ncpl.get_data()) == expected_ncpl


def test_the_lazy_disv_publishes_the_domain_it_resolved(tmp_path, boundary):
    """And `BuiltModel.context` must carry it, not the declared one.

    `ModelSpec.build` sets `context` BEFORE the package loop and passes that local to
    both the guard and `BuiltModel`. Without the re-read this assertion fails and, worse,
    the guard refuses a correct build.
    """

    built = _simulation(boundary, 12_000.0).build_flopy(tmp_path / "ws")

    model = built.built_model("m")
    assert np.asarray(model.context.domain).shape == (2, 92)
    assert np.asarray(model.model.myflopy_context.domain).shape == (2, 92)
    assert model.context.surfaces is not None, "the sampled elevations are published too"


def test_a_declared_domain_that_disagrees_is_refused(tmp_path, boundary):
    """Declaring one domain and resolving another describes two models."""

    grid = mf.GridSpec.voronoi(
        boundary=mf.GeoPackageSourceSpec(path=boundary), max_area=40_000.0
    ).resolve(workspace=tmp_path / "g")
    clashing = ModelContext(grid=grid, domain=np.zeros((2, grid.ncpl), dtype=int))
    model = mf.gwf(
        "m", context=clashing,
        packages=[mf.disv(_stack()), mf.npf(k=1.0), mf.ic(strt=90.0),
                  mf.oc(head_filerecord="m.hds")],
    )
    spec = mf.simulation(model, name="x",
                         packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])

    with pytest.raises(ValueError, match=r"domain declared on model 'm' is not the one"):
        spec.build_flopy(tmp_path / "ws")


def test_a_lazy_disv_with_no_grid_anywhere_says_where_a_grid_comes_from(tmp_path):
    """The old raise was `LayerStack.to_disv() needs a grid`; this one names the routes."""

    model = mf.gwf("m", packages=[mf.disv(_stack()), mf.npf(k=1.0), mf.ic(strt=1.0)])
    spec = mf.simulation(model, name="x",
                         packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])

    with pytest.raises(ValueError, match=r"this model's context carries none"):
        spec.build_flopy(tmp_path / "ws")


# --- the save, which is what could not be done at all ----------------------- #


def test_a_layered_model_saves_reloads_and_rebuilds(tmp_path, boundary):
    """Measured before this: `save()` persisted 0 simulations and `load()` was empty."""

    project = Project(tmp_path / "p", name="p")
    project.add_simulation(_simulation(boundary, 40_000.0))

    assert project.validate() == []
    report = project.save()
    assert report.simulations == 1 and report.skipped == ()

    reloaded = Project.load(tmp_path / "p")
    assert list(reloaded.simulations) == ["s"]
    run = reloaded.prepare_run("r", "s")
    run.build()
    disv = run.built.built_model("m").model.get_package("disv")
    assert int(disv.ncpl.get_data()) == 28


# --- the archive, which is the freeze decision ------------------------------ #


def test_a_lazy_run_rebuilds_with_its_sources_deleted(tmp_path, boundary):
    """The snapshot must hold ARRAYS, not the recipe that points back at the rasters.

    The first version froze the declared spec, so the thawed DISV carried a
    `LayerStack` and no `top`/`botm`/`idomain` at all -- ledger 251's failure mode
    restored. A raster top is used deliberately: a flat one would rebuild from the
    recipe just fine and prove nothing.
    """

    rasterio = pytest.importorskip("rasterio")
    from rasterio.transform import from_origin

    gis = tmp_path / "gis"
    gis.mkdir()
    dem = gis / "ground.tif"
    with rasterio.open(
        dem, "w", driver="GTiff", height=12, width=12, count=1, dtype="float32",
        crs="EPSG:2927", transform=from_origin(0, 600, 50, 50),
    ) as handle:
        handle.write(np.full((12, 12), 100.0, dtype="float32"), 1)

    project = Project(tmp_path / "p", name="p")
    project.add_simulation(_simulation(boundary, 40_000.0, stack=_stack(mf.Surface.raster(dem))))
    run = project.prepare_run("r", "s")
    run.build()
    before = run.built.built_model("m").model.get_package("disv").botm.get_data()

    shutil.rmtree(gis)

    thawed = thaw_simulation(run.workspace)
    rebuilt = thawed.build_flopy(workspace=tmp_path / "again")
    after = rebuilt.built_model("m").model.get_package("disv").botm.get_data()
    assert np.array_equal(before, after)


def test_the_snapshot_holds_the_arrays_not_the_stack(tmp_path, boundary):
    """The mechanism behind the test above, asserted directly."""

    project = Project(tmp_path / "p", name="p")
    project.add_simulation(_simulation(boundary, 40_000.0))
    run = project.prepare_run("r", "s")
    run.build()

    disv = next(p for p in thaw_simulation(run.workspace).models[0].packages
                if p.name == "disv")

    assert set(disv.options) >= {"top", "botm", "idomain", "ncpl"}
    assert not disv.inputs, "the recipe must not survive in place of what it built"


# --- the eight build kwargs ------------------------------------------------- #


def _thin_stack():
    """A unit thinner than `default_min_thickness`, so the build policy shows."""

    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    stack.add("thin", thickness=0.05)
    stack.add("thick", thickness=30.0)
    return stack


def _lazy_arrays(tmp_path, boundary, tag, **kwargs):
    """Build a lazy DISV with ``kwargs`` and return its idomain and botm."""

    model = mf.gwf(
        "m",
        packages=[
            mf.disv(_thin_stack(), **kwargs), mf.npf(k=1.0), mf.ic(strt=90.0),
            mf.oc(head_filerecord="m.hds"),
        ],
    ).with_grid(
        mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec(path=boundary),
                            max_area=40_000.0)
    )
    spec = mf.simulation(model, name="s",
                         packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    disv = spec.build_flopy(tmp_path / tag).built_model("m").model.get_package("disv")
    return np.asarray(disv.idomain.get_data()), np.asarray(disv.botm.get_data())


def test_a_forwarded_build_kwarg_changes_the_arrays(tmp_path, boundary):
    """`mf.disv` forwarded NONE of `to_disv`'s eight, so a stack silently used defaults.

    Asserting the arrays MOVE, not that the build succeeds: a test that built the same
    thing twice and checked it was not None would pass with the forwarding removed.

    This is why the forwarding matters once the DISV resolves lazily. Declare with
    `default_pinch="inactive"`, move to another mesh, and a disv that forwards nothing
    re-derives with `"floor"` -- the thin unit comes back ACTIVE instead of pinched
    out, and nothing says so.
    """

    idomain, botm = _lazy_arrays(tmp_path, boundary, "defaults")
    assert idomain[0].sum() == 28 and botm[0].sum() == pytest.approx(2772.0)

    pinched, _ = _lazy_arrays(tmp_path, boundary, "pinch", default_pinch="inactive")
    assert pinched[0].sum() == 0, "`default_pinch` must reach the resolution"

    _, thicker = _lazy_arrays(tmp_path, boundary, "thick", default_min_thickness=5.0)
    assert thicker[0].sum() == pytest.approx(2660.0), "`default_min_thickness` must reach it"

    _, unreconciled = _lazy_arrays(tmp_path, boundary, "recon", reconcile=False)
    assert unreconciled[0].sum() == pytest.approx(2798.6), "`reconcile` must reach it"


def test_a_build_kwarg_reaches_a_GRID_BOUND_stack_too(triangulated_grid):
    """`mf.disv(bound_stack, ...)` is a third path, and it forwarded nothing either.

    Found by mutation: removing the forwarding from the eager branch killed no test,
    because every kwarg test went through the lazy one.
    """

    def arrays(**kwargs):
        stack = mf.LayerStack(triangulated_grid, top=mf.Surface.flat(100.0))
        stack.add("thin", thickness=0.05)
        stack.add("thick", thickness=30.0)
        spec = mf.disv(stack, **kwargs)
        return np.asarray(spec.options["idomain"]), np.asarray(spec.options["botm"])

    idomain, botm = arrays()
    pinched, _ = arrays(default_pinch="inactive")
    _, thicker = arrays(default_min_thickness=5.0)

    assert idomain[0].sum() == triangulated_grid.ncpl
    assert pinched[0].sum() == 0, "`default_pinch` must reach a grid-bound stack"
    assert thicker[0].sum() != botm[0].sum(), "`default_min_thickness` must reach it"


def test_a_build_kwarg_is_refused_on_an_already_built_result(triangulated_grid):
    """It cannot act there, and a parameter that cannot act is a hidden defect.

    A `LayerBuildResult` already holds its arrays -- they were resolved when
    `stack.build(...)` ran -- so applying `default_pinch` here would describe a
    different model than the one `qc()` showed you.
    """

    stack = _stack()
    result = stack.build(triangulated_grid, attach=True)

    assert mf.disv(result) is not None
    with pytest.raises(TypeError, match=r"cannot apply default_pinch"):
        mf.disv(result, default_pinch="inactive")


def test_a_built_result_still_builds_eagerly(triangulated_grid):
    """`mf.disv(result)` is untouched: making it lazy costs 24 option-read assertions."""

    result = _stack().build(triangulated_grid, attach=True)

    spec = mf.disv(result)

    assert spec.options["ncpl"] == triangulated_grid.ncpl
    assert not spec.inputs, "an eager result carries no recipe"


def test_a_typo_on_a_lazy_disv_is_refused_at_the_call(triangulated_grid):
    """The declaration-time check every other package gets, which the lazy one lost.

    `_accepted_options` inspects a FloPy class, or a wrapper that DECLARES its target
    via `mf6_targets`. `_build_lazy_disv` is a plain function and declared none, so the
    check abstained and `mf.disv(stack, bogus=1)` was accepted at the call -- surfacing
    from inside FloPy at build, a whole `prepare_run` from the mistake. The eager form
    raised the whole time, which is what made the gap invisible.
    """

    with pytest.raises(TypeError, match=r"does not accept 'bogus_option'"):
        mf.disv(_stack(), bogus_option=1)

    # The near-miss path reaches a lazy spec too.
    with pytest.raises(TypeError, match=r"did you mean 'nlay'"):
        mf.disv(_stack(), nlays=3)

    # Controls, in both directions: a real FloPy option and a forwarded build kwarg
    # must still pass, or the check has been bought with a false refusal.
    assert mf.disv(_stack(), length_units="feet").options == {"length_units": "feet"}
    assert mf.disv(_stack(), default_pinch="inactive").inputs["build_kwargs"][
        "default_pinch"
    ] == "inactive"


def test_a_recipe_name_cannot_also_be_passed_as_an_option():
    """`build` splats inputs AND options, so a name in both is always fatal.

    Unchecked it read `_build_lazy_disv() got multiple values for keyword argument
    'layer_stack'` -- the private builder's name, not the call the user wrote. It is
    not a `_check_options` question: `layer_stack` is perfectly ACCEPTABLE to that
    builder; the defect is that the recipe already supplied it.
    """

    with pytest.raises(TypeError, match=r"'layer_stack' is part of this package's recipe"):
        mf.disv(_stack(), layer_stack="nonsense")

    with pytest.raises(TypeError, match=r"'disv_name' is part of this package's recipe"):
        mf.disv(_stack(), disv_name="other")


def test_a_source_in_inputs_may_shadow_the_option_it_supplies():
    """`inputs` carries TWO things, and only one of them can collide.

    A lazy recipe's entries are builder keyword arguments; a duplicate is a mistake.
    A `DataSourceSpec`/`LiteralSource` entry is an UNRESOLVED value for the option of
    the same name -- `inputs["k"]` says where `options["k"]` comes from -- so the
    names matching is the whole point, and `specs.py` says as much where `_json_value`
    handles it ("`PackageSpec.inputs` is exactly where source specs live").

    Refusing it would foreclose a design rather than catch a bug. It is not yet a
    working shape: nothing in `src/` constructs one and nothing resolves one, so a
    build still dies inside FloPy. That resolver is B2's, and this test is what stops
    the collision check from being "tightened" over the top of it in the meantime.
    """

    import dataclasses

    from myflopy.sources import DataSourceSpec, LiteralSource

    for source in (DataSourceSpec(path="x.tif"), LiteralSource({0: 1.0})):
        spec = dataclasses.replace(mf.npf(k=1.0), inputs={"k": source})
        assert spec.options["k"] == 1.0, "the option it supplies is left alone"
        assert spec.inputs["k"] is source

    # ...and the exemption is by TYPE, not by living in `inputs`: a plain value under
    # a name `options` also carries is still the collision this check exists for.
    with pytest.raises(TypeError, match=r"'k' is part of this package's recipe"):
        dataclasses.replace(mf.npf(k=1.0), inputs={"k": 2.0})


def test_a_lazy_disv_is_refused_on_a_simulation():
    """Declaring `mf6_targets` also makes the builder's SCOPE decidable.

    `_package_scope` reads the same signal, so before this a lazy DISV was undecidable
    and `mf.simulation(packages=[...])` waved it through -- to die inside FloPy as
    `Invalid package type`. The eager DISV was refused at the call all along.
    """

    from myflopy.package_api import _build_lazy_disv
    from myflopy.specs import _package_scope

    assert _package_scope(_build_lazy_disv) == "model"

    with pytest.raises(ValueError, match=r"was given model package"):
        mf.simulation(mf.gwf("m"), packages=[mf.disv(_stack())], name="s")


def test_voiding_a_recipe_is_refused_rather_than_given_a_verdict():
    """A spec with its recipe cleared cannot build on ANY mesh, so it must not exist.

    Measured before ledger 275, the failure was worse than an error: `package_move`
    called it `"unchanged"` -- the verdict that means MOVES UNTOUCHED -- so it sailed
    through a grid move and died later as `_build_lazy_disv() missing 2 required
    keyword-only arguments`, naming a private builder.

    Refused at construction rather than labelled, because a label cannot refuse
    anything and the state has no legitimate producer (the house pattern: ledger 208
    refuses a duplicate registration outright).
    """

    import dataclasses

    lazy = mf.disv(_stack())

    with pytest.raises(TypeError, match=r"recipe no longer supplies"):
        dataclasses.replace(lazy, inputs={})

    # A PARTIAL void is the same defect and is caught the same way.
    with pytest.raises(TypeError, match=r"'build_kwargs'"):
        dataclasses.replace(lazy, inputs={"layer_stack": lazy.inputs["layer_stack"]})

    # The message names the builder that actually requires them. `_builder_label_for_error`
    # answers with the four FloPy DISV classes here -- none of which requires these
    # names -- because ledger 272 gave the lazy builder `mf6_targets`.
    with pytest.raises(TypeError, match=r"_build_lazy_disv requires"):
        dataclasses.replace(lazy, inputs={})


def test_an_ordinary_package_pays_nothing_for_that_check():
    """The control: it must refuse nothing that exists.

    Measured: 108 `flopy.mf6.Modflow*` classes have zero required keyword-only
    parameters, and `_build_lazy_disv` is the only builder in the library with any.
    """

    from myflopy.specs import _required_keyword_only

    assert mf.npf(k=1.0) is not None
    assert mf.chd(stress_period_data={0: [[(0, 0), 1.0]]}) is not None
    assert _required_keyword_only(mf.npf(k=1.0).builder) == ()


def test_a_builder_can_ask_for_the_spec_it_came_from():
    """A builder is called as `builder(target, **inputs, **options)` and so learns
    nothing about its own spec -- not its name, not its metadata.

    Fine for a FloPy class, which is told everything through options. Wrong for a
    builder that RESOLVES a recipe: the package it produces has to be NAMED, and the
    only honest source of that name is the `PackageSpec`.

    Freezing the name in the recipe instead was measured worse than useless --
    `fork_package` renames through `_renamed_options`, which rewrites `options` and
    never touches `inputs`, so the frozen name survives a rename and two renamed forks
    collapse into ONE MF6 package under a FloPy "Replacing existing package" warning.
    """

    from myflopy.specs import PackageSpec, _wants_spec, fork_package

    def asks(target, *, _spec, **options):
        return _spec.name

    def plain(target, **options):
        return "no identity"

    assert _wants_spec(asks) and not _wants_spec(plain)

    spec = PackageSpec("chd", asks, {"save_flows": True})
    assert spec.build(None) == "chd"

    # The point of taking it from the SPEC: a rename reaches the built package.
    assert fork_package(spec, name="chd_west").build(None) == "chd_west"


def test_a_builder_that_does_not_ask_is_called_exactly_as_before():
    """Opt-in by signature, so nothing that exists today changes shape."""

    from myflopy.specs import _wants_spec

    assert not _wants_spec(mf.npf(k=1.0).builder), "a FloPy class must not be handed one"
    assert not _wants_spec(mf.chd(stress_period_data={0: [[(0, 0), 1.0]]}).builder)


def test_asking_for_a_spec_is_not_an_incomplete_recipe():
    """`_spec` is filled by `build`, not by the recipe.

    Without the exclusion, ledger 275's completeness check refuses a builder that asks
    for its own identity -- one guard firing on another's mechanism.
    """

    from myflopy.specs import PackageSpec

    def asks(target, *, _spec, layer_stack, **options):
        return _spec.name

    assert PackageSpec("x", asks, inputs={"layer_stack": object()}).build(None) == "x"

    with pytest.raises(TypeError, match=r"no longer supplies 'layer_stack'"):
        PackageSpec("x", asks, inputs={})
