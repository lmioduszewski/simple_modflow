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


# --- the package owns its stack (ledger 270, amended 2026-09-25) ------------ #
#
# A lazy `mf.disv(stack)` held the CALLER'S `LayerStack` whenever nothing in it named a
# file -- an all-inline stack -- so a later `stack.add(...)` (or the planned
# `stack.excavate(...)`) rewrote a package already declared, a library copy of it, and
# the digest `packages.digests()` had cached. A stack naming a file was rebuilt while
# being hash-stamped and did not follow; the two kinds of stack behaved differently.


def _dem(tmp_path):
    rasterio = pytest.importorskip("rasterio")
    from rasterio.transform import from_origin

    path = tmp_path / "ground.tif"
    with rasterio.open(
        path, "w", driver="GTiff", height=12, width=12, count=1, dtype="float32",
        crs="EPSG:2927", transform=from_origin(0, 600, 50, 50),
    ) as handle:
        handle.write(np.full((12, 12), 100.0, dtype="float32"), 1)
    return path


def _held(spec):
    return spec.inputs["layer_stack"]


def test_a_lazy_disv_does_not_follow_an_edit_to_an_inline_stack():
    from myflopy.specs import package_digest

    stack = _stack()
    spec = mf.disv(stack)
    before = package_digest(spec)

    stack.add("base", thickness=5.0)
    stack.replace("upper", thickness=1.0)

    assert _held(spec) is not stack
    assert [layer.name for layer in _held(spec)._layers] == ["upper", "lower"]
    assert _held(spec)._layers[0].surface == _stack()._layers[0].surface
    assert package_digest(spec) == before


def test_the_copy_is_rebuilt_from_the_payload_and_shares_no_record():
    """Rebuilt the way a reload rebuilds it, so it shares no mutable part with the caller.

    `_Layer` is a mutable dataclass. A structural copy (`for_grid`) would hold the same
    records -- isolated from every public edit, which replaces a record rather than
    changing one, but not from a write to one -- and it would drop anything a later
    field adds unless `for_grid` is taught it too, where `to_dict`/`from_dict` must be
    taught it anyway to be saved.
    """

    stack = _stack()

    held = _held(mf.disv(stack))

    assert not any(mine is theirs for mine, theirs in zip(held._layers, stack._layers, strict=True))
    assert held.to_dict() == stack.to_dict()


def test_a_lazy_disv_does_not_follow_an_edit_to_a_stack_naming_a_file(tmp_path):
    """The kind that was already copied, because stamping its hash rebuilt it."""

    stack = _stack(mf.Surface.raster(_dem(tmp_path)))
    spec = mf.disv(stack)

    stack.add("base", thickness=5.0)

    assert _held(spec) is not stack
    assert _held(spec)._top.source_shas, "the stamped branch is what this case is for"
    assert [layer.name for layer in _held(spec)._layers] == ["upper", "lower"]


def test_a_stack_holding_per_cell_arrays_is_copied_too():
    """It cannot be written down, so it is copied structurally rather than rebuilt.

    `for_grid(None)`: a new stack and a new record list, sharing the records -- which
    every public edit replaces rather than changes, so none of them reaches the copy.
    """

    stack = mf.LayerStack(top=mf.Surface.from_array(np.full(28, 100.0)))
    stack.add("upper", thickness=20.0)
    with pytest.raises(ValueError, match="per-cell arrays"):
        stack.to_dict()

    spec = mf.disv(stack)
    stack.add("lower", thickness=30.0)

    assert _held(spec) is not stack
    assert [layer.name for layer in _held(spec)._layers] == ["upper"]


def _as_it_was(stack, name="disv"):
    """The spec the old code built: the caller's own object in `inputs`."""

    from myflopy.package_api import _DISV_BUILD_DEFAULTS, _build_lazy_disv
    from myflopy.specs import PackageSpec

    return PackageSpec(name, _build_lazy_disv, {}, requires=(), inputs={
        "layer_stack": stack, "build_kwargs": dict(_DISV_BUILD_DEFAULTS), "disv_name": name,
    })


@pytest.mark.parametrize("make", [
    lambda: _stack(),
    lambda: mf.LayerStack(top=mf.Surface.flat(100.0) - 2.5).add("a", bottom=60.0),
    lambda: mf.LayerStack(top=mf.Surface.flat(100.0), length_units="meters")
    .add("sand", thickness=30.0, split=(0.3, 0.7), names=("sand_upper", "sand_lower"))
    .add("clay", thickness=10.0, pinch="inactive", min_thickness=2.0)
    .add("till", thickness=6.0, split=3),
    lambda: mf.LayerStack(top=mf.Surface.flat(100.0)).add(
        "fill", bottom=mf.Surface.where(box(0, 0, 250, 600), 100.0, 95.0), pinch="passthrough"
    ),
], ids=["two_units", "algebra", "split_pinch_units", "where_zone"])
def test_holding_a_copy_moves_no_digest(make):
    """The copy is rebuilt from the stack's own payload, so it digests as the original.

    Measured over the whole suite too: 124 declarations where the old code held the
    caller's object (121 stacks, 3 offsets), every copy canonically equal.
    """

    from myflopy.specs import package_digest

    stack = make()

    assert package_digest(mf.disv(stack)) == package_digest(_as_it_was(stack))


def test_a_library_package_does_not_follow_the_callers_stack(tmp_path):
    """`packages.digests()` caches; the cached value stayed while a fresh one moved."""

    from myflopy.specs import package_digest

    stack = _stack()
    project = Project(tmp_path / "p", name="p")
    project.add_package("disv/base", mf.disv(stack))
    cached = dict(project.packages.digests())

    stack.add("base", thickness=5.0)

    assert package_digest(project.packages["disv/base"]) == cached["disv/base"]
    assert dict(project.packages.digests()) == cached
    assert [layer.name for layer in _held(project.packages["disv/base"])._layers] == [
        "upper", "lower",
    ]


def test_a_payload_its_own_class_refuses_is_raised_not_swallowed():
    """It used to return the caller's object, unstamped and uncopied, in silence."""

    from myflopy.specs import stamp_source_hashes

    class _Refuses:
        _SPEC_TOKEN = "$refuses"

        def to_dict(self):
            return {"path": "nowhere.tif"}

        @classmethod
        def from_dict(cls, data):
            raise ValueError("cannot rebuild")

    with pytest.raises(ValueError, match="cannot rebuild"):
        stamp_source_hashes(_Refuses())


# --- ...and so does the EAGER form, and every payload funnel ------------------ #
#
# Review of the change above. `mf.disv(result)` (and `result.to_disv()`) held the
# RESULT'S OWN arrays -- `np.asarray` does not copy -- so an in-place edit of the result
# after the declaration (the shape an excavation of a built result takes) reached the
# declared package and moved its digest. And `_walk_surfaces` rebuilt a token object
# from its payload without taking off the `{"$tuple": ...}` tag the payload carries,
# which `_spec_value` does: `Spread(("a", "b"))` came back with a DICT for a field.


def test_an_eager_disv_does_not_follow_an_in_place_edit_of_its_result(triangulated_grid):
    from myflopy.specs import package_digest

    result = _stack().build(triangulated_grid, attach=False)
    spec = mf.disv(result)
    before = package_digest(spec)

    for key in ("top", "botm", "idomain"):
        assert not np.shares_memory(np.asarray(spec.options[key]), getattr(result, key)), key
    result.idomain[0, :5] = 0
    result.botm[1, :5] -= 3.0
    result.top[:5] += 1.0

    assert int(np.asarray(spec.options["idomain"])[0, :5].min()) == 1
    assert package_digest(spec) == before


def test_a_result_on_a_grid_of_another_size_is_refused_at_the_call(
    triangulated_grid, triangulated_grid_finer
):
    """It wrote the new `ncpl` over the old `top`; FloPy refused it only at build."""

    assert triangulated_grid.ncpl != triangulated_grid_finer.ncpl
    result = _stack().build(triangulated_grid, attach=False)

    with pytest.raises(ValueError, match=rf"{triangulated_grid.ncpl} cells per layer"):
        mf.disv(result, grid=triangulated_grid_finer)
    with pytest.raises(ValueError, match=rf"the grid given has {triangulated_grid_finer.ncpl}"):
        result.to_disv(triangulated_grid_finer)


def test_a_stamped_spread_keeps_its_tuple_and_equals_the_original():
    from myflopy.specs import package_digest, stamp_source_hashes

    spread = mf.Spread(("c0", "c1"))
    assert spread.to_dict()["field"] == {"$tuple": ["c0", "c1"]}, (
        "the payload carrying the tag is what makes this case worth testing"
    )

    held = stamp_source_hashes(spread)
    held_in_inputs = stamp_source_hashes({"conductance": spread})["conductance"]

    for copy in (held, held_in_inputs):
        assert type(copy.field) is tuple and copy.field == ("c0", "c1")
        assert copy == spread
    lazy = mf.PackageSpec("x", lambda target, **_: target, inputs={"conductance": spread})
    assert package_digest(lazy) == package_digest(
        mf.PackageSpec("x", lambda target, **_: target, inputs={"conductance": held})
    )


def test_a_repointed_payload_is_untagged_before_it_is_rebuilt():
    """The third payload-to-object funnel, pinned on a token class of its own.

    No shipped class writes a path AND a tagged tuple into one payload today, so this
    is the rule rather than a live case: a class that does must not get the tag back.
    """

    from myflopy.specs import repoint_source_paths

    received = []

    class _Carrier:
        _SPEC_TOKEN = "$carrier"

        def to_dict(self):
            return {"path": "gis/a.tif", "pair": {"$tuple": [1, 2]}}

        @classmethod
        def from_dict(cls, data):
            received.append(data)
            return cls()

    moved = repoint_source_paths(_Carrier(), lambda path, *, role, notes: "inputs/a.tif")

    assert isinstance(moved, _Carrier)
    assert received == [{"path": "inputs/a.tif", "pair": (1, 2)}]
    assert type(received[0]["pair"]) is tuple
