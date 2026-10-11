"""`mf.rch(...)` / `mf.evt(...)` declared without `context=` resolve at BUILD time.

The areal builders read two things from the model -- the idomain and a cell count --
and both are on the built discretization, so the declaration does not need a grid to
exist when it is written down. One `mf.rch(nper=1, recharge=6e-4)` therefore builds
correctly on any mesh, and `sim.with_grid(...)` reports it `rebuilt` where it used to
report `stale` (ledger 284).

**With `context=` nothing moves.** The eager branch is untouched, still resolves at
declaration, and still pins. The lazy branch is reached by exactly the declaration that
RAISED before -- which is what keeps every stored digest still, and is the same rule the
array form follows one commit along.

What a lazy recipe may NOT carry is cell identity the caller chose: an explicit `cells=`
list, `name_by_cell=`, a `{cellid: value}` mapping, or a per-cell sequence. Those are
refused at the call, naming `context=` as the way to resolve and pin them.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

import myflopy as mf
from myflopy.package_api import _lazy_areal_factory
from myflopy.specs import GRID_PIN, ModelContext, package_digest, package_move


def _stack():
    return mf.LayerStack(top=mf.Surface.flat(100.0)).add("sand", thickness=50.0).add(
        "clay", thickness=50.0
    )


def _eager_context(grid):
    built = _stack().build(grid)
    return ModelContext(grid=grid, domain=built.idomain, surfaces=built)


def _model(grid, *packages, context=None, discretization=None):
    return mf.gwf(
        "m",
        context=ModelContext(grid=grid) if context is None else context,
        packages=[
            mf.disv(_stack()) if discretization is None else discretization,
            mf.ic(strt=90.0),
            mf.npf(k=10.0),
            mf.oc(head_filerecord="m.hds"),
            *packages,
        ],
    )


def _build(tmp_path, grid, *packages, tag="w", **kwargs):
    model = _model(grid, *packages, **kwargs)
    sim = mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    return sim.build_flopy(tmp_path / tag).simulation.get_model("m")


def _norm(rows):
    """One period's records as comparable tuples, from a recarray or a list of lists.

    The shape has to be normalized because the two sides arrive differently -- FloPy
    hands back a recarray with empty aux fields, the eager spec holds plain lists -- and
    ledger 281 is what happens when a shape difference is allowed to mean something.
    """

    rows = rows.tolist() if hasattr(rows, "tolist") else rows
    return [
        (tuple(row[0]), tuple(float(v) for v in list(row)[1:] if isinstance(v, (int, float))))
        for row in rows
    ]


def _records(package):
    return _norm(package.stress_period_data.get_data()[0])


# --- the payoff ---------------------------------------------------------------


@pytest.mark.parametrize("package", ["rch", "evt"])
def test_one_declaration_equals_the_eager_form_on_two_different_meshes(
    tmp_path, triangulated_grid, triangulated_grid_finer, package
):
    """The commit's own test: lazy equals eager on A and on B.

    Not "builds on both" -- the same records. A lazy form that merely runs is the
    failure mode ledger 249 measured, where stale cell numbers built clean under
    Normal termination.
    """

    lazy = (
        mf.rch(nper=1, recharge=6.0e-4)
        if package == "rch"
        else mf.evt(nper=1, rate=2.0e-3, depth=2.5)
    )

    for tag, grid in (("coarse", triangulated_grid), ("fine", triangulated_grid_finer)):
        context = _eager_context(grid)
        eager = (
            mf.rch(context=context, nper=1, recharge=6.0e-4)
            if package == "rch"
            else mf.evt(context=context, nper=1, rate=2.0e-3, depth=2.5)
        )
        built = _build(tmp_path, grid, lazy, tag=f"{package}_{tag}")
        resolved = _records(built.get_package(package))

        # The eager side is read off its OWN spec, never built onto this model: building
        # it here would REPLACE the lazy package and compare the eager form to itself.
        assert resolved == _norm(eager.options["stress_period_data"][0])
        assert len(resolved) == int(grid.ncpl)


def test_a_lazy_areal_package_moves_with_the_grid(triangulated_grid, triangulated_grid_finer):
    """`stale` before, `rebuilt` after -- which is the whole point of B2.

    `with_grid` reports and stamps; the BUILD refuses. A package that reports `stale`
    cannot be moved without re-declaring it, so this row is the user-visible change.
    """

    sim = mf.simulation(
        _model(triangulated_grid, mf.rch(nper=1, recharge=6.0e-4), mf.evt(nper=1, rate=1e-3, depth=2.0)),
        name="s",
        packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])],
    )

    verdicts = {row["package"]: row["how"] for row in sim.with_grid(triangulated_grid_finer).move_report()}

    assert verdicts["rch"] == "rebuilt"
    assert verdicts["evt"] == "rebuilt"


def test_the_eager_form_still_resolves_and_pins_at_declaration(triangulated_grid):
    """The contract: passing `context=` changes nothing at all."""

    context = _eager_context(triangulated_grid)
    eager = mf.rch(context=context, nper=1, recharge=6.0e-4)

    assert eager.inputs == {}
    assert "stress_period_data" in eager.options
    assert eager.metadata[GRID_PIN]
    assert package_move(eager)[0] == "stale"


# --- the shapes the precondition must accept ----------------------------------


def test_it_resolves_from_the_discretization_when_the_context_has_no_grid(
    tmp_path, triangulated_grid
):
    """The areal builders need an idomain and a cell count, not a grid object.

    Routing through the `.gpkg` precondition verbatim would demand `context.grid` and so
    refuse two model shapes the EAGER form accepts -- a bare `ModelContext(domain=...)`,
    and a model declared with no context at all, which this repo's own examples use.
    A lazy form narrower than the eager one is the divergence ledger 279 made the
    contract against.
    """

    grid = triangulated_grid
    built = _stack().build(grid)
    ncpl = int(grid.ncpl)
    hand_written = mf.disv(
        nlay=2, top=built.top, botm=built.botm, idomain=built.idomain, **grid.get_disv_gridprops()
    )

    from_domain = _build(
        tmp_path, grid, mf.rch(nper=1, recharge=6.0e-4),
        context=ModelContext(domain=built.idomain), discretization=hand_written, tag="dom",
    )

    assert len(_records(from_domain.get_package("rch"))) == ncpl


def test_it_resolves_with_no_context_and_a_discretization_carrying_no_idomain(
    tmp_path, triangulated_grid
):
    """The case with neither route to a cell count.

    `_ArealBuilder` falls back to `context.grid.ncpl` when the domain is None, so with
    no context AND a DISV declared without `idomain=` both are None and it dies five
    frames down with a message naming neither the package nor the fix. An absent IDOMAIN
    means every cell ACTIVE in MF6, so the precondition says that outright rather than
    guessing.
    """

    grid = triangulated_grid
    ncpl = int(grid.ncpl)
    no_idomain = mf.disv(
        nlay=2,
        top=np.full(ncpl, 100.0),
        botm=np.stack([np.full(ncpl, 50.0), np.zeros(ncpl)]),
        **grid.get_disv_gridprops(),
    )

    built = _build(
        tmp_path, grid, mf.rch(nper=1, recharge=6.0e-4),
        context=ModelContext(), discretization=no_idomain, tag="noctx",
    )

    assert len(_records(built.get_package("rch"))) == ncpl


def test_the_et_surface_resolves_from_the_built_discretization(tmp_path, triangulated_grid):
    """`with_grid` DROPS `context.surfaces`, so a lazy EVT has to get the model top from
    the discretization it was built beside -- the same fallback the lazy `.gpkg` added
    for `CellSurfaceOffset` (ledger 279)."""

    grid = triangulated_grid
    built = _stack().build(grid)
    hand_written = mf.disv(
        nlay=2, top=built.top, botm=built.botm, idomain=built.idomain, **grid.get_disv_gridprops()
    )

    model = _build(
        tmp_path, grid,
        mf.evt(nper=1, rate=1e-3, depth=2.0, surface=mf.CellSurfaceOffset("model_top", offset=-2.0)),
        context=ModelContext(domain=built.idomain), discretization=hand_written, tag="sfc",
    )

    surfaces = {row[1][0] for row in _records(model.get_package("evt"))}
    assert surfaces == {98.0}, "model top 100.0 less the 2.0 offset"


def test_a_lazy_package_declared_before_the_discretization_is_refused_by_name(
    tmp_path, triangulated_grid
):
    """Otherwise it dies inside FloPy, several frames from the declaration that caused it."""

    model = mf.gwf(
        "m",
        context=ModelContext(grid=triangulated_grid),
        packages=[mf.rch(nper=1, recharge=6.0e-4), mf.disv(_stack()), mf.ic(strt=90.0), mf.npf(k=10.0)],
    )
    sim = mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])

    with pytest.raises(ValueError, match="declared BEFORE it|before it"):
        sim.build_flopy(tmp_path / "order")


# --- what a recipe may not carry ----------------------------------------------


@pytest.mark.parametrize(
    "declare, because",
    [
        (lambda n: mf.rch(nper=1, recharge=6e-4, cells=[1, 2, 3]), "cells"),
        (lambda n: mf.rch(nper=1, recharge=6e-4, name_by_cell={1: "a"}), "name_by_cell"),
        (lambda n: mf.rch(nper=1, recharge={(0, 1): 6e-4}), "recharge"),
        (lambda n: mf.rch(nper=1, recharge=[6e-4] * n), "recharge"),
        (lambda n: mf.rch(nper=1, recharge={0: [6e-4] * n}), "recharge"),
        (lambda n: mf.evt(nper=1, rate=1e-3, depth=2.0, surface=[10.0] * n), "surface"),
        # Array-likes that are not a `Sequence` were ACCEPTED, read `rebuilt` and made the
        # spec unsaveable (ledger 379).
        (lambda n: mf.rch(nper=1, recharge=np.full(n, 6e-4)), "recharge"),
        (lambda n: mf.rch(nper=1, recharge={0: np.full(n, 6e-4)}), "recharge"),
        (lambda n: mf.rch(nper=1, recharge=__import__("pandas").Series(np.full(n, 6e-4))), "recharge"),
        (lambda n: mf.rch(nper=1, recharge=(6e-4,) * n), "recharge"),
        (lambda n: mf.evt(nper=1, rate=1e-3, depth=np.full(n, 2.0)), "depth"),
    ],
)
def test_cell_identity_the_caller_chose_is_refused_in_a_lazy_recipe(
    triangulated_grid, declare, because
):
    """A recipe is a promise the package can be rebuilt on another mesh, and these
    cannot keep it.

    The explicit `cells=` list is the sharp one: measured before the refusal, a lazy
    `cells=[1, 2, 3]` built the SAME three records at ncpl 28 and at ncpl 92, reported
    `rebuilt`, and warned about nothing. `name_by_cell` does not even error on another
    mesh -- it MISSES and falls back to a generated name.

    All six are still available EAGERLY, so this refuses a combination rather than a
    capability, which is why the message names `context=`.
    """

    with pytest.raises(TypeError, match=because):
        declare(int(triangulated_grid.ncpl))


@pytest.mark.parametrize(
    "declare",
    [
        lambda: mf.rch(nper=1, recharge=6e-4),
        lambda: mf.rch(nper=2, recharge={0: 6e-4, 1: 3e-4}),
        lambda: mf.evt(nper=1, rate=1e-3, depth=2.0),
        lambda: mf.evt(nper=1, rate=1e-3, depth=2.0, surface=mf.CellSurfaceOffset("model_top", offset=-2.0)),
        # A numpy PERIOD key is a period, and a 0-d array a scalar (ledger 379).
        lambda: mf.rch(nper=2, recharge={np.int64(0): 6e-4, 1: 3e-4}),
        lambda: mf.rch(nper=1, recharge=np.array(6e-4)),
    ],
)
def test_the_movable_shapes_are_accepted(declare):
    """The refusing direction: a scalar, a per-period scalar and a `CellSurfaceOffset`
    all mean the same thing on any mesh, so refusing them would make the lazy form
    useless rather than safe."""

    spec = declare()

    assert spec.inputs
    assert package_move(spec)[0] == "rebuilt"


def test_an_incomplete_builder_declaration_still_says_what_is_missing():
    """`context=` became optional; `nper=`/`recharge=` did not."""

    with pytest.raises(TypeError, match="nper"):
        mf.rch(recharge=6e-4)
    with pytest.raises(TypeError, match="depth"):
        mf.evt(nper=1, rate=1e-3)


def test_segmented_et_is_refused_in_the_lazy_branch_too():
    """`mf.evt.array` refuses it, the eager builder refuses it (ledger 283), and a third
    branch that quietly wrote one segment would put the defect straight back."""

    with pytest.raises(ValueError, match="segment"):
        mf.evt(nper=1, rate=1e-3, depth=2.0, nseg=3)


# --- writing it down ----------------------------------------------------------


@pytest.mark.parametrize(
    "declare",
    [
        lambda: mf.rch(nper=2, recharge={0: 6e-4, 1: 3e-4}),
        lambda: mf.evt(nper=1, rate=1e-3, depth=2.0, surface=mf.CellSurfaceOffset("model_top", offset=-2.0)),
    ],
)
def test_a_recipe_round_trips_through_json(declare):
    """The `CellSurfaceOffset` is the one worth asserting: it survives ONLY through
    `inputs`, and elsewhere it comes back a plain dict with nothing raised (ledger 277)."""

    spec = declare()
    revived = type(spec).from_dict(json.loads(json.dumps(spec.to_dict())))

    assert revived.inputs == spec.inputs
    assert package_digest(revived) == package_digest(spec)


def test_the_builder_reference_carries_no_keywords():
    """A keyword on the partial moved the DIGEST, which nothing else here would catch.

    `package_digest` pickles `name` and the builder reference into ONE payload, so they
    share a memo table. A declared spec holds the same INTERNED string in both places
    and pickle emits a back-reference for the second; a spec revived from JSON holds two
    distinct equal strings and pickle writes it twice -- 168 bytes against 172, equal
    values, different digest. `runs.drift()` would then report every reloaded package as
    EDITED, permanently. Measured and checked in as
    `docs/reproducers/digest_pickle_memo.py` (ledger 284).

    `options` and `inputs` are pre-pickled separately, which is why their own strings
    are safe and why this had never bitten before.
    """

    for package in ("rch", "evt"):
        assert _lazy_areal_factory(package).keywords == {}


def test_the_recipe_records_only_what_the_caller_set():
    """`_explicit_options` discipline (ledger 258): a defaulted argument stays out, so
    adding another one later cannot move a stored digest."""

    spec = mf.rch(nper=1, recharge=6e-4)

    assert set(spec.inputs) == {"nper", "recharge"}
    assert mf.rch(nper=1, recharge=6e-4, boundnames=True).inputs["boundnames"] is True


def test_a_saved_project_keeps_a_lazy_areal_model(tmp_path, triangulated_grid):
    """A scalar list form always saved; this asserts the recipe does not break it."""

    sim = mf.simulation(
        _model(triangulated_grid, mf.rch(nper=1, recharge=6e-4), mf.evt(nper=1, rate=1e-3, depth=2.0)),
        name="sim",
        packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])],
    )
    project = mf.Project(tmp_path / "proj", name="p")
    project.add_simulation(sim)

    report = project.save()

    assert report.simulations == 1, report.skipped


def test_a_run_freezes_the_records_it_built_not_the_recipe(tmp_path, triangulated_grid):
    """Freezing the recipe sends a reopened run back to resolve itself again; the
    snapshot exists to remove that dependence (ledger 251/280)."""

    sim = mf.simulation(
        _model(triangulated_grid, mf.rch(nper=1, recharge=6e-4)),
        name="sim",
        packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])],
    )
    project = mf.Project(tmp_path / "proj", name="p")
    project.add_simulation(sim)
    run = project.prepare_run("r1", sim)
    run.write()

    rows = {row["package"]: row for row in run.resolution if row.get("package") == "rch"}
    assert rows["rch"]["how"] == "rebuilt"
    assert rows["rch"]["resolved_records"] == int(triangulated_grid.ncpl)
    assert (tmp_path / "proj" / "runs" / "r1" / ".myflopy" / "frozen.npz").exists()

    reopened = mf.load_run(tmp_path / "proj" / "runs" / "r1")
    assert len(_records(reopened.model("m").fl.rch)) == int(triangulated_grid.ncpl)


# --- the checks a lazy builder must not lose ----------------------------------


@pytest.mark.parametrize("package", ["rch", "evt"])
@pytest.mark.parametrize(
    "option",
    [
        {"print_input": True},
        {"print_flows": True},
        {"auxiliary": ["conc"]},
        {"save_flows": False},
        {"prnt_input": True},
        {"nonsense_xyz": 1},
    ],
)
def test_the_lazy_branch_accepts_exactly_what_the_eager_branch_accepts(package, option):
    """Both directions, which is the only version of this test that checks anything.

    Asserting "a typo is refused" passes with the whole mechanism deleted -- measured:
    with the `functools.partial` bound to `RCHBuilder` instead of to its FloPy class,
    `_accepted_options` narrowed to that builder's eight parameters and the lazy branch
    refused `print_input`, `print_flows` and `auxiliary`, all real MF6 options the eager
    branch accepts, while every typo test still passed. `lazy == eager` is the contract
    (ledger 279) and it is broken by narrowing just as surely as by widening.
    """

    # Surfaces included: EVT's default `CellSurfaceOffset("model_top")` needs them, and
    # an eager side that raises for an unrelated reason would make every row agree by
    # accident.
    eager_context = ModelContext(
        domain=np.ones((2, 6), dtype=int),
        surfaces={"top": np.full(6, 100.0), "botm": np.stack([np.full(6, 50.0), np.zeros(6)])},
    )
    value = {"recharge": 6e-4} if package == "rch" else {"rate": 1e-3, "depth": 2.0}
    helper = getattr(mf, package)

    def outcome(extra):
        try:
            helper(nper=1, **value, **extra)
        except TypeError:
            return "refused"
        return "accepted"

    lazy = outcome(option)
    eager = outcome({**option, "context": eager_context})

    assert lazy == eager, f"{option} is {lazy} lazily and {eager} eagerly"


def test_an_option_the_recipe_already_supplies_collides_by_name():
    """`boundnames` is both a recipe key and a real FloPy option, which is exactly why
    the recipe is FLAT: nested under one key it would be invisible to
    `_check_input_collisions` and would reach the builder twice, dying inside a
    `functools.partial`."""

    spec = mf.rch(nper=1, recharge=6e-4, boundnames=True)

    with pytest.raises(TypeError, match="boundnames"):
        spec.with_options(boundnames=False)


# --- the ARRAY form (ledger 285) ----------------------------------------------


@pytest.mark.parametrize("package", ["rch", "evt"])
def test_the_array_form_resolves_its_layer_index_at_build_time(
    tmp_path, triangulated_grid, triangulated_grid_finer, package
):
    """`irch="top_active"` is the documented default and the one that must stay one:
    MODFLOW 6 with no IRCH applies recharge to layer 1 unconditionally and silently
    skips any column whose layer 1 is inactive (ledger 163).

    Eagerly it resolves to an `(ncpl,)` array, which is why the declaration is mesh-bound
    -- and why `to_dict` refused it, making the whole model unsaveable. Lazily the recipe
    holds the STRING.
    """

    index = "irch" if package == "rch" else "ievt"
    lazy = (
        mf.rch.array(recharge=6.0e-4)
        if package == "rch"
        else mf.evt.array(rate=2.0e-3, depth=2.5)
    )

    assert index in lazy.inputs
    assert package_move(lazy)[0] == "rebuilt"

    for tag, grid in (("coarse", triangulated_grid), ("fine", triangulated_grid_finer)):
        context = _eager_context(grid)
        eager = (
            mf.rch.array(recharge=6.0e-4, context=context)
            if package == "rch"
            else mf.evt.array(rate=2.0e-3, depth=2.5, context=context)
        )
        built = _build(tmp_path, grid, lazy, tag=f"arr_{package}_{tag}")

        resolved = getattr(built.get_package(package), index).get_data()
        resolved = np.asarray(resolved[0] if isinstance(resolved, dict) else resolved)
        assert resolved.size == int(grid.ncpl)
        assert np.array_equal(resolved, np.asarray(eager.options[index]))


def test_the_array_form_is_what_made_a_model_unsaveable(tmp_path, triangulated_grid):
    """The payoff, stated as the defect it removes.

    With a lazy `mf.disv(stack)` so the discretization is not the blocker (ledger 270),
    the DOCUMENTED default array declaration was the next thing that made
    `Project.save()` persist zero simulations -- `irch="top_active"` resolves to an
    `(ncpl,)` array and `to_dict` refuses an ndarray (ledger 193). Same defect, one
    package over.
    """

    lazy = mf.simulation(
        _model(triangulated_grid, mf.rch.array(recharge=6.0e-4)),
        name="sim", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])],
    )
    eager = mf.simulation(
        _model(triangulated_grid, mf.rch.array(recharge=6.0e-4, context=_eager_context(triangulated_grid))),
        name="sim", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])],
    )

    lazy_project = mf.Project(tmp_path / "lazy", name="p")
    lazy_project.add_simulation(lazy)
    eager_project = mf.Project(tmp_path / "eager", name="p")
    eager_project.add_simulation(eager)

    assert lazy_project.save().simulations == 1
    assert eager_project.save().simulations == 0, "the eager array form still cannot be written down"


def test_an_omitted_et_surface_resolves_lazily_rather_than_raising(tmp_path, triangulated_grid):
    """Ledger 283 made an unresolvable surface RAISE rather than write the datum. With no
    context there is now something to resolve it against later, so it defers instead --
    and the set of declarations that raised is exactly the set that became a recipe."""

    lazy = mf.evt.array(rate=1e-3, depth=2.0, ievt=None)
    assert sorted(lazy.inputs) == ["surface"]

    built = _build(tmp_path, triangulated_grid, lazy, tag="evt_surface")
    surface = built.get_package("evt").surface.get_data()
    surface = np.asarray(surface[0] if isinstance(surface, dict) else surface)

    assert np.allclose(surface, 100.0), "the stack's model top, not MODFLOW's datum"


@pytest.mark.parametrize(
    "ievt, written",
    [
        (1, 1),
        (np.int64(1), 1),
        ("ones", 1),
        ({0: "ones"}, 1),
        (0, 0),
    ],
    ids=["int", "numpy int", "array", "per period", "zero"],
)
def test_a_typed_ievt_beside_an_omitted_surface_is_written(tmp_path, triangulated_grid, ievt, written):
    """With ``surface`` left to the build, a TYPED ``ievt`` rode in neither the deferred
    inputs nor the options, and was DROPPED: ``ievt=1`` wrote no IEVT at all and MF6 took
    the ET from layer 1 (ledger 324). It is carried now, and written as given -- a numpy
    integer too, which FloPy itself refused (it subclasses no ``int``)."""

    n = int(triangulated_grid.ncpl)
    ones = np.ones(n, dtype=int)
    value = {0: ones} if isinstance(ievt, dict) else ones if isinstance(ievt, str) else ievt
    spec = mf.evt.array(rate=1e-3, depth=2.0, ievt=value)
    assert sorted(spec.inputs) == ["surface"], "the surface is still left to the build"

    built = _build(tmp_path, triangulated_grid, spec, tag="typed")
    ievt_written = built.get_package("evt").ievt.get_data()
    assert ievt_written is not None, "the typed index was dropped"
    assert (np.asarray(ievt_written[0]).reshape(-1) == written).all()


@pytest.mark.parametrize("make", [
    lambda: mf.rch.array(recharge=6e-4, irch=np.int64(1)),
    lambda: mf.rch.array(recharge=6e-4, irch={0: np.int64(1)}),
    lambda: mf.evt.array(rate=1e-3, depth=2.0, surface=100.0, ievt=np.int64(1)),
], ids=["rch", "rch per period", "evt"])
def test_a_numpy_integer_index_reaches_flopy_as_an_int(tmp_path, triangulated_grid, make):
    """FloPy refused ``np.int64(1)`` as an index ("has no len()") on the EAGER form too:
    numpy's integer scalars subclass no ``int`` (ledger 265's class). Checked with
    ``type(x) is``, never isinstance."""

    spec = make()
    index = next(value for key, value in spec.options.items() if key in ("irch", "ievt"))
    for value in index.values() if isinstance(index, dict) else [index]:
        assert type(value) is int
    built = _build(tmp_path, triangulated_grid, spec, tag="numpy")
    kind = "rch" if "irch" in spec.options else "evt"
    written = getattr(built.get_package(kind), "irch" if kind == "rch" else "ievt").get_data()
    assert (np.asarray(written[0]).reshape(-1) == 1).all()


def test_an_omitted_ievt_beside_an_omitted_surface_is_left_as_it_was():
    """The digest of what was ALREADY right does not move: ``ievt=None`` and the
    ``"top_active"`` default carry nothing new."""

    assert "ievt" not in mf.evt.array(rate=1e-3, depth=2.0, ievt=None).options
    default = mf.evt.array(rate=1e-3, depth=2.0)
    assert "ievt" not in default.options and default.inputs["ievt"] == "top_active"


@pytest.mark.parametrize(
    "declare",
    [
        lambda n: mf.rch.array(recharge=6e-4, irch=0),
        lambda n: mf.rch.array(recharge=6e-4, irch=None),
        lambda n: mf.rch.array(recharge=6e-4, irch=np.zeros(n, dtype=int)),
        lambda n: mf.rch.array(recharge=np.full(n, 6e-4), irch=0),
        lambda n: mf.evt.array(rate=1e-3, depth=2.0, surface=95.0, ievt=0),
    ],
)
def test_a_declaration_that_already_worked_is_left_entirely_alone(triangulated_grid, declare):
    """**The reason the array trigger is narrower than the list form's.**

    A lazy `.gpkg` and a lazy `mf.rch(...)` can key on `context is None` because every
    such declaration RAISED before. These did not: they are legal, working and
    digestible today. Moving their payload into `inputs` would change `package_digest`,
    and `runs.drift()` would report untouched packages as EDITED the first time anyone
    re-ran their declaring script -- the exact failure ledger 250 kept the pin out of the
    digest to avoid.
    """

    spec = declare(int(triangulated_grid.ncpl))

    assert spec.inputs == {}, "this declaration must not become a recipe"
    assert package_digest(spec) is not None


def test_a_run_reports_what_an_array_package_resolved(tmp_path, triangulated_grid):
    """An array package has neither records nor an `ncpl`, so both existing columns said
    nothing about it -- the shape ledger 280 fixed for the list form."""

    sim = mf.simulation(
        _model(triangulated_grid, mf.rch.array(recharge=6e-4)),
        name="sim", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])],
    )
    project = mf.Project(tmp_path / "proj", name="p")
    project.add_simulation(sim)
    run = project.prepare_run("r1", sim)
    run.write()

    row = next(r for r in run.resolution if r.get("package") == "rch")
    assert row["how"] == "rebuilt"
    assert row["resolved_cells"] == int(triangulated_grid.ncpl)
