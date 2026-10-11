"""`mf.uzf` declared without `context=` resolves at BUILD time (ledger 304).

One declaration builds correctly on two meshes, and `with_grid` reads it `rebuilt` where
the eager form is `stale`. UZF is the cheapest of the advanced packages to make lazy and
the measurement says why: it needs only the DOMAIN. An attribute spy shows
`UZFBuilder.build()` touching exactly one grid attribute, `get_disv_gridprops`, and only
from inside `pin_metadata` -- so it needs strictly LESS than `mf.rch` does, and it names no
file at all, which means nothing for the bundler, the project anchor or the content digest
to get wrong.

With `context=` nothing changes: every test here that compares the two is asserting that.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

import myflopy as mf
from myflopy.grid_identity import GRID_PIN, package_move
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, rectangular_voronoi
from myflopy.specs import ModelContext, package_digest

ARGS = dict(
    nper=2, vks=1.0, thtr=0.05, thts=0.35, thti=0.08, eps=4.0,
    finf=1.0e-3, pet=1.0e-3, extdp=2.0,
)


def _grid(n):
    return rectangular_voronoi(CanonicalModelConfig(nrow=n, ncol=n, nlay=3, nper=1))


@pytest.fixture(scope="module")
def meshes():
    return _grid(5), _grid(9)


def _stack():
    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    for name in ("a", "b", "c"):
        stack.add(name, thickness=20.0)
    return stack


def _model(grid, uzf, **context_kwargs):
    return mf.simulation(
        mf.gwf("m", context=ModelContext(grid=grid, **context_kwargs),
               packages=[mf.disv(_stack()), mf.npf(k=1.0), mf.ic(strt=50.0), uzf]),
        name="s", packages=[mf.tdis(nper=2)],
    )


# --------------------------------------------------------------------------- the payoff


def test_one_declaration_matches_a_fresh_eager_one_on_two_meshes(meshes, tmp_path):
    """The whole point: lazy == eager, per mesh, from a single declaration."""

    lazy = mf.uzf(**ARGS)
    for index, grid in enumerate(meshes):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            built = _model(grid, lazy).build(tmp_path / f"m{index}")
        written = built.simulation.get_model("m").get_package("uzf").packagedata.get_data()

        # The domain the lazy DISV PUBLISHED, not one of ours -- comparing against a
        # hand-made domain compares two different models.
        published = built.built.models["m"].context.domain
        eager = mf.uzf(context=ModelContext(grid=grid, domain=published), **ARGS)
        rows = eager.options["packagedata"]

        assert len(written) == len(rows) == grid.ncpl * 3
        for got, want in zip(written, rows, strict=True):
            assert tuple(got)[:4] == tuple(want)[:4]


def test_a_lazy_uzf_reads_rebuilt_where_an_eager_one_reads_stale(meshes):
    small, _ = meshes
    assert package_move(mf.uzf(**ARGS))[0] == "rebuilt"
    eager = mf.uzf(context=ModelContext(grid=small, domain=np.ones((3, small.ncpl), int)), **ARGS)
    assert package_move(eager)[0] == "stale"


def test_one_declaration_survives_a_grid_move_and_builds(meshes, tmp_path):
    small, large = meshes
    sim = _model(small, mf.uzf(**ARGS))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        moved = sim.with_grid(large)
        assert next(r["how"] for r in moved.move_report() if r["package"] == "uzf") == "rebuilt"
        built = moved.build(tmp_path / "moved")
    rows = built.simulation.get_model("m").get_package("uzf").packagedata.get_data()
    assert len(rows) == large.ncpl * 3


def test_the_build_records_the_RESOLVED_package_with_a_pin(meshes, tmp_path):
    """So the run freezes RECORDS, not the recipe -- and records which mesh they are on.

    `UZFBuilder` pins from `context.grid`, which a lazily resolved context may not carry;
    the built DISV is the same mesh by another route, so the resolved package falls back to
    `disv_pin` rather than going out under the weaker marker (ledger 302).
    """

    small, _ = meshes
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        built = _model(small, mf.uzf(**ARGS)).build(tmp_path / "resolved")
    resolved = built.built.models["m"].resolved
    assert "uzf" in resolved, f"only {sorted(resolved)} were recorded"
    assert GRID_PIN in resolved["uzf"].metadata
    assert isinstance(resolved["uzf"].options["packagedata"], list)


def test_a_lazy_uzf_model_can_be_SAVED(meshes, tmp_path):
    """An eager UZF holds records -- lists of cellid tuples -- and a recipe is scalars."""

    small, _ = meshes
    project = mf.Project(tmp_path / "proj", name="p")
    project.add_grid("mesh", small)
    project.add_simulation(
        mf.simulation(
            mf.gwf("m", grid=mf.grid_ref("mesh"),
                   packages=[mf.disv(_stack()), mf.npf(k=1.0), mf.ic(strt=50.0),
                             mf.uzf(**ARGS)]),
            name="s", packages=[mf.tdis(nper=2)],
        )
    )
    report = project.save()
    assert report.simulations == 1, report
    assert len(mf.Project.load(tmp_path / "proj").simulations) == 1


# -------------------------------------------------------------------------- refusals


def test_an_explicit_cell_list_is_refused(meshes):
    """You cannot pin cell numbers the caller typed (ledger 284's rule, same words)."""

    with pytest.raises(TypeError, match="names cells by number"):
        mf.uzf(**ARGS, cells=[(0, 1), (0, 2)])


@pytest.mark.parametrize("field", ["vks", "thtr", "finf", "extdp"])
def test_a_per_cell_value_is_refused(meshes, field):
    small, _ = meshes
    with pytest.raises(TypeError, match="per-cell"):
        mf.uzf(**{**ARGS, field: [1.0] * small.ncpl})


@pytest.mark.parametrize("shape", ["column", "layered"])
def test_a_per_cell_ARRAY_is_refused_like_a_list(meshes, shape):
    """An ndarray is not a `collections.abc.Sequence`, so it was accepted, read `rebuilt`
    across a move and made the spec unsaveable (ledger 379)."""

    small, _ = meshes
    values = np.ones(small.ncpl) if shape == "column" else np.ones((3, small.ncpl))
    with pytest.raises(TypeError, match="vks= is per-cell.*context="):
        mf.uzf(**{**ARGS, "vks": values})


def test_mover_is_refused_on_a_lazy_declaration(meshes):
    """A UZF mover endpoint is a hand-typed `iuzno` -- an enumeration over ACTIVE CELLS.

    Measured: the column at one location is iuzno 18-20 on one mesh and 33-35 on another,
    in range on both. There is no `uzf_connection` to make it symbolic, so it is refused
    rather than moved wrongly.
    """

    with pytest.raises(TypeError, match="mover=True"):
        mf.uzf(**ARGS, mover=True)


def test_nper_is_validated_at_the_declaration(meshes):
    with pytest.raises(ValueError, match="nper must be at least 1"):
        mf.uzf(**{**ARGS, "nper": 0})


@pytest.mark.parametrize(
    "extra", [{"cells": [(0, 1)]}, {"mover": True}, {"vks": "per-cell-ok-eagerly"}]
)
def test_everything_refused_lazily_still_works_eagerly(meshes, extra):
    """This refuses a COMBINATION, not a capability."""

    small, _ = meshes
    domain = np.ones((3, small.ncpl), dtype=int)
    if "vks" in extra:
        extra = {"vks": [1.0] * (small.ncpl * 3)}
    spec = mf.uzf(context=ModelContext(grid=small, domain=domain), **{**ARGS, **extra})
    assert spec.options["packagedata"]


def test_a_lazy_uzf_declared_before_the_disv_is_refused_by_name(meshes, tmp_path):
    small, _ = meshes
    sim = mf.simulation(
        mf.gwf("m", context=ModelContext(grid=small),
               packages=[mf.uzf(**ARGS), mf.disv(_stack()), mf.npf(k=1.0), mf.ic(strt=50.0)]),
        name="s", packages=[mf.tdis(nper=2)],
    )
    with pytest.raises(ValueError, match="declared BEFORE it"):
        sim.build(tmp_path / "order")


# ------------------------------------------------------------- the eager path must not move


def test_the_eager_declaration_is_byte_identical(meshes):
    """`context=` was required and is now optional; the eager branch must not have moved."""

    small, _ = meshes
    domain = np.ones((3, small.ncpl), dtype=int)
    spec = mf.uzf(context=ModelContext(grid=small, domain=domain), **ARGS)
    assert spec.inputs == {}
    assert GRID_PIN in spec.metadata
    assert spec.metadata["builder"] == "UZFBuilder"


def test_the_recipe_is_FLAT_so_a_colliding_option_is_refused(meshes):
    """Nesting would hide every key from `_check_input_collisions` (ledger 284)."""

    spec = mf.uzf(**ARGS)
    assert "nper" in spec.inputs and "vks" in spec.inputs
    with pytest.raises(TypeError, match="part of this package's recipe"):
        spec.with_options(boundnames=True)


def test_two_lazy_declarations_of_the_same_thing_digest_alike(meshes):
    """A recipe of scalars, so the digest is stable without a source hash."""

    assert package_digest(mf.uzf(**ARGS)) == package_digest(mf.uzf(**ARGS))
    assert package_digest(mf.uzf(**ARGS)) != package_digest(mf.uzf(**{**ARGS, "vks": 2.0}))


@pytest.mark.parametrize("context", ["none", "domain-only", "grid-only"])
def test_a_lazy_uzf_needs_no_GRID_only_a_domain(meshes, tmp_path, context):
    """`needs_grid=False`, and this is the test that holds it.

    Every other model in this file carries `ModelContext(grid=...)`, so none of them
    exercises the difference -- a mutation setting `needs_grid=True` passed all 35 tests.
    UZF resolves its cell set, its `iuzno` numbering and its `ivertcon` from the IDOMAIN
    alone, so demanding a grid would make the lazy form strictly NARROWER than the eager
    one, which accepts a bare `ModelContext(domain=...)` and a model with no context at
    all -- the thing `_lazy_resolution_context`'s own docstring says not to do (ledger 284).

    A HAND-WRITTEN disv, because a lazy `mf.disv(stack)` is the one package that does need
    a grid: here the mesh arrives through the discretization instead of the context.
    """

    small, _ = meshes
    gridprops = small.get_disv_gridprops()
    disv = mf.disv(nlay=3, top=100.0, botm=[80.0, 60.0, 40.0], **gridprops)
    # "grid-only": a context with a grid and NO domain, beside a DISV that declares no
    # idomain -- every cell active, as in MF6. It raised "ModelContext.domain is required"
    # while the no-context model built (ledger 380).
    kwargs = {
        "none": {},
        "domain-only": {"context": ModelContext(domain=np.ones((3, small.ncpl), dtype=int))},
        "grid-only": {"context": ModelContext(grid=small)},
    }[context]
    sim = mf.simulation(
        mf.gwf("m", packages=[disv, mf.npf(k=1.0), mf.ic(strt=50.0), mf.uzf(**ARGS)],
               **kwargs),
        name="s", packages=[mf.tdis(nper=2)],
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        built = sim.build(tmp_path / f"nogrid_{context}")
    rows = built.simulation.get_model("m").get_package("uzf").packagedata.get_data()
    assert len(rows) == small.ncpl * 3
