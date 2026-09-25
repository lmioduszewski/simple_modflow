"""Array-form (``READASARRAYS``) recharge and ET (ledger 163).

`mf.rch` and `mf.evt` produced only the LIST form, so a whole-grid recharge field
was written one record per cell. `mf.rch.array` / `mf.evt.array` write one array
per period instead -- the same physics, and on a 72-period 9405-cell model 6.1x
faster end to end.

The default that matters is `irch="top_active"`. MODFLOW 6 given no IRCH array
applies recharge to layer 1 unconditionally, and a column whose layer 1 is
inactive is then outside the reduced node numbering and skipped WITHOUT A WORD.
"""

from __future__ import annotations

import numpy as np
import pytest

import myflopy as mf
from myflopy.specs import ModelContext


@pytest.fixture
def ctx(canonical_run):
    return ModelContext(grid=canonical_run.vor, domain=canonical_run.gwf.disv.idomain.array)


# --- the shape ----------------------------------------------------------------


def test_rch_array_targets_the_array_class(ctx, canonical_run):
    spec = mf.rch.array(recharge={0: np.full(canonical_run.vor.ncpl, 3.5e-4)}, context=ctx)
    assert "recharge" in spec.options
    assert "stress_period_data" not in spec.options


def test_evt_array_targets_the_array_class(ctx):
    spec = mf.evt.array(rate=1e-4, depth=5.0, context=ctx)
    assert {"rate", "depth"} <= set(spec.options)
    assert "stress_period_data" not in spec.options


def test_readasarrays_is_not_passed(ctx):
    """The CLASS is the choice; passing it again is a second place to get it wrong."""

    assert "readasarrays" not in mf.rch.array(recharge=1e-4, context=ctx).options


def test_the_list_form_is_untouched(canonical_run):
    """`mf.rch(...)` and `.flopy(...)` still produce list packages."""

    spec = mf.rch.flopy(stress_period_data={0: [[(0, 3), 6.0e-4]]})
    assert "stress_period_data" in spec.options


# --- irch, which is the whole safety story ------------------------------------


def test_irch_is_derived_from_the_domain(ctx, canonical_run):
    spec = mf.rch.array(recharge=1e-4, context=ctx)
    irch = np.asarray(spec.options["irch"])
    assert irch.shape == (canonical_run.vor.ncpl,)
    assert irch.min() >= 0


def test_irch_is_zero_based(canonical_run, ctx):
    """FloPy adds one on write, so the array handed to it must be a layer INDEX.

    Passing a 1-based layer produced `2` in the file and MODFLOW 6 died with
    `Invalid layer number: 3` at the inactive columns.
    """

    # ACTIVE is `> 0` (ledger 319): a pass-through -1 is not a cell MF6 keeps. The
    # canonical domain has none, so this also read `!= 0` without failing.
    active = np.asarray(canonical_run.gwf.disv.idomain.array) > 0
    expected = np.where(active.any(axis=0), np.argmax(active, axis=0), 0)
    assert np.array_equal(np.asarray(mf.rch.array(recharge=1e-4, context=ctx).options["irch"]), expected)


def test_irch_none_is_the_explicit_opt_out(ctx):
    assert "irch" not in mf.rch.array(recharge=1e-4, irch=None, context=ctx).options


def test_irch_top_active_needs_a_domain(canonical_run):
    """Silently defaulting to layer 1 is what this refuses to do.

    A `context=` that carries no domain still cannot resolve `top_active` and still
    raises. With NO context the declaration is now LAZY instead (ledger 285) -- the
    domain comes from the discretization at build time, so there is nothing to refuse.
    """

    from myflopy.specs import ModelContext

    with pytest.raises(ValueError, match="active-cell array"):
        mf.rch.array(recharge=1e-4, context=ModelContext(grid=canonical_run.vor))

    assert sorted(mf.rch.array(recharge=1e-4).inputs) == ["irch"]


def test_an_unknown_irch_word_says_so(ctx):
    with pytest.raises(ValueError, match="top_active"):
        mf.rch.array(recharge=1e-4, irch="highest", context=ctx)


def test_an_explicit_irch_passes_through(ctx, canonical_run):
    given = np.zeros(canonical_run.vor.ncpl, dtype=int)
    assert np.array_equal(
        np.asarray(mf.rch.array(recharge=1e-4, irch=given, context=ctx).options["irch"]), given
    )


# --- what the array form cannot do, said out loud -----------------------------


def test_segmented_et_is_refused_not_silently_flattened(ctx):
    """NSEG lives in DIMENSIONS, which MODFLOW 6 never reads under READASARRAYS.

    Dropping the segments would be a change to the physics, not the format.
    """

    with pytest.raises(ValueError, match="nseg"):
        mf.evt.array(rate=1e-4, depth=5.0, nseg=2, context=ctx)


def test_an_omitted_et_surface_is_the_model_top_not_the_datum(ctx, canonical_run, tmp_path):
    """MODFLOW 6 draws the FULL rate wherever the head is at or above the ET surface.

    `surface` was put into the FloPy call only when it was not None, so omitting it
    handed the file `ModflowGwfevta`'s own default -- `CONSTANT 0.00000000`, the datum
    -- while the parameter's documentation promised the model top. That makes `depth=`
    inert and runs maximum ET everywhere under Normal termination. Measured on two real
    MF6 runs before the fix: an EVT total of -151.13 against 0.0 for the same
    declaration naming the top (ledger 283).

    Asserted on the WRITTEN FILE, because the written file is where the two forms
    differed and where a reader would have to look to notice.
    """

    grid = canonical_run.vor
    ncpl, top = int(grid.ncpl), 10.0
    # Surfaces stated on the context, so the number the test expects and the number the
    # DISV declares are the same one rather than whatever the canonical grid's own
    # table happens to carry.
    flat = ModelContext(
        grid=grid,
        domain=np.ones((1, ncpl), dtype=int),
        surfaces={"top": np.full(ncpl, top), "botm": np.zeros((1, ncpl))},
    )

    omitted = _written_evt(tmp_path, "omitted", flat,
                           mf.evt.array(rate=1e-4, depth=5.0, ievt=None, context=flat))
    named = _written_evt(tmp_path, "named", flat,
                         mf.evt.array(rate=1e-4, depth=5.0, surface=top, ievt=None, context=flat))

    assert omitted["surface"] == named["surface"]
    assert float(omitted["surface"].split()[-1]) == pytest.approx(top)


def test_an_et_surface_that_cannot_be_resolved_is_refused(canonical_run):
    """Refusing beats substituting: writing a number the caller never gave is the defect.

    The refusal is for a `context=` that cannot answer -- one with no surfaces and a
    grid with no elevation table. With NO context the declaration DEFERS instead
    (ledger 285), and the model top comes off the built discretization.
    """

    from myflopy.specs import ModelContext

    with pytest.raises(ValueError, match="surface"):
        mf.evt.array(rate=1e-4, depth=5.0, ievt=0, context=ModelContext(domain=np.ones((1, 6), dtype=int)))

    assert sorted(mf.evt.array(rate=1e-4, depth=5.0, ievt=None).inputs) == ["surface"]


def test_a_flat_land_surface_resolves_to_a_scalar(ctx, canonical_run):
    """A constant top means the same thing on every mesh, so it must not name cells.

    Resolving it per column would write an `ncpl`-long array of one repeated number,
    which is both a larger file and -- since ledger 282 -- a package that classifies
    `unverifiable` and is refused on a grid move it should survive.
    """

    from myflopy.grid_identity import CELLWISE_OPTIONS
    from myflopy.specs import GRID_PIN

    spec = mf.evt.array(rate=1e-4, depth=5.0, ievt=None, context=ctx)

    tops = np.unique(np.asarray(canonical_run.gwf.disv.top.array))
    if tops.size == 1:
        assert isinstance(spec.options["surface"], float)
        assert GRID_PIN not in (spec.metadata or {})
        assert CELLWISE_OPTIONS not in (spec.metadata or {})
    else:
        assert np.asarray(spec.options["surface"]).size == int(canonical_run.vor.ncpl)


@pytest.mark.parametrize("lazily", [False, True], ids=["eager", "lazy"])
@pytest.mark.parametrize(
    "keyword, value",
    [
        ("irch", "spread"),
        ("auxmultname", "spread"),
        ("observations", "offset"),
    ],
)
def test_a_per_feature_value_is_refused_on_every_keyword(ctx, keyword, value, lazily):
    """The refusal saw `fields` only -- `{"recharge"}` -- so every other keyword passed.

    Measured: `recharge=Spread(...)` refused, `irch=Spread(...)` accepted AND digested,
    which is precisely the exposure the refusal's own docstring claims to close
    structurally, and what ledger 274's "no stored digest can move" was resting on.
    """

    payload = mf.Spread("x") if value == "spread" else mf.CellSurfaceOffset("model_top")
    # The LAZY branch returns before `_array_bc_values`, so the refusal has to be
    # repeated there or it fires at BUILD instead of at the call (ledger 285).
    index = {} if keyword == "irch" else ({} if lazily else {"irch": 0})

    with pytest.raises(TypeError, match=keyword):
        mf.rch.array(recharge=1e-4, **{keyword: payload}, **index)


def _written_evt(tmp_path, tag, ctx, spec):
    """Build a one-package model and return its EVT file's array blocks."""

    model = mf.gwf(
        "m",
        context=ctx,
        packages=[
            mf.disv(nlay=1, top=10.0, botm=[0.0], idomain=ctx.domain,
                    **ctx.grid.get_disv_gridprops()),
            mf.npf(k=1.0), mf.ic(strt=5.0), spec, mf.oc(head_filerecord="m.hds"),
        ],
    )
    sim = mf.simulation(model, name="evtw", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    workspace = tmp_path / f"ws_{tag}"
    built = sim.build_flopy(workspace)
    built.simulation.write_simulation(silent=True)
    written = [p for p in workspace.rglob("*") if p.is_file() and "evt" in p.suffix.lower()]
    lines = written[0].read_text().splitlines()
    return {
        line.strip().lower(): lines[i + 1].strip()
        for i, line in enumerate(lines)
        if line.strip().lower() in {"surface", "rate", "depth"}
    }


def test_boundnames_are_refused(ctx):
    with pytest.raises(ValueError, match="boundnames"):
        mf.rch.array(recharge=1e-4, boundnames=True, context=ctx)
    with pytest.raises(ValueError, match="boundnames"):
        mf.evt.array(rate=1e-4, depth=5.0, boundnames=True, context=ctx)


# --- the tiers still read it --------------------------------------------------


@pytest.fixture
def array_run(canonical_run_fresh):
    """The canonical model with its list RCH swapped for the array form."""

    model = canonical_run_fresh
    context = ModelContext(grid=model.vor, domain=model.gwf.disv.idomain.array)
    model.gwf.remove_package("rch")
    mf.rch.array(
        recharge={0: np.full(model.vor.ncpl, 3.5e-4)}, context=context
    ).build(model.gwf)
    model.sim.write_simulation(silent=True)
    success, report = model.run_simulation()
    assert success, "\n".join(report[-25:])
    return model


@pytest.mark.slow
def test_an_array_package_runs(array_run):
    import flopy

    assert isinstance(array_run.package("rch"), flopy.mf6.ModflowGwfrcha)


@pytest.mark.slow
def test_the_inputs_tier_reads_an_array_package(array_run):
    """It used to raise AttributeError: an array package has no stress_period_data."""

    table = array_run.packages.rch.inputs.get()
    assert not table.empty
    assert {"per", "layer", "cell", "recharge"} <= set(table.columns)
    assert len(table) == array_run.vor.ncpl
    assert array_run.packages.rch.inputs.summary() is not None


@pytest.mark.slow
def test_the_results_tier_reads_an_array_package(array_run):
    frame = array_run.packages.rch.results.q.get()
    assert not frame.empty


@pytest.mark.slow
def test_the_budget_record_is_matched_by_package_name(array_run):
    """FloPy matches `text=` by first-hit SUBSTRING, and "RCH" is inside "RCHA".

    It is also inside "UZF-GWRCH", so a model declaring UZF first would have had
    `packages.rch.results.q` quietly return UZF's recharge term.
    """

    frame = array_run.packages.rch.results.q.get()
    column = next(c for c in frame.columns if c.startswith("q"))
    assert frame[column].sum() != 0
