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

    active = np.asarray(canonical_run.gwf.disv.idomain.array) != 0
    expected = np.where(active.any(axis=0), np.argmax(active, axis=0), 0)
    assert np.array_equal(np.asarray(mf.rch.array(recharge=1e-4, context=ctx).options["irch"]), expected)


def test_irch_none_is_the_explicit_opt_out(ctx):
    assert "irch" not in mf.rch.array(recharge=1e-4, irch=None, context=ctx).options


def test_irch_top_active_needs_a_domain():
    """Silently defaulting to layer 1 is what this refuses to do."""

    with pytest.raises(ValueError, match="active-cell array"):
        mf.rch.array(recharge=1e-4)


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
