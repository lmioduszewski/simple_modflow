"""One verb per operation, at both tiers.

Adding or replacing a package used to have four spellings:
``ModelSpec.with_package``, ``ModelSpec.with_packages``,
``SimulationSpec.add_package`` and ``SimulationSpec.replace_package`` -- the last
two differing only in the absent-slot case, and ``add_package`` recording its own
lineage entry as ``'replace_package'`` when the slot happened to exist. Removing
had two, one of which existed only on the model tier: there was NO way to drop a
simulation-level package short of rebuilding the tuple by hand.

Now ``with_packages`` / ``without_packages``, one or many, with ``on=`` choosing
the tier. The must-exist assertion that ``replace_package`` carried moved onto
``slot=``, because that assertion is load-bearing: without it a swap that misses
silently ADDS, and MF6 sums two RCH packages on one model -- measured 0.8000
against an intended 0.4000, under normal termination.
"""

from __future__ import annotations

import pytest

import myflopy as mf

_DIS = dict(nlay=1, nrow=3, ncol=3, delr=10.0, delc=10.0, top=10.0, botm=0.0)


def _model():
    return mf.gwf(
        "flow", packages=[mf.dis(**_DIS), mf.ic(strt=5.0), mf.ref("npf/base")]
    )


def _sim():
    return mf.simulation(_model(), name="base")


# --- one or many, varargs or a list ----------------------------------------- #


@pytest.mark.parametrize(
    "call",
    [
        lambda m: m.with_packages(mf.oc(), mf.sto()),
        lambda m: m.with_packages([mf.oc(), mf.sto()]),
    ],
)
def test_with_packages_takes_varargs_or_a_list(call):
    assert call(_model()).package_names == ("dis", "ic", "npf/base", "oc", "sto")


def test_with_packages_takes_exactly_one():
    assert _model().with_packages(mf.oc()).package_names[-1] == "oc"


def test_without_packages_takes_one_or_many():
    assert _model().without_packages("dis", "ic").package_names == ("npf/base",)
    assert _model().without_packages(["dis", "ic"]).package_names == ("npf/base",)


def test_a_package_name_is_not_iterated_apart():
    """A str is Iterable; flattening it would ask for a package called 'd'."""

    assert _model().without_packages("dis").package_names == ("ic", "npf/base")


# --- slot= carries the assertion replace_package used to hold ---------------- #


def test_slot_targets_a_semantic_key_the_prefix_would_miss():
    swapped = _model().with_packages(mf.ref("k/calibrated"), slot="npf/base")
    assert swapped.package_names == ("dis", "ic", "k/calibrated")


def test_slot_refuses_a_slot_that_does_not_exist():
    """The double-RCH guard: without it, a missed swap silently ADDS."""

    with pytest.raises(KeyError, match="no 'sfr' package to replace"):
        _model().with_packages(mf.ref("k/calibrated"), slot="sfr")


def test_without_slot_an_unmatched_package_is_simply_added():
    assert _model().with_packages(mf.oc()).package_names[-1] == "oc"


def test_slot_with_several_packages_is_refused():
    with pytest.raises(ValueError, match="targets ONE package"):
        _model().with_packages(mf.oc(), mf.sto(), slot="npf/base")


# --- on= chooses the tier --------------------------------------------------- #


def test_without_on_the_packages_are_the_simulations_own():
    sim = _sim().with_packages(mf.tdis(nper=3, perioddata=[(1.0, 1, 1.0)] * 3))
    assert sim.package_names == ("tdis", "flow_ims")
    assert sim.model("flow").package_names == ("dis", "ic", "npf/base")


def test_on_targets_that_models_packages():
    sim = _sim().with_packages(mf.oc(), on="flow")
    assert sim.model("flow").package_names == ("dis", "ic", "npf/base", "oc")
    assert sim.package_names == ("tdis", "flow_ims")


def test_on_carries_slot_through_to_the_model():
    sim = _sim().with_packages(mf.ref("npf/high"), on="flow", slot="npf/base")
    assert sim.model("flow").package_names == ("dis", "ic", "npf/high")


def test_slot_without_on_is_refused_at_the_simulation_tier():
    """A simulation's own packages match by full name; there is no slot hint."""

    with pytest.raises(ValueError, match="pass on=<model name>"):
        _sim().with_packages(mf.tdis(), slot="tdis")


def test_dropping_a_simulation_package_is_now_possible():
    """There was no verb for this at all -- you rebuilt the tuple by hand."""

    assert _sim().without_packages("flow_ims").package_names == ("tdis",)


def test_dropping_a_model_package_through_the_simulation():
    sim = _sim().without_packages("ic", on="flow")
    assert sim.model("flow").package_names == ("dis", "npf/base")


def test_dropping_a_simulation_package_that_is_not_there_raises():
    with pytest.raises(KeyError, match="has no package 'nope' to remove"):
        _sim().without_packages("nope")


def test_the_simulation_removal_error_points_at_the_model_tier():
    with pytest.raises(KeyError, match="on=<model name>"):
        _sim().without_packages("ic")


# --- name= makes it a variant in one step ----------------------------------- #


def test_name_returns_a_derived_simulation():
    variant = _sim().with_packages(
        mf.ref("npf/high"), on="flow", slot="npf/base", name="high_k"
    )
    assert variant.name == "high_k"
    assert variant.derived_from == "base"


def test_name_records_both_the_edit_and_the_derivation():
    variant = _sim().with_packages(mf.oc(), on="flow", name="with_oc")
    assert [step["operation"] for step in variant.lineage] == [
        "with_packages",
        "derive",
    ]


def test_the_lineage_entry_lists_every_package_it_touched():
    variant = _sim().with_packages(mf.oc(), mf.sto(), on="flow")
    assert variant.lineage[-1] == {
        "operation": "with_packages",
        "model": "flow",
        "packages": ["oc", "sto"],
    }


def test_without_a_name_the_simulation_keeps_its_own():
    assert _sim().with_packages(mf.oc(), on="flow").name == "base"


# --- the deleted verbs ------------------------------------------------------ #


@pytest.mark.parametrize(
    "verb", ["with_package", "without_package", "add_package", "replace_package"]
)
def test_the_collapsed_verbs_are_gone_from_the_simulation_tier(verb):
    assert not hasattr(mf.SimulationSpec, verb)


@pytest.mark.parametrize("verb", ["with_package", "without_package"])
def test_the_singular_verbs_are_gone_from_the_model_tier(verb):
    assert not hasattr(mf.ModelSpec, verb)
