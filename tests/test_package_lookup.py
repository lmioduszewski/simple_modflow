"""One name, one answer: ``package`` / ``with_package`` / ``without_package``.

Three lookups had drifted apart. On a model carrying ``mf.ref("wel/north")`` and
``mf.ref("wel/south")``, asking each of them about ``"wel"`` gave three different
answers:

    package("wel")           -> wel/north
    with_package(x, "wel")   -> replaced wel/north
    without_package("wel")   -> deleted BOTH

Nobody chose any of that. ``_package_entry_index`` grew its label-first pass after
a swap deleted the wrong well package and MF6 reported normal termination on the
wrong model (ledger, `specs.py`); the fix reached that one lookup and neither of
the others. They now share a resolver, and an ambiguous name RAISES rather than
being settled by declaration order.

``without_package`` also used to be a silent no-op on a name it did not find,
while still returning a new object -- so a typo read exactly like success.
"""

from __future__ import annotations

import pytest

import myflopy as mf


def _ambiguous():
    """A model whose `wel` slot is claimed by two differently-keyed refs."""

    return mf.gwf(
        "flow",
        packages=[mf.npf(k=1.0), mf.ic(strt=1.0), mf.ref("wel/north"), mf.ref("wel/south")],
    )


def _unambiguous():
    return mf.gwf("flow", packages=[mf.npf(k=1.0), mf.ic(strt=1.0), mf.ref("wel/north")])


# --- an ambiguous name is refused, not resolved by list order --------------- #


def test_package_refuses_an_ambiguous_slot():
    with pytest.raises(KeyError, match="wel/north, wel/south"):
        _ambiguous().package("wel")


def test_with_package_refuses_an_ambiguous_slot():
    with pytest.raises(KeyError, match="wel/north, wel/south"):
        _ambiguous().with_package(mf.ref("wel/east"), slot="wel")


def test_without_package_refuses_an_ambiguous_slot_instead_of_deleting_both():
    model = _ambiguous()
    with pytest.raises(KeyError, match="wel/north, wel/south"):
        model.without_package("wel")
    assert model.package_names == ("npf", "ic", "wel/north", "wel/south")


def test_the_refusal_names_the_model_and_both_candidates():
    with pytest.raises(KeyError) as excinfo:
        _ambiguous().package("wel")
    message = str(excinfo.value)
    assert "Model 'flow'" in message
    assert "wel/north" in message and "wel/south" in message


# --- and an unambiguous one behaves identically across all three ------------ #


def test_an_exact_label_resolves_even_when_a_slot_hint_would_too():
    assert _ambiguous().package("wel/south").key == "wel/south"


def test_a_slot_hint_still_works_when_only_one_entry_claims_it():
    assert _unambiguous().package("wel").key == "wel/north"


def test_the_three_verbs_agree_on_the_entry_a_name_resolves_to():
    model = _unambiguous()
    replaced = model.with_package(mf.ref("wel/east"), slot="wel")
    removed = model.without_package("wel")

    assert model.package("wel").key == "wel/north"
    assert replaced.package_names == ("npf", "ic", "wel/east")
    assert removed.package_names == ("npf", "ic")


def test_without_package_removes_exactly_one_named_entry():
    model = _ambiguous().without_package("wel/north")
    assert model.package_names == ("npf", "ic", "wel/south")


# --- the silent no-op ------------------------------------------------------- #


def test_without_package_raises_on_a_name_it_cannot_find():
    """It used to return a new object having changed nothing."""

    with pytest.raises(KeyError, match="no package 'sfr' to remove"):
        _unambiguous().without_package("sfr")


def test_the_not_found_message_lists_what_the_model_carries():
    with pytest.raises(KeyError) as excinfo:
        _unambiguous().without_package("sfr")
    assert "npf, ic, wel/north" in str(excinfo.value)


def test_package_not_found_also_lists_what_the_model_carries():
    with pytest.raises(KeyError) as excinfo:
        _unambiguous().package("sfr")
    assert "npf, ic, wel/north" in str(excinfo.value)


def test_replace_package_reports_ambiguity_against_the_named_model():
    sim = mf.simulation(_ambiguous(), name="s")
    with pytest.raises(KeyError, match="wel/north, wel/south"):
        sim.replace_package("flow", mf.ref("wel/east"), name="wel")
