"""Coupling two models needed a name no IDE would offer.

The four ``build_*_exchange`` functions resolve on ``mf`` but appear in neither
``mf.__all__`` nor ``dir(mf)``, so autocomplete never offered the one name a
coupled model cannot be built without -- and the spelling shown in four
docstrings, ``mf.build_gwf_gwt_exchange("flow", "transport")``, raises
``ValueError``: the real builder takes ``(simulation, models, **options)`` and
has to be wrapped in an ``ExchangeSpec``.

Nothing wires an exchange for you. ``mf.simulation(flow, transport)`` builds both
models, moves nothing between them, and MF6 reports normal termination.
"""

from __future__ import annotations

import pytest

import myflopy as mf
from myflopy.builders import (
    build_gwf_gwe_exchange,
    build_gwf_gwf_exchange,
    build_gwf_gwt_exchange,
    build_gwf_prt_exchange,
)

_DIS = dict(nlay=1, nrow=3, ncol=3, delr=10.0, delc=10.0, top=10.0, botm=0.0)


def _flow(name="flow"):
    return mf.gwf(name, packages=[mf.dis(**_DIS), mf.ic(strt=5.0), mf.npf(k=1.0)])


def _transport():
    return mf.gwt("trans", packages=[mf.dis(**_DIS), mf.ic(strt=0.0)])


# --- discoverability, which was the whole complaint ------------------------- #


def test_exchange_is_offered_by_autocomplete():
    assert "exchange" in mf.__all__
    assert "exchange" in dir(mf)


def test_every_coupling_has_an_explicit_spelling():
    assert sorted(n for n in dir(mf.exchange) if not n.startswith("_")) == [
        "gwfgwe",
        "gwfgwf",
        "gwfgwt",
        "gwfprt",
    ]


# --- inference -------------------------------------------------------------- #


@pytest.mark.parametrize(
    "target, builder",
    [
        (mf.gwt("trans"), build_gwf_gwt_exchange),
        (mf.gwe("heat"), build_gwf_gwe_exchange),
        (mf.prt("part"), build_gwf_prt_exchange),
    ],
)
def test_the_pair_of_kinds_picks_the_mf6_class(target, builder):
    assert mf.exchange(_flow(), target).builder is builder


def test_the_gwf_model_is_put_first_whichever_way_round_you_ask():
    """MF6 wants (gwf, gwt); nobody gets that wrong on purpose."""

    assert mf.exchange(_transport(), _flow()).models == ("flow", "trans")


def test_the_exchange_is_named_after_the_pair():
    assert mf.exchange(_flow(), _transport()).name == "flow-trans"


def test_an_unsupported_pairing_is_refused_with_the_list():
    with pytest.raises(ValueError, match="Supported: gwf-gwe, gwf-gwf"):
        mf.exchange(_transport(), _transport())


def test_inference_refuses_bare_names_because_kinds_are_what_it_reads():
    with pytest.raises(TypeError, match="needs the model SPECS"):
        mf.exchange("flow", "trans")


# --- explicit spellings ----------------------------------------------------- #


def test_an_explicit_spelling_accepts_plain_names():
    assert mf.exchange.gwfgwt("flow", "trans").models == ("flow", "trans")


def test_gwfgwf_refuses_the_file_mf6_rejects():
    """A GWF-GWF with no options writes only an empty options block."""

    with pytest.raises(ValueError, match="nexg and exchangedata"):
        mf.exchange.gwfgwf(_flow("a"), _flow("b"))


def test_gwfgwf_accepts_the_dimensions_it_needs():
    spec = mf.exchange.gwfgwf(
        _flow("a"), _flow("b"), nexg=1, exchangedata=[[(0, 0), (0, 1), 1, 5.0, 5.0, 10.0]]
    )
    assert spec.builder is build_gwf_gwf_exchange
    assert spec.options["nexg"] == 1


def test_a_custom_name_survives():
    assert mf.exchange(_flow(), _transport(), name="coupling").name == "coupling"
