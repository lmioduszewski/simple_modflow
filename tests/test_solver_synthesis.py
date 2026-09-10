"""``mf.simulation`` has to pick a solver MF6 will accept, and it did not.

Two failures, both on one-liners the docstrings teach:

``mf.simulation(flow, transport)`` -- labelled "one IMS each (coupled-ready)" in
``simulation()``'s own docstring -- wrote ``COMPLEXITY simple`` onto the GWT
solver, which selects conjugate gradient. MF6 aborts: *"MODEL TRANS PRODUCES AN
ASYMMETRIC COEFFICIENT MATRIX, BUT THE CONJUGATE GRADIENT METHOD WAS SELECTED.
USE BICGSTAB INSTEAD."* Nothing in the suite built a spec-declared GWF-GWT and
ran it -- ``test_gwt_gwe_results.py`` hand-sets complexity AND a BICGSTAB
override -- which is exactly how it survived. These tests reach MF6.

``mf.simulation(flow)`` failed the same way whenever ``flow`` carried
``newtonoptions``, because Newton makes a GWF matrix asymmetric too. So the rule
is keyed on the MATRIX, not the model kind: a rule written against ``model_type``
would have fixed the coupled case and left every unconfined model broken.

Solver ORDER is the second trap: MF6 requires a GWF solution to be listed in
``mfsim.nam`` before the solution of any model it feeds, and nothing enforced it.
"""

from __future__ import annotations

import pytest

import myflopy as mf
from myflopy.package_api import _asymmetric_matrix, _ordered_solvers
from myflopy.specs import ExchangeSpec

_DIS = dict(nlay=1, nrow=3, ncol=3, delr=10.0, delc=10.0, top=10.0, botm=0.0)


def _flow(**model_options):
    """A GWF model whose CHD carries an aux concentration for SSM to read."""

    return mf.gwf(
        "flow",
        packages=[
            mf.dis(**_DIS),
            mf.ic(strt=5.0),
            mf.npf(k=1.0, icelltype=1),
            mf.chd(
                stress_period_data={0: [((0, 0, 0), 5.0, 1.0), ((0, 0, 2), 4.0, 0.0)]},
                auxiliary=["concentration"],
            ),
            mf.oc(
                head_filerecord="flow.hds",
                budget_filerecord="flow.cbc",
                saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
            ),
        ],
        **model_options,
    )


def _transport():
    return mf.gwt(
        "trans",
        packages=[
            mf.dis(**_DIS),
            mf.ic(strt=0.0),
            mf.adv(scheme="UPSTREAM"),
            mf.mst(porosity=0.3),
            mf.ssm(sources=[["chd", "AUX", "concentration"]]),
            mf.oc(
                concentration_filerecord="trans.ucn",
                saverecord=[("CONCENTRATION", "ALL")],
            ),
        ],
    )


def _gwf_gwt_exchange():
    from myflopy.builders import build_gwf_gwt_exchange

    return ExchangeSpec("gwfgwt", build_gwf_gwt_exchange, models=("flow", "trans"))


def _solvers(spec):
    return [p for p in spec.packages if p.name != "tdis"]


def _run(spec, tmp_path):
    built = spec.build_flopy(tmp_path)
    built.simulation.write_simulation(silent=True)
    success, _ = built.simulation.run_simulation(silent=True)
    return success


# --- the predicate ---------------------------------------------------------- #


@pytest.mark.parametrize(
    "model, expected",
    [
        (mf.gwt("t"), True),
        (mf.gwe("e"), True),
        (mf.gwf("f"), False),
        (mf.gwf("f", newtonoptions="NEWTON"), True),
        (mf.gwf("f", newtonoptions="NEWTON UNDER_RELAXATION"), True),
        (mf.prt("p"), False),
    ],
)
def test_asymmetric_matrix_is_keyed_on_the_matrix_not_the_kind(model, expected):
    assert _asymmetric_matrix(model) is expected


def test_the_synthesized_solver_only_asks_for_bicgstab_where_it_is_needed():
    spec = mf.simulation(_flow(), _transport(), name="c")
    acceleration = {p.name: p.options.get("linear_acceleration") for p in _solvers(spec)}
    assert acceleration == {"flow_ims": None, "trans_ims": "BICGSTAB"}


# --- solver ordering -------------------------------------------------------- #


def test_solvers_are_ordered_by_model_declaration_not_by_the_list_given():
    """MF6: "The IMS specified for GWF must be listed ... before the IMS for GWT"."""

    spec = mf.simulation(
        _flow(),
        _transport(),
        name="c",
        packages=[
            mf.ims(models=["trans"], name="trans_ims", linear_acceleration="BICGSTAB"),
            mf.ims(models=["flow"], name="flow_ims"),
        ],
    )
    assert [p.name for p in _solvers(spec)] == ["flow_ims", "trans_ims"]


def test_a_solver_naming_no_known_model_keeps_its_place_at_the_end():
    flow, transport = _flow(), _transport()
    stray = mf.ims(models=["nobody"], name="stray")
    ordered = _ordered_solvers(
        [stray, mf.ims(models=["trans"], name="t"), mf.ims(models=["flow"], name="f")],
        (flow, transport),
    )
    assert [p.name for p in ordered] == ["f", "t", "stray"]


# --- what MF6 says, which is the only thing that settles it ----------------- #


@pytest.mark.slow
def test_the_coupled_one_liner_actually_runs(tmp_path):
    """The exact line `simulation()`'s docstring calls "coupled-ready"."""

    spec = mf.simulation(
        _flow(), _transport(), name="c", exchanges=(_gwf_gwt_exchange(),)
    )
    assert _run(spec, tmp_path) is True


@pytest.mark.slow
def test_an_unconfined_newton_model_runs_from_the_one_liner(tmp_path):
    spec = mf.simulation(_flow(newtonoptions="NEWTON"), name="n")
    assert _run(spec, tmp_path) is True


@pytest.mark.slow
def test_a_backwards_solver_list_still_produces_a_runnable_name_file(tmp_path):
    spec = mf.simulation(
        _flow(),
        _transport(),
        name="c",
        exchanges=(_gwf_gwt_exchange(),),
        packages=[
            mf.ims(models=["trans"], name="trans_ims", linear_acceleration="BICGSTAB"),
            mf.ims(models=["flow"], name="flow_ims"),
        ],
    )
    assert _run(spec, tmp_path) is True


# --- assembling `packages=` ------------------------------------------------- #


def test_two_timing_packages_are_refused():
    """Both used to survive: FloPy warned "Replacing", and the last one won."""

    with pytest.raises(ValueError, match="2 timing packages"):
        mf.simulation(
            _flow(), name="s", packages=[mf.tdis(), mf.tdis(name="tdis2")]
        )


def test_two_solvers_with_no_models_are_refused():
    """Which one claims the leftovers has no defined answer."""

    with pytest.raises(ValueError, match="2 solvers with no models="):
        mf.simulation(_flow(), name="s", packages=[mf.ims(), mf.ims(name="ims2")])


def test_a_solver_with_no_models_claims_what_the_others_left():
    spec = mf.simulation(
        _flow(),
        _transport(),
        name="s",
        packages=[mf.ims(models=["trans"], name="trans_ims"), mf.ims(name="rest")],
    )
    claimed = {p.name: p.options.get("models") for p in _solvers(spec)}
    assert claimed == {"rest": ("flow",), "trans_ims": ("trans",)}


def test_prt_gets_its_own_ems_even_when_an_open_ims_is_present():
    """MF6 rejects an explicit model under IMS6, so PRT is never folded in."""

    particles = mf.prt("part", packages=[mf.dis(**_DIS), mf.mip(porosity=0.3)])
    spec = mf.simulation(
        _flow(), particles, name="s", packages=[mf.ims(complexity="COMPLEX")]
    )
    claimed = {p.name: p.options.get("models") for p in _solvers(spec)}
    assert claimed == {"ims": ("flow",), "part_ems": ("part",)}


def test_an_open_solver_with_nothing_left_to_claim_is_refused():
    with pytest.raises(ValueError, match="none left for it to claim"):
        mf.simulation(
            _flow(),
            name="s",
            packages=[mf.ims(models=["flow"], name="explicit"), mf.ims(name="open")],
        )


def test_a_supplied_timing_package_replaces_the_default():
    spec = mf.simulation(
        _flow(), name="s", packages=[mf.tdis(nper=3, perioddata=[(1.0, 1, 1.0)] * 3)]
    )
    timing = [p for p in spec.packages if p.name == "tdis"]
    assert len(timing) == 1
    assert timing[0].options["nper"] == 3


def test_a_ref_passes_through_and_does_not_suppress_the_timing_default():
    """A ref carries no builder, so it cannot be identified as timing or solver."""

    spec = mf.simulation(_flow(), name="s", packages=[mf.ref("tdis/monthly")])
    assert "tdis" in spec.package_names          # the default is still added
    assert "tdis/monthly" in spec.package_names  # and the ref is kept
