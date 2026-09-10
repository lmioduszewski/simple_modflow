"""A simulation carries timing and solvers. A `chd` handed to it is a mistake.

It used to be accepted at declaration and die much later inside FloPy as
``Invalid package type "chd". Unable to register package.`` -- a message that
names the type but not the mistake, and arrives a long way from the line that
made it.

The refusal reads the BUILDER's first constructor parameter, never the package
name. Measured across every ``flopy.mf6`` package class, that partitions them
cleanly (``model`` 78, ``simulation`` 15, ``parent_model_or_package`` 8,
``parent_package`` 7) while a name-prefix heuristic is wrong on 17: all six
exchange classes carry a model-kind prefix and are simulation-scoped, and eleven
``utl*`` classes read as simulation and are not.

The other half of the trap is ``functools.partial``. The seven list BCs all build
through ``partial(_build_named, ModflowGwfchd)``, so a check that does not unwrap
measures them as unknown -- and would wave through exactly the packages it exists
to catch, while still catching ``npf``.
"""

from __future__ import annotations

import flopy
import pytest

import myflopy as mf
from myflopy.specs import _package_scope

_DIS = dict(nlay=1, nrow=1, ncol=1, delr=1.0, delc=1.0, top=10.0, botm=0.0)


def _flow():
    return mf.gwf(
        "flow", packages=[mf.dis(**_DIS), mf.ic(strt=5.0), mf.npf(k=1.0)]
    )


# --- the signal ------------------------------------------------------------- #


@pytest.mark.parametrize(
    "label, spec, expected",
    [
        ("chd", mf.chd(stress_period_data={0: [((0, 0, 0), 1.0)]}), "model"),
        ("wel", mf.wel(stress_period_data={0: [((0, 0, 0), -1.0)]}), "model"),
        ("npf", mf.npf(k=1.0), "model"),
        ("oc", mf.oc(), "model"),
        ("tdis", mf.tdis(), "simulation"),
        ("ims", mf.ims(models=["flow"]), "simulation"),
        ("ems", mf.ems(models=["part"]), "simulation"),
    ],
)
def test_scope_is_read_from_the_builder(label, spec, expected):
    assert _package_scope(spec.builder) == expected


def test_a_partial_wrapped_list_bc_is_still_seen_as_a_model_package():
    """`mf.chd` builds through `partial(_build_named, ModflowGwfchd)`.

    Without unwrapping it measures as unknown -- so the check would pass every
    list BC and catch only the packages nobody hands to a simulation by mistake.
    """

    import functools

    builder = mf.chd(stress_period_data={0: [((0, 0, 0), 1.0)]}).builder
    assert isinstance(builder, functools.partial)
    assert _package_scope(builder) == "model"


@pytest.mark.parametrize(
    "cls, expected",
    [
        (flopy.mf6.ModflowGwfgwt, "simulation"),  # model-kind prefix, sim-scoped
        (flopy.mf6.ModflowGwfgwf, "simulation"),
        # "utl" prefix, and neither of them is simulation-scoped -- laktab hangs
        # off a MODEL, ats off a parent PACKAGE. Two different answers, one
        # prefix.
        (flopy.mf6.ModflowUtllaktab, "model"),
        (flopy.mf6.ModflowUtlats, None),
    ],
)
def test_a_name_prefix_heuristic_would_get_these_wrong(cls, expected):
    assert _package_scope(cls) == expected


def test_a_package_that_is_genuinely_either_is_not_refused():
    """MVR/GNC take `parent_model_or_package`; abstaining beats guessing."""

    assert _package_scope(flopy.mf6.ModflowGwfmvr) is None


# --- the refusal ------------------------------------------------------------ #


def test_a_model_package_at_the_simulation_tier_is_refused():
    with pytest.raises(ValueError, match="was given model package"):
        mf.simulation(
            _flow(),
            name="s",
            packages=[mf.chd(stress_period_data={0: [((0, 0, 0), 1.0)]})],
        )


def test_the_refusal_names_the_package_and_the_fix():
    with pytest.raises(ValueError) as excinfo:
        mf.simulation(_flow(), name="s", packages=[mf.npf(k=1.0)])
    message = str(excinfo.value)
    assert "npf" in message
    assert "mf.gwf(packages=[...])" in message


def test_with_package_on_a_simulation_refuses_a_model_package_too():
    """`with_package` rebuilds through `replace`, so `__post_init__` runs."""

    sim = mf.simulation(_flow(), name="s")
    with pytest.raises(ValueError, match="was given model package"):
        sim.with_packages(mf.oc(printrecord=[("BUDGET", "ALL")]))


def test_timing_and_solvers_are_accepted():
    sim = mf.simulation(_flow(), name="s")
    assert sim.package_names == ("tdis", "flow_ims")


def test_a_package_reference_is_not_refused():
    """A ref carries no builder; it cannot be classified before it resolves."""

    sim = mf.simulation(_flow(), name="s", packages=[mf.ref("ims/shared")])
    assert "ims/shared" in sim.package_names
