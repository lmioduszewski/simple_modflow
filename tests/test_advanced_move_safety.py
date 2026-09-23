"""The advanced packages must not move silently (ledgers 302, 303).

Two holes, both silent-wrong-model rather than wrong-error:

* a GRIDLESS ``mf.uzf`` carried no pin and no cell marker, so ``package_move`` answered
  "unchanged -- names no cells" and a move BUILT CLEAN with the old mesh's UZF objects
  sitting on a prefix of the new mesh's columns;
* ``mf.mvr`` names its endpoints as bare integers and nothing compared them against the
  providers, so a mover beside a moved SFR routed water to a reach that no longer means
  what it meant.

Both were invisible in the same way: the verdict said the move was fine.
"""

from __future__ import annotations

import numpy as np
import pytest

import myflopy as mf
from myflopy.grid_identity import CELLWISE_OPTIONS, GRID_PIN, package_move
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, rectangular_voronoi
from myflopy.specs import ModelContext


def _grid(n, nlay=3):
    return rectangular_voronoi(
        CanonicalModelConfig(nrow=n, ncol=n, nlay=nlay, nper=1)
    )


@pytest.fixture(scope="module")
def meshes():
    return _grid(5), _grid(9)


def _uzf(context):
    return mf.uzf(
        context=context, nper=1, finf=1e-3, pet=1e-3,
        vks=1.0, thtr=0.05, thts=0.35, thti=0.08, eps=4.0,
    )


def _stack():
    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    for name in ("a", "b", "c"):
        stack.add(name, thickness=20.0)
    return stack


# ------------------------------------------------------------------ C1: gridless uzf


def test_a_gridless_uzf_is_marked_as_addressing_cells(meshes):
    """It has no mesh to pin to, so it gets the weaker marker rather than nothing."""

    small, _ = meshes
    domain = np.ones((3, small.ncpl), dtype=int)
    spec = _uzf(ModelContext(domain=domain))

    assert GRID_PIN not in spec.metadata, "there was no grid to pin to"
    assert spec.metadata[CELLWISE_OPTIONS] == ["packagedata"]
    how, why = package_move(spec)
    assert how == "unverifiable", why


def test_a_grid_bearing_uzf_still_pins(meshes):
    """The eager, documented declaration must not move -- it was already safe."""

    small, _ = meshes
    domain = np.ones((3, small.ncpl), dtype=int)
    spec = _uzf(ModelContext(grid=small, domain=domain))

    assert GRID_PIN in spec.metadata
    assert CELLWISE_OPTIONS not in spec.metadata
    assert package_move(spec)[0] == "stale"


def test_the_marker_does_not_change_what_was_built(meshes):
    """A stamp is provenance; it must not touch the package's own records."""

    small, _ = meshes
    domain = np.ones((3, small.ncpl), dtype=int)
    gridless = _uzf(ModelContext(domain=domain))
    pinned = _uzf(ModelContext(grid=small, domain=domain))
    assert gridless.options["packagedata"] == pinned.options["packagedata"]
    assert gridless.options["perioddata"] == pinned.options["perioddata"]


def test_the_marker_moves_no_digest(meshes):
    """`package_digest` excludes metadata by contract -- assert it rather than trust it.

    Without this, closing the hole would report every stored run as edited on the day it
    landed, which is ledger 293's churn.
    """

    from dataclasses import replace

    from myflopy.specs import package_digest

    small, _ = meshes
    spec = _uzf(ModelContext(domain=np.ones((3, small.ncpl), dtype=int)))
    assert package_digest(spec) == package_digest(replace(spec, metadata={}))


def test_moving_a_gridless_uzf_is_now_REFUSED_at_build(meshes, tmp_path):
    """The defect itself: this used to build clean on a prefix of the new mesh.

    Measured before the fix: 75 UZF objects on columns 0..24 of an 81-cell mesh, zero
    warnings -- UZF over 31% of the domain.
    """

    small, large = meshes
    domain = np.ones((3, small.ncpl), dtype=int)
    sim = mf.simulation(
        mf.gwf("m", context=ModelContext(grid=small, domain=domain),
               packages=[mf.disv(_stack()), mf.npf(k=1.0), mf.ic(strt=50.0),
                         _uzf(ModelContext(domain=domain))]),
        name="s", packages=[mf.tdis(nper=1)],
    )
    moved = sim.with_grid(large)
    assert next(r["how"] for r in moved.move_report() if r["package"] == "uzf") == "unverifiable"
    with pytest.raises(ValueError):
        moved.build(tmp_path / "moved")


def test_a_gridless_uzf_still_builds_on_its_OWN_mesh(meshes, tmp_path):
    """The fix must not narrow the legitimate case.

    A gridless declaration is correct for a model that never moves, and
    `_lazy_resolution_context`'s own docstring records that refusing one would make the
    lazy form narrower than the eager one (ledger 284). Two tests in
    `test_uzf_builder.py` declare exactly this shape.
    """

    small, _ = meshes
    domain = np.ones((3, small.ncpl), dtype=int)
    sim = mf.simulation(
        mf.gwf("m", context=ModelContext(grid=small, domain=domain),
               packages=[mf.disv(_stack()), mf.npf(k=1.0), mf.ic(strt=50.0),
                         _uzf(ModelContext(domain=domain))]),
        name="s", packages=[mf.tdis(nper=1)],
    )
    built = sim.build(tmp_path / "own")
    uzf = built.simulation.get_model("m").get_package("uzf")
    assert len(uzf.packagedata.get_data()) == len(sim.model("m").packages[-1].options["packagedata"])


def test_the_limit_is_the_MOVE_not_every_wrong_grid(meshes, tmp_path):
    """What the marker does NOT close, pinned so the commit cannot overclaim.

    The stamp is read by `package_move`, which only `with_grid` and `move_report`
    consult. A model ASSEMBLED BY HAND against a different grid never goes through
    either, and still builds -- measured. Closing that is the build-time consistency
    guard's job (ledger 249), and it is a range check, not a proof.
    """

    small, large = meshes
    domain = np.ones((3, small.ncpl), dtype=int)
    hand_made = mf.simulation(
        mf.gwf("m", context=ModelContext(grid=large, domain=np.ones((3, large.ncpl), int)),
               packages=[mf.disv(_stack()), mf.npf(k=1.0), mf.ic(strt=50.0),
                         _uzf(ModelContext(domain=domain))]),
        name="s", packages=[mf.tdis(nper=1)],
    )
    built = hand_made.build(tmp_path / "byhand")
    rows = built.simulation.get_model("m").get_package("uzf").packagedata.get_data()
    columns = {int(r[1][1]) for r in rows}
    assert max(columns) < large.ncpl - 1, (
        "this test exists because the small mesh's cells are a PREFIX of the large "
        "one's, and nothing on this path objects"
    )
