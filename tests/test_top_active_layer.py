"""One rule for "the top active layer of a column": ``idomain > 0`` (ledger 319).

Five places asked the question and two predicates answered it. The areal RCH/EVT
builders, UZF and SFR read ``> 0``; the array-form helpers (``irch``/``ievt=
"top_active"``) and the MODFLOW-USG importer read ``!= 0``. For the array form that
was a defect: an MF6 ``idomain`` of ``-1`` is a vertical PASS-THROUGH cell, removed
from the model, so a column whose top cell was one had its recharge written to a cell
MF6 then skips without a word -- 40 of 92 columns, 42.2% of the model's recharge,
under Normal termination (``docs/reproducers/top_active_passthrough.py``). For the
USG importer ``!= 0`` was neutral: its ``idomain`` is already 0/1, because a
NEGATIVE USG IBOUND is a constant head, an ACTIVE cell.

Every caller now asks :func:`~myflopy.modflow.mf6.grid.geometry.top_active_layer`.
"""

from __future__ import annotations

import contextlib
import io
import pathlib
import types
import warnings

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import LineString, box

import myflopy as mf
from myflopy.modflow.mf6.grid.geometry import top_active_layer
from myflopy.specs import ModelContext


def _passthrough_domain(grid):
    """Three layers on ``grid``: -1 on top in the west, a 0 under part of it, a hole."""

    x = np.asarray(grid.centroids_x, dtype=float)
    y = np.asarray(grid.centroids_y, dtype=float)
    domain = np.ones((3, int(grid.ncpl)), dtype=int)
    domain[0, x < 250] = -1
    domain[1, (x < 250) & (y < 300)] = 0
    domain[:, (x > 450) & (y > 450)] = 0
    return domain


# --- the rule itself ----------------------------------------------------------


def test_a_passthrough_or_inactive_cell_is_not_the_top():
    domain = np.array([
        [1, -1, 0, -1, 0],
        [1, 1, 1, -1, 0],
        [1, 1, 1, 1, 0],
    ])

    assert top_active_layer(domain).tolist() == [0, 1, 1, 2, -1]


def test_a_one_dimensional_domain_is_one_layer():
    """The reading `domain_active` gives it: the same activity in every layer."""

    assert top_active_layer([1, 0, -1, 2]).tolist() == [0, -1, -1, 0]


def test_what_asarray_accepts_is_accepted():
    """A nested list and a float array read the same as an int array."""

    as_ints = top_active_layer(np.array([[0, 1], [1, 1]]))
    assert top_active_layer([[0, 1], [1, 1]]).tolist() == as_ints.tolist() == [1, 0]
    assert top_active_layer(np.array([[0.0, 1.0], [1.0, 1.0]])).tolist() == [1, 0]
    assert as_ints.dtype.kind == "i"


def test_a_structured_domain_is_refused_by_shape():
    with pytest.raises(ValueError, match=r"\(nlay, nrow \* ncol\)"):
        top_active_layer(np.ones((2, 3, 4), dtype=int))


# --- the array form, eager and lazy -------------------------------------------


@pytest.mark.parametrize("package", ["rch", "evt"])
def test_the_array_form_skips_a_passthrough_cell(triangulated_grid, package):
    """The defect itself: the index named layer 0 -- the removed cell -- in those columns."""

    domain = _passthrough_domain(triangulated_grid)
    context = ModelContext(grid=triangulated_grid, domain=domain)
    spec = (
        mf.rch.array(recharge=1e-3, context=context)
        if package == "rch"
        else mf.evt.array(rate=1e-3, depth=2.0, surface=100.0, context=context)
    )
    index = np.asarray(spec.options["irch" if package == "rch" else "ievt"])

    passthrough = domain[0] == -1
    assert passthrough.sum() >= 3, "a domain with no pass-through proves nothing here"
    empty = ~(domain > 0).any(axis=0)
    assert np.all(index[passthrough & ~empty] >= 1)
    assert np.array_equal(index, np.where(empty, 0, top_active_layer(domain)))
    assert np.all(domain[index[~empty], np.flatnonzero(~empty)] > 0), (
        "every column with an active cell must be indexed AT one"
    )


def _passthrough_stack():
    """A stack whose top unit pinches out in the west under `pinch="passthrough"`."""

    return (
        mf.LayerStack(top=mf.Surface.flat(100.0))
        .add("fill", bottom=mf.Surface.where(box(0, 0, 250, 600), 100.0, 95.0),
             pinch="passthrough")
        .add("sand", bottom=80.0)
        .add("base", bottom=50.0)
    )


def test_the_lazy_array_form_resolves_past_a_passthrough_cell(tmp_path, triangulated_grid):
    """The lazy recipe holds the STRING, so its digest cannot move -- its BUILD did."""

    model = mf.gwf("m", context=ModelContext(grid=triangulated_grid), packages=[
        mf.disv(_passthrough_stack(), attach=False),
        mf.npf(k=10.0), mf.ic(strt=95.0), mf.oc(head_filerecord="m.hds"),
        mf.rch.array(recharge=1e-3),
    ])
    simulation = mf.simulation(
        model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])]
    )

    built = simulation.build_flopy(tmp_path / "w").simulation.get_model("m")

    idomain = np.asarray(built.disv.idomain.array).reshape(3, -1)
    irch = built.get_package("rch").irch.get_data()
    irch = np.asarray(irch[0] if isinstance(irch, dict) else irch).reshape(-1)
    assert (idomain[0] == -1).sum() >= 3
    assert np.all(irch[idomain[0] == -1] == 1)
    assert np.array_equal(irch, top_active_layer(idomain))


# --- every caller agrees ------------------------------------------------------


def test_every_top_active_caller_gives_the_same_layer(triangulated_grid):
    """Areal RCH/EVT, UZF, SFR's default reach layer and the array form: one answer.

    Each is asked of the SAME pass-through domain. They were already `> 0` except the
    array form; this pins that they stay one rule rather than five copies.
    """

    from myflopy.modflow.mf6.evapotranspiration import EVTBuilder
    from myflopy.modflow.mf6.recharge import RCHBuilder
    from myflopy.modflow.mf6.sfr import SFRBuilder
    from myflopy.modflow.mf6.uzf import UZFBuilder

    domain = _passthrough_domain(triangulated_grid)
    context = ModelContext(grid=triangulated_grid, domain=domain)
    top = top_active_layer(domain)
    expected = {(int(layer), cell) for cell, layer in enumerate(top) if layer >= 0}
    assert any(domain[0, cell] == -1 for _, cell in expected)

    rch = set(RCHBuilder(context=context, nper=1, recharge=1e-3).selected_cells)
    evt = set(EVTBuilder(context=context, nper=1, surface=100.0, rate=1e-3,
                         depth=2.0).selected_cells)
    uzf = set(UZFBuilder(context=context, nper=1, vks=1.0, thtr=0.1, thts=0.3, thti=0.2,
                         cells="surface_only").uzf_cells)
    sfr = SFRBuilder(
        context=context, nper=1,
        streams=gpd.GeoDataFrame({"id": ["a"]}, geometry=[LineString([(10, 300), (590, 300)])],
                                 crs=2927),
    )
    array = np.asarray(mf.rch.array(recharge=1e-3, context=context).options["irch"])

    assert rch == evt == uzf == expected
    for cell, layer in enumerate(top):
        if layer >= 0:
            assert sfr._layer_for_cell("a", cell) == layer
            assert array[cell] == layer
        else:
            with pytest.raises(ValueError, match="inactive grid column"):
                sfr._layer_for_cell("a", cell)


def test_a_usg_constant_head_is_the_top_active_cell():
    """USG's negative IBOUND is a constant head -- ACTIVE, unlike MF6's pass-through -1.

    `UsgModel.idomain` maps it to 1 before `uppermost_active` reads it, which is what
    makes the shared `> 0` rule the right one there. A change that let -1 through to
    `idomain` would move the ETS/RCH target of every constant-head column down a layer.
    """

    from myflopy.modflow.usg.model import UsgModel

    stub = types.SimpleNamespace(
        bas=types.SimpleNamespace(ibound=np.array([-1, 0, 1, 0, 1, 1, 1, 0])),
        nlay=2, ncpl=4,
    )
    stub._layered = types.MethodType(UsgModel._layered, stub)
    stub.idomain = UsgModel.idomain.fget(stub)

    assert UsgModel.uppermost_active.fget(stub).tolist() == [0, 1, 0, 0]


# --- what MF6 does with it ----------------------------------------------------


@pytest.mark.slow
@pytest.mark.parametrize("fixed_cell", [False, True])
def test_a_passthrough_column_gets_its_recharge(tmp_path, triangulated_grid_finer, fixed_cell):
    """Run MF6: the recharge the pass-through columns were declared is the recharge applied.

    Asserted on the budget, and on MF6's percent discrepancy -- never on Normal
    termination alone, which the defect also reached.

    With `fixed_cell=True` the pass-through columns still get theirs (their IRCH cell is
    the wet sand below), and what is lost is exactly the DRY fill cells' share: the
    steady head sits below the fill's 95 ft bottom, so every active fill cell is dry.
    The `rch.array` docstring states that; it first stated the OLD `irch`'s number.
    """

    import flopy

    grid = triangulated_grid_finer
    rate = 1e-3
    east = [int(c) for c in np.flatnonzero(grid.gdf_vorPolys.geometry.bounds.maxx >= 600 - 1e-6)]
    model = mf.gwf("m", context=ModelContext(grid=grid), packages=[
        mf.disv(_passthrough_stack(), attach=False),
        mf.npf(k=10.0, icelltype=1),
        mf.ic(strt=95.0),
        mf.chd(stress_period_data={0: [[(2, c), 90.0] for c in east]}),
        mf.rch.array(recharge=rate, fixed_cell=fixed_cell),
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]), mf.ims(models=["m"])])
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        run = mf.Project(tmp_path / "p", name="p").run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])

    workspace = pathlib.Path(run.workspace)
    listing = (workspace / "m.lst").read_text()
    discrepancy = float(
        [line for line in listing.splitlines() if "PERCENT DISCREPANCY" in line][-1].split()[-1]
    )
    record = flopy.utils.CellBudgetFile(str(workspace / "m.cbc"), precision="double").get_data(
        text="RCH")[-1]
    idomain = _passthrough_stack().build(grid, attach=False).idomain
    passthrough = np.flatnonzero(idomain[0] == -1)
    area = np.asarray(grid.gdf_vorPolys.area, dtype=float)
    column = (np.asarray(record["node"]) - 1) % int(grid.ncpl)
    applied = np.asarray(record["q"], dtype=float)

    assert len(passthrough) >= 10
    assert abs(discrepancy) < 0.1
    assert applied[np.isin(column, passthrough)].sum() == pytest.approx(
        rate * area[passthrough].sum(), rel=1e-6
    )
    if not fixed_cell:
        assert applied.sum() == pytest.approx(rate * area.sum(), rel=1e-6)
        return
    head = np.asarray(flopy.utils.HeadFile(str(workspace / "m.hds")).get_data())
    fill = idomain[0] > 0
    dry = fill & (head.reshape(idomain.shape)[0] <= -1e29)
    assert fill.sum() >= 10 and np.array_equal(dry, fill), "the premise: every fill cell dry"
    assert applied.sum() == pytest.approx(rate * (area.sum() - area[dry].sum()), rel=1e-6)
