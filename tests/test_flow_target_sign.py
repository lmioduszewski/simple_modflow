"""A flow target is compared in MF6's own sign, on both routes (ledger 350, the user's rule).

``SfrFlowTargets`` read one reach with two signs. The observation route handed on
MF6's ``DOWNSTREAM-FLOW``, which is NEGATIVE (water leaving the reach), while the
budget-file route built the same quantity from the SFR budget file and took its
MAGNITUDE -- so the name a target carried chose the sign, and a gauged flow ``Q``
matched exactly on the observation route scored a residual of ``-2Q`` (phi 327,106
on three exact observations of the scoping toy). ``DrnFlowTargets`` read MF6's
negative sign on both routes but said nowhere what to enter, and a positive spring
discharge scored ``-2Q`` on either (phi 7.9e6). The user's decisions: both flows are
MF6's own number on both routes; a gauged ``Q`` is entered as ``-Q``; a target of
the sign MF6 can never write (above zero) is REFUSED when the targets are built,
naming ``-Q``, and zero is accepted.

One real MF6 run of the toy of ``test_flow_targets_period_end.py`` -- the 6 x 6
irregular Voronoi grid, one confined layer, CHD at 102 ft on the west and east
edges, heads filling from 90 ft -- with a MOVER: a four-reach stream whose reach 1
passes 30% of its flow to a lake (``FACTOR 0.3``) and whose outlet, reach 3, passes
up to 20 (``UPTO 20``), so the four reaches are a plain interior reach (0, 2), an
interior reach through a mover (1) and an outlet through a mover (3). SFR
``DOWNSTREAM-FLOW`` is declared for every reach under ``r0`` .. ``r3``, DRN for a
three-cell zone under ``seep``: those names take the observation route, any other
name the budget files. Three 10-day periods of three steps, ``SAVE ALL``.

The expected values are read from MF6's own files -- the observation CSVs with
pandas, the budget files with FloPy -- not through the tables under test.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from flopy.utils import CellBudgetFile

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.pest.forward_run import write_named_series_targets

NPER, NSTP, PERLEN = 3, 3, 10.0
REACHES = [13, 14, 15, 16]  # one interior row, west to east
LAST = len(REACHES) - 1
DRN_CELLS = [20, 21, 27]
LAKE_CELLS = [8, 9]
#: The names MF6 observes each reach under -- the observation route.
OBSERVED = {f"r{reach}": reach for reach in range(len(REACHES))}
#: The same reaches under names nothing declares -- the budget-file route.
GAUGED = {f"g{reach}": reach for reach in range(len(REACHES))}


def _build_and_run(root: Path):
    """The toy in the module docstring, run; returns the live run model."""

    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    stream = mf.sfr.flopy(
        packagedata=[[i, (0, c), 100.0, 5.0, 0.001, 95.0 - 0.1 * i, 1.0, 1.0, 0.03,
                      1 if i in (0, LAST) else 2, 1.0, 0] for i, c in enumerate(REACHES)],
        connectiondata=[[i, *([i - 1] if i > 0 else []), *([-(i + 1)] if i < LAST else [])]
                        for i in range(len(REACHES))],
        perioddata={0: [[0, "INFLOW", 50.0]], 1: [[0, "INFLOW", 200.0]],
                    2: [[0, "INFLOW", 20.0]]},
        stage_filerecord="m.sfr.stg", budget_filerecord="m.sfr.bud", mover=True,
        observations=mf.SfrFlowTargets(locations=OBSERVED).to_flopy_obs(
            csv_name="sfr_flow_targets.csv"),
    )
    drain = mf.drn(
        stress_period_data={0: [[(0, c), 96.0, 50.0] for c in DRN_CELLS]},
        # In its layer: with no model there are no layers to read (a zone given none
        # reads every layer, ledger 550).
        observations=mf.DrnFlowTargets(
            locations=[{"name": "seep", "layer": 0, "cells": DRN_CELLS}]
        ).to_flopy_obs(None, csv_name="drn_flow_targets.csv"),
    )
    lake = mf.lak.flopy(
        nlakes=1, packagedata=[[0, 104.0, len(LAKE_CELLS), "pond"]],
        connectiondata=[[0, i, (0, c), "VERTICAL", 0.5, 0.0, 0.0, 0.0, 0.0]
                        for i, c in enumerate(LAKE_CELLS)],
        perioddata={0: [[0, "WITHDRAWAL", 100.0]]}, boundnames=True, mover=True,
        stage_filerecord="m.lak.stg",
    )
    mover = mf.mvr(nper=1, moves=(
        mf.Move(mf.MoverConnection("sfr", 1), mf.MoverConnection("lak", 0), "FACTOR", 0.3),
        mf.Move(mf.MoverConnection("sfr", LAST), mf.MoverConnection("lak", 0), "UPTO", 20.0),
    ))
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=110.0, botm=[0.0]),
        mf.npf(k=10.0, icelltype=0),
        mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=90.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), 102.0] for c in edge]}),
        stream, drain, lake, mover,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=NPER, perioddata=[(PERLEN, NSTP, 1.0)] * NPER)])
    with warnings.catch_warnings():
        # A built grid object cannot be frozen into the run's snapshot, and this run
        # is never reopened.
        warnings.filterwarnings("ignore", "Run 'r' built but its frozen snapshot", UserWarning)
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    return run.model()


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    """The run in the module docstring."""

    return _build_and_run(tmp_path_factory.mktemp("flow_target_sign"))


def _at_period_ends(model, csv: str) -> pd.DataFrame:
    """MF6's observation file at t = 10, 20, 30 -- found by time -- to the last digit."""

    raw = pd.read_csv(Path(model.workspace) / csv, float_precision="round_trip")
    assert len(raw) == NPER * NSTP
    rows = [int(np.flatnonzero(np.isclose(raw["time"], PERLEN * (per + 1), rtol=1e-9))[0])
            for per in range(NPER)]
    return raw.iloc[rows].reset_index(drop=True)


def _budget_terms(model) -> dict[str, dict[int, np.ndarray]]:
    """``{term: {reach: [q at each period end]}}`` from the SFR budget file, with FloPy.

    ``FLOW-JA-FACE`` is the reach's OUTBOUND rows added; every value as MF6 wrote it.
    """

    reader = CellBudgetFile(Path(model.workspace) / "m.sfr.bud", precision="double")
    terms = {"FLOW-JA-FACE": {}, "TO-MVR": {}, "EXT-OUTFLOW": {}}
    for reach in range(len(REACHES)):
        for term in terms:
            values = []
            for per in range(NPER):
                data = reader.get_data(text=term, kstpkper=(NSTP - 1, per))[0]
                q = data["q"][data["node"] == reach + 1]
                values.append(float(q[q < 0].sum()) if term == "FLOW-JA-FACE" else float(q.sum()))
            terms[term][reach] = np.array(values)
    reader.close()
    return terms


# -- one number, on both routes ---------------------------------------------------------


@pytest.mark.slow
def test_both_routes_read_mf6s_downstream_flow_at_every_period_end(toy):
    """The same reach under an observed and an undeclared name: the SAME number now,
    MF6's negative ``DOWNSTREAM-FLOW``, on a plain reach, a reach through a mover and
    an outlet through a mover. The budget route read each with the opposite sign."""

    terms = _budget_terms(toy)
    # The premise: reach 1 sends water on AND to the mover, the outlet out of the
    # network AND to the mover, so every piece of the budget route is exercised.
    assert (terms["TO-MVR"][1] < 0).all() and (terms["FLOW-JA-FACE"][1] < 0).all()
    assert (terms["TO-MVR"][LAST] < 0).all() and (terms["EXT-OUTFLOW"][LAST] < 0).all()
    assert (terms["TO-MVR"][0] == 0).all() and (terms["TO-MVR"][2] == 0).all()

    written = _at_period_ends(toy, "sfr_flow_targets.csv")
    observed = mf.SfrFlowTargets(locations=OBSERVED).simulated_series(toy)
    gauged = mf.SfrFlowTargets(locations=GAUGED).simulated_series(toy)
    assert list(observed["time"]) == list(gauged["time"]) == list(range(NPER))
    for reach in range(len(REACHES)):
        mf6 = written[f"R{reach}"].to_numpy()
        assert (mf6 < 0).all(), reach  # MF6's frame: water leaving the reach
        np.testing.assert_array_equal(observed[f"r{reach}"], mf6, err_msg=str(reach))
        np.testing.assert_allclose(gauged[f"g{reach}"], mf6, rtol=1e-10, err_msg=str(reach))
        # ... and it IS the budget file's pieces, each leaving the reach, added.
        pieces = terms["FLOW-JA-FACE"][reach] + terms["TO-MVR"][reach] + terms["EXT-OUTFLOW"][reach]
        np.testing.assert_allclose(gauged[f"g{reach}"], pieces, rtol=1e-12, err_msg=str(reach))


@pytest.mark.slow
def test_a_gauged_flow_entered_as_minus_q_reads_a_zero_residual_on_both_routes(toy):
    """A flow measured exactly as the model makes it, entered as ``-Q``: residual 0 on
    every route. Entered as ``+Q`` it scored ``-2Q`` on the observation route and 0 on
    the budget route (SFR), ``-2Q`` on both (DRN); it is refused now."""

    stream = mf.SfrFlowTargets(locations=OBSERVED).simulated_series(toy)
    seep = mf.DrnFlowTargets(locations={"seep": DRN_CELLS}).simulated_series(toy)
    for route, names in (("observed", OBSERVED), ("budget", GAUGED)):
        values = pd.DataFrame({"time": range(NPER)})
        for name, observed in zip(names, OBSERVED, strict=True):
            values[name] = stream[observed].to_numpy()
        compare = mf.SfrFlowTargets(locations=names, values=values).compare(toy)
        assert len(compare) == NPER * len(REACHES), route
        np.testing.assert_allclose(compare["residual"], 0.0, atol=1e-10, err_msg=route)
    for name in ("seep", "zone"):  # "zone" is declared nowhere: the budget file
        compare = mf.DrnFlowTargets(
            locations={name: DRN_CELLS},
            values=pd.DataFrame({"time": range(NPER), name: seep["seep"].to_numpy()}),
        ).compare(toy)
        assert (compare["flow_target"] < 0).all()
        np.testing.assert_allclose(compare["residual"], 0.0, atol=1e-10, err_msg=name)


@pytest.mark.slow
def test_pest_is_handed_mf6s_sign_on_both_routes(toy, tmp_path):
    """The build writes each target as ``obsval`` beside a base-run table in MF6's sign,
    so a target the model reproduces carries a zero residual on every route -- phi
    was ``sum (2Q)^2`` for an observation-route SFR target -- and every forward run
    writes the base run's table again."""

    pyemu = pytest.importorskip("pyemu")
    stream = mf.SfrFlowTargets(locations=OBSERVED).simulated_series(toy)
    seep = mf.DrnFlowTargets(locations={"seep": DRN_CELLS}).simulated_series(toy)
    sets = {
        "sfr_obs": mf.SfrFlowTargets(locations={"r1": 1, "r3": LAST}, values=pd.DataFrame(
            {"time": range(NPER), "r1": stream["r1"], "r3": stream["r3"]})),
        "sfr_bud": mf.SfrFlowTargets(locations={"g1": 1, "g3": LAST}, values=pd.DataFrame(
            {"time": range(NPER), "g1": stream["r1"], "g3": stream["r3"]})),
        "drn_obs": mf.DrnFlowTargets(locations={"seep": DRN_CELLS}, values=pd.DataFrame(
            {"time": range(NPER), "seep": seep["seep"]})),
        "drn_bud": mf.DrnFlowTargets(locations={"zone": DRN_CELLS}, values=pd.DataFrame(
            {"time": range(NPER), "zone": seep["seep"]})),
    }
    cal = toy.pest("sign", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    for prefix, targets in sets.items():
        cal.observe(targets, prefix=prefix)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pst = cal.build("sign.pst", noptmax=0)
    template = Path(cal.template_workspace)
    obs = pst.observation_data
    weighted = obs.loc[obs["weight"] > 0]
    assert len(weighted) == 2 * 2 * NPER + 2 * 1 * NPER  # two SFR sets of two, two DRN of one

    read = {}
    for prefix in sets:
        kind = "sfr_flow" if prefix.startswith("sfr") else "drn_flow"
        table = template / f"{prefix}_simulated_{kind}.csv"
        read.update(pyemu.pst_utils.InstructionFile(f"{table}.ins").read_output_file(
            str(table))["obsval"].to_dict())
        # A forward run writes the base run's table again, sign and all.
        again = tmp_path / f"{prefix}_forward.csv"
        write_named_series_targets(
            [{"kind": kind, "locations_file": template / f"{prefix}_target_locations.csv",
              "output_csv": again}],
            sim_ws=toy.workspace, model_name="m")
        assert again.read_text() == table.read_text(), prefix
    assert (weighted["obsval"] < 0).all()
    simulated = np.array([read[name] for name in weighted.index], dtype=float)
    np.testing.assert_allclose(simulated, weighted["obsval"].to_numpy(dtype=float), rtol=1e-12)


# -- the impossible sign is refused -------------------------------------------------------


@pytest.mark.parametrize(
    ("make", "locations", "refusal"),
    [
        (mf.SfrFlowTargets, {"gage": 3},
         r"SfrFlowTargets values are MF6's DOWNSTREAM-FLOW, in MF6's own sign: MF6 writes "
         r"water leaving the reach as negative, so the flow is never above 0, and 2 value\(s\) "
         r"are: 'gage' at time 0 holds 165\.7; 'gage' at time 2 holds 0\.5\. Enter a gauged "
         r"streamflow Q as -Q \(0, a dry reading, is accepted\)\."),
        (mf.DrnFlowTargets, {"gage": [20, 21]},
         r"DrnFlowTargets values are MF6's DRN flow, in MF6's own sign: MF6 writes water "
         r"leaving the aquifer through a drain as negative, so the flow is never above 0, and "
         r"2 value\(s\) are: 'gage' at time 0 holds 165\.7; 'gage' at time 2 holds 0\.5\. "
         r"Enter a measured discharge Q as -Q \(0, a dry reading, is accepted\)\."),
    ],
    ids=["sfr", "drn"],
)
def test_a_flow_target_above_zero_is_refused_naming_minus_q(make, locations, refusal):
    """Both flows are water LEAVING -- the reach, the aquifer -- which MF6 writes
    negative, so a value above zero can never be matched: it is a sign error, and it
    scored ``-2Q``. Refused in every shape the values take; zero and a blank read."""

    wide = pd.DataFrame({"time": [0, 1, 2, 3], "gage": [165.7, -3.0, 0.5, np.nan]})
    long = wide.melt(id_vars="time", var_name="name", value_name="flow")
    series = wide.set_index("time")["gage"]
    for values in (wide, long, series, wide.to_dict(orient="list")):
        with pytest.raises(ValueError, match=refusal):
            make(locations=locations, values=values)

    accepted = make(locations=locations, values=pd.DataFrame(
        {"time": [0, 1, 2], "gage": [0.0, -3.0, np.nan]}))
    assert accepted.to_wide()["gage"].tolist()[:2] == [0.0, -3.0]
    # Six or more are named five at a time, with the rest counted.
    with pytest.raises(ValueError, match=r"7 value\(s\) are: .*'gage' at time 4 holds 1\.0; and 2 more\."):
        make(locations=locations, values=pd.DataFrame({"time": range(7), "gage": 1.0}))


def test_a_model_without_outputs_is_not_read_through_its_exchange():
    """A branch for doubles without ``outputs`` read ``sfr.results.q`` -- the
    stream-AQUIFER exchange -- as the reach's streamflow, and took its magnitude. The
    budget route reads the SFR budget file or nothing."""

    exchange = pd.DataFrame({"per": [0, 1], "reach": [0, 0], "q_gwf": [-123.0, 456.0]})
    double = SimpleNamespace(
        sfr=SimpleNamespace(results=SimpleNamespace(q=SimpleNamespace(get=lambda: exchange))))
    with pytest.raises(AttributeError, match="outputs"):
        mf.SfrFlowTargets(locations={"gage": 0}).simulated_series(double)
