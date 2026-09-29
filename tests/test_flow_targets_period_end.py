"""A period's flow is its END-of-period flow (ledger 336, the user's rule).

``SfrFlowTargets`` and ``DrnFlowTargets`` ADDED every time step a period saved on
their budget-file route: at ``nstp=3`` under ``SAVE BUDGET ALL`` the canonical
model's outlet reach read 113,453.29 where each step held 37,817.76, and a
one-cell DRN zone -292.51 against -97.50. Summing a zone's CELLS within one step
is right -- that is what a seepage zone is -- and summing its time steps is not.
Both flow families, and the head and stage families beside them, now read the
record that ENDS each period through one rule (``myflopy.modflow.mf6.period_end``),
on the clock of the file the value comes from, and a period whose end that file
does not hold is refused.

Five small real MF6 runs of ONE model -- the 6 x 6 irregular Voronoi grid of
``test_stage_targets_one_stage.py``, one confined layer, CHD at 102 ft on the
west and east edges, heads filling from 90 ft so every flow moves within a
period; a four-reach stream (``up`` = reach 0, ``down`` = reach 3, the outlet,
which has no downstream reach), a three-cell DRN zone ``seep`` and a two-cell
lake; three 10-day periods of THREE steps each. SFR ``DOWNSTREAM-FLOW`` and
per-cell DRN observations are written for the names ``up``/``down``/``seep``, so
those names take the observation route and any other name the budget file.

``all``           ``SAVE HEAD ALL`` + ``SAVE BUDGET ALL``.
``budget_last``   ``SAVE HEAD ALL`` + ``SAVE BUDGET LAST``: 9 head records against
                  3 budget records, the two clocks disagreeing (ledger 301).
``budget_first``  ``SAVE HEAD LAST`` + ``SAVE BUDGET FIRST``: no period's END in
                  the budget file.
``head_first``    ``SAVE HEAD FIRST`` + ``SAVE BUDGET LAST``: no period's END in
                  the head or stage files.
``ats``           ``SAVE ALL`` under ATS, one declared step per period: MF6 decides
                  the step count at run time, so the END is found by time, not by
                  ``nstp``.

A sixth run, ``hourly``, is ``all`` with three ONE-HOUR periods (``1/24`` d), read
live: FloPy writes that length at 8 decimals, so the period end MF6 records is not
the length the live model holds.

The expected values are read from the budget files with FloPy directly, not
through the myflopy tables under test, and each fixture asserts what its files
hold before it is used.
"""

from __future__ import annotations

import ast
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from flopy.utils import CellBudgetFile, HeadFile

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.observations._shared import _aggregate_named_series
from myflopy.modflow.mf6.period_end import (
    at_period_end,
    period_end_steps,
    period_lengths,
    saved_steps,
)
from myflopy.modflow.mf6.pest.forward_run import (
    _write_conc_target_csv,
    _write_head_target_csv,
    write_named_series_targets,
)
from myflopy.modflow.mf6.pest.observations import (
    prepare_drn_flow_observations,
    prepare_sfr_flow_observations,
)
from myflopy.modflow.mf6.pest.specs import DrnFlowObservationSpec, SfrFlowObservationSpec

NPER, NSTP, PERLEN = 3, 3, 10.0
REACHES = [13, 14, 15, 16]  # one interior row, west to east
DRN_CELLS = [20, 21, 27]
LAKE_CELLS = [8, 9]
SAVES = {
    "all": [("HEAD", "ALL"), ("BUDGET", "ALL")],
    "budget_last": [("HEAD", "ALL"), ("BUDGET", "LAST")],
    "budget_first": [("HEAD", "LAST"), ("BUDGET", "FIRST")],
    "head_first": [("HEAD", "FIRST"), ("BUDGET", "LAST")],
    "ats": [("HEAD", "ALL"), ("BUDGET", "ALL")],
}


def _build_and_run(root: Path, saverecord, *, ats: bool = False, perlen: float = PERLEN):
    """The model in the module docstring, run under ``saverecord``, periods ``perlen`` long."""

    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    last = len(REACHES) - 1
    packagedata = [
        [i, (0, c), 100.0, 5.0, 0.001, 95.0 - 0.1 * i, 1.0, 1.0, 0.03,
         1 if i in (0, last) else 2, 1.0, 0]
        for i, c in enumerate(REACHES)
    ]
    connectiondata = [
        [i, *([i - 1] if i > 0 else []), *([-(i + 1)] if i < last else [])]
        for i in range(len(REACHES))
    ]
    stream = mf.sfr.flopy(
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata={0: [[0, "INFLOW", 50.0]], 1: [[0, "INFLOW", 200.0]],
                    2: [[0, "INFLOW", 20.0]]},
        stage_filerecord="m.sfr.stg",
        budget_filerecord="m.sfr.bud",
        observations=mf.SfrFlowTargets(locations={"up": 0, "down": last}).to_flopy_obs(
            csv_name="sfr_flow_targets.csv"),
    )
    drain = mf.drn(
        stress_period_data={0: [[(0, c), 96.0, 50.0] for c in DRN_CELLS]},
        observations=mf.DrnFlowTargets(locations={"seep": DRN_CELLS}).to_flopy_obs(
            None, csv_name="drn_flow_targets.csv"),
    )
    lake = mf.lak.flopy(
        nlakes=1,
        packagedata=[[0, 104.0, len(LAKE_CELLS), "pond"]],
        connectiondata=[[0, i, (0, c), "VERTICAL", 0.5, 0.0, 0.0, 0.0, 0.0]
                        for i, c in enumerate(LAKE_CELLS)],
        perioddata={0: [[0, "WITHDRAWAL", 100.0]]},
        boundnames=True,
        stage_filerecord="m.lak.stg",
    )
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=110.0, botm=[0.0]),
        mf.npf(k=10.0, icelltype=0),
        mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=90.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), 102.0] for c in edge]}),
        stream,
        drain,
        lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc", saverecord=saverecord),
    ])
    timing = {"ats": True, "ats_dt0": 1.0, "ats_dtmin": 0.01, "ats_dtmax": 4.0} if ats else {}
    steps = 1 if ats else NSTP
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=NPER, perioddata=[(perlen, steps, 1.0)] * NPER, **timing)])
    with warnings.catch_warnings():
        # A built grid object cannot be frozen into the run's snapshot, and this run
        # is never reopened.
        warnings.filterwarnings("ignore", "Run 'r' built but its frozen snapshot", UserWarning)
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    return run.model()


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    """The five runs in the module docstring, by name."""

    root = tmp_path_factory.mktemp("flow_period_end")
    return {name: _build_and_run(root / name, save, ats=name == "ats")
            for name, save in SAVES.items()}


#: One hour, in days. FloPy writes it ``0.04166667`` -- 8 decimals, 3.3e-9 short of
#: the ``1/24`` a live model holds -- and MF6 ends the period at what it READ.
HOUR = 1.0 / 24.0


@pytest.fixture(scope="module")
def hourly(tmp_path_factory):
    """The ``all`` run with three ONE-HOUR periods, read live (``run.model()``)."""

    return _build_and_run(tmp_path_factory.mktemp("hourly"), SAVES["all"], perlen=HOUR)


def _budget_flows(model):
    """``{(kstp, kper): {"up", "down", "seep": flow}}`` read with FloPy, zero-based.

    A reach's flow is what it sends downstream -- outbound FLOW-JA-FACE plus
    EXT-OUTFLOW, the quantity MF6's DOWNSTREAM-FLOW reports -- as a magnitude;
    the zone's is its cells' DRN rows added, in MF6's sign.
    """

    workspace = Path(model.workspace)
    sfr = CellBudgetFile(workspace / "m.sfr.bud", precision="double")
    cbc = CellBudgetFile(workspace / "m.cbc", precision="double")
    flows = {}
    for kstp, kper in sfr.get_kstpkper():
        key = (int(kstp), int(kper))
        routed = sfr.get_data(text="FLOW-JA-FACE", kstpkper=key)[0]
        outflow = sfr.get_data(text="EXT-OUTFLOW", kstpkper=key)[0]
        by_reach = []
        for reach in (0, len(REACHES) - 1):
            rows = routed[routed["node"] == reach + 1]
            out = -rows["q"][rows["q"] < 0].sum() - outflow["q"][outflow["node"] == reach + 1].sum()
            by_reach.append(float(out))
        drains = cbc.get_data(text="DRN", kstpkper=key)[0]
        seep = float(drains["q"][np.isin(drains["node"], np.array(DRN_CELLS) + 1)].sum())
        flows[key] = {"up": by_reach[0], "down": by_reach[1], "seep": seep}
    sfr.close()
    cbc.close()
    return flows


def _period_end(flows, name):
    """``name``'s flow at each period's last step, from :func:`_budget_flows`."""

    last = {kper: max(k for k, p in flows if p == kper) for _, kper in flows}
    return [flows[(last[per], per)][name] for per in range(NPER)]


def _targets(values=None):
    """Flow targets under names the observation CSVs do not hold: the budget route."""

    sfr = mf.SfrFlowTargets(
        locations={"a": 0, "b": len(REACHES) - 1},
        values=None if values is None else pd.DataFrame(
            {"time": range(NPER), "a": values, "b": values}),
    )
    drn = mf.DrnFlowTargets(
        locations={"zone": DRN_CELLS},
        values=None if values is None else pd.DataFrame(
            {"time": range(NPER), "zone": values}),
    )
    return sfr, drn


# -- the rule -------------------------------------------------------------------


def test_the_rule_takes_the_record_that_ends_each_period():
    """By time, not by ``nstp``; a period that saved nothing is absent, not wrong."""

    saved = [(0, 0, 3.0), (1, 0, 7.0), (2, 0, 10.0), (6, 1, 10.0), (0, 3, 2.5), (1, 3, 10.0)]
    assert period_end_steps(saved, perlen=[10.0] * 4) == {0: (2, 0), 1: (6, 1), 3: (1, 3)}
    # Without lengths there is nothing to check: the last saved step, as given.
    assert period_end_steps([(0, 0), (2, 0), (1, 0)]) == {0: (2, 0)}
    table = pd.DataFrame({"kstpkper": [(0, 0), (2, 0), (6, 1)], "q": [1.0, 2.0, 3.0]})
    assert at_period_end(table, {0: (2, 0), 1: (6, 1)})["q"].tolist() == [2.0, 3.0]


def test_a_period_that_saved_no_end_is_refused_naming_it_and_the_fix():
    saved = [(0, 0, 10.0), (0, 1, 2.0), (2, 1, 6.0), (0, 2, 10.0)]
    with pytest.raises(ValueError, match=(
            r"END is missing from the SFR budget file: period 1 saved steps 1, 3, the "
            r"last ending 6 into a period 10 long.*SAVE BUDGET LAST \(or ALL\)")):
        period_end_steps(saved, perlen=[10.0] * 3, file="the SFR budget file", saves="BUDGET")
    with pytest.raises(ValueError, match="TDIS declares 2 periods"):
        period_end_steps(saved, perlen=[10.0] * 2)


def test_a_period_end_within_flopys_write_rounding_is_the_end():
    """MF6 ends a period at the length it READ, which FloPy wrote at 8 decimals; a
    caller holds the length FloPy was GIVEN. At a tolerance of 1e-9 an hourly, a
    1/3- or a 4/3-unit period was refused while the file held its end (ledger 336)."""

    for length in (HOUR, 1 / 3, 4 / 3, 365 / 12):
        written = float(f"{length:15.8f}")  # FloPy's fixed format, 0.001 to 1e5
        assert period_end_steps([(2, 0, written)], perlen=[length]) == {0: (2, 0)}
    written = float(f"{100000.00049:.8E}")  # and its E format above 1e5
    assert period_end_steps([(0, 0, written)], perlen=[100000.00049]) == {0: (0, 0)}
    # A record short of the end by more than that is still refused, and the refusal
    # prints digits enough to tell the two numbers apart: at six figures this one
    # read "ending 10 into a period 10 long".
    with pytest.raises(ValueError, match=(
            r"period 0 saved step 3, the last ending 9\.9999998 into a period 10 long")):
        period_end_steps([(2, 0, 9.9999998)], perlen=[10.0])


def test_time_steps_never_add_even_for_an_extensive_value():
    """``kind="extensive"`` adds the CELLS of one step; the steps are read at the end."""

    table = pd.DataFrame({
        "kstpkper": [(0, 0), (0, 0), (1, 0), (1, 0), (0, 1)],
        "per": [0, 0, 0, 0, 1], "zone": [0] * 5, "q": [-1.0, -2.0, -4.0, -8.0, -16.0]})
    reduced = _aggregate_named_series(table, id_column="zone", value_column="q", kind="extensive")
    assert reduced["q"].tolist() == [-12.0, -16.0]


# -- the budget-file route --------------------------------------------------------


@pytest.mark.slow
def test_a_flow_target_reads_each_periods_end_not_the_sum_of_its_steps(runs):
    """The 3.0x of ledger 331: SFR and DRN now read the budget file's period-end step."""

    model = runs["all"]
    flows = _budget_flows(model)
    assert len(flows) == NPER * NSTP
    for name in ("up", "down", "seep"):
        per_step = [flows[(k, 0)][name] for k in range(NSTP)]
        assert len(set(np.round(per_step, 6))) == NSTP, "the flow must move within a period"

    sfr, drn = _targets()
    stream = sfr.simulated_series(model)
    seep = drn.simulated_series(model)
    assert list(stream["time"]) == list(range(NPER)) == list(seep["time"])
    np.testing.assert_allclose(stream["a"], _period_end(flows, "up"), rtol=1e-12)
    np.testing.assert_allclose(stream["b"], _period_end(flows, "down"), rtol=1e-12)
    np.testing.assert_allclose(seep["zone"], _period_end(flows, "seep"), rtol=1e-12)
    summed = sum(flows[(k, 0)]["seep"] for k in range(NSTP))
    assert not np.isclose(seep["zone"].iloc[0], summed)


@pytest.mark.slow
def test_the_outlet_reach_reads_what_leaves_the_network(runs):
    """Reach 3 has no downstream reach: its water leaves as EXT-OUTFLOW, which the
    route did not read, so the outlet -- the reach a gage sits on -- had no column."""

    series = mf.SfrFlowTargets(locations={"outlet": len(REACHES) - 1}).simulated_series(runs["all"])
    np.testing.assert_allclose(series["outlet"], _period_end(_budget_flows(runs["all"]), "down"),
                               rtol=1e-12)
    assert (series["outlet"] > 0).all()


@pytest.mark.slow
def test_a_flow_is_read_on_the_budget_clock_not_the_head_clock(runs):
    """``SAVE HEAD ALL`` + ``SAVE BUDGET LAST``: 9 head records against 3 budget
    records, and each file's period end is its own."""

    model = runs["budget_last"]
    head_file = HeadFile(Path(model.workspace) / "m.hds", precision="double")
    assert len(head_file.get_kstpkper()) == NPER * NSTP
    head_file.close()
    flows = _budget_flows(model)
    assert sorted(flows) == [(NSTP - 1, per) for per in range(NPER)]

    sfr, drn = _targets()
    np.testing.assert_allclose(sfr.simulated_series(model)["a"], _period_end(flows, "up"),
                               rtol=1e-12)
    np.testing.assert_allclose(drn.simulated_series(model)["zone"], _period_end(flows, "seep"),
                               rtol=1e-12)
    # The same periods' ends as the all-saved run: the schedule chose which records
    # exist, not which value a period has.
    np.testing.assert_allclose(drn.simulated_series(model)["zone"],
                               drn.simulated_series(runs["all"])["zone"], rtol=1e-12)


@pytest.mark.slow
def test_compare_and_the_files_pest_reads_carry_the_period_end_flow(runs, tmp_path):
    """``compare()``, the base run's CSV and a forward run's, together."""

    model = runs["all"]
    flows = _budget_flows(model)
    sfr, drn = _targets(values=[1.0] * NPER)
    compare = sfr.compare(model)
    by_name = {name: frame.sort_values("per") for name, frame in compare.groupby("name")}
    np.testing.assert_allclose(by_name["a"]["sim_flow"], _period_end(flows, "up"), rtol=1e-12)
    np.testing.assert_allclose(
        drn.compare(model).sort_values("per")["sim_flow"], _period_end(flows, "seep"), rtol=1e-12)

    project = SimpleNamespace(model=model, template_workspace=tmp_path,
                              pf=SimpleNamespace(add_observations=lambda *a, **k: None))
    for kind, targets, prepare, spec, column, name in (
        ("sfr_flow", sfr, prepare_sfr_flow_observations, SfrFlowObservationSpec, "a", "up"),
        ("drn_flow", drn, prepare_drn_flow_observations, DrnFlowObservationSpec, "zone", "seep"),
    ):
        prepared = prepare(project, spec(targets=targets, prefix=kind))
        config = prepared["named_series_forward_run_config"]
        base = tmp_path / config["output_csv"]
        forward = tmp_path / f"{kind}_forward.csv"
        write_named_series_targets(
            kind, tmp_path / config["locations_file"], forward, sim_ws=model.workspace)
        assert forward.read_text() == base.read_text()
        np.testing.assert_allclose(pd.read_csv(base)[column], _period_end(flows, name),
                                   rtol=1e-12)


@pytest.mark.slow
def test_a_period_whose_end_the_budget_did_not_save_is_refused(runs, tmp_path):
    """``SAVE BUDGET FIRST``: each period's flow was read at step 1 of 3, silently."""

    model = runs["budget_first"]
    assert sorted(_budget_flows(model)) == [(0, per) for per in range(NPER)]
    sfr, drn = _targets(values=[1.0] * NPER)
    refusal = (r"END is missing from the (SFR budget|model's cell-budget) file: period 0 "
               r"saved step 1, the last ending 3\.33333 into a period 10 long.*"
               r"SAVE BUDGET LAST")
    for call in (sfr.simulated_series, sfr.compare, drn.simulated_series, drn.compare):
        with pytest.raises(ValueError, match=refusal):
            call(model)
    project = SimpleNamespace(model=model, template_workspace=tmp_path,
                              pf=SimpleNamespace(add_observations=lambda *a, **k: None))
    with pytest.raises(ValueError, match=refusal):
        prepare_sfr_flow_observations(project, SfrFlowObservationSpec(targets=sfr, prefix="q"))

    # Its head file saved each period's END (SAVE HEAD LAST), so heads read.
    heads = mf.HeadTargets(
        locations=pd.DataFrame({"name": ["h"], "cell": [22], "layer": [0]}),
        values=pd.DataFrame({"time": range(NPER), "name": "h", "head": 100.0}))
    assert list(heads.simulated_heads(model)["per"]) == list(range(NPER))


#: ``SAVE HEAD FIRST``: the same decision for every family (ledger 336). Heads, SFR
#: stage and LAK stage are written when MF6 saves HEAD, so none holds a period end.
HEAD_REFUSAL = r"END is missing from {}: period 0 saved step 1.*SAVE HEAD LAST"


@pytest.mark.slow
def test_a_head_target_is_refused_when_the_head_file_saved_no_period_end(runs):
    heads = mf.HeadTargets(
        locations=pd.DataFrame({"name": ["h"], "cell": [22], "layer": [0]}),
        values=pd.DataFrame({"time": range(NPER), "name": "h", "head": 100.0}))
    with pytest.raises(ValueError, match=HEAD_REFUSAL.format("the head file")):
        heads.simulated_heads(runs["head_first"])
    with pytest.raises(ValueError, match=HEAD_REFUSAL.format("the head file")):
        heads.compare(runs["head_first"])


@pytest.mark.slow
def test_the_head_noun_reads_each_periods_end_and_refuses_a_file_without_one(runs):
    """``model.hds.get()`` documents period-end values, and under ``SAVE HEAD
    FIRST`` it returned step 1 of 3 as each period's head (97.61 where period 0
    ends at 101.18 on the ``all`` run). It reads the one rule now (ledger 336)."""

    model = runs["all"]
    table = model.all_heads.reset_index()
    table = table[(table["layer"] == 0) & (table["cell"] == 22)]
    at = {tuple(int(v) for v in key): float(value)
          for key, value in zip(table["kstpkper"], table["elev"], strict=True)}
    got = model.hds.get(cells=[22], layer=0)
    assert got["per"].tolist() == list(range(NPER))
    np.testing.assert_allclose(got["head"], [at[(NSTP - 1, per)] for per in range(NPER)])
    for call in (runs["head_first"].hds.get, runs["head_first"].hds.summary):
        with pytest.raises(ValueError, match=HEAD_REFUSAL.format("the head file")):
            call()


@pytest.mark.slow
def test_an_sfr_stage_target_is_refused_when_the_stage_file_saved_no_period_end(runs):
    targets = mf.SfrStageTargets(locations={"a": 0},
                                 values=pd.DataFrame({"time": range(NPER), "a": 95.0}))
    for call in (targets.simulated_series, targets.compare):
        with pytest.raises(ValueError, match=HEAD_REFUSAL.format("the SFR stage file")):
            call(runs["head_first"])


@pytest.mark.slow
def test_a_lake_stage_target_is_refused_when_the_stage_file_saved_no_period_end(runs):
    """One record per period made the record index look like the period, with every
    value read at step 1 of 3. The record clock itself stays (ledger 331's limit)."""

    targets = mf.LakeStageTargets(locations={"pond": 0},
                                  values=pd.DataFrame({"time": range(NPER), "pond": 100.0}))
    for call in (targets.simulated_series, targets.compare):
        with pytest.raises(ValueError, match=HEAD_REFUSAL.format("the LAK stage file")):
            call(runs["head_first"])


@pytest.mark.slow
def test_the_flows_beside_them_read_their_own_clock(runs):
    """The same run's BUDGET saved each period's end (``SAVE BUDGET LAST``)."""

    model = runs["head_first"]
    sfr, drn = _targets()
    flows = _budget_flows(model)
    np.testing.assert_allclose(drn.simulated_series(model)["zone"], _period_end(flows, "seep"),
                               rtol=1e-12)
    np.testing.assert_allclose(sfr.simulated_series(model)["a"], _period_end(flows, "up"),
                               rtol=1e-12)


@pytest.mark.slow
def test_the_forward_run_refuses_a_head_file_with_no_period_end(runs, tmp_path, monkeypatch):
    """PEST's forward run reads the same rule: with the period lengths the build
    records it refuses too; a forward run generated before them reads unchecked."""

    mapping = tmp_path / "map.csv"
    pd.DataFrame({"obsname": ["h"], "name": ["h"], "layer": [0], "cell": [22]}).to_csv(
        mapping, index=False)
    monkeypatch.chdir(runs["head_first"].workspace)  # it reads <model>.hds from the cwd
    with pytest.raises(ValueError, match=HEAD_REFUSAL.format(r"m\.hds")):
        _write_head_target_csv("m", str(mapping), str(tmp_path / "h.csv"),
                               perlen=[PERLEN] * NPER)
    _write_head_target_csv("m", str(mapping), str(tmp_path / "h.csv"))
    assert list(pd.read_csv(tmp_path / "h.csv")["per"]) == list(range(NPER))


@pytest.mark.parametrize(
    ("writer", "file"),
    [(_write_head_target_csv, r"flow\.hds"), (_write_conc_target_csv, r"flow\.ucn")],
    ids=["heads", "concentration"],
)
def test_both_forward_run_writers_read_the_period_end_through_the_rule(
    writer, file, tmp_path, monkeypatch
):
    """Both kept a copy of "the last saved step"; both read the one rule now."""

    import flopy

    read = []

    class _FakeFile:
        recordarray = None  # set per case, one-based as on disk

        def __init__(self, *_args, **_kwargs):
            pass

        def get_data(self, kstpkper=None):
            read.append(kstpkper)
            return np.full((1, 4), 7.0)

        def close(self):
            pass

    monkeypatch.setattr(flopy.utils, "HeadFile", _FakeFile)
    monkeypatch.chdir(tmp_path)
    mapping = tmp_path / "map.csv"
    pd.DataFrame({"obsname": ["w"], "name": ["w"], "layer": [0], "cell": [2]}).to_csv(
        mapping, index=False)
    dtype = [("kstp", int), ("kper", int), ("pertim", float)]

    _FakeFile.recordarray = np.array([(1, 1, 4.0), (3, 1, 10.0)], dtype=dtype)
    writer("flow", str(mapping), str(tmp_path / "out.csv"), perlen=[10.0])
    assert read == [(2, 0)]

    _FakeFile.recordarray = np.array([(1, 1, 4.0)], dtype=dtype)
    with pytest.raises(ValueError, match=f"END is missing from {file}: period 0 saved step 1"):
        writer("flow", str(mapping), str(tmp_path / "out.csv"), perlen=[10.0])


@pytest.mark.slow
def test_the_forward_run_reads_the_period_end_head_the_build_did(runs, tmp_path, monkeypatch):
    """The forward run kept its own copy of the rule; it reads the shared one now."""

    model = runs["budget_last"]
    heads = mf.HeadTargets(
        locations=pd.DataFrame({"name": ["h"], "cell": [22], "layer": [0]}),
        values=pd.DataFrame({"time": range(NPER), "name": "h", "head": 100.0}))
    built = heads.simulated_heads(model)
    mapping = tmp_path / "map.csv"
    heads.observation_labels(model).loc[:, ["obsname", "name", "layer", "cell"]].to_csv(
        mapping, index=False)
    monkeypatch.chdir(model.workspace)
    _write_head_target_csv("m", str(mapping), str(tmp_path / "h.csv"), perlen=[PERLEN] * NPER)
    written = pd.read_csv(tmp_path / "h.csv")
    np.testing.assert_allclose(written["h"], built["h"].astype(float), rtol=1e-12)
    table = model.all_heads.reset_index()
    table = table[(table["layer"] == 0) & (table["cell"] == 22)]
    at = {tuple(int(v) for v in key): float(value)
          for key, value in zip(table["kstpkper"], table["elev"], strict=True)}
    np.testing.assert_allclose(written["h"], [at[(NSTP - 1, per)] for per in range(NPER)])


@pytest.mark.slow
def test_an_ats_period_ends_by_time_not_by_its_declared_step_count(runs):
    """One declared step per period, and MF6 takes several: ``kstp == nstp - 1`` would
    refuse every period, and the last saved step is the one whose ``pertim`` is
    the period's length."""

    model = runs["ats"]
    reader = CellBudgetFile(Path(model.workspace) / "m.cbc", precision="double")
    saved = saved_steps(reader)
    reader.close()
    counts = {per: sum(1 for _, kper, _ in saved if kper == per) for per in range(NPER)}
    assert all(count > 1 for count in counts.values()), counts
    ends = {kper: (kstp, pertim) for kstp, kper, pertim in saved}
    assert all(pertim == PERLEN for _, pertim in ends.values())

    flows = _budget_flows(model)
    sfr, drn = _targets()
    np.testing.assert_allclose(drn.simulated_series(model)["zone"], _period_end(flows, "seep"),
                               rtol=1e-12)
    np.testing.assert_allclose(sfr.simulated_series(model)["b"], _period_end(flows, "down"),
                               rtol=1e-12)


# -- a live model's period lengths ---------------------------------------------------


@pytest.mark.slow
def test_a_live_model_with_hourly_periods_reads_each_periods_end(hourly, tmp_path, monkeypatch):
    """The live model holds ``1/24``, the file ends each period at ``0.04166667``:
    every family and the forward run given the lengths the PEST build records
    refused all three periods at a tolerance of 1e-9, while the same run reloaded
    read fine (ledger 336)."""

    model = hourly
    lengths = period_lengths(model)
    assert lengths == [HOUR] * NPER  # what FloPy was given ...
    heads_file = HeadFile(Path(model.workspace) / "m.hds", precision="double")
    ends = [pertim for kstp, _, pertim in saved_steps(heads_file) if kstp == NSTP - 1]
    assert ends == [0.04166667] * NPER  # ... is not what MF6 read
    assert ends != lengths
    expected = [float(heads_file.get_data(kstpkper=(NSTP - 1, per)).ravel()[22])
                for per in range(NPER)]
    heads_file.close()

    heads = mf.HeadTargets(
        locations=pd.DataFrame({"name": ["h"], "cell": [22], "layer": [0]}),
        values=pd.DataFrame({"time": range(NPER), "name": "h", "head": 100.0}))
    np.testing.assert_allclose(heads.simulated_heads(model)["h"].astype(float), expected,
                               rtol=1e-12)
    np.testing.assert_allclose(model.hds.get(cells=[22], layer=0)["head"], expected,
                               rtol=1e-12)
    sfr, drn = _targets()
    flows = _budget_flows(model)
    np.testing.assert_allclose(drn.simulated_series(model)["zone"], _period_end(flows, "seep"),
                               rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(sfr.simulated_series(model)["a"], _period_end(flows, "up"),
                               rtol=1e-12)

    # The forward run, given the lengths exactly as the build writes them into its call.
    mapping = tmp_path / "map.csv"
    heads.observation_labels(model).loc[:, ["obsname", "name", "layer", "cell"]].to_csv(
        mapping, index=False)
    monkeypatch.chdir(model.workspace)
    recorded = ast.literal_eval(repr(lengths))  # the literal the build writes
    _write_head_target_csv("m", str(mapping), str(tmp_path / "h.csv"), perlen=recorded)
    np.testing.assert_allclose(pd.read_csv(tmp_path / "h.csv")["h"], expected, rtol=1e-12)


# -- the observation route, keyed by period since ledger 346 -------------------------


@pytest.mark.slow
def test_the_observation_route_is_keyed_by_period_and_agrees_at_each_period_end(runs):
    """MF6 writes a flow observation at every time step, keyed by simulated time; the
    route hands on the row that ENDS each period, keyed by the period (ledger 346), and
    it equals the budget route -- the DRN zone to rounding, each of its cells to the bit,
    the SFR reach with MF6's own NEGATIVE sign, where the budget route reports the magnitude. It handed on every
    step, so a target keyed by PERIOD matched nothing; one keyed by MF6's time now
    names the period to key it by."""

    model = runs["all"]
    observed_sfr = mf.SfrFlowTargets(locations={"up": 0, "down": len(REACHES) - 1})
    observed_drn = mf.DrnFlowTargets(locations={"seep": DRN_CELLS})
    stream = observed_sfr.simulated_series(model)
    seep = observed_drn.simulated_series(model)
    assert list(stream["time"]) == list(seep["time"]) == list(range(NPER))
    sfr, drn = _targets()
    np.testing.assert_allclose(-stream["up"], sfr.simulated_series(model)["a"], rtol=1e-10)
    np.testing.assert_allclose(-stream["down"], sfr.simulated_series(model)["b"], rtol=1e-10)
    # Each cell is MF6's number on both routes; the zone's sum is added in another order.
    np.testing.assert_allclose(seep["seep"], drn.simulated_series(model)["zone"], rtol=1e-14)

    by_period = mf.DrnFlowTargets(
        locations={"seep": DRN_CELLS},
        values=pd.DataFrame({"time": range(NPER), "seep": -800.0})).compare(model)
    np.testing.assert_array_equal(by_period["sim_flow"], seep["seep"])
    by_time = mf.DrnFlowTargets(
        locations={"seep": DRN_CELLS},
        values=pd.DataFrame({"time": [10, 20, 30], "seep": -800.0}))
    with pytest.raises(ValueError, match=(
            r"DrnFlowTargets.compare\(...\): the observation route is keyed by zero-based "
            r"stress period.*periods are 0 \.\. 2; 3 target row\(s\) name none of them.*"
            r"10 is the end of period 0 in simulated time -- key it 0; 20 is the end of "
            r"period 1 in simulated time -- key it 1")):
        by_time.compare(model)


@pytest.mark.slow
def test_a_zone_given_as_several_rows_sums_the_same_cells_on_both_routes(runs):
    """A zone of rows ``[20]`` and ``[21, 27]`` is ONE zone of three cells. The
    budget route united the rows; the observation route kept the LAST row's
    labels, summing cells 21 and 27 only (-474.37 where the zone ends period 0 at
    -709.47, ledger 336)."""

    model = runs["all"]
    expected = _period_end(_budget_flows(model), "seep")
    observed = mf.DrnFlowTargets(
        locations=pd.DataFrame({"name": ["seep", "seep"], "cells": [[20], [21, 27]]}),
    ).simulated_series(model)
    assert list(observed["time"]) == list(range(NPER))  # one row per period (ledger 346)
    np.testing.assert_allclose(observed["seep"], expected, rtol=1e-10)
    budget = mf.DrnFlowTargets(
        locations=pd.DataFrame({"name": ["zone", "zone"], "cells": [[20], [21, 27]]}),
    ).simulated_series(model)
    np.testing.assert_allclose(budget["zone"], expected, rtol=1e-12)
