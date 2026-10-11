"""The observation route is one row per stress period, at its end (ledger 346).

MF6 writes a continuous observation at the end of EVERY time step, labelled by the
simulated time ``totim`` and carrying no step or period number, and
``LakeStageTargets``, ``SfrStageTargets``, ``SfrFlowTargets`` and ``DrnFlowTargets``
handed every one of those rows on, keyed by that time. So on this route alone:

* a target keyed by stress period -- the key every other route reads -- matched
  nothing: ``compare()`` answered ``NaN``, and the PEST build refused it (ledger 340;
  before that it dropped it and anchored PEST to the base run);
* a PEST base run and its forward runs wrote ``nper x nstp`` rows, and under ATS the
  count moved with the parameters (12, 13, 14 rows at K = 10, 1, 2), while PEST's
  instruction file reads the BASE run's rows by position -- a K = 2 realization read
  t = 23 as the end of the last period, silently.

Each family reads one row per period now, the row that ENDS it, keyed ``0 .. nper - 1``
(``period_end.period_end_rows``, checked by the rule every route reads), and a target
keyed by MF6's time is refused naming the period to key it by.

Real MF6 runs of ONE toy -- the 6 x 6 irregular Voronoi grid of
``test_flow_targets_period_end.py``, one confined layer, CHD at 102 ft on the west and
east edges, heads filling from 90 ft; a four-reach stream, a three-cell DRN zone and a
two-cell lake holding water above its bed (115 ft over a 110 ft bed) and drained by
leakage, so every observed quantity moves within a period. OBS declares all four
families: LAK ``STAGE`` ``pond``, SFR ``STAGE`` ``gage`` (reach 1), SFR
``DOWNSTREAM-FLOW`` ``up``/``down`` (reaches 0 and 3), DRN ``seep``. Three 10-day
periods:

``steps``    three steps each, ``SAVE ALL``: 9 observation rows for 3 periods.
``ats_k10``  ATS, one declared step per period, K = 10 ft/d.
``ats_k2``   the same at K = 2 ft/d: MF6 takes another number of steps.
``hourly``   ``steps`` with SIX one-HOUR periods, read live: FloPy writes 1/24 at 8
             decimals, so MF6's period-end time runs 3.3e-9 further ahead of the
             lengths the live model holds every period -- past ``END_RTOL`` from the
             fourth.
``twelve``   ``steps`` with TWELVE two-day periods: more periods than the time of the
             early period ends, so a key of MF6's time there is also a period number
             (the review of ledger 346, ledger 347).

The expected values are read from MF6's own files -- the observation CSVs with pandas,
the budget and stage files with FloPy -- not through the tables under test.
"""

from __future__ import annotations

import shutil
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from flopy.utils import CellBudgetFile

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.period_end import END_RTOL, period_end_rows, period_lengths
from myflopy.modflow.mf6.pest.forward_run import write_named_series_targets

NPER, NSTP, PERLEN = 3, 3, 10.0
HOURLY_NPER = 6
HOUR = 1.0 / 24.0
REACHES = [13, 14, 15, 16]
DRN_CELLS = [20, 21, 27]
LAKE_CELLS = [8, 9]
#: A value no simulated table holds, written into one cell of a probe copy.
PROBE = 123456.789


def _build_and_run(root: Path, *, ats=False, k=10.0, perlen=PERLEN, nper=NPER):
    """The toy in the module docstring, run; returns the live run model."""

    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    last = len(REACHES) - 1
    stream_obs = mf.SfrFlowTargets(locations={"up": 0, "down": last}).to_flopy_obs(
        csv_name="sfr_flow_targets.csv")
    stream_obs.update(mf.SfrStageTargets(locations={"gage": 1}).to_flopy_obs(
        csv_name="sfr_stage_targets.csv"))
    stream = mf.sfr.flopy(
        packagedata=[[i, (0, c), 100.0, 5.0, 0.001, 95.0 - 0.1 * i, 1.0, 1.0, 0.03,
                      1 if i in (0, last) else 2, 1.0, 0] for i, c in enumerate(REACHES)],
        connectiondata=[[i, *([i - 1] if i > 0 else []), *([-(i + 1)] if i < last else [])]
                        for i in range(len(REACHES))],
        perioddata={0: [[0, "INFLOW", 50.0]], 1: [[0, "INFLOW", 200.0]],
                    2: [[0, "INFLOW", 20.0]]},
        stage_filerecord="m.sfr.stg", budget_filerecord="m.sfr.bud",
        observations=stream_obs,
    )
    drain = mf.drn(
        stress_period_data={0: [[(0, c), 96.0, 50.0] for c in DRN_CELLS]},
        observations=mf.DrnFlowTargets(locations={"seep": DRN_CELLS}).to_flopy_obs(
            None, csv_name="drn_flow_targets.csv"),
    )
    lake = mf.lak.flopy(
        nlakes=1, packagedata=[[0, 115.0, len(LAKE_CELLS), "pond"]],
        connectiondata=[[0, i, (0, c), "VERTICAL", 0.001, 0.0, 0.0, 0.0, 0.0]
                        for i, c in enumerate(LAKE_CELLS)],
        perioddata={0: [[0, "WITHDRAWAL", 100.0]]}, boundnames=True,
        stage_filerecord="m.lak.stg",
        observations=mf.LakeStageTargets(locations={"pond": 0}).to_flopy_obs(
            csv_name="lake_stage_targets.csv"),
    )
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=110.0, botm=[0.0]),
        mf.npf(k=k, icelltype=0),
        mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=90.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), 102.0] for c in edge]}),
        stream, drain, lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    timing = {"ats": True, "ats_dt0": 1.0, "ats_dtmin": 0.01, "ats_dtmax": 4.0} if ats else {}
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=nper, perioddata=[(perlen, 1 if ats else NSTP, 1.0)] * nper, **timing)])
    with warnings.catch_warnings():
        # A built grid object cannot be frozen into the run's snapshot, and this run
        # is never reopened.
        warnings.filterwarnings("ignore", "Run 'r' built but its frozen snapshot", UserWarning)
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    return run.model()


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    """The four runs in the module docstring, by name."""

    root = tmp_path_factory.mktemp("obs_route_period_end")
    return {
        "steps": _build_and_run(root / "steps"),
        "ats_k10": _build_and_run(root / "ats_k10", ats=True, k=10.0),
        "ats_k2": _build_and_run(root / "ats_k2", ats=True, k=2.0),
        "hourly": _build_and_run(root / "hourly", perlen=HOUR, nper=HOURLY_NPER),
        "twelve": _build_and_run(root / "twelve", perlen=2.0, nper=12),
    }


def _raw(model, csv: str) -> pd.DataFrame:
    """An observation file as MF6 wrote it -- one row per time step -- read to the last digit."""

    return pd.read_csv(Path(model.workspace) / csv, float_precision="round_trip")


def _at_period_ends(raw: pd.DataFrame, perlen: float = PERLEN, nper: int = NPER) -> pd.DataFrame:
    """The rows of ``raw`` at t = perlen, 2 perlen, ... -- found by time, not position."""

    rows = [int(np.flatnonzero(np.isclose(raw["time"], perlen * (per + 1), rtol=1e-6))[0])
            for per in range(nper)]
    return raw.iloc[rows].reset_index(drop=True)


def _families():
    """Each family on the observation route: its class, locations, file and raw columns."""

    last = len(REACHES) - 1
    return {
        "lake": (mf.LakeStageTargets, {"pond": 0}, "lake_stage_targets.csv",
                 {"pond": ["POND"]}),
        "sfr_stage": (mf.SfrStageTargets, {"gage": 1}, "sfr_stage_targets.csv",
                      {"gage": ["GAGE"]}),
        "sfr_flow": (mf.SfrFlowTargets, {"up": 0, "down": last}, "sfr_flow_targets.csv",
                     {"up": ["UP"], "down": ["DOWN"]}),
        "drn": (mf.DrnFlowTargets, {"seep": DRN_CELLS}, "drn_flow_targets.csv",
                {"seep": [f"SEEP_C{c}" for c in DRN_CELLS]}),
    }


# -- one row per period, at its end ---------------------------------------------------


@pytest.mark.slow
def test_every_family_reads_one_row_per_period_at_its_end(runs):
    """9 rows MF6 wrote, 3 periods: each family hands on the period-end rows, keyed by
    the period, where it handed on all nine keyed by t = 3.33 .. 30."""

    model = runs["steps"]
    for family, (make, locations, csv, columns) in _families().items():
        raw = _raw(model, csv)
        assert len(raw) == NPER * NSTP, family
        ends = _at_period_ends(raw)
        series = make(locations=locations).simulated_series(model)
        assert list(series["time"]) == list(range(NPER)), family
        for name, labels in columns.items():
            expected = ends.loc[:, labels].sum(axis=1)
            if len(labels) == 1:  # MF6's own number, to the bit
                np.testing.assert_array_equal(series[name], expected, err_msg=family)
            else:  # a zone: its cells added, in another order than here
                np.testing.assert_allclose(series[name], expected, rtol=1e-14, err_msg=family)
            # Every quantity moves within a period, so which row is read is visible.
            assert not np.isclose(raw.loc[:, labels].sum(axis=1).iloc[NSTP - 2],
                                  expected.iloc[0]), family


@pytest.mark.slow
def test_the_period_end_rows_are_the_file_routes_values(runs):
    """The two routes agree at each period's end: the lake with the stage file's
    period-end record, the reach with the stage file, the flows with the budget files
    (SFR in MF6's NEGATIVE sign on both routes since ledger 350, where the budget route
    reported the magnitude), and the DRN zone cell by cell to the bit."""

    model = runs["steps"]
    lake = mf.LakeStageTargets(locations={"pond": 0}).simulated_series(model)
    stage_file = model.outputs.lak.stage.get()
    assert stage_file.shape == (NPER * NSTP, 1)  # the stage file's record clock
    np.testing.assert_array_equal(lake["pond"], stage_file[NSTP - 1::NSTP, 0])

    gage = mf.SfrStageTargets(locations={"gage": 1}).simulated_series(model)
    np.testing.assert_allclose(
        gage["gage"], mf.SfrStageTargets(locations={"a": 1}).simulated_series(model)["a"],
        rtol=1e-14)
    flow = mf.SfrFlowTargets(locations={"up": 0}).simulated_series(model)
    np.testing.assert_allclose(
        flow["up"], mf.SfrFlowTargets(locations={"a": 0}).simulated_series(model)["a"],
        rtol=1e-10)

    # Each cell to the bit: a one-cell zone ``seep`` of cell c reads MF6's ``seep_c<c>``.
    cbc = CellBudgetFile(Path(model.workspace) / "m.cbc", precision="double")
    for cell in DRN_CELLS:
        one = mf.DrnFlowTargets(locations={"seep": [cell]}).simulated_series(model)["seep"]
        for per in range(NPER):
            drains = cbc.get_data(text="DRN", kstpkper=(NSTP - 1, per))[0]
            assert drains["q"][drains["node"] == cell + 1][0] == one.iloc[per], (cell, per)
    cbc.close()
    seep = mf.DrnFlowTargets(locations={"seep": DRN_CELLS}).simulated_series(model)
    budget = mf.DrnFlowTargets(locations={"zone": DRN_CELLS}).simulated_series(model)
    np.testing.assert_allclose(seep["seep"], budget["zone"], rtol=1e-14)


@pytest.mark.slow
def test_compare_reads_a_period_key_and_refuses_a_time_key_naming_the_period(runs):
    """Keyed by period, every family's ``compare()`` read ``NaN``; keyed by MF6's time
    (10, 20, 30) it read the period ends. Now the period reads, and the time is
    refused, naming the period each one ends."""

    model = runs["steps"]
    for family, (make, locations, _csv, columns) in _families().items():
        name = next(iter(columns))
        locations = {name: locations[name]}
        series = make(locations=locations).simulated_series(model)
        # A flow is entered in MF6's sign, negative (ledger 350); a stage as it reads.
        value = -1.0 if make in (mf.SfrFlowTargets, mf.DrnFlowTargets) else 1.0
        by_period = make(locations=locations,
                         values=pd.DataFrame({"time": range(NPER), name: value})).compare(model)
        simulated = "sim_stage" if "sim_stage" in by_period.columns else "sim_flow"
        np.testing.assert_array_equal(by_period.sort_values("per")[simulated], series[name],
                                      err_msg=family)
        by_time = make(locations=locations,
                       values=pd.DataFrame({"time": [10, 20, 30], name: value}))
        with pytest.raises(ValueError, match=(
                r"\.compare\(\.\.\.\): the observation route is keyed by zero-based stress "
                r"period, one row per period at its end \(ledger 346\), and this model's "
                r"periods are 0 \.\. 2; 3 target row\(s\) name none of them.*10 is the end "
                r"of period 0 in simulated time -- key it 0; 20 is the end of period 1.*"
                r"30 is the end of period 2 in simulated time -- key it 2")):
            by_time.compare(model)


@pytest.mark.slow
@pytest.mark.parametrize("route", ["observation", "file"])
@pytest.mark.parametrize("family", ["lake", "sfr_stage", "sfr_flow", "drn"])
def test_the_residuals_of_a_set_given_no_values_are_one_row_per_period(runs, family, route):
    """A capture template -- a set given no values, a forecast set -- answered
    ``residuals()`` with ``KeyError: 'per'`` for lakes, reach stages and reach flows on
    BOTH routes: ``compare()`` renamed ``per`` to ``time`` (on the observation route into
    a SECOND ``time`` column). DRN zones were fixed in ledger 342; all four now answer
    through one helper, ``_shared._compare_without_values``. The ``file`` route uses a
    name the OBS does not declare, so the stage or budget file is read."""

    model = runs["steps"]
    make, locations, _csv, _columns = _families()[family]
    if route == "file":
        locations = {f"{name}_file": where for name, where in locations.items()}
    targets = make(locations=locations)
    series = targets.simulated_series(model)
    residuals = targets.residuals(model)
    simulated = "sim_stage" if "sim_stage" in residuals.columns else "sim_flow"
    target = "stage_target" if simulated == "sim_stage" else "flow_target"
    assert list(residuals.columns).count("time") == 1
    assert len(residuals) == NPER * len(locations)
    for name in locations:
        rows = residuals[residuals["name"] == name].sort_values("per")
        assert rows["per"].tolist() == rows["time"].tolist() == list(range(NPER)), name
        np.testing.assert_array_equal(rows[simulated], series[name])
    assert residuals[target].isna().all() and residuals["residual"].isna().all()
    assert targets.stats(model)["n"].item() == 0


# -- what PEST is handed ------------------------------------------------------------


def _build(model, workspace: Path, observe):
    cal = model.pest("obsper", workspace=workspace, start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    for targets in observe:
        cal.observe(targets)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pst = cal.build("obsper.pst", noptmax=0)
    return cal, pst


def _read_by_pest(template: Path, table: str, output: Path | None = None) -> pd.Series:
    """What PEST's instruction file reads from ``output`` (default: the table itself)."""

    pyemu = pytest.importorskip("pyemu")
    ins = pyemu.pst_utils.InstructionFile(str(template / f"{table}.ins"))
    return ins.read_output_file(str(output or template / table))["obsval"]


def _observation_at(template: Path, table: str, period: int, column: str) -> str:
    """The observation PEST reads from one cell of a simulated table, found by what it reads."""

    frame = pd.read_csv(template / table)
    probe = frame.astype({column: float})
    probe.loc[frame.iloc[:, 0].astype(float) == period, column] = PROBE
    path = template / f"probe_{table}"
    probe.to_csv(path, index=False)
    read = _read_by_pest(template, table, path)
    hits = read.index[np.isclose(read.to_numpy(dtype=float), PROBE)]
    assert len(hits) == 1, (table, period, column, list(hits))
    return hits[0]


@pytest.mark.slow
def test_the_pest_build_puts_each_period_target_on_its_period_end_row(runs, tmp_path):
    """Every family keyed by period, observation route: 0 of 9 rows carried a target
    and all nine were weighted at the base run's own value (F1, measured); ledger 340
    refused the targets instead. Each lands on its period's row now, the table holds
    one row per period, and nothing else is weighted."""

    pytest.importorskip("pyemu")
    model = runs["steps"]
    values = {"pond": [114.0, 113.0, 112.0], "gage": [100.0, 101.0, 102.0],
              "up": [-150.0, -90.0, -210.0], "seep": [-700.0, -850.0, -860.0]}
    observe = [
        mf.LakeStageTargets(locations={"pond": 0},
                            values=pd.DataFrame({"time": range(NPER), "pond": values["pond"]})),
        mf.SfrStageTargets(locations={"gage": 1},
                           values=pd.DataFrame({"time": range(NPER), "gage": values["gage"]})),
        mf.SfrFlowTargets(locations={"up": 0},
                          values=pd.DataFrame({"time": range(NPER), "up": values["up"]})),
        mf.DrnFlowTargets(locations={"seep": DRN_CELLS},
                          values=pd.DataFrame({"time": range(NPER), "seep": values["seep"]})),
    ]
    cal, pst = _build(model, tmp_path / "pest", observe)
    template = cal.template_workspace
    obs = pst.observation_data
    tables = {"pond": "lak_stage_simulated_lake_stage.csv", "gage": "sfr_stage_simulated_sfr_stage.csv",
              "up": "sfr_flow_simulated_sfr_flow.csv", "seep": "drn_flow_simulated_drn_flow.csv"}
    for column, table in tables.items():
        frame = pd.read_csv(template / table)
        assert list(frame["time"]) == list(range(NPER)), table
        for period, value in enumerate(values[column]):
            name = _observation_at(template, table, period, column)
            assert obs.loc[name, "obsval"] == value, (table, period)
            assert obs.loc[name, "weight"] == 1.0, (table, period)
    assert int((obs["weight"] > 0).sum()) == 4 * NPER
    assert len(obs) == 4 * NPER  # nothing per step left in the control file

    late = mf.DrnFlowTargets(locations={"seep": DRN_CELLS},
                             values=pd.DataFrame({"time": [10, 30], "seep": -700.0}))
    with pytest.raises(ValueError, match=(
            r"DrnFlowTargets \(PEST build, prefix 'drn_flow'\): 2 target\(s\) match nothing.*"
            r"'seep' at time 10 \(period 10\): no row with time = 10.*Its rows are time = 0, "
            r"1, 2.*10 is the end of period 0 in simulated time -- key it 0\. 30 is the end "
            r"of period 2 in simulated time -- key it 2\.")):
        _build(model, tmp_path / "late", [late])


@pytest.mark.slow
def test_an_ats_realization_writes_the_rows_the_base_instruction_file_reads(runs, tmp_path):
    """ATS decides the number of steps at run time, so the rows MF6 writes move with
    the parameters -- and PEST reads each realization by the BASE run's line positions.
    The base run (K = 10) and a realization (K = 2) each write one row per period, and
    the base instruction file reads the realization's period-end values."""

    pytest.importorskip("pyemu")
    base, other = runs["ats_k10"], runs["ats_k2"]
    raw_base, raw_other = (_raw(m, "drn_flow_targets.csv") for m in (base, other))
    assert len(raw_base) != len(raw_other)  # ATS took another number of steps
    seep = mf.DrnFlowTargets(locations={"seep": DRN_CELLS},
                             values=pd.DataFrame({"time": range(NPER), "seep": -700.0}))
    cal, _pst = _build(base, tmp_path / "pest", [seep])
    template = cal.template_workspace
    table = "drn_flow_simulated_drn_flow.csv"
    assert list(pd.read_csv(template / table)["time"]) == list(range(NPER))

    realization = tmp_path / "realization.csv"
    write_named_series_targets(
        "drn_flow", template / "drn_flow_target_locations.csv", realization,
        sim_ws=other.workspace, model_name="m")
    read = _read_by_pest(template, table, realization)
    expected = _at_period_ends(raw_other).loc[:, [f"SEEP_C{c}" for c in DRN_CELLS]].sum(axis=1)
    np.testing.assert_allclose(read.to_numpy(dtype=float), expected, rtol=1e-12)


@pytest.mark.slow
def test_hourly_periods_read_live_find_every_period_end(runs):
    """MF6 ends each hourly period at the 1/24 FloPy WROTE (0.04166667), so its period
    ends run ahead of the lengths the live model holds by 3.3e-9 more every period --
    past ``END_RTOL`` from the fourth. Measured from MF6's own previous period end, the
    time into each period is the written length every period."""

    model = runs["hourly"]
    held = np.cumsum(period_lengths(model))
    raw = _raw(model, "drn_flow_targets.csv")
    ends = raw["time"].to_numpy()[NSTP - 1::NSTP]
    drift = np.abs(ends - held)
    assert (drift[3:] > END_RTOL * np.maximum(held[3:], 1.0)).all()  # the premise
    series = mf.DrnFlowTargets(locations={"seep": DRN_CELLS}).simulated_series(model)
    assert list(series["time"]) == list(range(HOURLY_NPER))
    np.testing.assert_allclose(
        series["seep"], raw.iloc[NSTP - 1::NSTP, 1:].sum(axis=1), rtol=1e-14)


@pytest.mark.slow
def test_mf6s_thirteen_digit_observation_time_stays_inside_end_rtol(runs):
    """MF6 writes a text observation file's VALUES with every digit but its TIME at 13
    significant digits (``0.1388889000000E-01``) -- a rounding ``END_RTOL`` was never
    calibrated for. Rounded as MF6 rounds it, ten years of hourly periods (87,600
    ends, up to t = 3,650) still read every period, each within 4.4e-9 of its length
    (ledger 347)."""

    lines = (Path(runs["hourly"].workspace) / "drn_flow_targets.csv").read_text().splitlines()
    for line in lines[1:]:
        token = line.split(",")[0]
        digits = token.upper().split("E")[0].lstrip("+-").replace(".", "").lstrip("0")
        assert len(digits) == 13, token

    nper = 24 * 3650
    # MF6 ends each period at the last end plus the length it READ, FloPy's 8 decimals.
    ends = np.cumsum(np.full(nper, float(f"{HOUR:15.8f}")))
    times = [float(f"{end:.13G}") for end in ends]
    assert period_end_rows(times, [HOUR] * nper) == {per: per for per in range(nper)}
    assert np.abs(np.diff(np.concatenate([[0.0], times])) - HOUR).max() < END_RTOL


@pytest.mark.slow
def test_a_time_key_below_the_last_period_reads_as_that_period(runs):
    """The limit, pinned (ledger 347): a set keyed by MF6's time -- this route's
    spelling until ledger 346 -- is refused only when a key names no period. With
    twelve 2-day periods, t = 2, 4, 6 END periods 0, 1, 2 and ARE periods 2, 4, 6, and
    they read there in silence, each another flow; a key past the last period (t = 24,
    the end of period 11) refuses the set, naming the period it ends."""

    model = runs["twelve"]
    series = mf.DrnFlowTargets(locations={"seep": DRN_CELLS}).simulated_series(model)
    old = mf.DrnFlowTargets(locations={"seep": DRN_CELLS},
                            values=pd.DataFrame({"time": [2, 4, 6], "seep": -700.0}))
    read = old.compare(model).sort_values("per")["sim_flow"].to_numpy()
    np.testing.assert_array_equal(read, series["seep"].to_numpy()[[2, 4, 6]])
    assert not np.isclose(read, series["seep"].to_numpy()[[0, 1, 2]]).any()

    reaching = mf.DrnFlowTargets(locations={"seep": DRN_CELLS},
                                 values=pd.DataFrame({"time": [2, 4, 24], "seep": -700.0}))
    with pytest.raises(ValueError, match=(
            r"periods are 0 \.\. 11; 1 target row\(s\) name none of them: 'seep' at time 24\. "
            r"24 is the end of period 11 in simulated time -- key it 11")):
        reaching.compare(model)


@pytest.mark.slow
def test_a_zone_declared_over_a_text_and_a_binary_file_reads_as_one(runs, tmp_path):
    """MF6 writes a text observation file's time at 13 significant digits and a BINARY
    one's in full, so a zone declared over one of each -- cell 20 in ``drn_a.csv``,
    cells 21 and 27 in ``drn_b.bsv`` -- was refused as two files holding "different
    times" (3.333333333333 against 3.3333333333333335). They are MF6's one clock: the
    zone reads what the one-file declaration reads, to the bit (ledger 347)."""

    import subprocess

    from myflopy.modflow.mf6.observations._shared import _observation_output_table

    copy = tmp_path / "ws"
    shutil.copytree(runs["steps"].workspace, copy)
    records = [f"  seep_c{cell}  DRN  1 {cell + 1}" for cell in DRN_CELLS]
    (copy / "m.drn.obs").write_text(
        "BEGIN options\nEND options\n\n"
        f"BEGIN continuous  FILEOUT  drn_a.csv\n{records[0]}\nEND continuous\n\n"
        "BEGIN continuous  FILEOUT  drn_b.bsv  BINARY\n"
        + "\n".join(records[1:]) + "\nEND continuous\n")
    (copy / "drn_flow_targets.csv").unlink()
    run = subprocess.run(["mf6"], cwd=copy, capture_output=True, text=True)
    assert run.returncode == 0, run.stdout[-2000:]
    text_time = _observation_output_table(copy / "drn_a.csv", owner="t")[1][0].to_numpy()
    binary_time = _observation_output_table(copy / "drn_b.bsv", owner="t")[1][0].to_numpy()
    assert not np.array_equal(text_time, binary_time)  # the premise: one clock, two roundings
    np.testing.assert_allclose(text_time, binary_time, rtol=1e-12)

    reloaded = mf.load_mf6_run(copy, verbosity_level=0)
    split = mf.DrnFlowTargets(locations={"seep": DRN_CELLS}).simulated_series(reloaded)
    whole = mf.DrnFlowTargets(locations={"seep": DRN_CELLS}).simulated_series(runs["steps"])
    assert list(split["time"]) == list(range(NPER))
    np.testing.assert_array_equal(split["seep"], whole["seep"])

    # Still a check: one step 2e-10 of its time off (past MF6's 13 digits) is refused.
    text = pd.read_csv(copy / "drn_a.csv", dtype=str)
    text.loc[0, "time"] = f"{float(text.loc[0, 'time']) * (1 + 2e-10):.13G}"
    text.to_csv(copy / "drn_a.csv", index=False)
    with pytest.raises(ValueError, match=r"drn_a\.csv and drn_b\.bsv hold different times"):
        mf.DrnFlowTargets(locations={"seep": DRN_CELLS}).simulated_series(reloaded)


# -- what is refused ------------------------------------------------------------------


@pytest.mark.slow
def test_a_file_short_of_a_period_end_is_refused_naming_the_period(runs, tmp_path):
    """A run that stopped, or a copy cut short: the last period's end missing, or the
    whole last period -- each read a value from before the end, or one row short of
    the base run's instruction file."""

    copy = tmp_path / "ws"
    shutil.copytree(runs["steps"].workspace, copy)
    raw = pd.read_csv(copy / "drn_flow_targets.csv")
    reloaded = mf.load_mf6_run(copy, verbosity_level=0)
    zone = mf.DrnFlowTargets(locations={"seep": DRN_CELLS})

    raw.iloc[:-1].to_csv(copy / "drn_flow_targets.csv", index=False)
    with pytest.raises(ValueError, match=(
            r"A stress period's END is missing from the observation file "
            r"drn_flow_targets\.csv \(package 'drn'\): period 2 saved steps 1, 2, the last "
            r"ending 6\.66667 into a period 10 long.*MF6 writes a continuous observation at "
            r"the end of EVERY time step.*Run the model to its end")):
        zone.simulated_series(reloaded)

    raw.iloc[:-NSTP].to_csv(copy / "drn_flow_targets.csv", index=False)
    with pytest.raises(ValueError, match=(
            r"drn_flow_targets\.csv \(package 'drn'\) holds no row in stress period\(s\) 2 "
            r"\(zero-based\) of the 3 its TDIS declares: its last row is t = 20\.0, in "
            r"period 1")):
        zone.simulated_series(reloaded)


def test_the_rule_measures_each_period_from_mf6s_own_previous_end():
    """Six hourly periods as MF6 writes them -- each end the sum of the 8-decimal
    lengths it read. Comparing with the cumulative lengths the model holds refuses
    from the fourth period; differencing reads every one."""

    held = [HOUR] * HOURLY_NPER
    written = np.cumsum([float(f"{HOUR:15.8f}")] * HOURLY_NPER)
    times = []
    for per in range(HOURLY_NPER):
        start = written[per - 1] if per else 0.0
        times += [start + (written[per] - start) * (step + 1) / NSTP for step in range(NSTP)]
    assert period_end_rows(times, held) == {per: per * NSTP + NSTP - 1 for per in range(HOURLY_NPER)}
    assert (np.abs(written - np.cumsum(held))[3:] > END_RTOL).all()

    # ATS: another number of steps in each period; and a zero-length period.
    assert period_end_rows([4.0, 10.0, 11.0, 15.0, 20.0], [10.0, 10.0]) == {0: 1, 1: 4}
    assert period_end_rows([1.0, 1.0, 2.0], [1.0, 0.0, 1.0]) == {0: 0, 1: 1, 2: 2}


def test_a_file_that_is_not_this_simulations_is_refused():
    with pytest.raises(ValueError, match=r"holds 1 row\(s\) past the end.*t = 25\.0"):
        period_end_rows([10.0, 20.0, 25.0], [10.0, 10.0], file="x.csv")
    with pytest.raises(ValueError, match=r"x\.csv's times go back at row 2 \(10\.0 then 5\.0\)"):
        period_end_rows([10.0, 5.0], [10.0, 10.0], file="x.csv")
    with pytest.raises(ValueError, match=r"not a number at row 2"):
        period_end_rows(["10.0", "n/a"], [10.0, 10.0], file="x.csv")
    with pytest.raises(ValueError, match=r"holds no row at all"):
        period_end_rows([], [10.0], file="x.csv")


def test_a_model_with_no_tdis_is_refused_on_the_observation_route(tmp_path):
    """MF6 labels an observation row by time alone, so a double with a declared OBS
    and no TDIS cannot tell which row ends a period."""

    pd.DataFrame({"time": [10.0, 20.0], "POND": [101.0, 102.0]}).to_csv(
        tmp_path / "lake.csv", index=False)
    package = SimpleNamespace(
        package_name="lak",
        obs=[SimpleNamespace(continuous=SimpleNamespace(
            get_data=lambda: {"lake.csv": [("pond", "STAGE", 1, None)]}))],
        packagedata=SimpleNamespace(get_data=lambda: None),
    )
    model = SimpleNamespace(workspace=tmp_path, package=lambda name: package)
    with pytest.raises(ValueError, match=(
            r"LakeStageTargets: cannot tell which row of the observation file lake\.csv "
            r"\(package 'lak'\) ends each stress period.*carries no TDIS")):
        mf.LakeStageTargets(locations={"pond": 0}).simulated_series(model)
