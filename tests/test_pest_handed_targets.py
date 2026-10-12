"""What PEST is handed matches the targets (ledger 340).

The PEST build rebuilt each target's observation NAME from the spelling of its key
and skipped any name the control file did not hold, and it left every simulated row
no target named at pyEMU's defaults -- weight 1, the BASE run's own value. So:

* a target keyed ``1.0`` (pyEMU wrote ``per:1``) reached nothing, while
  ``compare()`` paired it; a key of ``1.5`` was truncated to period 1 in
  ``compare()`` and dropped in the build;
* a target at a period the table does not hold vanished without a word;
* two readings of one well in one period: the pst kept the LAST;
* every unmeasured period -- and every period of a forecast set but the forecast
  one -- history-matched the calibration back to the uncalibrated model;
* a blank in a wide table wrote a NaN ``obsval`` and pyEMU refused the file, while
  a typo such as ``"12,5"`` became that same NaN; a set given no values died on
  ``KeyError: 'time'``;
* heads were written by two writers -- the base copy a pivot, whose columns pandas
  sorts, every forward run in mapping order -- so PEST's positional instruction
  file read W1's head as W10's in every realization;
* the lake/SFR/DRN forward run reloaded the simulation's FIRST GWF model, not the
  calibrated one.

Every test here builds a real control file on a real MF6 run of one small model --
the 6 x 6 irregular Voronoi grid of ``test_flow_targets_period_end.py``, one confined
layer, CHD at 102 ft on the west and east edges, heads filling from 90 ft, a
four-reach stream, a three-cell DRN zone and a two-cell lake; three 10-day periods
of three steps, ``SAVE HEAD LAST`` and ``SAVE BUDGET LAST``. It carries no MF6 OBS
package, so every family reads its FILE route, keyed by stress period. Which
observation holds which cell of a simulated table is found the way PEST finds it:
by reading a probe copy of the table through pyEMU's own ``InstructionFile``.
"""

from __future__ import annotations

import re
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid

pyemu = pytest.importorskip("pyemu")

NPER, NSTP, PERLEN = 3, 3, 10.0
REACHES = [13, 14, 15, 16]
DRN_CELLS = [20, 21, 27]
LAKE_CELLS = [8, 9]
#: A value no simulated table holds, written into one cell of a probe copy.
PROBE = 123456.789


def _grid_and_edges():
    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    return grid, props, edge


def _disv(grid, props):
    return mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                   vertices=props["vertices"], cell2d=props["cell2d"], top=110.0, botm=[0.0])


def _run(root: Path, simulation):
    with warnings.catch_warnings():
        # A built grid object cannot be frozen into the run's snapshot, and this run
        # is never reopened.
        warnings.filterwarnings("ignore", "Run 'r' built but its frozen snapshot", UserWarning)
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    return run


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    """The model in the module docstring, run."""

    return run_toy(tmp_path_factory.mktemp("pest_handed"))


def run_toy(root: Path):
    """The model in the module docstring, run under ``root``; its run model.

    A function as well as the fixture, so another module's fixture can run the same
    model (``tests/test_pest_observation_groups.py``).
    """

    grid, props, edge = _grid_and_edges()
    last = len(REACHES) - 1
    stream = mf.sfr.flopy(
        packagedata=[[i, (0, c), 100.0, 5.0, 0.001, 95.0 - 0.1 * i, 1.0, 1.0, 0.03,
                      1 if i in (0, last) else 2, 1.0, 0] for i, c in enumerate(REACHES)],
        connectiondata=[[i, *([i - 1] if i > 0 else []), *([-(i + 1)] if i < last else [])]
                        for i in range(len(REACHES))],
        perioddata={0: [[0, "INFLOW", 50.0]], 1: [[0, "INFLOW", 200.0]],
                    2: [[0, "INFLOW", 20.0]]},
        stage_filerecord="m.sfr.stg", budget_filerecord="m.sfr.bud",
    )
    lake = mf.lak.flopy(
        nlakes=1, packagedata=[[0, 104.0, len(LAKE_CELLS), "pond"]],
        connectiondata=[[0, i, (0, c), "VERTICAL", 0.5, 0.0, 0.0, 0.0, 0.0]
                        for i, c in enumerate(LAKE_CELLS)],
        perioddata={0: [[0, "WITHDRAWAL", 10.0]]}, boundnames=True,
        stage_filerecord="m.lak.stg",
    )
    model = mf.gwf("m", grid=grid, packages=[
        _disv(grid, props),
        mf.npf(k=10.0, icelltype=0),
        mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=90.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), 102.0] for c in edge]}),
        stream,
        mf.drn(stress_period_data={0: [[(0, c), 96.0, 50.0] for c in DRN_CELLS]}),
        lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=NPER, perioddata=[(PERLEN, NSTP, 1.0)] * NPER)])
    return _run(root, simulation).model()


def _build(model, workspace: Path, *, observe=(), forecast=()):
    """A calibration of ``k`` against ``observe``, forecasting ``forecast``, built."""

    cal = model.pest("handed", workspace=workspace, start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    for targets in observe:
        cal.observe(targets)
    for targets in forecast:
        cal.forecast(targets)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pst = cal.build("handed.pst", noptmax=0)
    return cal, pst


def _read_by_pest(template: Path, table: str, output: Path | None = None) -> pd.Series:
    """What PEST's instruction file reads from ``output`` (default: the table itself)."""

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


def _base_value(template: Path, table: str, period: int, column: str) -> float:
    frame = pd.read_csv(template / table)
    return float(frame.loc[frame.iloc[:, 0].astype(float) == period, column].iloc[0])


def _head_targets(names_cells, values, **kwargs):
    locations = pd.DataFrame({"name": [n for n, _ in names_cells],
                              "cell": [c for _, c in names_cells], "layer": 0})
    return mf.HeadTargets(locations=locations, values=values, **kwargs)


# -- (a) and (b): matched by number; refused when nothing holds a target --------


def test_a_period_key_written_as_a_float_reaches_its_row_in_every_family(toy, tmp_path):
    """``1.0`` is period 1: ``compare()`` always read it so, the build dropped it --
    pyEMU wrote ``per:1``, the build looked for ``per:1.0`` and skipped it."""

    heads = _head_targets([("W", 22)], pd.DataFrame({"time": [0.0, 1.0], "W": [50.0, 51.0]}))
    flow = mf.SfrFlowTargets(locations={"a": 0}, values=pd.DataFrame({"time": [1.0], "a": [-7.0]}))
    seep = mf.DrnFlowTargets(locations={"zone": DRN_CELLS},
                             values=pd.DataFrame({"time": [2.0], "zone": [-3.0]}))
    pond = mf.LakeStageTargets(locations={"pond": 0},
                               values=pd.DataFrame({"time": [1.0], "pond": [99.0]}))
    cal, pst = _build(toy, tmp_path / "pest", observe=(heads, flow, seep, pond))
    template = cal.template_workspace
    obs = pst.observation_data

    for table, period, column, value in (
        ("hds_simulated_heads.csv", 0, "W", 50.0),
        ("hds_simulated_heads.csv", 1, "W", 51.0),
        ("sfr_flow_simulated_sfr_flow.csv", 1, "a", -7.0),
        ("drn_flow_simulated_drn_flow.csv", 2, "zone", -3.0),
        ("lak_stage_simulated_lake_stage.csv", 1, "pond", 99.0),
    ):
        name = _observation_at(template, table, period, column)
        assert obs.loc[name, "obsval"] == value, (table, period)
        assert obs.loc[name, "weight"] == 1.0, (table, period)
    # compare() and the build agree on what 1.0 means.
    assert heads.compare(toy).set_index("per")["head_target"].to_dict() == {0: 50.0, 1: 51.0}


def test_a_target_nothing_holds_is_refused_by_the_build_naming_it(toy, tmp_path):
    """A period the table does not hold, and a name no target defines: both vanished
    from the pst without a word, leaving the base run's value weighted. (A lake the
    model does not have is refused before the build reaches the table now, ledger
    503 -- ``tests/test_target_hygiene.py``.)"""

    late = _head_targets([("W", 22)], pd.DataFrame({"time": [5], "W": [50.0]}))
    with pytest.raises(ValueError, match=(
        r"HeadTargets \(PEST build, prefix 'hds'\): 1 target\(s\) match nothing in the "
        r"table PEST reads, hds_simulated_heads\.csv.*'W' at time 5 \(period 5\): no row "
        r"with per = 5.*Its rows are per = 0, 1, 2")):
        _build(toy, tmp_path / "late", observe=(late,))

    ghost = mf.LakeStageTargets(locations={"pond": 0},
                                values=pd.DataFrame({"time": [0], "ghost": [99.0]}))
    with pytest.raises(ValueError, match=r"'ghost' at time 0 \(period 0\): no column 'ghost'"):
        _build(toy, tmp_path / "ghost", observe=(ghost,))


def test_a_fractional_period_is_refused_by_compare_and_by_the_build(toy, tmp_path):
    """``1.5`` was truncated to period 1 by ``compare()`` and dropped by the build."""

    heads = _head_targets([("W", 22)], pd.DataFrame({"time": [1.5], "W": [50.0]}))
    flow = mf.SfrFlowTargets(locations={"a": 0}, values=pd.DataFrame({"time": [10.9], "a": [-7.0]}))
    refusal = r"reads `time` as a zero-based stress period, and {} is not a whole number"
    with pytest.raises(ValueError, match=refusal.format(r"1\.5")):
        heads.compare(toy)
    with pytest.raises(ValueError, match=refusal.format(r"10\.9")):
        flow.compare(toy)
    with pytest.raises(ValueError, match=refusal.format(r"1\.5")):
        _build(toy, tmp_path / "pest", observe=(heads,))


def test_two_measurements_of_one_target_in_one_period_are_refused(toy, tmp_path):
    """The pst kept the LAST of two readings while ``compare()`` counted both."""

    heads = mf.HeadTargets(
        locations=pd.DataFrame({"name": ["A"], "cell": [22], "layer": [0]}),
        values=pd.DataFrame({"time": [0, 0, 1], "name": ["A", "A", "A"],
                             "head": [99.0, 95.0, 97.0]}))
    flow = mf.SfrFlowTargets(
        locations={"r": 3},
        values=pd.DataFrame({"time": [0, 0], "name": ["r", "r"], "flow": [-100.0, -140.0]}))
    refusal = r"two or more measurements of one target in one stress period -- {} at period 0: 2 values"
    with pytest.raises(ValueError, match=refusal.format("'A'")):
        heads.compare(toy)
    with pytest.raises(ValueError, match=refusal.format("'r'")):
        flow.compare(toy)
    with pytest.raises(ValueError, match=refusal.format("'A'")):
        _build(toy, tmp_path / "pest", observe=(heads,))
    # A blank is not a measurement, so it collides with nothing.
    blank = mf.HeadTargets(
        locations=pd.DataFrame({"name": ["A"], "cell": [22], "layer": [0]}),
        values=pd.DataFrame({"time": [0, 0], "name": ["A", "A"], "head": [99.0, None]}))
    assert blank.compare(toy)["head_target"].notna().sum() == 1


def _wide_family(family, values):
    """One target ``w`` of ``family``, given the long ``time``/``name``/``value`` rows."""

    frame = pd.DataFrame(values, columns=["time", "name", "value"])
    sign = -1.0 if family in ("sfr_flow", "drn") else 1.0
    frame["value"] = frame["value"] * sign
    column = {"head": "head", "conc": "conc", "lake": "stage", "sfr_stage": "stage",
              "sfr_flow": "flow", "drn": "flow"}[family]
    frame = frame.rename(columns={"value": column})
    make = {
        "head": lambda: mf.HeadTargets(locations={"w": 22}, values=frame),
        "conc": lambda: mf.ConcTargets(locations={"w": 22}, values=frame),
        "lake": lambda: mf.LakeStageTargets(locations={"w": 0}, values=frame),
        "sfr_stage": lambda: mf.SfrStageTargets(locations={"w": 1}, values=frame),
        "sfr_flow": lambda: mf.SfrFlowTargets(locations={"w": 3}, values=frame),
        "drn": lambda: mf.DrnFlowTargets(locations={"w": [20]}, values=frame),
    }[family]
    return make(), sign


@pytest.mark.parametrize("family", ["head", "conc", "lake", "sfr_stage", "sfr_flow", "drn"])
def test_to_wide_refuses_two_measurements_of_one_target_in_one_period_naming_both(family):
    """``to_wide()`` pivoted with ``aggfunc="first"`` and kept the first of two readings
    in silence, while ``compare()`` and the PEST build refuse the same set (ledger 340).
    It reads the same check now (ledger 371), whose message names each row."""

    twice, sign = _wide_family(family, [(0, "w", 10.0), (1, "w", 11.0), (1, "w", 12.0)])
    first, second = _shown_values(11.0 * sign, 12.0 * sign)
    with pytest.raises(ValueError, match=(
            rf"\.to_wide\(\): two or more measurements of one target in one stress period "
            rf"-- 'w' at period 1: 2 values \({first} at time 1, {second} at time 1\)")):
        twice.to_wide()
    # A blank row is not a measurement: it collides with nothing, and the value reads.
    blank, sign = _wide_family(family, [(0, "w", 10.0), (1, "w", None), (1, "w", 12.0)])
    assert blank.to_wide()["w"].tolist() == [10.0 * sign, 12.0 * sign]
    # A key that is no stress period is refused here as compare() refuses it.
    fraction, _ = _wide_family(family, [(0.5, "w", 10.0)])
    with pytest.raises(ValueError, match=r"\.to_wide\(\) reads `time` as a zero-based stress "
                                         r"period, and 0\.5 is not a whole number"):
        fraction.to_wide()


def _shown_values(*values):
    return [re.escape(repr(float(value))) for value in values]


def test_to_wide_keeps_one_column_for_a_two_screen_well():
    """``get()`` repeats each value row for each screen of a well; that is one value,
    not two measurements, so the check reads the values as given and passes."""

    well = mf.HeadTargets(
        locations=pd.DataFrame({"name": ["w", "w"], "cell": [22, 22], "layer": [0, 1]}),
        values=pd.DataFrame({"time": [0, 1], "w": [10.0, 11.0]}))
    wide = well.to_wide()
    assert wide.columns.tolist() == ["time", "w"]
    assert wide["w"].tolist() == [10.0, 11.0]


_FAMILY_COLUMN = {"head": "head", "conc": "conc", "lake": "stage", "sfr_stage": "stage",
                  "sfr_flow": "flow", "drn": "flow"}


def _two_target_family(family, values, names=("w", "v")):
    """Targets ``names`` of ``family``, given the long ``time``/``name``/``value`` rows."""

    frame = pd.DataFrame(values, columns=["time", "name", _FAMILY_COLUMN[family]])
    ids = {"head": [22, 23], "conc": [22, 23], "lake": [0, 0], "sfr_stage": [1, 2],
           "sfr_flow": [3, 4], "drn": [[20], [21]]}[family]
    locations = dict(zip(names, ids, strict=False))
    return {
        "head": mf.HeadTargets, "conc": mf.ConcTargets, "lake": mf.LakeStageTargets,
        "sfr_stage": mf.SfrStageTargets, "sfr_flow": mf.SfrFlowTargets,
        "drn": mf.DrnFlowTargets,
    }[family](locations=locations, values=frame)


@pytest.mark.parametrize("family", ["head", "conc", "lake", "sfr_stage", "sfr_flow", "drn"])
def test_to_wide_has_a_row_per_period_and_a_column_per_target(family):
    """``to_wide()`` pivoted on the ``time`` label as given, so ``"1"`` and ``1`` were
    two rows of period 1, each blank where the other held the value; and pandas'
    pivot dropped a target whose every row is blank, which ``compare()`` keeps. It
    is keyed by the period now, ``int``, and keeps every target (ledger 422)."""

    sign = -1.0 if family in ("sfr_flow", "drn") else 1.0
    targets = _two_target_family(family, [
        (0, "w", 10.0 * sign), ("1", "w", None), (1, "w", 11.0 * sign),
        (1.0, "v", None), (2, "v", None)])
    wide = targets.to_wide()
    assert wide.columns.tolist() == ["time", "v", "w"]
    assert wide["time"].dtype == np.int64 and wide["time"].tolist() == [0, 1, 2]
    assert wide["v"].isna().all()
    assert wide["w"].tolist()[:2] == [10.0 * sign, 11.0 * sign] and np.isnan(wide["w"].iloc[2])


@pytest.mark.parametrize("family", ["head", "conc", "lake", "sfr_stage", "sfr_flow", "drn"])
@pytest.mark.parametrize("name", ["per", " PER ", "time"])
def test_a_target_named_like_a_key_column_is_refused_naming_it(family, name):
    """A target named ``per`` collided with the period column of the simulated table
    PEST reads (``cannot insert per, already exists``; the build blamed the columns'
    order), and one named ``time`` with the values' own (ledger 423)."""

    with pytest.raises(ValueError, match=(
            rf"{re.escape(repr(name))} cannot name a target: 'per' and 'time' \(in any "
            r"case\) are the key columns")):
        _two_target_family(family, [(0, name, -1.0)], names=(name, "v"))


_OWN_COLUMNS = {"head": ["head_target"], "conc": ["head_target"],
                "lake": ["stage", "sim_stage"], "sfr_stage": ["stage_target", "sim_stage"],
                "sfr_flow": ["flow_target", "sim_flow"], "drn": ["flow_target", "sim_flow"]}


@pytest.mark.parametrize("family, name", [
    (family, name) for family, names in _OWN_COLUMNS.items() for name in names])
def test_a_target_named_like_its_family_s_value_column_is_refused_naming_it(family, name):
    """A target named like the column a wide table is melted into (``stage`` on a lake
    set) or ``compare()`` melts the simulated one into (``sim_stage``) raised pandas'
    ``value_name (stage) cannot match an element in the DataFrame columns`` (review of
    ledger 423)."""

    with pytest.raises(ValueError, match=(
            rf"{re.escape(repr(name.upper()))} cannot name a target: .*; and .*"
            rf"{re.escape(repr(name))}.* are the value columns")):
        _two_target_family(family, [(0, name.upper(), -1.0)], names=(name.upper(), "v"))


@pytest.mark.parametrize("family", ["head", "conc", "lake", "sfr_stage", "sfr_flow", "drn"])
def test_to_wide_of_a_set_given_no_values_is_typed_as_a_full_one(family):
    """A set given no values has ``time`` and one empty column per target, typed and
    sorted as the full table is (they were ``object`` columns in the order given;
    review of ledger 422)."""

    frame = pd.DataFrame(columns=["time", "name", _FAMILY_COLUMN[family]])
    ids = {"head": [22, 23], "conc": [22, 23], "lake": [0, 0], "sfr_stage": [1, 2],
           "sfr_flow": [3, 4], "drn": [[20], [21]]}[family]
    cls = {"head": mf.HeadTargets, "conc": mf.ConcTargets, "lake": mf.LakeStageTargets,
           "sfr_stage": mf.SfrStageTargets, "sfr_flow": mf.SfrFlowTargets,
           "drn": mf.DrnFlowTargets}[family]
    wide = cls(locations=dict(zip(("w", "v"), ids, strict=True)), values=frame).to_wide()
    assert wide.columns.tolist() == ["time", "v", "w"] and wide.empty
    assert wide["time"].dtype == np.int64
    assert all(dtype == np.float64 for dtype in wide.loc[:, ["v", "w"]].dtypes)


# -- (d): nothing is history-matched that nobody measured -----------------------


def test_every_row_no_target_names_carries_weight_zero_forecasts_included(toy, tmp_path):
    """Measured at period 1 only; forecast at period 2 only. Every other row kept
    pyEMU's weight 1 and the BASE run's value -- anchoring the calibration to the
    uncalibrated model -- and the forecast's other periods were history-matched."""

    heads = _head_targets([("W", 22)], pd.DataFrame({"time": [1], "W": [51.0]}))
    predicted = _head_targets([("P", 23)], pd.DataFrame({"time": [2], "P": [0.0]}))
    cal, pst = _build(toy, tmp_path / "pest", observe=(heads,), forecast=(predicted,))
    template = cal.template_workspace
    obs = pst.observation_data

    measured = _observation_at(template, "hds_simulated_heads.csv", 1, "W")
    assert (obs.loc[measured, "obsval"], obs.loc[measured, "weight"]) == (51.0, 1.0)
    for period in (0, 2):
        name = _observation_at(template, "hds_simulated_heads.csv", period, "W")
        assert obs.loc[name, "weight"] == 0.0
        assert obs.loc[name, "obsval"] == pytest.approx(
            _base_value(template, "hds_simulated_heads.csv", period, "W"), rel=1e-12)

    forecast = _observation_at(template, "fore1_simulated_heads.csv", 2, "P")
    assert pst.pestpp_options["forecasts"] == forecast
    for period in (0, 1):
        name = _observation_at(template, "fore1_simulated_heads.csv", period, "P")
        assert obs.loc[name, "weight"] == 0.0
    assert pst.nnz_obs == 1


# -- (e): a blank is unmeasured; a typo is refused; no values still registers ---


def test_a_blank_is_unmeasured_and_a_value_that_is_not_a_number_is_refused(toy, tmp_path):
    """A wide table with a gap -- real well data -- could not be built: its blank
    became a NaN ``obsval`` and pyEMU refused the control file. And a typo became
    the same NaN, so a measurement could vanish unnoticed."""

    gappy = _head_targets(
        [("W1", 22), ("W2", 23)],
        pd.DataFrame({"time": [0, 1], "W1": [50.0, None], "W2": [52.0, 53.0]}))
    cal, pst = _build(toy, tmp_path / "pest", observe=(gappy,))
    obs = pst.observation_data
    assert obs["obsval"].notna().all()
    gap = _observation_at(cal.template_workspace, "hds_simulated_heads.csv", 1, "W1")
    assert obs.loc[gap, "weight"] == 0.0
    assert pst.nnz_obs == 3

    with pytest.raises(ValueError, match=(
            r"HeadTargets values must be finite numbers, or blank where a target was not "
            r"measured, and 1 are neither: 'W1' at time 1 holds '12,5' \(column 'W1'\)")):
        _head_targets([("W1", 22)], pd.DataFrame({"time": [0, 1], "W1": ["50", "12,5"]}))
    with pytest.raises(ValueError, match=r"'r' at time 0 holds 'dry' \(column 'flow'\)"):
        mf.SfrFlowTargets(locations={"r": 3}, values=pd.DataFrame(
            {"time": [0], "name": ["r"], "flow": ["dry"]}))


def test_a_set_given_no_values_registers_and_a_forecast_of_one_names_every_row(
    toy, tmp_path
):
    """A capture template or a forecast placeholder died on ``KeyError: 'time'``."""

    captured = mf.DrnFlowTargets(locations={"zone": DRN_CELLS})
    outlet = mf.SfrFlowTargets(locations={"out": len(REACHES) - 1})
    cal, pst = _build(toy, tmp_path / "pest", observe=(captured,), forecast=(outlet,))
    template = cal.template_workspace
    obs = pst.observation_data

    zone = [_observation_at(template, "drn_flow_simulated_drn_flow.csv", p, "zone")
            for p in range(NPER)]
    out = [_observation_at(template, "fore1_simulated_sfr_flow.csv", p, "out")
           for p in range(NPER)]
    assert (obs.loc[zone + out, "weight"] == 0.0).all()
    assert pst.pestpp_options["forecasts"].split(",") == out


# -- (f): one writer, one order -------------------------------------------------


@pytest.mark.slow
def test_pest_reads_each_well_its_own_head_whatever_order_the_names_are_in(toy, tmp_path):
    """Names out of sorted order -- W2, W1, W10. The base table was a pivot, whose
    columns pandas sorts, and the forward run wrote mapping order; PEST's instruction
    file reads by position, so every realization read W1's head as W10's."""

    wells = _head_targets([("W2", 1), ("W1", 2), ("W10", 34)],
                          pd.DataFrame({"time": [0], "W2": [1.0], "W1": [2.0], "W10": [3.0]}))
    cal, _pst = _build(toy, tmp_path / "pest", observe=(wells,))
    template = cal.template_workspace
    table = "hds_simulated_heads.csv"
    assert (template / table).read_text().splitlines()[0] == "per,W2,W1,W10"

    subprocess.run([sys.executable, "forward_run.py"], cwd=template, check=True,
                   capture_output=True, text=True)
    written = pd.read_csv(template / table)
    assert list(written.columns) == ["per", "W2", "W1", "W10"]
    read = _read_by_pest(template, table)
    assert len(read) == 3 * NPER
    for name, value in read.items():
        label, per = re.search(r"usecol:(.+)_per:(\d+)$", name).groups()
        assert value == pytest.approx(written.loc[int(per), label.upper()], rel=1e-12), name


def test_a_supplied_simulated_table_in_another_order_is_refused(toy, tmp_path):
    """The one writer cannot vouch for a table the caller supplies, so its header is
    checked against what every forward run writes."""

    wells = _head_targets([("W2", 1), ("W1", 2)],
                          pd.DataFrame({"time": [0], "W2": [1.0], "W1": [2.0]}))
    sorted_table = pd.DataFrame({"per": range(NPER), "W1": 1.0, "W2": 2.0})
    spec = mf.HeadTargetObservationSpec(targets=wells, simulated_values=sorted_table)
    with pytest.raises(ValueError, match=(
            r"has columns \['per', 'W1', 'W2'\], and every forward run writes "
            r"\['per', 'W2', 'W1'\]")):
        _build(toy, tmp_path / "pest", observe=(spec,))


def test_two_sets_under_one_prefix_are_refused(toy, tmp_path):
    """The prefix names every file a set writes: the second set overwrote the
    first's table, instruction file and map, and PEST++ refused the control file."""

    first = _head_targets([("W1", 22)], pd.DataFrame({"time": [0], "W1": [1.0]}))
    second = _head_targets([("W2", 23)], pd.DataFrame({"time": [0], "W2": [1.0]}))
    with pytest.raises(ValueError, match="Two observation or forecast sets share the prefix 'hds'"):
        _build(toy, tmp_path / "pest", observe=(first, second))


# -- (g): the forward run reloads the calibrated model --------------------------


@pytest.mark.slow
def test_the_named_series_forward_run_reads_the_calibrated_model_not_the_first(tmp_path):
    """Two GWF models, calibrating the SECOND. The forward run reloaded the first
    one: its DRN zone's flow in every realization, the second's in the base run."""

    grid, props, edge = _grid_and_edges()

    def flow_model(name, cond):
        return mf.gwf(name, grid=grid, packages=[
            _disv(grid, props),
            mf.npf(k=10.0, icelltype=0),
            mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}),
            mf.ic(strt=100.0),
            mf.chd(stress_period_data={0: [[(0, int(c)), 102.0] for c in edge]}),
            mf.drn(stress_period_data={0: [[(0, c), 96.0, cond] for c in DRN_CELLS]}),
            mf.oc(head_filerecord=f"{name}.hds", budget_filerecord=f"{name}.cbc",
                  saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")]),
        ])

    simulation = mf.simulation(
        flow_model("m", 1.0), flow_model("n", 50.0), name="s",
        packages=[mf.tdis(nper=2, perioddata=[(PERLEN, 1, 1.0)] * 2), mf.ims()])
    run = _run(tmp_path, simulation)
    first, second = run.model("m"), run.model("n")
    targets = mf.DrnFlowTargets(locations={"seep": DRN_CELLS},
                                values=pd.DataFrame({"time": [0, 1], "seep": [-500.0, -500.0]}))
    wanted = targets.simulated_series(second)["seep"].to_numpy()
    other = targets.simulated_series(first)["seep"].to_numpy()
    assert not np.allclose(wanted, other)

    cal = second.pest("two", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    cal.parameterize("drn.cond", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 1000.0))
    cal.observe(targets)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build("two.pst", noptmax=0)
    template = cal.template_workspace
    table = template / "drn_flow_simulated_drn_flow.csv"
    base = table.read_text()
    np.testing.assert_allclose(pd.read_csv(table)["seep"], wanted, rtol=1e-12)

    subprocess.run([sys.executable, "forward_run.py"], cwd=template, check=True,
                   capture_output=True, text=True)
    assert table.read_text() == base


# -- the review of ledger 340 (ledger 341) --------------------------------------


def test_a_head_target_on_no_active_cell_is_refused_by_the_build_naming_it(toy, tmp_path):
    """``compare()`` leaves a target on no active cell out with a warning -- one file may
    serve several grids (ledger 211) -- and that warning goes to a logger nobody reads.
    The build joined the values to the matched locations and dropped the rest, so the
    measurement vanished from the control file in silence: the one silent drop 340
    left. A DRN zone that lands nowhere was already refused (no column)."""

    import geopandas as gpd
    from shapely.geometry import Point

    xmin, ymin, xmax, ymax = toy.vor.gdf_vorPolys.total_bounds
    inside = toy.vor.gdf_vorPolys.geometry.iloc[22].centroid
    wells = gpd.GeoDataFrame(
        {"name": ["IN", "OUT"], "layer": [0, 0]},
        geometry=[inside, Point(xmax + 5000.0, ymax + 5000.0)], crs=toy.vor.crs)
    targets = mf.HeadTargets(
        locations=wells,
        values=pd.DataFrame({"time": [0, 1], "IN": [50.0, 51.0], "OUT": [60.0, 61.0]}))

    # compare() is unchanged: it answers, the location left out holding no head.
    compared = targets.compare(toy).set_index("name")["sim_head"]
    assert compared.loc["IN"].notna().all() and compared.loc["OUT"].isna().all()
    with pytest.raises(ValueError, match=(
            r"HeadTargets \(PEST build, prefix 'hds'\): 1 target\(s\) name a location "
            r"that lands on no active cell of this model, so PEST would never see "
            r"them: 'OUT' at time\(s\) 0, 1")):
        _build(toy, tmp_path / "pest", observe=(targets,))


def test_a_name_pandas_would_read_as_a_number_or_a_blank_reaches_pest_verbatim(
    toy, tmp_path
):
    """Wells named ``007``, ``0012`` and ``NA``: the forward-run writer read the target
    map back with ``pd.read_csv``, which reads ``007`` as the integer 7 and ``NA`` as a
    missing value, so it wrote columns ``7``, ``12`` and ``nan`` -- and the base table,
    written by that same writer since 340, failed the build's own header check,
    naming columns the user never wrote. A lake's forward run read its snapshot the
    same way and wrote ``7`` under the base run's ``007``."""

    from myflopy.modflow.mf6.pest.forward_run import (
        _write_head_target_csv,
        write_named_series_targets,
    )

    wells = _head_targets(
        [("007", 22), ("0012", 23), ("NA", 24)],
        pd.DataFrame({"time": [0], "007": [50.0], "0012": [51.0], "NA": [52.0]}))
    pond = mf.LakeStageTargets(locations={"007": 0},
                               values=pd.DataFrame({"time": [0], "007": [99.0]}))
    cal, pst = _build(toy, tmp_path / "pest", observe=(wells, pond))
    template = cal.template_workspace
    heads = template / "hds_simulated_heads.csv"
    stage = template / "lak_stage_simulated_lake_stage.csv"
    assert heads.read_text().splitlines()[0] == "per,007,0012,NA"
    assert stage.read_text().splitlines()[0] == "time,007"
    assert pst.nnz_obs == 4

    # Each forward-run writer, as forward_run.py calls it, writes the same header.
    rewritten = tmp_path / "heads_again.csv"
    _write_head_target_csv("m", template / "hds_head_target_map.csv", rewritten,
                           sim_ws=toy.workspace)
    assert rewritten.read_text() == heads.read_text()
    lake_again = tmp_path / "stage_again.csv"
    write_named_series_targets("lake_stage", template / "lak_stage_target_locations.csv",
                               lake_again, sim_ws=toy.workspace, model_name="m")
    assert lake_again.read_text().splitlines()[0] == "time,007"


def test_the_ies_review_reads_a_name_pandas_would_read_as_a_number_or_a_blank_as_text(
    toy, tmp_path
):
    """The forward-run writers read the build's snapshots as TEXT since ledger 341, and
    the IES review read the same files with a plain ``read_csv``: a well named ``007``
    joined to no cell on the residual map (its map row said ``7``), ``NA`` neither (it
    said ``nan``), a DRN zone named ``007`` was dropped from the map outright, and the
    review panel titled ``NA`` as ``nan``. All of them read through
    ``forward_run.read_target_table`` now. The review is driven over the REAL
    template, with one realization standing in for PESTPP-IES's posterior."""

    import json
    from types import SimpleNamespace

    from myflopy.modflow.mf6.pest.ies import IesResults

    wells = _head_targets(
        [("007", 22), ("0012", 23), ("NA", 24)],
        pd.DataFrame({"time": [0], "007": [50.0], "0012": [51.0], "NA": [52.0]}))
    zone = mf.DrnFlowTargets(locations={"007": DRN_CELLS},
                             values=pd.DataFrame({"time": [0], "007": [-5.0]}))
    cal, pst = _build(toy, tmp_path / "pest", observe=(wells, zone))
    template = cal.template_workspace
    posterior = pd.DataFrame([pst.observation_data["obsval"] + 1.0], index=["base"])

    class _Review(IesResults):
        def __init__(self):
            self.workspace = template
            self.model = toy
            self.pst = pst
            self.__dict__["_metadata"] = json.loads(
                (template / "myflopy_pest_metadata.json").read_text())

        @property
        def posterior(self):
            return SimpleNamespace(_df=posterior)

    review = _Review()
    residuals = review.obs_residuals().set_index(["prefix", "location"])
    assert residuals.loc[("hds", "007"), "cell"] == 22
    assert residuals.loc[("hds", "0012"), "cell"] == 23
    assert residuals.loc[("hds", "na"), "cell"] == 24
    assert residuals.loc[("drn_flow", "007"), "cells"] == DRN_CELLS
    titles = [panel["title"] for panel in review._location_panels(None, 20)]
    assert "hds: NA [hds]" in titles and "hds: 007 [hds]" in titles


def test_an_infinite_value_is_refused_naming_it(toy, tmp_path):
    """``inf`` is a float, so it passed the "is it a number" test and reached the
    control file as an ``obsval``; pyEMU refused the file as holding a NaN and wrote
    ``pst.obs_data.nans.csv`` into the caller's working directory."""

    with pytest.raises(ValueError, match=(
            r"HeadTargets values must be finite numbers, or blank where a target was "
            r"not measured, and 1 are neither: 'W' at time 0 holds inf \(column 'W'\)")):
        _head_targets([("W", 22)], pd.DataFrame({"time": [0], "W": [np.inf]}))
    with pytest.raises(ValueError, match=r"'r' at time 1 holds '-inf' \(column 'flow'\)"):
        mf.SfrFlowTargets(locations={"r": 3}, values=pd.DataFrame(
            {"time": [1], "name": ["r"], "flow": ["-inf"]}))
