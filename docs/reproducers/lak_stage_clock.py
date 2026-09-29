"""A lake's stage file read on its own clock, and a period read at its last saved step (ledger 348).

``build_lak_stage_result_table`` labelled ``per`` by the INDEX of the saved stage
record (ledger 301's stage clock), and every view that reduced a period to one value
combined the steps the period saved -- adding a flow's, averaging a stage's, or
painting its first. This script measures each on real MF6 runs:

    /home/lukem/python/envs/gw/bin/python docs/reproducers/lak_stage_clock.py

It uses only public spellings present before and after the fix, and reads the
expected numbers straight from the binary files with FloPy, so it measures whichever
tree is on ``PYTHONPATH`` -- the "before" is ``12b3a57``'s ``src`` (``git archive
12b3a57 src | tar -x -C <dir>``, then ``PYTHONPATH=<dir>/src``). Needs ``mf6``.

THE SETUP. The toy of ``tests/test_lak_stage_clock.py``: a 6 x 6 irregular Voronoi
grid of ~100 ft cells, one confined layer, CHD at 50 ft west and east, ten-day
periods.

  steps    one lake of three vertical connections and a 100 ft3/d withdrawal,
           three periods of THREE steps, SAVE HEAD ALL + SAVE BUDGET ALL
  uneven   the same, periods of 1, 3 and 2 steps
  asym     the same at nstp=3, SAVE HEAD LAST + SAVE BUDGET ALL
  stream   a four-reach stream at nstp=3, everything saved (heads filling from
           90 ft toward 102, so the stage moves within a period)

  1. The LAK stage table's labels, and the stage-file target route: series, compare()
     keyed by period 0/1/2, and a PEST build (rows, obsval, weight).
  2. The lake's stage views at per=1: map, plot(), stage_change.
  3. The lake's exchange views at per=1: q.map (value and hover), q.plot(),
     q.budget, the combined surface-water table; and CHD's q map and series.
  4. The same under SAVE HEAD LAST + SAVE BUDGET ALL: the two files' clocks.
  5. SFR: stage map / plot, the long profile's exchange, q map; the group SFR table.
  6. The group LAK table and diff, and region_heads(per=1).
  7. The review of ledger 348 (ledger 349): the TABLES keyed by period -- a cell-budget
     noun's wide(), the LAK/SFR budget namespaces' summary()/wide() -- and the SFR
     field profiles. Its "before" is a196f2e's src (348 itself).

MEASURED (mf6 6.8.0.dev0, 2026-09-29), before (12b3a57) -> after. Period ends of the
lake (steps): 101.228206, 100.366551, 100.098692 (records 2, 5, 8 of 9).
  1  table per labels                0..8 -> 0, 1, 2 (and a kstpkper column)
     simulated_series                nine rows, the records -> three, the period ends
     compare() keyed 0/1/2           102.702536, 101.823609, 101.228206 (steps 1-3 of
                                     period 0) -> the period ends
     uneven 1/3/2 series             six rows, records 0..5 -> 101.634036, 100.492709,
                                     100.153053 (records 0, 3, 5)
     PEST build                      nine rows: the targets on time:0..2, the base run's
                                     stage on time:3..8 at weight 0 -> three rows, the
                                     targets
  2  stage.map(per=1)                101.823609 (record 1, step 2 of period 0) -> 100.366551
     stage.plot()                    nine points -> the three period ends
     stage_change                    eight transitions between steps -> (0, 1), (1, 2)
                                     between period ends
  3  q.map(per=1) at cell 1          z -0.08300292, hover q -2793.32 and area 33653.24
                                     (three steps ADDED), stage 101.823609 (record 1)
                                     -> -0.05236442, -587.41, 11217.75, 100.366551
     q.plot()                        -25992.61, -7873.47, -2240.85 (sums of steps) ->
                                     -5547.85, -1655.72, -445.79 (the period ends)
     q.budget.get(per=1)             -7873.47 over 9 records -> -1655.72 over 3
     surface-water table, per=1      -0.08300292 -> -0.05236442
     chd q.map(per=1), q.plot()      -2152.678 (-1018.687 - 681.279 - 452.712) ->
                                     -452.712; period sums -> period ends
  4  asym table (kstpkper, per)      no kstpkper -> (2, 0), (2, 1), (2, 2): the stage
                                     file's three records against the budget's nine
     q.map(per=1) hover q            -2793.32 -> -587.41 (the stage was right: record 1
                                     of three IS period 1 there)
  5  sfr stage.map(per=1), reach 0   102.026218 (step 1) -> 102.212593 (step 3)
     sfr stage.plot(), reach 0       100.033536, 102.131971, 101.536691 (step means) ->
                                     101.143706, 102.212593, 101.540962
     profile.get(per=1) q_gwf        364.231 (128.274 + 119.799 + 116.158) -> 116.158
     sfr q.map(per=1)                1.2141 -> 1.16158
     group sfr get(per=1, reach=3)   101.294695 (reach 1 at (0, 1): "per" was the
                                     reach, "reach" the record) -> 101.830363
     group sfr get() per / reach     per 0..3, reach 0..8 -> per 0..2, reach 0..3
  6  group lak get(per=1)            101.823609 / 101.099786 (each member's record 1,
                                     two different instants) -> 100.366551 / 100.492709
     diff summary n                  6 (records paired by index across nstp=3 and
                                     1/3/2) -> 3
     region_heads(per=1)             9 rows, steps (0..2, 1) -> 3 rows, (2, 1)
  7  (before is a196f2e, 348)
     chd q.wide(per=1) at cell 0     -2152.678 (three steps added, beside a map(per=1)
                                     of -452.712 since 348) -> -452.712
     lak.budget.gwf.summary() per 1  q -7873.47, flow_area 94857.75, 9 records,
                                     q/area -0.08300292 -> -1655.72, 31619.25, 3,
                                     -0.05236442 (q.budget.get(per=1) since 348)
     lak.budget.wide(term="GWF")     -7873.47 -> -1655.72
     sfr.budget.gwf.summary(per=1)   reach 0 -364.31 (-128.305 - 119.818 - 116.188)
                                     -> -116.188
     sfr stage.profile.get(per=1)    12 rows, reach 0 at three steps 102.026218,
                                     102.157103, 102.212593 -> 4 rows, 102.212593
     sfr q.profile.get(per=1)        12 rows, reach 0 128.274, 119.799, 116.158 ->
                                     4 rows, 116.158

Why each is so, and what was decided: ledger 348 (sections 1-6) and 349 (section 7).
"""

from __future__ import annotations

import pathlib
import shutil
import sys
import tempfile
import warnings

import numpy as np
import pandas as pd

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid

NPER = 3
REACH_CELLS = [13, 14, 15, 16]
SAVE_ALL = [("HEAD", "ALL"), ("BUDGET", "ALL")]


def _grid():
    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    disv = mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                   vertices=props["vertices"], cell2d=props["cell2d"], top=100.0, botm=[0.0])
    return grid, disv, edge


def _run(root, simulation):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success
    return run.model()


def _lake(root, nstp, saverecord):
    grid, disv, edge = _grid()
    interior = [c for c in range(grid.ncpl) if c not in set(edge.tolist())]
    lake = mf.lak.flopy(
        nlakes=1, packagedata=[[0, 104.0, 3, "pond"]],
        connectiondata=[[0, i, (0, c), "VERTICAL", 0.5, 0.0, 0.0, 0.0, 0.0]
                        for i, c in enumerate(interior[:3])],
        perioddata={0: [[0, "WITHDRAWAL", 100.0]]}, boundnames=True,
        stage_filerecord="m.lak.stg", budget_filerecord="m.lak.bud")
    model = mf.gwf("m", grid=grid, packages=[
        disv, mf.npf(k=10.0, icelltype=0),
        mf.sto(ss=1e-5, sy=0.2, iconvert=0, transient={0: True}), mf.ic(strt=50.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), 50.0] for c in edge]}), lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc", saverecord=saverecord)])
    return _run(root, mf.simulation(model, name="s", packages=[
        mf.tdis(nper=len(nstp), perioddata=[(10.0, n, 1.0) for n in nstp])]))


def _stream(root):
    grid, _, edge = _grid()
    props = grid.get_disv_gridprops()
    last = len(REACH_CELLS) - 1
    stream = mf.sfr.flopy(
        packagedata=[[i, (0, c), 100.0, 5.0, 0.001, 95.0 - 0.1 * i, 1.0, 1.0, 0.03,
                      1 if i in (0, last) else 2, 1.0, 0] for i, c in enumerate(REACH_CELLS)],
        connectiondata=[[i, *([i - 1] if i > 0 else []), *([-(i + 1)] if i < last else [])]
                        for i in range(len(REACH_CELLS))],
        perioddata={0: [[0, "INFLOW", 50.0]], 1: [[0, "INFLOW", 200.0]], 2: [[0, "INFLOW", 20.0]]},
        stage_filerecord="m.sfr.stg", budget_filerecord="m.sfr.bud")
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=110.0, botm=[0.0]),
        mf.npf(k=10.0, icelltype=0), mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=90.0), mf.chd(stress_period_data={0: [[(0, int(c)), 102.0] for c in edge]}),
        stream, mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc", saverecord=SAVE_ALL)])
    return _run(root, mf.simulation(model, name="s", packages=[
        mf.tdis(nper=NPER, perioddata=[(10.0, 3, 1.0)] * NPER)]))


def _file(reader):
    """``{(kstp, kper): values}`` of a stage file, from its own record headers."""

    clock = [(int(k) - 1, int(p) - 1) for k, p in
             dict.fromkeys(zip(reader.recordarray["kstp"], reader.recordarray["kper"],
                               strict=True))]
    data = np.asarray(reader.get_alldata(), dtype=float).reshape(len(clock), -1)
    return dict(zip(clock, data, strict=True))


def _step(frame, key):
    return frame[frame["kstpkper"].map(lambda v: tuple(map(int, v))) == key]


def _show(label, fn):
    try:
        value = fn()
    except Exception as exc:  # noqa: BLE001 -- a measurement, either outcome is printed
        value = f"RAISES {type(exc).__name__}: {str(exc)[:160]}"
    print(f"   {label:<44} {value}")


def _r(values, digits=6):
    return [round(float(v), digits) for v in values]


def main() -> None:
    root = pathlib.Path(tempfile.mkdtemp(prefix="lak_stage_clock_"))
    try:
        _measure(root)
    finally:
        shutil.rmtree(root, ignore_errors=True)


def _measure(root: pathlib.Path) -> None:
    warnings.simplefilter("ignore")
    steps = _lake(root / "steps", [3, 3, 3], SAVE_ALL)
    uneven = _lake(root / "uneven", [1, 3, 2], SAVE_ALL)
    asym = _lake(root / "asym", [3, 3, 3], [("HEAD", "LAST"), ("BUDGET", "ALL")])
    stream = _stream(root / "stream")

    at = {k: v[0] for k, v in _file(steps.fl.lak.output.stage()).items()}
    ends = [at[(2, p)] for p in range(NPER)]
    print("1. the LAK stage table and the stage-file target route (steps)")
    print(f"   stage file: {len(at)} records; period ends (2, p) {_r(ends)}")
    table = steps.lak.results.stage.get()
    _show("table per labels", lambda: sorted(set(table["per"])))
    _show("table has kstpkper", lambda: "kstpkper" in table.columns)
    series = mf.LakeStageTargets(locations={"a": 0}).simulated_series(steps)
    _show("simulated_series time", lambda: list(series["time"]))
    _show("simulated_series values", lambda: _r(series["a"]))
    targets = mf.LakeStageTargets(
        locations={"a": 0}, values=pd.DataFrame({"time": [0, 1, 2], "a": [111.0, 112.0, 113.0]}))
    _show("compare() keyed 0/1/2", lambda: _r(targets.compare(steps).sort_values("per")["sim_stage"]))
    ua = {k: v[0] for k, v in _file(uneven.fl.lak.output.stage()).items()}
    _show("uneven 1/3/2: series", lambda: _r(
        mf.LakeStageTargets(locations={"a": 0}).simulated_series(uneven)["a"]))
    print(f"   {'uneven: period ends (records 0, 3, 5)':<44} "
          f"{_r([ua[(0, 0)], ua[(2, 1)], ua[(1, 2)]])}")

    def _pest():
        cal = steps.pest("lk", workspace=root / "pest", start_datetime="2024-01-01")
        cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
        cal.observe(targets)
        obs = cal.build("lk.pst", noptmax=0).observation_data
        lake = obs[obs.obsnme.str.contains("usecol:a")]
        return [(name.split("usecol:a_")[-1], round(float(v), 4), float(w))
                for name, v, w in zip(lake.obsnme, lake.obsval, lake.weight, strict=True)]
    _show("PEST build (row, obsval, weight)", _pest)

    print("\n2. the lake's stage views at per=1 (steps)")
    cells = sorted(set(table["cell"]))
    print(f"   {'the stage at (2, 1)':<44} {round(float(at[(2, 1)]), 6)}")
    _show("stage.map(per=1) at the lake's cells", lambda: _r(
        np.asarray(steps.lak.results.stage.map(per=1).fig.data[0].z, dtype=float)[cells]))
    _show("stage.plot() x / y", lambda: (list(map(int, steps.lak.results.stage.plot().data[0].x)),
                                         _r(steps.lak.results.stage.plot().data[0].y)))
    _show("stage_change (per0, per1, stage0, stage1)", lambda: [
        (int(a), int(b), round(float(c), 6), round(float(d), 6)) for a, b, c, d in
        steps.lak.results.stage_change.get()[["per0", "per1", "stage0", "stage1"]].itertuples(
            index=False)])

    print("\n3. the lake's exchange views at per=1 (steps)")
    q = steps.lak.results.q.get()
    cell = int(q["cell"].iloc[0])
    by_step = {k: _step(q, (k, 1)) for k in range(3)}
    one = [by_step[k][by_step[k]["cell"] == cell] for k in range(3)]
    print(f"   {f'q_lake at cell {cell}, steps 1/2/3 of period 1':<44} "
          f"{_r([f['q_lake'].iloc[0] for f in one], 2)}, area {round(float(one[2]['flow_area'].iloc[0]), 2)}")

    def _qmap():
        picture = steps.lak.results.q.map(per=1)
        hover = picture._custom_hover
        return (round(float(np.asarray(picture.fig.data[0].z)[cell]), 8),
                round(float(hover["q"][cell]), 2), round(float(hover["flow_area"][cell]), 2),
                round(float(hover["stage"][cell]), 6))
    _show("q.map(per=1): z, hover q, area, stage", _qmap)
    _show("q.plot() y", lambda: _r(steps.lak.results.q.plot().data[0].y, 2))
    _show("lake total q_lake at (2, p)", lambda: _r(
        [_step(q, (2, p))["q_lake"].sum() for p in range(NPER)], 2))
    _show("q.budget.get(per=1): q_lake, record_count", lambda: (
        round(float(steps.lak.results.q.budget.get(per=1)["q_lake"].sum()), 2),
        int(steps.lak.results.q.budget.get(per=1)["record_count"].sum())))
    from myflopy.modflow.mf6.package_budget import build_surface_water_exchange_cell_table
    _show("surface-water table per=1 at the cell", lambda: round(float(
        build_surface_water_exchange_cell_table(steps, per=1, include="lak")
        .set_index("cell")["exchange_intensity"][cell]), 8))
    chd = steps.chd.results.q.get()
    edge_cell = int(chd["cell"].iloc[0])
    _show(f"chd q at cell {edge_cell}, steps 1/2/3 of per 1", lambda: _r(
        [_step(chd, (k, 1)).set_index("cell")["q_gwf"][edge_cell] for k in range(3)], 3))
    _show("chd q.map(per=1) at that cell", lambda: round(float(
        np.asarray(steps.chd.results.q.map(per=1).fig.data[0].z)[edge_cell]), 3))
    _show("chd q.plot() y", lambda: _r(steps.chd.results.q.plot().data[0].y, 2))
    _show("chd total q_gwf at (2, p)", lambda: _r(
        [_step(chd, (2, p))["q_gwf"].sum() for p in range(NPER)], 2))

    print("\n4. SAVE HEAD LAST + SAVE BUDGET ALL (asym)")
    aa = _file(asym.fl.lak.output.stage())
    budget_steps = len(dict.fromkeys(zip(asym.fl.lak.output.budget().recordarray["kstp"],
                                         asym.fl.lak.output.budget().recordarray["kper"],
                                         strict=True)))
    print(f"   {'stage file clock / LAK budget steps':<44} {list(aa)} / {budget_steps}")
    _show("table (kstpkper, per)", lambda: sorted(set(zip(
        asym.lak.results.stage.get().get("kstpkper", pd.Series([None] * 9)).map(str),
        asym.lak.results.stage.get()["per"], strict=True))))

    def _asym_hover():
        hover = asym.lak.results.q.map(per=1)._custom_hover
        return round(float(hover["stage"][cell]), 6), round(float(hover["q"][cell]), 2)
    _show("q.map(per=1) hover stage, q", _asym_hover)
    aq = asym.lak.results.q.get()
    print(f"   {'stage at (2, 1); q_lake at (2, 1)':<44} {round(float(aa[(2, 1)][0]), 6)}; "
          f"{round(float(_step(aq, (2, 1)).set_index('cell')['q_lake'][cell]), 2)}")

    print("\n5. SFR (stream)")
    sa = _file(stream.fl.sfr.output.stage())
    print(f"   {'reach 0 at steps 1/2/3 of period 1':<44} {_r([sa[(k, 1)][0] for k in range(3)])}")
    print(f"   {'reach 0 period ends':<44} {_r([sa[(2, p)][0] for p in range(NPER)])}")
    _show("stage.map(per=1) at reach 0", lambda: round(float(
        np.asarray(stream.sfr.results.stage.map(per=1).fig.data[0].z)[REACH_CELLS[0]]), 6))
    _show("stage.plot() reach 0", lambda: _r(
        {t.name: t.y for t in stream.sfr.results.stage.plot().data}["Reach 0"]))
    sq = stream.sfr.results.q.get()
    _show("reach 0 q_gwf, steps 1/2/3 of period 1", lambda: _r(
        [_step(sq, (k, 1)).set_index("reach")["q_gwf"][0] for k in range(3)], 3))
    _show("profile.get(per=1) reach 0 q_gwf", lambda: round(float(
        stream.sfr.results.profile.get(per=1).set_index("reach")["q_gwf"][0]), 3))
    _show("q.map(per=1) at reach 0 (q / rlen)", lambda: round(float(
        np.asarray(stream.sfr.results.q.map(per=1).fig.data[0].z)[REACH_CELLS[0]]), 5))
    sfr_group = mf.ModelGroup({"a": stream, "b": stream}, reference="a").packages.sfr.results.stage
    _show("group get(per=1, reach=3) stage", lambda: _r(
        sfr_group.get(model_name="a", per=1, reach=3)["stage"]))
    print(f"   {'reach 3 at (2, 1); reach 1 at (0, 1)':<44} "
          f"{round(float(sa[(2, 1)][3]), 6)}; {round(float(sa[(0, 1)][1]), 6)}")
    _show("group get() per / reach values", lambda: (
        sorted(set(sfr_group.get(model_name="a")["per"])),
        sorted(set(sfr_group.get(model_name="a")["reach"]))))

    print("\n6. the group LAK tier and region_heads")
    group = mf.ModelGroup({"steps": steps, "uneven": uneven}, reference="steps")
    _show("group lak stage get(per=1)", lambda: [
        (m, round(float(s), 6)) for m, s in
        group.packages.lak.results.stage.get(per=1)[["model", "stage"]].itertuples(index=False)])
    print(f"   {'steps (2, 1); uneven (2, 1)':<44} {round(float(at[(2, 1)]), 6)}; "
          f"{round(float(ua[(2, 1)]), 6)}")
    _show("diff summary n, argmax_per", lambda: tuple(
        group.diff().packages.lak.results.stage.summary().iloc[0][["n", "argmax_per"]]))
    steps.add_region("west", cellids=[0, 1, 2], layer=0, overwrite=True)
    _show("region_heads('west', per=1) rows, steps", lambda: (
        len(steps.region_heads("west", per=1)),
        sorted(set(map(str, steps.region_heads("west", per=1)["kstpkper"])))))

    print("\n7. the tables keyed by period, and the SFR field profiles (the review, 349)")
    _show(f"chd q.wide(per=1) at cell {edge_cell}", lambda: round(float(
        steps.chd.results.q.wide(per=1).set_index("cell")["per_1"][edge_cell]), 3))
    at_end = _step(q, (2, 1))
    print(f"   {'lake q_lake, flow_area summed at (2, 1)':<44} "
          f"{round(float(at_end['q_lake'].sum()), 2)}, {round(float(at_end['flow_area'].sum()), 2)}")

    def _lak_summary():
        row = steps.lak.budget.gwf.summary().set_index("per").loc[1]
        return (round(float(row["q"]), 2), round(float(row["flow_area"]), 2),
                int(row["record_count"]), round(float(row["q_per_area"]), 8))
    _show("lak.budget.gwf.summary() per 1: q, area, n, q/a", _lak_summary)
    _show("lak.budget.wide(term='GWF') per 1", lambda: round(float(
        steps.lak.budget.wide(term="GWF").set_index("per")["GWF"][1]), 2))
    sb = stream.sfr.budget.get(term="GWF")
    print(f"   {'sfr budget GWF reach 0, steps 1/2/3 of per 1':<44} "
          f"{_r([_step(sb, (k, 1)).set_index('reach')['q'][0] for k in range(3)], 3)}")
    _show("sfr.budget.gwf.summary(per=1) reach 0", lambda: round(float(
        stream.sfr.budget.gwf.summary(per=1).set_index("reach")["q"][0]), 3))
    _show("sfr stage.profile.get(per=1) rows, reach 0", lambda: (
        len(stream.sfr.results.stage.profile.get(per=1)),
        _r(stream.sfr.results.stage.profile.get(per=1).query("reach == 0")["stage"])))
    _show("sfr q.profile.get(per=1) rows, reach 0", lambda: (
        len(stream.sfr.results.q.profile.get(per=1)),
        _r(stream.sfr.results.q.profile.get(per=1).query("reach == 0")["q_gwf"], 3)))


if __name__ == "__main__":
    sys.exit(main())
