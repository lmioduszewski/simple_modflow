"""What MF6 does with three SFR inputs ``mf.sfr``/``SFRBuilder`` write.

Ledger 363, ledger 364, ledger 365 and ledger 366.

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/sfr_stream_inputs.py

Run on d71e421's ``src`` it prints the "before": sections A1 and A3 refused by MF6,
section C's stream-keyed rainfall on 2 of 6 reaches. Sections B and D use ``mf.sfr.flopy``
only, so they measure MF6 on any tree (D's "validate" line needs this tree's builder).

THE SETUP. One small model: a 6 x 6 irregular Voronoi grid of ~100 ft cells, one
confined layer (top 110, bottom 0, K 10), CHD at 90 ft on the west and east edges, one
steady period of one day, and two straight streams, beds at 100 ft and streambed K
1e-6 so they barely leak: ``source`` (y = 150, 4 reaches, 100 ft3/d into its
headwater) and ``canal`` (y = 420, 2 reaches), nothing connecting them. The tests are
``tests/test_sfr_stream_inputs.py``.

WHAT IT MEASURES
  A. A ``StreamDiversion`` of FRACTION 0.25 from source's third reach to the canal's
     first: (1) declared alone, (2) beside a ``StreamConnection`` naming the same two
     reaches -- the old workaround -- and (3) with the canal first in the stream table,
     so the receiving reach is reach 0. MF6 diverts only along a reach connection:
     without one it stops reading the package with "Reach 3 is not a upstream reach
     for reach 5 ... diversion 1 from reach 3 to reach 5 is not possible". A downstream
     connection to reach 0 has no negative integer; FloPy writes -0.0 as "-1".
  B. Which SFR ``packagedata`` values MF6 reads as a time-series name: RHK and MAN
     run, RWID, RGRD, RTP and RBTH stop with 'Error converting "kts" to a real number'
     (FloPy's ``gwf-sfr.dfn`` marks only rhk and man ``time_series true``).
  C. Rainfall 0.01 ft/d as ``{0: {"source": r, "canal": r}}`` and as ``r``: the
     RAINFALL MF6 writes per reach to the SFR budget file, in ft3/d and divided by the
     reach's width x length. RAINFALL is a rate per unit of streambed area -- MF6 rains
     rate x width x length on a reach -- so each reach takes the rate unchanged and the
     two spellings must agree.
  D. Networks a diversion makes that MF6 cannot run (ledger 366), written straight
     from the builder's tables -- ``SFRBuilder.validate`` refuses each of them now. MF6's
     ``sfr_check_ustrf`` wants every reach's downstream reaches' USTRF to sum to one,
     diversion receivers included, and a receiving reach's USTRF to be 0: (1) the
     DEFAULT ``source_location="downstream"`` on a stream with no outlet connection --
     its last reach's only downstream reach is the receiver -- stops on "Upstream
     fractions ... not equal to one (0.0)"; (2) the same with the receiver's USTRF set
     to 1, which MF6 tolerates with a warning when the receiver is the source's only
     downstream reach, CRASHES MF6 6.8 in that warning ("End of record",
     ``gwf-sfr.f90``); (3) a diversion into the canal's SECOND reach, whose first reach
     then has no USTRF to flow by; (4) a ``StreamConnection`` leaving ``source`` at its
     third reach (no diversion): that reach then has two downstream reaches of USTRF 1.
"""

from __future__ import annotations

import pathlib
import shutil
import subprocess
import tempfile
import warnings

import flopy
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid

LINES = {
    "source": LineString([(120.0, 150.0), (480.0, 150.0)]),
    "canal": LineString([(320.0, 420.0), (480.0, 420.0)]),
}
DIVERSION = mf.StreamDiversion(
    "source", "canal", amount=0.25, source_location=Point(300.0, 150.0), receiver_location="upstream"
)


def _grid():
    return irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)


def _builder(order=("source", "canal"), **values):
    grid = _grid()
    streams = gpd.GeoDataFrame({"sid": list(order)}, geometry=[LINES[name] for name in order], crs=grid.crs)
    settings = dict(
        connection_mode="explicit",
        connections=(mf.StreamConnection("source", None), mf.StreamConnection("canal", None)),
        reach_top=100.0, streambed_k=1.0e-6, inflow={0: {"source": 100.0}},
    )
    settings.update(values)
    return mf.SFRBuilder(
        context=mf.ModelContext(grid=grid, domain=np.ones((1, grid.ncpl), dtype=int)),
        nper=1, streams=streams, stream_id="sid", **settings,
    )


def _run(root, grid, sfr):
    """Run the setup's model with ``sfr``; (success, workspace, MF6's error lines)."""

    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    edge = np.flatnonzero((xy[:, 0] < 100.0) | (xy[:, 0] > 500.0))
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]), vertices=props["vertices"],
                cell2d=props["cell2d"], top=110.0, botm=[0.0]),
        mf.npf(k=10.0, icelltype=0), mf.ic(strt=90.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), 90.0] for c in edge]}),
        sfr,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = mf.Project(root, name="p").run("r", simulation)
    workspace = pathlib.Path(run.model().workspace)
    listing = (workspace / "mfsim.lst").read_text().splitlines()
    start = next((i for i, line in enumerate(listing) if "ERROR REPORT" in line), None)
    errors = [] if start is None else [line.strip() for line in listing[start + 1 : start + 4] if line.strip()]
    return run.success, workspace, errors


def _observed(builder, kinds):
    observations = [
        (f"{kind.lower().replace('-', '_')}_{reach}", kind, reach + 1)  # MF6 ids are one-based
        for kind in kinds
        for reach in range(builder.total_nreaches)
    ]
    return builder.with_updates(options={"observations": {"sfr.obs.csv": observations}}).build()


def _section_a(root):
    print("A. a FRACTION 0.25 diversion from source to canal")
    both = (
        mf.StreamConnection("source", "canal", source_location=Point(300.0, 150.0), receiver_location="upstream"),
        mf.StreamConnection("canal", None),
    )
    forms = {
        "1 alone": _builder(diversions=(DIVERSION,)),
        "2 with its StreamConnection": _builder(diversions=(DIVERSION,), connections=both),
        "3 into reach 0": _builder(order=("canal", "source"), diversions=(DIVERSION,)),
    }
    for label, builder in forms.items():
        (_, source, receiver, _), = builder.resolved_diversions
        print(f"   {label}: reach {source} -> reach {receiver}; connectiondata {builder.connectiondata}")
        ok, workspace, errors = _run(root / label.split()[0], builder.grid,
                                     _observed(builder, ("UPSTREAM-FLOW", "DOWNSTREAM-FLOW")))
        if not ok:
            print(f"      MF6 refused: {' '.join(errors)}")
            continue
        flows = pd.read_csv(workspace / "sfr.obs.csv").iloc[-1].rename(str.lower)
        after = builder.stream_reaches["source"][builder.stream_reaches["source"].index(source) + 1]
        print(f"      ran: {-flows[f'downstream_flow_{source}']:.3f} ft3/d leaves reach {source}; "
              f"{flows[f'upstream_flow_{receiver}']:.3f} enters canal reach {receiver}, "
              f"{flows[f'upstream_flow_{after}']:.3f} stays in source (reach {after})")


def _section_b(root):
    print("\nB. which SFR packagedata values MF6 reads as a time-series name")
    grid = _grid()
    numbers = dict(rwid=5.0, rgrd=0.001, rtp=100.0, rbth=1.0, rhk=0.001, man=0.03)
    series = {"filename": "sfr.ts", "timeseries": [(0.0, 0.002, 0.04), (1.0, 0.002, 0.04)],
              "time_series_namerecord": [("kts", "mts")], "interpolation_methodrecord": ["linear", "linear"]}
    trials = {"rhk and man": {"rhk": "kts", "man": "mts"}, **{name: {name: "kts"} for name in ("rwid", "rgrd", "rtp", "rbth")}}
    for label, swap in trials.items():
        values = {**numbers, **swap}
        rows = [[reach, (0, cell), 100.0, values["rwid"], values["rgrd"], values["rtp"], values["rbth"],
                 values["rhk"], values["man"], 1 if reach != 1 else 2, 1.0, 0]
                for reach, cell in enumerate((13, 14, 15))]
        sfr = mf.sfr.flopy(packagedata=rows, connectiondata=[[0, -1], [1, 0, -2], [2, 1]],
                           perioddata={0: [[0, "INFLOW", 50.0]]}, timeseries=series)
        ok, _, errors = _run(root / label.replace(" ", "_"), grid, sfr)
        print(f"   {label:>11} as a time series: {'runs' if ok else 'refused: ' + errors[0]}")


def _section_c(root):
    print("\nC. rainfall 0.01 ft/d per stream and as a scalar: MF6's RAINFALL per reach")
    rate = 0.01
    for label, rainfall in (("{stream: r}", {0: {"source": rate, "canal": rate}}), ("r", rate)):
        builder = _builder(rainfall=rainfall)
        rows = [row[0] for row in builder.perioddata[0] if row[1] == "RAINFALL"]
        ok, workspace, errors = _run(root / f"rain_{len(label)}", builder.grid, builder.build())
        assert ok, errors
        budget = flopy.utils.CellBudgetFile(workspace / "m_sfr_budget.sfr", precision="double")
        rain = budget.get_data(text="RAINFALL")[-1]
        budget.close()
        area = 10.0 * builder.reaches["rlen"].to_numpy(float)  # the default width, 10 ft
        print(f"   {label:>11}: rows on reaches {rows}")
        print(f"   {'':>11}  ft3/d {[round(float(q), 3) for q in rain['q']]}")
        print(f"   {'':>11}  / (width x length) {[round(float(q) / float(a), 5) for q, a in zip(rain['q'], area, strict=True)]}")


def _section_d(root):
    print("\nD. networks MF6 cannot run, written without SFRBuilder.validate")
    mid = Point(300.0, 150.0)
    forms = {
        "1 default source_location (last reach)": (_builder(diversions=(mf.StreamDiversion("source", "canal", amount=0.25),)), None),
        "2 the same, receiving USTRF 1": (_builder(diversions=(mf.StreamDiversion("source", "canal", amount=0.25),)), 1.0),
        "3 into the canal's second reach": (_builder(diversions=(mf.StreamDiversion(
            "source", "canal", amount=0.25, source_location=mid, receiver_location="downstream"),)), None),
        "4 a StreamConnection from the third reach": (_builder(connections=(
            mf.StreamConnection("source", "canal", source_location=mid, receiver_location="upstream"),
            mf.StreamConnection("canal", None))), None),
    }
    for label, (builder, receiving_ustrf) in forms.items():
        packagedata = builder.packagedata
        for _, _, receiver, _ in builder.resolved_diversions:
            if receiving_ustrf is not None:
                packagedata[receiver][10] = receiving_ustrf
        try:
            builder.validate()
            refusal = "none"
        except ValueError as error:
            refusal = str(error).split(":")[0]
        sfr = mf.sfr.flopy(packagedata=packagedata, connectiondata=builder.connectiondata,
                           perioddata=builder.perioddata, diversions=builder.diversiondata or None)
        ok, workspace, errors = _run(root / f"d{label.split()[0]}", builder.grid, sfr)
        crash = ""
        if not ok and not errors:  # no ERROR REPORT: rerun MF6 to read how it stopped
            result = subprocess.run(["mf6"], cwd=workspace, capture_output=True, text=True, check=False)
            crash = next((line.strip() for line in (result.stdout + result.stderr).splitlines()
                          if "Fortran runtime error" in line), "no error report")
        print(f"   {label}: {'runs' if ok else 'MF6 refused: ' + (' '.join(errors) or crash)}")
        print(f"      validate: {refusal}")


def main() -> None:
    root = pathlib.Path(tempfile.mkdtemp(prefix="sfr_stream_inputs_"))
    try:
        _section_a(root)
        _section_b(root)
        _section_c(root)
        _section_d(root)
    finally:
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    main()
