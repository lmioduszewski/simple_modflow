"""The review of ledger 340: three ways a target still failed to reach PEST as written.

Ledger 341. Ledger 340 refused every target the PEST build could not place -- except a
head or concentration target whose LOCATION lands on no active cell, which the build's
inner join still dropped in silence; it wrote the base head table with the forward
run's own writer, which reads the target map back with ``pd.read_csv`` and so turns a
well named ``007`` into ``7`` and one named ``NA`` into a missing value; and it let an
infinite value through as a number. This script measures each, and one thing 340's
docs claimed that is not so on the observation route.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/pest_handed_targets_review.py

It uses only the public API, so it measures whichever tree is on ``PYTHONPATH`` -- the
"before" column below is this script run against ledger 340's commit's ``src``
(``git archive 58988d7 src | tar -x -C <dir>``, then ``PYTHONPATH=<dir>/src``).
Needs ``mf6``.

THE SETUP. ``pest_handed_targets.py``'s model (imported from it): a 6 x 6 irregular
Voronoi grid, one confined layer, three 10-day periods of three steps, ``SAVE ...
LAST``, no MF6 OBS, so every family reads its FILE route. Section D's model is the
same grid with a DRN zone and an MF6 DRN observation, run as three ONE-day periods of
one step, so its DRN targets read the OBSERVATION route, whose rows are MF6's simulated
times 1, 2, 3.

  A. Head targets IN (cell 22's centroid) and OUT (5,000 ft beyond the grid's corner),
     each with values at periods 0 and 1.
  B. Head targets named ``007``, ``0012`` and ``NA``, built alone; then a lake target
     named ``007``, built alone, and its forward-run writer called once.
  C. A head value of ``inf``.
  D. DRN targets keyed by PERIOD 1 and 2 on the observation route.

MEASURED (mf6 6.8.0.dev0, pyEMU 1.4.0, 2026-09-27)          before           after
  A  build with OUT in the targets              built, nnz_obs 2 of 3     REFUSED, naming
                                                (IN's column only):       'OUT' at times 0, 1
                                                OUT's 2 values dropped,
                                                no Python warning
  B  build with heads 007, 0012, NA             REFUSED by 340's header   built, nnz_obs 3;
                                                check, naming columns     header per,007,0012,NA
                                                ['per', '7.0', '12.0',
                                                'nan'] as the base's
     lake 007: forward-run header (base         time,7                    time,007
     time,007)
  C  a head of inf                              constructed; the build    REFUSED when the
                                                died in pyEMU ("NaNs in   targets are built
                                                obs_data") and wrote
                                                pst.obs_data.nans.csv
                                                into the caller's cwd
  D  period keys 1, 2 on the observation route  built, weight 1, reading  unchanged here; since
                                                t = 1 and t = 2: the ENDS ledger 346 the rows are
                                                of periods 0 and 1        periods 0, 1, 2, and 1
                                                                          and 2 read periods 1, 2

"before" is ledger 340's own commit. Its parent (``b46366b``) built the three head
names -- header ``per,0012,007,NA``, the sorted pivot 340 replaced -- because PEST read
the forward run's ``7`` by position under the base table's ``007``. D is the
interim ledger 340 described as "refused": a period key matches nothing on that route
only while no simulated time equals it.
"""

from __future__ import annotations

import pathlib
import sys
import tempfile
import warnings

import numpy as np
import pandas as pd

import myflopy as mf

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from pest_handed_targets import (  # noqa: E402
    DRN_CELLS,
    _build,
    _grid,
    _heads,
    _run,
    _toy,
)


def _attempt(section, label, action):
    try:
        outcome = action()
    except Exception as error:  # noqa: BLE001 - reports whatever each tree does, crash included
        text = str(error).splitlines()[0]
        print(f"{section}  {label}: {type(error).__name__}: {text[:220]}")
        return
    print(f"{section}  {label}: {outcome}")


def section_a(model, root):
    import geopandas as gpd
    from shapely.geometry import Point

    _xmin, _ymin, xmax, ymax = model.vor.gdf_vorPolys.total_bounds
    wells = gpd.GeoDataFrame(
        {"name": ["IN", "OUT"], "layer": [0, 0]},
        geometry=[model.vor.gdf_vorPolys.geometry.iloc[22].centroid,
                  Point(xmax + 5000.0, ymax + 5000.0)],
        crs=model.vor.crs)
    targets = mf.HeadTargets(locations=wells, values=pd.DataFrame(
        {"time": [0, 1], "IN": [50.0, 51.0], "OUT": [60.0, 61.0]}))

    def build():
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            _cal, pst = _build(model, root / "a", observe=(targets,))
        said = [w for w in caught if "OUT" in str(w.message)]
        return (f"built, nnz_obs {pst.nnz_obs} of {pst.nobs}; Python warnings naming "
                f"OUT: {len(said)}")

    _attempt("A", "build with OUT in the targets", build)


def section_b(model, root):
    wells = _heads([("007", 22), ("0012", 23), ("NA", 24)], pd.DataFrame(
        {"time": [0], "007": [50.0], "0012": [51.0], "NA": [52.0]}))
    pond = mf.LakeStageTargets(locations={"007": 0},
                               values=pd.DataFrame({"time": [0], "007": [99.0]}))

    def heads():
        cal, pst = _build(model, root / "b1", observe=(wells,))
        header = (pathlib.Path(cal.template_workspace) / "hds_simulated_heads.csv")
        return f"built, nnz_obs {pst.nnz_obs}; header {header.read_text().splitlines()[0]}"

    def lake():
        from myflopy.modflow.mf6.pest.forward_run import write_named_series_targets

        cal, _pst = _build(model, root / "b2", observe=(pond,))
        template = pathlib.Path(cal.template_workspace)
        again = root / "stage_again.csv"
        # Found by suffix: the lake set's default prefix is ``lak_stage`` since ledger
        # 372 and was ``stage`` before, and this runs against both trees.
        (locations,) = template.glob("*_target_locations.csv")
        (simulated,) = template.glob("*_simulated_lake_stage.csv")
        write_named_series_targets(
            [{"kind": "lake_stage", "locations_file": locations,
              "output_csv": again}],
            sim_ws=model.workspace, model_name="m")
        base = simulated.read_text().splitlines()[0]
        return f"base header {base}, forward-run header {again.read_text().splitlines()[0]}"

    _attempt("B", "build with heads 007, 0012, NA", heads)
    _attempt("B", "lake 007", lake)


def section_c(model, root):
    def build():
        targets = _heads([("W", 22)], pd.DataFrame({"time": [0], "W": [np.inf]}))
        _build(model, root / "c", observe=(targets,))
        return "built"

    _attempt("C", "a head of inf", build)


def section_d(root):
    grid, disv, edge = _grid()
    zone = mf.DrnFlowTargets(locations={"seep": DRN_CELLS})
    model = mf.gwf("m", grid=grid, packages=[
        disv, mf.npf(k=10.0, icelltype=0),
        mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}), mf.ic(strt=90.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), 102.0] for c in edge]}),
        mf.drn(stress_period_data={0: [[(0, c), 96.0, 50.0] for c in DRN_CELLS]},
               observations=zone.to_flopy_obs(None, csv_name="drn_flow_targets.csv")),
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")])])
    run = _run(root / "d", mf.simulation(model, name="s", packages=[
        mf.tdis(nper=3, perioddata=[(1.0, 1, 1.0)] * 3)]))
    daily = run.model()
    targets = mf.DrnFlowTargets(locations={"seep": DRN_CELLS}, values=pd.DataFrame(
        {"time": [1, 2], "seep": [-500.0, -500.0]}))
    rows = list(targets.simulated_series(daily)["time"])

    def build():
        _cal, pst = _build(daily, root / "d" / "pest", observe=(targets,))
        obs = pst.observation_data
        weighted = obs.loc[obs["weight"] > 0].index.str.extract(r"(time:[\d.]+)$")[0].tolist()
        # Before ledger 346 the rows were MF6's times (1, 2, 3), so keys 1 and 2 read the
        # ENDS of periods 0 and 1; since, they are the periods (0, 1, 2).
        return f"rows at time = {rows}; built, weight > 0 on {weighted}"

    _attempt("D", "period keys 1, 2 on the observation route", build)


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as scratch:
        root = pathlib.Path(scratch)
        toy = _toy(root / "toy")
        section_a(toy, root)
        section_b(toy, root)
        section_c(toy, root)
        section_d(root)
