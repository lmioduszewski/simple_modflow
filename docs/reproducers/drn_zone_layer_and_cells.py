"""A DRN zone is read by layer and cell, once per zone: seven ways it was not.

Ledger 342. ``DrnFlowTargets`` summed a zone's drains on the budget-file route by
matching MF6's cell-budget NODE -- the model-wide user node, ``layer * ncpl + cell``
on a DISV -- against the zone's in-layer cell numbers, so a zone below layer 0 read
the LAYER-0 drain that shares a cell number, or nothing; the PEST zone snapshot
dropped ``layer`` besides. A mixed ``GeoDataFrame`` never resolved its geometry
rows; ``to_flopy_obs`` wrote one record per locations ROW; a dict carrying a
``geometry`` key lost it; ``residuals`` on a value-less set raised ``KeyError``;
and a partial miss was said only to a logger nobody hears (ledger 288).

    /home/lukem/python/envs/gw/bin/python docs/reproducers/drn_zone_layer_and_cells.py

It uses only the public API (plus the forward run's own writer, which a PEST
template calls), so it measures whichever tree is on ``PYTHONPATH`` -- the "before"
column below is this script run against the parent commit's ``src``
(``git archive 4fbcc03 src | tar -x -C <dir>``, then ``PYTHONPATH=<dir>/src``);
for H, added by ledger 343's review, "before" is ledger 342's own commit, 6aba608.
Needs ``mf6``.

THE SETUP. The 6 x 6 irregular Voronoi grid of ``tests/test_flow_targets_period_end.py``
(``ncpl`` 36), TWO confined layers (top 110 ft, bottoms 50 and 0 ft), CHD at 102 ft on
the west and east edges of layer 0, heads filling from 100 ft; two 10-day periods of
one step, ``SAVE ... LAST``. Cell 3 of layer 0 is INACTIVE (``idomain`` 0), so every
later node's REDUCED number is one less than its USER number. Drains: the zone
``deep`` at cells 20, 21 and 27 of LAYER 1 (elevation 60 ft, conductance 50), and a
DECOY at cell 20 of LAYER 0 (elevation 99 ft, conductance 5). An MF6 DRN observation
records ``deep``'s three cells, so a target set named ``deep`` reads the observation
route and one named ``basin`` the budget file.

  A. What MF6's cell-budget file holds, read with FloPy: the DRN nodes (1-based) and
     the layer-1 zone's sum; then ``basin`` (budget route) and ``deep`` (observation
     route) through ``simulated_series``.
  B. A PEST build observing ``basin``, and its forward-run writer called once on the
     same workspace: the base table and the forward run's.
  C. A GeoDataFrame of zone ``a`` (explicit cells ``[20]``, layer 1) and zone ``b``
     (NO cells, drawn over cell 27, layer 1): ``b``'s resolved cells.
  D. A geometry zone ``g`` given as TWO rows (cells 20 and 21, layer 1):
     ``to_flopy_obs`` records and ``summary()['n_cells']``; then MF6 run with those
     records as the DRN's observations.
  E. ``{"name": ["g"], "layer": [1], "geometry": [<cell 27>]}`` and the same as a
     list of one record dict.
  F. ``residuals`` of ``basin`` given no values.
  G. Zones ``in`` (cell 27) and ``out`` (5,000 ft beyond the grid); head targets
     ``IN``/``OUT`` the same: Python warnings raised by ``zone_definitions`` and
     ``match_to_model``.
  H. (ledger 343) Zone ``basin`` given cells ``[20, 99]`` in layer 1 -- cell 99 is
     past the 36 cells of a layer -- and cells ``[20]`` in layer 5, on the budget
     route; then the same pair as a DRN observation MF6 is asked to run.
  I. (ledger 343, found, not fixed) A zone with a 36-character name reaching two
     layers -- the decoy cell in layer 0, cells 20 and 21 in layer 1 -- written by
     ``to_flopy_obs`` as three 45-character labels (``<name>_lay<l>_c<c>``): the column
     names MF6 writes, and which route ``simulated_series`` then reads.

MEASURED (mf6 6.8.0.dev0, pyEMU 1.4.0, 2026-09-27; H, I 2026-09-29)  before    after
  A  CBC DRN nodes (1-based)                    [21, 57, 58, 64]    the same: USER nodes
                                                                    (reduced would read 20,
                                                                    56, 57, 63 past the
                                                                    inactive cell 3)
     layer-1 zone, raw CBC (period 0 / 1)       -5975.73 / -6003.62 the same
     basin (budget route)                       -6.77 / -7.79: the  -5975.73 / -6003.62
                                                LAYER-0 decoy
     deep (observation route)                   -5975.73 / -6003.62 the same
  B  base table and forward run                 -6.77 / -7.79 both  -5975.73 / -6003.62
                                                (the decoy); the    both; the snapshot
                                                snapshot's columns  carries layer
                                                name, group,
                                                weight, cells
  C  b's cells                                  []                  [27]
  D  to_flopy_obs records                       4: g_c20, g_c21,    2: g_c20, g_c21
                                                g_c20, g_c21
     summary n_cells                            4                   2
     MF6 given the 4 records                    runs, Normal        (no duplicate to
                                                termination: a      write)
                                                duplicate name is
                                                accepted
  E  dict with a geometry key                   ValueError: polygon [[27]]
                                                zones require a
                                                GeoDataFrame
     list of record dicts                       the same            [[27]]
  F  residuals, no values                       KeyError: "['per']  one row per period,
                                                not in index"       per == time
  G  Python warnings, 1 of 2 DRN zones missed   0                   1, naming 'out'
     Python warnings, 1 of 2 head targets       0                   1, naming 'OUT'
  H  [20, 99], layer 1 (budget route)           -1992.36 / -2002.16 ValueError naming
                                                -- cell 20 ALONE,   'basin' [(1, 99)] and
                                                with no word        "2 layer(s) of 36 cells"
     [20], layer 5 (budget route)               an empty series     the same, [(5, 20)]
     the pair (1, 99) as an MF6 observation     MF6 refuses it on both trees: "Cell2d
                                                number in list (100) is outside of the grid"
  I  the three 45-character labels, as MF6      the same on both trees: all three
     writes them, and the route read            truncated to ONE 40-character name,
                                                `<name>_LAY` (MF6 accepts it, Normal
                                                termination); the reader finds none of
                                                its labels and reads the budget file
                                                (times 0, 1, periods). Since ledger
                                                344 (f0fe9c7 -> its commit) the
                                                observation route reads each DECLARED
                                                record by its column: times 10, 20
"""

from __future__ import annotations

import pathlib
import tempfile
import warnings

import geopandas as gpd
import numpy as np
import pandas as pd
from flopy.utils import CellBudgetFile
from shapely.geometry import Point

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid

ZONE = [20, 21, 27]
DECOY = 20
INACTIVE = 3


def _grid():
    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    return grid, props, [int(c) for c in edge if int(c) != INACTIVE]


def _toy(root, *, observations=None):
    """The model in the module docstring, run; ``observations`` replaces the DRN OBS."""

    grid, props, edge = _grid()
    idomain = np.ones((2, grid.ncpl), dtype=int)
    idomain[0, INACTIVE] = 0
    if observations is None:
        observations = mf.DrnFlowTargets(
            locations=[{"name": "deep", "layer": 1, "cells": ZONE}]
        ).to_flopy_obs(None, csv_name="drn_flow_targets.csv")
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=2, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=110.0,
                botm=[50.0, 0.0], idomain=idomain),
        mf.npf(k=10.0, icelltype=0),
        mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=100.0),
        mf.chd(stress_period_data={0: [[(0, c), 102.0] for c in edge]}),
        mf.drn(stress_period_data={0: [[(1, c), 60.0, 50.0] for c in ZONE]
                                   + [[(0, DECOY), 99.0, 5.0]]},
               observations=observations),
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=2, perioddata=[(10.0, 1, 1.0)] * 2)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = mf.Project(root / "p", name="p").run("r", simulation)
    return run


def _attempt(section, label, action):
    try:
        outcome = action()
    except Exception as error:  # noqa: BLE001 - reports whatever each tree does, crash included
        text = str(error).splitlines()[0] if str(error) else ""
        print(f"{section}  {label}: {type(error).__name__}: {text[:200]}")
        return
    print(f"{section}  {label}: {outcome}")


def _two(values):
    return " / ".join(f"{float(v):.2f}" for v in values)


def section_a(model):
    cbc = CellBudgetFile(pathlib.Path(model.workspace) / "m.cbc", precision="double")
    ncpl = model.vor.ncpl
    sums = []
    for kstpkper in cbc.get_kstpkper():
        record = cbc.get_data(text="DRN", kstpkper=kstpkper)[0]
        nodes = sorted(int(n) for n in record["node"])
        sums.append(record["q"][np.isin(record["node"], [ncpl + c + 1 for c in ZONE])].sum())
    cbc.close()
    print(f"A  CBC DRN nodes (1-based): {nodes}")
    print(f"A  layer-1 zone, raw CBC: {_two(sums)}")
    basin = mf.DrnFlowTargets(locations=[{"name": "basin", "layer": 1, "cells": ZONE}])
    _attempt("A", "basin (budget route)",
             lambda: _two(basin.simulated_series(model)["basin"]))
    deep = mf.DrnFlowTargets(locations=[{"name": "deep", "layer": 1, "cells": ZONE}])
    _attempt("A", "deep (observation route, period ends)",
             lambda: _two(deep.simulated_series(model)["deep"]))


def section_b(model, root):
    from myflopy.modflow.mf6.pest.forward_run import write_named_series_targets

    targets = mf.DrnFlowTargets(
        locations=[{"name": "basin", "layer": 1, "cells": ZONE}],
        values=pd.DataFrame({"time": [0, 1], "basin": [-2000.0, -2000.0]}))

    def build():
        cal = model.pest("zone", workspace=root / "b", start_datetime="2024-01-01")
        cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
        cal.observe(targets, prefix="drn")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cal.build("zone.pst")
        template = pathlib.Path(cal.template_workspace)
        base = pd.read_csv(template / "drn_simulated_drn_flow.csv")["basin"]
        again = root / "forward.csv"
        locations = next(template.glob("drn_target_locations.*"))
        write_named_series_targets(
            [{"kind": "drn_flow", "locations_file": locations,
              "output_csv": again}],
            sim_ws=model.workspace, model_name="m")
        forward = pd.read_csv(again)["basin"]
        snapshot = (pd.read_csv(locations) if locations.suffix == ".csv"
                    else gpd.read_file(locations))
        return (f"base {_two(base)}, forward {_two(forward)}; snapshot columns "
                f"{list(snapshot.columns)}")

    _attempt("B", "base table / forward run", build)


def section_c(model):
    poly = model.vor.gdf_vorPolys.geometry
    zones = gpd.GeoDataFrame(
        {"name": ["a", "b"], "layer": [1, 1], "cells": [[20], np.nan]},
        geometry=[None, poly.loc[27]], crs=model.vor.crs)
    _attempt("C", "b's cells", lambda: mf.DrnFlowTargets(locations=zones)
             .zone_definitions(model).set_index("name").loc["b", "cells"])


def section_d(model, root):
    poly = model.vor.gdf_vorPolys.geometry
    zone = gpd.GeoDataFrame({"name": ["g", "g"], "layer": [1, 1]},
                            geometry=[poly.loc[20], poly.loc[21]], crs=model.vor.crs)
    targets = mf.DrnFlowTargets(locations=zone)
    records = targets.to_flopy_obs(model)["drn_flow_targets.csv"]
    print(f"D  to_flopy_obs records: {len(records)}: {[r[0] for r in records]}")
    _attempt("D", "summary n_cells", lambda: int(targets.summary(model).loc[0, "n_cells"]))
    if len({r[0] for r in records}) == len(records):
        print("D  MF6 with duplicate labels: not written (no duplicate labels)")
        return

    def run():
        result = _toy(root / "d", observations={"drn_flow_targets.csv": records})
        return f"run.success {result.success}; last lines {result.report[-3:]}"

    _attempt("D", "MF6 with duplicate labels", run)


def section_e(model):
    cell27 = model.vor.gdf_vorPolys.geometry.loc[27]
    _attempt("E", "dict with a geometry key", lambda: mf.DrnFlowTargets(
        locations={"name": ["g"], "layer": [1], "geometry": [cell27]})
        .zone_definitions(model)["cells"].tolist())
    _attempt("E", "list of record dicts", lambda: mf.DrnFlowTargets(
        locations=[{"name": "g", "layer": 1, "geometry": cell27}])
        .zone_definitions(model)["cells"].tolist())


def section_f(model):
    basin = mf.DrnFlowTargets(locations=[{"name": "basin", "layer": 1, "cells": ZONE}])
    _attempt("F", "residuals, no values", lambda: basin.residuals(model)
             .loc[:, ["time", "per", "sim_flow"]].to_dict("records"))


def section_g(model):
    polys = model.vor.gdf_vorPolys
    _xmin, _ymin, xmax, ymax = polys.total_bounds
    far = Point(xmax + 5000.0, ymax + 5000.0)
    zones = gpd.GeoDataFrame({"name": ["in", "out"], "layer": [1, 1]},
                             geometry=[polys.geometry.loc[27], far.buffer(10.0)],
                             crs=model.vor.crs)
    wells = gpd.GeoDataFrame({"name": ["IN", "OUT"], "layer": [1, 1]},
                             geometry=[polys.geometry.loc[27].centroid, far],
                             crs=model.vor.crs)
    for label, call in (
        ("1 of 2 DRN zones missed", lambda: mf.DrnFlowTargets(locations=zones)
         .zone_definitions(model)),
        ("1 of 2 head targets missed", lambda: mf.HeadTargets(
            locations=wells, values=pd.DataFrame({"time": [0], "IN": [100.0], "OUT": [100.0]}))
         .match_to_model(model)),
    ):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            call()
        said = [str(w.message)[:90] for w in caught if issubclass(w.category, UserWarning)]
        print(f"G  Python warnings, {label}: {len(said)} {said}")


def section_h(model, root):
    for label, layer, cells in (("[20, 99], layer 1", 1, [20, 99]), ("[20], layer 5", 5, [20])):
        zone = mf.DrnFlowTargets(locations=[{"name": "basin", "layer": layer, "cells": cells}])
        _attempt("H", f"{label} (budget route)",
                 lambda zone=zone: _two(zone.simulated_series(model)["basin"])
                 if len(zone.simulated_series(model)) else "an empty series")

    def run():
        result = _toy(root / "h", observations={"drn_flow_targets.csv": [
            ("basin_c99", "DRN", (1, 99))]})
        errors = [line.strip() for line in result.report if "outside of the grid" in line]
        return f"run.success {result.success}; {errors[:1]}"

    _attempt("H", "the pair (1, 99) as an MF6 observation", run)


def section_i(root):
    name = "a_zone_name_thirty_six_characters_x"
    name += "x" * (36 - len(name))
    zone = mf.DrnFlowTargets(locations=[{"name": name, "layer": 0, "cells": [DECOY]},
                                        {"name": name, "layer": 1, "cells": [20, 21]}])
    records = zone.to_flopy_obs(None)["drn_flow_targets.csv"]

    def run():
        result = _toy(root / "i", observations={"drn_flow_targets.csv": records})
        header = (pathlib.Path(result.model().workspace) / "drn_flow_targets.csv"
                  ).read_text().splitlines()[0].split(",")
        times = zone.simulated_series(result.model())["time"].tolist()
        # Before ledger 346 the times told the routes apart -- 0, 1 the budget file's
        # periods, 10, 20 the observation route's MF6 times; since, both are periods.
        return (f"run.success {result.success}; asked {[len(r[0]) for r in records]} "
                f"characters, MF6 wrote {header[1:]}; the reader's times {times}")

    _attempt("I", "45-character labels", run)


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as scratch:
        root = pathlib.Path(scratch)
        toy = _toy(root / "toy").model()
        section_a(toy)
        section_b(toy, root)
        section_c(toy)
        section_d(toy, root)
        section_e(toy)
        section_f(toy)
        section_g(toy)
        section_h(toy, root)
        section_i(root)
