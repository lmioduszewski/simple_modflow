"""What the lake-family docstrings say MF6 and myflopy do, measured.

Ledger 318. The measurements behind the ``LakeConnection``, ``LakeOutlet``,
``LakeTable``, ``LakeTableBuilder`` and ``mf.lak`` docstrings, and behind the list of
things ``mf.lak`` and ``LakeTableBuilder`` accepted without checking -- each now
refused or read correctly, so sections 2 and 6-9 print the refusal next to what MF6
did with the same input, which is the reason for it. Run it::

    python docs/reproducers/lak_docstring_claims.py [workspace]

It needs ``mf6`` on PATH (or ``MF6=/path/to/mf6``) and writes nothing outside the
workspace it is given (a temporary directory by default). About 10 s.

**The MF6 bench.** DIS, one layer, 3 x 3 cells of 10 x 10 ft, top 8.0, bottom -10.0,
``k = 1``, a CHD at 5.0 in one corner, and one lake connected to the centre cell
(100 ft2). Sections 1-5 and 9 run it; ``aquifer`` moves the start head and the CHD.

**The myflopy bench.** A 4 x 4 DISV grid of unit cells, two layers (tops 12, 8, 0),
and a 2 x 2 lake in the middle, for what ``mf.lak`` WRITES (sections 6-7, 9), plus a
2 x 2 raster of 1 ft pixels (bed 8, 9 / 9, 10) for ``LakeTableBuilder`` (section 8).

WHAT IT MEASURES
  1. MF6 ignores BELEV/TELEV on a VERTICAL connection: 3.0, 6.0 and 0.0 all hold
     VOLUME 400.0 = 100 ft2 x (12 - 8), counted from the cell top. A HORIZONTAL face
     on the SAME cell from the bed (3.0) to the top -- what ``"bathy"`` writes -- is
     read: with the lake held at 6.0, below the top, the VERTICAL connection carries
     0.0 and holds VOLUME 0.0 while in a steady period the face carries 3.4641 out of
     the lake (0.0 when sealed). In a TRANSIENT period MF6 caps what leaves at what is
     available: 0.0 out, 2.0 with an INFLOW of 2.0, and -5.1475 (in) with the aquifer
     at 7.0.
  2. MF6's table checks (``gwf-lak.f90``, ``lak_read_table``) against ``LakeTable``'s:
     an equal stage, an equal volume and a falling area stop MF6 and are refused by
     ``LakeTable``; an equal area passes both. MF6 RUNS a negative volume or area
     (VOLUME -2.0, SURFACE-AREA -3.0), and a NaN reaches it only as a row FloPy cut
     short, which it cannot read; ``LakeTable`` refuses all three.
  3. Beyond the table MF6 grows the volume at the FULL last-row area, so a
     ``storage_coefficient`` table stopped at 11 gives 10.0 where one tabulated to 13
     gives 7.0 at stage 12.5.
  4. A SPECIFIED outlet's RATE: +5 raises the stage, -5 lowers it. A period left out
     keeps the previous period's value (rainfall and rate both).
  5. STATUS in lower case is accepted; an unknown one, and a negative rainfall,
     evaporation, runoff, inflow or withdrawal, are refused by name. LakeConnection's
     defaults: CONNLEN 0 is refused, CONNWIDTH 0 carries nothing, BELEV == TELEV is the
     whole cell, a TELEV above the cell top is refused, and a BEDLEAK of 0 carries
     nothing.
  6. What each ``connection_modes`` writes, and that a rectangular lake whose
     ``lake_bottom`` is outside ``only_layer`` is refused as a ``"bathy"`` one is (it
     used to lose its VERTICAL connection silently).
  7. What ``mf.lak`` refuses or reads now: an integer key past the last period (it
     wrote nothing), ``nlakes=`` in ``**options`` (it replaced the lake count), and a
     numpy array, now split by period -- where MF6, given the array whole, dies
     reading it (SIGSEGV, for a RATE and a RAINFALL). A Series keyed by integer lake
     ids out of lake order gives each lake its own value, as the dict does (it was
     read by position); a whole-grid ``lake_bottom`` Series over lakes numbered 1 and
     2 stays per cell, where a label reading would give each lake one bed.
  8. What ``LakeTableBuilder``/``LakeTable`` accept or refuse: a ``str`` footprint, a
     footprint in another CRS (reprojected), one off the raster or narrower than a
     pixel (refused, where it built a table of zeros), a GeoPackage's second layer
     (``footprint_layer``), a DataFrame, a ``storage_coefficient`` out of (0, 1], an
     integer ``tables`` key (it finds its lake).
  9. ``mf.lak.flopy``: an unknown option is a ``TypeError``, and a count that
     disagrees with its block is refused at the call. Beside each refusal, what MF6
     does with the same records written unchecked (through ``lak_spec``): FloPy
     RAISES ``NOUTLETS``/``NTABLES`` to the rows given but never lowers them, and
     writes ``nlakes`` as given, so a count above its rows (or with no block) and an
     ``nlakes`` either side of ``packagedata`` reach MF6, which refuses each by name;
     a boundname column without ``boundnames=True`` fails the FloPy build. The
     auxiliary values before a boundname are counted by NAME: ``[("a", "b")]``,
     ``["auxiliary", "a", "b"]`` and FloPy's own ``lak.auxiliary.get_data()``
     recarray (one row, ``("auxiliary", "a", "b")``) each write ``auxiliary  a  b``
     and run, so a row carrying two values is accepted under each (``len()`` counted
     two of them as one name and refused the row). And a
     ``pname`` through ``mf.lak``'s ``**options`` is refused -- it renamed only
     FloPy's package, while the spec, the ``lak_index`` a mover reads and the table
     children kept ``name``.
 10. The DEM mode keeps ONE zero-volume row of several: MF6 cannot read the table
     the requested stages make (two zero volumes), and reads the kept table as the
     requested curve at every stage, below the kept first row included.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path

import flopy
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import box

import myflopy as mf
from myflopy import LAKBuilder, LakeOutlet, LakeTable, LakeTableBuilder, ModelContext
from myflopy.advanced import lak_spec
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus

MF6 = os.environ.get("MF6", "mf6")


def _bench(ws: Path, *, connections, perioddata, stage=12.0, nper=1, transient=False,
           noutlets=0, outlets=None, table=None, obs=("STAGE", "VOLUME", "SURFACE-AREA"),
           aquifer=5.0):
    """Write and run the MF6 bench; return (ok, last obs row or None, first error line).

    ``aquifer`` is both the starting head and the corner CHD.
    """
    shutil.rmtree(ws, ignore_errors=True)
    ws.mkdir(parents=True)
    sim = flopy.mf6.MFSimulation(sim_name="t", sim_ws=str(ws), exe_name=MF6)
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=[(1.0, 1, 1.0)] * nper)
    flopy.mf6.ModflowIms(sim, complexity="SIMPLE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=3, ncol=3, delr=10, delc=10, top=8.0, botm=-10.0)
    flopy.mf6.ModflowGwfic(gwf, strt=aquifer)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0, icelltype=1)
    if transient:
        flopy.mf6.ModflowGwfsto(gwf, sy=0.1, ss=1e-5, iconvert=1, transient={0: True})
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 0), aquifer]])
    kw = {}
    if noutlets:
        kw.update(noutlets=noutlets, outlets=outlets)
    if table is not None:
        kw.update(ntables=1, tables=[[0, "t.lak.tab"]])
    # a name is a lake-wide observation; ("q", "LAK", 1, 1) is connection 1's flow
    records = [o if isinstance(o, tuple) else (o.lower().replace("-", ""), o, 1) for o in obs]
    flopy.mf6.ModflowGwflak(
        gwf, nlakes=1, packagedata=[[0, stage, len(connections)]],
        connectiondata=[[0, i, *c] for i, c in enumerate(connections)],
        perioddata=perioddata, observations={"lak.obs.csv": records}, **kw)
    if table is not None:
        flopy.mf6.ModflowUtllaktab(gwf, nrow=len(table), ncol=3, table=[list(r) for r in table],
                                   filename="t.lak.tab", pname="tab")
    sim.write_simulation(silent=True)
    proc = subprocess.run([MF6], cwd=ws, capture_output=True, text=True)
    lst = (ws / "mfsim.lst").read_text() if (ws / "mfsim.lst").exists() else ""
    error = ""
    if "ERROR REPORT:" in lst:
        error = " ".join(line.strip() for line in lst.split("ERROR REPORT:")[1].splitlines()[:3]).strip()
    elif proc.returncode < 0:
        error = f"killed by signal {-proc.returncode}"
    ok = proc.returncode == 0 and "Normal termination" in lst
    rows = pd.read_csv(ws / "lak.obs.csv").round(4).to_dict("records") if ok else None
    return ok, rows, error


VERTICAL = [(0, 1, 1), "VERTICAL", 0.001, 0.0, 0.0, 0.0, 0.0]
CONSTANT_12 = {0: [[0, "STATUS", "CONSTANT"], [0, "STAGE", 12.0]]}


def section_1(ws: Path) -> None:
    print("1. BELEV/TELEV written on a VERTICAL connection (cell top 8.0, lake at 12.0)")
    for bed in (3.0, 6.0, 0.0):
        conn = [(0, 1, 1), "VERTICAL", 0.001, bed, bed, 0.0, 0.0]
        ok, rows, _ = _bench(ws / f"s1_{bed}", connections=[conn], perioddata=CONSTANT_12)
        print(f"   written {bed:4}: VOLUME {rows[-1]['VOLUME']}  SURFACE-AREA {rows[-1]['SURFACEAREA']}")
    # What "bathy" writes on a lake cell whose bed (3.0) lies inside it: the VERTICAL
    # connection plus a HORIZONTAL face on the SAME cell from the bed to the top.
    print("   ... and a HORIZONTAL face on the same cell, belev 3.0 (the bed) to telev 8.0,")
    print("       both bedleak 1.0, lake held at 6.0 -- BELOW the cell top (positive = out of the lake)")
    held = [[0, "STATUS", "CONSTANT"], [0, "STAGE", 6.0]]
    cases = (
        ("steady, aquifer ~5", 1.0, False, 5.0, []),
        ("steady, face sealed (bedleak 1e-30)", 1e-30, False, 5.0, []),
        ("TRANSIENT, aquifer ~5", 1.0, True, 5.0, []),
        ("TRANSIENT, INFLOW 2.0", 1.0, True, 5.0, [[0, "INFLOW", 2.0]]),
        ("TRANSIENT, aquifer 7.0 (above the stage)", 1.0, True, 7.0, []),
    )
    for label, horizontal_leak, transient, aquifer, extra in cases:
        vertical = [(0, 1, 1), "VERTICAL", 1.0, 3.0, 3.0, 0.0, 0.0]
        horizontal = [(0, 1, 1), "HORIZONTAL", horizontal_leak, 3.0, 8.0, 5.0, 10.0]
        tag = "".join(ch for ch in label if ch.isalnum())
        ok, rows, error = _bench(
            ws / f"s1_h_{tag}", connections=[vertical, horizontal], stage=6.0,
            transient=transient, aquifer=aquifer, perioddata={0: held + extra},
            obs=("VOLUME", ("qv", "LAK", 1, 1), ("qh", "LAK", 1, 2)))
        print(f"   {label}: VERTICAL {rows[-1]['QV']}  HORIZONTAL {rows[-1]['QH']}  "
              f"VOLUME {rows[-1]['VOLUME']}")


def section_2(ws: Path) -> None:
    print("2. MF6's table checks against LakeTable's (lake held at 8.5 then 9.5)")
    cases = {
        "equal stage": [(8.0, 0.0, 0.0), (8.0, 1.0, 4.0), (10.0, 5.0, 4.0)],
        "equal volume": [(8.0, 0.0, 0.0), (9.0, 0.0, 0.0), (10.0, 5.0, 4.0)],
        "falling area": [(8.0, 0.0, 5.0), (9.0, 1.0, 4.0), (10.0, 5.0, 4.0)],
        "equal area": [(8.0, 0.0, 4.0), (9.0, 1.0, 4.0), (10.0, 5.0, 4.0)],
        "negative volume": [(8.0, -5.0, 0.0), (9.0, 1.0, 4.0), (10.0, 5.0, 4.0)],
        "negative area": [(8.0, 0.0, -4.0), (9.0, 1.0, -2.0), (10.0, 5.0, 4.0)],
        "NaN volume": [(8.0, 0.0, 0.0), (9.0, float("nan"), 4.0), (10.0, 5.0, 4.0)],
    }
    held = {i: [[0, "STATUS", "CONSTANT"], [0, "STAGE", s]] for i, s in enumerate((8.5, 9.5))}
    for label, rows in cases.items():
        try:
            LakeTable(rows)
            verdict = "accepts"
        except ValueError as error:
            verdict = f"refuses ({str(error)[:60]}...)"
        ok, obs, error = _bench(ws / f"s2_{label.replace(' ', '_')}", connections=[VERTICAL],
                                perioddata=held, nper=2, stage=8.5, table=rows)
        mf6 = (f"ran, VOLUME/SURFACE-AREA {[(r['VOLUME'], r['SURFACEAREA']) for r in obs]}" if ok
               else f"stopped: {error[:100]}")
        print(f"   {label}: LakeTable {verdict}; MF6 {mf6}")


def section_3(ws: Path, bed: Path) -> None:
    print("3. storage_coefficient=0.5 beyond the table, lake held at 12.5")
    for label, stages in (("to 13", [8, 10, 11, 12, 13]), ("stopped at 11", [8, 10, 11])):
        rows = LakeTableBuilder(dem=bed, footprint=box(0, 0, 2, 2), stages=stages,
                                storage_coefficient=0.5).build().rows
        ok, obs, _ = _bench(ws / f"s3_{len(stages)}", connections=[VERTICAL], stage=12.5,
                            perioddata={0: [[0, "STATUS", "CONSTANT"], [0, "STAGE", 12.5]]},
                            table=rows)
        print(f"   table {label}: rows {rows[-2:]}  VOLUME {obs[-1]['VOLUME']}  (normal={ok})")


def section_4(ws: Path) -> None:
    print("4. Outlet RATE sign, and a period left out (transient, 2 periods)")
    outlet = [[0, 0, -1, "SPECIFIED", 0.0, 0.0, 0.0, 0.0]]
    leaky = [(0, 1, 1), "VERTICAL", 1e-9, 0.0, 0.0, 0.0, 0.0]
    for rate in (5.0, -5.0):
        ok, rows, _ = _bench(ws / f"s4_{rate}", connections=[leaky], stage=10.0, nper=2,
                             transient=True, perioddata={0: [[0, "RATE", rate]], 1: []},
                             noutlets=1, outlets=outlet, obs=("STAGE", "EXT-OUTFLOW"))
        print(f"   RATE {rate:+}: stage {[r['STAGE'] for r in rows]}, "
              f"EXT-OUTFLOW {[r['EXTOUTFLOW'] for r in rows]} (period 2 block empty)")
    ok, rows, _ = _bench(ws / "s4_rain", connections=[leaky], stage=10.0, nper=2, transient=True,
                         perioddata={0: [[0, "RAINFALL", 1.0]], 1: []}, obs=("STAGE", "RAINFALL"))
    print(f"   RAINFALL 1.0 in period 1 only: RAINFALL {[r['RAINFALL'] for r in rows]}")


def section_5(ws: Path) -> None:
    print("5. What MF6 accepts and refuses by name")
    for label, pdata in (("STATUS 'constant'", {0: [[0, "STATUS", "constant"], [0, "STAGE", 12.0]]}),
                         ("STATUS 'BOGUS'", {0: [[0, "STATUS", "BOGUS"]]})):
        ok, _, error = _bench(ws / f"s5_{label[8:12]}", connections=[VERTICAL], perioddata=pdata)
        print(f"   {label}: ran={ok} {error[:70]}")
    for keyword in ("RAINFALL", "EVAPORATION", "RUNOFF", "INFLOW", "WITHDRAWAL"):
        ok, _, error = _bench(ws / f"s5_{keyword}", connections=[VERTICAL],
                              perioddata={0: [[0, keyword, -1.0]]})
        print(f"   {keyword} -1: ran={ok} {error[:80]}")
    for label, conn in (
        ("HORIZONTAL, CONNLEN 0", [(0, 1, 1), "HORIZONTAL", 0.1, 0.0, 0.0, 0.0, 0.0]),
        ("HORIZONTAL, CONNWIDTH 0", [(0, 1, 1), "HORIZONTAL", 0.1, 0.0, 0.0, 5.0, 0.0]),
        ("HORIZONTAL, BELEV == TELEV", [(0, 1, 1), "HORIZONTAL", 0.1, 0.0, 0.0, 5.0, 10.0]),
        ("HORIZONTAL, whole cell", [(0, 1, 1), "HORIZONTAL", 0.1, -10.0, 8.0, 5.0, 10.0]),
        ("HORIZONTAL, TELEV above top", [(0, 1, 1), "HORIZONTAL", 0.1, 0.0, 9.0, 5.0, 10.0]),
        ("VERTICAL, BEDLEAK 0", [(0, 1, 1), "VERTICAL", 0.0, 0.0, 0.0, 0.0, 0.0]),
    ):
        tag = "".join(ch for ch in label if ch.isalnum())
        ok, rows, error = _bench(ws / f"s5_{tag}", connections=[conn], perioddata=CONSTANT_12,
                                 obs=[("q", "LAK", 1, 1)])
        result = f"connection flow {rows[-1]['Q']:.4f}" if ok else error[:90]
        print(f"   {label}: ran={ok} {result}")


def _toy_grid() -> VoronoiGridPlus:
    n = 4
    verts = np.array([[x, y] for y in range(n + 1) for x in range(n + 1)], float)
    iverts, centers = [], []
    for r in range(n):
        for c in range(n):
            v0 = r * (n + 1) + c
            iverts.append([v0, v0 + 1, v0 + n + 2, v0 + n + 1])
            centers.append([c + 0.5, r + 0.5])
    grid = VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=np.array(centers))
    grid.gdf_topbtm = gpd.GeoDataFrame({0: [12.0] * 16, 1: [8.0] * 16, 2: [0.0] * 16},
                                       geometry=grid.gdf_vorPolys.geometry, crs=grid.crs)
    grid.get_disu_connectivity(validate=True)
    return grid


def section_6() -> None:
    print("6. What each connection mode WRITES (4 x 4 unit cells, 2 x 2 lake, layer tops 12/8/0)")
    grid = _toy_grid()
    lakes = gpd.GeoDataFrame({"name": ["pond"]}, geometry=[box(1, 1, 3, 3)], crs=grid.crs)
    ctx = ModelContext(grid=grid, domain=np.ones((2, 16), dtype=int))
    bed = pd.Series({c: 12.0 for c in range(16)})
    bed[5], bed[6], bed[9], bed[10] = 9.0, 10.0, 10.0, 11.0
    base = dict(context=ctx, nper=1, lakes=lakes, lake_id_field="name", starting_stage=11.5,
                bed_leakance=0.1)
    cases = {
        "bathy, per-cell bed": dict(lake_bottom=bed),
        "rectangular": dict(connection_modes="rectangular", lake_bottom=9.0, lake_top=11.0),
        "rectangular, interior": dict(connection_modes="rectangular", lake_bottom=9.0,
                                      lake_top=11.0, rectangular_interior=True),
        "rectangular, bottom in layer 1, only_layer=0": dict(
            connection_modes="rectangular", lake_bottom=6.0, lake_top=11.0, only_layer=0),
    }
    for label, kw in cases.items():
        with warnings.catch_warnings():
            warnings.simplefilter("error")  # a silent drop must not warn
            try:
                rows = LAKBuilder(**{**base, **kw}).connectiondata
            except ValueError as error:
                print(f"   {label}: ValueError: {str(error)[:120]}")
                continue
        kinds = pd.Series([r[3] for r in rows]).value_counts().to_dict()
        own_cells = all(r[2][1] in (5, 6, 9, 10) for r in rows)
        print(f"   {label}: {kinds}; every connection on a lake cell: {own_cells}")
        print(f"      first rows: {[(r[2], r[3], r[5], r[6], round(r[7], 2), round(r[8], 2)) for r in rows[:3]]}")


def section_7(ws: Path) -> None:
    print("7. What mf.lak refuses or reads now, and what MF6 does with an array written whole")
    grid = _toy_grid()
    lakes = gpd.GeoDataFrame({"num": [101]}, geometry=[box(1, 1, 3, 3)], crs=grid.crs)
    ctx = ModelContext(grid=grid, domain=np.ones((2, 16), dtype=int))
    base = dict(context=ctx, nper=2, lakes=lakes, lake_id_field="num", starting_stage=11.5,
                lake_bottom=9.0, connection_modes="rectangular", lake_top=11.0, bed_leakance=0.1,
                unexcavated="ignore")

    def attempt(label, fn):
        try:
            print(f"   {label}: {fn()}")
        except (ValueError, TypeError) as error:
            print(f"   {label}: {type(error).__name__}: {str(error)[:130]}")

    attempt("rainfall={101: 0.5}", lambda: mf.lak(**base, rainfall={101: 0.5}).options["perioddata"])
    attempt("rainfall={'101': 0.5}", lambda: mf.lak(**base, rainfall={"101": 0.5}).options["perioddata"])
    attempt("**options nlakes=5 with one lake", lambda: mf.lak(**base, nlakes=5).options["nlakes"])
    attempt("LakeOutlet(rate=np.array([-1, -2]))", lambda: mf.lak(
        **base, outlets=[LakeOutlet("101", rate=np.array([-1.0, -2.0]))]).options["perioddata"])
    ok, _, error = _bench(ws / "s7_array", connections=[VERTICAL], stage=10.0,
                          perioddata={0: [[0, "RATE", np.array([-1.0, -2.0])]]}, noutlets=1,
                          outlets=[[0, 0, -1, "SPECIFIED", 0.0, 0.0, 0.0, 0.0]], obs=())
    written = [line.strip() for line in (ws / "s7_array" / "m.lak").read_text().splitlines() if "RATE" in line]
    print(f"   MF6 on {written}: ran={ok} {error}")
    ok, _, error = _bench(ws / "s7_rain_array", connections=[VERTICAL], stage=10.0,
                          perioddata={0: [[0, "RAINFALL", np.array([0.1, 0.2])]]}, obs=())
    written = [line.strip() for line in (ws / "s7_rain_array" / "m.lak").read_text().splitlines()
               if "RAINFALL" in line]
    print(f"   MF6 on {written}: ran={ok} {error}")

    # A Series keyed by integer lake ids is read by lake, as the dict is; to lake_bottom an
    # integer is a cell number too, so a whole-grid bed stays per cell.
    pair = gpd.GeoDataFrame({"num": [101, 102]}, geometry=[box(0, 0, 2, 2), box(2, 2, 4, 4)], crs=grid.crs)
    two = {**base, "lakes": pair}
    for label, stage in (("dict {102: 12.0, 101: 11.5}", {102: 12.0, 101: 11.5}),
                         ("pd.Series({102: 12.0, 101: 11.5})", pd.Series({102: 12.0, 101: 11.5}))):
        attempt(f"lakes 101, 102; starting_stage={label}", lambda stage=stage: [
            (row[3], row[1]) for row in mf.lak(**{**two, "starting_stage": stage}).options["packagedata"]])
    numbered = gpd.GeoDataFrame({"num": [1, 2]}, geometry=[box(0, 0, 2, 2), box(2, 2, 4, 4)], crs=grid.crs)
    whole = pd.Series({cell: 9.0 + 0.1 * cell for cell in range(16)})
    rows = mf.lak(**{**base, "lakes": numbered, "lake_bottom": whole, "connection_modes": "bathy",
                     "nper": 1}).options["connectiondata"]
    floors = {lake: sorted(round(r[5], 3) for r in rows if r[3] == "VERTICAL" and r[0] == lake) for lake in (0, 1)}
    print(f"   lakes 1, 2; a whole-grid lake_bottom Series: VERTICAL beds {floors}; read by label, "
          f"each lake would take ONE bed, {whole[1]:.1f} and {whole[2]:.1f}")


def section_8(ws: Path, bed: Path) -> None:
    print("8. What LakeTableBuilder / LakeTable accept or refuse")
    full = box(0, 0, 2, 2)
    fp = ws / "fp.gpkg"
    gpd.GeoDataFrame(geometry=[full], crs=2927).to_file(fp)
    fp4326 = ws / "fp4326.gpkg"
    gpd.GeoDataFrame(geometry=[full], crs=2927).to_crs(4326).to_file(fp4326)
    two = ws / "two_layers.gpkg"
    gpd.GeoDataFrame(geometry=[box(0, 0, 1, 2)], crs=2927).to_file(two, layer="a")
    gpd.GeoDataFrame(geometry=[box(1, 0, 2, 2)], crs=2927).to_file(two, layer="b")

    def attempt(label, fn):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                print(f"   {label}: {fn()}")
            except (ValueError, TypeError) as error:
                print(f"   {label}: {type(error).__name__}: {str(error)[:150]}")
        for warning in caught:
            if "More than one layer" not in str(warning.message):
                print(f"      warns: {str(warning.message)[:120]}")

    stages = [9, 10, 11]
    attempt("footprint=Path", lambda: LakeTableBuilder(dem=bed, footprint=fp, stages=stages).build().rows)
    attempt("footprint=str", lambda: LakeTableBuilder(dem=bed, footprint=str(fp), stages=stages).build().rows)
    attempt("footprint file in EPSG:4326 (reprojected)", lambda: LakeTableBuilder(
        dem=bed, footprint=fp4326, stages=stages).build().rows)
    attempt("footprint geometry in EPSG:4326 coordinates, no CRS", lambda: LakeTableBuilder(
        dem=bed, footprint=box(-122.5, 47.0, -122.4, 47.1), stages=stages).build().rows)
    attempt("footprint off the raster", lambda: LakeTableBuilder(
        dem=bed, footprint=box(100, 100, 101, 101), stages=stages).build().rows)
    attempt("footprint narrower than a pixel centre", lambda: LakeTableBuilder(
        dem=bed, footprint=box(0, 0, 0.4, 2), stages=stages).build().rows)
    attempt("two-layer GeoPackage, default (a = left column)", lambda: LakeTableBuilder(
        dem=bed, footprint=two, stages=stages).build().rows)
    attempt("two-layer GeoPackage, footprint_layer='b' (right column)", lambda: LakeTableBuilder(
        dem=bed, footprint=two, footprint_layer="b", stages=stages).build().rows)
    frame = pd.DataFrame({"stage": [96.0, 98.0], "volume": [0.0, 1.0], "area": [0.0, 2.0]})
    attempt("rows=DataFrame (named)", lambda: LakeTableBuilder(rows=frame).build().rows)
    attempt("rows=DataFrame (labels 0, 1, 2)", lambda: LakeTableBuilder(
        rows=pd.DataFrame([[1, 0, 0], [2, 1, 1]])).build().rows)
    for coefficient in (-1, 0, 1.5):
        attempt(f"storage_coefficient={coefficient} (prism)", lambda c=coefficient: LakeTableBuilder(
            area=1, bottom=9, top=11, storage_coefficient=c).build().rows)
    attempt("LakeTable with a NaN and a negative volume", lambda: LakeTable(
        [(1, -5.0, 0.0), (2, float("nan"), 1.0)]).rows)

    grid = _toy_grid()
    lakes = gpd.GeoDataFrame({"num": [101]}, geometry=[box(1, 1, 3, 3)], crs=grid.crs)
    ctx = ModelContext(grid=grid, domain=np.ones((2, 16), dtype=int))
    table = LakeTableBuilder(area=4.0, bottom=9.0, top=12.0)
    base = dict(context=ctx, nper=1, lakes=lakes, lake_id_field="num", starting_stage=11.5,
                lake_bottom=9.0, connection_modes="rectangular", lake_top=11.0, unexcavated="ignore")
    attempt("tables={'101': ...}", lambda: mf.lak(**base, tables={"101": table}).options["tables"])
    attempt("tables={101: ...}", lambda: mf.lak(**base, tables={101: table}).options["tables"])


def _flopy_bench(ws: Path, **overrides) -> str:
    """The MF6 bench with its lake from ``mf.lak.flopy(**overrides)``; one summary line."""
    shutil.rmtree(ws, ignore_errors=True)
    ws.mkdir(parents=True)
    sim = flopy.mf6.MFSimulation(sim_name="t", sim_ws=str(ws), exe_name=MF6)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="SIMPLE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=3, ncol=3, delr=10, delc=10, top=8.0, botm=-10.0)
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0, icelltype=1)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 0), 5.0]])
    kw = dict(packagedata=[[0, 10.0, 1]], connectiondata=[[0, 0, *VERTICAL]],
              perioddata={0: [[0, "STATUS", "ACTIVE"]]})
    kw.update(overrides)
    try:
        mf.lak.flopy(**kw)
        verdict = "mf.lak.flopy accepts"
    except ValueError as error:
        verdict = f"mf.lak.flopy refuses ({str(error).splitlines()[1].strip()[:70]}...)"
    # The same records written UNCHECKED -- lak_spec is what mf.lak.flopy wraps -- so what
    # MF6 does with them is still measured.
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            lak_spec(kw.pop("packagedata"), kw.pop("connectiondata"), kw.pop("perioddata"), **kw).build(gwf)
    except flopy.mf6.mfbase.MFDataException as error:
        return f"{verdict}; unchecked, the build raised MFDataException ({str(error).splitlines()[0][:60]}...)"
    sim.write_simulation(silent=True)
    (ws / "t.lak.tab").write_text("BEGIN DIMENSIONS\n  NROW 3\n  NCOL 3\nEND DIMENSIONS\nBEGIN TABLE\n"
                                  "  8.0 0.0 0.0\n  10.0 200.0 100.0\n  14.0 600.0 100.0\nEND TABLE\n")
    proc = subprocess.run([MF6], cwd=ws, capture_output=True, text=True)
    lst = (ws / "mfsim.lst").read_text()
    written = [line.strip() for line in (ws / "m.lak").read_text().splitlines()
               if line.strip().startswith(("NLAKES", "NOUTLETS", "NTABLES", "auxiliary"))]
    error = ""
    if "ERROR REPORT:" in lst:
        error = " ".join(line.strip() for line in lst.split("ERROR REPORT:")[1].splitlines()[:4]).strip()
    return f"{verdict}; unchecked, writes {written}; ran={proc.returncode == 0} {error[:90]}"


def section_9(ws: Path) -> None:
    print("9. mf.lak.flopy: the counts against the blocks, boundnames; a pname through mf.lak's **options")
    try:
        mf.lak.flopy(packagedata=[[0, 10.0, 1]], connectiondata=[], perioddata={}, print_stagee=True)
    except TypeError as error:
        print(f"   print_stagee=True: TypeError: {str(error)[:80]}")
    one = [[0, 0, -1, "SPECIFIED", 0.0, 0.0, 0.0, 0.0]]
    two = one + [[1, 0, -1, "SPECIFIED", 0.0, 0.0, 0.0, 0.0]]
    two_lakes = dict(packagedata=[[0, 10.0, 1], [1, 10.0, 1]],
                     connectiondata=[[0, 0, *VERTICAL], [1, 0, (0, 1, 2), *VERTICAL[1:]]])
    cases = {
        "nlakes=3, one packagedata row": dict(nlakes=3),
        "nlakes=1, two packagedata rows": dict(nlakes=1, **two_lakes),
        "two outlet rows, noutlets=0 (default)": dict(outlets=two),
        "two outlet rows, noutlets=1": dict(outlets=two, noutlets=1),
        "one outlet row, noutlets=2": dict(outlets=one, noutlets=2),
        "noutlets=1, no outlets": dict(noutlets=1),
        "one table, ntables=0 (default)": dict(tables=[[0, "t.lak.tab"]]),
        "one table, ntables=2": dict(tables=[[0, "t.lak.tab"]], ntables=2),
        "ntables=1, no tables": dict(ntables=1),
        "a boundname column, boundnames unset": dict(packagedata=[[0, 10.0, 1, "pond"]]),
        "a boundname column, boundnames=True": dict(packagedata=[[0, 10.0, 1, "pond"]], boundnames=True),
    }
    # FloPy's own shape for the auxiliary names, as `lak.auxiliary.get_data()` returns it.
    probe = flopy.mf6.MFSimulation(sim_name="p", sim_ws=str(ws))
    flopy.mf6.ModflowTdis(probe, nper=1)
    probe_gwf = flopy.mf6.ModflowGwf(probe, modelname="p")
    flopy.mf6.ModflowGwfdis(probe_gwf, nlay=1, nrow=3, ncol=3, top=8.0, botm=-10.0)
    flopys_own = flopy.mf6.ModflowGwflak(
        probe_gwf, nlakes=1, auxiliary=["a", "b"], packagedata=[[0, 10.0, 1, 5.0, 6.0]],
        connectiondata=[[0, 0, *VERTICAL]], perioddata={0: []}).auxiliary.get_data()
    print(f"   FloPy's own auxiliary.get_data(): {flopys_own.tolist()} (len {len(flopys_own)})")
    two_aux = [[0, 10.0, 1, 5.0, 6.0]]
    for label, auxiliary in (("['a', 'b']", ["a", "b"]), ("[('a', 'b')]", [("a", "b")]),
                             ("['auxiliary', 'a', 'b']", ["auxiliary", "a", "b"]),
                             ("FloPy's own get_data()", flopys_own)):
        cases[f"two auxiliary values, auxiliary={label}"] = dict(packagedata=two_aux, auxiliary=auxiliary)
    for label, overrides in cases.items():
        print(f"   {label}: {_flopy_bench(ws / ('s9_' + ''.join(c for c in label if c.isalnum())), **overrides)}")

    grid = _toy_grid()
    lakes = gpd.GeoDataFrame({"name": ["pond"]}, geometry=[box(1, 1, 3, 3)], crs=grid.crs)
    ctx = ModelContext(grid=grid, domain=np.ones((2, 16), dtype=int))
    kwargs = dict(context=ctx, nper=1, lakes=lakes, lake_id_field="name", starting_stage=11.5,
                  lake_bottom=9.0, connection_modes="rectangular", lake_top=11.0, unexcavated="ignore",
                  tables={"pond": LakeTableBuilder(area=4.0, bottom=9.0, top=12.0)})
    try:
        mf.lak(**kwargs, pname="zz")
    except TypeError as error:
        print(f"   mf.lak(..., pname='zz'): TypeError: {str(error)[:110]}")
    spec = mf.lak(**kwargs, name="zz")
    print(f"   mf.lak(..., name='zz'): FloPy pname {spec.options['pname']!r}, spec name {spec.name!r}, "
          f"lak_index package {spec.metadata['lak_index']['package']!r}, "
          f"table child {spec.options['lake_tables'][0]['pname']!r}")


def section_10(ws: Path, bed: Path) -> None:
    print("10. The DEM mode's zero-volume rows: MF6 on the requested list, and on the kept table")
    requested = [6.0, 7.0, 7.5, 8.0, 9.0, 10.0, 11.0]
    full = [(6.0, 0.0, 0.0), (7.0, 0.0, 0.0), (7.5, 0.0, 0.0), (8.0, 0.0, 0.0),
            (9.0, 1.0, 1.0), (10.0, 4.0, 3.0), (11.0, 8.0, 4.0)]
    kept = LakeTableBuilder(dem=bed, footprint=box(0, 0, 2, 2), stages=requested).build().rows
    print(f"   requested stages {requested} -> kept rows {kept}")
    held = [6.5, 7.25, 8.0, 8.5, 9.5, 10.75]
    pdata = {i: [[0, "STATUS", "CONSTANT"], [0, "STAGE", s]] for i, s in enumerate(held)}
    ok, _, error = _bench(ws / "s10_full", connections=[VERTICAL], perioddata=pdata, nper=len(held),
                          stage=held[0], table=full)
    print(f"   the requested list: MF6 ran={ok} {error[:90]}")
    ok, obs, error = _bench(ws / "s10_kept", connections=[VERTICAL], perioddata=pdata,
                            nper=len(held), stage=held[0], table=kept)
    stages, volumes, areas = zip(*full, strict=True)
    for stage, row in zip(held, obs, strict=True):
        print(f"   held at {stage}: MF6 VOLUME {row['VOLUME']} SURFACE-AREA {row['SURFACEAREA']}; "
              f"the requested curve {np.interp(stage, stages, volumes):.4f} "
              f"{np.interp(stage, stages, areas):.4f}")


def main() -> None:
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(tempfile.mkdtemp(prefix="lak_claims_"))
    root.mkdir(parents=True, exist_ok=True)
    bed = root / "bed.tif"
    with rasterio.open(bed, "w", driver="GTiff", height=2, width=2, count=1, dtype="float64",
                       crs="EPSG:2927", transform=from_origin(0, 2, 1, 1)) as dst:
        dst.write(np.array([[8.0, 9.0], [9.0, 10.0]]), 1)
    section_1(root)
    section_2(root)
    section_3(root, bed)
    section_4(root)
    section_5(root)
    section_6()
    section_7(root)
    section_8(root, bed)
    section_9(root)
    section_10(root, bed)
    print(f"workspace: {root}")


if __name__ == "__main__":
    main()
