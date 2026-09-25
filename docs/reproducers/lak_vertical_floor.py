"""Does a VERTICAL LAK connection honour the bed elevation written for it? (No -- and the cut fixes it.)

The measurement behind ledgers 320, 325 and 327: why `mf.lak` on an un-excavated layer stack
leaves most of a pond's floor exchanging nothing, and why `LayerStack.excavate` is the fix
rather than a different connection record. Run it::

    /home/lukem/python/envs/gw/bin/python docs/reproducers/lak_vertical_floor.py [workspace]

Needs ``mf6`` on PATH (measured with 6.8.0.dev0+0418efe4; case g7 needs a build with LAK
``IMPLICIT``). About 10 s. Writes nothing outside the workspace it is given (a temporary
directory when none is).

THE CLAIM. A VERTICAL lake-GWF connection exchanges through the TOP face of its GWF cell.
MF6 overwrites the BELEV/TELEV written for it with ``belev = cell top`` and ``telev = top +
SURFDEP`` (gwf-lak.f90, ``lak_read_lake_connections``: ``this%belev(ipos) = top``). So a
vertical connection placed in the layer that CONTAINS the lakebed -- the bed partway down
the cell, which is what ``"bathy"`` writes on a stack whose top is the ground -- leaks
nothing until the stage passes that cell's TOP. A lake table changes volume and area, not
this.

PART 1 -- the claim, on MF6 alone (raw FloPy, no myflopy). A 2-layer 3 x 3 DIS of 10 ft
cells: top 10, botm [5, 0], K = K33 = 10, confined unless stated, the layer-2 OUTER RING at
CHD 1.0 so the head under the lake stays far below every bed tested and a wet lake should
LOSE water. One lake, BEDLEAK 1.0, every connection on the centre column. Each row prints
the measured lake-GWF flow beside two hand predictions from MF6's own exchange formula
``q = C (max(h, belev) - max(stage, belev))`` (``lak_calculate_conn_exchange``; C from
``lak_calculate_satcond``): ``claim`` forces belev to the cell TOP for a vertical
connection, ``honour`` takes the BELEV written. Measured (every case matches ``claim``):

  - bed 7 in a cell topped at 10, stage 8: q = 0.000 (``honour`` says -80.000), and MF6
    prints the stage as 10.000 and writes DRY to the stage file: its lake bottom is the
    cell top, 10, not the bed (a, c5; the same at bed 9.9 and at 0.5, below the cell's own
    bottom -- c1, c3). Stage 12: -160.000 in all three, whatever bed was written (b, c2,
    c4).
  - a TAB6 table does not change it (e1, e3); nor does NEWTON (g1, g2), icelltype 1 (g3)
    or LAK IMPLICIT (g7). An ACTIVE lake with 80 ft3/d of inflow stands at 11.000, the
    formula from the cell top, not the 8.000 from the bed (e2). SURFDEP 1.0 counts from
    the cell top too (g6, NEWTON), and under Picard it carried nothing even at stage 12
    (g5).
  - a vertical to LAYER 2 counts from layer 2's top, 5 (d1-d4), with the layer-1 cell
    above active or not; a vertical to an INACTIVE cell is refused by MF6 (g8).
  - a HORIZONTAL connection on the same cell DOES honour its elevations: 7 to 10 carries
    -20.000 at stage 8 (f), and beside a vertical (what ``"bathy"`` writes at a step) it is
    the ONLY one that moves water (h).

PART 2 -- the same bench through myflopy, on a DISV of the same nine 10 ft squares
(``LayerStack(top=mf.Flat(10)).add("upper", bottom=mf.Flat(5)).add("lower",
bottom=mf.Flat(0))``), the lake the centre cell, bed 7, bedleak 1, CONSTANT stage 8 and 12.

  L  the un-excavated stack and ``mf.lak(lake_bottom=pd.Series({4: 7.0}))`` -- the bed
     per cell, as a real pond's bathymetry is, the land around it at the model top --
     under the old rule (which the ``unexcavated=`` guard now reports): a vertical into
     layer 1, whose top is 10, and four walls on the lake's OWN cell from 7 to 10.
  E  ``stack.excavate(bottom=7.0, within=<the centre cell>)`` and ``mf.lak`` with no
     ``lake_bottom`` (the record's): the centre column's top IS 7, the layer-1 cell under
     the bed runs 7 to 5, and the four walls go INTO the neighbours' layer-1 cells.

Measured (ft3/d, + into the lake; GWF discrepancy 0.00 in all four):

  stage 8   L  vertical    0.000 (its cell top is 10)   walls  -60.000 on the lake cell
            E  vertical  -90.909 (predicted -90.909)    walls  -80.000 into cells 1/3/5/7
  stage 12  L  vertical -160.000 (from 10)              walls -300.000
            E  vertical -454.545 (from 7, predicted)    walls -400.000

E's vertical is MF6's formula from the BED exactly, ``C = 100 / (1/1 + 0.5 * 2 / 10)``
through the 2 ft cell the cut leaves; L's counts from the cell top, so at stage 8 the
floor carries nothing and the lake's only exchange is through faces written on its own
column -- which, on a real pond, is lake water (ledger 327: 411 of 612 vertical
connections carried exactly 0 all run). MF6's lake bottom is 7.000 in both, because L's
walls reach down to the bed: the lake is wet, and its floor still exchanges nothing.
"""

from __future__ import annotations

import contextlib
import io
import re
import shutil
import sys
import tempfile
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import flopy
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

import myflopy as mf
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.lak_implicit import insert_implicit

# ---- the aquifer ------------------------------------------------------------------------
NLAY, NROW, NCOL = 2, 3, 3
DELR = DELC = 10.0
TOP = 10.0
BOTM = [5.0, 0.0]
K = K33 = 10.0
BEDLEAK = 1.0
CHD_HEAD = 1.0  # layer-2 outer ring; keeps heads far below the beds tested
CENTER = (1, 1)  # (row, col), zero-based
AREA = DELR * DELC
DRY = -1.0e30  # MF6 DHDRY, written to the stage file for a lake at/below its bottom


@dataclass
class Conn:
    layer: int  # zero-based
    ctype: str  # "vertical" | "horizontal"
    belev: float
    telev: float
    connlen: float = 0.0
    connwidth: float = 0.0

    @property
    def cellid(self):
        return (self.layer, *CENTER)

    @property
    def tag(self):
        return f"{'V' if self.ctype == 'vertical' else 'H'}L{self.layer + 1}"


@dataclass
class Case:
    name: str
    conns: list[Conn]
    stage: float | None = None  # CONSTANT lake stage; None -> ACTIVE lake
    inflow: float | None = None  # ACTIVE lake: specified INFLOW (ft3/d)
    table: bool = False  # add a TAB6 lake table, stages 7..12, prism on 100 ft2
    surfdep: float | None = None
    newton: bool = False
    icelltype: int = 0
    l1_center_active: bool = True
    implicit: bool = False
    chd_head: float = CHD_HEAD  # raise it to put the aquifer head ABOVE the cell top
    plain_formula: bool = True  # False where saturation terms make it inapplicable
    note: str = ""
    result: dict = field(default_factory=dict)


# ---- hand predictions from MF6's own formula --------------------------------------------
def _cell_top_bot(layer):
    tops = [TOP, *BOTM]
    return tops[layer], tops[layer + 1]


def satcond(conn: Conn, thickness: float | None = None) -> float:
    """MF6 ``lak_calculate_satcond`` for a confined cell with a BEDLEAK."""

    top, bot = _cell_top_bot(conn.layer)
    thick = top - bot if thickness is None else thickness
    clb = 1.0 / BEDLEAK
    if conn.ctype == "vertical":
        return AREA / (clb + 0.5 * thick / K33)
    area = (conn.telev - conn.belev) * conn.connwidth
    return area / (clb + conn.connlen / K)


def _belev(conn: Conn, honour_written: bool) -> float:
    if conn.ctype == "vertical" and not honour_written:
        return _cell_top_bot(conn.layer)[0]
    return conn.belev


def predict_q(conn: Conn, stage: float, head: float, honour_written: bool) -> float:
    belev = _belev(conn, honour_written)
    return satcond(conn) * (max(head, belev) - max(stage, belev))


def predict_stage(conn: Conn, inflow: float, honour_written: bool) -> float:
    return _belev(conn, honour_written) + inflow / satcond(conn)


def _listing_numbers(listing: str) -> tuple[float, float]:
    """The stage MF6 used and MF6's own lake bottom, from its printed tables."""

    # a lake with a boundname prints its NAME before its number
    m = re.search(r"STAGES FOR EACH CONTROL VOLUME.*?-+\n.*?-+\n\s*(?:[A-Za-z]\S*\s+)?1\s+([-\d.E+]+)",
                  listing, re.DOTALL)
    stage = float(m.group(1))
    m = re.search(r"STAGE/VOLUME RELATION FOR LAKE\s+1.*?-+\n\s*([-\d.E+]+)", listing, re.DOTALL)
    return stage, float(m.group(1))


# ---- PART 1: build + run one raw-FloPy case -----------------------------------------------
def build_and_run(root: Path, case: Case) -> dict:
    ws = root / "part1" / case.name
    if ws.exists():
        shutil.rmtree(ws)
    ws.mkdir(parents=True)

    sim = flopy.mf6.MFSimulation(sim_name=case.name, sim_ws=str(ws), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)], time_units="days")
    nonlinear = case.newton or case.icelltype == 1 or case.implicit
    flopy.mf6.ModflowIms(
        sim,
        complexity="moderate" if nonlinear else "simple",
        linear_acceleration="bicgstab" if (case.newton or case.implicit) else "cg",
        outer_dvclose=1e-8, inner_dvclose=1e-10, outer_maximum=500, inner_maximum=300,
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf", save_flows=True,
                               newtonoptions="newton" if case.newton else None)
    idomain = np.ones((NLAY, NROW, NCOL), dtype=int)
    if not case.l1_center_active:
        idomain[0, CENTER[0], CENTER[1]] = 0
    # idomain is applied AFTER the LAK package is declared (below): FloPy refuses a
    # connection to an idomain-0 cell at declaration, and case g8 wants MF6's own
    # answer to that, not FloPy's.
    dis = flopy.mf6.ModflowGwfdis(gwf, nlay=NLAY, nrow=NROW, ncol=NCOL, delr=DELR, delc=DELC,
                                  top=TOP, botm=BOTM, idomain=np.ones_like(idomain))
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=case.icelltype, k=K, k33=K33)
    flopy.mf6.ModflowGwfic(gwf, strt=case.chd_head)
    ring = [(1, i, j) for i in range(NROW) for j in range(NCOL) if (i, j) != CENTER]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: [[c, case.chd_head] for c in ring]})
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord="gwf.hds", budget_filerecord="gwf.cbc",
                           saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])

    if case.stage is not None:
        strt, period = case.stage, [[0, "status", "constant"], [0, "stage", case.stage]]
    else:
        strt, period = 10.5, [[0, "status", "active"], [0, "inflow", case.inflow]]
    lak_kwargs = {}
    if case.table:
        stages = np.arange(7.0, 12.0 + 1e-9, 0.25)
        rows = "\n".join(f"  {s:.4f} {AREA * (s - 7.0):.4f} {AREA:.4f}" for s in stages)
        (ws / "lake.tab").write_text(
            f"BEGIN dimensions\n  nrow {len(stages)}\n  ncol 3\nEND dimensions\n\n"
            f"BEGIN table\n{rows}\nEND table\n")
        lak_kwargs.update(ntables=1, tables=[[0, "lake.tab"]])
    flopy.mf6.ModflowGwflak(
        gwf,
        print_input=True, print_stage=True, print_flows=True, save_flows=True,
        surfdep=case.surfdep,
        stage_filerecord="lak.stg", budget_filerecord="lak.bud",
        nlakes=1, noutlets=0,
        packagedata=[[0, strt, len(case.conns)]],
        connectiondata=[[0, i, c.cellid, c.ctype, BEDLEAK, c.belev, c.telev, c.connlen, c.connwidth]
                        for i, c in enumerate(case.conns)],
        perioddata={0: period},
        **lak_kwargs,
    )
    dis.idomain.set_data(idomain)
    sim.write_simulation(silent=True)
    if case.implicit:  # MF6 6.8 dev option; FloPy 3.10 cannot write it (ledger 326)
        insert_implicit(ws / "gwf.lak")
    ok, _buff = sim.run_simulation(silent=True, report=True)

    out = {"ok": ok}
    if not ok:
        mfsim = (ws / "mfsim.lst").read_text() if (ws / "mfsim.lst").exists() else ""
        errs = [ln.strip() for ln in mfsim.splitlines()
                if re.match(r"\s*\d+\.\s", ln) and "Error occurred" not in ln]
        out["error"] = " | ".join(dict.fromkeys(errs))
        return out

    listing = (ws / "gwf.lst").read_text()
    # per-connection lake-frame flow (+ into lake), in connection order
    gwf_rec = flopy.utils.CellBudgetFile(ws / "lak.bud", precision="double").get_data(text="GWF")[-1]
    out["q"] = [float(v) for v in gwf_rec["q"]]
    out["node_used"] = [int(v) for v in gwf_rec["node2"]]  # 1-based GWF node MF6 exchanged with
    # aquifer-frame cross-check from the GWF model's own budget (+ into aquifer)
    lak_rec = flopy.utils.CellBudgetFile(ws / "gwf.cbc", precision="double").get_data(text="LAK")[-1]
    assert abs(sum(out["q"]) + float(lak_rec["q"].sum())) < 1e-6, case.name
    heads = flopy.utils.HeadFile(ws / "gwf.hds").get_data().ravel()
    out["h_used"] = [float(heads[n - 1]) for n in out["node_used"]]
    live = heads[np.abs(heads) < 1e29]
    out["hmax"] = float(live.max())
    stg = float(flopy.utils.HeadFile(ws / "lak.stg", text="STAGE", precision="double")
                .get_data().ravel()[0])
    out["stage_file_dry"] = stg <= DRY
    out["mf6_stage"], out["mf6_lakebot"] = _listing_numbers(listing)
    return out


def V(layer, b, t=None):
    """Vertical connection; BELEV = b, TELEV = t (defaults to b, as 'bathy' writes it)."""
    return Conn(layer, "vertical", b, b if t is None else t)


H_SIDE = Conn(0, "horizontal", 7.0, 10.0, connlen=5.0, connwidth=10.0)

CASES = [
    Case("a", [V(0, 7.0)], stage=8.0, note="bed mid-L1, stage above bed, below L1 top"),
    Case("b", [V(0, 7.0)], stage=12.0, note="stage above L1 top"),
    Case("c1", [V(0, 9.9)], stage=8.0, note="only the written bed changes"),
    Case("c2", [V(0, 9.9)], stage=12.0),
    Case("c3", [V(0, 0.5)], stage=8.0, note="written bed BELOW the L1 cell bottom"),
    Case("c4", [V(0, 0.5)], stage=12.0),
    Case("c5", [V(0, 7.0, 10.0)], stage=8.0, note="belev 7 / telev 10 written"),
    Case("d1", [V(1, 5.0)], stage=8.0, note="vertical to LAYER 2 (top 5), L1 cell active"),
    Case("d2", [V(1, 5.0)], stage=8.0, l1_center_active=False,
         note="same, L1 cell above made INACTIVE"),
    Case("d3", [V(1, 3.0)], stage=4.0, note="bed mid-L2, stage above bed, below L2 top"),
    Case("d4", [V(1, 3.0)], stage=8.0, note="same bed, stage above L2 top"),
    Case("e1", [V(0, 7.0)], stage=8.0, table=True, note="case a + TAB6 table 7..12"),
    Case("e2", [V(0, 7.0)], inflow=80.0, note="ACTIVE lake, INFLOW 80: steady stage"),
    Case("e3", [V(0, 7.0)], inflow=80.0, table=True, note="e2 + TAB6 table 7..12"),
    Case("f", [H_SIDE], stage=8.0, note="HORIZONTAL, same L1 cell"),
    Case("h", [V(0, 7.0), H_SIDE], stage=8.0,
         note="both in one lake (what 'bathy' writes at a step)"),
    Case("i1", [V(0, 7.0)], stage=8.0, chd_head=12.0,
         note="aquifer head ABOVE L1 top (CHD 12), stage 8"),
    Case("i2", [V(0, 7.0)], stage=12.0, chd_head=13.0,
         note="aquifer head (CHD 13) AND stage above L1 top"),
    Case("g1", [V(0, 7.0)], stage=8.0, newton=True, note="NEWTON"),
    Case("g2", [V(0, 7.0)], stage=12.0, newton=True, note="NEWTON"),
    Case("g3", [V(0, 7.0)], stage=8.0, newton=True, icelltype=1, plain_formula=False,
         note="icelltype=1 + NEWTON"),
    Case("g4", [V(0, 7.0)], stage=12.0, icelltype=1,
         note="icelltype=1, Picard: L1 dries, conn moves to L2"),
    Case("g5", [V(0, 7.0)], stage=12.0, surfdep=1.0, plain_formula=False,
         note="SURFDEP 1.0, Picard"),
    Case("g6", [V(0, 7.0)], stage=12.0, surfdep=1.0, newton=True, plain_formula=False,
         note="SURFDEP 1.0 + NEWTON"),
    Case("g7", [V(0, 7.0)], inflow=80.0, implicit=True, note="e2 with LAK IMPLICIT"),
    Case("g8", [V(0, 7.0)], stage=8.0, l1_center_active=False,
         note="vertical to an INACTIVE L1 cell"),
]


def part_one(root: Path) -> None:
    print("PART 1 -- MF6 alone (raw FloPy, 3 x 3 DIS of 10 ft cells, top 10, botm 5 / 0)\n")
    hdr = (f"{'case':4s} {'conn':4s} {'belev/telev':>11s} {'stage in':>8s} {'MF6 stg':>8s} "
           f"{'q meas':>9s} {'claim':>9s} {'honour':>9s} {'h used':>6s} {'node':>4s} "
           f"{'MF6 bot':>7s}  note")
    print(hdr)
    print("-" * len(hdr))
    for case in CASES:
        r = case.result = build_and_run(root, case)
        if not r["ok"]:
            c = case.conns[0]
            print(f"{case.name:4s} {c.tag:4s} {f'{c.belev:g}/{c.telev:g}':>11s} "
                  f"{case.stage:8.2f}  MF6 REFUSED THE INPUT: {r['error']}   [{case.note}]")
            continue
        for i, c in enumerate(case.conns):
            be = f"{c.belev:g}/{c.telev:g}"
            h = r["h_used"][i]
            if case.stage is None:  # ACTIVE lake: the steady stage is the answer
                sin = f"Q={case.inflow:g}"
                cl = f"s={predict_stage(c, case.inflow, False):.3f}"
                ho = f"s={predict_stage(c, case.inflow, True):.3f}"
            else:
                sin = f"{case.stage:.2f}"
                cl = f"{predict_q(c, case.stage, h, False):9.3f}" if case.plain_formula else "-"
                ho = f"{predict_q(c, case.stage, h, True):9.3f}" if case.plain_formula else "-"
            stg = f"{r['mf6_stage']:.3f}" + ("*" if r["stage_file_dry"] else "")
            name = case.name if i == 0 else ""
            note = case.note if i == 0 else ""
            print(f"{name:4s} {c.tag:4s} {be:>11s} {sin:>8s} {stg:>8s} {r['q'][i]:9.3f} "
                  f"{cl:>9s} {ho:>9s} {h:6.3f} {r['node_used'][i]:4d} {r['mf6_lakebot']:7g}  {note}")
    print(f"""
q meas  : LAK budget 'GWF' term per connection, ft3/d, + INTO the lake (asserted equal and
          opposite to the GWF model's own 'LAK' cell-budget term in every run).
claim   : MF6's exchange formula with a VERTICAL connection's belev forced to its cell TOP.
honour  : the same formula with belev = the BELEV written. '-' = saturation terms apply.
MF6 stg : the stage MF6 used, from its printed stage table. '*' = MF6 wrote DRY (-1e30) to
          the binary stage file, i.e. it considers the lake at or below its own bottom.
h used  : head in the GWF cell MF6 actually exchanged with ('node', 1-based: 5 = L1 centre,
          14 = L2 centre). MF6 bot: first row of MF6's STAGE/VOLUME RELATION table.
Max head in any active cell over the low-water-table runs (all but i1/i2):
{max(c.result['hmax'] for c in CASES if c.result['ok'] and c.chd_head == CHD_HEAD):.3f} ft
(CHD 1.0 on the L2 outer ring; every bed tested is at or above 3.0 ft, so every one of those
lakes is perched above the water table and should lose water wherever it is wet).""")


# ---- PART 2: the same bench through myflopy, un-excavated against excavated -------------
def square_grid() -> VoronoiGridPlus:
    """The nine 10 ft squares of PART 1 as a DISV grid; cell ``r * 3 + c``, centre cell 4."""

    n = NROW
    verts = np.array([[x * DELR, y * DELC] for y in range(n + 1) for x in range(n + 1)], float)
    iverts, centres = [], []
    for r in range(n):
        for c in range(n):
            v0 = r * (n + 1) + c
            iverts.append([v0, v0 + n + 1, v0 + n + 2, v0 + 1])  # CLOCKWISE, or MF6 refuses
            centres.append([(c + 0.5) * DELR, (r + 0.5) * DELC])
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=np.array(centres))


def run_part_two(root: Path, tag: str, stage: float, excavated: bool) -> dict:
    grid = square_grid()
    pond = gpd.GeoDataFrame({"name": ["pond"]}, geometry=[box(DELR, DELC, 2 * DELR, 2 * DELC)],
                            crs=grid.crs)
    stack = mf.LayerStack(top=mf.Flat(TOP)).add("upper", bottom=mf.Flat(BOTM[0])).add(
        "lower", bottom=mf.Flat(BOTM[1]))
    lake_kw = dict(nper=1, lakes=pond, lake_id_field="name", starting_stage=stage,
                   bed_leakance=BEDLEAK, status="CONSTANT", stage=stage, print_input=True,
                   print_stage=True)
    if excavated:
        result = stack.excavate(bottom=7.0, within=pond, name="pond").build(grid, attach=True)
        lake_kw.update(unexcavated="raise")  # the bed is the record's: no lake_bottom
    else:
        result = stack.build(grid, attach=True)
        # per CELL, as a real pond's bathymetry is: the land around it sits at the model
        # top, so "bathy" writes a wall at each of the four steps (a scalar has no steps)
        lake_kw.update(lake_bottom=pd.Series({4: 7.0}), unexcavated="ignore")
    context = mf.ModelContext(grid=grid, domain=result.idomain)
    ring = [[(1, cell), CHD_HEAD] for cell in range(grid.ncpl) if cell != 4]
    lake = mf.lak(context=context, **lake_kw)
    model = mf.gwf("m", context=context, packages=[
        mf.disv(result), mf.npf(k=K, k33=K33, icelltype=0), mf.ic(strt=CHD_HEAD),
        mf.chd(stress_period_data={0: ring}), lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
        mf.ims(models=["m"], complexity="MODERATE", outer_dvclose=1e-8, inner_dvclose=1e-10,
               outer_maximum=500, inner_maximum=300, linear_acceleration="BICGSTAB")])
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        # the icelltype-0 wall and no-NEWTON warnings: this bench is confined on purpose,
        # so MF6's plain formula applies and the prediction can be made by hand
        warnings.simplefilter("ignore")
        run = mf.Project(root / "part2" / tag, name="p").run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    ws = Path(run.workspace)
    listing = (ws / "m.lst").read_text()
    discrepancy = [float(line.split()[-1]) for line in listing.splitlines() if "PERCENT DISCREPANCY" in line]
    record = flopy.utils.CellBudgetFile(ws / "m_lak_budget.lak", precision="double").get_data(text="GWF")[-1]
    rows = lake.options["connectiondata"]
    mf6_stage, mf6_bottom = _listing_numbers(listing)
    heads = np.asarray(flopy.utils.HeadFile(ws / "m.hds").get_data()).reshape(NLAY, -1)
    return dict(rows=rows, q=[float(v) for v in record["q"]], discrepancy=discrepancy,
                mf6_stage=mf6_stage, mf6_bottom=mf6_bottom, rules=lake.metadata.get("lake_rules"),
                top=float(result.top[4]), heads=heads)


def part_two(root: Path) -> None:
    print("\nPART 2 -- the same bench through myflopy (DISV of the same nine cells, lake = the"
          " centre cell,\nbed 7, bedleak 1, CONSTANT stage): the old rule on the un-excavated"
          " stack (L) against\nLayerStack.excavate(bottom=7.0) + mf.lak with no lake_bottom (E)\n")
    print(f"{'run':3s} {'stage':>5s} {'rule':>9s} {'col top':>7s} {'MF6 bot':>7s} {'MF6 stg':>7s} "
          f"{'vertical q':>10s} {'predicted':>9s} {'walls q':>8s} {'walls on':>9s}  discrepancy")
    for stage in (8.0, 12.0):
        for excavated in (False, True):
            tag = f"{'E' if excavated else 'L'}{stage:g}"
            r = run_part_two(root, tag, stage, excavated)
            vertical = [q for row, q in zip(r["rows"], r["q"], strict=True) if row[3] == "VERTICAL"]
            walls = [(row, q) for row, q in zip(r["rows"], r["q"], strict=True) if row[3] == "HORIZONTAL"]
            on = sorted({row[2][1] for row, _ in walls})
            where = "lake cell" if on == [4] else f"cells {on}"
            (layer, cell) = next(row[2] for row in r["rows"] if row[3] == "VERTICAL")
            cell_top = r["top"] if layer == 0 else BOTM[0]
            thickness = cell_top - (BOTM[0] if layer == 0 else BOTM[1])
            head = float(r["heads"][layer, cell])
            conn = Conn(layer, "vertical", cell_top, cell_top)
            predicted = satcond(conn, thickness) * (max(head, cell_top) - max(stage, cell_top))
            print(f"{tag[0]:3s} {stage:5g} {r['rules']['pond']:>9s} {r['top']:7.3f} "
                  f"{r['mf6_bottom']:7.3f} {r['mf6_stage']:7.3f} {sum(vertical):10.3f} "
                  f"{predicted:9.3f} {sum(q for _, q in walls):8.3f} {where:>9s}  {r['discrepancy']}")
    print("""
rule      : metadata['lake_rules'] -- 'legacy' is the old rule, 'excavated' the record's.
col top   : the model top of the lake's column in the DISV (the ground, or the cut bed).
MF6 bot   : MF6's own lake bottom (first row of its STAGE/VOLUME RELATION table).
predicted : MF6's formula C (max(h, belev) - max(stage, belev)) with belev the vertical
            cell's TOP, h its head, and C = area / (1/bedleak + 0.5 thickness / K33).
walls q   : the HORIZONTAL connections together; 'walls on' names the cells they are on.
discrepancy: the GWF model's percent discrepancy, as MF6 reports it.""")


def main(root: Path) -> None:
    if shutil.which("mf6") is None:
        raise SystemExit("mf6 is not on PATH")
    part_one(root)
    part_two(root)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(Path(sys.argv[1]))
    else:
        with tempfile.TemporaryDirectory() as tmp:
            main(Path(tmp))
