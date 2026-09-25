"""What MF6 requires of LAK ``IMPLICIT``, and what FloPy does with it (ledger 326).

    /home/lukem/python/envs/gw/bin/python docs/reproducers/lak_implicit_requirements.py

Needs an ``mf6`` on PATH that has the option -- MF6's develop branch (commit
``e2dfc84c``, PR #2860) or the 6.8 line; measured with 6.8.0.dev0+0418efe4. About 5 s.
Each MF6 claim is shown by MF6 itself: the verdict is ``Normal termination`` in
``mfsim.lst`` or the ERROR lines MF6 printed, never myflopy's refusal, so the premise of
every refusal in ``grid_identity.implicit_solver_complaints`` and ``LAKBuilder`` can be
checked against the executable in hand.

PART 1 -- the linear acceleration. A two-cell model (one lake on cell 1, a CHD on cell
2) written as plain MF6 text, its IMS varied. Measured: with ``IMPLICIT``, BICGSTAB runs;
``LINEAR_ACCELERATION CG`` STOPS ("... PRODUCES AN ASYMMETRIC COEFFICIENT MATRIX, BUT THE
CONJUGATE GRADIENT METHOD WAS SELECTED. USE BICGSTAB INSTEAD."), and so do COMPLEXITY
SIMPLE with no LINEAR_ACCELERATION and an IMS naming neither -- MF6's default there is
CG; COMPLEXITY MODERATE and COMPLEX with none run, their default being BICGSTAB. Control:
the same CG IMS runs without ``IMPLICIT``.

PART 2 -- a lake with NO connection. Two lakes, the second with ``nlakeconn`` 0.
Measured: with ``IMPLICIT`` MF6 stops ("lake 2 has no connections; the IMPLICIT option
requires each lake to have at least one GWF connection"); without it, Normal
termination. ``mf.lak(implicit=True)`` refuses such a lake at declaration.

PART 3 -- FloPy's loader. ``flopy.mf6.MFSimulation.load`` of PART 1's BICGSTAB model and
a rewrite: the rewritten LAK file carries no ``IMPLICIT`` and nothing warned (measured).
``lak_implicit.keep_implicit`` -- what ``Run.open`` and ``mf.load_mf6_run`` call --
re-classes the loaded package, and a rewrite then carries it once.

PART 4 -- the converged answer. A pit excavated into a 12 x 12 Voronoi grid of 20 ft
cells (``LayerStack`` 100 / 82 / 40 ft; ``stack.excavate(bottom=82, within=<30 ft
circle>)``), NPF k = k33 = 10, icelltype 1, NEWTON; CHD 94 ft on the west and east
edges; three 1-day transient periods, the pit filled at 2,000 ft3/d from 95 ft; IMS
COMPLEX, BICGSTAB, outer_dvclose 1e-8, inner_dvclose 1e-10. Built through ``mf.lak`` with
and without ``implicit=True``. Measured: both Normal termination; the implicit listing
says the stage is solved in the matrix ("IMPLICIT FORMULATION") and the default one does
not; stages agree to 1.0e-8 ft; GWF and LAK percent discrepancy 0.00 in every period.
"""

from __future__ import annotations

import contextlib
import io
import pathlib
import shutil
import subprocess
import tempfile
import warnings

import flopy
import geopandas as gpd
import numpy as np
from shapely.geometry import Point

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.lak_implicit import keep_implicit

MFSIM = (
    "BEGIN options\nEND options\nBEGIN timing\n  TDIS6 p.tdis\nEND timing\n"
    "BEGIN models\n  gwf6 p.nam p\nEND models\nBEGIN exchanges\nEND exchanges\n"
    "BEGIN solutiongroup 1\n  ims6 p.ims p\nEND solutiongroup\n"
)
FILES = {
    "p.tdis": "BEGIN options\nEND options\nBEGIN dimensions\n  NPER 1\nEND dimensions\n"
              "BEGIN perioddata\n  1.0 1 1.0\nEND perioddata\n",
    "p.nam": "BEGIN options\nEND options\nBEGIN packages\n  DIS6 p.dis\n  NPF6 p.npf\n"
             "  IC6 p.ic\n  CHD6 p.chd\n  LAK6 p.lak\nEND packages\n",
    "p.dis": "BEGIN options\nEND options\nBEGIN dimensions\n  NLAY 1\n  NROW 1\n  NCOL 2\n"
             "END dimensions\nBEGIN griddata\n  delr\n    CONSTANT 10.0\n  delc\n    CONSTANT 10.0\n"
             "  top\n    CONSTANT 10.0\n  botm\n    CONSTANT 0.0\nEND griddata\n",
    "p.npf": "BEGIN options\nEND options\nBEGIN griddata\n  icelltype\n    CONSTANT 0\n"
             "  k\n    CONSTANT 1.0\nEND griddata\n",
    "p.ic": "BEGIN options\nEND options\nBEGIN griddata\n  strt\n    CONSTANT 9.0\nEND griddata\n",
    "p.chd": "BEGIN options\nEND options\nBEGIN dimensions\n  MAXBOUND 1\nEND dimensions\n"
             "BEGIN period 1\n  1 1 2 9.0\nEND period\n",
}
ONE_LAKE = ("BEGIN dimensions\n  NLAKES 1\n  NOUTLETS 0\n  NTABLES 0\nEND dimensions\n"
            "BEGIN packagedata\n  1 9.5 1\nEND packagedata\n"
            "BEGIN connectiondata\n  1 1 1 1 1 VERTICAL 1.0 0.0 0.0 0.0 0.0\nEND connectiondata\n")
TWO_LAKES = ("BEGIN dimensions\n  NLAKES 2\n  NOUTLETS 0\n  NTABLES 0\nEND dimensions\n"
             "BEGIN packagedata\n  1 9.5 1\n  2 9.5 0\nEND packagedata\n"
             "BEGIN connectiondata\n  1 1 1 1 1 VERTICAL 1.0 0.0 0.0 0.0 0.0\nEND connectiondata\n")


def _ims(linear: str) -> str:
    """An IMS file; ``linear`` is ``"<OPTIONS lines>|<LINEAR lines>"``."""

    return (f"BEGIN options\n{linear.split('|')[0]}END options\nBEGIN nonlinear\n  OUTER_MAXIMUM 50\n"
            f"END nonlinear\nBEGIN linear\n  INNER_MAXIMUM 50\n{linear.split('|')[1]}END linear\n")


def _write(root: pathlib.Path, *, ims: str, implicit: bool, lakes: str = ONE_LAKE) -> pathlib.Path:
    root.mkdir(parents=True)
    (root / "mfsim.nam").write_text(MFSIM)
    for name, text in FILES.items():
        (root / name).write_text(text)
    (root / "p.ims").write_text(_ims(ims))
    (root / "p.lak").write_text(
        ("BEGIN options\n  IMPLICIT\nEND options\n" if implicit else "BEGIN options\nEND options\n")
        + lakes
    )
    return root


def _verdict(root: pathlib.Path) -> str:
    subprocess.run(["mf6"], cwd=root, capture_output=True, text=True, timeout=120)
    listing = (root / "mfsim.lst").read_text() if (root / "mfsim.lst").exists() else ""
    if "Normal termination" in listing:
        return "Normal termination"
    # MF6's own words: the numbered items of its ERROR REPORT, each joined onto one line.
    report = listing.split("\nERROR REPORT:", 1)[-1].split("UNIT ERROR REPORT:", 1)[0]
    return "STOPPED: " + " ".join(report.split())


def part1(root: pathlib.Path) -> None:
    print("PART 1 -- the linear acceleration (two-cell model, IMPLICIT unless said)")
    cases = {
        "LINEAR_ACCELERATION BICGSTAB": "|  LINEAR_ACCELERATION BICGSTAB\n",
        "LINEAR_ACCELERATION CG": "|  LINEAR_ACCELERATION CG\n",
        "COMPLEXITY SIMPLE, no LINEAR_ACCELERATION": "  COMPLEXITY SIMPLE\n|",
        "neither option": "|",
        "COMPLEXITY MODERATE, no LINEAR_ACCELERATION": "  COMPLEXITY MODERATE\n|",
        "COMPLEXITY COMPLEX, no LINEAR_ACCELERATION": "  COMPLEXITY COMPLEX\n|",
    }
    for index, (label, ims) in enumerate(cases.items()):
        print(f"  {label:48s} {_verdict(_write(root / f'p1_{index}', ims=ims, implicit=True))}")
    control = _write(root / "p1_control", ims="|  LINEAR_ACCELERATION CG\n", implicit=False)
    print(f"  {'LINEAR_ACCELERATION CG, WITHOUT IMPLICIT':48s} {_verdict(control)}")


def part2(root: pathlib.Path) -> None:
    print("PART 2 -- a lake with no connection (two lakes, the second nlakeconn 0)")
    for implicit in (True, False):
        workspace = _write(root / f"p2_{implicit}", ims="|  LINEAR_ACCELERATION BICGSTAB\n",
                           implicit=implicit, lakes=TWO_LAKES)
        print(f"  {'IMPLICIT' if implicit else 'default formulation':48s} {_verdict(workspace)}")


def part3(root: pathlib.Path) -> None:
    print("PART 3 -- FloPy's loader")
    source = _write(root / "p3", ims="|  LINEAR_ACCELERATION BICGSTAB\n", implicit=True)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        simulation = flopy.mf6.MFSimulation.load(sim_ws=str(source), verbosity_level=0)
    lak = simulation.get_model("p").get_package("lak")
    print(f"  loaded class: {type(lak).__name__}; warnings while loading: {len(caught)}")
    simulation.set_sim_path(str(root / "p3_rewritten"))
    simulation.write_simulation(silent=True)
    text = (root / "p3_rewritten" / "p.lak").read_text().upper()
    print(f"  a FloPy rewrite carries IMPLICIT: {'IMPLICIT' in text}")
    simulation = flopy.mf6.MFSimulation.load(sim_ws=str(source), verbosity_level=0)
    print(f"  keep_implicit re-classed: {keep_implicit(simulation)}")
    simulation.set_sim_path(str(root / "p3_kept"))
    simulation.write_simulation(silent=True)
    text = (root / "p3_kept" / "p.lak").read_text().upper()
    print(f"  then a rewrite carries IMPLICIT {text.count('IMPLICIT')} time(s); "
          f"MF6 on it: {_verdict(root / 'p3_kept')}")


def _pit_simulation(grid, pond, result, *, implicit: bool):
    context = mf.ModelContext(grid=grid, domain=result.idomain, surfaces=result)
    xy = np.array([(c[1], c[2]) for c in grid.get_disv_gridprops()["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 20) | (xy[:, 0] > xmax - 20))
    lake = mf.lak(context=context, nper=3, lakes=pond, lake_id_field="name",
                  starting_stage=95.0, bed_leakance=3000.0, implicit=implicit,
                  unexcavated="raise", inflow=[2000.0, 2000.0, 2000.0])
    model = mf.gwf("m", context=context, newtonoptions="NEWTON", packages=[
        mf.disv(result), mf.npf(k=10.0, k33=10.0, icelltype=1), mf.ic(strt=95.0),
        mf.chd(stress_period_data={0: [[(1, int(c)), 94.0] for c in edge]}),
        mf.sto(ss=1e-5, sy=0.2, iconvert=1, transient={0: True}), lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    return mf.simulation(model, name="s", packages=[
        mf.tdis(nper=3, perioddata=[(1.0, 1, 1.0)] * 3),
        mf.ims(models=["m"], complexity="COMPLEX", outer_dvclose=1e-8, inner_dvclose=1e-10,
               linear_acceleration="BICGSTAB"),
    ])


def part4(root: pathlib.Path) -> None:
    print("PART 4 -- the converged answer (an excavated pit filled at 2,000 ft3/d, 3 days)")
    grid = irregular_voronoi_grid(nrow=12, ncol=12, cell_size=20.0, refined=False)
    xmin, ymin, xmax, ymax = grid.gdf_vorPolys.total_bounds
    pond = gpd.GeoDataFrame({"name": ["pit"]}, crs=grid.crs, geometry=[
        Point((xmin + xmax) / 2, (ymin + ymax) / 2).buffer(30.0, 32)])
    result = (mf.LayerStack(top=mf.Flat(100.0)).add("upper", bottom=mf.Flat(82.0))
              .add("lower", bottom=mf.Flat(40.0))
              .excavate(bottom=82.0, within=pond, name="pit").build(grid))
    stages = {}
    for implicit in (False, True):
        with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
            warnings.simplefilter("ignore")
            run = mf.Project(root / f"p4_{implicit}", name="p").run(
                "r", _pit_simulation(grid, pond, result, implicit=implicit))
        workspace = pathlib.Path(run.workspace)
        listing = (workspace / "m.lst").read_text()
        stage = flopy.utils.HeadFile(str(workspace / "m_lak_stage.lak"), text="STAGE",
                                     precision="double").get_alldata().ravel()
        discrepancy = [float(line.split()[-1]) for line in listing.splitlines()
                       if "PERCENT DISCREPANCY" in line]
        stages[implicit] = stage
        print(f"  {'IMPLICIT' if implicit else 'default formulation':20s} success={run.success} "
              f"'IMPLICIT FORMULATION' in listing={'IMPLICIT FORMULATION' in listing} "
              f"stages={np.round(stage, 6).tolist()} "
              f"max |discrepancy| (GWF and LAK, every period)={max(map(abs, discrepancy)):.2f}%")
    print(f"  max |stage difference| = {np.max(np.abs(stages[True] - stages[False])):.1e} ft")


def main() -> None:
    if shutil.which("mf6") is None:
        raise SystemExit("mf6 is not on PATH")
    print(subprocess.run(["mf6", "--version"], capture_output=True, text=True).stdout.strip())
    with tempfile.TemporaryDirectory(prefix="lak_implicit_") as scratch:
        root = pathlib.Path(scratch)
        part1(root)
        part2(root)
        part3(root)
        part4(root)


if __name__ == "__main__":
    main()
