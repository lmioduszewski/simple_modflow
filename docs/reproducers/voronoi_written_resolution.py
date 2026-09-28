"""A corner split closer than the DISV file records: a face written with zero length, NaN heads.

Ledger 316, 2026-09-28.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/voronoi_written_resolution.py [run_dir ...]

WHAT HAPPENED. A 20,323-cell ``GridSpec.voronoi`` model at state plane (EPSG:2927,
x ~ 1.26e6 ft) built, ran to "Normal termination" and wrote NaN for every head: MF6's
first outer iteration already had a NaN residual. The repair (from bb56668) had rewritten
9 cells and appended 3 vertices. One appended vertex sat 7.4e-5 ft from one of FloPy's,
and the two were the ends of a face between two cells. FloPy writes a coordinate above
1e5 to nine significant digits -- 0.01 ft at x = 1.26e6 -- so both printed as
``1.25858105E+06 7.19399294E+05`` and MF6 read a face of zero length. Rewriting that ONE
vertex line at full precision, nothing else changed, the same model solved in 32 s with a
0.00% budget discrepancy and every head finite (91-703 ft).

Why the two vertices were that close: four of Triangle's nodes (8692, 8694, 8695, 8878)
sat almost on one circle -- a near-rectangle about 10 x 15 ft, the fourth node 5.7e-5 ft
off the circle through the other three (radius 7.9 ft), inside a 100 ft2 refinement the
smoothing loop had made near-regular. The two triangles on their shared edge then have
circumcentres 7.4e-5 ft apart: a true Voronoi face that short. FloPy (Qhull, in raw
state-plane coordinates) drew ONE corner there, so its cells had no face and ran. The
repair drew both, as the exact dual. Of 60,646 faces in that mesh's dual, 8 are shorter
than 0.01 ft; a second pair the repair split, 3.5e-4 ft either side of FloPy's corner,
happened to print one digit apart and ran -- by the luck of the rounding.

The rule now: dual vertices joined by a face shorter than one written unit on both axes
are ONE corner -- FloPy's where it drew one -- so the repair never writes a face the file
cannot hold. On that model it rewrites 2 cells (the reflex-corner fix it exists for),
appends nothing, and keeps 41,813 vertices, FloPy's count. A model whose DISV would still
be written with a zero-length face -- a grid rebuilt from a snapshot is the grid that ran
-- is refused at build, naming the cells.

WHAT THIS RUNS. The same geometry reduced to four quads around one corner at
(1.2e6 + 100, 6.5e5 + 100), the corner split in two 7.4e-5 ft apart with a face between
cells 0 and 2, the two MF6 solves for (constant heads in 1 and 3; between two
constant-head cells the same face would reach only the budget):
  1. MF6 on the split corner at FloPy's default number format: NaN heads, all four.
  2. The same split written at ``float_precision = 16``: 5, 10, 5, 0 -- it is the FILE.
  3. The corner joined, as the repair now draws it: 5, 10, 5, 0.
  4. What the repair does with that corner now, and what it did before this rule.
  5. FloPy's OWN cells, no repair involved: a plain four-corner boundary at (1.2e6, 6.5e5),
     max_area 2000, default profile, no points. FloPy drew a face 3.5e-4 ft long whose ends
     print as one point (NaN heads in MF6); the repair joins it and rewrites 3 cells, and the
     same boundary on the "clean" profile, which has no such face, is untouched.
Pass run directories to check their written ``.disv`` files for zero-length faces.
"""

from __future__ import annotations

import pathlib
import sys
import tempfile
import warnings

import flopy
import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon

import myflopy as mf
from myflopy.grid_identity import grid_pin, zero_length_faces
from myflopy.modflow.mf6.grid.tessellation import _onto_flopy

X0, Y0 = 1.2e6 + 100.0, 6.5e5 + 100.0
QUADS = np.array([(X0 - 100, Y0 - 100), (X0, Y0 - 100), (X0 + 100, Y0 - 100), (X0 + 100, Y0),
                  (X0 + 100, Y0 + 100), (X0, Y0 + 100), (X0 - 100, Y0 + 100), (X0 - 100, Y0), (X0, Y0)])
JOINED = [[0, 7, 8, 1], [1, 8, 3, 2], [8, 5, 4, 3], [7, 6, 5, 8]]
SPLIT = [[0, 7, 8, 9, 1], [1, 9, 3, 2], [8, 5, 4, 3, 9], [7, 6, 5, 8]]
CENTRES = [(X0 - 50, Y0 - 50), (X0 + 50, Y0 - 50), (X0 + 50, Y0 + 50), (X0 - 50, Y0 + 50)]
GAP = 7.4e-5


def run(verts, rings, precision=None):
    """MF6 on four quads, CHD 10 ft in cell 1 and 0 ft in cell 3 -- the split face is between
    cells 0 and 2, which MF6 solves for: ``(heads, zero faces)``."""
    with tempfile.TemporaryDirectory() as ws:
        sim = flopy.mf6.MFSimulation(sim_ws=ws, exe_name="mf6")
        if precision is not None:
            sim.simulation_data.float_precision = precision
            sim.simulation_data.float_characters = precision + 8
        flopy.mf6.ModflowTdis(sim)
        flopy.mf6.ModflowIms(sim, complexity="simple")
        gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
        disv = flopy.mf6.ModflowGwfdisv(
            gwf, nlay=1, ncpl=4, nvert=len(verts), top=10.0, botm=[0.0],
            vertices=[[i, x, y] for i, (x, y) in enumerate(verts)],
            cell2d=[[k, *CENTRES[k], len(r), *r] for k, r in enumerate(rings)])
        flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
        flopy.mf6.ModflowGwfic(gwf, strt=5.0)
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 1), 10.0], [(0, 3), 0.0]])
        flopy.mf6.ModflowGwfoc(gwf, head_filerecord="m.hds", saverecord=[("HEAD", "ALL")])
        faces = zero_length_faces(disv, sim.simulation_data)
        sim.write_simulation(silent=True)
        ok, _ = sim.run_simulation(silent=True)
        heads = gwf.output.head().get_data().ravel() if ok else None
    return heads, faces


def check_run(folder: pathlib.Path) -> None:
    """Zero-length faces in a written ``.disv``, from the file's own text."""
    for path in sorted(folder.glob("*.disv")):
        lines = path.read_text().splitlines()

        def block(name):
            s = next(i for i, line in enumerate(lines) if line.strip().upper() == f"BEGIN {name}")
            e = next(i for i, line in enumerate(lines) if line.strip().upper() == f"END {name}")
            return [line.split() for line in lines[s + 1:e]]

        written = [tuple(t[1:3]) for t in block("VERTICES")]
        owners: dict = {}
        for t in block("CELL2D"):
            ring = [int(v) - 1 for v in t[4:4 + int(t[3])]]
            for a, b in zip(ring, ring[1:] + ring[:1], strict=True):
                if written[a] == written[b]:
                    owners.setdefault(frozenset((a, b)), []).append(int(t[0]) - 1)
        shared = [(cells, written[min(e)]) for e, cells in owners.items() if len(cells) > 1]
        print(f"  {path}: {len(shared)} face(s) written with zero length"
              + "".join(f"\n    cells {c} at {xy}" for c, xy in shared))


def main() -> None:
    warnings.simplefilter("ignore")
    split = np.vstack([QUADS, [(X0 + GAP, Y0)]])

    heads, faces = run(split, SPLIT)
    print(f"1. split corner, FloPy's default format: zero-length faces {faces}; heads {heads}")
    heads, faces = run(split, SPLIT, precision=16)
    print(f"2. the same split at float_precision=16: zero-length faces {faces}; heads {np.round(heads, 4)}")
    heads, faces = run(QUADS, JOINED)
    print(f"3. the corner joined: zero-length faces {faces}; heads {np.round(heads, 4)}")

    verts, iverts, rewritten, appended, unsplit, was = _onto_flopy(split, SPLIT, QUADS, JOINED, np.zeros(4, bool))
    print(f"4. the repair on that corner: {unsplit} vertex joined, {appended} appended, {len(rewritten)} cells "
          f"rewritten, rings {'FloPy' if [list(r) for r in iverts] == JOINED else 'rewritten'}'s; before this "
          f"rule it appended {len(was[0]) - len(QUADS)} and wrote {len(was[0])} vertices")

    corners = [(0, 0), (1100.3, 130.7), (1250.9, 940.2), (310.4, 1180.6)]
    quad = Polygon([(1.2e6 + x, 6.5e5 + y) for x, y in corners])
    with tempfile.TemporaryDirectory() as root:
        root = pathlib.Path(root)
        gpd.GeoDataFrame({"n": ["b"]}, geometry=[quad], crs=2927).to_file(root / "b.gpkg")
        for profile in ("balanced", "clean"):
            for tessellation in ("flopy", "dual"):
                spec = mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec(path=root / "b.gpkg"),
                                           max_area=2000.0, crs="EPSG:2927", profile=profile,
                                           tessellation=tessellation)
                vor = spec.resolve(workspace=root / f"{profile}-{tessellation}")
                props = vor.get_disv_gridprops()
                sim = flopy.mf6.MFSimulation(sim_ws=str(root / "sim"))
                gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
                disv = flopy.mf6.ModflowGwfdisv(gwf, nlay=1, ncpl=props["ncpl"], nvert=props["nvert"], top=1.0,
                                                botm=[0.0], vertices=props["vertices"], cell2d=props["cell2d"])
                report = vor.tessellation_report
                print(f"5. {profile:8s} {tessellation:5s}: {vor.ncpl} cells, pin {grid_pin(vor)}, zero-length faces "
                      f"{zero_length_faces(disv, sim.simulation_data)}"
                      + (f"; {report.summary()}" if report is not None else ""))

    for folder in sys.argv[1:]:
        check_run(pathlib.Path(folder))


if __name__ == "__main__":
    main()
