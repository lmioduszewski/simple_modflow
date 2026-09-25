"""Array-form recharge written into a pass-through cell, which MF6 drops without a word.

Ledger 319. Found while scoping `LayerStack.excavate` (its D4): `mf.rch.array` and
`mf.evt.array` resolved `irch`/`ievt="top_active"` with `idomain != 0`, so a
vertical pass-through cell (`idomain == -1`) counted as the column's top ACTIVE cell.
MF6 has removed that cell, and a flux addressed to it is skipped -- the column gets
nothing, and the run ends in Normal termination.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/top_active_passthrough.py

Needs `mf6` on PATH (measured with 6.8.0.dev0). About 3 s.

THE SETUP. A real triangulated grid -- `GridSpec.voronoi`, `max_area=12000` -- on a
600 x 600 ft square: 92 cells. A `LayerStack` with a flat top at 100 ft and three
units: "fill" to 95 ft, except west of x = 250 ft where its bottom IS the top, so it
pinches out there under `pinch="passthrough"`; "sand" to 80 ft; "base" to 50 ft. That
is 40 pass-through cells in layer 1 (MF6's numbering), every one the top of its
column. Steady state, CHD at 90 ft on the east edge, NPF k = 10 ft/d, `icelltype` 1,
array recharge 0.001 ft/d everywhere.

Four runs on it:

  OLD    `irch` resolved by the old rule, applied here by hand to the stack's own
         idomain -- so this script still shows the "before" once the library is fixed.
  LAZY   `mf.rch.array(recharge=...)` with no context: `irch="top_active"` resolves
         at build against the built discretization (ledger 285).
  EAGER  the same with `context=ModelContext(grid=, domain=)`, resolved at
         declaration. Its DIGEST is what moves; the lazy recipe holds the string.
  FIXED  LAZY with `fixed_cell=True`: MF6 keeps recharge in the IRCH cell even when
         that cell is DRY. The steady head sits below the fill's 95 ft bottom, so
         every ACTIVE fill cell -- the top of each of the other 52 columns -- is dry.

WHAT IT MEASURES, per run: the recharge MF6 applied in total and in the pass-through
columns, against the rate times the area, and MF6's own percent discrepancy. Then how
many active fill cells are dry and what falls on them, the eager declaration's digest
under the old rule and the new one.

What `fixed_cell=True` loses is the DRY cells' share, 207.98 of 360.00 ft3/d; it
applies 152.02, all of it in the pass-through columns, whose IRCH cell is the wet sand
below. (Measured first against the OLD `irch`, it read 207.98 -> 0.0: that pair is the
old rule, which also dropped those 152.02, and it reached the docstring and ledger 319
before a review caught it -- which is why the case is run here.)
"""

from __future__ import annotations

import contextlib
import io
import pathlib
import tempfile
import warnings

import flopy
import geopandas as gpd
import numpy as np
from shapely.geometry import box

import myflopy as mf
from myflopy.specs import ModelContext, package_digest

RATE = 0.001
SIDE = 600.0


def _grid(root: pathlib.Path):
    boundary = root / "boundary.gpkg"
    gpd.GeoDataFrame({"n": ["d"]}, geometry=[box(0, 0, SIDE, SIDE)], crs=2927).to_file(boundary)
    return mf.GridSpec.voronoi(
        boundary=mf.GeoPackageSourceSpec(path=boundary), max_area=12_000.0
    ).resolve(workspace=root / "grid")


def _stack():
    return (
        mf.LayerStack(top=mf.Surface.flat(100.0))
        .add("fill", bottom=mf.Surface.where(box(0, 0, 250, SIDE), 100.0, 95.0),
             pinch="passthrough")
        .add("sand", bottom=80.0)
        .add("base", bottom=50.0)
    )


def _old_rule(idomain):
    """What `_array_bc_values` computed before ledger 319: `!= 0` is "active"."""

    active = np.asarray(idomain) != 0
    return np.where(active.any(axis=0), np.argmax(active, axis=0), 0).astype(int)


def _run(root, vor, tag, recharge):
    east = [int(c) for c in np.flatnonzero(vor.gdf_vorPolys.geometry.bounds.maxx >= SIDE - 1e-6)]
    model = mf.gwf("m", context=ModelContext(grid=vor), packages=[
        mf.disv(_stack(), attach=False),
        mf.npf(k=10.0, icelltype=1),
        mf.ic(strt=95.0),
        mf.chd(stress_period_data={0: [[(2, c), 90.0] for c in east]}),
        recharge,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]), mf.ims(models=["m"])])
    # `Run.execute(silent=True)` still lets FloPy's "writing package ..." lines through
    # (its `write()` is called without `silent`), so they are caught here.
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        run = mf.Project(root / f"p_{tag}", name=tag).run(tag, simulation)
    listing = (pathlib.Path(run.workspace) / "m.lst").read_text()
    discrepancy = [line for line in listing.splitlines() if "PERCENT DISCREPANCY" in line][-1]
    budget = flopy.utils.CellBudgetFile(str(pathlib.Path(run.workspace) / "m.cbc"),
                                        precision="double")
    record = budget.get_data(text="RCH")[-1]
    head = flopy.utils.HeadFile(str(pathlib.Path(run.workspace) / "m.hds")).get_data()
    return run.success, record, float(discrepancy.split()[-1]), np.asarray(head)


def main() -> None:
    root = pathlib.Path(tempfile.mkdtemp(prefix="top_active_"))
    vor = _grid(root)
    ncpl = int(vor.ncpl)
    area = np.asarray(vor.gdf_vorPolys.area, dtype=float)
    idomain = _stack().build(vor, attach=False).idomain
    passthrough = np.flatnonzero(idomain[0] == -1)
    print(f"{ncpl} columns, {len(passthrough)} with a pass-through cell on top "
          f"(idomain values {sorted(set(np.unique(idomain).tolist()))})")
    print(f"expected recharge: {RATE * area.sum():.2f} ft3/d in total, "
          f"{RATE * area[passthrough].sum():.2f} in the pass-through columns\n")

    context = ModelContext(grid=vor, domain=idomain)
    old = _old_rule(idomain)
    cases = {
        "OLD": mf.rch.array(recharge=RATE, irch=old),
        "LAZY": mf.rch.array(recharge=RATE),
        "EAGER": mf.rch.array(recharge=RATE, context=context),
        "FIXED": mf.rch.array(recharge=RATE, fixed_cell=True),
    }
    heads = {}
    for tag, spec in cases.items():
        success, record, discrepancy, heads[tag] = _run(root, vor, tag, spec)
        column = (np.asarray(record["node"]) - 1) % ncpl
        applied = np.asarray(record["q"], dtype=float)
        into = applied[np.isin(column, passthrough)].sum()
        print(f"{tag:6s} normal termination={success}  applied {applied.sum():8.2f} "
              f"({100 * applied.sum() / (RATE * area.sum()):5.1f}%)  in the pass-through "
              f"columns {into:7.2f}  percent discrepancy {discrepancy:+.2f}")

    # MF6 writes HDRY (-1e30) for a dry convertible cell; `idomain > 0` is active.
    fill = idomain[0] > 0
    dry = fill & (heads["LAZY"].reshape(idomain.shape)[0] <= -1e29)
    print(f"\nfill layer: {int(fill.sum())} active cells, {int(dry.sum())} dry "
          f"(LAZY run); recharge falling on the dry ones {RATE * area[dry].sum():.2f} ft3/d "
          "-- passed to the cell below by default, lost with fixed_cell=True")

    eager_old = mf.rch.array(recharge=RATE, irch=old, context=context)
    print(f"\neager irch in the pass-through columns: old rule layer "
          f"{sorted(set(old[passthrough].tolist()))}, new rule layer "
          f"{sorted(set(np.asarray(cases['EAGER'].options['irch'])[passthrough].tolist()))}"
          " (zero-based)")
    print("eager digest moved:", package_digest(eager_old) != package_digest(cases["EAGER"]))
    print("lazy recipe inputs:", cases["LAZY"].inputs, "(the digest holds the string)")


if __name__ == "__main__":
    main()
