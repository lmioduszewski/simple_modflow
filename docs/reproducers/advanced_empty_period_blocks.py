"""The measurement behind ledger 468: an empty SFR or LAK period block changes nothing.

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/advanced_empty_period_blocks.py

``SFRBuilder.perioddata`` and ``LAKBuilder.perioddata`` gave EVERY period an entry,
``[]`` where a period sets nothing, and FloPy writes each as an empty ``BEGIN period``
block. FloPy 3.10's ``MFTransientList`` reads the FIRST period's records for its dtype,
so a package whose period 0 sets nothing made ``gwf.plot()`` and the package's own
``.plot()`` raise ``AttributeError: 'list' object has no attribute 'dtype'``.

THE SETUP. A 4 x 4 grid of 25 ft cells, one convertible layer (top 10, bottom 0, K 1),
CHD at 9.0 and 7.0 in two corners, three transient periods of two steps:

  SFR  one stream across row 1, 50 ft3/d of INFLOW set in period 1 ONLY (or in no
       period), so periods 0 and 2 -- or all three -- set nothing.
  LAK  a 2 x 2-cell pond on the same grid in two layers (bottoms 5 and 0), bed 6 ft,
       RAINFALL set in period 1 only.

Each model is run twice: as the builder writes it now (no empty block) and with an
empty block put back for every period that sets nothing. MEASURED: every output file
-- the model budget, the package budget, heads, the package stage file -- is
BYTE-IDENTICAL between the two, and FloPy's plot raises only with the empty blocks.
MF6 keeps an advanced package's period settings until a block changes them, so an
empty block is a no-op. (An empty MVR block is NOT: it turns the movers off. That
builder is not touched.)
"""

from __future__ import annotations

import dataclasses
import filecmp
import itertools
import pathlib
import shutil
import sys
import tempfile
import warnings

import geopandas as gpd
import matplotlib
import numpy as np
from shapely.geometry import LineString, box

matplotlib.use("Agg")
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))

import myflopy as mf  # noqa: E402
from myflopy.specs import ModelContext  # noqa: E402
from tests.test_gpkg_on_cell_faces import _context  # noqa: E402


def _with_empty_blocks(spec, nper):
    data = spec.options["perioddata"]
    full = {period: list(data.get(period, [])) for period in range(nper)}
    return dataclasses.replace(spec, options={**spec.options, "perioddata": full})


def _sfr_model(empty, inflow):
    context = _context()
    disv = mf.disv(nlay=1, top=10.0, botm=[0.0], **context.grid.get_disv_gridprops())
    streams = gpd.GeoDataFrame(
        {"name": ["s0"]}, geometry=[LineString([(0, 37.5), (100, 37.5)])], crs="EPSG:2927"
    )
    sfr = mf.sfr(
        context=context, nper=3, streams=streams, stream_id="name", reach_top=8.0, width=5.0,
        inflow=inflow, save_flows=True, budget_filerecord="m.sfr.cbc",
        stage_filerecord="m.sfr.stg",
    )
    if empty:
        sfr = _with_empty_blocks(sfr, 3)
    return [disv, mf.npf(k=1.0, icelltype=1), mf.ic(strt=9.0),
            mf.sto(transient={0: True}, sy=0.2, ss=1e-5, iconvert=1),
            mf.chd(stress_period_data={0: [[(0, 0), 9.0], [(0, 15), 7.0]]}), sfr], "sfr"


def _lak_model(empty, _inflow):
    context = _context()
    grid = context.grid
    grid.gdf_topbtm = gpd.GeoDataFrame(
        {0: [10.0] * 16, 1: [5.0] * 16, 2: [0.0] * 16}, geometry=grid.gdf_vorPolys.geometry,
        crs=grid.crs,
    )
    context = ModelContext(grid=grid, domain=np.ones((2, 16), dtype=int))
    lakes = gpd.GeoDataFrame({"name": ["pond"]}, geometry=[box(25, 25, 75, 75)], crs=grid.crs)
    disv = mf.disv(nlay=2, top=10.0, botm=[5.0, 0.0], **grid.get_disv_gridprops())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lak = mf.lak(
            context=context, nper=3, lakes=lakes, lake_id_field="name", starting_stage=9.0,
            lake_bottom=6.0, connection_modes="rectangular", lake_top=10.0, bed_leakance=0.1,
            rainfall={1: 0.01}, budget_filerecord="m.lak.cbc", stage_filerecord="m.lak.stg",
        )
    if empty:
        lak = _with_empty_blocks(lak, 3)
    return [disv, mf.npf(k=1.0, icelltype=1), mf.ic(strt=9.0),
            mf.sto(transient={0: True}, sy=0.2, ss=1e-5, iconvert=1),
            mf.chd(stress_period_data={0: [[(1, 0), 9.0], [(1, 15), 8.0]]}), lak], "lak"


#: One numbered workspace per run, all under the one scratch root `main` removes.
_RUNS = itertools.count()


def run(factory, empty, inflow, root):
    packages, kind = factory(empty, inflow)
    model = mf.gwf("m", newtonoptions="NEWTON", packages=[
        *packages,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    sim = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=3, perioddata=[(1.0, 2, 1.0)] * 3),
        mf.ims(complexity="COMPLEX", linear_acceleration="BICGSTAB"),
    ])
    workspace = root / f"run_{next(_RUNS)}"
    workspace.mkdir()
    simulation = sim.build(str(workspace)).simulation
    simulation.write_simulation(silent=True)
    ok, report = simulation.run_simulation(silent=True)
    assert ok, report[-10:]
    gwf = simulation.get_model("m")
    plots = {}
    for what, call in (("gwf.plot()", gwf.plot), (f"{kind}.plot()", gwf.get_package(kind).plot)):
        try:
            call()
            plots[what] = "ok"
        except AttributeError as error:
            plots[what] = f"AttributeError: {error}"
    blocks = (workspace / f"m.{kind}").read_text().count("BEGIN period")
    return workspace, plots, blocks, kind


def main() -> None:
    root = pathlib.Path(tempfile.mkdtemp())
    try:
        cases = [
            ("SFR, inflow in period 1 only", _sfr_model, {1: {"s0": 50.0}}),
            ("SFR, no setting in any period", _sfr_model, None),
            ("LAK, rainfall in period 1 only", _lak_model, None),
        ]
        for label, factory, inflow in cases:
            print(f"== {label}")
            out = {}
            for empty in (True, False):
                workspace, plots, blocks, kind = run(factory, empty, inflow, root)
                out[empty] = workspace
                print(f"   {'with' if empty else 'without'} empty blocks: {blocks} period "
                      f"block(s) written; {plots}")
            for name in ("m.cbc", f"m.{kind}.cbc", "m.hds", f"m.{kind}.stg"):
                same = filecmp.cmp(out[True] / name, out[False] / name, shallow=False)
                print(f"   {name}: byte-identical = {same}")
    finally:
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    main()
