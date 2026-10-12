"""What the single-model read views' docstrings say they do, measured.

Ledger 441. The parameter entries written for ``modflow/mf6/{package_inputs,
headsplus,budget}.py``, ``modflow/utils/inputs.py``, ``project/model_config.py``,
``library.py`` and ``modflow/mf6/observations/registry.py`` (ledger 338's guard,
extended to the read grammar by ledger 436) state what each parameter ACCEPTS.
Several of those statements describe behaviour a reader would not guess -- some of it
defective, and documented as it is rather than fixed in the same pass. This measures
each of them. Run it::

    python docs/reproducers/view_docstring_claims.py [workspace]

It needs ``mf6`` on PATH and writes nothing outside the workspace it is given (a
temporary directory, removed at the end, by default). Under a minute.

**The bench.** ``CanonicalModelConfig.testing()`` (21 x 21, four layers, six periods,
every package), built and run once; a second, unrun copy is edited in memory where a
measurement needs a package written differently.

WHAT IT MEASURES
  1. With GHB's periods 1-5 removed (MF6 runs period 0's records in all six),
     ``ghb.inputs.get(per=3)`` had no rows, ``ghb.inputs.cond.map(per=3)`` coloured no
     cell, and the read ADDED an empty block 3, so the model written afterwards had no
     GHB records from period 3 on against 84. Fixed by ledger 449 and printed here
     as now: 84 rows at ``per=3``, no block added, 84 records in every period.
  2. A record noun's ``**trace_kwargs`` tail carried the per-layer names the noun does
     not name to ``model.plot.map``, unrefused: ``show_mounding=True`` kept every colour
     and changed every tooltip, on five nouns. Refused since ledger 465, and printed
     as now.
  3. ``npf.transmissivity_max.get(per=...)`` ignores ``per``.
  4. ``model.hds.long()``/``wide()`` accept only ``values="elev"`` (``"head"`` raises)
     and kept MF6's no-value sentinels where ``get()`` masks them -- they mask them
     too since ledger 484, printed as now (``NaN`` beside ``get()``'s).
  5. ``model.hds.get_obs_cells(crs=...)`` returned the same cells whatever ``crs`` named,
     a points file with no CRS raised, and a polygon over several cells got the first
     of them; ``get_obs_heads`` read the constructor's ``obs_path`` when given a ``str``
     path, and refused a second ``Path`` on one reader. Fixed by ledger 476 and
     printed as now: ``crs`` labels a file carrying none (naming another raises), the
     polygon raises, and each path given is read.
  6. ``model.bud("drn").plot_budget_obs``: a polygon file in another CRS read 0.0 for
     every area; a multi-part polygon raised; a ``q`` the table lacks returned an empty
     table; a second record on one cell was dropped; the default y range was ``[0, 4]``;
     a file with no CRS raised. Fixed by ledger 476 and printed as now: the same
     sums in either CRS, records added, a missing ``q`` refused, the y range fitted.
  7. ``model.outputs.lak.bud.get(..., return_df=0)`` returns a DataFrame (only ``False``
     returns the raw list).
  8. ``model.outputs.sfr.bud.plot_flows()`` raised ``KeyError: 'rno'`` (FloPy 3.10 names
     the column ``ifno``); it draws since ledger 476.
  9. ``model.inputs.uzf.finf``: a NumPy integer ``per`` raised ``UnboundLocalError``
     (reads since ledger 449); two land-surface UZF objects under one cell raise
     ``ValueError`` naming the cells (ledger 476 draws a stacked column's
     land-surface object). ``model.inputs.rch.map()`` raised ``KeyError: None`` (ledger
     449), and on an RCH carrying boundnames ``multiplier=2.5`` raised ``TypeError`` and
     ``2`` doubled each name; the multiplier scales only the numbers since ledger
     476.
 10. ``AoiLibrary.setdefault``: a non-string key was stored and then raised
     ``TypeError``; it is refused before anything is stored since ledger 476.
     ``setdefault("new")`` raises ``TypeError``.
 11. ``TargetRegistry``: a set registered as ``registry[1]`` was ``in`` the registry and
     ``registry[1]`` raised ``KeyError``; a non-string name raises ``TypeError`` since
     ledger 476.
 12. ``ModelConfig.from_model(model.gwf)`` was empty, and ``config.section("zzz")`` was
     ``{}``, with no error; both raise since ledger 476. ``ModelConfig.model_options`` reported ``"NEWTON"`` for a
     model whose ``newtonoptions`` is ``""`` or ``[]``, which FloPy writes as no option;
     it leaves the option out since ledger 450.
"""

from __future__ import annotations

import contextlib
import io
import shutil
import sys
import tempfile
import warnings
from pathlib import Path

import flopy
import geopandas as gpd
import numpy as np
import pandas as pd
import plotly.io as pio

import myflopy as mf
from myflopy.library import AoiLibrary
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, build_canonical_model
from myflopy.modflow.mf6.observations.registry import TargetRegistry
from myflopy.project.model_config import ModelConfig

pio.renderers.default = "json"  # `fig.show()` prints JSON (swallowed) instead of a browser


def _attempt(label, call):
    """Print what ``call`` returned or the exception it raised, and return the result."""

    try:
        with contextlib.redirect_stdout(io.StringIO()):
            result = call()
    except Exception as error:  # noqa: BLE001 - a reproducer reports every outcome
        print(f"  {label}: {type(error).__name__}: {str(error)[:110]}")
        return None
    shown = (
        f"{type(result).__name__} {getattr(result, 'shape', '')}"
        if hasattr(result, "shape") else repr(result)[:110]
    )
    print(f"  {label}: {shown}")
    return result


def _build(root: Path, name: str, *, run: bool):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = build_canonical_model(root / name, config=CanonicalModelConfig.testing())
    if run:
        with contextlib.redirect_stdout(io.StringIO()):
            ok, _ = model.run_simulation()
        assert ok, "the canonical model did not run"
    return model


def _ghb_records_per_period(model) -> list[int]:
    """Write and run ``model``, then count the GHB records MF6 applied in each period."""

    with contextlib.redirect_stdout(io.StringIO()):
        model.sim.write_simulation()
        ok, _ = model.sim.run_simulation(silent=True)
    assert ok, f"{model.name} did not run"
    budget = flopy.utils.CellBudgetFile(str(Path(model.sim.sim_path) / f"{model.name}.cbc"))
    return [
        sum(len(record) for record in budget.get_data(text="GHB", kstpkper=key))
        for key in budget.get_kstpkper()
    ]


def _without_ghb_blocks_after_0(model):
    """Remove GHB's period blocks 1..nper-1, so MF6 runs period 0's records in all."""

    package = model.package("ghb")
    for period in range(1, model.nper):
        package.stress_period_data.remove_transient_key(period)
    return package


def _trace(picture):
    trace = picture.fig.data[0]
    return np.asarray(trace.z, float), (trace.hovertemplate, str(trace.customdata))


def main(root: Path) -> None:
    model = _build(root, "run", run=True)
    edited = _build(root, "edited", run=False)

    print("1. input tables hold the period BLOCKS written")
    layer = int(edited.ghb.inputs.get(per=0).layer.iloc[0])
    ghb = _without_ghb_blocks_after_0(edited)
    print(f"  GHB blocks: {sorted(ghb.stress_period_data.get_data())}")
    _attempt("ghb.inputs.get(per=0) rows", lambda: len(edited.ghb.inputs.get(per=0)))
    _attempt("ghb.inputs.get(per=3) rows", lambda: len(edited.ghb.inputs.get(per=3)))
    for per in (0, 3):
        z, _ = _trace(edited.ghb.inputs.cond.map(per=per, layer=layer))
        print(f"  ghb.inputs.cond.map(per={per}): {int((z != 0).sum())} cells coloured")
    print(f"  GHB blocks after those reads: {sorted(ghb.stress_period_data.get_data())}")
    unread = _build(root, "unread", run=False)
    _without_ghb_blocks_after_0(unread)
    print(f"  written without the reads, GHB records per period: "
          f"{_ghb_records_per_period(unread)}")
    print(f"  written after the reads, GHB records per period: "
          f"{_ghb_records_per_period(edited)}")

    print("2. a record noun's tail forwards show_mounding")
    nouns = {
        "ghb.inputs.cond": lambda **k: model.ghb.inputs.cond.map(**k),
        "npf.k": lambda **k: model.npf.k.map(**k),
        "uzf.inputs.finf": lambda **k: model.uzf.inputs.finf.map(**k),
        "inputs.uzf": lambda **k: model.inputs.uzf.map(**k),
        "npf.transmissivity": lambda **k: model.npf.transmissivity.map(**k),
    }
    for label, draw in nouns.items():
        _attempt(f"{label}.map(show_mounding=True)", lambda d=draw: d(show_mounding=True))

    print("3. transmissivity_max ignores per")
    tmax = model.npf.transmissivity_max
    print(f"  get(per=3).equals(get()): {tmax.get(per=3).equals(tmax.get())}")

    print("4. the heads value table")
    hds = model.hds
    _attempt("hds.long(values='head')", lambda: hds.long(values="head"))
    table = hds.all_values
    (kstp, kper), lay, cell = table.index[0]
    table.iloc[0, 0] = -1e30  # what MF6 writes for a DRY cell
    print(f"  after marking (({kstp}, {kper}), {lay}, {cell}) dry: long() {hds.long().iloc[0]}, "
          f"get() {hds.get(per=kper, layer=lay, cells=cell)['head'].iloc[0]}")

    print("5. head observation points")
    centroids = model.vor.gdf_vorPolys.geometry.centroid
    points = root / "points.gpkg"
    gpd.GeoDataFrame({"ExploName": ["b", "a"]}, geometry=[centroids.iloc[10], centroids.iloc[200]],
                     crs=hds.crs).to_file(points)
    other = root / "other.gpkg"
    gpd.GeoDataFrame({"ExploName": ["z"]}, geometry=[centroids.iloc[300]], crs=hds.crs).to_file(other)
    naive = root / "naive.gpkg"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gpd.GeoDataFrame({"ExploName": ["n"]}, geometry=[centroids.iloc[5]]).to_file(naive)
    _attempt("get_obs_cells(points)", lambda: model.hds.get_obs_cells(points))
    _attempt("get_obs_cells(points, crs='EPSG:4326')",
             lambda: model.hds.get_obs_cells(points, crs="EPSG:4326"))
    _attempt("get_obs_cells(no-CRS file)", lambda: model.hds.get_obs_cells(naive))
    spanning = root / "spanning.gpkg"
    over = model.vor.gdf_vorPolys.loc[[0, 1, 22]].union_all().convex_hull.buffer(-1)
    gpd.GeoDataFrame({"ExploName": ["poly"]}, geometry=[over], crs=hds.crs).to_file(spanning)
    _attempt("get_obs_cells(a polygon over cells 0, 1 and 22)",
             lambda: model.hds.get_obs_cells(spanning))
    reader = type(hds)(model=model, obs_path=points)
    heads = _attempt("get_obs_heads(str(other)) on a reader opened with obs_path=points",
                     lambda: reader.get_obs_heads(str(other)))
    if heads is not None:
        print(f"    columns read: {list(heads.columns)} (other.gpkg holds 'z')")
    reader = model.hds
    _attempt("get_obs_heads(points)", lambda: reader.get_obs_heads(points))
    _attempt("get_obs_heads(other) on the same reader", lambda: reader.get_obs_heads(other))

    print("6. Budget.plot_budget_obs")
    bud = model.bud("drn")
    cells = sorted(set(bud.df.index.get_level_values(0)))
    polys = model.vor.gdf_vorPolys
    hull = polys.loc[cells[:3]].union_all().convex_hull.buffer(-1)
    areas = gpd.GeoDataFrame({"name": ["A"]}, geometry=[hull], crs=model.vor.crs)
    same = root / "areas.gpkg"
    areas.to_file(same)
    wgs = root / "areas_4326.gpkg"
    areas.to_crs(4326).to_file(wgs)
    split = root / "areas_multi.gpkg"
    gpd.GeoDataFrame({"name": ["M"]}, geometry=[polys.loc[cells[:3]].union_all().buffer(-1)],
                     crs=model.vor.crs).to_file(split)
    _attempt("areas in the grid's CRS, row 0",
             lambda: bud.plot_budget_obs(shp_gpkg=same, plot_fig=False).iloc[0].tolist())
    _attempt("areas in EPSG:4326, row 0",
             lambda: bud.plot_budget_obs(shp_gpkg=wgs, plot_fig=False).iloc[0].tolist())
    _attempt("a multi-part area", lambda: bud.plot_budget_obs(shp_gpkg=split, plot_fig=False))
    split_shp = root / "areas_multi.shp"
    gpd.read_file(split).to_file(split_shp)
    _attempt("the same multi-part area in a .shp, row 0",
             lambda: bud.plot_budget_obs(shp_gpkg=split_shp, plot_fig=False).iloc[0].tolist())
    naive_areas = root / "areas_no_crs.gpkg"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gpd.GeoDataFrame({"name": ["A"]}, geometry=[hull]).to_file(naive_areas)
    _attempt("a file with no CRS",
             lambda: bud.plot_budget_obs(shp_gpkg=naive_areas, plot_fig=False))
    _attempt("q='nope'", lambda: bud.plot_budget_obs(shp_gpkg=same, plot_fig=False, q="nope"))
    _attempt("a str path", lambda: bud.plot_budget_obs(shp_gpkg=str(same), plot_fig=False))
    figure = bud.plot_budget_obs(shp_gpkg=same, plot_fig=False, return_fig=True)
    print(f"  default y range: {figure.layout.yaxis.range}")
    single = bud.plot_budget_obs(shp_gpkg=same, plot_fig=False).iloc[0].tolist()
    bud._df = pd.concat([bud.df, bud.df])  # every cell now carries a second, equal record
    doubled = bud.plot_budget_obs(shp_gpkg=same, plot_fig=False).iloc[0].tolist()
    print(f"  one record per cell {single} / two records per cell {doubled}")

    print("7. LakBudget.get(return_df=...)")
    _attempt("return_df=0", lambda: model.outputs.lak.bud.get("GWF", return_df=0))
    _attempt("return_df=False", lambda: type(model.outputs.lak.bud.get("GWF", return_df=False)))

    print("8. SFRBudget.plot_flows")
    _attempt("plot_flows()", lambda: model.outputs.sfr.bud.plot_flows())

    print("9. the legacy inputs helpers")
    _attempt("inputs.uzf.finf(np.int64(0))", lambda: model.inputs.uzf.finf(np.int64(0)))
    _attempt("inputs.rch.map()", lambda: model.inputs.rch.map())
    uzf = edited.package("uzf")
    packagedata = uzf.packagedata.get_data().copy()
    packagedata["cellid"][1] = (1, packagedata["cellid"][0][1])  # ifno 1 under ifno 0's cell
    uzf.packagedata.set_data(packagedata)
    _attempt("inputs.uzf.finf(0), two objects under one cell", lambda: edited.inputs.uzf.finf(0))
    rch = edited.package("rch")
    rch.boundnames.set_data(True)
    rch.stress_period_data.set_data({
        period: [
            (cell, rate, "zone")
            for cell, rate in zip(records["cellid"], records["recharge"], strict=True)
        ]
        for period, records in rch.stress_period_data.get_data().items()
    })
    for multiplier in (2, 2.5):
        _attempt(f"inputs.rch.map(0, multiplier={multiplier}) with boundnames, first name",
                 lambda m=multiplier: edited.inputs.rch.map(0, multiplier=m)._custom_hover["boundname"][0])

    print("10. AoiLibrary.setdefault")
    library = AoiLibrary()
    _attempt("setdefault(1, box)", lambda: library.setdefault(1, (0, 0, 1, 1)))
    print(f"  stored keys: {list(library)}")
    _attempt("setdefault('new')", lambda: library.setdefault("new"))

    print("11. TargetRegistry keys")
    registry = TargetRegistry(None)
    conc = mf.ConcTargets(
        locations=pd.DataFrame({"name": ["w"], "cell": [3], "layer": [0]}),
        values=pd.DataFrame({"name": ["w"], "time": [0], "conc": [1.0]}),
    )
    _attempt("registry[1] = <ConcTargets>", lambda: registry.__setitem__(1, conc))
    print(f"  1 in registry: {1 in registry}")
    _attempt("registry[1]", lambda: registry[1])

    print("12. ModelConfig")
    _attempt("from_model(model.gwf).sections", lambda: ModelConfig.from_model(model.gwf).sections)
    _attempt("model.config.section('zzz')", lambda: model.config.section("zzz"))
    for index, newton in enumerate(("", [], "NEWTON")):
        workspace = root / f"newton_{index}"
        workspace.mkdir()
        sim = flopy.mf6.MFSimulation(sim_ws=str(workspace), verbosity_level=0)
        flopy.mf6.ModflowTdis(sim)
        gwf = flopy.mf6.ModflowGwf(sim, modelname="m", newtonoptions=newton)
        sim.write_simulation(silent=True)
        written = [line.strip() for line in (workspace / "m.nam").read_text().splitlines()
                   if "NEWTON" in line.upper()]
        print(f"  newtonoptions={newton!r}: name file writes {written}, "
              f"model_options reads {ModelConfig.model_options(gwf)}")


if __name__ == "__main__":
    given = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    workspace = given or Path(tempfile.mkdtemp(prefix="view_docstring_claims_"))
    try:
        main(workspace)
    finally:
        if given is None:
            shutil.rmtree(workspace, ignore_errors=True)
