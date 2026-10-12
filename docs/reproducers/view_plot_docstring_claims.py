"""What the plot views' and composers' docstrings say they do, measured.

Ledger 438. The parameter entries written for ``modflow/mf6/package_plotting.py``,
``layers.py``, ``modflow/mf6/grid/plotting.py``, ``modflow/mf6/interactive_plotting.py``
and ``modflow/utils/datatypes/xsections.py`` (ledger 338's guard, extended to the read
grammar) state what each parameter does. Several statements describe behaviour a
reader would not guess -- defects documented as they are rather than fixed in the same
pass. This measures each of them. Run it::

    python docs/reproducers/view_plot_docstring_claims.py [workspace]

It needs ``mf6`` on PATH and writes nothing outside the workspace it is given (a
temporary directory, removed at the end, by default). About 30 s.

Items 10-15 were added by ledger 436, from the second (verification) pass.

Since ledgers 461, 462 and 463 (and a 2026-10-11 fix of item 14 and
15), items 1-3, 5, 7-9 and 12-15 print what the code does NOW -- each answered or
refused as its docstring says; the WHAT IT MEASURES list below is the defect as
first measured.

**The benches.**

* ``canon``: ``CanonicalModelConfig.testing()`` (21 x 21, every package), run once.
* ``stack``: a 400 x 300 ft ``TriangleGrid`` (angle 30) under a ``LayerStack`` of three
  flat 10 ft units, ``a``/``b``/``c``, top 50 ft.
* ``gwt``: a one-layer GWF ``flow`` + GWT ``trans`` pair on a finer copy of that grid
  (``max_area=1500``; the default leaves 4 cells, too few to section), CHD
  west (head 10, concentration 1) and east (head 8, concentration 0), one period.

WHAT IT MEASURES
  1. ``SpatialView.mosaic`` drops ``agg`` on ``kind="map"``: ``q.mosaic(agg="count")``
     draws the same cells as ``agg="sum"``, where ``q.map(agg="count")`` does not.
  2. ``mosaic(kind="section", cells=...)``: the named ``cells`` never reaches the
     section, which fails with ``object of type 'NoneType' has no len()``.
  3. ``mosaic(kind="plot", <anything>=...)`` accepts the keyword and draws.
  4. ``SpatialView.plot(layer=3)`` on a table with no ``layer`` column keeps every row.
  5. ``FieldMappable.section`` raises ``NotImplementedError`` on every field of every
     package namespace of the canonical model (50 of 50, ``surface_water`` included),
     and of a group of it (14 of 14).
  6. ``FieldMappable`` matches ``field=`` exactly (``"Stage"`` raises).
  7. ``LayerSection``: ``plot_mpl(color_by="bogus")`` draws the layer picture; under
     ``color_by="thickness"`` the fill is ``viridis`` whatever ``cmap`` says, with no
     legend; ``cmap="viridis"`` gives three layers legend colours within 0.02 of each
     other.
  8. ``stack.plot.map(1)`` raises ``ValueError: list.index(x): x not in list``; a name
     works.
  9. ``stack.plot.map(bounds_padding=0.1)`` without ``extent=`` raises ``TypeError``.
 10. ``model.drn.inputs.elev.plot()`` draws the SUM of the drains' elevations per
     period (``agg="sum"``), where the field's ``map`` takes its kind's rule (ledger 297).
 11. A diff's LAK ``results`` lists ``'stage'`` among its fields, and
     ``map(field='stage')``/``section(field='stage')`` raise ``AttributeError``.
 12. ``mosaic(kind="section", <unknown>=...)``'s ``TypeError`` lists ``cells`` as a
     parameter the section takes, which the composer's own ``cells`` captures.
 13. ``stack.plot.section(legend="bottom")`` puts the legend at the upper right.
 14. ``myflopy.plot.grid(vor, backend="mpl")`` draws onto pyplot's CURRENT figure.
 15. The interpolated profile of a GWT model's concentration reads ``<name>.hds`` and
     raises ``FileNotFoundError`` (the model writes ``trans.ucn``).
"""

from __future__ import annotations

import shutil
import sys
import tempfile
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from shapely.geometry import LineString  # noqa: E402

import myflopy as mf  # noqa: E402
from myflopy.layers import Flat, LayerStack  # noqa: E402
from myflopy.modflow.mf6.canonical_example import (  # noqa: E402
    CanonicalModelConfig,
    build_canonical_model,
)
from myflopy.modflow.mf6.package_plotting import FieldMappable, SpatialView  # noqa: E402


def _z(figure_or_choro):
    fig = getattr(figure_or_choro, "fig", figure_or_choro)
    return np.asarray(fig.data[0].z, dtype=float)


def canon(root: Path) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = build_canonical_model(root / "canonical", config=CanonicalModelConfig.testing())
    ok, _ = model.run_simulation()
    assert ok, "the canonical model did not run"
    q = model.drn.results.q

    same = np.allclose(_z(q.mosaic(agg="sum", layer=[0])), _z(q.mosaic(agg="count", layer=[0])),
                       equal_nan=True)
    differs = not np.allclose(_z(q.map(layer=0, agg="sum")), _z(q.map(layer=0, agg="count")),
                              equal_nan=True)
    print(f"1. mosaic agg='count' == agg='sum': {same}; map agg='count' != 'sum': {differs}")

    try:
        model.hds.mosaic(kind="section", cells=[0, 1, 2], by="period", per=[0, 1])
        print("2. mosaic(kind='section', cells=...): drew")
    except TypeError as error:
        print(f"2. mosaic(kind='section', cells=...): TypeError: {error}")

    try:
        figure = q.mosaic(kind="plot", by="layer", not_a_parameter=1)
        print(f"3. mosaic(kind='plot', not_a_parameter=1): {type(figure).__name__}")
    except TypeError as error:
        print(f"3. mosaic(kind='plot', not_a_parameter=1): TypeError: {error}")

    xmin, ymin, xmax, ymax = model.vor.gdf_vorPolys.total_bounds
    line = LineString([(xmin, (ymin + ymax) / 2), (xmax, (ymin + ymax) / 2)])
    raised, total = 0, 0
    for pkg in ("chd", "drn", "ghb", "riv", "wel", "rch", "evt", "lak", "sfr", "uzf",
                "surface_water"):
        for side in ("inputs", "results"):
            namespace = getattr(getattr(model, pkg), side, None)
            if not isinstance(namespace, FieldMappable):
                continue
            for field in namespace.field_names():
                total += 1
                node = getattr(namespace, field)
                if not hasattr(node, "section"):
                    continue
                try:
                    node.section(line=line)
                except NotImplementedError:
                    raised += 1
    print(f"5. FieldMappable.section raised NotImplementedError on {raised} of {total} fields "
          f"(offered by the namespace: {hasattr(FieldMappable, 'section')})")
    group = mf.ModelGroup({"a": model, "b": model}, reference="a")
    raised, total = 0, 0
    for pkg in dir(type(group.packages)):
        if pkg.startswith("_"):
            continue
        for side in ("inputs", "results"):
            namespace = getattr(getattr(group.packages, pkg, None), side, None)
            if not isinstance(namespace, FieldMappable):
                continue
            for field in namespace.field_names():
                total += 1
                node = getattr(namespace, field)
                if not hasattr(node, "section"):
                    continue
                try:
                    node.section(line=line)
                except NotImplementedError:
                    raised += 1
    print(f"   ... and on {raised} of {total} fields of a group of it")

    try:
        model.lak.results.map(field="Stage")
        print("6. field='Stage': drew")
    except ValueError as error:
        print(f"6. field='Stage': ValueError: {error}")

    elev = model.drn.inputs.elev
    table = elev.get()
    drawn = [float(y) for y in elev.plot().data[0].y]
    means = table.groupby("per")["elev"].mean().round(2).tolist()
    sums = table.groupby("per")["elev"].sum().round(2).tolist()
    print(f"10. drn.inputs.elev.plot() draws {[round(y, 2) for y in drawn][:3]}; "
          f"per-period sums {sums[:3]}, means {means[:3]}")

    diff_lak = group.diff().packages.lak.results
    listed = "stage" in diff_lak.field_names()
    for verb in ("map", "section"):
        try:
            getattr(diff_lak, verb)(field="stage")
            print(f"11. diff lak results {verb}(field='stage'): drew (listed: {listed})")
        except AttributeError as error:
            print(f"11. diff lak results {verb}(field='stage'): AttributeError: {error} "
                  f"(listed: {listed})")

    try:
        model.hds.mosaic(kind="section", by="period", per=[0, 1], not_a_parameter=1)
        print("12. mosaic(kind='section', not_a_parameter=1): drew")
    except TypeError as error:
        print(f"12. mosaic(kind='section', not_a_parameter=1): TypeError: {error}")


def series_stub() -> None:
    class Stub(SpatialView):
        value_name = "q"

        def get(self):
            return pd.DataFrame({"per": [0, 0, 1, 1], "cell": [1, 2, 1, 2],
                                 "lake": [1, 2, 1, 2], "q": [1.0, 2.0, 3.0, 4.0]})

    lines = [len(trace.x) for trace in Stub().plot(layer=3).data]
    print(f"4. plot(layer=3) on a table with no layer column: lines of {lines} points")


def stack(root: Path) -> None:
    tri = mf.TriangleGrid(model_ws=str(root / "tri"), angle=30)
    tri.set_domain_rectangle(x_dist=400, y_dist=300, origin=(0, 0))
    tri.build(verbose=False)
    result = (LayerStack(mf.VoronoiGridPlus(tri), top=Flat(50))
              .add("a", bottom=Flat(40)).add("b", bottom=Flat(30)).add("c", bottom=Flat(20))
              .build())

    try:
        bogus = result.plot.section(y=150).plot_mpl(color_by="bogus")
        print(f"7. plot_mpl(color_by='bogus') drew: {bogus.get_title()!r}")
    except ValueError as error:
        print(f"7. plot_mpl(color_by='bogus'): ValueError: {error}")
    cmaps = set()
    for cmap in ("Reds", "Blues"):
        axes = result.plot.section(y=150, color_by="thickness", cmap=cmap, legend=True).axes
        cmaps |= {c.get_cmap().name for c in axes.collections if c.get_array() is not None}
        legend = axes.get_legend() is not None or len(axes.figure.axes) > 1
    print(f"   thickness fill colormaps for cmap='Reds'/'Blues': {sorted(cmaps)}; "
          f"legend or colorbar: {legend}")
    handles = result.plot.section(y=150, cmap="viridis").axes.get_legend().legend_handles
    colours = np.array([h.get_facecolor()[:3] for h in handles])
    spread = float(np.abs(colours - colours[0]).max())
    print(f"   cmap='viridis' legend colours differ by at most {spread:.3f}")

    print(f"8. stack.plot.map('b'): {result.plot.map('b').title!r}")
    try:
        print(f"   stack.plot.map(1): {result.plot.map(1).title!r}")
    except ValueError as error:
        print(f"   stack.plot.map(1): ValueError: {error}")
    try:
        limits = result.plot.map(bounds_padding=0.1).axes.get_xlim()
        print(f"9. stack.plot.map(bounds_padding=0.1): drew, x limits {limits}")
    except TypeError as error:
        print(f"9. stack.plot.map(bounds_padding=0.1): TypeError: {error}")

    import matplotlib.pyplot as plt

    for where in ("bottom", "upper left"):
        legend = result.plot.section(y=150, legend=where).axes.get_legend()
        print(f"13. stack.plot.section(legend={where!r}): legend loc code {legend._loc} "
              "(1 is 'upper right')")
        plt.close("all")

    mine = plt.figure()
    mine.add_subplot(1, 1, 1, aspect="equal").plot([0, 1], [0, 1])
    drawn = mf.plot.grid(result.vor, backend="mpl")
    print(f"14. plot.grid(backend='mpl') returned the caller's current figure: {drawn is mine}")
    plt.close("all")


def gwt(root: Path) -> None:
    tri = mf.TriangleGrid(model_ws=str(root / "tri_gwt"), angle=30)
    tri.set_domain_rectangle(x_dist=400, y_dist=300, origin=(0, 0), max_area=1500)
    tri.build(verbose=False)
    vor = mf.VoronoiGridPlus(tri)
    gp = vor.get_disv_gridprops()
    ctx = mf.ModelContext(grid=vor, domain=np.ones((1, vor.ncpl), dtype=int))
    cx = np.asarray(vor.centroids_x, dtype=float)
    west = [int(c) for c in np.where(cx <= cx.min() + 60.0)[0]]
    east = [int(c) for c in np.where(cx >= cx.max() - 60.0)[0]]

    def disv():
        return mf.disv(nlay=1, ncpl=gp["ncpl"], nvert=len(gp["vertices"]),
                       vertices=gp["vertices"], cell2d=gp["cell2d"], top=10.0, botm=0.0)

    chd = [[(0, c), 10.0, 1.0] for c in west] + [[(0, c), 8.0, 0.0] for c in east]
    flow = mf.gwf("flow", context=ctx, packages=[
        disv(), mf.ic(strt=9.0), mf.npf(k=1.0, save_specific_discharge=True),
        mf.chd(stress_period_data={0: chd}, auxiliary="concentration"),
        mf.oc(head_filerecord="flow.hds", budget_filerecord="flow.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    trans = mf.gwt("trans", context=ctx, packages=[
        disv(), mf.ic(strt=0.0), mf.adv(scheme="UPSTREAM"), mf.mst(porosity=0.2),
        mf.ssm(sources=[["chd", "AUX", "concentration"]]),
        mf.oc(concentration_filerecord="trans.ucn", budget_filerecord="trans.cbc",
              saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")]),
    ])
    sim = mf.SimulationSpec("gwtcpl", models=[flow, trans], packages=[
        mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
        mf.ims(name="ims_flow", models=["flow"], complexity="SIMPLE"),
        mf.ims(name="ims_trans", models=["trans"], complexity="SIMPLE",
               linear_acceleration="BICGSTAB"),
    ], exchanges=[mf.exchange(flow, trans)])
    run = mf.Project(root / "gwtproj", name="p").prepare_run("r", sim, overwrite=True)
    ok, _ = run.execute()
    assert ok, "the GWT bench did not run"
    model = run.model("trans")
    written = sorted(p.name for p in Path(model.model_output_folder_path).iterdir()
                     if p.suffix in (".hds", ".ucn"))
    try:
        model.conc._sections(cells=[0, 5, 10], interpolate=True)[model.name].to_frame()
        series = model.conc._sections(cells=[0, 5, 10])[model.name].to_frame()["series"]
        print(f"15. interpolated GWT profile: drew (files {written}); stepped series "
              f"{sorted(set(series))}")
    except FileNotFoundError as error:
        print(f"15. interpolated GWT profile: FileNotFoundError: {error} (files {written})")


def main(argv: list[str]) -> None:
    owned = not argv
    root = Path(tempfile.mkdtemp(prefix="view_plot_claims_")) if owned else Path(argv[0])
    root.mkdir(parents=True, exist_ok=True)
    try:
        canon(root)
        series_stub()
        stack(root)
        gwt(root)
    finally:
        if owned:
            shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    main(sys.argv[1:])
