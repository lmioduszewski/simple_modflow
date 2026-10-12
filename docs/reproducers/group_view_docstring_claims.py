"""What the GROUP read views' docstrings say they do, measured.

Ledger 439. The parameter entries written for ``project/group/{budget,inputs,results,sfr,uzf,lak,
spatial,surface_water}.py`` (ledger 338's guard, extended to the read grammar) state
what each selector ACCEPTS and how each table aligns two members. Several of those
statements describe behaviour a reader would not guess -- some of it defective, and
documented as it is rather than fixed in the same pass. This measures each of them.
Run it::

    python docs/reproducers/group_view_docstring_claims.py [workspace]

It needs ``mf6`` on PATH and writes nothing outside the workspace it is given (a
temporary directory, removed at the end, by default). About 40 s.

**The benches.**

* ``canon``: ``CanonicalModelConfig.testing()`` (21 x 21, every package), and
  ``canon+1``, the same with one more time step per period -- two members of
  DIFFERENT ``nstp``, the usual reason a group table has to align by period.
* ``ghb2``: a 600 ft square on a ``GridSpec.voronoi`` grid (one layer, top 20, bottom 0,
  K 10, CHD 10.0 on the west edge) with TWO GHB records on one east-edge cell, ``bhead``
  5.0 and 7.0 (member ``a``) and 6.0 and 8.0 (member ``b``), conductance 100. Each record
  changes by exactly 1.0.
* ``dry``: the sentinel bench of ``tests/conftest.py::sentinel_run`` (600 ft square,
  ~92 cells, cell bottoms 10,050 ft in the north-east corner) run twice on ONE
  discretization: CHD at 10,080/10,090 ft (``wet``: every cell wet) and at
  9,990/10,010 ft (``dry``: the north-east corner DRY, -1e30).
* ``inherit``: three 1-day periods, CHD declared in period 0 ONLY (MF6 runs periods 1-2
  on that block), no LAK, SFR or UZF.

WHAT IT MEASURES
  1. ``lak.results.q.compare()`` and ``group.bud(...).compare()`` pair members on the
     EXACT ``kstpkper``: at nstp 1 against 2 they diff the reference's period end
     against the other member's FIRST step, and LAK's also keys on the wetted
     ``FLOW-AREA`` and drops rows where it differs. ``sfr.results.q.compare()`` and
     ``lak.results.q.compare_map()`` reduce each member to its own period end first.
  2. On a cell holding two records, the group inputs ``map()`` ADDED ``bhead`` (12.0
     where the single-model map draws 6.0), ``compare()`` cross-paired the records
     (four rows) and ``compare_map()`` drew 4.0 for a change of 1.0 -- through the
     public ``group.diff()`` route too. Fixed by ledger 472, and printed as now:
     6.0, two rows, 1.0.
  3. ``group.hds.compare()`` differences MF6's -1e30 at a cell dry in one member only,
     and ``compare_map()``'s default colour range becomes plus/minus 1e30.
  4. A member with no LAK or UZF package ends the group LAK and UZF views with
     ``AttributeError``; the inputs, cell-result, SFR and surface-water views leave it
     out.
  5. CHD declared in period 0 only: the budget holds the same cells in every period.
     The input table had no rows at ``per=1`` (single-model and group alike) until
     ledger 449; it now reads the block MF6 has in force, so the counts agree.
  6. ``values=`` in a ``compare_map`` tail REPLACES the drawn differences.
  7. ``lak.connections.map``'s ``multiplier`` acts AFTER ``agg`` and scales the filled
     cells too.
  8. What the selectors accept: a tuple ``layer`` and a bare-int ``cells`` on
     ``group.hds.get``, a NumPy ``layer`` on the surface-water table, an unknown
     connection type, an unknown ``include``, an unknown ``model_name``.
  9. ``compare_map`` refuses a member whose discretization file differs from the
     reference's in ``botm`` alone.
 10. The UZF ``sat`` accessor's table holds TWO ``sat`` columns (NPF's ``DATA-SAT``
     record's 0.0 flow renamed beside the saturation), so its ``compare``,
     ``compare_map`` and ``map`` raise.
"""

from __future__ import annotations

import dataclasses
import shutil
import sys
import tempfile
import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
from shapely.geometry import box

import myflopy as mf
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, build_canonical_model
from myflopy.project.group._shared import _reduce_to_period_end, _stable_compare_keys


def _attempt(label, call):
    """Print what ``call`` returned (a frame's row count) or the exception it raised."""

    try:
        result = call()
    except Exception as exc:  # noqa: BLE001 - a reproducer reports whatever is raised
        print(f"    {label}: {type(exc).__name__}: {str(exc)[:110]}")
        return None
    rows = len(result) if hasattr(result, "__len__") else result
    print(f"    {label}: {rows} rows" if hasattr(result, "columns") else f"    {label}: {rows}")
    return result


def _z(picture, ncpl):
    """The per-cell values a map picture draws."""

    for trace in picture.fig.data:
        values = getattr(trace, "z", None)
        if values is not None and len(values) == ncpl:
            return np.asarray(values, dtype=float)
    raise AssertionError("no per-cell trace")


def _square(root: Path, max_area: float):
    root.mkdir(parents=True, exist_ok=True)
    boundary = root / "boundary.gpkg"
    gpd.GeoDataFrame({"n": ["d"]}, geometry=[box(0, 0, 600, 600)], crs=2927).to_file(boundary)
    return mf.GridSpec.voronoi(
        boundary=mf.GeoPackageSourceSpec(path=boundary), max_area=max_area
    ).resolve(workspace=root / "grid")


def _run(root: Path, name: str, grid, packages, *, nper: int = 1):
    props = grid.get_disv_gridprops()
    disv, *rest = packages
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=props["ncpl"], nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], **disv),
        *rest,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=nper, perioddata=[(1.0, 1, 1.0)] * nper), mf.ims(models=["m"])])
    run = mf.Project(root / name, name=name).run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    return mf.load_mf6_run(run.model().workspace)


def section_1(canon, canon_plus):
    print("\n1. members of different nstp: which steps a compare pairs")
    group = mf.ModelGroup({"a": canon, "b": canon_plus}, reference="a")
    print(f"    a saves {canon.kstpkper[:2]}...; b saves {canon_plus.kstpkper[:4]}...")
    lak = group.packages.lak.results.q
    compared = lak.compare()
    print(f"    lak compare(): {len(compared)} rows of a's {len(lak.get(model_name='a'))};"
          f" b's steps used: {sorted(set(compared.kstpkper))}")
    data = _reduce_to_period_end(lak.get())
    keys = _stable_compare_keys(data, ["q_lake", "q_per_area"])
    reference = data[data.model == "a"][keys + ["q_lake"]].rename(columns={"q_lake": "ref"})
    at_end = data[data.model == "b"].merge(reference, on=keys)
    paired = compared.merge(
        at_end[["per", "iconn", "lake", "q_lake"]].rename(columns={"q_lake": "b_end"}),
        on=["per", "iconn", "lake"],
    )
    print(f"    period-end alignment pairs {len(at_end)} rows; on the rows compare() kept,"
          f" b's first step differs from its period end by up to"
          f" {np.nanmax(np.abs(paired.q_lake - paired.b_end)):.2f} ft3/d")
    budget = group.bud("ghb").compare()
    print(f"    bud('ghb').compare(): b's steps used: {sorted(set(budget.kstpkper))[:3]}...")
    sfr = group.packages.sfr.results.q.compare()
    print(f"    sfr compare(): {len(sfr)} rows, aligned by per (no kstpkper column:"
          f" {'kstpkper' not in sfr.columns})")
    lak_map = lak.compare_map(model_name="b", per=5)
    print(f"    lak compare_map(per=5) draws {np.count_nonzero(_z(lak_map, canon.vor.ncpl))}"
          " non-zero cells, each member at its own period end")


def section_2(root: Path):
    print("\n2. two GHB records on one cell: bhead 5, 7 (a) and 6, 8 (b)")
    grid = _square(root / "ghb2", 40_000.0)
    edges = grid.gdf_vorPolys.geometry.bounds
    west = [[(0, int(c)), 10.0] for c in np.flatnonzero(edges.minx <= 1e-6)]
    cell = int(np.flatnonzero(edges.maxx >= 600 - 1e-6)[0])

    def member(name, shift):
        ghb = [[(0, cell), 5.0 + shift, 100.0], [(0, cell), 7.0 + shift, 100.0]]
        return _run(root / "ghb2", name, grid, [
            {"top": 20.0, "botm": [0.0]}, mf.npf(k=10.0), mf.ic(strt=8.0),
            mf.chd(stress_period_data={0: west}), mf.ghb(stress_period_data={0: ghb}),
        ])

    a, b = member("pa", 0.0), member("pb", 1.0)
    group = mf.ModelGroup({"a": a, "b": b}, reference="a")
    ncpl = a.vor.ncpl
    print(f"    single-model a.ghb.inputs.bhead.map(): {_z(a.ghb.inputs.bhead.map(), ncpl)[cell]}")
    print(f"    group.packages.ghb.inputs.bhead.map(): "
          f"{_z(group.packages.ghb.inputs.bhead.map(), ncpl)[cell]}")
    rows = group.packages.ghb.inputs.compare(cells=[cell])
    print(f"    group compare(): {len(rows)} rows, bhead_diff {rows.bhead_diff.tolist()}")
    drawn = group.packages.ghb.inputs.compare_map(model_name="b", value_column="bhead")
    print(f"    compare_map(value_column='bhead'): {_z(drawn, ncpl)[cell]} (each record moved 1.0)")
    public = group.diff().packages.ghb.inputs.bhead.map("b")
    print(f"    group.diff().packages.ghb.inputs.bhead.map('b'): {_z(public, ncpl)[cell]}")
    return group


def section_3(root: Path):
    print("\n3. a cell dry in one member only")
    grid = _square(root / "dry", 12_000.0)
    x, y = np.asarray(grid.centroids_x), np.asarray(grid.centroids_y)
    raised = (x > 400) & (y > 400)
    edges = grid.gdf_vorPolys.geometry.bounds

    def member(name, dry):
        west, east = (9_990.0, 10_010.0) if dry else (10_080.0, 10_090.0)
        chd = [[(0, int(c)), west] for c in np.flatnonzero(edges.minx <= 1e-6)]
        chd += [[(0, int(c)), east] for c in np.flatnonzero(edges.maxx >= 600 - 1e-6)
                if not (dry and raised[c])]
        return _run(root / "dry", name, grid, [
            {"top": 10_100.0, "botm": [np.where(raised, 10_050.0, 9_900.0)]},
            mf.npf(k=10.0, icelltype=1), mf.ic(strt=10_000.0), mf.chd(stress_period_data={0: chd}),
        ])

    group = mf.ModelGroup({"wet": member("wet", False), "dry": member("dry", True)},
                          reference="wet")
    print(f"    max |diff| in group.hds.compare(): {group.hds.compare()['diff'].abs().max():.3g}")
    trace = group.hds.compare_map(model_name="dry").fig.data[0]
    print(f"    compare_map colour range: {trace.zmin:.3g} .. {trace.zmax:.3g}")


def section_4(canon, inherit):
    print("\n4. one member without LAK, SFR or UZF")
    group = mf.ModelGroup({"a": canon, "o": inherit}, reference="a")
    for label, call in [
        ("lak.results.q.get()", lambda: group.packages.lak.results.q.get()),
        ("lak.connections.get()", lambda: group.packages.lak.connections.get()),
        ("uzf.inputs.finf.get()", lambda: group.packages.uzf.inputs.finf.get()),
        ("sfr.results.q.get()", lambda: group.packages.sfr.results.q.get()),
        ("ghb.inputs.get()", lambda: group.packages.ghb.inputs.get()),
        ("surface_water.results.q.get()", lambda: group.packages.surface_water.results.q.get()),
    ]:
        _attempt(label, call)


def section_5(inherit):
    print("\n5. CHD declared in period 0 only (MF6 runs periods 1-2 on that block)")
    group = mf.ModelGroup({"a": inherit, "b": inherit}, reference="a")
    print(f"    model.chd.inputs.get(per=k): {[len(inherit.chd.inputs.get(per=k)) for k in range(3)]}")
    inputs = group.packages.chd.inputs
    print(f"    group inputs.get(per=k): {[len(inputs.get(model_name='a', per=k)) for k in range(3)]}")
    results = group.packages.chd.results.q
    print(f"    group results.q.get(per=k): {[len(results.get(model_name='a', per=k)) for k in range(3)]}")


def section_6(ghb_group):
    print("\n6. values= in a compare_map tail")
    ncpl = ghb_group.models["a"].vor.ncpl
    drawn = ghb_group.packages.ghb.inputs.compare_map(
        model_name="b", values=np.arange(ncpl, dtype=float))
    print(f"    draws the array passed: {np.allclose(_z(drawn, ncpl), np.arange(ncpl))}")


def section_7(canon):
    print("\n7. lak.connections.map: multiplier against agg and fill_value")
    group = mf.ModelGroup({"a": canon, "b": canon}, reference="a")
    connections = group.packages.lak.connections
    ncpl = canon.vor.ncpl
    raw = connections.get(model_name="a")
    plain = _z(connections.map(value_column="connlen", agg="max"), ncpl)
    flipped = _z(connections.map(value_column="connlen", agg="max", multiplier=-1.0), ncpl)
    cell = int(np.nanargmax(plain))
    in_cell = raw.loc[raw.cell == cell, "connlen"]
    print(f"    cell {cell}: connlen {sorted(in_cell.round(2).tolist())}; agg='max' draws"
          f" {plain[cell]:.2f}, with multiplier=-1 {flipped[cell]:.2f}"
          f" (max of the negated values would be {(-in_cell).max():.2f})")
    filled = _z(connections.map(fill_value=-1.0, multiplier=2.0), ncpl)
    print(f"    fill_value=-1, multiplier=2: unconnected cells draw {sorted(set(filled[filled < 0]))}")


def section_8(canon):
    print("\n8. what the selectors accept")
    group = mf.ModelGroup({"a": canon, "b": canon}, reference="a")
    water = group.packages.surface_water.results.q
    for label, call in [
        ("hds.get(layer=(0, 1))", lambda: group.hds.get(layer=(0, 1))),
        ("hds.get(layer=[0, 1])", lambda: group.hds.get(layer=[0, 1])),
        ("hds.get(cells=5)", lambda: group.hds.get(cells=5)),
        ("hds.get(per=99)", lambda: group.hds.get(per=99)),
        ("surface_water get(layer=np.array([0, 1]))",
         lambda: water.get(layer=np.array([0, 1]))),
        ("surface_water get(include='lak', lak_connection_type='horiz')",
         lambda: water.get(include="lak", lak_connection_type="horiz")),
        ("surface_water get(include='uzf')", lambda: water.get(include="uzf")),
        ("lak.results.q.get(connection_type='horiz')",
         lambda: group.packages.lak.results.q.get(connection_type="horiz")),
        ("ghb.inputs.get(model_name='nope')", lambda: group.packages.ghb.inputs.get(model_name="nope")),
    ]:
        _attempt(label, call)


def section_9(root: Path):
    print("\n9. compare_map against a member whose botm alone differs")
    grid = _square(root / "botm", 40_000.0)
    edges = grid.gdf_vorPolys.geometry.bounds
    chd = [[(0, int(c)), 10.0] for c in np.flatnonzero(edges.minx <= 1e-6)]
    chd += [[(0, int(c)), 5.0] for c in np.flatnonzero(edges.maxx >= 600 - 1e-6)]
    members = {
        name: _run(root / "botm", name, grid, [
            {"top": 20.0, "botm": [bottom]}, mf.npf(k=10.0), mf.ic(strt=8.0),
            mf.chd(stress_period_data={0: chd}),
        ])
        for name, bottom in (("a", 0.0), ("b", -1.0))
    }
    group = mf.ModelGroup(members, reference="a")
    _attempt("group.hds.compare() rows", lambda: group.hds.compare())
    _attempt("group.hds.compare_map(model_name='b')", lambda: group.hds.compare_map(model_name="b"))


def section_10(canon):
    print("\n10. the UZF sat accessor")
    group = mf.ModelGroup({"a": canon, "b": canon}, reference="a")
    sat = group.packages.uzf.results.sat
    print(f"    get() columns: {list(sat.get().columns)}")
    for label, call in [
        ("compare()", lambda: sat.compare()),
        ("compare_map(model_name='b')", lambda: sat.compare_map(model_name="b")),
        ("map()", lambda: sat.map()),
    ]:
        _attempt(label, call)


def main() -> None:
    given = len(sys.argv) > 1
    root = Path(sys.argv[1]) if given else Path(tempfile.mkdtemp(prefix="group_claims_"))
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            config = CanonicalModelConfig.testing()
            built = {}
            for label, cfg in (("canon", config),
                               ("canon_plus", dataclasses.replace(
                                   config, steps_per_period=config.steps_per_period + 1))):
                model = build_canonical_model(root / label, config=cfg)
                assert model.run_simulation()[0]
                built[label] = mf.load_mf6_run(model.workspace)
            grid = _square(root / "inherit", 40_000.0)
            edges = grid.gdf_vorPolys.geometry.bounds
            chd = [[(0, int(c)), 10.0] for c in np.flatnonzero(edges.minx <= 1e-6)]
            chd += [[(0, int(c)), 5.0] for c in np.flatnonzero(edges.maxx >= 600 - 1e-6)]
            inherit = _run(root / "inherit", "p", grid, [
                {"top": 20.0, "botm": [0.0]}, mf.npf(k=10.0), mf.ic(strt=8.0),
                mf.sto(transient={0: True}, ss=1e-5, sy=0.1), mf.chd(stress_period_data={0: chd}),
            ], nper=3)

            section_1(built["canon"], built["canon_plus"])
            ghb_group = section_2(root)
            section_3(root)
            section_4(built["canon"], inherit)
            section_5(inherit)
            section_6(ghb_group)
            section_7(built["canon"])
            section_8(built["canon"])
            section_9(root)
            section_10(built["canon"])
    finally:
        if not given:
            shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    main()
