"""What the LAK / SFR / surface-water / budget / HFB read views' docstrings say, measured.

Ledger 440. The parameter entries written for
``modflow/mf6/{package_surface_water,package_results,package_model}.py`` (ledger 338's
guard, extended to the read grammar by ledger 436) state what each selector ACCEPTS
and what each picture draws. Several of those statements describe behaviour a reader
would not guess -- some of it defective, and documented as it is rather than fixed in
the same pass. This measures each of them. Run it::

    python docs/reproducers/surface_water_hfb_view_claims.py [workspace]

It needs ``mf6`` on PATH and writes nothing outside the workspace it is given (a
temporary directory, removed at the end, by default). About 30 s.

**The benches.** ``canon``: ``CanonicalModelConfig.testing()`` (21 x 21, every package,
six one-step periods; its SFR puts reaches 22, 37 and 38 in cell 195). ``hfb``: the same
model with a wall of barriers straight across the middle, declared in period 0
(``hydchr`` 1e-8), and a period-3 block naming only three of those faces (1e-6) --
MODFLOW 6 REPLACES a period's barriers with a later block's.

WHAT IT MEASURES
  1. ``model.sfr.results.q.get()`` joins each ``SFR`` record to the reach table on the
     CELL, not on the record's own reach (``node2``): a cell of three reaches carries
     nine rows, and ``model.sfr.results.profile`` gives each of the three the cell's
     total.
  2. ``model.sfr.budget.<term>.wide()``'s default ``index=("per", "lake")`` raises
     ``KeyError: 'lake'`` on every SFR term.
  3. ``model.<lak|sfr>.budget.get(term="ext_inflow")`` matches nothing, silently.
  4. Both SFR profile views draw by reach number for ANY ``x`` but ``"distance"``.
  5. ``model.lak.connections.map(per=1)`` raises a ``TypeError`` naming
     ``ModelPlots.map`` while its (spliced) docstring documents ``per``.
  6. ``multiplier`` on a record map is applied after ``agg`` and scales the
     ``fill_value`` cells -- but not when the selection holds no record at all.
  7. A map's unknown ``**trace_kwargs`` key raises at ``.fig`` on Plotly and is
     ignored under ``backend="mpl"``.
  8. ``model.budget[key]`` returns a namespace attribute for ``"model"``/``"types"``
     and raises ``AttributeError`` (not ``KeyError``) for an unknown term.
  9. ``exchange_intensity`` is length²/time on an SFR cell and length/time on a LAK
     cell.
 10. The HFB maps colour their backdrop by CELL NUMBER; the results map does not colour
     a barrier by its flow; ``inputs.map(per=None)`` draws every block; ``segments(per)``
     reads the block DECLARED in a period, not the one in force; ``results.q.get()``
     lists every block's barriers at every step; ``results.q.map`` with a ``per`` and a
     ``kstpkper`` that disagree raises ``KeyError: 'layer'``.
 11. ``model.lak.results.q``'s ``flow_area`` is MODFLOW 6's wetted ``FLOW-AREA`` (a
     vertical connection reads its cell's whole area or 0.0), not the connection
     table's ``connection_area``; a package budget's ``FLOW-JA-FACE`` is written from
     both ends, ``q`` positive INTO ``reach``; and ``model.surface_water.results.q``'s
     SFR ``q``/``rlen`` are tripled on cell 195 while ``exchange_intensity`` is right.
"""

from __future__ import annotations

import shutil
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np

from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, build_canonical_model


def _attempt(label, call):
    """Print what ``call`` returned, or the exception it raised."""

    try:
        result = call()
    except Exception as exc:  # noqa: BLE001 - a reproducer reports whatever is raised
        print(f"    {label}: {type(exc).__name__}: {str(exc)[:110]}")
        return None
    rows = f"{len(result)} rows" if hasattr(result, "columns") else repr(result)
    print(f"    {label}: {rows}")
    return result


def _z(picture):
    """The per-cell values a Plotly map draws (its choropleth trace)."""

    return np.asarray(picture.fig.data[0].z, dtype=float)


def _canonical(path: Path):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=r"mf\.lak: lake")
        model = build_canonical_model(path, config=CanonicalModelConfig.testing())
    success, report = model.run_simulation()
    assert success, "\n".join(report[-20:])
    return model


def _with_barriers(path: Path):
    import flopy
    import shapely

    from myflopy.modflow.mf6.grid.barriers import barrier_faces

    model = _canonical(path)
    xmin, ymin, xmax, ymax = model.vor.gdf_vorPolys.total_bounds
    middle = (ymin + ymax) / 2
    faces = barrier_faces(model.vor, shapely.LineString([(xmin, middle), (xmax, middle)]))
    flopy.mf6.ModflowGwfhfb(model.gwf, stress_period_data={
        0: [[[0, a], [0, b], 1e-8] for a, b in faces],
        3: [[[0, a], [0, b], 1e-6] for a, b in faces[:3]],
    })
    model.sim.write_simulation(silent=True)
    success, report = model.run_simulation()
    assert success, "\n".join(report[-20:])
    return model


def section_1(m):
    print("\n1. SFR records joined to reaches on the CELL")
    records = m._get_budget_reader().get_data(text="SFR", kstpkper=(0, 3))[0]
    at = records[records["node"] == 196]  # cell 195, one-based node
    print(f"    the budget file at cell 195: reaches {(at['node2'] - 1).tolist()},"
          f" q {np.round(at['q'], 3).tolist()}")
    rows = m.sfr.results.q.get(per=3, cells=[195])
    print(f"    sfr.results.q.get(per=3, cells=[195]): {len(rows)} rows for {len(at)} records")
    profile = m.sfr.results.profile.get(per=3)
    print("    sfr.results.profile, reaches 22/37/38 q_gwf:",
          np.round(profile.set_index("reach").loc[[22, 37, 38], "q_gwf"], 3).tolist())


def section_2_to_5(m):
    print("\n2. an SFR budget term's wide() default")
    _attempt("sfr.budget.gwf.wide()", lambda: m.sfr.budget.gwf.wide())
    _attempt("sfr.budget.gwf.wide(index=('per', 'reach'))",
             lambda: m.sfr.budget.gwf.wide(index=("per", "reach")))
    print("\n3. a term in its attribute spelling")
    _attempt("sfr.budget.get(term='ext_inflow')", lambda: m.sfr.budget.get(term="ext_inflow"))
    _attempt("sfr.budget.get(term='EXT-INFLOW')", lambda: m.sfr.budget.get(term="EXT-INFLOW"))
    print("\n4. the profile axis")
    for view in (m.sfr.results.profile, m.sfr.results.q.profile):
        title = view.plot(x="Reach-typo").layout.xaxis.title.text
        print(f"    {type(view).__name__}.plot(x='Reach-typo') x axis: {title!r}")
    print("\n5. lak.connections.map(per=...)")
    _attempt("lak.connections.map(per=1)", lambda: m.lak.connections.map(per=1))


def section_6_to_9(m):
    print("\n6. multiplier and fill_value")
    q = m.drn.results.q
    drawn = _z(q.map(multiplier=2.0, fill_value=1.0))
    holds = set(q.get(per=0, layer=0).cell)
    empty = next(c for c in range(len(drawn)) if c not in holds)
    print(f"    drn q.map(multiplier=2, fill_value=1): a cell with no record draws {drawn[empty]}")
    print(f"    the same in layer 2, which holds no drain: "
          f"{sorted(set(_z(q.map(layer=2, multiplier=2.0, fill_value=1.0)).tolist()))}")
    print("\n7. an unknown trace property")
    picture = q.map(bogus_property=1)
    _attempt("plotly: .fig", lambda: picture.fig)
    _attempt("backend='mpl'", lambda: type(q.map(backend="mpl", bogus_property=1)).__name__)
    print("\n8. model.budget[key]")
    _attempt("budget['model'] is the model", lambda: m.budget["model"] is m)
    _attempt("budget['types']", lambda: m.budget["types"][:3])
    _attempt("budget['nope']", lambda: m.budget["nope"])
    print("\n9. exchange_intensity units")
    table = m.surface_water.results.q.get(per=0)
    sfr = table[table.source == "sfr"].iloc[0]
    lak = table[table.source == "lak"].iloc[0]
    print(f"    SFR cell: -q/rlen = {-sfr.q / sfr.rlen:.4f} == {sfr.exchange_intensity:.4f}")
    print(f"    LAK cell: q/flow_area = {lak.q / lak.flow_area:.4f}"
          f" == {lak.exchange_intensity:.4f}")


def section_10(m):
    print("\n10. the HFB views")
    print(f"    barriers per declared block: {m.hfb.get().groupby('per').size().to_dict()}")
    backdrop = _z(m.hfb.map())
    print(f"    hfb.map() backdrop z[:5] = {backdrop[:5].tolist()} (cell numbers)")
    lines = [t for t in m.hfb.results.q.map().fig.data if t.type == "scattermap"]
    print(f"    hfb.results.q.map(): {len(lines)} barrier traces, colours set:"
          f" {sorted({str(t.line.color) for t in lines})}")
    overlay = [t for t in m.hfb.map().fig.data if t.type == "scattermap"][0]
    print(f"    hfb.map() draws {sum(v is None for v in overlay.lon)} segments;"
          f" segments(per=1): {len(m.hfb.segments(per=1))}")
    summary = m.hfb.results.q.summary()
    print("    results.q.summary() barriers / q_abs_total by period:",
          list(zip(summary.per, summary.barriers, np.round(summary.q_abs_total, 4), strict=True)))
    _attempt("hfb.results.q.map(per=2, kstpkper=(0, 3))",
             lambda: m.hfb.results.q.map(per=2, kstpkper=(0, 3)))


def section_11(m):
    print("\n11. LAK flow_area, FLOW-JA-FACE, surface-water q on a multi-reach cell")
    rows = m.lak.results.q.get().merge(
        m.lak.connections.get()[["lake", "iconn", "connection_area"]], on=["lake", "iconn"]
    )
    vertical = rows[rows.claktype == "VERTICAL"]
    print("    vertical flow_area / connection_area:",
          sorted(set((vertical.flow_area / vertical.connection_area).round(4))),
          "; zero per period:",
          vertical.assign(z=vertical.flow_area == 0).groupby("per").z.sum().to_dict(),
          f"of {int(vertical.groupby('per').size().iloc[0])}")
    face = m.sfr.budget.get(term="FLOW-JA-FACE", per=1)
    pair = face[face.reach.isin([0, 1]) & face.reach_to.isin([0, 1])]
    print("    sfr FLOW-JA-FACE reach/reach_to/q at period 1:",
          [(int(r.reach), int(r.reach_to), round(float(r.q), 2)) for r in pair.itertuples()])
    table = m.surface_water.results.q.get(per=3)
    row = table[(table.source == "sfr") & (table.cell == 195)].iloc[0]
    print(f"    surface_water cell 195 at period 3: q {row.q:.2f}, rlen {row.rlen:.2f},"
          f" exchange_intensity {row.exchange_intensity:.4f}")


def main(argv: list[str]) -> None:
    workspace = Path(argv[1]) if len(argv) > 1 else Path(tempfile.mkdtemp())
    try:
        canon = _canonical(workspace / "canon")
        section_1(canon)
        section_2_to_5(canon)
        section_6_to_9(canon)
        section_11(canon)
        section_10(_with_barriers(workspace / "hfb"))
    finally:
        if len(argv) <= 1:
            shutil.rmtree(workspace, ignore_errors=True)


if __name__ == "__main__":
    main(sys.argv)
