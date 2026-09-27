"""`model.hds.section(fill="layer")` returned a plain head profile, byte-identical.

Found 2026-09-24 while restating ledger 315 (recorded there under "Not done").

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/section_noun_kwargs_tail.py

THE MECHANISM. The noun's `section` (`SpatialView.section`) took `**kwargs` and handed
them to `_sections`, which built `XSection(model=, line=, cells=, per=, layer=, **kwargs)`.
`XSection` keeps its own `**kwargs` ONLY to pass to `InterpolatedSurface`, and builds one
of those only when interpolating. So a keyword the noun did not name went one of three
ways, and nothing told the caller which:

  * it was one of `XSection`'s own parameters, and it reached the class (it may act);
  * it was not, and there was no `interpolator=`: stored, never read -- SWALLOWED;
  * it was not, and there was one: `TypeError` from `InterpolatedSurface.__init__`, a
    class the caller never named.

THE SETUP. The model `tests/conftest.py::sentinel_run` builds, with one change: TWO
periods instead of one, so that a keyword choosing a TIME (`kstpkper=`) can show whether
it acts. Real `GridSpec.voronoi` grid on a 600 ft square (92 cells), one layer,
`icelltype=1`, no Newton, steady; CHD 9,990 ft west / 10,010 ft east in period 0 and
9,995 / 10,005 in period 1; the north-east corner is DRY and an interior block is
idomain 0. The section runs corner to corner. Each keyword is passed ALONE on top of the
same base call, with a value that differs from its default, and the result is compared
with the base call's:

  plotly  `fig.to_json()`, byte for byte;
  mpl     the Agg PNG's SHA-256.

  SAME     the picture is byte-identical to the base -- the keyword did nothing;
  changed  it acted;
  raises   the first line of the exception.

Run it at a69ac87 for the defect, and at the fixing commit (ledger 334) for the fix.
Section 3 names `fill=None` and `backend=` on the verb, so it draws the same picture --
the line profile -- at both commits (the verb's defaults moved in ledger 333).

WHAT IT FOUND at a69ac87: 24 of the 28 keywords came back byte-identical from the noun,
on both backends and on the group and diff nouns; with an interpolator, every keyword
`XSection` does not name became a `TypeError` from `InterpolatedSurface.__init__`;
`section_name=` crashed the group and diff nouns ("got multiple values"); the section
composers swallowed the same set, `kstpkper` included; and the verb swallowed every
typo too. After the fix every keyword the noun does not name raises a `TypeError` that
names it, and the cross-section's own point at `model.plot.section(fill=...)`.
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import pathlib
import tempfile
import warnings

import geopandas as gpd
import matplotlib
import numpy as np
from shapely.geometry import LineString, box

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

import myflopy as mf  # noqa: E402
from myflopy.project.model_group import ModelGroup  # noqa: E402
from myflopy.viz import Picture  # noqa: E402

SIDE = 600.0
LINE = LineString([(10.0, 10.0), (590.0, 590.0)])
HALF = box(0.0, 0.0, 300.0, 600.0)

#: One non-default value per keyword the noun's tail could reach: every parameter
#: the free verb names and the noun did not, the names `InterpolatedSurface` takes
#: that `XSection` does not already pass, a map option, and two typos.
KEYWORDS = {
    # the free verb's line-profile parameters
    "kstpkper": (0, 1),                 # the base reads the FIRST saved time, (0, 0)
    "x_or_y": "y",
    "spacing": 50,
    "num_points": 20,
    "interpolate": True,
    "use_rbf": False,
    "interpolator": "griddata",
    "extrapolate_beyond_section_ends": True,
    "show_model_top": False,
    "show_model_btm": True,
    "surf_type": "lyr",
    "section_name": "renamed",
    "clip": HALF,
    "animation_kstpkpers": [(0, 0)],
    # the free verb's filled-section parameters
    "fill": "layer",
    "fill_cmap": "magma",
    "fill_label": "label",
    "layers": [0],
    "head_layers": None,
    "layer_labels": ["sand"],
    "show_grid": False,
    "legend": False,
    "reverse": True,
    # reach `InterpolatedSurface` through `XSection._kwargs`
    "resolution": 50,
    "crs": "EPSG:2927",
    # a map option, and two typos
    "extent": (0.0, 0.0, 300.0, 300.0),
    "interpolater": "griddata",
    "fil": "layer",
}


def build(root: pathlib.Path):
    """The sentinel model over two periods (see the module docstring)."""

    boundary = root / "boundary.gpkg"
    gpd.GeoDataFrame({"n": ["d"]}, geometry=[box(0, 0, SIDE, SIDE)], crs=2927).to_file(boundary)
    vor = mf.GridSpec.voronoi(
        boundary=mf.GeoPackageSourceSpec(path=boundary), max_area=12_000.0
    ).resolve(workspace=root / "grid")
    x = np.asarray(vor.centroids_x, dtype=float)
    y = np.asarray(vor.centroids_y, dtype=float)
    raised = (x > 400) & (y > 400)
    botm = np.where(raised, 10_050.0, 9_900.0)
    idomain = np.where((x > 225) & (x < 375) & (y > 225) & (y < 375), 0, 1)
    edges = vor.gdf_vorPolys.geometry.bounds
    west = np.flatnonzero(edges.minx <= 1e-6)
    east = [c for c in np.flatnonzero(edges.maxx >= SIDE - 1e-6) if not raised[c]]

    def chd(w, e):
        return [[(0, int(c)), w] for c in west] + [[(0, int(c)), e] for c in east]

    props = vor.get_disv_gridprops()
    model = mf.gwf("m", grid=vor, packages=[
        mf.disv(nlay=1, ncpl=props["ncpl"], nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"],
                top=10_100.0, botm=[botm], idomain=[idomain]),
        mf.npf(k=10.0, icelltype=1), mf.ic(strt=10_000.0),
        mf.chd(stress_period_data={0: chd(9_990.0, 10_010.0), 1: chd(9_995.0, 10_005.0)}),
        mf.oc(head_filerecord="m.hds", saverecord=[("HEAD", "ALL")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=2, perioddata=[(1.0, 1, 1.0), (1.0, 1, 1.0)]),
        mf.ims(models=["m"]),
    ])
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "Run 'r' built but its frozen snapshot", UserWarning)
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    return run.model()


def fingerprint(result) -> str:
    """What the caller SEES: plotly's JSON, or the rendered PNG."""

    if isinstance(result, Picture):
        # The verb's XSection. Not `hasattr(result, "fig")`: building `.fig` can
        # raise AttributeError, which `hasattr` would swallow as "no such attribute".
        result = result.fig
    if hasattr(result, "to_json"):
        return hashlib.sha256(result.to_json().encode()).hexdigest()
    buffer = io.BytesIO()
    result.savefig(buffer, format="png")
    plt.close(result)
    return hashlib.sha256(buffer.getvalue()).hexdigest()


def outcome(call, base: str | None) -> str:
    """``SAME`` / ``changed`` / ``raises <first line>``, against the base call."""

    try:
        # `XSection.model_top` reloads the simulation through FloPy at its
        # default verbosity, which prints; that is not what is measured here.
        with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
            warnings.simplefilter("ignore")
            print_ = fingerprint(call())
    except Exception as error:  # noqa: BLE001 - the report IS every failure mode
        head = str(error).splitlines()[0] if str(error) else ""
        return f"raises {type(error).__name__}: {head[:90]}"
    if base is None:
        return "ok"
    return "SAME" if print_ == base else "changed"


def sweep(label: str, draw) -> None:
    """One column of the table: every keyword through ``draw(**kw)``."""

    print(f"\n{label}")
    base_outcome = outcome(draw, None)
    if base_outcome != "ok":                   # e.g. a keyword the base call relies on
        print(f"  the base call itself {base_outcome}")
        return
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        base = fingerprint(draw())
    for name, value in KEYWORDS.items():
        print(f"  {name:32s} {outcome(lambda: draw(**{name: value}), base)}")


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        model = build(pathlib.Path(tmp))
        group = ModelGroup({"a": model, "b": mf.load_mf6_run(model.workspace)}, reference="a")
        diff = group.diff()

        print("1. THE DEFECT")
        with contextlib.redirect_stdout(io.StringIO()):
            plain = model.hds.section(line=LINE)
            try:
                filled = model.hds.section(line=LINE, fill="layer").to_json()
            except (TypeError, ValueError) as error:
                filled = f"raises {type(error).__name__}: {error}".splitlines()[0]
        print("   hds.section(line=L).to_json() == hds.section(line=L, fill='layer').to_json():",
              plain.to_json() == filled if filled.startswith("{") else filled)
        print("   hds.section(line=L, fill='layer', interpolator='griddata'):",
              outcome(lambda: model.hds.section(line=LINE, fill="layer",
                                                interpolator="griddata"), None))
        print("   plot.section(line=L, fill='layer', backend='mpl'):",
              outcome(lambda: model.plot.section(line=LINE, fill="layer", backend="mpl"), None))
        print("   hds.section(line=L, fill='layer', backend='mpl'):",
              outcome(lambda: model.hds.section(line=LINE, fill="layer", backend="mpl"), None))

        print("\n2. EVERY KEYWORD THE TAIL REACHES, one at a time")
        sweep("model.hds.section(line=L), plotly",
              lambda **kw: model.hds.section(line=LINE, **kw))
        sweep("model.hds.section(line=L), mpl",
              lambda **kw: model.hds.section(line=LINE, backend="mpl", **kw))
        sweep("model.hds.section(line=L, interpolator='griddata'), plotly",
              lambda **kw: model.hds.section(line=LINE, **{"interpolator": "griddata", **kw}))
        sweep("group.hds.section(line=L), plotly",
              lambda **kw: group.hds.section(line=LINE, **kw))
        sweep("diff.hds.section(line=L), plotly",
              lambda **kw: diff.hds.section(line=LINE, **kw))
        sweep("model.hds.mosaic(kind='section', line=L), plotly",
              lambda **kw: model.hds.mosaic(kind="section", line=LINE, **kw))
        sweep("model.hds.animate(kind='section', line=L), plotly",
              lambda **kw: model.hds.animate(kind="section", line=LINE, **kw))

        print("\n3. THE VERB's line profile, for comparison")
        sweep("model.plot.section(line=L, fill=None), plotly",
              lambda **kw: model.plot.section(
                  line=LINE, **{"fill": None, "backend": "plotly", **kw}))
        sweep("model.plot.section(line=L, fill=None), mpl",
              lambda **kw: model.plot.section(
                  line=LINE, **{"fill": None, "backend": "mpl", **kw}))

        print("\n4. A SINGLE-CELL PATH, where x_or_y is documented to act")
        cell = int(model.vor.gdf_vorPolys.index[0])
        for scope, draw in (
            ("noun", lambda **kw: model.hds.section(**kw)),
            ("verb", lambda **kw: model.plot.section(fill=None, backend="plotly", **kw)),
        ):
            # `XSection`'s docstring shows `cells=69`; its constructor skips the
            # setter that wraps an int, so a bare int fails at `len()`.
            print(f"   {scope} cells={cell}: {outcome(lambda: draw(cells=cell), None)}")
            with contextlib.redirect_stdout(io.StringIO()):
                base = fingerprint(draw(cells=[cell]))
            print(f"   {scope} cells=[{cell}], x_or_y='y': "
                  f"{outcome(lambda: draw(cells=[cell], x_or_y='y'), base)}")


if __name__ == "__main__":
    main()
