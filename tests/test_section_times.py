"""`per=N` is the last saved step of stress period N, on every picture (ledger 337).

The user's rule: "i want per to equal the last saved step of that stress period
(0-based numbering). so per=0 is the last saved step of the 1st stress period". It
held on the filled section and the map, and NOT on the line profile, the noun, the
group or ``surface``, which read ``per`` as an index into the SAVED records -- at
``nstp=3``, ``per=1`` drew (2, 1) on one picture and (1, 0) on the next. One helper
decides it now, ``period_end.last_saved_step``, and every reader defers to it.

The same pass: both sections open on the FIRST saved time (the filled one opened on
the last); ``per`` with ``kstpkper`` raises; the profile's interpolation-only
arguments raise without interpolation, naming ``interpolate=True``; ``spacing`` and
``extrapolate_beyond_section_ends`` are gone; ``surf_type="lyr"`` reads the model's
own discretization; the profile reads its model top from the model it draws,
without printing; and ``.ani`` opens on its first frame.

The fixture: ONE simulation holding two uncoupled GWF models on the same 92-cell
Voronoi grid, ``m`` (top 10,100 ft) declared before ``n`` (top 10,300 ft), both
transient over two periods of THREE steps, every step saved -- so a saved-record
index and a period differ, and ``sim.gwf[0]`` is not ``n``.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import geopandas as gpd  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402
from shapely.geometry import LineString, box  # noqa: E402

import myflopy as mf  # noqa: E402
from myflopy import plot  # noqa: E402
from myflopy.modflow.mf6.period_end import (  # noqa: E402
    last_saved_step,
    last_saved_steps,
    period_end_steps,
)
from myflopy.project.group._shared import _reduce_to_period_end  # noqa: E402

SIDE = 600.0
LINE = LineString([(10.0, 10.0), (590.0, 590.0)])
TIMES = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1)]


@pytest.fixture(scope="module")
def two_models(tmp_path_factory):
    """``{"m": model, "n": model}`` from one run: 2 periods x 3 saved steps each."""

    root = tmp_path_factory.mktemp("section_times")
    boundary = root / "boundary.gpkg"
    gpd.GeoDataFrame({"n": ["d"]}, geometry=[box(0, 0, SIDE, SIDE)], crs=2927).to_file(boundary)
    vor = mf.GridSpec.voronoi(
        boundary=mf.GeoPackageSourceSpec(path=boundary), max_area=12_000.0
    ).resolve(workspace=root / "grid")
    x = np.asarray(vor.centroids_x, dtype=float)
    edges = vor.gdf_vorPolys.geometry.bounds
    west = np.flatnonzero(edges.minx <= 1e-6)
    east = np.flatnonzero(edges.maxx >= SIDE - 1e-6)
    props = vor.get_disv_gridprops()
    assert len(x) == props["ncpl"]

    def chd(w, e):
        return [[(0, int(c)), w] for c in west] + [[(0, int(c)), e] for c in east]

    def gwf(name, top):
        return mf.gwf(name, grid=vor, packages=[
            mf.disv(nlay=1, ncpl=props["ncpl"], nvert=len(props["vertices"]),
                    vertices=props["vertices"], cell2d=props["cell2d"],
                    top=top, botm=[9_900.0]),
            mf.npf(k=10.0, icelltype=1), mf.ic(strt=10_000.0),
            mf.sto(iconvert=1, sy=0.1, ss=1e-5, transient={0: True}),
            mf.chd(stress_period_data={0: chd(9_990.0, 10_010.0), 1: chd(9_995.0, 10_005.0)}),
            mf.oc(head_filerecord=f"{name}.hds", saverecord=[("HEAD", "ALL")]),
        ])

    simulation = mf.simulation(gwf("m", 10_100.0), gwf("n", 10_300.0), name="s", packages=[
        mf.tdis(nper=2, perioddata=[(3.0, 3, 1.0), (3.0, 3, 1.0)]), mf.ims(),
    ])
    run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    models = {"m": run.model("m"), "n": run.model("n")}
    assert [tuple(t) for t in models["m"].kstpkper] == TIMES
    return models


def _json(picture) -> str:
    return (picture.fig if hasattr(picture, "ani") or hasattr(picture, "_assemble") else picture).to_json()


def _profile(model, **kwargs):
    return plot.section(model, LINE, fill=None, backend="plotly", **kwargs)


# -- the one rule ---------------------------------------------------------------


def test_the_helper_is_the_last_saved_step_of_each_period():
    saved = [(0, 0), (2, 0), (1, 0), (0, 1), (0, 3), (4, 3)]
    assert last_saved_steps(saved) == {0: (2, 0), 1: (0, 1), 3: (4, 3)}
    assert last_saved_step(saved, 0) == (2, 0)
    assert last_saved_step(saved, 3) == (4, 3)
    for missing in (2, 4, -1):
        with pytest.raises(ValueError, match=r"Available periods: \[0, 1, 3\]"):
            last_saved_step(saved, missing)
    # ledger 336's period-end rule is built on it, and its no-perlen route is it
    assert period_end_steps(saved) == last_saved_steps(saved)


def test_a_group_reduces_each_model_to_its_own_last_saved_step():
    frame = pd.DataFrame({
        "model": ["a", "a", "a", "a", "b", "b"],
        "cell": [0, 0, 1, 1, 0, 0],
        "kstpkper": [(0, 0), (2, 0), (0, 0), (1, 0), (0, 0), (1, 0)],
        "value": [1.0, 3.0, 5.0, 6.0, 7.0, 8.0],
    })
    reduced = _reduce_to_period_end(frame)
    # model a's last saved step is (2, 0): cell 1 has no value there, and an
    # earlier step's value is a different number; model b's own last is (1, 0)
    assert reduced[["model", "cell", "value"]].values.tolist() == [["a", 0, 3.0], ["b", 0, 8.0]]


# -- per is a period on every picture -----------------------------------------------


@pytest.mark.parametrize("per", [0, 1])
def test_per_draws_the_periods_last_saved_step_on_the_profile(two_models, per):
    model = two_models["m"]
    expected = last_saved_step(TIMES, per)
    assert _profile(model, per=per).kstpkper == expected
    assert _json(_profile(model, per=per)) == _json(_profile(model, kstpkper=expected))


@pytest.mark.parametrize("per", [0, 1])
def test_per_draws_the_periods_last_saved_step_on_the_noun(two_models, per):
    model = two_models["m"]
    expected = last_saved_step(TIMES, per)
    assert (model.hds.section(line=LINE, per=per).to_json()
            == model.hds.section(line=LINE, kstpkper=expected).to_json())
    assert model.hds.section(line=LINE, per=per).to_json() != (
        model.hds.section(line=LINE, kstpkper=TIMES[per]).to_json()
    ), "the fixture cannot tell a period from a saved-record index"


def test_per_agrees_across_the_pictures(two_models):
    model = two_models["m"]
    filled = plot.section(model, LINE, backend="plotly", per=1)
    assert filled.fig.to_json() == plot.section(
        model, LINE, backend="plotly", kstpkper=(2, 1)).fig.to_json()
    assert model.plot.map(per=1).kstpkper == (2, 1)
    surface = plot.surface(model, per=1, resolution=40)
    assert tuple(surface.kstpkper) == (2, 1)


def test_per_on_a_group_is_each_members_period(two_models):
    from myflopy.project.model_group import ModelGroup

    group = ModelGroup(dict(two_models), reference="m")
    table = group.hds.get(per=1)
    assert set(table["kstpkper"]) == {(2, 1)}
    assert (group.hds.section(line=LINE, per=1).to_json()
            == group.hds.section(line=LINE, kstpkper=(2, 1)).to_json())


def test_a_period_the_output_did_not_save_raises(two_models):
    with pytest.raises(ValueError, match=r"Stress period 2 is unavailable"):
        _profile(two_models["m"], per=2)
    with pytest.raises(ValueError, match=r"kstpkper \(5, 0\) is unavailable"):
        _profile(two_models["m"], kstpkper=(5, 0))


@pytest.mark.parametrize("draw", [
    lambda model: plot.section(model, LINE, per=1, kstpkper=(0, 0)),
    lambda model: plot.section(model, LINE, fill=None, per=1, kstpkper=(0, 0)),
    lambda model: model.hds.section(line=LINE, per=1, kstpkper=(0, 0)),
    lambda model: plot.surface(model, per=1, kstpkper=(0, 0), resolution=40),
], ids=["cross-section", "profile", "noun", "surface"])
def test_per_with_kstpkper_raises(two_models, draw):
    # Documented as exclusive; the profile let per win and the cross-section kstpkper.
    with pytest.raises(ValueError, match="per= and kstpkper= both choose the time"):
        draw(two_models["m"])


# -- both sections open on the first saved time -----------------------------------


@pytest.mark.parametrize("fill", ["layer", "results"])
def test_the_cross_section_opens_on_the_first_saved_time(two_models, fill):
    model = two_models["m"]
    default = plot.section(model, LINE, fill=fill, backend="plotly")
    first = plot.section(model, LINE, fill=fill, backend="plotly", kstpkper=TIMES[0])
    last = plot.section(model, LINE, fill=fill, backend="plotly", kstpkper=TIMES[-1])
    assert default.fig.to_json() == first.fig.to_json()
    assert default.fig.to_json() != last.fig.to_json()


def test_the_profile_opens_on_the_first_saved_time(two_models):
    model = two_models["m"]
    assert _profile(model).kstpkper == TIMES[0]


def test_ani_opens_on_its_first_frame(two_models):
    fig = _profile(two_models["m"]).ani
    assert fig.layout.sliders[0].active in (0, None)
    still = {trace.name: list(trace.y) for trace in fig.data}
    first = {trace.name: list(trace.y) for trace in fig.frames[0].data}
    head = "Lyr 0 hds - m"
    assert still[head] == first[head]
    assert still[head] != {t.name: list(t.y) for t in fig.frames[-1].data}[head]


# -- the profile's own arguments ---------------------------------------------------


@pytest.mark.parametrize("name, value", [
    ("num_points", 20), ("use_rbf", False), ("clip", box(0, 0, 300, 600)),
    ("surf_type", "lyr"), ("resolution", 40),
])
def test_an_interpolation_argument_without_interpolation_raises(two_models, name, value):
    model = two_models["m"]
    with pytest.raises(ValueError, match=rf"{name} shape the INTERPOLATED profile.*interpolate=True"):
        _profile(model, **{name: value})
    with pytest.raises(ValueError, match="interpolate=True"):
        plot.section(model, LINE, fill=None, backend="mpl", **{name: value})
    if name in ("num_points", "use_rbf", "clip"):    # what the noun names
        with pytest.raises(ValueError, match="interpolate=True"):
            model.hds.section(line=LINE, **{name: value})
    # ...and it is taken when interpolating
    kwargs = {"resolution": 40} if name != "resolution" else {}
    _profile(model, interpolate=True, **{name: value}, **kwargs).fig


def test_an_interpolated_profile_refuses_several_layers(two_models):
    with pytest.raises(ValueError, match="draws ONE layer's surface"):
        _profile(two_models["m"], interpolate=True, layer=[0, 0], resolution=40)


@pytest.mark.parametrize("name, value", [
    ("spacing", 50), ("extrapolate_beyond_section_ends", True),
])
@pytest.mark.parametrize("fill", [None, "layer"])
def test_the_removed_arguments_raise(two_models, name, value, fill):
    # Nothing ever read `spacing`; `extrapolate_beyond_section_ends` only raised.
    with pytest.raises(TypeError, match=name):
        plot.section(two_models["m"], LINE, fill=fill, **{name: value})
    with pytest.raises(TypeError, match=name):
        two_models["m"].hds.section(line=LINE, **{name: value})


def test_surf_type_lyr_draws_the_layer_top_from_the_models_discretization(two_models):
    model = two_models["m"]
    # A hand-written DISV's grid carries no gdf_topbtm; this was the crash.
    assert model.vor.gdf_topbtm is None
    picture = _profile(model, interpolate=True, surf_type="lyr", resolution=40)
    traces = {trace.name: np.asarray(trace.y, dtype=float) for trace in picture.fig.data}
    assert set(traces) == {"Lyr 0 top - m", "model top"}
    assert np.allclose(traces["Lyr 0 top - m"], 10_100.0)
    surface = plot.surface(model, surf_type="lyr", resolution=40)
    assert np.nanmax(np.abs(surface.clipped_surface - 10_100.0)) < 1e-6
    with pytest.raises(ValueError, match="surf_type must be"):
        plot.surface(model, surf_type="layer", resolution=40)


# -- the model top: the drawn model's own, read without printing ---------------------


def test_the_profile_draws_its_own_models_top(two_models):
    # `sim.gwf[0]` is m: n's profile drew m's 10,100 ft top.
    for name, top in (("m", 10_100.0), ("n", 10_300.0)):
        traces = {t.name: np.asarray(t.y, dtype=float) for t in _profile(two_models[name]).fig.data}
        assert np.allclose(traces["model top"], top), name


@pytest.mark.parametrize("draw", [
    lambda model: _profile(model).fig,
    lambda model: plot.section(model, LINE, fill=None),
    lambda model: model.hds.section(line=LINE),
    lambda model: model.hds.mosaic(kind="section", line=LINE),
], ids=["verb-plotly", "verb-mpl", "noun", "composer"])
def test_a_line_profile_prints_nothing(two_models, draw, capfd):
    import matplotlib.pyplot as plt

    capfd.readouterr()
    draw(two_models["m"])
    plt.close("all")
    out, _ = capfd.readouterr()
    assert out == "", f"printed: {out.splitlines()[:3]}"
