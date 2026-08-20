"""One contract for every drawable thing (plan 8.1).

Before this, four picture classes had four different names for their figure --
`Choro.choropleth`, `XSection.fig`, `GridSection.figure`, and
`InterpolatedSurface` with none at all -- and `.plot()` meant three incompatible
things depending on which one you held:

* `Choro.plot()`   returned the figure
* `GridSection.plot()` showed it and returned None
* `InterpolatedSurface.plot()` opened a BROWSER WINDOW and returned None

So you learned each class separately, and `model.cor().plot().show()` was the
price of not knowing which. Now every picture answers `.fig`, `.show()`,
`.save()`, `.html()`, and renders itself in a notebook.
"""

from __future__ import annotations

import pytest

from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface
from myflopy.modflow.mf6.grid.plotting import GridSection
from myflopy.modflow.utils.datatypes.choros import Choro
from myflopy.modflow.utils.datatypes.xsections import XSection
from myflopy.viz import Fig, Picture

PICTURES = [Choro, XSection, GridSection, InterpolatedSurface]
IDS = [cls.__name__ for cls in PICTURES]


@pytest.mark.parametrize("cls", PICTURES, ids=IDS)
def test_every_picture_answers_the_contract(cls):
    assert issubclass(cls, Picture)
    for verb in ("fig", "show", "save", "html", "_repr_mimebundle_"):
        assert hasattr(cls, verb), f"{cls.__name__} is missing {verb}"


@pytest.mark.parametrize("cls", PICTURES, ids=IDS)
def test_the_old_spellings_are_gone(cls):
    """`.plot()` especially: on a picture it meant "give me the figure", while in
    the view grammar it means TIME SERIES. One name, two meanings, and the
    grammar's is the one worth keeping."""

    for dead in ("plot", "choropleth", "figure", "get_figure"):
        assert not hasattr(cls, dead), (
            f"{cls.__name__}.{dead} still exists; the contract is `.fig`"
        )


def test_get_choropleth_survives_because_it_returns_a_trace():
    """Not a fourth spelling of the figure -- it builds the choropleth TRACE, and
    `viz.mosaic` composes panels out of exactly that. Deleting it as a duplicate
    (which this plan's first draft nearly did) would have broken mosaics."""

    assert hasattr(Choro, "get_choropleth")


# --- the bug the contract had to fix to exist -------------------------------
class _TwoTraceStub(Picture):
    """A picture whose assembly appends, like `Choro.add_choropleth` does."""

    def __init__(self):
        self._fig = Fig()
        self._assembled = False

    @property
    def fig(self):
        if not self._assembled:
            self._fig.add_scatter(x=[0, 1], y=[0, 1])
            self._assembled = True
        return self._fig


def test_assembling_twice_does_not_draw_twice():
    """`Choro.add_choropleth()` calls `fig.add_trace(...)` unconditionally, so
    the old `.choropleth` property was NOT idempotent -- and `.plot()` just
    returned it. Touching the figure twice silently drew every trace twice.
    The contract requires assemble-once, which is what makes
    `picture.fig.update_layout(...)` followed by `picture.show()` safe.
    """

    picture = _TwoTraceStub()
    assert len(picture.fig.data) == 1
    assert len(picture.fig.data) == 1          # the second touch must be free
    picture.fig.update_layout(title="edited")
    assert picture.fig.layout.title.text == "edited"
    assert len(picture.fig.data) == 1


def test_choro_assembly_flag_survives_object_new():
    """`Choro` is built via `object.__new__` by several test doubles, so the
    assembly flag is a CLASS attribute -- a `.fig` that raises AttributeError on
    those is a contract that only half-holds."""

    bare = object.__new__(Choro)
    assert bare._assembled is False


def test_save_names_kaleido_rather_than_leaking_plotlys_error(tmp_path, monkeypatch):
    """`kaleido` is not a declared dependency, so static export can fail on a
    clean install. The message must name it, and point at the .html route that
    needs nothing extra."""

    import builtins

    real_import = builtins.__import__

    def _no_kaleido(name, *args, **kwargs):
        if name == "kaleido":
            raise ImportError("No module named 'kaleido'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _no_kaleido)
    with pytest.raises(ImportError, match="kaleido is required"):
        _TwoTraceStub().save(tmp_path / "x.png")


def test_save_dispatches_html_without_needing_kaleido(tmp_path):
    out = _TwoTraceStub().save(tmp_path / "picture.html")
    assert out.exists() and "plotly" in out.read_text(encoding="utf-8").lower()


def test_html_writes_a_standalone_file(tmp_path):
    out = _TwoTraceStub().html(tmp_path / "standalone.html")
    assert out.exists() and out.stat().st_size > 0


def test_a_picture_without_fig_says_so():
    class _Undefined(Picture):
        pass

    with pytest.raises(NotImplementedError, match="does not define `fig`"):
        _Undefined().show()


@pytest.mark.canonical
@pytest.mark.slow
def test_an_animation_is_the_choros_figure_from_then_on(canonical_run):
    """`.ani` swaps in the animation figure; the static traces must not come back.

    `.ani` builds a frames figure and assigns it to `_fig`. It used to leave
    `_assembled` False, so the next `.fig` access -- which `.show()`, `.html()`
    and `.save()` ALL go through -- re-ran the assembly and appended the static
    choropleth, contours, locs and overlays on top of the animation. The bug was
    invisible to the shipping exporters because they read `.ani` and never touch
    `.fig`.
    """

    choro = canonical_run.plot.map(layer=0)
    animation = choro.ani

    assert animation.frames, "no frames on the animation figure"
    before = len(animation.data)

    # The Picture output path. Idempotent, and still the animation.
    assert choro.fig is animation
    assert len(choro.fig.data) == before
    assert choro.fig.frames, "accessing .fig dropped the frames"


# --- 8.5a: a Picture whose renderer is not Plotly -----------------------------
def _demo_mpl_picture():
    """A minimal MplPicture, so these test the CONTRACT and not layers.py."""

    from myflopy.viz import MplPicture, mpl_axes

    class Demo(MplPicture):
        title = "demo"

        def __init__(self):
            self.draws = 0

        def draw(self, ax=None, **kwargs):
            self.draws += 1
            if ax is None:
                _, ax = mpl_axes()
            ax.plot([0, 1], [0, 1])
            return ax

    return Demo()


def test_an_mpl_picture_answers_the_same_verbs(tmp_path):
    """The point of `MplPicture`: a filled cross-section has no Plotly form, and
    exempting it from the grammar would mean callers learning which pictures are
    "real" ones. It answers show/save/html/inline like everything else."""

    import matplotlib

    matplotlib.use("Agg")
    picture = _demo_mpl_picture()

    assert picture.save(tmp_path / "a.png").exists()
    assert picture.save(tmp_path / "a.pdf").exists()
    assert "image/png" in picture._repr_mimebundle_()


def test_an_mpl_picture_draws_once(tmp_path):
    """`Picture` requires idempotence -- the defect 8.1 and ledger 130 both fixed
    for Plotly. The Axes is cached, so repeated output does not redraw."""

    import matplotlib

    matplotlib.use("Agg")
    picture = _demo_mpl_picture()

    picture.axes
    picture.axes
    picture.save(tmp_path / "a.png")
    assert picture.draws == 1


def test_an_mpl_pictures_html_is_self_contained(tmp_path):
    """No plotly.js to fetch, because there is no Plotly figure. This is the
    artifact you email someone -- it must not need a network."""

    import matplotlib

    matplotlib.use("Agg")
    out = _demo_mpl_picture().html(tmp_path / "a.html")
    text = out.read_text(encoding="utf-8")

    assert "data:image/png;base64," in text
    assert "http://" not in text and "https://" not in text


def test_asking_an_mpl_picture_for_fig_says_what_to_use_instead():
    """`fig` is documented package-wide as the PLOTLY figure. Returning an
    `mpl.Figure` would satisfy the letter and break every caller reaching for
    `.add_trace`/`.update_layout`, so it raises and names `.axes`."""

    import pytest

    picture = _demo_mpl_picture()
    with pytest.raises(TypeError, match=r"\.axes"):
        picture.fig


@pytest.mark.slow
def test_the_stack_verbs_are_all_pictures(tmp_path):
    """Imports prove nothing; each verb must actually draw."""

    import matplotlib

    matplotlib.use("Agg")
    import myflopy as mf
    from myflopy.layers import Flat, LayerStack
    from myflopy.viz import Fig, MplPicture, Picture

    tri = mf.TriangleGrid(model_ws=str(tmp_path), angle=30)
    tri.set_domain_rectangle(x_dist=400, y_dist=300, origin=(0, 0))
    tri.build(verbose=False)
    vor = mf.VoronoiGridPlus(tri)
    result = LayerStack(vor, top=Flat(50)).add("a", bottom=Flat(20)).build()

    thickness = result.plot.map()
    assert isinstance(thickness, MplPicture)
    assert hasattr(thickness.axes, "set_title")

    section = result.plot.section(y=150)
    assert isinstance(section, MplPicture)
    assert hasattr(section.axes, "set_title")

    surface = result.plot.surface("top", resolution=20)
    assert isinstance(surface, Picture)
    assert isinstance(surface.fig, Fig)   # a house Fig, not a bare go.Figure


# --- the contract, applied to the REAL classes (8.6a) --------------------------
@pytest.mark.canonical
@pytest.mark.slow
def test_every_real_picture_is_idempotent(canonical_run):
    """`test_assembling_twice_does_not_draw_twice` uses a STUB, so it proved the
    contract about a class written to satisfy it -- and every real subclass went
    unchecked. `XSection.fig` rebuilt from scratch on every access for exactly
    that reason: repeated access returned equivalent but DIFFERENT figures, so
    `xs.fig.update_layout(...)` then `xs.show()` silently dropped the edit.

    Parametrized over the live classes so a new Picture cannot quietly opt out.
    """

    import matplotlib
    from shapely.geometry import LineString

    matplotlib.use("Agg")
    model, vor = canonical_run, canonical_run.vor
    xmin, ymin, xmax, ymax = vor.gdf_vorPolys.total_bounds
    mid = (ymin + ymax) / 2
    pictures = {
        "Choro": model.plot.map(layer=0),
        "XSection": model.plot.section(cells=[0, 1, 2]),
        "GridMesh": vor.plot.grid(),
        "GridSection": vor.plot.section(LineString([(xmin, mid), (xmax, mid)])),
        "InterpolatedSurface": model.plot.surface(layer=0),
    }

    # Every PLOTLY Picture subclass must be represented -- the point is that a
    # new one cannot quietly opt out. MplPicture/VtkScene are the two whose
    # `.fig` deliberately raises, and LayerSurface needs a built stack.
    plotly_pictures = {
        cls.__name__
        for cls in Picture.__subclasses__()
        if cls.__module__.startswith("myflopy.")          # not this file's stubs
        and cls.__name__ not in {"MplPicture", "VtkScene", "LayerSurface"}
    }
    assert plotly_pictures <= set(pictures), (
        f"unchecked Picture subclasses: {sorted(plotly_pictures - set(pictures))}"
    )

    for name, picture in pictures.items():
        first = picture.fig
        assert picture.fig is first, f"{name}.fig rebuilt on second access"
        before = len(first.data)
        picture.fig.update_layout(title="edited")
        assert picture.fig.layout.title.text == "edited", f"{name} lost a layout edit"
        assert len(picture.fig.data) == before, f"{name} drew its traces twice"


@pytest.mark.canonical
@pytest.mark.slow
def test_a_section_animation_is_the_sections_figure_from_then_on(canonical_run):
    """The `XSection` half of ledger 130.

    `Choro.ani` was fixed to persist its animation figure; `XSection.ani` built
    one, returned it, and never assigned it -- so `plot.section(...).ani` then
    `.show()` or `.html()` wrote the STATIC section. One sibling was fixed and
    the other was not.
    """

    section = canonical_run.plot.section(cells=[0, 1, 2])
    animation = section.ani

    assert animation.frames, "no frames on the section animation"
    assert section.fig is animation
    assert section.fig.frames, "accessing .fig dropped the frames"
