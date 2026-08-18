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
