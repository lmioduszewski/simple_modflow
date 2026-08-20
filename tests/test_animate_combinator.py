"""`animate` as a Layer-2 combinator over pictures (plan 8.6a).

Before this stage `plot.animate(model)` returned a config bag -- an object with
`.sliders` and `.updatemenus` and nothing else. It was in `__all__`, documented
as a verb, and could not be shown, saved, or written to HTML. These tests pin
the shape that replaced it.

The heterogeneous-frames test is the one the plan names as the acceptance
criterion for the stage.
"""

from __future__ import annotations

import pytest

from myflopy import plot, viz


def test_animate_is_a_combinator_not_a_subject_verb():
    """`animate` takes FRAMES first, like `mosaic` -- not a subject like `map`.

    The model-bound `animate(model, periods=)` is gone: it had no callers
    outside its own docstring, and it could not answer any Picture verb.
    """

    import inspect

    parameters = list(inspect.signature(plot.animate).parameters)
    assert parameters[0] == "frames"
    assert "periods" not in parameters
    assert "backend" in parameters


def test_animate_takes_the_shapes_mosaic_takes():
    """Bare pictures or `(label, picture)` pairs, normalized identically."""

    from myflopy.modflow.mf6.interactive_plotting import _normalize_frames

    a, b = object(), object()
    assert _normalize_frames([("first", a), ("second", b)]) == [("first", a), ("second", b)]
    assert _normalize_frames([a, b]) == [("Frame 1", a), ("Frame 2", b)]
    with pytest.raises(ValueError, match="at least one frame"):
        _normalize_frames([])


def test_the_png_animation_has_no_single_figure():
    """Its frames are rendered images, so `.fig` raises and names the way out --
    the pattern `MplPicture` and `VtkScene` already use."""

    animation = plot.animate([object(), object()], backend="png")
    assert isinstance(animation, viz.Picture)
    with pytest.raises(TypeError, match=r"\.frames"):
        animation.fig


def test_an_animation_cannot_be_saved_as_a_still():
    """A still image holds one frame; saying so beats writing frame 1 silently."""

    animation = plot.animate([object()], backend="png")
    with pytest.raises(ValueError, match="one frame"):
        animation.save("out.png")


def test_animate_rejects_an_unknown_backend():
    with pytest.raises(ValueError, match="'plotly' or 'png'"):
        plot.animate([object()], backend="gif")


# --- the acceptance criterion -------------------------------------------------
@pytest.mark.canonical
@pytest.mark.slow
def test_both_backends_draw_the_canonical_model(canonical_run, tmp_path):
    """Imports prove nothing; these must actually write a file."""

    import matplotlib

    matplotlib.use("Agg")
    model = canonical_run
    count = min(3, len(model.kstpkper))
    frames = [(f"per {i}", model.plot.map(per=i, layer=0)) for i in range(count)]

    live = plot.animate(frames)
    assert isinstance(live, viz.Picture)
    assert len(live.fig.frames) == count
    assert live.labels == [f"per {i}" for i in range(count)]
    assert live.html(tmp_path / "live.html").stat().st_size > 0

    paged = plot.animate(frames, backend="png")
    written = paged.html(tmp_path / "paged.html")
    assert written.stat().st_size > 0
    assert "data:image/png;base64," in written.read_text(encoding="utf-8")

    # The richer handle is still reachable when you want the frame manifest.
    handle = paged.export(tmp_path / "handle.html")
    assert handle.frame_count == count


@pytest.mark.canonical
@pytest.mark.slow
def test_an_animation_can_mix_picture_kinds(canonical_run, tmp_path):
    """The plan's stated acceptance criterion for 8.6.

    A map and a cross-section share no trace structure, so no single plotly
    figure can hold both. Rasterizing each frame independently can -- which is
    why the png backend is the GENERAL one and plotly is the fast special case.
    """

    import matplotlib

    matplotlib.use("Agg")
    model = canonical_run
    mixed = [
        ("a map", model.plot.map(layer=0)),
        ("a section", model.plot.section(cells=[0, 1, 2])),
    ]

    written = plot.animate(mixed, backend="png").html(tmp_path / "mixed.html")
    assert written.read_text(encoding="utf-8").count("data:image/png;base64,") == 2

    # And the plotly backend must SAY it cannot, rather than fall back quietly:
    # swapping an interactive figure for a raster page changes what you get.
    with pytest.raises(ValueError, match="same trace structure"):
        plot.animate(mixed, backend="plotly").fig
