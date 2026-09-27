"""A field noun's `section` names what acts on it and refuses the rest (2026-09-26).

Ledger 315's "Not done", measured 2026-09-24: `model.hds.section(line=...,
fill="layer")` returned a plain head profile, byte-identical to the call without
`fill=`, and with `interpolator=` added it raised a `TypeError` from
`InterpolatedSurface.__init__`, a class the caller never named. The noun took
`**kwargs` and handed them to `XSection`, which keeps its tail only for the surface
it builds when interpolating: 24 of 28 keywords came back byte-identical
(`docs/reproducers/section_noun_kwargs_tail.py`).

The user's call: the noun stays its field's LINE PROFILE, and the cross-section's
arguments raise, naming `model.plot.section(fill=...)`, which draws it.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import io  # noqa: E402

import pytest  # noqa: E402
import shapely as shp  # noqa: E402

import myflopy as mf  # noqa: E402
from myflopy.project.model_group import ModelGroup  # noqa: E402

#: The model's cross-section's own arguments, each with a value it would act on.
CROSS_SECTION = {
    "fill": "layer", "fill_cmap": "magma", "fill_label": "head", "layers": [0],
    "head_layers": None, "layer_labels": ["unit"], "show_grid": False, "legend": False,
}


def _line(model):
    xmin, ymin, xmax, ymax = model.vor.gdf_vorPolys.total_bounds
    return shp.LineString([(xmin + 10, ymin + 10), (xmax - 10, ymax - 10)])


def _png(figure) -> bytes:
    import matplotlib.pyplot as plt

    buffer = io.BytesIO()
    figure.savefig(buffer, format="png")
    plt.close(figure)
    return buffer.getvalue()


@pytest.fixture(scope="module")
def sentinel_group(sentinel_run):
    """The sentinel run beside a reloaded copy of itself: a group and its diff."""

    return ModelGroup({"a": sentinel_run, "b": mf.load_mf6_run(sentinel_run.workspace)},
                      reference="a")


@pytest.mark.parametrize("fill", ["layer", "results"])
def test_the_noun_and_the_verb_agree_about_fill(fill, sentinel_run):
    """Asserted as AGREEMENT, not as one outcome.

    Whichever way the noun goes, it may not do what it did: take `fill=` and hand
    back a different picture. Refusing is agreement if the refusal names the verb
    that draws the picture; drawing is agreement if it is the verb's picture.
    """

    line = _line(sentinel_run)
    drawn_by_the_verb = _png(sentinel_run.plot.section(line=line, fill=fill, backend="mpl"))
    try:
        noun = sentinel_run.hds.section(line=line, fill=fill, backend="mpl")
    except TypeError as refusal:
        assert "model.plot.section(" in str(refusal), (
            f"the noun refused fill={fill!r} without naming the verb that draws it: {refusal}"
        )
    else:
        assert _png(noun) == drawn_by_the_verb, (
            f"the noun took fill={fill!r} and drew something other than the verb's picture"
        )


@pytest.mark.parametrize("typo", ["fil", "interpolater", "extent", "sectoin_name"])
def test_a_misspelled_keyword_raises_on_the_noun_and_the_verb(typo, sentinel_run):
    """Both pictures of the verb, and the noun, name the keyword they did not take."""

    line = _line(sentinel_run)
    calls = {
        "model.hds.section": lambda: sentinel_run.hds.section(line=line, **{typo: 1}),
        "model.plot.section": lambda: sentinel_run.plot.section(line=line, **{typo: 1}),
        "model.plot.section(fill=None)": lambda: sentinel_run.plot.section(
            line=line, fill=None, **{typo: 1}
        ),
    }
    for label, call in calls.items():
        with pytest.raises(TypeError) as refusal:
            call()
        assert typo in str(refusal.value), f"{label} did not name {typo!r}: {refusal.value}"


@pytest.mark.parametrize("name", CROSS_SECTION)
@pytest.mark.parametrize("scope", ["model", "group", "diff"])
def test_every_heads_noun_refuses_the_cross_sections_arguments(scope, name, sentinel_run,
                                                                sentinel_group):
    noun = {
        "model": sentinel_run.hds,
        "group": sentinel_group.hds,
        "diff": sentinel_group.diff().hds,
    }[scope]
    with pytest.raises(TypeError, match=r"model\.plot\.section\(") as refusal:
        noun.section(line=_line(sentinel_run), **{name: CROSS_SECTION[name]})
    assert f"{name}=" in str(refusal.value)


def test_section_name_is_refused_rather_than_crashing_a_group(sentinel_run, sentinel_group):
    """It collided with the name each member's line carries: `XSection() got
    multiple values for keyword argument 'section_name'`."""

    for noun in (sentinel_group.hds, sentinel_group.diff().hds):
        with pytest.raises(TypeError, match="unexpected keyword argument") as refusal:
            noun.section(line=_line(sentinel_run), section_name="x")
        assert "multiple values" not in str(refusal.value)


def test_surf_type_is_refused_for_what_it_is(sentinel_run):
    with pytest.raises(TypeError, match="surf_type= would section the layer elevations"):
        sentinel_run.hds.section(line=_line(sentinel_run), surf_type="lyr")


def _profile(model, **options):
    return model.hds.section(line=_line(model), **options).to_json()


@pytest.mark.parametrize(
    ("name", "value", "base"),
    [
        ("interpolate", True, {}),
        ("interpolator", "griddata", {}),
        ("use_rbf", False, {"interpolate": True}),
        ("num_points", 20, {"interpolator": "griddata"}),
        ("show_model_top", False, {}),
        ("show_model_btm", True, {}),
        ("layer", [0, 1], {}),
    ],
)
def test_what_the_noun_names_acts(name, value, base, canonical_run):
    """Each named option changes the picture -- nothing on the signature is inert."""

    assert _profile(canonical_run, **base, **{name: value}) != _profile(canonical_run, **base)


def test_kstpkper_and_clip_act_too(canonical_run):
    last = tuple(canonical_run.kstpkper[-1])
    assert _profile(canonical_run, kstpkper=last) != _profile(canonical_run)
    xmin, ymin, xmax, ymax = canonical_run.vor.gdf_vorPolys.total_bounds
    half = shp.box(xmin, ymin, (xmin + xmax) / 2, ymax)
    assert (_profile(canonical_run, interpolator="griddata", clip=half)
            != _profile(canonical_run, interpolator="griddata"))


def test_a_single_cell_path_may_be_an_int_and_x_or_y_turns_it(sentinel_run):
    """`XSection`'s own docstring showed `cells=69`, and both scopes document
    `int or list of int`; a bare int died at `len()`."""

    cell = int(sentinel_run.vor.gdf_vorPolys.index[0])
    as_int = sentinel_run.hds.section(cells=cell).to_json()
    assert as_int == sentinel_run.hds.section(cells=[cell]).to_json()
    assert as_int != sentinel_run.hds.section(cells=cell, x_or_y="y").to_json()
    profile = sentinel_run.plot.section(cells=cell, fill=None, backend="plotly")
    assert profile.fig.to_json() == sentinel_run.plot.section(
        cells=[cell], fill=None, backend="plotly"
    ).fig.to_json()


@pytest.mark.parametrize("verb", ["mosaic", "animate"])
@pytest.mark.parametrize(
    ("options", "message"),
    [
        ({"fill": "layer"}, r"model\.plot\.section\("),
        ({"fil": "layer"}, "fil"),
        ({"kstpkper": (0, 0)}, "per="),
    ],
)
def test_a_section_composer_refuses_what_the_noun_refuses(verb, options, message, sentinel_run):
    """`kind="section"` hands its tail to the same sections; it swallowed the same
    keywords, and `kstpkper` too, since each panel sets its own time."""

    with pytest.raises(TypeError, match=message):
        getattr(sentinel_run.hds, verb)(kind="section", line=_line(sentinel_run), **options)


def test_a_section_composer_still_takes_the_profiles_options(sentinel_run):
    line = _line(sentinel_run)
    plain = sentinel_run.hds.mosaic(kind="section", line=line).to_json()
    assert sentinel_run.hds.mosaic(kind="section", line=line, interpolate=True).to_json() != plain
