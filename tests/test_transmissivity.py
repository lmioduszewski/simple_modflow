"""Transmissivity maps and nouns (ledger 308).

`model.plot.map(type="trans")` / `type="trans_max"` and
`model.packages.npf.transmissivity` / `.transmissivity_max` -- one computation,
`hydraulics.model_transmissivity`, so the four cannot disagree.

**Saturated is the default** (the user's decision): K times the thickness MF6
actually conducts through at a time. The rule is MF6's ICELLTYPE rule, not "head
minus bottom" everywhere, and the headline test holds it to MF6's own `DATA-SAT`
at every saved time.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from myflopy import plot
from myflopy.modflow.mf6.hydraulics import (
    flow_model,
    model_transmissivity,
    transmissivity,
    unit_labels,
)


@pytest.fixture(scope="module")
def unrun_model(tmp_path_factory):
    from myflopy.modflow.mf6.canonical_example import (
        CanonicalModelConfig,
        build_canonical_model,
    )

    return build_canonical_model(
        tmp_path_factory.mktemp("unrun_t") / "canonical",
        config=CanonicalModelConfig.testing(),
    )


def _heads(model, kstpkper):
    nlay = model.gwf.modelgrid.nlay
    return np.vstack([
        model.hds.array(layer=layer, kstpkper=kstpkper, masked=False) for layer in range(nlay)
    ])


# --- the rule, on numbers you can check by hand ----------------------------------

TOP = [10.0, 10.0]
BOTM = [[6.0, 6.0], [0.0, 0.0]]  # layer 0 is 4 thick, layer 1 is 6


def test_the_maximum_is_k_times_the_full_thickness():
    result = transmissivity([[2.0, 3.0], [1.0, 1.0]], TOP, BOTM)
    assert result.transmissivity.tolist() == [[8.0, 12.0], [6.0, 6.0]]
    assert not result.saturated


def test_a_convertible_layer_conducts_through_its_saturated_thickness():
    heads = [[8.0, 12.0], [5.0, -1.0]]  # half full / above its top | 5 of 6 / below
    result = transmissivity(1.0, TOP, BOTM, heads=heads, icelltype=[1, 1])
    assert result.thickness.tolist() == [[2.0, 4.0], [5.0, 0.0]]
    assert np.array_equal(result.transmissivity, result.k * result.thickness)


def test_a_confined_layer_conducts_through_its_full_thickness_whatever_the_head():
    """MF6 treats ICELLTYPE 0 as fully saturated -- head minus bottom would be wrong."""

    heads = [[8.0, 7.0], [3.0, 3.0]]
    result = transmissivity(1.0, TOP, BOTM, heads=heads, icelltype=[1, 0])
    assert result.thickness.tolist() == [[2.0, 1.0], [6.0, 6.0]]


def test_a_dry_cell_conducts_nothing():
    """MF6 writes HDRY (-1e30); that is not a water level 1e30 ft down."""

    result = transmissivity(1.0, TOP, BOTM, heads=[[-1e30, 9.0], [5.0, 5.0]], icelltype=1)
    assert result.transmissivity[0].tolist() == [0.0, 3.0]


def test_inactive_cells_are_gaps_including_vertical_pass_through():
    """MF6 removes idomain -1 from horizontal flow as well as 0."""

    result = transmissivity(1.0, TOP, BOTM, idomain=[[1, 0], [-1, 1]])
    assert np.isnan(result.transmissivity[0, 1]) and np.isnan(result.transmissivity[1, 0])
    assert result.transmissivity[0, 0] == 4.0 and result.transmissivity[1, 1] == 6.0


def test_thickstrt_takes_a_negative_celltype_layer_from_the_starting_head():
    result = transmissivity(
        1.0, TOP, BOTM, heads=[[7.0, 7.0], [5.0, 5.0]], icelltype=[[-1, 1], [0, 0]],
        strt=[[9.0, 9.0], [5.0, 5.0]], thickstrt=True,
    )
    assert result.thickness[0].tolist() == [3.0, 1.0]  # strt-based | head-based


def test_thickstrt_without_a_starting_head_refuses():
    with pytest.raises(ValueError, match="strt"):
        transmissivity(1.0, TOP, BOTM, heads=[[7.0] * 2] * 2, icelltype=-1, thickstrt=True)


# --- the rule, against MODFLOW itself -------------------------------------------


def test_the_saturated_thickness_is_mf6s_own_at_every_saved_time(canonical_run):
    """The headline measurement: `DATA-SAT` is what MF6 conducted through."""

    model = canonical_run
    nlay = model.gwf.modelgrid.nlay
    full = model_transmissivity(model).thickness
    saturation = model.budget.data_sat.get()
    stamps = saturation["kstpkper"].apply(lambda key: tuple(int(v) for v in key))
    assert len(model.hds.kstpkper) > 1
    for key in model.hds.kstpkper:
        key = tuple(int(v) for v in key)
        ours = model_transmissivity(model, heads=_heads(model, key)).thickness / full
        rows = saturation[stamps == key]
        for layer in range(nlay):
            at = rows[rows["layer"] == layer].sort_values("cell")
            assert np.allclose(ours[layer][at["cell"].to_numpy()], at["sat"].to_numpy(), atol=1e-6)


def test_saturated_and_maximum_differ_where_the_water_table_is(canonical_run):
    """Why the choice mattered: the measured median in the water-table layer."""

    key = tuple(canonical_run.hds.kstpkper[-1])
    saturated = model_transmissivity(canonical_run, heads=_heads(canonical_run, key))
    maximum = model_transmissivity(canonical_run)
    ratio = saturated.transmissivity[0] / maximum.transmissivity[0]
    assert np.nanmedian(ratio) < 0.5
    confined = np.asarray(canonical_run.gwf.npf.icelltype.array)[:, 0] == 0
    assert np.allclose(saturated.transmissivity[confined], maximum.transmissivity[confined])


# --- the free map ----------------------------------------------------------------


@pytest.mark.parametrize("spelling", ["bound", "free"])
def test_trans_max_draws_on_an_unrun_model(unrun_model, spelling):
    picture = (
        unrun_model.plot.map(type="trans_max", layer=1)
        if spelling == "bound"
        else plot.map(unrun_model, type="trans_max", layer=1)
    )
    expected = model_transmissivity(unrun_model).transmissivity[1]
    assert np.allclose(picture.linear_zs, expected, equal_nan=True)


def test_saturated_trans_on_an_unrun_model_names_the_maximum(unrun_model):
    with pytest.raises(ValueError, match="trans_max"):
        unrun_model.plot.map(type="trans")


def test_trans_reads_the_heads_of_the_period_asked_for(canonical_run):
    key = next(tuple(k) for k in reversed(canonical_run.hds.kstpkper) if int(k[1]) == 2)
    picture = canonical_run.plot.map(type="trans", per=2, layer=0)
    assert picture.kstpkper == key
    expected = model_transmissivity(canonical_run, heads=_heads(canonical_run, key))
    assert np.allclose(picture.linear_zs, expected.transmissivity[0], equal_nan=True)


def test_the_hover_is_t_in_the_models_units_with_k_and_thickness(canonical_run):
    assert unit_labels(canonical_run) == ("ft", "d")
    picture = canonical_run.plot.map(type="trans", layer=0)
    spec = picture.hover_spec
    assert (spec.primary, spec.title) == ("T", "Transmissivity")
    assert spec.units == {"T": "ft²/d", "K": "ft/d", "thickness": "ft"}
    trace = picture.fig.data[0]
    rendered = " ".join(str(value) for value in trace.customdata[0])
    assert "ft²/d" in rendered and "saturated thickness" in trace.hovertemplate
    assert "Head" not in trace.hovertemplate  # T's hover, not the heads rows

    maximum = canonical_run.plot.map(type="trans_max", layer=0).hover_spec
    assert maximum.title == "Max transmissivity" and maximum.footer == ()


def test_a_bare_grid_or_a_transport_model_refuses_by_name(canonical_run):
    with pytest.raises(ValueError, match="bare grid"):
        plot.map(canonical_run.vor, type="trans_max")
    transport = SimpleNamespace(gwf=SimpleNamespace(), model_type="gwt6")
    with pytest.raises(ValueError, match="NPF"):
        flow_model(transport)


# --- the nouns -------------------------------------------------------------------


def test_the_nouns_draw_the_free_maps(canonical_run):
    npf = canonical_run.packages.npf
    for noun, spelling in ((npf.transmissivity, "trans"), (npf.transmissivity_max, "trans_max")):
        picture = noun.map(per=3, layer=1)
        free = canonical_run.plot.map(type=spelling, per=3, layer=1)
        assert np.allclose(picture.linear_zs, free.linear_zs, equal_nan=True)
        assert picture.hover_spec.title == free.hover_spec.title


def test_the_noun_table_is_k_times_thickness_by_period(canonical_run):
    table = canonical_run.packages.npf.transmissivity.get(per=[1, 4], layer=0)
    assert sorted(table["per"].unique()) == [1, 4]
    assert np.allclose(table["transmissivity"], table["k"] * table["thickness"], equal_nan=True)
    maximum = canonical_run.packages.npf.transmissivity_max.get()
    assert "per" not in maximum.columns


def test_a_noun_mosaic_is_one_panel_per_layer_on_both_backends(canonical_run):
    noun = canonical_run.packages.npf.transmissivity
    figure = noun.mosaic(per=2, logscale=True)
    nlay = canonical_run.gwf.modelgrid.nlay
    maps = {trace.subplot for trace in figure.data if trace.type == "choroplethmap"}
    assert len(maps) == nlay
    static = noun.mosaic(per=2, backend="mpl")
    assert len([axis for axis in static.axes if axis.collections and axis.get_label() != "<colorbar>"]) == nlay


def test_the_saturated_noun_animates_over_the_saved_periods(canonical_run):
    periods = sorted({int(k[1]) for k in canonical_run.hds.kstpkper})
    assert len(canonical_run.packages.npf.transmissivity.animate().frames) == len(periods)


def test_the_saturated_noun_refuses_an_unrun_model_by_name(unrun_model):
    with pytest.raises(ValueError, match="transmissivity_max"):
        unrun_model.packages.npf.transmissivity.get()
    assert unrun_model.packages.npf.transmissivity_max.get().shape[0] > 0


def test_derived_fields_are_not_declared_package_arrays(canonical_run):
    """`summary()`/`ModelPackages.mosaic` enumerate what the package STORES."""

    assert canonical_run.packages.npf.fields["field"].tolist() == ["k", "k22", "k33"]


# --- review fixes (ledger 311) ---------------------------------------------------


def test_a_saturated_map_redraws_when_its_time_moves(canonical_run):
    """The cache was unkeyed: every `.ani` frame and a reassigned `per` drew the
    FIRST time's transmissivity while the footer advanced."""

    picture = canonical_run.plot.map(type="trans", layer=0)
    first = np.asarray(picture.linear_zs, dtype=float)
    picture.per = 5
    later = np.asarray(picture.linear_zs, dtype=float)
    expected = model_transmissivity(canonical_run, heads=_heads(canonical_run, picture.kstpkper))
    assert np.allclose(later, expected.transmissivity[0], equal_nan=True)
    assert not np.allclose(later, first, equal_nan=True)

    frames = canonical_run.plot.map(type="trans", layer=0).ani.frames
    sums = {round(float(np.nansum(np.asarray(frame.data[0].z, dtype=float))), 3) for frame in frames}
    assert len(sums) > 1


def test_the_saturated_table_keys_on_the_period(canonical_run):
    noun = canonical_run.packages.npf.transmissivity
    series = noun.long(layer=0)
    assert series.index.names == ["per", "layer", "cell"] and series.index.is_unique
    wide = noun.wide(per=[1, 4])
    assert sorted(wide["per"].unique()) == [1, 4]
    assert noun.map(per=None, layer=0) is not None  # the default time, like every noun
