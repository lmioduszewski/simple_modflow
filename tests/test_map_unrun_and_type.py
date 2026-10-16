"""A map of a model that has not been run, and a map `type` that does not exist
(ledger 305).

Both used to draw NODE NUMBERS -- 0..ncpl-1 under a colorbar -- with no error:

* the dispatch bound a model only when it had RESULTS, so an unrun model was
  drawn as a bare grid: `type='ks'` could not read its K, and a region name in
  `select=` fell through to being opened as a FILE;
* `type` was never checked, so a typo (`type='transmissivity'`) and the documented
  `'rch'` (which had no branch) all fell through to the cell-index branch.

The unrun model is built privately: the session `canonical_model` is the same
object `canonical_run` runs, so whether it has results depends on test order.
"""

from __future__ import annotations

import re

import numpy as np
import pytest

from myflopy import plot
from myflopy.modflow.utils.datatypes.choros import MAP_TYPES, check_map_type


@pytest.fixture(scope="module")
def unrun_model(tmp_path_factory):
    """A canonical model that is built and never run (0.6 s)."""

    from myflopy.modflow.mf6.canonical_example import (
        CanonicalModelConfig,
        build_canonical_model,
    )

    return build_canonical_model(
        tmp_path_factory.mktemp("unrun") / "canonical",
        config=CanonicalModelConfig.testing(),
    )


def _is_node_numbers(picture) -> bool:
    values = np.asarray(picture.linear_zs, dtype=float)
    return np.array_equal(values, np.arange(values.size))


# --- an unrun model is still a model ------------------------------------------


@pytest.mark.parametrize("spelling", ["bound", "free"])
def test_an_unrun_model_maps_its_conductivity(unrun_model, spelling):
    """`type='ks'` reads K, an INPUT, so it must not need a run."""

    if spelling == "bound":
        picture = unrun_model.plot.map(type="ks", layer=1)
    else:
        picture = plot.map(unrun_model, type="ks", layer=1)

    expected = np.asarray(unrun_model.gwf.npf.k.array)[1]
    assert np.allclose(np.asarray(picture.linear_zs, dtype=float), expected)
    assert not _is_node_numbers(picture)
    # It is NOT pretending to have results: the results binding stays empty.
    assert picture.model is None


def test_an_unrun_models_region_names_resolve(unrun_model):
    """`select="all_streams"` raised `DataSourceError: No such file` before a run."""

    ncpl = unrun_model.vor.ncpl
    expected = sorted(int(c) for c in unrun_model.get_region_cells("all_streams"))

    picture = unrun_model.plot.map(values=np.arange(ncpl, dtype=float), select="all_streams")
    assert picture._select_cells == expected

    # ...and through a noun, which calls the free verb underneath.
    noun = unrun_model.npf.k.map(select="all_streams")
    assert noun._select_cells == expected


def test_an_unrun_model_refuses_a_results_map_by_name(unrun_model):
    """Heads of an unrun model were node numbers under a head colorbar."""

    with pytest.raises(ValueError, match="has not been run"):
        unrun_model.plot.map()
    with pytest.raises(ValueError, match="has not been run"):
        unrun_model.plot.map(type="rch")


def test_values_still_draw_on_an_unrun_model(unrun_model):
    """`values=` names its own numbers, so the default `type` asks nothing of a run."""

    values = np.linspace(0.0, 1.0, unrun_model.vor.ncpl)
    picture = unrun_model.plot.map(values=values)
    assert np.allclose(picture.linear_zs, values)


def test_a_bare_grid_still_draws_its_node_numbers(canonical_model):
    """The grid map's subject IS the node numbering; that must not move."""

    assert _is_node_numbers(canonical_model.vor.plot.map())


# --- type is checked -------------------------------------------------------------


@pytest.mark.parametrize("value", ["transmissivity", "bogus", "Ks_", ""])
def test_an_unknown_type_raises_naming_the_set(unrun_model, value):
    with pytest.raises(ValueError, match="not a map type") as caught:
        unrun_model.plot.map(type=value)
    assert "'ks'" in str(caught.value)


def test_input_rch_is_retired_by_name():
    """It read stress period 37 whatever was asked, and KeyError'd on smaller models."""

    with pytest.raises(ValueError, match=r"model\.rch\.inputs"):
        check_map_type("input_rch")


def test_type_is_case_insensitive(canonical_run):
    """`'HDS'` drew node numbers under a Head hover: one lookup lowered, one did not."""

    upper = canonical_run.plot.map(type="HDS", layer=0)
    lower = canonical_run.plot.map(type="hds", layer=0)
    assert not _is_node_numbers(upper)
    assert np.allclose(
        np.asarray(upper.linear_zs, dtype=float), np.asarray(lower.linear_zs, dtype=float)
    )


def _recharge_truth(model, kstpkper):
    """Recharge RATE per cell straight from FloPy's reader of the budget file."""

    import flopy

    cbc = next(model.workspace.glob("*.cbc"))
    record = flopy.utils.CellBudgetFile(str(cbc)).get_data(text="RCH", kstpkper=kstpkper)[0]
    ncpl = model.vor.ncpl
    cells = (np.asarray(record["node"], dtype=int) - 1) % ncpl  # FloPy's node is 1-based
    flux = np.bincount(cells, weights=np.asarray(record["q"], dtype=float), minlength=ncpl)
    return flux / np.asarray(model.vor.area_list, dtype=float)


def test_the_documented_rch_type_draws_recharge_on_the_right_cells(canonical_run):
    """`'rch'` was documented and drew node numbers; routed onto the old reader it
    drew every cell's rate on the PREVIOUS cell (map[i] == truth[i+1]) -- and the
    test that shipped with it compared against that same reader, so it passed
    (ledger 311). The truth here is FloPy's own budget-file reader."""

    for per in (0, 2):
        picture = canonical_run.plot.map(type="rch", per=per, layer=0)
        truth = _recharge_truth(canonical_run, picture.kstpkper)
        assert truth.max() > 0.0
        assert np.allclose(np.asarray(picture.linear_zs, dtype=float), truth)


def test_a_bare_grid_refuses_a_field_type_it_cannot_read(canonical_model):
    """A bare grid has no model: `type='ks'` drew node numbers under a K colorbar."""

    with pytest.raises(TypeError, match="bare grid"):
        canonical_model.vor.plot.map(type="ks")
    with pytest.raises(ValueError, match="bare grid"):
        plot.map(canonical_model.vor, type="rch")
    assert _is_node_numbers(plot.map(canonical_model.vor))  # the grid's own map


def test_the_documented_type_set_is_the_accepted_one():
    """Both directions: nothing documented is refused, nothing accepted is hidden
    except the one legacy alias `'output_rch'`."""

    match = re.search(r"^\s*type : \{([^}]*)\}", plot.map.__doc__, flags=re.M)
    assert match, "plot.map no longer documents its type set"
    documented = {item.strip().strip("'\"") for item in match.group(1).split(",")}
    assert documented == set(MAP_TYPES) - {"output_rch"}
    for value in documented:
        assert check_map_type(value) == value


# --- layer is checked at the call ------------------------------------------------


@pytest.mark.parametrize(
    ("layer", "error", "match"),
    [
        (4, ValueError, r"outside this model's 4 layers, 0\.\.3"),
        (None, TypeError, "zero-based integer"),
        ("all", TypeError, "zero-based integer"),
        (True, TypeError, "zero-based integer"),
        (-1, ValueError, "show_mounding"),
    ],
)
def test_a_layer_that_cannot_be_drawn_raises_at_the_call(unrun_model, layer, error, match):
    """All of these constructed a picture that died at `.fig` -- or, for `-1`,
    silently drew the BOTTOM layer through numpy's negative indexing."""

    with pytest.raises(error, match=match):
        unrun_model.plot.map(type="ks", layer=layer)


def test_a_numpy_layer_is_accepted(unrun_model):
    picture = unrun_model.plot.map(type="ks", layer=np.int64(2))
    assert type(picture.layer) is int and picture.layer == 2


def test_layer_minus_one_keeps_its_mounding_meaning(canonical_run):
    picture = canonical_run.plot.map(layer=-1, show_mounding=True)
    assert picture._show_mounding_above_ground
    assert picture.layer == 0


# --- values= is ONE layer's field, checked at the call (ledger 601) -------
#
# `ies.calibrated().get("npf", "k")` is `(nlay, ncpl)`, and `map(values=)` took it:
# both verbs did `custom_zs = list(values)`, which made it nlay ROWS, and Plotly pairs
# z with locations by position, so it drew z of length 4 against 441 cells, a
# flattened field's LAYER 0 under `layer=1`, a dict's keys and a DataFrame's column
# names -- while Matplotlib raised pandas' "Length of values" error naming neither.
# The check lives in `Choro.__init__`, the funnel every map reaches.


def _k(model):
    return np.asarray(model.gwf.npf.k.array, dtype=float)


# The four spellings that take values=: the free verb on a model and on a grid, the
# bound model namespace, and the grid namespace (its own path -- it calls
# `_choropleth_factory` directly, and `model.hfb.inputs.map` delegates to it).
VALUES_SCOPES = {
    "free_model": lambda m, v, **kw: plot.map(m, v, layer=1, **kw),
    "bound_model": lambda m, v, **kw: m.plot.map(values=v, layer=1, **kw),
    "free_grid": lambda m, v, **kw: plot.map(m.vor, v, **kw),
    "grid": lambda m, v, **kw: m.vor.plot.map(values=v, **kw),
}
MODEL_SCOPES = {"free_model", "bound_model"}


def _render(picture):
    """Force the picture to draw (a Choro builds its trace lazily)."""

    return picture.fig if hasattr(picture, "fig") else picture


def _mpl_drawn(figure):
    """The per-cell values a static map coloured its cells by."""

    return np.ma.filled(np.ma.asarray(figure.axes[0].collections[0].get_array(), dtype=float), np.nan)


@pytest.mark.parametrize("backend", ["plotly", "mpl"])
@pytest.mark.parametrize("scope", sorted(VALUES_SCOPES))
@pytest.mark.parametrize("flatten", [False, True], ids=["stacked", "flattened"])
def test_a_field_for_every_layer_raises_naming_one_layer(unrun_model, scope, backend, flatten):
    """Plotly drew z of length 4 against 441 cells, and a flattened field painted
    LAYER 0 under `layer=1`; the refusal names the one-layer spelling."""

    k = _k(unrun_model)
    values = k.reshape(-1) if flatten else k
    row = r"array\.reshape\(4, 441\)" if flatten else "array"
    expected = rf"values={row}\[1\], layer=1" if scope in MODEL_SCOPES else rf"values={row}\[i\]"
    with pytest.raises(ValueError, match=r"ONE layer.* 441 numbers in cell order") as caught:
        _render(VALUES_SCOPES[scope](unrun_model, values, backend=backend))
    assert caught.match(expected)


@pytest.mark.parametrize("backend", ["plotly", "mpl"])
@pytest.mark.parametrize("scope", sorted(VALUES_SCOPES))
@pytest.mark.parametrize("size", [100, 442])
def test_a_wrong_length_raises_naming_the_cell_count(unrun_model, scope, backend, size):
    with pytest.raises(ValueError, match=rf"ONE layer.*shape \(441,\).*got shape \({size},\).*441 cells"):
        _render(VALUES_SCOPES[scope](unrun_model, np.ones(size), backend=backend))


@pytest.mark.parametrize("scope", sorted(VALUES_SCOPES))
def test_a_per_period_dict_raises_naming_one_period(unrun_model, scope):
    """`get("rch", "recharge")` is `{period: array}`; plotly drew its KEYS."""

    k = _k(unrun_model)
    with pytest.raises(TypeError, match=r"not a dict.*values=data\[period\]"):
        _render(VALUES_SCOPES[scope](unrun_model, {0: k[0], 1: k[1]}))


def test_a_table_raises(unrun_model):
    """Plotly drew a one-column DataFrame as its column NAME."""

    import pandas as pd

    table = pd.DataFrame({"k": _k(unrun_model)[1]})
    with pytest.raises(TypeError, match=r"not a table \(441 rows x 1 columns\)"):
        _render(unrun_model.plot.map(values=table, layer=1))


def test_a_single_number_raises(unrun_model):
    with pytest.raises(TypeError, match="not the single value 3.0"):
        _render(unrun_model.plot.map(values=3.0, layer=1))


@pytest.mark.parametrize("spelling", ["list", "ndarray", "series", "row"])
@pytest.mark.parametrize("backend", ["plotly", "mpl"])
def test_one_layer_draws_that_layer(unrun_model, spelling, backend):
    """Every one-layer spelling draws exactly that layer -- a single `(1, ncpl)` row
    too, which is what a ONE-layer model's `calibrated().get("npf", "k")` returns
    (plotly drew it as z of length 1; mpl raised)."""

    import pandas as pd

    row = _k(unrun_model)[1]
    values = {
        "list": list(row),
        "ndarray": row,
        "series": pd.Series(row, index=np.arange(row.size)[::-1]),  # read by POSITION
        "row": row[None, :],
    }[spelling]
    picture = unrun_model.plot.map(values=values, layer=1, backend=backend)
    drawn = (
        np.asarray(picture.fig.data[0].z, dtype=float) if backend == "plotly" else _mpl_drawn(picture)
    )
    assert drawn.shape == row.shape  # allclose alone broadcasts a (1, 441) row
    assert np.array_equal(drawn, row)


def test_a_masked_cell_draws_as_a_gap(unrun_model):
    """A mask means "no value": plotly painted masked cells 0.0 (measured, 63 of 63)."""

    row = _k(unrun_model)[1]
    mask = np.arange(row.size) % 7 == 0
    picture = unrun_model.plot.map(values=np.ma.masked_array(row, mask=mask), layer=1)
    drawn = np.asarray(picture.fig.data[0].z, dtype=float)  # the TRACE: 0.0 before
    assert np.isnan(drawn[mask]).all()
    assert np.array_equal(drawn[~mask], row[~mask])


def test_the_choro_funnel_refuses_a_mis_shaped_custom_zs(unrun_model):
    """Every map reaches `Choro.__init__`, the nouns' own `custom_zs` included, so the
    check lives there; the setter goes through the same rule, not an `assert`."""

    from myflopy.modflow.mf6.grid.plotting import _choropleth_factory

    with pytest.raises(ValueError, match="441"):
        _choropleth_factory(unrun_model.vor, custom_zs=[1.0] * 5)
    picture = _choropleth_factory(unrun_model.vor, custom_zs=np.ones(441))
    with pytest.raises(ValueError, match="441"):
        picture.custom_zs = [1.0] * 5
    picture.custom_zs = np.full((1, 441), 2.0)
    assert picture.custom_zs == [2.0] * 441


def test_a_mosaic_panel_refuses_too(unrun_model):
    """`model.plot.mosaic([...])` composed a garbled panel without a word."""

    with pytest.raises(ValueError, match="ONE layer"):
        unrun_model.plot.mosaic([unrun_model.plot.map(values=_k(unrun_model), layer=1)])


@pytest.mark.parametrize("spelling", ["generator", "dict_values", "map_object", "range", "tuple"])
def test_every_iterable_list_read_still_draws(unrun_model, spelling):
    """`list(values)` read any iterable; `np.asarray` makes a generator a 0-d object
    array, so the check materializes one first. A guard: these drew before the fix."""

    row = _k(unrun_model)[1]
    values = {
        "generator": lambda: (x for x in row),
        "dict_values": lambda: dict(enumerate(row)).values(),
        "map_object": lambda: map(float, row),
        "range": lambda: range(row.size),
        "tuple": lambda: tuple(row),
    }[spelling]()
    expected = np.arange(row.size, dtype=float) if spelling == "range" else row
    drawn = np.asarray(unrun_model.plot.map(values=values, layer=1).fig.data[0].z, dtype=float)
    assert np.array_equal(drawn, expected)
