"""A map of a model that has not been run, and a map `type` that does not exist
(ledger 305).

Both used to draw NODE NUMBERS -- 0..ncpl-1 under a colorbar -- with no error:

* the dispatch bound a model only when it had RESULTS, so an unrun model was
  drawn as a bare grid: `type='ks'` could not read its K, and a region name in
  `select=` fell through to being opened as a FILE;
* `type` was never checked, so `type='trans'`, a typo, and the documented
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
    noun = unrun_model.packages.npf.k.map(select="all_streams")
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


@pytest.mark.parametrize("value", ["trans", "bogus", "Ks_", ""])
def test_an_unknown_type_raises_naming_the_set(unrun_model, value):
    with pytest.raises(ValueError, match="not a map type") as caught:
        unrun_model.plot.map(type=value)
    assert "'ks'" in str(caught.value)


def test_input_rch_is_retired_by_name():
    """It read stress period 37 whatever was asked, and KeyError'd on smaller models."""

    with pytest.raises(ValueError, match=r"packages\.rch\.inputs"):
        check_map_type("input_rch")


def test_type_is_case_insensitive(canonical_run):
    """`'HDS'` drew node numbers under a Head hover: one lookup lowered, one did not."""

    upper = canonical_run.plot.map(type="HDS", layer=0)
    lower = canonical_run.plot.map(type="hds", layer=0)
    assert not _is_node_numbers(upper)
    assert np.allclose(
        np.asarray(upper.linear_zs, dtype=float), np.asarray(lower.linear_zs, dtype=float)
    )


def test_the_documented_rch_type_draws_recharge(canonical_run):
    """`'rch'` was documented, carried a Recharge hover row, and drew node numbers."""

    picture = canonical_run.plot.map(type="rch", layer=0)
    assert not _is_node_numbers(picture)
    values = np.asarray(picture.linear_zs, dtype=float)
    assert values.max() > 0.0
    # It IS the recharge rate: the same numbers the older spelling draws.
    alias = np.asarray(canonical_run.plot.map(type="output_rch", layer=0).linear_zs, dtype=float)
    assert np.array_equal(values, alias)


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
        (4, ValueError, r"outside this model's layers 0\.\.3"),
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
