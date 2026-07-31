"""Defects that a broad `except Exception` was hiding (plan 7.3).

Narrowing an exception handler is only half of 7.3. The other half is what the
narrowing turned up: a handler written to catch "anything" tends to be written
once and never re-read, and the code inside it drifts. Each test here pins one
concrete wrong answer the old handlers were returning, so the fix cannot quietly
regress.

These are cheap unit tests on purpose. Every defect below was reachable from a
supported public path, but reproducing it end-to-end would mean building a
broken model -- which is exactly the state where nobody is watching closely.
"""

from __future__ import annotations

import logging
from types import SimpleNamespace

import numpy as np
import pytest

from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.modflow.utils.datatypes.readers import read_gpkg


def _grid(**kwargs):
    """A stand-in modelgrid; attributes absent from kwargs raise AttributeError."""

    return SimpleNamespace(**kwargs)


# --- nlay was reporting a cell count ------------------------------------------
def test_nlay_falls_back_to_the_layer_count_not_the_number_of_botm_values():
    """`_model_nlay`'s fallback read `np.asarray(botm).reshape(-1).size` -- the
    total number of bottom-elevation ENTRIES. On a 3-layer, 4-cell grid that is
    12, reported as the layer count, and the broad handler around it meant the
    absurd answer was returned rather than raised."""

    botm = np.array([[10.0, 10.0, 10.0, 10.0],
                     [5.0, 5.0, 5.0, 5.0],
                     [0.0, 0.0, 0.0, 0.0]])
    model = SimpleNamespace(modelgrid=_grid(botm=botm))   # no `nlay` attribute

    assert SimulationBase._model_nlay(model) == 3


def test_nlay_prefers_the_grid_over_the_fallback():
    model = SimpleNamespace(modelgrid=_grid(nlay=7, botm=np.zeros((3, 4))))
    assert SimulationBase._model_nlay(model) == 7


def test_a_one_dimensional_botm_is_a_single_layer():
    """DISU stores one bottom per node, so its first axis is the NODE count."""

    model = SimpleNamespace(modelgrid=_grid(botm=np.zeros(25)))
    assert SimulationBase._model_nlay(model) == 1


@pytest.mark.parametrize(
    "modelgrid",
    [_grid(), _grid(botm=None), _grid(botm=np.array([]))],
    ids=["no-botm", "botm-is-None", "botm-is-empty"],
)
def test_nlay_is_unknown_rather_than_wrong_when_there_is_no_geometry(modelgrid):
    assert SimulationBase._model_nlay(SimpleNamespace(modelgrid=modelgrid)) is None


def test_the_nlay_fallback_says_so_at_debug(caplog):
    """The swallow convention: a fallback that fires leaves a trace."""

    model = SimpleNamespace(modelgrid=_grid(botm=np.zeros((2, 4))))
    with caplog.at_level(logging.DEBUG, logger="myflopy"):
        SimulationBase._model_nlay(model)
    assert "deriving it from botm" in caplog.text


# --- a partial GeoPackage read reported success -------------------------------
def test_an_unsupported_geometry_type_is_raised_not_swallowed(tmp_path):
    """`read_gpkg` raises `TypeError: Unexpected geometry type` for a geometry it
    cannot place -- and the bare `except:` four lines below CAUGHT ITS OWN
    RAISE, ending the read and returning whatever it had collected so far. A
    user got a GeoDataFrame that was silently missing features.
    """

    import geopandas as gpd
    import shapely as shp

    path = tmp_path / "mixed.gpkg"
    gpd.GeoDataFrame(
        {"name": ["ok", "unsupported"]},
        geometry=[
            shp.Point(0, 0),
            shp.GeometryCollection([shp.Point(1, 1), shp.Point(2, 2)]),
        ],
        crs="EPSG:2927",
    ).to_file(path, driver="GPKG")

    with pytest.raises(TypeError, match="Unexpected geometry type"):
        read_gpkg(path)


def test_a_geopackage_with_no_usable_geometry_says_which_file(tmp_path):
    """The old code only produced its "Could not read gpkg file" message when an
    exception happened to fire; an empty layer fell through to a downstream
    error about a missing geometry column instead."""

    import geopandas as gpd

    path = tmp_path / "empty.gpkg"
    gpd.GeoDataFrame({"name": []}, geometry=[], crs="EPSG:2927").to_file(
        path, driver="GPKG"
    )

    with pytest.raises(ValueError, match="no usable geometries"):
        read_gpkg(path)
