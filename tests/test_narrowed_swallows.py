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


# --- a misspelled contour method drew nothing ---------------------------------
def _square_grid():
    """A stand-in grid exposing only the centroids the contour code reads."""

    return SimpleNamespace(
        centroids_x=np.array([0.0, 1.0, 0.0, 1.0, 0.5]),
        centroids_y=np.array([0.0, 0.0, 1.0, 1.0, 0.5]),
    )


def test_an_unknown_contour_method_is_named_not_silently_empty():
    """`contour_line_segments` raises `ValueError: contour method must be...`
    for an unrecognized `method=` -- and the `except Exception` that wrapped the
    whole body caught its own raise and returned `[]`. A typo produced a map
    with no contours and no explanation.
    """

    from myflopy.modflow.mf6.contour_plotting import contour_line_segments

    with pytest.raises(ValueError, match="must be 'linear' or 'cubic'"):
        contour_line_segments(_square_grid(), [0.0, 1.0, 1.0, 2.0, 1.0], method="lienar")


def test_the_method_is_checked_before_the_data_is_given_up_on():
    """The validation sits above the `return []` early exits, so a typo is
    reported even for values too flat to contour -- otherwise the message
    depends on the data, which is the least helpful possible behaviour."""

    from myflopy.modflow.mf6.contour_plotting import contour_line_segments

    with pytest.raises(ValueError, match="must be 'linear' or 'cubic'"):
        contour_line_segments(_square_grid(), [1.0] * 5, method="lienar")


def test_a_valid_contour_method_still_draws():
    """The guard must not have cost the working path."""

    from myflopy.modflow.mf6.contour_plotting import contour_line_segments

    segments = contour_line_segments(
        _square_grid(), [0.0, 1.0, 1.0, 2.0, 1.0], levels=[0.5, 1.5], method="linear"
    )
    assert segments and {"level", "x", "y"} <= set(segments[0])


# --- resample handed the caller None ------------------------------------------
def test_a_bad_resample_rule_raises_instead_of_returning_none():
    """`_resample_timeseries_df` did `return print(...)` on failure: it printed
    a message and returned None, so the real symptom arrived later as an
    AttributeError on None, far from the bad argument."""

    import pandas as pd

    from myflopy.modflow.mf6.mf2Dplots import WaterLevelPlot

    frame = pd.DataFrame(
        {"when": pd.date_range("2024-01-01", periods=3, freq="D"), "wl": [1.0, 2.0, 3.0]}
    )
    with pytest.raises(ValueError, match="valid pandas resample rule"):
        WaterLevelPlot._resample_timeseries_df(frame, time_step="nonsense", time_col_name="when")


# --- a masked value is "not int-able", not a crash -----------------------------
def test_nested_int_conversion_leaves_masked_values_alone():
    """`convert_nested_to_int`'s contract is "return obj unchanged if it is not
    int-able", and `np.ma.MaskError` subclasses Exception DIRECTLY -- so a
    narrowing to the obvious (TypeError, ValueError) would have let it escape.
    Masked values are the normal shape of flopy head output at nodata cells."""

    from myflopy.modflow.utils.datatypes.datalists import convert_nested_to_int

    result = convert_nested_to_int([1, "2", None, np.ma.masked])
    assert result[:3] == [1, 2, None]
    assert result[3] is np.ma.masked
