"""Calibration plots draw what the table holds and refuse what they cannot draw (ledger 477).

Each test is one of ledger 436's calibration defects: ``type="heads"`` raised
``KeyError`` unless every location had rows in every layer; a table with neither
``per`` nor ``time`` raised a bare ``KeyError('time')``; a lone ``target_column`` was
ignored; a typo'd ``type`` drew the calibration plot; ``title``/``yaxis_title`` were
ignored for heads and residuals; the constructor's calibration plot labelled its axes
the wrong way round for the points it drew; ``from_timeseries`` drew an empty baseline
and one line zig-zagging between a well's two layers; a bound ``ConcTargets``' pictures
were titled for heads; ``to_flopy_obs`` wrote any ``kind`` and gave a ``ConcTargets``
``HEAD``.
"""

from __future__ import annotations

from types import SimpleNamespace

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import Point

import myflopy as mf
from myflopy.modflow.calcs.calibration import CalibrationPlot
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.observations.heads import BoundHeadTargetPlots, HeadTargets


def _compare(rows):
    """A heads compare table from (name, layer, per, target, simulated) rows."""

    return pd.DataFrame(rows, columns=["name", "layer", "per", "head_target", "sim_head"]).assign(
        time=lambda frame: frame["per"],
        abs_residual=lambda frame: (frame["sim_head"] - frame["head_target"]).abs(),
    )


_TWO_WELLS = _compare([
    ("W1", 0, 0, 10.0, 10.5), ("W1", 0, 1, 11.0, 11.5),
    ("W2", 1, 0, 20.0, 19.0), ("W2", 1, 1, 21.0, 20.0),
])
_TWO_LAYERS = _compare([
    ("W1", 0, 0, 10.0, 10.5), ("W1", 1, 0, 30.0, 29.0),
    ("W1", 0, 1, 11.0, 11.5), ("W1", 1, 1, 31.0, 30.0),
])


def test_heads_draws_a_well_in_each_layer_it_has_rows_in():
    figure = CalibrationPlot.from_compare(_TWO_WELLS, type="heads")
    assert [trace.name for trace in figure.data] == [
        "W1 observed", "W1 simulated", "W2 observed", "W2 simulated",
    ]
    layered = CalibrationPlot.from_compare(_TWO_LAYERS, type="heads")
    assert [trace.name for trace in layered.data] == [
        "W1 layer 0 observed", "W1 layer 0 simulated", "W1 layer 1 observed", "W1 layer 1 simulated",
    ]


def test_heads_without_a_period_column_is_refused_by_name():
    with pytest.raises(ValueError, match="neither a 'per' nor a 'time' column"):
        CalibrationPlot.from_compare(_TWO_WELLS.drop(columns=["per", "time"]), type="heads")


@pytest.mark.parametrize("given", [{"target_column": "head_target"}, {"simulated_column": "sim_head"}])
def test_one_column_of_the_pair_is_refused(given):
    with pytest.raises(ValueError, match="without"):
        CalibrationPlot.from_compare(_TWO_WELLS, **given)


def test_a_misspelt_type_is_refused():
    with pytest.raises(ValueError, match="is not a plot it draws"):
        CalibrationPlot.from_compare(_TWO_WELLS, type="calibraton")
    with pytest.raises(ValueError, match="draws 'calibration' or 'heads'"):
        CalibrationPlot(observed=[1.0, 2.0], simulated=[1.5, 2.5], type="head")


def test_titles_reach_the_heads_and_residual_plots():
    heads = CalibrationPlot.from_compare(_TWO_WELLS, type="heads", title="T", yaxis_title="Y")
    assert (heads.layout.title.text, heads.layout.yaxis.title.text) == ("T", "Y")
    residuals = CalibrationPlot.from_compare(
        _TWO_WELLS[["time", "abs_residual"]], type="residuals_by_period", title="R", yaxis_title="E")
    assert (residuals.layout.title.text, residuals.layout.yaxis.title.text) == ("R", "E")


def test_the_calibration_cross_plot_puts_observed_on_its_observed_axis():
    observed = pd.Series([1.0, 2.0, 3.0])
    simulated = pd.Series([1.5, 2.5, 3.5])
    figure = CalibrationPlot(observed=observed, simulated=simulated)
    assert figure.layout.xaxis.title.text == "Observed"
    assert figure.layout.yaxis.title.text == "Simulated"
    assert list(figure.data[0].x) == observed.tolist()


def test_a_timeseries_draws_one_line_per_layer_and_refuses_an_absent_baseline():
    figure = CalibrationPlot.from_timeseries(_TWO_LAYERS)
    assert [trace.name for trace in figure.data] == [
        "Target (layer 0)", "Model (layer 0)", "Target (layer 1)", "Model (layer 1)",
    ]
    assert list(figure.data[1].y) == [10.5, 11.5]
    with pytest.raises(ValueError, match="baseline_compare has no rows for observation name 'W1'"):
        CalibrationPlot.from_timeseries(
            _TWO_WELLS, name="W1", baseline_compare=_TWO_WELLS.query("name == 'W2'"))


# --- the bound targets' pictures and their OBS records -----------------------------


def _grid(n: int = 4, length: float = 100.0) -> VoronoiGridPlus:
    d = length / n
    verts = np.array([[i * d, j * d] for j in range(n + 1) for i in range(n + 1)], dtype=float)
    iverts = [
        [j * (n + 1) + i, (j + 1) * (n + 1) + i, (j + 1) * (n + 1) + i + 1, j * (n + 1) + i + 1]
        for j in range(n)
        for i in range(n)
    ]
    xcyc = np.array([[(i + 0.5) * d, (j + 0.5) * d] for j in range(n) for i in range(n)])
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)


def _model(grid, nlay: int = 1):
    ncpl = len(grid.gdf_vorPolys)
    return SimpleNamespace(
        vor=grid,
        gwf=SimpleNamespace(
            disv=SimpleNamespace(idomain=SimpleNamespace(array=np.ones((nlay, ncpl), dtype=int)))
        ),
    )


def _locations():
    return gpd.GeoDataFrame({"name": ["W1"], "layer": [0]}, geometry=[Point(37.5, 37.5)],
                            crs="EPSG:2927")


def test_obs_records_are_the_familys_own_type():
    model = _model(_grid())
    heads = HeadTargets(locations=_locations(), values=pd.DataFrame({"per": [0], "W1": [1.0]}),
                        time_column="per")
    conc = mf.ConcTargets(locations=_locations(), values=pd.DataFrame({"per": [0], "W1": [1.0]}),
                          time_column="per")
    assert heads.to_flopy_obs(model)["head_targets.csv"][0][1] == "HEAD"
    assert conc.to_flopy_obs(model)["head_targets.csv"][0][1] == "CONCENTRATION"
    assert heads.to_flopy_obs(model, kind="head")["head_targets.csv"][0][1] == "HEAD"
    with pytest.raises(ValueError, match="writes HEAD observations"):
        heads.to_flopy_obs(model, kind="bogus")
    with pytest.raises(ValueError, match="writes CONCENTRATION observations"):
        conc.to_flopy_obs(model, kind="HEAD")


def test_a_bound_conc_targets_pictures_are_titled_for_concentration():
    conc = mf.ConcTargets(locations=_locations(), values=pd.DataFrame({"per": [0], "W1": [1.0]}),
                          time_column="per")
    table = _compare([("W1", 0, 0, 1.0, 1.2)])
    stub = SimpleNamespace(compare=lambda model: table.rename(columns={}),
                           _FIELD_LABEL=conc._FIELD_LABEL)
    plots = BoundHeadTargetPlots(SimpleNamespace(model=None, targets=stub))
    assert plots.obs_vs_sim().layout.title.text == "Observed vs simulated concentrations"
    assert plots.timeseries().layout.yaxis.title.text == "Concentration"
    with pytest.raises(TypeError, match="ax"):
        plots.obs_vs_sim(ax=object())


def test_calibration_plots_read_the_backend_through_the_one_resolver():
    """``"go"``/``"sns"`` were this module's own aliases (ledger 480)."""

    table = _TWO_WELLS
    assert CalibrationPlot.from_obs_vs_sim(table, backend=" Interactive ").data
    with pytest.raises(ValueError, match="backend must be"):
        CalibrationPlot.from_obs_vs_sim(table, backend="go")
