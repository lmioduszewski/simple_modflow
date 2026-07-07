"""Tests for the unified grammar's ``plot`` (series) and ``xs`` (cross-section)
panel verbs -- Phase 1 of the panel-verbs/composer redesign.

Fast tests exercise the generic series engine and the section-line
normalization on fakes; the slow test renders the real verbs on the canonical
model across the three surfaces (model / group / diff).
"""

from __future__ import annotations

import matplotlib
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pytest
from matplotlib.figure import Figure
from shapely.geometry import LineString

matplotlib.use("Agg")

from myflopy.modflow.mf6.package_plotting import FieldMappable, SpatialView
from myflopy.modflow.utils.datatypes.xsections import _normalize_section_line


# --- section-line normalization ------------------------------------------------
def test_normalize_section_line_accepts_all_specs():
    line = LineString([(0.0, 0.0), (10.0, 10.0)])
    assert _normalize_section_line(line) is line

    from_pairs = _normalize_section_line([(0, 0), (5, 5), (10, 0)])
    assert isinstance(from_pairs, LineString)
    assert len(from_pairs.coords) == 3

    from_dict = _normalize_section_line({"line": [(0, 0), (10, 10)]})
    assert isinstance(from_dict, LineString)


def test_normalize_section_line_rejects_bad_specs():
    with pytest.raises(TypeError):
        _normalize_section_line(42)
    with pytest.raises(ValueError):
        _normalize_section_line([(0, 0)])  # one point is not a line


# --- generic series engine (SpatialView.plot) ----------------------------------
class _FakeSeriesHost(SpatialView):
    value_name = "q"

    def __init__(self, frame, models=None):
        self._frame = frame
        self._models = models

    def get(self):
        return self._frame

    def _spatial_models(self):
        return self._models


def _series_frame():
    return pd.DataFrame(
        {
            "model": ["a"] * 4 + ["b"] * 4,
            "per": [0, 0, 1, 1] * 2,
            "layer": [0] * 8,
            "cell": [1, 2, 1, 2] * 2,
            "lake": [0, 1, 0, 1] * 2,
            "q": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
        }
    )


def test_plot_draws_one_line_per_model_and_entity():
    host = _FakeSeriesHost(_series_frame(), models=["a", "b"])
    fig = host.plot()
    assert isinstance(fig, go.Figure)
    assert len(fig.data) == 4  # 2 models x 2 lakes
    names = {trace.name for trace in fig.data}
    assert names == {"a / Lake 0", "a / Lake 1", "b / Lake 0", "b / Lake 1"}
    # values aggregate (sum) over cells inside each line: a/Lake0 = q at cell 1
    lake0_a = next(t for t in fig.data if t.name == "a / Lake 0")
    assert list(lake0_a.x) == [0, 1] and list(lake0_a.y) == [1.0, 3.0]


def test_plot_cells_selection_adds_per_cell_lines_and_filters():
    host = _FakeSeriesHost(_series_frame(), models=["a", "b"])
    fig = host.plot(cells=[1])
    # cell 1 belongs to lake 0 only -> one line per model
    assert {trace.name for trace in fig.data} == {"a / Lake 0 / C1", "b / Lake 0 / C1"}


def test_plot_model_filter_and_backends():
    host = _FakeSeriesHost(_series_frame(), models=["a", "b"])
    fig = host.plot("b")
    assert {trace.name for trace in fig.data} == {"b / Lake 0", "b / Lake 1"}

    mpl_fig = host.plot(backend="mpl")
    assert isinstance(mpl_fig, Figure)

    with pytest.raises(KeyError):
        host.plot("ghost")
    with pytest.raises(ValueError):
        host.plot(backend="crayon")


def test_plot_single_model_host_rejects_model_and_aggregates_cells():
    frame = _series_frame().drop(columns=["model", "lake"])
    host = _FakeSeriesHost(frame)  # single-model: no model axis
    with pytest.raises(ValueError):
        host.plot("a")
    fig = host.plot()  # cells aggregate -> a single total line
    assert len(fig.data) == 1
    total = fig.data[0]
    assert list(total.y) == [pytest.approx(1 + 2 + 5 + 6), pytest.approx(3 + 4 + 7 + 8)]


def test_plot_requires_period_column():
    host = _FakeSeriesHost(pd.DataFrame({"cell": [1], "q": [1.0]}))
    with pytest.raises(ValueError):
        host.plot()


def test_plot_agg_override_and_layer_split():
    frame = _series_frame().drop(columns=["model", "lake"])
    frame["layer"] = [0, 1] * 4  # two layers -> split into layer lines
    host = _FakeSeriesHost(frame)
    fig = host.plot(agg="mean")
    assert {trace.name for trace in fig.data} == {"L0", "L1"}


# --- field= sugar dispatches plot like map/mosaic/animate ----------------------
def test_field_mappable_plot_dispatch():
    host = _FakeSeriesHost(_series_frame(), models=["a", "b"])

    class _Namespace(FieldMappable):
        _default_field = "q"

        def _field_names(self):
            return ["q"]

        @property
        def q(self):
            return host

    namespace = _Namespace()
    fig = namespace.plot()
    assert isinstance(fig, go.Figure) and len(fig.data) == 4
    with pytest.raises(ValueError):
        namespace.plot(field="stage")


# --- heads leaf gains the grammar verbs at the class level ---------------------
def test_headsplus_is_a_spatial_view_leaf():
    from myflopy.modflow.mf6.headsplus import HeadsPlus

    assert issubclass(HeadsPlus, SpatialView)
    # the grammar's plot deliberately shadows flopy's legacy LayerFile.plot
    assert HeadsPlus.plot is SpatialView.plot
    for verb in ("get", "summary", "map", "plot", "xs", "mosaic", "animate"):
        assert callable(getattr(HeadsPlus, verb))


def test_group_heads_is_a_group_spatial_view_leaf():
    from myflopy.project.model_group import GroupHeads, _GroupSpatialView

    assert issubclass(GroupHeads, _GroupSpatialView)
    for verb in ("get", "compare", "map", "plot", "xs", "mosaic", "animate"):
        assert callable(getattr(GroupHeads, verb))


def test_heads_result_diff_has_plot_and_xs():
    from myflopy.project.model_results_diff import HeadsResultDiff

    for verb in ("get", "summary", "map", "plot", "xs", "mosaic", "animate"):
        assert callable(getattr(HeadsResultDiff, verb))


# --- real renders on the canonical model ---------------------------------------
@pytest.mark.slow
def test_plot_and_xs_verbs_render_on_canonical(canonical_run):
    """The plot/xs panel verbs render real figures on all three surfaces."""

    import myflopy as mf
    from myflopy.project.model_group import ModelGroup

    model = canonical_run

    # -- single model: heads leaf ------------------------------------------
    heads = model.hds.get()
    assert list(heads.columns) == ["per", "layer", "cell", "head"]
    assert heads["per"].nunique() >= 1 and not heads["head"].isna().all()

    summary = model.hds.summary()
    assert int(summary.loc[0, "records"]) == len(heads)

    hydrograph = model.hds.plot(cells=[0, 1])
    assert isinstance(hydrograph, go.Figure) and len(hydrograph.data) >= 2
    assert isinstance(model.hds.plot(backend="mpl"), Figure)

    # section line across the model interior
    bounds = model.vor.gdf_vorPolys.total_bounds  # xmin, ymin, xmax, ymax
    line = LineString(
        [
            (bounds[0] + 0.25 * (bounds[2] - bounds[0]), bounds[1] + 0.25 * (bounds[3] - bounds[1])),
            (bounds[0] + 0.75 * (bounds[2] - bounds[0]), bounds[1] + 0.75 * (bounds[3] - bounds[1])),
        ]
    )
    section = model.hds.xs(line=line)
    assert isinstance(section, go.Figure) and len(section.data) >= 2  # head + top
    assert isinstance(model.hds.xs(line=line, backend="mpl"), Figure)

    # -- single model: package leaves --------------------------------------
    ghb_series = model.packages.ghb.results.plot()  # flux by period (per layer)
    assert isinstance(ghb_series, go.Figure) and len(ghb_series.data) >= 1
    # no redundant model-name prefix on a single-model surface
    assert all(model.name not in (trace.name or "") for trace in ghb_series.data)

    lak_stage_series = model.packages.lak.results.plot(field="stage")
    assert isinstance(lak_stage_series, go.Figure) and len(lak_stage_series.data) >= 1

    stage_change = model.packages.lak.results.stage_change.plot()
    assert isinstance(stage_change, go.Figure)

    # -- group surface -------------------------------------------------------
    reloaded = mf.load_mf6_run(model.workspace)
    group = ModelGroup({"a": model, "b": reloaded}, reference="a")

    assert group.hds.map() is not None  # reference heads choropleth
    nlay = int(model.gwf.modelgrid.nlay)
    group_series = group.hds.plot()  # one line per model x layer
    assert isinstance(group_series, go.Figure)
    assert len(group_series.data) == 2 * nlay
    group_layer0 = group.hds.plot(layer=0)  # layer filter -> one line per model
    assert len(group_layer0.data) == 2
    group_section = group.hds.xs(line=line)
    assert isinstance(group_section, go.Figure)
    assert len(group_section.data) >= 3  # two member profiles + model top

    group_stage = group.packages.lak.results.plot(field="stage")
    assert isinstance(group_stage, go.Figure) and len(group_stage.data) >= 2

    # -- diff surface ---------------------------------------------------------
    diff = group.diff()
    diff_series = diff.hds.plot(layer=0)  # mean dHead line per model
    assert isinstance(diff_series, go.Figure) and len(diff_series.data) == 1
    diff_section = diff.hds.xs(line=line)
    assert isinstance(diff_section, go.Figure) and len(diff_section.data) >= 3
    diff_q_series = diff.packages.ghb.results.q.plot()
    assert isinstance(diff_q_series, go.Figure) and len(diff_q_series.data) >= 1
