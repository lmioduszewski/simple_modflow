"""Streambed K is a calibration target: one multiplier, one per reach, per named stream or per zone.

Ledger 558. ``sfr.rhk`` (alias ``rhk``) scales the SFR PACKAGEDATA ``rhk``
column. ``"constant"`` (the default) scales every reach; ``"grid"`` is one parameter per
reach, keyed ``(ifno, cellid)`` so pyEMU locates it at its cell; ``"feature"`` one per
BOUNDNAME (a named stream); ``"zone"`` one per ``zones=`` label of the reach's cell. A
reach whose cellid is ``NONE`` is refused: its row is one column short, and pyEMU,
reading by position, multiplied its Manning's n (measured). The model is
``tests/pest_feature_models.py``'s: reaches 0-2 ``creek_a``, 3-4 ``Creek_B``, reaches 1
and 2 in one cell.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6.pest import native_parameters as np_mod
from tests.pest_feature_models import (
    REACHES,
    build,
    build_model,
    calibration,
    forward_run,
    records,
)

pyemu = pytest.importorskip("pyemu")
flopy = pytest.importorskip("flopy")

SFR_FILE = "m.sfr_packagedata.txt"
RHK = [row[2] for row in REACHES]
MAN = [row[3] for row in REACHES]


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    root = tmp_path_factory.mktemp("sfr_rhk")
    return {
        "named": build_model(root / "named"),
        "unnamed": build_model(root / "unnamed", sfr_boundnames=False),
        "none": build_model(root / "none", none_reach=True),
        "partly": build_model(root / "partly", unnamed_reach=True),
    }


def _columns(worker):
    rows = records(worker / SFR_FILE)
    return [float(value) for value in rows.iloc[:, 8]], [float(value) for value in rows.iloc[:, 9]]


def test_a_constant_doubles_every_reach_and_leaves_mannings_n_alone(runs, tmp_path):
    cal = calibration(runs["named"], tmp_path)
    spec = cal.parameterize("rhk")
    assert spec.style == "constant" and spec.recipe.canonical == "sfr.rhk"
    pst = build(cal)
    assert len(pst.parameter_data) == 1
    worker = forward_run(cal, pst, tmp_path, {"sfrrhk": 2.0})
    rhk, man = _columns(worker)
    assert np.allclose(rhk, [2 * value for value in RHK]) and np.allclose(man, MAN)


def test_grid_is_one_parameter_per_reach_and_each_reaches_its_own(runs, tmp_path):
    """Reaches 1 and 2 share a cell: keyed by the cell alone they would be one record."""

    cal = calibration(runs["named"], tmp_path)
    cal.parameterize("rhk", style="grid")
    pst = build(cal)
    par = pst.parameter_data
    assert len(par) == 5
    reach_two = [name for name in par.index if "idx0:2_" in name]
    assert len(reach_two) == 1
    par.loc[reach_two, "parval1"] = 3.0
    worker = forward_run(cal, pst, tmp_path, {})
    rhk, _man = _columns(worker)
    assert np.allclose(rhk, [0.5, 0.6, 2.1, 0.8, 0.9])


def test_feature_is_one_parameter_per_named_stream(runs, tmp_path):
    cal = calibration(runs["named"], tmp_path)
    spec = cal.parameterize("rhk", style="feature", bounds=(0.1, 10.0))
    pst = build(cal)
    assert sorted(set(pst.parameter_data["pargp"])) == ["sfrrhk_creek_a", "sfrrhk_creek_b"]
    assert [(item["label"], item["features"]) for item in spec.feature_sets] == [
        ("creek_a", [0, 1, 2]), ("creek_b", [3, 4])]
    worker = forward_run(cal, pst, tmp_path, {"sfrrhk_creek_a": 2.0, "sfrrhk_creek_b": 0.5})
    rhk, man = _columns(worker)
    assert np.allclose(rhk, [1.0, 1.2, 1.4, 0.4, 0.45]) and np.allclose(man, MAN)


def test_features_name_streams_and_reaches(runs, tmp_path):
    cal = calibration(runs["named"], tmp_path)
    cal.parameterize("rhk", style="feature", features=["CREEK_B", 0], name="pick")
    cal.parameterize("rhk", style="grid", features=["creek_a"], name="cells")
    pst = build(cal)
    groups = set(pst.parameter_data["pargp"])
    assert {"pick_creek_b", "pick_0"} <= groups
    assert int((pst.parameter_data["pargp"] == "cells").sum()) == 3
    worker = forward_run(cal, pst, tmp_path, {"pick_creek_b": 2.0, "pick_0": 3.0})
    rhk, _man = _columns(worker)
    assert np.allclose(rhk, [1.5, 0.6, 0.7, 1.6, 1.8])
    with pytest.raises(ValueError, match=r"names reach 3 twice -- through 'creek_b' and '3'"):
        cal.parameterize("rhk", style="feature", features=["creek_b", 3], name="twice")
    with pytest.raises(ValueError, match=r"reach 5 \(one past the last.*reaches 0 \.\. 4"):
        cal.parameterize("rhk", features=[5], name="past")


def test_zone_is_one_parameter_per_label_of_the_reachs_cell(runs, tmp_path):
    """Reaches 1 and 2 share cell 11 (zone 2) and take one multiplier; pyEMU's cell key
    holds them as one record, which applies to both."""

    cal = calibration(runs["named"], tmp_path)
    zones = np.ones(20, dtype=int)
    zones[11] = 2
    cal.parameterize("rhk", style="zone", zones=zones)
    pst = build(cal)
    par = pst.parameter_data
    assert len(par) == 2
    par.loc[[name for name in par.index if "zone:2" in name], "parval1"] = 3.0
    worker = forward_run(cal, pst, tmp_path, {})
    rhk, _man = _columns(worker)
    assert np.allclose(rhk, [0.5, 1.8, 2.1, 0.8, 0.9])


def test_every_style_stacks_on_one_file_and_each_multiplies(runs, tmp_path):
    """A constant, a per-reach grid, a zone and a per-stream feature on one SFR file --
    four keys (``ifno``; ``ifno, cell``; the cell; ``ifno``) -- and a constant under the
    per-lake leakance: every multiplier on a row multiplies."""

    cal = calibration(runs["named"], tmp_path)
    cal.parameterize("rhk", name="rc")
    cal.parameterize("rhk", style="grid", name="rg")
    zones = np.ones(20, dtype=int)
    zones[11] = 2
    cal.parameterize("rhk", style="zone", zones=zones, name="rz")
    cal.parameterize("rhk", style="feature", name="rf")
    cal.parameterize("bedleak", style="constant", name="lc")
    cal.parameterize("bedleak", name="lf")
    pst = build(cal)
    par = pst.parameter_data
    par.loc[[name for name in par.index if "rg_" in name and "idx0:2_" in name], "parval1"] = 3.0
    par.loc[[name for name in par.index if "rz_" in name and "zone:2" in name], "parval1"] = 5.0
    worker = forward_run(cal, pst, tmp_path, {"rc": 2.0, "rf_creek_b": 7.0, "lc": 2.0,
                                              "lf_marjorie": 3.0})
    rhk, man = _columns(worker)
    assert np.allclose(rhk, [1.0, 6.0, 21.0, 11.2, 12.6]) and np.allclose(man, MAN)
    bedleak = records(worker / "m.lak_connectiondata.txt").iloc[:, 5].astype(float)
    assert np.allclose(bedleak, [0.2, 0.4, 0.5, 1.8, 2.1])


# -- refusals ---------------------------------------------------------------------------


def test_a_reach_with_no_cell_is_refused_where_pyemu_multiplied_its_mannings_n(runs, tmp_path):
    """The regression: before the refusal pyEMU built and scaled column 8 of every row,
    which on the ``NONE`` reach, one column shorter, is Manning's n (0.034 became 0.068
    at x2, measured)."""

    cal = calibration(runs["none"], tmp_path)
    cal.parameterize("rhk")
    with pytest.raises(ValueError, match=r"reach\(es\) \[4\] \(zero-based\) have cellid NONE.*"
                                         r"Manning's n"):
        build(cal)


def test_feature_on_streams_with_no_boundnames_is_refused_naming_grid(runs, tmp_path):
    cal = calibration(runs["unnamed"], tmp_path)
    with pytest.raises(ValueError, match=r"declares no BOUNDNAMES.*style='grid'"):
        cal.parameterize("rhk", style="feature")
    spec = cal.parameterize("rhk", style="feature", features=[0, 3])
    pst = build(cal)
    assert sorted(set(pst.parameter_data["pargp"])) == ["sfrrhk_0", "sfrrhk_3"]
    assert [item["rows"] for item in spec.feature_sets] == [1, 1]


def test_a_reach_no_named_stream_covers_is_left_alone_and_said_to_be(runs, tmp_path):
    """Reach 4 carries no boundname: one parameter per named stream covers reaches 0-3,
    and the build says, as a Python warning (ledger 288), that reach 4 keeps its rhk."""

    cal = calibration(runs["partly"], tmp_path)
    cal.parameterize("rhk", style="feature")
    with pytest.warns(UserWarning, match=r"1 reach\(es\) of SFR package 'sfr' carry no "
                                         r"boundname -- \[4\]"):
        pst = cal.build("feat.pst", noptmax=0)
    worker = forward_run(cal, pst, tmp_path, {"sfrrhk_creek_a": 2.0, "sfrrhk_creek_b": 2.0})
    rhk, _man = _columns(worker)
    assert np.allclose(rhk, [1.0, 1.2, 1.4, 1.6, 0.9])


@pytest.mark.parametrize("kwargs,message", [
    ({"style": "pilotpoints"}, r"style='pilotpoints' on target 'sfr.rhk': it takes"),
    ({"style": "feature", "correlation": 50.0}, r"style='feature'.*not located at a cell"),
    ({"style": "grid", "zones": np.ones(20, dtype=int)}, r"zones= on 'sfr.rhk'.*style='zone'"),
    ({"style": "zone", "zones": np.ones(20, dtype=int), "features": [0]},
     r"features= with style='zone'"),
    ({"package": ("sfr", "creeks"), "features": [0]}, r"numbered and named within ONE package"),
])
def test_an_argument_that_cannot_act_on_a_reach_is_refused(kwargs, message):
    with pytest.raises(ValueError, match=message):
        np_mod.NativeParameterSpec(target="rhk", **kwargs)


def test_a_correlated_grid_over_reaches_sharing_a_cell_is_refused(runs, tmp_path):
    """pyEMU places reaches 1 and 2 at one point and then draws EVERY reach uncorrelated
    (off-diagonal within 0.03 of 0, measured); without the shared cell the same draw
    correlates neighbours at 0.82. ``docs/reproducers/sfr_grid_shared_cell_prior.py``."""

    cal = calibration(runs["named"], tmp_path)
    cal.parameterize("rhk", style="grid", correlation=50.0)
    with pytest.raises(ValueError, match=r"reaches \[1, 2\] \(zero-based\) lie in one model column"):
        build(cal)

    cal = calibration(runs["named"], tmp_path, name="apart")
    cal.parameterize("rhk", style="grid", correlation=50.0, features=[0, 1, 3, 4])
    build(cal, "apart")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        drawn = cal.draw_prior(500, use_specsim=False)
    frame = drawn._df if hasattr(drawn, "_df") else drawn
    logs = np.log10(frame[[name for name in frame.columns if "sfrrhk" in name]])
    assert logs.corr().to_numpy()[0, 1] > 0.6


def test_the_canonical_models_sfr_scales_its_rhk_column(canonical_run, tmp_path):
    """On the library's own SFR builder output (the canonical model) the column pyEMU
    scales is the package's ``rhk``: the template's original file holds, in that column,
    exactly what FloPy reads as ``rhk``, so the column (8 on DISV) is the builder's rhk.
    The shared run is only read."""

    import geopandas as gpd

    grid = canonical_run.gwf.modelgrid
    point = gpd.GeoDataFrame({"name": ["W"], "layer": [0]}, geometry=gpd.points_from_xy(
        [float(grid.xcellcenters[0])], [float(grid.ycellcenters[0])]))
    cal = canonical_run.pest("canon", workspace=tmp_path / "canon", start_datetime="2024-01-01")
    cal.observe(mf.HeadTargets(locations=point, values=pd.DataFrame({"time": [0], "W": [1.0]})))
    spec = cal.parameterize("rhk")
    build(cal, "canon")
    # pyEMU keeps its copy of the original comma-separated.
    org = pd.read_csv(cal.template_workspace / "org" / spec.resolved_files[0], header=None,
                      sep=r"[\s,]+", engine="python")
    flopy_rhk = pd.DataFrame(canonical_run.gwf.sfr.packagedata.get_data())["rhk"]
    assert len(org) == len(flopy_rhk) > 1
    assert np.allclose(org.iloc[:, 8].astype(float).to_numpy(), flopy_rhk.to_numpy(dtype=float))
