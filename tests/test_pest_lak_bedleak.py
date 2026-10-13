"""Lake-bed leakance is a calibration target, one multiplier per LAKE by default.

Ledger 558. ``lak.bedleak`` (alias ``bedleak``) scales the LAK
CONNECTIONDATA ``bedleak`` column, keyed by ``(lakeno, iconn)`` -- unique per row, where
a lake connects to one cell twice -- with ``use_rows`` computed from the staged package,
never pyEMU's all-rows fallback. ``style="feature"`` (the default) makes one parameter per lake, named after
its boundname; ``features=`` picks lakes by boundname (any case) or zero-based number.
The model is ``tests/pest_feature_models.py``'s.
"""

from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6.pest import native_parameters as np_mod
from tests.pest_feature_models import (
    LAKE_CONNECTIONS,
    build,
    build_model,
    calibration,
    forward_run,
    records,
)

pyemu = pytest.importorskip("pyemu")
flopy = pytest.importorskip("flopy")

LAK_FILE = "m.lak_connectiondata.txt"


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    root = tmp_path_factory.mktemp("lak_bedleak")
    return {
        "named": build_model(root / "named"),
        "unnamed": build_model(root / "unnamed", lake_boundnames=False),
        "none": build_model(root / "none", none_bedleak=True),
        "two": build_model(root / "two", two_lak=True),
    }


def _bedleak(worker):
    return [float(value) for value in records(worker / LAK_FILE).iloc[:, 5]]


BASE = [row[4] for row in LAKE_CONNECTIONS]


# -- one parameter per lake ----------------------------------------------------------


def test_the_default_is_one_parameter_per_lake_and_it_reaches_only_its_lake(runs, tmp_path):
    """Keevie x10 multiplies Keevie's three connections -- two of them to ONE cell --
    and leaves Marjorie's alone; MF6 runs. Each row is one ``(lakeno, iconn)`` key, so
    the duplicated cell is two records."""

    cal = calibration(runs["named"], tmp_path)
    spec = cal.parameterize("bedleak", bounds=(0.1, 10.0))
    assert spec.style == "feature" and spec.name == "lakbedleak"
    pst = build(cal)
    groups = sorted(set(pst.parameter_data["pargp"]))
    assert groups == ["lakbedleak_keevie", "lakbedleak_marjorie"]
    assert len(pst.parameter_data) == 2
    assert [record["features"] for record in spec.feature_sets] == [[0], [1]]
    assert [record["rows"] for record in spec.feature_sets] == [3, 2]

    worker = forward_run(cal, pst, tmp_path, {"lakbedleak_keevie": 10.0})
    assert np.allclose(_bedleak(worker), [1.0, 2.0, 2.5, 0.3, 0.35])
    worker = forward_run(cal, pst, tmp_path, {"lakbedleak_marjorie": 0.5})
    assert np.allclose(_bedleak(worker), [0.1, 0.2, 0.25, 0.15, 0.175])
    # The claktype column, text, survives the rewrite.
    assert list(records(worker / LAK_FILE).iloc[:, 4]) == [row[3] for row in LAKE_CONNECTIONS]


def test_features_pick_lakes_by_boundname_any_case_or_by_zero_based_number(runs, tmp_path):
    cal = calibration(runs["named"], tmp_path)
    keevie = cal.parameterize("bedleak", features=["KEEVIE"], name="kv")
    marjorie = cal.parameterize("lak.bedleak", features=1, name="mj")
    shared = cal.parameterize("bedleak", style="constant", features=[0, "marjorie"],
                              name="both")
    assert keevie.features == ("KEEVIE",) and marjorie.features == (1,)
    pst = build(cal)
    assert sorted(set(pst.parameter_data["pargp"])) == ["both", "kv_keevie", "mj_marjorie"]
    assert shared.feature_sets == [{"pargp": "both", "label": None, "package": "lak",
                                    "features": [0, 1], "rows": 5}]
    worker = forward_run(cal, pst, tmp_path, {"kv_keevie": 2.0, "mj_marjorie": 3.0,
                                              "both": 10.0})
    assert np.allclose(_bedleak(worker), [2.0, 4.0, 5.0, 9.0, 10.5])


def test_a_constant_scales_every_lake_together(runs, tmp_path):
    cal = calibration(runs["named"], tmp_path)
    spec = cal.parameterize("bedleak", style="constant")
    pst = build(cal)
    assert len(pst.parameter_data) == 1 and spec.feature_sets[0]["rows"] == 5
    worker = forward_run(cal, pst, tmp_path, {"lakbedleak": 2.0})
    assert np.allclose(_bedleak(worker), [2 * value for value in BASE])


def test_a_lake_with_no_boundname_is_named_by_its_number(runs, tmp_path):
    cal = calibration(runs["unnamed"], tmp_path)
    cal.parameterize("bedleak")
    pst = build(cal)
    assert sorted(set(pst.parameter_data["pargp"])) == ["lakbedleak_0", "lakbedleak_1"]
    with pytest.raises(ValueError, match=r"boundname 'keevie'.*lakes 0 \.\. 1.*boundnames: none"):
        cal.parameterize("bedleak", features=["keevie"], name="kv")


# -- refusals -------------------------------------------------------------------------


def test_a_lake_the_package_does_not_have_is_refused_at_the_call(runs, tmp_path):
    cal = calibration(runs["named"], tmp_path)
    with pytest.raises(ValueError, match=(
            r"boundname 'kevie'.*lakes 0 \.\. 1 \(numbered from 0\); boundnames: "
            r"'keevie', 'marjorie'")):
        cal.parameterize("bedleak", features=["kevie"])
    with pytest.raises(ValueError, match=r"lake 2 \(one past the last"):
        cal.parameterize("bedleak", features=[2])
    with pytest.raises(ValueError, match="names lake 0 twice"):
        cal.parameterize("bedleak", features=["keevie", 0])
    assert cal._native_parameter_specs == []


def test_a_none_bedleak_on_a_selected_lake_is_refused_and_elsewhere_is_left_alone(
        runs, tmp_path):
    """MF6 reads ``NONE`` as no lakebed -- nothing to multiply. On a lake not selected it
    is written back as it was and MF6 runs."""

    cal = calibration(runs["none"], tmp_path)
    cal.parameterize("bedleak")
    with pytest.raises(ValueError, match=r"lake 1 connection 1 \(both zero-based\) has bedleak NONE"):
        build(cal)

    cal = calibration(runs["none"], tmp_path, name="kv")
    cal.parameterize("bedleak", features=["keevie"])
    pst = build(cal, "kv")
    worker = forward_run(cal, pst, tmp_path, {"lakbedleak_keevie": 10.0})
    written = list(records(worker / LAK_FILE).iloc[:, 5])
    assert np.allclose([float(value) for value in written[:4]], [1.0, 2.0, 2.5, 0.3])
    assert written[4].lower() == "none"


@pytest.mark.parametrize("kwargs,message", [
    ({"style": "grid"}, r"style='grid' on target 'lak.bedleak': it takes 'feature', 'constant'"),
    ({"style": "zone"}, r"style='zone' on target 'lak.bedleak'.*'feature' is one multiplier"),
    ({"style": "pilotpoints"}, r"style='pilotpoints' on target 'lak.bedleak'"),
    ({"correlation": 500.0}, r"correlation=500.0 with style='feature'.*not located at a cell"),
    ({"zones": np.ones(20, dtype=int)}, r"zones= on 'lak.bedleak'.*features="),
    ({"capture": True}, r"capture=True on 'lak.bedleak'.*deferred"),
])
def test_an_argument_that_cannot_act_on_a_lake_is_refused(kwargs, message):
    with pytest.raises(ValueError, match=message):
        np_mod.NativeParameterSpec(target="bedleak", **kwargs)


def test_features_and_the_feature_style_are_refused_where_there_are_no_features():
    with pytest.raises(ValueError, match=r"features=\[0\] on target 'k'.*\['lak.bedleak', 'sfr.rhk'\]"):
        np_mod.NativeParameterSpec(target="k", features=[0])
    with pytest.raises(ValueError, match=r"style='feature' on target 'riv.cond'"):
        np_mod.NativeParameterSpec(target="riv.cond", style="feature")
    for bad, why in ((["  "], "is blank"), ([0.5], "not a whole number"),
                     ([-1], "negative"), ([True], "not a number"), ([], "names no lake")):
        with pytest.raises(ValueError, match=why):
            np_mod.NativeParameterSpec(target="bedleak", features=bad)
    assert np_mod.NativeParameterSpec(target="bedleak", features=" 3 ").features == (3,)


def test_two_lakes_whose_names_make_one_pest_name_are_refused():
    data = np.rec.fromrecords([(0, 101.0, 1, "North Lake"), (1, 101.0, 1, "north-lake "),
                               (2, 101.0, 1, "NORTH  LAKE")],
                              names="ifno,strt,nlakeconn,boundname")
    package = SimpleNamespace(package_name="lak",
                              packagedata=SimpleNamespace(get_data=lambda: data))
    spec = np_mod.NativeParameterSpec(target="bedleak", features=["north lake", "north-lake"])
    with pytest.raises(ValueError, match="names lake 0 twice"):
        # "North Lake" and "NORTH  LAKE" are two lakes; "north lake" names lake 0 alone.
        np_mod.resolve_features(
            np_mod.NativeParameterSpec(target="bedleak", features=["north lake", 0]), package)
    sets = np_mod.resolve_features(spec, package)
    assert [(item.label, item.features) for item in sets] == [
        ("north_lake", (0,)), ("north-lake", (1,))]
    with pytest.raises(ValueError, match=r"'North Lake' and 'NORTH  LAKE' both make the PEST "
                                         r"name lakbedleak_north_lake"):
        np_mod.resolve_features(np_mod.NativeParameterSpec(target="bedleak"), package)


def test_several_lake_packages_need_package_and_one_is_scaled(runs, tmp_path):
    cal = calibration(runs["two"], tmp_path)
    with pytest.raises(mf.SeveralPackagesError, match=r"'lak', 'ponds'"):
        cal.parameterize("bedleak")
    with pytest.raises(ValueError, match=r"package=\('lak', 'ponds'\) with style='feature'"):
        cal.parameterize("bedleak", package=("lak", "ponds"))
    spec = cal.parameterize("bedleak", package="ponds")
    pst = build(cal)
    assert spec.name == "lakbedleak_ponds" and spec.resolved_files == ["m.ponds_connectiondata.txt"]
    assert sorted(set(pst.parameter_data["pargp"])) == ["lakbedleak_ponds_0"]
    worker = forward_run(cal, pst, tmp_path, {"lakbedleak_ponds_0": 4.0})
    assert np.allclose(_bedleak(worker), BASE)
    assert float(records(worker / "m.ponds_connectiondata.txt").iloc[0, 5]) == pytest.approx(2.0)


def test_the_metadata_records_the_features_and_their_parameter_sets(runs, tmp_path):
    cal = calibration(runs["named"], tmp_path)
    cal.parameterize("bedleak", features=["marjorie", "Keevie"])
    build(cal)
    meta = json.loads((cal.template_workspace / "myflopy_pest_metadata.json").read_text())
    record = next(item for item in meta["parameters"] if item["name"] == "lakbedleak")
    assert record["features"] == ["marjorie", "Keevie"]
    assert [(item["pargp"], item["features"]) for item in record["feature_sets"]] == [
        ("lakbedleak_marjorie", [1]), ("lakbedleak_keevie", [0])]
    assert record["style"] == "feature" and record["files"][0]["filename"] == LAK_FILE
    text = str(cal.settings())
    assert "features=['marjorie', 'Keevie']" in text
    assert "lakbedleak_keevie: 3 rows of package 'lak' (features [0])" in text
    assert cal.settings().parameter_frame()["features"].tolist() == [["marjorie", "Keevie"]]


def test_a_template_scaling_other_rows_than_selected_is_refused(tmp_path):
    """pyEMU sets up EVERY row when no ``use_rows`` matches (measured); the build reads
    the template pyEMU wrote back and refuses one that is not the selection."""

    tpl = tmp_path / "x.csv.tpl"
    tpl.write_text('ptf ~\n,0,1,sidx,idx_strs,pargp5,parval1_5,5\n'
                   '0,0,0,"(0, 0)",a,k,1.0,~ k ~\n1,1,0,"(1, 0)",b,k,1.0,~ k ~\n')
    project = SimpleNamespace(template_workspace=tmp_path)
    frame = pd.DataFrame({"tpl_filename": [str(tpl)]})
    spec = np_mod.NativeParameterSpec(target="bedleak")
    np_mod._refuse_rows_pyemu_did_not_select(project, frame, {(0, 0), (1, 0)}, spec)
    with pytest.raises(RuntimeError, match=r"scales 2 row\(s\) where 1 were selected"):
        np_mod._refuse_rows_pyemu_did_not_select(project, frame, {(0, 0)}, spec)


# -- a real PESTPP-IES run ------------------------------------------------------------


@pytest.mark.slow
def test_pestpp_ies_runs_per_lake_leakance_and_per_stream_rhk(runs, tmp_path):
    """PESTPP-IES runs a template carrying one leakance multiplier per lake and one rhk
    multiplier per named stream: every realization's forward run applies them and MF6
    runs, and the posterior holds all four."""

    cal = calibration(runs["named"], tmp_path, name="ies")
    cal.parameterize("bedleak", bounds=(0.1, 10.0))
    cal.parameterize("rhk", style="feature", bounds=(0.2, 5.0))
    build(cal, "ies")
    ies = cal.run_ies(reals=4, iterations=1, workers=2, noise=False,
                      ies_lambda_mults=1.0, lambda_scale_fac=1.0)
    groups = set(cal.pst.parameter_data["pargp"])
    assert groups == {"lakbedleak_keevie", "lakbedleak_marjorie", "sfrrhk_creek_a",
                      "sfrrhk_creek_b"}
    posterior = ies.posterior_parameters._df
    assert set(cal.pst.adj_par_names) <= {str(name).lower() for name in posterior.columns}
    assert ies.posterior_iteration >= 1
    # Each realization drew its own multipliers, so they differ across realizations.
    assert (posterior.std(axis=0) > 0).all()
