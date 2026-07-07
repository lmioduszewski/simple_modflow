"""Tests for the ModelDiff connection-geometry tier (Phase 3: LAK / SFR).

The connection/reach table builders are patched so the multiset-difference
engine runs on small, controlled geometry tables without a real MF6 model.
"""

from __future__ import annotations

import pandas as pd
import pytest

from myflopy.project import model_diff as md
from myflopy.project.model_config import ModelConfig
from myflopy.project.model_diff import LakConnectionDiff, SfrReachDiff
from myflopy.project.model_group import ModelGroup

_LAK_COLS = ["lake", "layer", "cell", "claktype", "belev", "telev", "connlen", "connwidth"]
_SFR_COLS = ["reach", "layer", "cell", "rlen"]


def _lak(name, rows):
    return pd.DataFrame(rows, columns=_LAK_COLS).assign(model=name, package="lak")


def _sfr(name, rows):
    return pd.DataFrame(rows, columns=_SFR_COLS).assign(model=name)


def _empty_config():
    return ModelConfig(pd.DataFrame(columns=["section", "setting", "value"]))


class _FakeModel:
    def __init__(self, name, package_names=("LAK", "SFR")):
        self.name = name
        self.package_names = list(package_names)
        self.config = _empty_config()


@pytest.fixture
def connection_tables(monkeypatch):
    lak: dict[str, pd.DataFrame] = {}
    sfr: dict[str, pd.DataFrame] = {}

    def lak_stub(model):
        if model.name not in lak:
            raise KeyError(f"lak not attached to {model.name}")
        return lak[model.name].copy()

    def sfr_stub(model):
        if model.name not in sfr:
            raise KeyError(f"sfr not attached to {model.name}")
        return sfr[model.name].copy()

    monkeypatch.setattr(md, "build_lak_connection_table", lak_stub)
    monkeypatch.setattr(md, "build_sfr_reach_table", sfr_stub)
    return lak, sfr


def _group(names, *, reference, package_names=("LAK", "SFR")):
    models = {name: _FakeModel(name, package_names) for name in names}
    return ModelGroup(models, reference=reference)


# --- LAK connections ---------------------------------------------------------
def test_identical_lake_networks_have_no_connection_difference(connection_tables):
    lak, _ = connection_tables
    rows = [
        (0, 0, 10, "VERTICAL", 0.0, 0.0, 0.0, 0.0),
        (0, 0, 11, "HORIZONTAL", 5.0, 8.0, 50.0, 100.0),
    ]
    lak["ref"] = _lak("ref", rows)
    lak["twin"] = _lak("twin", list(rows))

    diff = _group(["ref", "twin"], reference="ref").diff()
    assert diff.packages.lak.inputs.connections().empty
    assert bool(diff.packages.lak.inputs.summary().iloc[0]["identical"])


def test_lake_connection_geometry_change_shows_as_remove_and_add(connection_tables):
    lak, _ = connection_tables
    # Lake 0 identical; lake 1 has a horizontal connection whose connlen changed.
    lak["ref"] = _lak(
        "ref",
        [
            (0, 0, 10, "VERTICAL", 0.0, 0.0, 0.0, 0.0),
            (1, 0, 20, "HORIZONTAL", 5.0, 8.0, 50.0, 100.0),
        ],
    )
    lak["variant"] = _lak(
        "variant",
        [
            (0, 0, 10, "VERTICAL", 0.0, 0.0, 0.0, 0.0),
            (1, 0, 20, "HORIZONTAL", 5.0, 8.0, 75.0, 100.0),  # connlen 50 -> 75
        ],
    )

    diff = _group(["ref", "variant"], reference="ref").diff()
    connections = diff.packages.lak.inputs.connections()
    # Lake 0 unaffected; only lake 1 cell 20 differs, as one remove + one add.
    assert set(connections["lake"]) == {1}
    memberships = dict(zip(connections["connlen"], connections["membership"]))
    assert memberships[50.0] == "only_in_reference"
    assert memberships[75.0] == "only_in_model"

    row = diff.packages.lak.inputs.summary().iloc[0]
    assert row["only_in_reference"] == 1
    assert row["only_in_model"] == 1
    assert not row["identical"]


def test_lak_absent_in_one_model(connection_tables):
    lak, _ = connection_tables
    lak["ref"] = _lak("ref", [(0, 0, 10, "VERTICAL", 0.0, 0.0, 0.0, 0.0)])
    # variant has no LAK
    models = {
        "ref": _FakeModel("ref", ["LAK"]),
        "variant": _FakeModel("variant", []),
    }
    diff = ModelGroup(models, reference="ref").diff()
    row = diff.packages.lak.inputs.summary().iloc[0]
    assert row["present_in_reference"]
    assert not row["present_in_model"]
    assert row["only_in_reference"] == 1
    assert not row["identical"]


def test_float_geometry_noise_below_round_to_is_ignored(connection_tables):
    lak, _ = connection_tables
    lak["ref"] = _lak("ref", [(0, 0, 10, "HORIZONTAL", 5.0, 8.0, 50.0000001, 100.0)])
    lak["variant"] = _lak("variant", [(0, 0, 10, "HORIZONTAL", 5.0, 8.0, 50.0000002, 100.0)])
    diff = _group(["ref", "variant"], reference="ref").diff()
    assert diff.packages.lak.inputs.connections(round_to=3).empty


# --- SFR reaches -------------------------------------------------------------
def test_sfr_reach_length_change_detected(connection_tables):
    _, sfr = connection_tables
    sfr["ref"] = _sfr("ref", [(0, 0, 5, 100.0), (1, 0, 6, 100.0)])
    sfr["variant"] = _sfr("variant", [(0, 0, 5, 100.0), (1, 0, 6, 140.0)])  # rlen change

    diff = _group(["ref", "variant"], reference="ref").diff()
    reaches = diff.packages.sfr.inputs.reaches()
    assert set(reaches["reach"]) == {1}
    assert set(reaches["membership"]) == {"only_in_reference", "only_in_model"}
    # .connections() is an alias for SFR
    assert not diff.packages.sfr.inputs.connections().empty
    assert not bool(diff.packages.sfr.inputs.summary().iloc[0]["identical"])


# --- namespace, report, errors ----------------------------------------------
def test_namespace_returns_connection_diff_types(connection_tables):
    lak, sfr = connection_tables
    lak["ref"] = _lak("ref", [(0, 0, 10, "VERTICAL", 0.0, 0.0, 0.0, 0.0)])
    lak["variant"] = _lak("variant", [(0, 0, 10, "VERTICAL", 0.0, 0.0, 0.0, 0.0)])
    sfr["ref"] = _sfr("ref", [(0, 0, 5, 100.0)])
    sfr["variant"] = _sfr("variant", [(0, 0, 5, 100.0)])
    diff = _group(["ref", "variant"], reference="ref").diff()
    # the node's inputs tier carries the connection-geometry diff
    assert isinstance(diff.packages.lak.inputs, LakConnectionDiff)
    assert isinstance(diff.packages.sfr.inputs, SfrReachDiff)
    assert "lak" in diff.connection_package_names
    assert "sfr" in diff.connection_package_names


def test_connection_difference_in_report(connection_tables):
    lak, _ = connection_tables
    lak["ref"] = _lak("ref", [(1, 0, 20, "HORIZONTAL", 5.0, 8.0, 50.0, 100.0)])
    lak["variant"] = _lak("variant", [(1, 0, 20, "HORIZONTAL", 5.0, 8.0, 75.0, 100.0)])
    report = _group(["ref", "variant"], reference="ref", package_names=["LAK"]).diff().report()
    assert "Connection differences" in report
    assert "differs from reference" in report


def test_connection_self_reference_errors(connection_tables):
    lak, _ = connection_tables
    lak["ref"] = _lak("ref", [(0, 0, 10, "VERTICAL", 0.0, 0.0, 0.0, 0.0)])
    lak["variant"] = _lak("variant", [(0, 0, 10, "VERTICAL", 0.0, 0.0, 0.0, 0.0)])
    diff = _group(["ref", "variant"], reference="ref").diff()
    with pytest.raises(ValueError):
        diff.packages.lak.inputs.connections(model_name="ref")
    with pytest.raises(KeyError):
        diff.packages.lak.inputs.connections(model_name="ghost")
