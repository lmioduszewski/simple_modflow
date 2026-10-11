"""One reader of a TDIS ``TIME_UNITS``, and every label made from it (ledger 442).

``prt._tdis_time_units`` (the PRT TDIS and the pathline hover), ``package_budget.
_declared_time_unit`` (budget labels) and ``hydraulics.unit_labels`` (the map hovers)
each read the declared unit, each with its own abbreviation table. They are one
function now, ``hydraulics.declared_time_unit``, and one table, ``MF6_TIME_UNITS`` --
MF6's six words, measured in ``docs/reproducers/tdis_time_units.py``.
"""

from __future__ import annotations

import ast
import logging
import pathlib
import types

import flopy
import pytest

from myflopy.modflow.mf6.hydraulics import MF6_TIME_UNITS, declared_time_unit, unit_labels
from myflopy.modflow.mf6.package_budget import budget_value_units
from myflopy.modflow.mf6.prt_maps import _time_unit

SRC = pathlib.Path(__file__).resolve().parents[1] / "src" / "myflopy"


def _model(time_units=None, *, set_it=True):
    """A FloPy GWF model in feet whose TDIS declares ``time_units`` (or nothing)."""

    sim = flopy.mf6.MFSimulation(sim_ws=".")
    flopy.mf6.ModflowTdis(sim, **({"time_units": time_units} if set_it else {}))
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
    flopy.mf6.ModflowGwfdis(gwf, nrow=1, ncol=2, length_units="feet")
    return types.SimpleNamespace(sim=sim, gwf=gwf, model_type="gwf6")


# (TDIS word or None for unset, the one reader's answer, the abbreviation)
CASES = [
    pytest.param(None, None, None, id="unset"),
    pytest.param("unknown", None, None, id="unknown"),
    pytest.param("UNKNOWN", None, None, id="UNKNOWN"),
    pytest.param("d", None, None, id="a-word-mf6-refuses"),
    pytest.param("seconds", "seconds", "s", id="seconds"),
    pytest.param("MINUTES", "minutes", "min", id="MINUTES"),
    pytest.param(" Hours ", "hours", "h", id="Hours-padded"),
    pytest.param("days", "days", "d", id="days"),
    pytest.param("Years", "years", "yr", id="Years"),
]


@pytest.mark.parametrize(("word", "unit", "label"), CASES)
def test_every_time_unit_label_reads_the_one_reader(word, unit, label):
    model = _model(word, set_it=word is not None)

    assert declared_time_unit(model) == unit
    # the map hovers (head, transmissivity, footer)
    assert unit_labels(model) == ("ft", label)
    # a budget value
    assert budget_value_units(model) == f"ft³/{label or 'T'}"
    # a pathline's travel time
    assert _time_unit(model) == (label or "")
    # a raw FloPy model, reached through gwf.simulation
    assert declared_time_unit(types.SimpleNamespace(gwf=model.gwf)) == unit


@pytest.mark.parametrize("blank", ["", "   "])
def test_a_blank_time_unit_is_no_unit(blank, caplog):
    """FloPy refuses a blank unit at construction, but an object can still hold one."""

    model = types.SimpleNamespace(sim=types.SimpleNamespace(tdis=types.SimpleNamespace(time_units=blank)))
    with caplog.at_level(logging.WARNING, logger="myflopy"):
        assert declared_time_unit(model) is None
    assert not caplog.records


def test_no_tdis_is_no_unit():
    assert declared_time_unit(types.SimpleNamespace()) is None
    assert declared_time_unit(types.SimpleNamespace(sim=None, gwf=None)) is None


def test_a_word_mf6_refuses_is_logged_by_name(caplog):
    with caplog.at_level(logging.WARNING, logger="myflopy"):
        assert declared_time_unit(_model("day")) is None
    assert "'day'" in caplog.text and "MF6 stops on it" in caplog.text


def test_the_table_is_mf6s_six_words():
    """The words MF6 reads (any case) -- every other one stops the run (reproducer)."""

    assert set(MF6_TIME_UNITS) == {"unknown", "seconds", "minutes", "hours", "days", "years"}
    assert MF6_TIME_UNITS["unknown"] is None


#: Readers kept out of the one rule ON PURPOSE. ``legacy_prt`` is the deprecated
#: tier: it copies FloPy's ``modeltime.time_units``, which answers ``"unknown"`` for a
#: TDIS naming no unit, and MF6 reads that (``docs/reproducers/tdis_time_units.py``).
_EXEMPT = {"modflow/mp3du/legacy_prt.py"}


def _time_unit_reads(path: pathlib.Path) -> list[str]:
    """``<x>.time_units`` / ``getattr(<x>, "time_units")`` where ``<x>`` names a TDIS."""

    source = path.read_text()
    found = []
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Attribute) and node.attr == "time_units":
            owner = ast.get_source_segment(source, node.value) or ""
        elif (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id in {"getattr", "_option", "_mfdata_value"}
            and len(node.args) >= 2
            and isinstance(node.args[1], ast.Constant)
            and node.args[1].value == "time_units"
        ):
            owner = ast.get_source_segment(source, node.args[0]) or ""
        else:
            continue
        if "tdis" in owner.lower() or "modeltime" in owner.lower():
            found.append(f"{path.relative_to(SRC).as_posix()}:{node.lineno}")
    return found


def test_only_the_one_reader_reads_a_tdis_time_unit():
    """A second reader of TDIS ``TIME_UNITS`` fails here, naming the line."""

    reads = [
        read
        for path in sorted(SRC.rglob("*.py"))
        if path.relative_to(SRC).as_posix() not in _EXEMPT
        for read in _time_unit_reads(path)
    ]
    assert [read.rsplit(":", 1)[0] for read in reads] == ["modflow/mf6/hydraulics.py"], reads


def test_the_exempt_legacy_reader_still_exists():
    """An exemption for a file that no longer reads the unit is a stale one."""

    for exempt in _EXEMPT:
        assert _time_unit_reads(SRC / exempt), exempt


def test_a_prt_tdis_declares_no_unit_for_a_flow_tdis_saying_unknown(tmp_path):
    """``unknown`` is MF6's undefined: the PRT TDIS mirrors it by naming none."""

    from myflopy.modflow.mf6.prt import PRTProject, PRTReleasePoints
    from tests.test_mf6_prt import _flow_model

    model = _flow_model("prt_unknown_units", tmp_path)
    model.sim.tdis.time_units = "unknown"
    project = PRTProject(
        model,
        workspace=tmp_path / "prt",
        release_points=PRTReleasePoints.from_cells(model, [0], local_z=0.5),
    )
    assert project.tdis.time_units.get_data() is None
    model.sim.tdis.time_units = "HOURS"
    hours = PRTProject(
        model,
        workspace=tmp_path / "prt_hours",
        release_points=PRTReleasePoints.from_cells(model, [0], local_z=0.5),
    )
    assert hours.tdis.time_units.get_data().upper() == "HOURS"
