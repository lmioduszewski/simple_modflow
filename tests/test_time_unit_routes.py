"""Time and length units each have one table, checked where they enter (ledger 456).

- ``mf.tdis(time_units=)`` and every edit of a TDIS spec took any word, FloPy wrote
  it, and MF6 STOPPED on anything but its six words
  (``docs/reproducers/tdis_time_units.py``);
- ``ModelContext(time_units=)`` and ``LayerStack(time_units=)`` took any string,
  refused only where SFR/LAK derived a ``TIME_CONVERSION``;
- the build warning's flux label and the map hovers abbreviated through two other
  tables, and the PRT model copied its flow DIS's ``length_units`` word as given.
"""

from __future__ import annotations

import ast
import shutil
import subprocess
from pathlib import Path

import flopy
import pytest

import myflopy as mf
from myflopy._length_units import (
    _MF6_LENGTH_UNITS,
    length_abbreviation,
    mf6_length_abbreviation,
    mf6_length_mirror,
)
from myflopy._time_units import (
    _TIME_SPELLINGS,
    MF6_TIME_UNITS,
    check_tdis_time_units,
    mf6_time_conversion,
    time_abbreviation,
)
from myflopy.grid_identity import _flux_units
from myflopy.specs import ModelContext, fork_package

SRC = Path(__file__).resolve().parents[1] / "src" / "myflopy"


# --- a TDIS takes MF6's six words, at every route ---------------------------- #


@pytest.mark.parametrize("word", ["d", "day", "Days ", "hr", "undefined", "weeks", 1])
def test_mf_tdis_refuses_a_word_mf6_stops_on(word):
    with pytest.raises(ValueError, match="not a word MF6 reads"):
        mf.tdis(time_units=word)


def test_mf_tdis_names_mf6s_word_for_a_short_spelling():
    with pytest.raises(ValueError, match=r"time_units='d'.*Write time_units='days'"):
        mf.tdis(time_units="d")
    with pytest.raises(ValueError, match=r"Write time_units='hours'"):
        mf.tdis(time_units="hr")


@pytest.mark.parametrize("word", [None, *MF6_TIME_UNITS, "DAYS", "Hours"])
def test_mf_tdis_takes_mf6s_words_in_any_case_as_given(word):
    spec = mf.tdis(time_units=word)
    assert spec.options.get("time_units") == word


def _library(route, tmp_path):
    def edit(spec, units):
        project = mf.Project(tmp_path / "project", name="p")
        project.add_package("tdis/base", spec)
        if route == "packages.fork":
            return project.packages.fork("tdis/base", "tdis/variant", time_units=units)
        return project.packages.edit("tdis/base", time_units=units)

    return edit


@pytest.mark.parametrize("route", ["with_options", "fork_package", "packages.fork",
                                   "packages.edit"])
def test_every_tdis_edit_route_refuses_a_word_mf6_stops_on(tmp_path, route):
    spec = mf.tdis(time_units="days")
    edit = {
        "with_options": lambda s, u: s.with_options(time_units=u),
        "fork_package": lambda s, u: fork_package(s, time_units=u),
        "packages.fork": _library("packages.fork", tmp_path),
        "packages.edit": _library("packages.edit", tmp_path),
    }[route]
    with pytest.raises(ValueError, match=r"time_units='d'.*Write time_units='days'"):
        edit(spec, "d")
    assert edit(spec, "HOURS").options["time_units"] == "HOURS"


@pytest.mark.skipif(shutil.which("mf6") is None, reason="needs mf6 on PATH")
@pytest.mark.parametrize("word", ["DAYS", "unknown", None])
def test_mf6_runs_every_tdis_time_unit_the_declaration_accepts(tmp_path, word):
    """What ``mf.tdis`` lets through, MF6 reads; ``d`` -- refused above -- it stops on."""

    sim = flopy.mf6.MFSimulation(sim_ws=str(tmp_path), exe_name="mf6", verbosity_level=0)
    mf.tdis(time_units=word).build(sim)
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
    flopy.mf6.ModflowGwfdis(gwf, nrow=1, ncol=2, length_units="feet")
    flopy.mf6.ModflowGwfnpf(gwf)
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[((0, 0, 0), 1.0)])
    sim.write_simulation(silent=True)
    subprocess.run(["mf6"], cwd=tmp_path, capture_output=True, check=False)
    assert "Normal termination" in (tmp_path / "mfsim.lst").read_text()


# --- a ModelContext / LayerStack takes any spelling of a unit --------------- #


@pytest.mark.parametrize("units", ["weeks", "fortnight", "daily", 86400])
def test_model_context_refuses_a_time_unit_nothing_reads(units):
    with pytest.raises(ValueError, match=r"ModelContext\(time_units="):
        ModelContext(time_units=units)


@pytest.mark.parametrize("units", [None, "", "unknown", " Days ", "d", "HR", "years"])
def test_model_context_takes_every_spelling_its_consumers_read(units):
    assert ModelContext(time_units=units).time_units == units


def test_layer_stack_refuses_a_time_unit_nothing_reads():
    with pytest.raises(ValueError, match=r"^LayerStack\(time_units='weeks'\) is not a time unit"):
        mf.LayerStack(top=10.0, time_units="weeks")
    assert mf.LayerStack(top=10.0, time_units="d").time_units == "d"


def test_every_spelling_converts_as_its_mf6_word_does():
    for spelling, word in _TIME_SPELLINGS.items():
        assert mf6_time_conversion(spelling.upper()) == mf6_time_conversion(word)
    assert mf6_time_conversion("d") == 86400.0
    assert mf6_time_conversion(None) == mf6_time_conversion("unknown") == 1.0


# --- one table per dimension, and every label agrees ------------------------- #


@pytest.mark.parametrize("length", sorted(_MF6_LENGTH_UNITS))
@pytest.mark.parametrize("time", sorted(_TIME_SPELLINGS))
def test_the_build_warning_flux_label_reads_the_one_table(length, time):
    """``metres``/``us-ft``/``centimeter`` were spelled out in full (``metres^3/days``)."""

    label = _flux_units(ModelContext(length_units=f" {length.upper()} ", time_units=time))
    word = _MF6_LENGTH_UNITS[length]
    assert label == f"{mf6_length_abbreviation(word)}3/{MF6_TIME_UNITS[_TIME_SPELLINGS[time]]}"


def test_a_context_declaring_no_unit_says_so():
    assert _flux_units(ModelContext(length_units=None, time_units="unknown")) == "length3/time"
    assert _flux_units(ModelContext(length_units="feet", time_units=None)) == "ft3/time"


def test_a_dis_word_mf6_reads_as_undefined_labels_no_unit():
    """``FT`` is unit 0 to MF6, so a label read back from such a DIS names none."""

    assert mf6_length_abbreviation("FEET") == mf6_length_abbreviation("feet") == "ft"
    assert mf6_length_abbreviation("FT") is None
    assert mf6_length_abbreviation("UNKNOWN") is None
    assert length_abbreviation("FT") == "ft"
    assert time_abbreviation("hr") == "h" and time_abbreviation("weeks") is None


_UNIT_WORDS = {"seconds", "minutes", "hours", "days", "years", "feet", "meters",
               "metres", "centimeters"}
_TABLE_HOMES = {"_length_units.py", "_time_units.py"}


def test_no_module_but_the_two_holds_a_unit_spelling_table():
    """A dict literal keyed by two or more unit names outside the two homes fails here.

    ``grid_identity._UNIT_ABBREVIATIONS``, ``hydraulics._LENGTH_ABBREVIATIONS`` and
    ``specs._MF6_TIME_TO_MODEL`` were three such tables, disagreeing about
    ``metres``, ``us-ft`` and spaces.
    """

    found = []
    for path in sorted(SRC.rglob("*.py")):
        if path.name in _TABLE_HOMES or "_vendor" in path.parts:
            continue
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.Dict):
                keys = {k.value.lower() for k in node.keys
                        if isinstance(k, ast.Constant) and isinstance(k.value, str)}
                if len(keys & _UNIT_WORDS) >= 2:
                    found.append(f"{path.relative_to(SRC)}:{node.lineno}")
    assert found == []


# --- PRT mirrors its flow model's units ------------------------------------- #


def test_check_tdis_time_units_is_what_prt_project_asks(tmp_path):
    from myflopy.modflow.mf6.prt import PRTProject, PRTReleasePoints
    from tests.test_mf6_prt import _flow_model

    model = _flow_model("prt_time_units", tmp_path)
    with pytest.raises(ValueError, match=r"PRTProject\(time_units='d'\).*'days'"):
        PRTProject(
            model,
            workspace=tmp_path / "prt",
            release_points=PRTReleasePoints.from_cells(model, [0], local_z=0.5),
            time_units="d",
        )
    assert not (tmp_path / "prt").exists()
    assert check_tdis_time_units("YEARS", helper="x") == "YEARS"


def test_the_legacy_temporal_discretization_refuses_a_word_mf6_stops_on(tmp_path):
    """The legacy OO ``TemporalDiscretization`` wrote ``"d"`` and MF6 stopped on it
    (review of ledger 456); it is the same TDIS, so it keeps ``mf.tdis``'s rule."""

    from myflopy.modflow.mf6.simulation.discretization import TemporalDiscretization
    from tests.test_mf6_prt import _flow_model

    model = _flow_model("legacy_tdis_units", tmp_path)
    with pytest.raises(
        ValueError, match=r"TemporalDiscretization\(time_units='d'\).*Write time_units='days'"
    ):
        TemporalDiscretization(model=model, time_units="d", per_len=1.0, num_steps=1)
    assert model.sim.tdis.time_units.get_data().upper() == "DAYS"   # left as it was
    TemporalDiscretization(model=model, time_units="Hours", per_len=1.0, num_steps=1)
    assert model.sim.tdis.time_units.get_data().lower() == "hours"


@pytest.mark.parametrize(("flow_word", "prt_word"),
                         [("ft", "feet"), ("METERS", "meters"), ("cubits", "cubits")])
def test_a_prt_dis_writes_mf6s_word_for_its_flow_dis(tmp_path, flow_word, prt_word):
    """A loaded run whose DIS says ``ft`` gives PRT ``feet``; a word MF6 read as unit 0 stays."""

    from myflopy.modflow.mf6.prt import PRTProject, PRTReleasePoints
    from tests.test_mf6_prt import _flow_model

    model = _flow_model(f"prt_len_{flow_word}", tmp_path)
    model.gwf.disv.length_units = flow_word
    project = PRTProject(
        model,
        workspace=tmp_path / "prt",
        release_points=PRTReleasePoints.from_cells(model, [0], local_z=0.5),
    )
    assert str(project.discretization.length_units.get_data()).lower() == prt_word


def test_the_mirror_maps_only_what_it_knows():
    assert mf6_length_mirror(None) is None
    assert mf6_length_mirror("feet") == "feet"
    assert mf6_length_mirror(" Ft ") == "FEET"
    assert mf6_length_mirror("cubits") == "cubits"
