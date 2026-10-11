"""A length unit is checked wherever it enters a model (ledgers 432, 433).

Two entries ledger 404 did not reach:

- a DIS/DISV/DISU spec's ``length_units`` given through ``with_options``,
  ``packages.fork`` or ``packages.edit``, which FloPy wrote as given -- ``ft`` as
  ``FT``, read by MF6 as unit 0 -- and the builder, which every route (a saved
  project, a frozen run) reaches;
- ``ModelContext(length_units=)``, never checked, and a ``Surface`` declaring
  ``units=`` sampled with NO model unit, which was read as model units in silence.
"""

from __future__ import annotations

import filecmp
import re
import shutil
import subprocess

import flopy
import numpy as np
import pytest
from shapely.geometry import box

import myflopy as mf
from myflopy._length_units import mf6_length_conversion
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.specs import ModelContext, PackageSpec, fork_package
from myflopy.surfaces import LayerSurfaces, Surface, sample_bed

UNIT_LINE = re.compile(r"Model length unit \[0=UND, 1=FEET, 2=METERS, 3=CENTIMETERS\] set as\s+(\d)")


@pytest.fixture(scope="module")
def vor():
    return irregular_voronoi_grid(nrow=2, ncol=2, cell_size=100.0)


def _dis(units=None):
    extra = {} if units is None else {"length_units": units}
    return mf.dis(nlay=1, nrow=1, ncol=2, delr=1.0, delc=1.0, top=10.0, botm=0.0, **extra)


def _disu(units=None):
    extra = {} if units is None else {"length_units": units}
    return mf.disu(nodes=1, nja=1, top=1.0, bot=0.0, area=1.0, iac=[1], ja=[0], **extra)


def _disv(units=None):
    extra = {} if units is None else {"length_units": units}
    return mf.disv(nlay=1, ncpl=2, nvert=6, top=10.0, botm=[0.0],
                   vertices=[[0, 0.0, 0.0], [1, 1.0, 0.0], [2, 2.0, 0.0],
                             [3, 0.0, 1.0], [4, 1.0, 1.0], [5, 2.0, 1.0]],
                   cell2d=[[0, 0.5, 0.5, 4, 0, 3, 4, 1], [1, 1.5, 0.5, 4, 1, 4, 5, 2]],
                   **extra)


# --- (a) the edit routes refuse what the helpers refuse --------------------- #


def _with_options(spec, units):
    return spec.with_options(length_units=units)


def _fork(spec, units):
    return fork_package(spec, length_units=units)


def _library(route):
    def edit(spec, units, tmp_path):
        project = mf.Project(tmp_path / "project", name="p")
        project.add_package("dis/base", spec)
        if route == "fork":
            return project.packages.fork("dis/base", "dis/variant", length_units=units)
        return project.packages.edit("dis/base", length_units=units)

    return edit


@pytest.mark.parametrize("make", [_dis, _disv, _disu], ids=["dis", "disv", "disu"])
@pytest.mark.parametrize("route", ["with_options", "fork_package", "packages.fork",
                                   "packages.edit"])
def test_every_edit_route_refuses_a_word_mf6_reads_as_undefined(tmp_path, make, route):
    spec = make()
    edit = {
        "with_options": lambda s, u: _with_options(s, u),
        "fork_package": lambda s, u: _fork(s, u),
        "packages.fork": lambda s, u: _library("fork")(s, u, tmp_path),
        "packages.edit": lambda s, u: _library("edit")(s, u, tmp_path),
    }[route]
    with pytest.raises(ValueError, match=r"length_units='ft'.*Write length_units='feet'"):
        edit(spec, "ft")
    with pytest.raises(ValueError, match=r"length_units='inches'.*not a word MF6 reads"):
        edit(spec, "inches")
    assert edit(spec, "METERS").options["length_units"] == "METERS"


def test_a_package_without_length_units_is_not_asked():
    """The check keys on the builder taking ``length_units``; an NPF edit is untouched."""

    assert mf.npf(k=1.0).with_options(k=2.0).options["k"] == 2.0


# --- (a) the builder writes MF6's word, whatever route the spec took -------- #


def _built_word(spec):
    sim = flopy.mf6.MFSimulation(sim_ws=".", exe_name="mf6")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
    return spec.build(gwf).length_units.get_data().lower()  # FloPy reads it back lower case


@pytest.mark.parametrize("make", [_dis, _disv, _disu], ids=["dis", "disv", "disu"])
def test_a_spec_holding_an_abbreviation_is_written_as_mf6s_word(make):
    """A saved project or a frozen run holding ``ft`` rebuilds declaring FEET.

    Built straight through ``PackageSpec`` and through ``from_dict(to_dict())`` --
    the two routes the declaration check does not see. The spec keeps what it held,
    so its digest does not move.
    """

    declared = make()
    held = PackageSpec(declared.name, declared.builder, {**declared.options, "length_units": "ft"})
    reloaded = PackageSpec.from_dict(held.to_dict())
    for spec in (held, reloaded):
        assert _built_word(spec) == "feet"
        assert spec.options["length_units"] == "ft"
    for word in ("feet", "Meters", "unknown"):
        assert _built_word(declared.with_options(length_units=word)) == word.lower()
    assert _built_word(PackageSpec(declared.name, declared.builder,
                                   {**declared.options, "length_units": "cm"})) == "centimeters"
    unknown = PackageSpec(declared.name, declared.builder,
                          {**declared.options, "length_units": "furlongs"})
    with pytest.raises(ValueError, match=r"length_units='furlongs' is not a length unit"):
        _built_word(unknown)


@pytest.mark.skipif(shutil.which("mf6") is None, reason="needs the mf6 executable")
def test_mf6_reads_a_rebuilt_ft_run_as_feet_and_computes_the_same_numbers(tmp_path):
    """Measured on MF6: the rebuilt ``ft`` run's listing says unit 1 (it said 0), and
    its outputs are byte-identical to the run under ``FT`` -- MF6 reads the unit for
    its labels only (``docs/reproducers/mf6_length_units.py``)."""

    declared = _disv()
    held = PackageSpec(declared.name, declared.builder, {**declared.options, "length_units": "ft"})

    def run(workspace, disv_spec, word=None):
        model = mf.gwf("m", packages=[disv_spec, mf.npf(k=1.0), mf.ic(strt=5.0),
                                      mf.chd(stress_period_data={0: [[(0, 0), 5.0]]}),
                                      mf.rch(stress_period_data={0: [[(0, 1), 0.01]]}),
                                      mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
                                            saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])])
        sim = mf.simulation(model, name="s",
                            packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
        sim.build_flopy(workspace).simulation.write_simulation(silent=True)
        if word is not None:  # what FloPy wrote for `ft` before: the word as given
            path = workspace / "m.disv"
            path.write_text(re.sub(r"LENGTH_UNITS\s+\S+", f"LENGTH_UNITS  {word}", path.read_text()))
        out = subprocess.run(["mf6"], cwd=workspace, capture_output=True, text=True, check=False)
        assert "Normal termination" in out.stdout, out.stdout[-1500:]
        return int(UNIT_LINE.search((workspace / "m.lst").read_text()).group(1))

    assert run(tmp_path / "rebuilt", held) == 1
    assert run(tmp_path / "before", held, word="FT") == 0
    for output in ("m.hds", "m.cbc"):
        assert filecmp.cmp(tmp_path / "rebuilt" / output, tmp_path / "before" / output,
                           shallow=False), output


# --- (b) ModelContext's unit is checked where it enters --------------------- #


def test_model_context_refuses_a_unit_nothing_reads():
    with pytest.raises(ValueError, match=r"ModelContext\(length_units='inches'\).*'ft'"):
        ModelContext(length_units="inches")
    with pytest.raises(ValueError, match=r"ModelContext\(length_units=3\.28"):
        ModelContext(length_units=3.28)
    for spelling in ("feet", " FT ", "metres", "us-ft", "Centimeter", None, "", "UNKNOWN"):
        assert ModelContext(length_units=spelling).length_units == spelling


def test_the_conversion_reads_the_one_table():
    """``"metres"``/``"us-ft"`` converted a surface but raised in ``mf.sfr``."""

    assert mf6_length_conversion("metres") == 1.0
    assert mf6_length_conversion("us-ft") == pytest.approx(3.28081)
    assert mf6_length_conversion(" Centimeter ") == 100.0
    for none in (None, "", "unknown"):
        assert mf6_length_conversion(none) == 1.0


# --- (b) a declared unit with no model unit is refused ---------------------- #


def test_every_sampler_refuses_a_declared_unit_with_no_model_unit(vor):
    metres = Surface.flat(30.0, units="meters")
    ground = np.full(vor.ncpl, 200.0)
    for no_unit in (None, "unknown"):
        with pytest.raises(ValueError, match=r"units='meters'.*length_units="):
            sample_bed(metres, vor, ground, length_units=no_unit)
    with pytest.raises(ValueError, match=r"units='meters'.*length_units="):
        LayerSurfaces([Surface.flat(100.0), metres]).sample(vor)
    stack = mf.LayerStack(top=100.0, length_units=None).add("a", bottom=metres)
    with pytest.raises(ValueError, match=r"units='meters'.*length_units="):
        stack.build(vor)
    # An operand inside algebra is asked too.
    with pytest.raises(ValueError, match=r"units='meters'"):
        sample_bed(metres.capped_at(1000.0), vor, ground)
    # Nothing declared: sampled as it is, as before.
    np.testing.assert_allclose(sample_bed(Surface.flat(30.0), vor, ground), 30.0)
    built = mf.LayerStack(top=100.0, length_units=None).add("a", bottom=Surface.flat(30.0))
    np.testing.assert_allclose(built.build(vor).botm[0], 30.0)


def test_mf_lak_refuses_a_metre_bed_on_a_context_with_no_unit(vor):
    lakes = __import__("geopandas").GeoDataFrame(
        {"name": ["pond"]}, geometry=[box(0.0, 0.0, 100.0, 100.0)], crs=vor.crs
    )
    bed = Surface.flat(3.0, units="meters")
    layers = mf.LayerStack(vor, top=20.0).add("a", thickness=10.0).build()
    for no_unit in (None, "unknown"):
        ctx = ModelContext(grid=vor, surfaces=layers, domain=layers.idomain,
                           length_units=no_unit)
        with pytest.raises(ValueError, match=r"units='meters'.*ModelContext\(length_units=\)"):
            mf.lak(context=ctx, nper=1, lakes=lakes, lake_id_field="name",
                   starting_stage=11.5, lake_bottom=bed, unexcavated="ignore")
