"""Layer-surface units: what the DISV says, what every operand is in, and what a refusal leaves.

Four defects of `surfaces.py`/`layers.py`, all silent:

- **LENGTH_UNITS** (ledger 404): a stack converts `"ft"`/`"m"`/`"cm"` but wrote
  them to the DISV upper-cased, `FT`/`M`/`CM`, which MF6 reads as length unit 0 --
  undefined -- under Normal termination. One helper now maps every spelling to MF6's
  word; a raw `mf.dis*(length_units=...)` must already BE MF6's word.
- **Operand units** (ledger 405): only the surface a stack was handed converted
  its `units=`; inside `capped_at`/`shift`/`where`/... a metre DEM was read as feet.
- **A refused `to_disv`** (ledger 407) had already written `vor.gdf_topbtm`,
  whose setter clears `vor.excavations`.
- **A repeated label** (ledger 406) kept the last surface's values under every
  position naming it.
"""

from __future__ import annotations

import re
import shutil
import subprocess

import numpy as np
import pytest
from shapely.geometry import box

import myflopy as mf
from myflopy._length_units import _LENGTH_TO_FEET, _MF6_LENGTH_UNITS, _mf6_length_units
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.surfaces import (
    LayerSurfaces,
    Surface,
    sample_bed,
)

M_TO_FT = 3.280839895013123
_WEST = box(-1000.0, -1000.0, 100.0, 1000.0)   # the 2x2 test grid spans 0..200


@pytest.fixture
def vor():
    return irregular_voronoi_grid(nrow=2, ncol=2, cell_size=100.0)


def _stack(length_units):
    return mf.LayerStack(top=Surface.flat(100.0), length_units=length_units).add(
        "a", thickness=10.0
    )


# --- LENGTH_UNITS: MF6's word for every spelling ---------------------------- #


def test_every_spelling_the_conversion_takes_has_the_mf6_word_of_its_own_factor():
    """The two tables are one rule: a spelling that converts must also be written."""

    assert set(_MF6_LENGTH_UNITS) == set(_LENGTH_TO_FEET)
    per_token = {"FEET": 1.0, "METERS": M_TO_FT, "CENTIMETERS": M_TO_FT / 100}
    for spelling, token in _MF6_LENGTH_UNITS.items():
        assert _LENGTH_TO_FEET[spelling] == pytest.approx(per_token[token]), spelling


@pytest.mark.parametrize(
    "units, word",
    [("ft", "FEET"), ("Foot", "FEET"), ("usft", "FEET"), ("feet", "FEET"),
     ("m", "METERS"), ("metres", "METERS"), ("cm", "CENTIMETERS")],
)
def test_every_to_disv_writes_mf6s_word(vor, units, word):
    ls = LayerSurfaces([Surface.flat(100.0), Surface.flat(90.0)])
    assert ls.to_disv(vor, length_units=units, attach=False).options["length_units"] == word
    assert _stack(units).to_disv(vor, attach=False).options["length_units"] == word
    assert _stack(units).build(vor).to_disv().options["length_units"] == word
    assert mf.disv(_stack(units).build(vor)).options["length_units"] == word


def test_an_unknown_length_unit_is_refused_naming_the_spellings(vor):
    with pytest.raises(ValueError, match=r"Unknown length unit 'furlongs'.*'ft'"):
        _mf6_length_units("furlongs")
    # At the declaration, so a lazy `mf.disv(stack)` cannot carry it to a build.
    with pytest.raises(ValueError, match="furlongs"):
        mf.LayerStack(top=100.0, length_units="furlongs")
    ls = LayerSurfaces([Surface.flat(100.0), Surface.flat(90.0)])
    with pytest.raises(ValueError, match="furlongs"):
        ls.to_disv(vor, length_units="furlongs")
    # `None` still declares none.
    assert "length_units" not in ls.to_disv(vor, attach=False).options
    assert "length_units" not in _stack(None).build(vor).to_disv().options


@pytest.mark.parametrize(
    "helper",
    [
        lambda units: mf.dis(nlay=1, nrow=1, ncol=1, delr=1.0, delc=1.0, top=1.0,
                             botm=0.0, length_units=units),
        lambda units: mf.disu(nodes=1, nja=1, top=1.0, bot=0.0, area=1.0, iac=[1],
                              ja=[0], length_units=units),
        lambda units: mf.disv(nlay=1, ncpl=1, nvert=4,
                              vertices=[[0, 0.0, 0.0], [1, 0.0, 1.0], [2, 1.0, 1.0],
                                        [3, 1.0, 0.0]],
                              cell2d=[[0, 0.5, 0.5, 4, 0, 1, 2, 3]], top=1.0, botm=[0.0],
                              length_units=units),
        lambda units: mf.disv(_stack("feet"), length_units=units),
    ],
    ids=["dis", "disu", "disv", "disv-lazy"],
)
def test_a_raw_length_units_option_must_be_a_word_mf6_reads(helper):
    """A FloPy keyword is written as given: 'ft' is refused naming 'feet', not mapped."""

    for good in ("feet", "FEET", "Meters", "centimeters", "unknown"):
        assert helper(good).options["length_units"] == good
    with pytest.raises(ValueError, match=r"length_units='ft'.*Write length_units='feet'"):
        helper("ft")
    with pytest.raises(ValueError, match=r"length_units='M'.*Write length_units='meters'"):
        helper("M")
    with pytest.raises(ValueError, match=r"not a word MF6 reads"):
        helper("furlongs")


@pytest.mark.skipif(shutil.which("mf6") is None, reason="needs the mf6 executable")
@pytest.mark.parametrize("units, unit_number", [("ft", 1), ("m", 2), ("cm", 3)])
def test_mf6_reads_the_unit_a_stack_declared_by_its_abbreviation(
    tmp_path, vor, units, unit_number
):
    """Measured on MF6 itself: its listing names the unit; `FT` read as 0 (undefined)."""

    model = mf.gwf(
        "m",
        packages=[
            _stack(units).to_disv(vor, attach=False),
            mf.npf(k=1.0), mf.ic(strt=95.0),
            mf.chd(stress_period_data={0: [[(0, 0), 95.0]]}),
        ],
    )
    sim = mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    built = sim.build_flopy(tmp_path / units)
    built.simulation.write_simulation(silent=True)
    run = subprocess.run(["mf6"], cwd=tmp_path / units, capture_output=True, text=True, check=False)
    assert "Normal termination" in run.stdout, run.stdout[-2000:]
    listing = (tmp_path / units / "m.lst").read_text()
    found = re.search(r"Model length unit \[0=UND, 1=FEET, 2=METERS, 3=CENTIMETERS\] "
                      r"set as\s+(\d)", listing)
    assert found is not None and int(found.group(1)) == unit_number


# --- units= on every operand ------------------------------------------------- #


def _meters(value):
    return Surface.flat(value, units="meters")


@pytest.mark.parametrize(
    "surface, expected_m",
    [
        (Surface.minimum(_meters(30.0), 1000.0), 30.0 * M_TO_FT),
        (Surface.maximum(_meters(30.0), 0.0), 30.0 * M_TO_FT),
        (_meters(30.0).capped_at(1000.0), 30.0 * M_TO_FT),
        (_meters(30.0).floored_at(_meters(10.0)), 30.0 * M_TO_FT),
        (Surface.clamp(_meters(30.0), lower=_meters(40.0)), 40.0 * M_TO_FT),
        (_meters(30.0).between(upper=_meters(20.0)), 20.0 * M_TO_FT),
        (_meters(30.0).below(5.0), 30.0 * M_TO_FT - 5.0),      # distance in model units
        (_meters(30.0) + 5.0, 30.0 * M_TO_FT + 5.0),
        (Surface.where(_WEST, _meters(1.0), _meters(2.0)), None),
        (Surface.isopach(_meters(10.0)), 200.0 - 10.0 * M_TO_FT),  # a thickness in m
        (Surface.toward(_meters(10.0), 0.5), 100.0 + 5.0 * M_TO_FT),  # halfway from 200 ft
        (Surface.minimum(Surface.flat(30.0, units="feet"), _meters(5.0)), 5.0 * M_TO_FT),
    ],
    ids=["minimum", "maximum", "capped_at", "floored_at", "clamp", "between", "below",
         "add", "where", "isopach", "toward", "mixed"],
)
def test_units_on_every_operand_are_converted_before_they_combine(vor, surface, expected_m):
    """Every sampling path: the bed sampler, a LayerSurfaces column and a stack's top."""

    if surface.kind == "where":
        west = np.asarray(vor.centroids[0]) < 100.0
        assert 0 < west.sum() < vor.ncpl
        expected = np.where(west, 1.0 * M_TO_FT, 2.0 * M_TO_FT)
    else:
        expected = np.full(vor.ncpl, expected_m)
    ground = np.full(vor.ncpl, 200.0)
    np.testing.assert_allclose(sample_bed(surface, vor, ground, length_units="feet"), expected)
    ls = LayerSurfaces([Surface.flat(200.0), surface])
    sampled = ls.sample(vor, reconcile=False, length_units="ft")
    np.testing.assert_allclose(sampled[1].to_numpy(dtype=float), expected)
    if not surface.is_relative and surface.kind not in ("isopach", "toward"):
        result = mf.LayerStack(top=surface, length_units="feet").add(
            "a", bottom=-1000.0).build(vor, reconcile=False)
        np.testing.assert_allclose(result.top, expected)


def test_a_converted_raster_fills_from_the_surface_above_as_it_is(tmp_path, vor):
    """`fill="propagate"` in metres: the inherited cells are already in feet."""

    rasterio = pytest.importorskip("rasterio")
    from rasterio.transform import from_origin

    xmin, ymin, xmax, ymax = vor.gdf_vorPolys.total_bounds
    pixels = np.full((2, 2), 10.0, dtype="float32")
    pixels[:, 1] = -9999.0                                  # the east half has no data
    path = tmp_path / "half.tif"
    with rasterio.open(path, "w", driver="GTiff", height=2, width=2, count=1,
                       dtype="float32", crs=vor.crs, nodata=-9999.0,
                       transform=from_origin(xmin, ymax, (xmax - xmin) / 2,
                                             (ymax - ymin) / 2)) as dst:
        dst.write(pixels, 1)
    east = vor.centroids[0] > (xmin + xmax) / 2
    ls = LayerSurfaces([Surface.flat(100.0),
                        Surface.raster(path, units="meters", fill="propagate")])
    bottom = ls.sample(vor, reconcile=False, length_units="feet")[1].to_numpy(dtype=float)
    np.testing.assert_allclose(bottom[east], 100.0)          # was 328.08: 100 ft x 3.28
    np.testing.assert_allclose(bottom[~east], 10.0 * M_TO_FT)


def test_no_model_unit_refuses_a_declared_unit_and_units_on_a_non_leaf_are_refused(vor):
    ground = np.full(vor.ncpl, 200.0)
    # With no model unit there is nothing to convert into: it read 30.0 m as 30.0 in
    # model units, silently (ledger 433).
    with pytest.raises(ValueError, match=r"units='meters'.*length_units="):
        sample_bed(_meters(30.0).capped_at(1000.0), vor, ground)
    import dataclasses

    on_algebra = dataclasses.replace(Surface.maximum(_meters(30.0), 0.0), units="meters")
    with pytest.raises(ValueError, match=r"'max' surface.*operands"):
        sample_bed(on_algebra, vor, ground, length_units="feet")
    on_relative = dataclasses.replace(Surface.offset_below(5.0), units="meters")
    with pytest.raises(ValueError, match="absolute surfaces"):
        sample_bed(_meters(30.0).capped_at(on_relative), vor, ground, length_units="feet")
    with pytest.raises(ValueError, match="Unknown length unit 'parsecs'"):
        sample_bed(Surface.maximum(Surface.flat(1.0, units="parsecs"), 0.0), vor, ground,
                   length_units="feet")


def test_a_surface_digest_does_not_move_with_the_conversion():
    """What moved is the sampled values; the declaration is the same payload."""

    surface = Surface.raster("dem.tif", units="meters").capped_at(1000.0)
    assert surface.to_dict() == {
        "kind": "min",
        "operands": [{"kind": "raster", "path": "dem.tif", "units": "meters"},
                     {"kind": "flat", "value": 1000.0}],
    }


# --- a refused to_disv leaves the grid alone --------------------------------- #


@pytest.mark.parametrize(
    "kwargs",
    [
        {"pinch_out": True, "idomain": "ones"},
        {"length_units": "furlongs"},
        {"pinch_out": True, "pinch": "evaporate"},
        {"pinch_out": True, "min_sep": 5.0, "minimum_thickness": 1.0},
    ],
    ids=["idomain-and-pinch-out", "unknown-unit", "unknown-pinch", "min-sep"],
)
def test_a_refused_to_disv_has_not_attached_anything(vor, kwargs):
    if kwargs.get("idomain") == "ones":
        kwargs = {**kwargs, "idomain": np.ones((1, vor.ncpl), dtype=int)}
    sentinel = object()
    vor._gdf_topbtm = None
    vor._excavations = {"pond": sentinel}
    ls = LayerSurfaces([Surface.flat(100.0), Surface.flat(90.0)])
    with pytest.raises(ValueError):
        ls.to_disv(vor, **kwargs)
    assert vor._gdf_topbtm is None
    assert vor._excavations == {"pond": sentinel}
    # A call that succeeds still attaches, and drops the records with the old table.
    ls.to_disv(vor)
    assert vor._gdf_topbtm is not None and not vor.excavations


# --- one column per surface --------------------------------------------------- #


def test_a_repeated_label_is_refused_naming_it(vor):
    flats = [Surface.flat(100.0), Surface.flat(90.0), Surface.flat(80.0)]
    with pytest.raises(ValueError, match=r"'a' at positions \[1, 2\]"):
        LayerSurfaces(flats, labels=["t", "a", "a"])
    with pytest.raises(ValueError, match=r"0 at positions \[0, 2\]"):
        LayerSurfaces(flats, labels=[0, 1, 0])
    sampled = LayerSurfaces(flats, labels=["t", "a", "b"]).sample(vor, reconcile=False)
    assert sampled.iloc[0][["t", "a", "b"]].tolist() == [100.0, 90.0, 80.0]


def test_mf6s_unknown_is_refused_in_a_stack_naming_none():
    """MF6 KNOWS `UNKNOWN`, so the generic "a word MF6 does not know" message was false.

    The stack's spelling of "no unit" is `None` (no LENGTH_UNITS: MF6 reads it as
    undefined, as it reads UNKNOWN), and the refusal names it (ledger 404).
    """

    for spelling in ("unknown", "UNKNOWN"):
        with pytest.raises(ValueError, match=r"length_units=None") as raised:
            mf.LayerStack(top=100.0, length_units=spelling)
        assert "does not know" not in str(raised.value)
        with pytest.raises(ValueError, match=r"length_units=None"):
            _mf6_length_units(spelling)
