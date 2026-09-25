"""The lake-family docstrings: every parameter documented, and what they claim, pinned.

``mf.lak`` documented 13 of its 29 parameters and ``LakeOutlet`` called its rate an
"outflow rate" -- a positive one ADDS water. The rewrite (ledger 318) gave every
parameter its own ``name : type`` entry in signature order and measured each claim;
this file keeps both from rotting -- the entries' existence, the corrected TEXT (the
outlet sign, what MF6 does with a VERTICAL and a HORIZONTAL connection's elevations),
and ``LAKBuilder``'s copy of ``mf.lak``'s entries. It reads the docstrings with a
small parser rather than numpydoc, which is not a dependency here.

The MF6 claims run the model (``slow``). What the builders accept without checking is
pinned too, under ledger 318: those tests describe CURRENT behaviour the docstrings now
state, so a fix must update the docstring and the test together.
"""

from __future__ import annotations

import dataclasses
import inspect
import re
import shutil
import subprocess
import warnings
from pathlib import Path

import flopy
import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import box

import myflopy as mf
from myflopy import (
    LAKBuilder,
    LakeConnection,
    LakeOutlet,
    LakeTable,
    LakeTableBuilder,
    ModelContext,
)
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus

ENTRY = re.compile(r"^(\*{0,2}[A-Za-z_]\w*) : (\S.*)$")


def _parameters(doc: str) -> list[tuple[str, str]]:
    """``(name, type)`` for every entry of a numpydoc Parameters section."""

    lines = inspect.cleandoc(doc).splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip() == "Parameters")
    assert set(lines[start + 1].strip()) == {"-"}, "Parameters has no underline"
    entries = []
    body = lines[start + 2:]
    for index, line in enumerate(body):
        following = body[index + 1].strip() if index + 1 < len(body) else ""
        if following and set(following) == {"-"}:
            break  # this line is the NEXT section's title
        if line.startswith(" ") or not line.strip():
            continue
        match = ENTRY.match(line)
        assert match, f"a Parameters line is not 'name : type': {line!r}"
        entries.append((match.group(1), match.group(2)))
    return entries


def _signature_names(function) -> list[str]:
    names = []
    for name, parameter in inspect.signature(function).parameters.items():
        if name == "self":
            continue
        names.append("**" + name if parameter.kind is parameter.VAR_KEYWORD else name)
    return names


@pytest.mark.parametrize(
    "subject, expected",
    [
        (LakeConnection, None),
        (LakeOutlet, None),
        (LakeTable, None),
        (LakeTableBuilder, None),
        (LAKBuilder, None),
        (type(mf.lak).__call__, "signature"),
        (type(mf.lak).flopy, "signature"),
        (mf.lak_connection, "signature"),
    ],
    ids=["LakeConnection", "LakeOutlet", "LakeTable", "LakeTableBuilder", "LAKBuilder", "mf.lak",
         "mf.lak.flopy", "mf.lak_connection"],
)
def test_every_parameter_has_its_own_typed_entry_in_signature_order(subject, expected):
    if expected == "signature":
        wanted = _signature_names(subject)
    else:
        wanted = [field.name for field in dataclasses.fields(subject) if field.init]
    documented = _parameters(subject.__doc__)
    assert [name for name, _ in documented] == wanted
    assert all(kind.strip() for _, kind in documented)


def test_the_text_an_editor_shows_for_mf_lak_points_at_the_call():
    """Signature help on ``mf.lak(`` shows the helper CLASS's docstring (jedi, measured),
    which holds none of the call's parameters -- so it must say where they are."""

    doc = type(mf.lak).__doc__
    assert "help(mf.lak)" in doc and "mf.lak.__call__" in doc
    assert "IGNORES the bed elevation" in doc
    # ...and that the bed still reaches MF6 through the HORIZONTAL faces (review of D0:
    # "the bed only picks the layer" was false for "bathy").
    assert "sets the bottom of every HORIZONTAL face written on the lake cell" in _flat(doc)
    assert "exchanges water through those faces from the bed up" in _flat(doc)
    assert "in a TRANSIENT period it lets no more out of such a lake than arrives" in _flat(doc)


def _flat(text: str) -> str:
    """``text`` with its line breaks and indentation collapsed to single spaces."""

    return " ".join(text.split())


def _entry(doc: str, name: str) -> str:
    """The description of one Parameters entry, whitespace collapsed."""

    lines = inspect.cleandoc(doc).splitlines()
    start = next(i for i, line in enumerate(lines) if ENTRY.match(line) and ENTRY.match(line).group(1) == name)
    body = []
    for line in lines[start + 1:]:
        if line and not line.startswith(" "):
            break
        body.append(line)
    return _flat("\n".join(body))


def test_the_docstrings_state_the_outlet_sign_and_the_vertical_rule():
    """The corrections themselves, as TEXT: revert one and this names it.

    The parity test only sees that an entry exists, and the MF6 tests pin MF6's
    physics, not what a docstring says about it.
    """

    call = type(mf.lak).__call__.__doc__
    rate = _entry(LakeOutlet.__doc__, "rate")
    assert "A NEGATIVE rate removes water from the lake and a POSITIVE rate ADDS water" in rate
    assert "outflow rate" not in _flat(LakeOutlet.__doc__).lower()
    outlets = _entry(call, "outlets")
    assert "a NEGATIVE rate removes water from the lake and a POSITIVE rate ADDS water" in outlets

    modes = _entry(call, "connection_modes")
    assert "MF6 IGNORES the BELEV and TELEV written for a VERTICAL connection" in modes
    assert "does not lower the VERTICAL connection" in modes
    assert "writes the lake cell's own HORIZONTAL faces from the bed up" in modes
    assert "MF6 exchanges water from the bed, below the cell top" in modes
    assert "in a TRANSIENT period MF6 caps the water leaving a lake" in modes
    assert "only picks the layer" not in _flat(call) + _flat(type(mf.lak).__doc__)
    bottom = _entry(call, "lake_bottom")
    assert "through which MF6 exchanges water from the bed up" in bottom


def test_lakbuilder_repeats_mf_lak_entry_for_entry():
    """``LAKBuilder`` is public (``from myflopy import LAKBuilder``) and ``mf.lak`` is it
    plus a ``lak_index``, so its Parameters are ``mf.lak``'s TEXT. Four entries differ
    on purpose: ``tables`` (``None`` raises here), ``name``/``mover`` (the endpoint is
    ``builder.connection``) and ``options`` (a mapping read at ``build``). If this
    fails, copy ``mf.lak``'s entry into ``lakes.py``."""

    call = type(mf.lak).__call__.__doc__
    own = {"tables", "name", "mover", "options"}
    for name in [field.name for field in dataclasses.fields(LAKBuilder) if field.init]:
        if name in own:
            continue
        assert _entry(LAKBuilder.__doc__, name) == _entry(call, name), name
    tail = "Stage / volume / surface-area tables"
    ours, theirs = _entry(LAKBuilder.__doc__, "tables"), _entry(call, "tables")
    assert ours[ours.index(tail):] == theirs[theirs.index(tail):]
    for paragraph in ("**Per-lake values.**", "**Per-period values.**"):
        mine, facade = _flat(LAKBuilder.__doc__), _flat(call)
        end = "(measured, all five)." if "period" in paragraph else "(``LAKBuilder.lake_ids``)."
        assert mine[mine.index(paragraph):mine.index(end, mine.index(paragraph))] == \
            facade[facade.index(paragraph):facade.index(end, facade.index(paragraph))]


def test_mf_lak_flopy_states_what_it_requires_and_raises():
    doc = type(mf.lak).flopy.__doc__
    for block in ("packagedata", "connectiondata", "perioddata"):
        assert _entry(doc, block).startswith("Required."), block
    assert "Raises" in doc and "TypeError" in doc
    assert "must equal the number of ``outlets`` rows" in _entry(doc, "noutlets")
    assert "must equal the number of ``tables`` rows" in _entry(doc, "ntables")
    assert "REPLACES" in _entry(doc, "nlakes")
    records = dict(packagedata=[[0, 101.0, 1]], perioddata={0: []},
                   connectiondata=[[0, 0, (0, 42), "VERTICAL", 0.1, 0.0, 0.0, 0.0, 0.0]])
    assert mf.lak.flopy(**records, nlakes=3).options["nlakes"] == 3, "written as given"
    with pytest.raises(TypeError, match="did you mean 'print_stage'"):
        mf.lak.flopy(**records, print_stagee=True)
    with pytest.raises(TypeError, match="packagedata"):
        mf.lak.flopy(connectiondata=records["connectiondata"], perioddata={})


# --------------------------------------------------------------------------- toy lake


def _toy_grid() -> VoronoiGridPlus:
    """4 x 4 unit cells, two layers (tops 12, 8, 0)."""

    n = 4
    verts = np.array([[x, y] for y in range(n + 1) for x in range(n + 1)], float)
    iverts, centers = [], []
    for r in range(n):
        for c in range(n):
            v0 = r * (n + 1) + c
            iverts.append([v0, v0 + 1, v0 + n + 2, v0 + n + 1])
            centers.append([c + 0.5, r + 0.5])
    grid = VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=np.array(centers))
    grid.gdf_topbtm = gpd.GeoDataFrame(
        {0: [12.0] * 16, 1: [8.0] * 16, 2: [0.0] * 16}, geometry=grid.gdf_vorPolys.geometry, crs=grid.crs
    )
    grid.get_disu_connectivity(validate=True)
    return grid


LAKE_CELLS = {5, 6, 9, 10}


@pytest.fixture(scope="module")
def toy():
    grid = _toy_grid()
    lakes = gpd.GeoDataFrame({"name": ["pond"], "num": [101]}, geometry=[box(1, 1, 3, 3)], crs=grid.crs)
    ctx = ModelContext(grid=grid, domain=np.ones((2, 16), dtype=int))
    return {"context": ctx, "lakes": lakes}


def _connections(toy, **kw):
    base = dict(context=toy["context"], nper=1, lakes=toy["lakes"], lake_id_field="name",
                starting_stage=11.5, bed_leakance=0.1)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        return LAKBuilder(**{**base, **kw}).connectiondata


def test_bathy_writes_what_connection_modes_says(toy):
    bed = pd.Series({cell: 12.0 for cell in range(16)})
    bed[5], bed[6], bed[9], bed[10] = 9.0, 10.0, 10.0, 11.0
    rows = _connections(toy, lake_bottom=bed)

    vertical = [row for row in rows if row[3] == "VERTICAL"]
    assert [(row[2], row[5], row[6]) for row in vertical] == [
        ((0, 5), 9.0, 9.0), ((0, 6), 10.0, 10.0), ((0, 9), 10.0, 10.0), ((0, 10), 11.0, 11.0)
    ], "one VERTICAL per lake cell, belev = telev = the bed"
    # Every HORIZONTAL is written on the lake cell ITSELF, never on the neighbour.
    assert {row[2][1] for row in rows} == LAKE_CELLS
    # Cell 5 (bed 9) steps up to two land cells (12) and two lake cells (10).
    assert sorted(row[6] for row in rows if row[2] == (0, 5) and row[3] == "HORIZONTAL") == [
        10.0, 10.0, 12.0, 12.0
    ]
    # Cell 10, the highest bed (11), has no lower neighbour inside the lake.
    assert sum(1 for row in rows if row[2] == (0, 10) and row[3] == "HORIZONTAL") == 2
    # CONNLEN the lake cell's OWN side of the face (0.5 on a unit square), not the 1.0
    # between the centres that grid.cl12 held until ledger 328; CONNWIDTH the face.
    assert {(row[7], row[8]) for row in rows if row[3] == "HORIZONTAL"} == {(0.5, 1.0)}


def test_rectangular_writes_what_connection_modes_says(toy):
    kw = dict(connection_modes="rectangular", lake_bottom=9.0, lake_top=11.0)
    rows = _connections(toy, **kw)
    assert {(row[3], row[5], row[6]) for row in rows} == {
        ("VERTICAL", 9.0, 11.0), ("HORIZONTAL", 9.0, 11.0)
    }
    assert sum(row[3] == "HORIZONTAL" for row in rows) == 8, "edges only: 2 outside faces per cell"
    assert {row[2][1] for row in rows} == LAKE_CELLS
    assert sum(row[3] == "HORIZONTAL" for row in _connections(toy, rectangular_interior=True, **kw)) == 16
    assert all(row[3] == "VERTICAL" for row in _connections(toy, only_vertical=True, **kw))


def test_a_rectangular_lake_drops_its_vertical_connection_silently(toy):
    """Ledger 318: documented, not fixed. ``"bathy"`` raises in the same place."""

    rows = _connections(toy, connection_modes="rectangular", lake_bottom=6.0, lake_top=11.0,
                        only_layer=0)
    assert rows and not any(row[3] == "VERTICAL" for row in rows)
    with pytest.raises(ValueError, match="does not intersect active cell"):
        _connections(toy, lake_bottom=6.0, only_layer=0)


def _period(toy, key="name", **kw):
    ids = {"name": "pond", "num": "101"}
    base = dict(context=toy["context"], nper=3, lakes=toy["lakes"], lake_id_field=key,
                starting_stage=11.5, lake_bottom=9.0, connection_modes="rectangular",
                lake_top=11.0, bed_leakance=0.1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return mf.lak(**base, **kw).options["perioddata"], ids[key]


def test_a_per_period_value_takes_the_documented_forms(toy):
    def amounts(value):
        data, _ = _period(toy, rainfall=value)
        return [[row[2] for row in data[period]] for period in range(3)]

    assert amounts(0.1) == [[0.1], [0.1], [0.1]]
    assert amounts([0.1, 0.2]) == [[0.1], [0.2], []], "position 0 = the first period"
    assert amounts({1: 0.2}) == [[], [0.2], []], "keys are ZERO-based periods"
    assert amounts({1: {"pond": 0.2}}) == [[], [0.2], []]
    assert amounts({"pond": 0.1}) == [[0.1], [0.1], [0.1]]
    assert amounts({"pond": [0.1, 0.2]}) == [[0.1], [0.2], []]
    assert amounts({"pond": {2: 0.3}}) == [[], [], [0.3]]


def test_an_outlet_rate_is_keyed_by_zero_based_period(toy):
    data, _ = _period(toy, outlets=[LakeOutlet("pond", rate={1: -2.0})])
    assert data == {0: [], 1: [[0, "RATE", -2.0]], 2: []}


def test_what_mf_lak_accepts_without_checking(toy):
    """Ledger 318: documented, not fixed -- a fix updates the docstrings with this."""

    silent, lake_id = _period(toy, key="num", rainfall={101: 0.5})
    assert silent == {0: [], 1: [], 2: []}, "all-int keys are read as periods"
    keyed, _ = _period(toy, key="num", rainfall={lake_id: 0.5})
    assert keyed[0] == [[0, "RAINFALL", 0.5]]

    unsplit, _ = _period(toy, rainfall=np.array([0.1, 0.2, 0.3]))
    assert isinstance(unsplit[0][0][2], np.ndarray), "an array is not split by period"

    base = dict(context=toy["context"], nper=1, lakes=toy["lakes"], lake_id_field="num",
                starting_stage=11.5, lake_bottom=9.0, connection_modes="rectangular", lake_top=11.0)
    table = LakeTableBuilder(area=4.0, bottom=9.0, top=12.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert mf.lak(**base, tables={"101": table}).options["tables"] == [[0, "lak_101.lak.tab"]]
        with pytest.raises(ValueError, match=r"unknown lake '101'"):
            mf.lak(**base, tables={101: table})
        assert mf.lak(**base, nlakes=5).options["nlakes"] == 5


def test_a_pname_through_options_renames_only_what_flopy_writes(toy):
    """Ledger 318: documented, not fixed -- ``name=`` is the way to rename."""

    spec = mf.lak(context=toy["context"], nper=1, lakes=toy["lakes"], lake_id_field="name",
                  starting_stage=11.5, lake_bottom=9.0, connection_modes="rectangular",
                  lake_top=11.0, pname="zz",
                  tables={"pond": LakeTableBuilder(area=4.0, bottom=9.0, top=12.0)})
    assert spec.options["pname"] == "zz"
    assert (spec.name, spec.metadata["lak_index"]["package"]) == ("lak", "lak")
    assert spec.options["lake_tables"][0]["pname"] == "lak_pond_table"
    assert "NOT with ``pname=``/``filename=``" in _entry(type(mf.lak).__call__.__doc__, "**options")


def test_the_conversions_read_the_unit_names_the_docstring_lists(toy):
    call = type(mf.lak).__call__.__doc__
    for units, length, time in (((" FT ", "D"), 3.28081, 86400.0), (("M", "Hr"), 1.0, 3600.0),
                                (("Centimeter", "s"), 100.0, 1.0), (("", "yr"), 1.0, 31557600.0)):
        ctx = ModelContext(grid=toy["context"].grid, domain=toy["context"].domain,
                           length_units=units[0], time_units=units[1])
        options = mf.lak(context=ctx, nper=1, lakes=toy["lakes"], lake_id_field="name",
                         starting_stage=11.5, lake_bottom=9.0, connection_modes="rectangular",
                         lake_top=11.0).options
        assert (options["length_conversion"], options["time_conversion"]) == (length, time), units
    for word in ('``"ft"``', '``"cm"``', '``"m"``'):
        assert word in _entry(call, "length_conversion")
    for word in ('``"d"``', '``"hr"``', '``"min"``', '``"s"``', '``"yr"``'):
        assert word in _entry(call, "time_conversion")
    ctx = ModelContext(grid=toy["context"].grid, domain=toy["context"].domain, length_units="inches")
    with pytest.raises(ValueError, match="Unknown length unit 'inches'"):
        mf.lak(context=ctx, nper=1, lakes=toy["lakes"], lake_id_field="name", starting_stage=11.5,
               lake_bottom=9.0, connection_modes="rectangular", lake_top=11.0)


def test_a_flat_bathy_bed_warns_as_the_first_example_says(toy):
    """The first ``mf.lak`` example used a flat ``lake_bottom`` under ``"bathy"``,
    the configuration the same docstring calls degenerate. It now says why not."""

    with pytest.warns(UserWarning, match="no horizontal"):
        rows = LAKBuilder(context=toy["context"], nper=1, lakes=toy["lakes"], lake_id_field="name",
                          starting_stage=11.5, lake_bottom=9.0, bed_leakance=0.1).connectiondata
    assert not any(row[3] == "HORIZONTAL" for row in rows)
    examples = type(mf.lak).__call__.__doc__.split("Examples")[1]
    first = examples.split(">>> mf.lak(")[1].split(">>>")[0]
    assert "mf.Raster(" in first, "the first example's bed varies by cell"


# ------------------------------------------------------------------ LakeTableBuilder


@pytest.fixture(scope="module")
def bed(tmp_path_factory) -> Path:
    """A 2 x 2 raster of 1 ft pixels: bed 8, 9 / 9, 10."""

    path = tmp_path_factory.mktemp("lakebed") / "bed.tif"
    with rasterio.open(path, "w", driver="GTiff", height=2, width=2, count=1, dtype="float64",
                       crs="EPSG:2927", transform=from_origin(0, 2, 1, 1)) as dst:
        dst.write(np.array([[8.0, 9.0], [9.0, 10.0]]), 1)
    return path


def test_past_the_highest_bed_the_volume_grows_by_area_times_the_coefficient(bed):
    rows = LakeTableBuilder(dem=bed, footprint=box(0, 0, 2, 2), stages=[8, 10, 11, 12],
                            storage_coefficient=0.5).build().rows
    assert rows == ((8.0, 0.0, 0.0), (10.0, 2.0, 3.0), (11.0, 4.0, 4.0), (12.0, 6.0, 4.0))


def test_the_footprint_rules_the_docstring_states(bed, tmp_path):
    stages = [9, 10, 11]
    two = tmp_path / "two.gpkg"
    gpd.GeoDataFrame(geometry=[box(0, 0, 1, 2)], crs=2927).to_file(two, layer="a")
    gpd.GeoDataFrame(geometry=[box(1, 0, 2, 2)], crs=2927).to_file(two, layer="b")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        first_layer = LakeTableBuilder(dem=bed, footprint=two, stages=stages).build().rows
    left_column = LakeTableBuilder(dem=bed, footprint=box(0, 0, 1, 2), stages=stages).build().rows
    assert first_layer == left_column, "only the default (first) layer is read"

    zeros = ((9.0, 0.0, 0.0), (10.0, 0.0, 0.0), (11.0, 0.0, 0.0))
    for footprint in (box(-122.5, 47.0, -122.4, 47.1), box(100, 100, 101, 101), box(0, 0, 0.4, 2)):
        assert LakeTableBuilder(dem=bed, footprint=footprint, stages=stages).build().rows == zeros

    with pytest.raises(ValueError, match="footprint geometry or path"):
        LakeTableBuilder(dem=bed, footprint=str(two), stages=stages).build()
    with pytest.raises(Exception, match="No such file") as missing:
        LakeTableBuilder(dem=bed, footprint=tmp_path / "nope.gpkg", stages=stages).build()
    assert type(missing.value).__name__ == "DataSourceError"


def test_the_refusals_the_rows_and_coefficient_entries_state():
    frame = pd.DataFrame({"stage": [96.0, 98.0], "volume": [0.0, 1.0], "area": [0.0, 2.0]})
    with pytest.raises(ValueError, match="could not convert string to float"):
        LakeTableBuilder(rows=frame).build()
    with pytest.raises(TypeError, match="'int' object is not iterable"):
        LakeTableBuilder(rows=pd.DataFrame([[1, 0, 0], [2, 1, 1]])).build()
    assert LakeTableBuilder(rows=frame.to_numpy()).build().rows == ((96.0, 0.0, 0.0), (98.0, 1.0, 2.0))
    with pytest.raises(ValueError, match="volumes must be nondecreasing"):
        LakeTableBuilder(area=1.0, bottom=9.0, top=11.0, storage_coefficient=-1.0).build()
    # Tables MF6 refuses pass LakeTable (section 2 of the reproducer runs MF6 on them).
    LakeTable([(8.0, 0.0, 0.0), (9.0, 0.0, 0.0), (10.0, 5.0, 4.0)])
    LakeTable([(8.0, 0.0, 5.0), (9.0, 1.0, 4.0)])


# ------------------------------------------------------------------------------ MF6


def _bench(ws: Path, *, connection, perioddata, stage=12.0, nper=1, outlets=None, obs=(),
           extra_connections=(), transient=True, aquifer=5.0):
    """DIS 3 x 3 of 10 ft cells, top 8, one lake on the centre cell (100 ft2).

    A name in ``obs`` is a lake-wide observation; a tuple is passed as given, so
    ``("q", "LAK", 1, 2)`` is connection 2's flow (positive out of the lake).
    ``aquifer`` is both the starting head and the corner CHD; ``transient=False``
    leaves STO out, so every period is steady.
    """

    shutil.rmtree(ws, ignore_errors=True)
    ws.mkdir(parents=True)
    sim = flopy.mf6.MFSimulation(sim_name="t", sim_ws=str(ws), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=[(1.0, 1, 1.0)] * nper)
    flopy.mf6.ModflowIms(sim, complexity="SIMPLE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=3, ncol=3, delr=10, delc=10, top=8.0, botm=-10.0)
    flopy.mf6.ModflowGwfic(gwf, strt=aquifer)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0, icelltype=1)
    if transient:
        flopy.mf6.ModflowGwfsto(gwf, sy=0.1, ss=1e-5, iconvert=1, transient={0: True})
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 0), aquifer]])
    extra = {"noutlets": 1, "outlets": outlets} if outlets else {}
    connections = [connection, *extra_connections]
    flopy.mf6.ModflowGwflak(
        gwf, nlakes=1, packagedata=[[0, stage, len(connections)]],
        connectiondata=[[0, index, *row] for index, row in enumerate(connections)],
        perioddata=perioddata,
        observations={"lak.obs.csv": [
            name if isinstance(name, tuple) else (name.lower().replace("-", ""), name, 1) for name in obs
        ]},
        **extra,
    )
    sim.write_simulation(silent=True)
    proc = subprocess.run(["mf6"], cwd=ws, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout[-2000:]
    return pd.read_csv(ws / "lak.obs.csv")


@pytest.mark.slow
def test_mf6_ignores_the_elevations_written_on_a_vertical_connection(tmp_path):
    """The rule under ``connection_modes``: the lake counts from the CELL TOP (8.0)."""

    constant = {0: [[0, "STATUS", "CONSTANT"], [0, "STAGE", 12.0]]}
    volumes = {
        written: _bench(tmp_path / str(written), obs=("VOLUME",), perioddata=constant,
                        connection=[(0, 1, 1), "VERTICAL", 0.001, written, written, 0.0, 0.0])["VOLUME"].iloc[-1]
        for written in (3.0, 6.0, 0.0)
    }
    assert volumes == {3.0: 400.0, 6.0: 400.0, 0.0: 400.0}


@pytest.mark.slow
def test_mf6_exchanges_through_a_horizontal_face_from_the_bed_below_the_cell_top(tmp_path):
    """What "bathy" writes on a lake cell whose bed (3.0) lies inside it: a VERTICAL
    connection plus a HORIZONTAL face on the SAME cell from the bed to the top (8.0).
    Held at 6.0, below the top, the VERTICAL connection carries nothing and the lake
    holds nothing, while the face exchanges water -- MF6 reads a HORIZONTAL belev. In a
    TRANSIENT period what LEAVES is capped at what is available (nothing is stored
    below the top), while water from the aquifer still comes in."""

    def run(tag, *, leak=1.0, extra=(), **kw):
        held = {0: [[0, "STATUS", "CONSTANT"], [0, "STAGE", 6.0], *extra]}
        return _bench(tmp_path / tag, stage=6.0, perioddata=held,
                      connection=[(0, 1, 1), "VERTICAL", 1.0, 3.0, 3.0, 0.0, 0.0],
                      extra_connections=[[(0, 1, 1), "HORIZONTAL", leak, 3.0, 8.0, 5.0, 10.0]],
                      obs=("VOLUME", ("qv", "LAK", 1, 1), ("qh", "LAK", 1, 2)), **kw).iloc[-1]

    steady = run("steady", transient=False)
    assert steady["QV"] == 0.0 and steady["VOLUME"] == 0.0
    assert steady["QH"] == pytest.approx(3.4641, abs=1e-3), "positive = out of the lake"
    sealed = run("sealed", transient=False, leak=1e-30)
    assert sealed["QH"] == pytest.approx(0.0, abs=1e-20), "the face is what carries it"
    assert run("transient")["QH"] == 0.0, "nothing stored below the top, nothing arrives"
    assert run("inflow", extra=[[0, "INFLOW", 2.0]])["QH"] == pytest.approx(2.0), "only what arrives"
    assert run("aquifer_above", aquifer=7.0)["QH"] == pytest.approx(-5.1475, abs=1e-3), "IN is not capped"


@pytest.mark.slow
def test_a_positive_specified_rate_adds_water_and_an_empty_period_keeps_it(tmp_path):
    connection = [(0, 1, 1), "VERTICAL", 1e-9, 0.0, 0.0, 0.0, 0.0]
    outlet = [[0, 0, -1, "SPECIFIED", 0.0, 0.0, 0.0, 0.0]]
    for rate in (5.0, -5.0):
        obs = _bench(tmp_path / f"rate{rate}", connection=connection, stage=10.0, nper=2,
                     perioddata={0: [[0, "RATE", rate]], 1: []}, outlets=outlet, obs=("STAGE",))
        change = np.diff([10.0, *obs["STAGE"]])
        assert np.allclose(change, rate / 100.0), (rate, list(obs["STAGE"]))
    rain = _bench(tmp_path / "rain", connection=connection, stage=10.0, nper=2,
                  perioddata={0: [[0, "RAINFALL", 1.0]], 1: []}, obs=("RAINFALL",))
    assert list(rain["RAINFALL"]) == [100.0, 100.0], "the empty second period kept the rainfall"
