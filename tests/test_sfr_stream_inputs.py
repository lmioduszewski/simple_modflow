"""What ``mf.sfr``/``SFRBuilder`` accept, read the way MF6 reads it.

Four defects the 2026-09-27 docstring pass documented as current behaviour:

* a numpy integer ``reach_layer`` was read as a COLUMN NAME (``isinstance(x, int)`` is
  False for ``numpy.int64``, ledger 265's class), and so were numpy period and reach
  keys in a period setting;
* ``rainfall``/``evaporation``/``status`` keyed by STREAM landed on that stream's
  FIRST reach only, while the scalar form reaches every reach (ledger 363).
  RAINFALL and EVAPORATION are MF6 rates per unit of streambed area (L/T), so each
  reach takes the stream's value unchanged -- the MF6 run below rains rate x width
  x length on every reach under both spellings;
* a ``StreamDiversion`` between streams not otherwise connected built, and MF6
  refused the package ("diversion 1 from reach 3 to reach 5 is not possible"); the
  diversion is a reach connection now, written once beside a ``StreamConnection``
  naming the same pair, and a diversion INTO reach 0 is written as FloPy's ``-0.0``
  (ledger 364). The MF6 runs read the diverted water arriving;
* a per-reach value that is not a number -- a string naming no column, a mapping
  that matches no reach -- died as a bare ``TypeError`` in ``validate``, or for
  ``reach_top`` went to FloPy as written. It is refused naming the argument, the
  reach and its stream; ``streambed_k``/``roughness`` keep taking a time-series
  name, the only two MF6 reads one for (ledger 365, measured in
  ``docs/reproducers/sfr_stream_inputs.py``).

Every MF6 run is one small model: a 6 x 6 irregular Voronoi grid of ~100 ft cells, one
confined layer (top 110, bottom 0, K 10), CHD at 90 ft on the west and east edges, one
1-day steady period, and two straight streams with beds at 100 ft and a streambed K of
1e-6, so the streams barely leak: ``source`` (4 reaches, 100 ft3/d entering its
headwater) and ``canal`` (2 reaches), 270 ft apart.
"""

from __future__ import annotations

import re
import subprocess
import warnings

import flopy
import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import LineString, Point

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.surface_water_validation import validate_sfr_configuration
from tests.test_sfr_builder import _builder

LINES = {
    "source": LineString([(120.0, 150.0), (480.0, 150.0)]),
    "canal": LineString([(320.0, 420.0), (480.0, 420.0)]),
}
#: Leaves ``source`` at its third reach and enters ``canal`` at its first.
DIVERSION = mf.StreamDiversion(
    "source", "canal", amount=0.25, source_location=Point(300.0, 150.0), receiver_location="upstream"
)


def _sfr_builder(order=("source", "canal"), **values) -> mf.SFRBuilder:
    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    streams = gpd.GeoDataFrame({"sid": list(order)}, geometry=[LINES[name] for name in order], crs=grid.crs)
    settings = dict(
        connection_mode="explicit",
        connections=(mf.StreamConnection("source", None), mf.StreamConnection("canal", None)),
        reach_top=100.0,
        streambed_k=1.0e-6,
        inflow={0: {"source": 100.0}},
    )
    settings.update(values)
    return mf.SFRBuilder(
        context=mf.ModelContext(grid=grid, domain=np.ones((1, grid.ncpl), dtype=int)),
        nper=1, streams=streams, stream_id="sid", **settings,
    )


def _run(root, builder: mf.SFRBuilder, kinds=("UPSTREAM-FLOW",), **options) -> pd.Series:
    """Run the module docstring's model on ``builder``'s SFR; MF6's observations, by name.

    Each ``kinds`` entry is observed at every reach under ``<kind>_<reach>``; MF6
    observation reach ids are one-based, and FloPy writes them as given. The SFR
    budget file's ``RAINFALL`` (volumetric, ft3/d) rides along as ``rain_<reach>``.
    """

    grid = builder.grid
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    edge = np.flatnonzero((xy[:, 0] < 100.0) | (xy[:, 0] > 500.0))
    observations = [
        (f"{kind.lower().replace('-', '_')}_{reach}", kind, reach + 1)
        for kind in kinds
        for reach in range(builder.total_nreaches)
    ]
    sfr = builder.with_updates(
        options={"observations": {"sfr.obs.csv": observations}, **options}
    ).build()
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]), vertices=props["vertices"],
                cell2d=props["cell2d"], top=110.0, botm=[0.0]),
        mf.npf(k=10.0, icelltype=0), mf.ic(strt=90.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), 90.0] for c in edge]}),
        sfr,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success, (root / "p").as_posix()
    workspace = run.model().workspace
    observed = pd.read_csv(f"{workspace}/sfr.obs.csv").iloc[-1].rename(str.lower)
    budget = flopy.utils.CellBudgetFile(f"{workspace}/m_sfr_budget.sfr", precision="double")
    rain = budget.get_data(text="RAINFALL")[-1]
    budget.close()
    return pd.concat([observed, pd.Series({f"rain_{int(node) - 1}": q for node, q in zip(rain["node"], rain["q"], strict=True)})])


# -- [0] a numpy integer is an integer -----------------------------------------------


@pytest.mark.parametrize("layer", [np.int64(0), np.int32(0), np.uint8(0)])
def test_a_numpy_integer_layer_is_a_layer(layer):
    assert _builder(reach_layer=layer).reaches["cellid"].tolist() == _builder(reach_layer=0).reaches["cellid"].tolist()


def test_numpy_period_and_reach_keys_are_periods_and_reaches():
    plain = _builder(rainfall={0: {3: 0.1}}, runoff={0: {4: 2.0}}, inflow={1: {5: 7.0}}).perioddata
    numpy = _builder(
        rainfall={np.int64(0): {np.int64(3): 0.1}},
        runoff={np.int64(0): {np.int64(4): 2.0}},
        inflow={np.int64(1): {np.int64(5): 7.0}},
    ).perioddata
    assert numpy == plain
    assert [3, "RAINFALL", 0.1] in numpy[0] and [5, "INFLOW", 7.0] in numpy[1]


# -- [1] a stream's rate holds on every reach of the stream ---------------------------


@pytest.mark.parametrize(
    ("argument", "keyword", "value"),
    [("rainfall", "RAINFALL", 0.001), ("evaporation", "EVAPORATION", 0.002), ("status", "STATUS", "INACTIVE")],
)
def test_a_stream_keyed_setting_reaches_every_reach_of_the_stream(argument, keyword, value):
    builder = _builder(**{argument: {0: {"main": value}}})
    rows = [row for row in builder.perioddata[0] if row[1] == keyword]
    assert rows == [[rno, keyword, value] for rno in builder.stream_reaches["main"]]
    assert len(builder.stream_reaches["main"]) > 1  # or the test proves nothing


@pytest.mark.parametrize("reach_first", [False, True])
def test_a_reach_key_beats_its_streams_key_in_either_order(reach_first):
    main = _builder().stream_reaches["main"]
    pairs = [("main", 0.001), (main[1], 0.005)]
    builder = _builder(rainfall={0: dict(pairs[::-1] if reach_first else pairs)})
    rows = [row for row in builder.perioddata[0] if row[1] == "RAINFALL"]
    assert rows == [[rno, "RAINFALL", 0.005 if rno == main[1] else 0.001] for rno in main]


def test_inflow_keyed_by_stream_still_enters_at_its_first_reach():
    builder = _builder(inflow={0: {"main": 50.0}})
    assert [row for row in builder.perioddata[0] if row[1] == "INFLOW"] == [
        [builder.stream_reaches["main"][0], "INFLOW", 50.0]
    ]


@pytest.mark.parametrize("argument", ["rainfall", "evaporation", "status", "inflow", "runoff"])
def test_an_unknown_stream_is_refused_naming_the_streams(argument):
    with pytest.raises(KeyError, match=rf"{argument} names stream 'mian'; the streams are \['main', 'tributary'\]"):
        _ = _builder(**{argument: {0: {"mian": 1.0}}}).perioddata


def test_mf6_reads_the_same_rainfall_per_reach_from_a_stream_keyed_rate(tmp_path):
    """Each reach takes the stream's rate, not a share of it: MF6 rains rate x width x
    length on every reach of both streams, exactly as for the scalar form."""

    rate = 0.01
    builder = _sfr_builder(rainfall={0: {"source": rate, "canal": rate}})
    keyed = _run(tmp_path / "keyed", builder)
    scalar = _run(tmp_path / "scalar", _sfr_builder(rainfall=rate))
    rain = [f"rain_{reach}" for reach in range(builder.total_nreaches)]
    area = 10.0 * builder.reaches["rlen"].to_numpy(float)  # the default width, 10 ft
    np.testing.assert_allclose(keyed[rain].to_numpy(float), rate * area, rtol=1e-9)
    np.testing.assert_allclose(keyed[rain].to_numpy(float), scalar[rain].to_numpy(float), rtol=1e-12)


# -- [2] a diversion is a reach connection --------------------------------------------


def test_a_diversion_is_a_reach_connection_written_once():
    bare = _sfr_builder(diversions=(DIVERSION,))
    (_, source, receiver, _), = bare.resolved_diversions
    assert -receiver in bare.connectiondata[source] and source in bare.connectiondata[receiver][1:]
    both = _sfr_builder(
        diversions=(DIVERSION,),
        connections=(
            mf.StreamConnection("source", "canal", source_location=Point(300.0, 150.0), receiver_location="upstream"),
            mf.StreamConnection("canal", None),
        ),
    )
    assert both.resolved_connections == ((source, receiver),)
    assert both.connectiondata == bare.connectiondata
    assert [row[9] for row in both.packagedata] == [len(row) - 1 for row in bare.connectiondata]


def test_a_diversion_into_reach_0_is_written_downstream():
    """Reach 0 has no negative integer; FloPy reads -0.0 as downstream to reach 0, and 0
    as upstream FROM it -- which MF6 refused as the same broken diversion."""

    builder = _sfr_builder(order=("canal", "source"), diversions=(DIVERSION,))
    (_, source, receiver, _), = builder.resolved_diversions
    assert receiver == 0
    downstream = builder.connectiondata[source][-2]
    assert downstream == 0 and np.signbit(downstream)


@pytest.mark.parametrize("form", ["diversion alone", "with its StreamConnection", "into reach 0"])
def test_mf6_runs_a_diversion_and_the_water_arrives(tmp_path, form):
    """MF6 refused the first and last forms before; the second was the workaround.
    A FRACTION of 0.25 at the third reach of ``source`` puts a quarter of its flow
    into the canal's first reach and leaves three quarters in ``source``."""

    values = {"diversions": (DIVERSION,)}
    if form == "with its StreamConnection":
        values["connections"] = (
            mf.StreamConnection("source", "canal", source_location=Point(300.0, 150.0), receiver_location="upstream"),
            mf.StreamConnection("canal", None),
        )
    if form == "into reach 0":
        values["order"] = ("canal", "source")
    builder = _sfr_builder(**values)
    flows = _run(tmp_path, builder, kinds=("UPSTREAM-FLOW", "DOWNSTREAM-FLOW"))
    (_, source, receiver, _), = builder.resolved_diversions
    leaving = -flows[f"downstream_flow_{source}"]
    assert leaving == pytest.approx(100.0, rel=1e-3)
    assert flows[f"upstream_flow_{receiver}"] == pytest.approx(0.25 * leaving, rel=1e-4)
    after = builder.stream_reaches["source"][builder.stream_reaches["source"].index(source) + 1]
    assert flows[f"upstream_flow_{after}"] == pytest.approx(0.75 * leaving, rel=1e-4)


#: Networks MF6's ``sfr_check_ustrf`` stops on (ledger 366), and what
#: ``validate`` says for each. ``mid`` is ``source``'s third reach.
_MID = Point(300.0, 150.0)
UNRUNNABLE = {
    "default source_location: the last reach": (
        {"diversions": (mf.StreamDiversion("source", "canal", amount=0.25),)},
        r"SFR diversion from reach 3 \(stream 'source'\) to reach 5 \(stream 'canal'\): MF6 diverts only "
        r"from a reach that has a downstream reach of its own .* source_location above",
        r"Upstream fractions for reach +4 is not equal to one \(0\.0",
    ),
    "into a reach another reach flows into": (
        {"diversions": (mf.StreamDiversion("source", "canal", amount=0.25, source_location=_MID, receiver_location="downstream"),)},
        r"SFR reach 5 \(stream 'canal'\) receives a diversion from reach 2 \(stream 'source'\) and is also "
        r"the only reach that reach 4 \(stream 'canal'\) flows into: .* receiver_location='upstream'",
        r"Upstream fractions for reach +5 is not equal to one \(0\.0",
    ),
    "a StreamConnection above the last reach": (
        {"connections": (
            mf.StreamConnection("source", "canal", source_location=_MID, receiver_location="upstream"),
            mf.StreamConnection("canal", None),
        )},
        r"SFR reach 2 \(stream 'source'\) flows into reach 3 \(stream 'source'\), reach 4 \(stream "
        r"'canal'\), whose upstream fractions \(USTRF\) sum to 2; MF6 requires one",
        r"Upstream fractions for reach +3 is not equal to one \(2\.0",
    ),
    # A diversion is a reach connection (ledger 364), so one leaving and entering
    # the same reach connected it to itself; nothing refused it before.
    "a diversion onto its own headwater reach": (
        {"diversions": (mf.StreamDiversion(
            "source", "source", amount=0.25, source_location="upstream", receiver_location="upstream"),)},
        r"SFR StreamDiversion\('source', 'source'\) leaves and enters the same reach, reach 0 "
        r"\(stream 'source'\).*cannot connect a reach to itself",
        r"Could not fill SFR package IA and JA connection data",
    ),
    "a diversion onto the canal's first reach": (
        {"diversions": (mf.StreamDiversion(
            "canal", "canal", amount=0.25, source_location="upstream", receiver_location="upstream"),)},
        r"SFR StreamDiversion\('canal', 'canal'\) leaves and enters the same reach, reach 4 "
        r"\(stream 'canal'\)",
        r"Could not fill SFR package IA and JA connection data",
    ),
}


def test_a_diversion_onto_its_own_reach_is_named_before_the_upstream_fractions():
    """Mid-stream, reach 1 flows into the self-diverting reach alone, which the USTRF
    check refused blaming reach 1; the refusal names the diversion now."""

    diversion = mf.StreamDiversion("source", "source", amount=0.25, source_location=_MID, receiver_location=_MID)
    with pytest.raises(ValueError, match=r"leaves and enters the same reach, reach 2 \(stream 'source'\)"):
        _sfr_builder(diversions=(diversion,)).build()


@pytest.mark.parametrize("form", list(UNRUNNABLE))
def test_a_network_mf6_cannot_run_is_refused_naming_the_reaches(form):
    values, refusal, _ = UNRUNNABLE[form]
    with pytest.raises(ValueError, match=refusal):
        _sfr_builder(**values).build()


@pytest.mark.parametrize("form", list(UNRUNNABLE))
def test_mf6_stops_on_each_network_validate_refuses(tmp_path, form):
    """The same tables written past ``validate`` (``mf.sfr.flopy``): MF6 stops on
    each, which is what makes the refusal a description of MF6 and not a guess."""

    values, _, mf6_says = UNRUNNABLE[form]
    builder = _sfr_builder(**values)
    grid = builder.grid
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    edge = np.flatnonzero((xy[:, 0] < 100.0) | (xy[:, 0] > 500.0))
    sfr = mf.sfr.flopy(
        packagedata=builder.packagedata, connectiondata=builder.connectiondata,
        perioddata=builder.perioddata, diversions=builder.diversiondata or None,
    )
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]), vertices=props["vertices"],
                cell2d=props["cell2d"], top=110.0, botm=[0.0]),
        mf.npf(k=10.0, icelltype=0), mf.ic(strt=90.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), 90.0] for c in edge]}),
        sfr,
        mf.oc(head_filerecord="m.hds", saverecord=[("HEAD", "ALL")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = mf.Project(tmp_path / "p", name="p").run("r", simulation)
    assert not run.success
    listing_path = next((tmp_path / "p").rglob("mfsim.lst"))
    for _ in range(5):
        listing = " ".join(listing_path.read_text().split())
        if "ERROR REPORT" in listing:
            break
        # MF6 6.8 sometimes SEGFAULTS on a reach connected to itself before it writes
        # its report (docs/reproducers/sfr_stream_inputs.py E counts how often); run
        # it again to read what it objects to.
        subprocess.run(["mf6"], cwd=listing_path.parent, capture_output=True, check=False)
    assert re.search(mf6_says, listing), listing[-600:]


# -- [47] a per-reach value is a number, or a time-series name where MF6 takes one ----

NUMERIC_ONLY = [("width", "RWID"), ("gradient", "RGRD"), ("reach_top", "RTP"), ("streambed_thickness", "RBTH")]
SERIES_TOO = [("streambed_k", "RHK"), ("roughness", "MAN")]


@pytest.mark.parametrize(("argument", "mf6_name"), NUMERIC_ONLY)
def test_a_string_naming_no_column_is_refused_where_mf6_reads_a_number(argument, mf6_name):
    with pytest.raises(
        ValueError,
        match=rf"SFR {argument} \({mf6_name}\) for reach 0 \(stream 'tributary'\): 'ts1' names no column "
        rf".*MF6 reads {mf6_name} as a number",
    ):
        _builder(**{argument: "ts1"}).build()


@pytest.mark.parametrize(("argument", "mf6_name"), NUMERIC_ONLY + SERIES_TOO)
def test_a_mapping_matching_no_reach_is_refused_naming_the_stream(argument, mf6_name):
    with pytest.raises(
        ValueError,
        match=rf"SFR {argument} \({mf6_name}\) for reach 0 \(stream 'tributary'\): the {argument}= mapping names none",
    ):
        _builder(**{argument: {"source": 4.0, "canal": 8.0}}).build()


@pytest.mark.parametrize(("argument", "mf6_name"), NUMERIC_ONLY + SERIES_TOO)
def test_an_empty_column_is_refused_naming_the_column(argument, mf6_name):
    streams = _builder().streams.assign(value=[np.nan, 2.0])
    with pytest.raises(ValueError, match=r"for reach 0 \(stream 'tributary'\): column 'value' of the stream table holds nan"):
        _builder(streams=streams, **{argument: "value"}).build()


@pytest.mark.parametrize(("argument", "mf6_name"), SERIES_TOO)
def test_a_time_series_name_is_kept_where_mf6_takes_one(argument, mf6_name):
    spec = _builder(**{argument: "ts1"}).build()
    position = {"streambed_k": 7, "roughness": 8}[argument]
    assert {row[position] for row in spec.options["packagedata"]} == {"ts1"}


CODES = {
    "width": "sfr_reach_width", "gradient": "sfr_reach_gradient", "reach_top": "sfr_reach_top",
    "streambed_thickness": "sfr_streambed_thickness", "streambed_k": "sfr_streambed_k",
    "roughness": "sfr_mannings",
}


@pytest.mark.parametrize(("argument", "mf6_name"), SERIES_TOO)
def test_validate_sfr_configuration_takes_a_time_series_name_where_mf6_does(argument, mf6_name):
    # It reported 'must be non-negative' / 'must be positive' for every reach,
    # though MF6 runs RHK and MAN given as a time series (ledger 365).
    report = validate_sfr_configuration(_builder(**{argument: "ts1"}))
    assert report.ok, report.to_frame().to_string()


@pytest.mark.parametrize(("argument", "mf6_name"), NUMERIC_ONLY)
def test_validate_sfr_configuration_refuses_a_string_where_mf6_reads_a_number(argument, mf6_name):
    report = validate_sfr_configuration(_builder(**{argument: "ts1"}))
    assert {issue.code for issue in report.errors} == {CODES[argument]}
    assert all("not a time-series name" in issue.message for issue in report.errors)


@pytest.mark.parametrize(("argument", "mf6_name"), NUMERIC_ONLY + SERIES_TOO)
@pytest.mark.parametrize("value", ["ts1", np.nan, None, 0.0, -1.0, np.float64(2.0), True])
def test_validate_sfr_configuration_and_the_builder_agree(argument, mf6_name, value):
    """One rule: the report finds a fault in a per-reach value exactly when ``validate`` raises."""

    builder = _builder(streams=_builder().streams.assign(value=value), **{argument: "value"})
    try:
        builder.validate()
        refused = False
    except ValueError:
        refused = True
    codes = {issue.code for issue in validate_sfr_configuration(builder).errors}
    assert (CODES[argument] in codes) == refused


def test_the_bounds_name_the_stream_and_the_value():
    with pytest.raises(ValueError, match=r"SFR reach 3 \(stream 'main'\) has nonpositive width \(RWID\): 0.0"):
        _builder(width={"tributary": 4.0, "main": 0.0}).build()


def test_mf6_runs_streambed_k_and_roughness_given_as_time_series(tmp_path):
    builder = _sfr_builder(streambed_k="kts", roughness="mts")
    series = {
        "filename": "sfr.ts",
        "timeseries": [(0.0, 1.0e-6, 0.03), (1.0, 1.0e-6, 0.03)],
        "time_series_namerecord": [("kts", "mts")],
        "interpolation_methodrecord": ["linear", "linear"],
    }
    flows = _run(tmp_path, builder, kinds=("DOWNSTREAM-FLOW",), timeseries=series)
    assert -flows["downstream_flow_3"] == pytest.approx(100.0, rel=1e-3)
