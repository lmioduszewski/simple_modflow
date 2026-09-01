"""The grid-independent GIS export, and the CLN fields it depends on.

Like ``test_usg_import.py`` every fixture here is SYNTHESIZED -- the model this
was built against is proprietary and 140 MB. The synthetic network is chosen to
carry the three things that actually broke: a feature name with a digit in it
(``Wet217``, which a regex-based number scan turns into a phantom column), a
tributary that *touches* its trunk (so grouping by connected component merges
two named features into one), and a centerline digitized mouth-first (so a
station projected onto it decreases downstream).
"""

from __future__ import annotations

import numpy as np
import pytest

flopy = pytest.importorskip("flopy")
gpd = pytest.importorskip("geopandas")
pytest.importorskip("rasterio")

from shapely.geometry import LineString  # noqa: E402

from myflopy.modflow.usg import read_usg  # noqa: E402
from myflopy.modflow.usg._io import free_floats, free_row  # noqa: E402

# A 3 x 3 plan mesh of 10 ft cells, two layers thick. Big enough for a stream to
# run along one row and for a lake to be a genuine 2-D patch -- a two-node pond
# is a chain, and the classifier is right to say so.
NCPL, NLAY, NODES = 9, 2, 18
PLAN = [(x, y) for y in (0.0, 10.0, 20.0, 30.0) for x in (0.0, 10.0, 20.0, 30.0)]
CORNERS = [
    [r * 4 + c, r * 4 + c + 1, (r + 1) * 4 + c + 1, (r + 1) * 4 + c]
    for r in range(3)
    for c in range(3)
]
CENTERS = [(c * 10.0 + 5.0, r * 10.0 + 5.0) for r in range(3) for c in range(3)]
LEVELS = [10.0, 5.0, 0.0]
PLAN_NEIGHBOURS = {
    i: [
        j
        for j in (i - 3, i + 3, i - 1 if i % 3 else None, i + 1 if i % 3 != 2 else None)
        if j is not None and 0 <= j < 9
    ]
    for i in range(9)
}


def _connectivity():
    iac, ja = [], []
    for node in range(NODES):
        layer, cell = divmod(node, NCPL)
        neighbours = [layer * NCPL + n for n in PLAN_NEIGHBOURS[cell]]
        neighbours += [node + NCPL] if layer == 0 else [node - NCPL]
        iac.append(1 + len(neighbours))
        ja += [node + 1] + [n + 1 for n in neighbours]
    return np.array(iac), np.array(ja)


def _write_gsf(path):
    lines = ["UNSTRUCTURED", f"{NODES} 1 1 1", str(len(PLAN) * len(LEVELS))]
    for z in LEVELS:
        lines += [f"{x:.6f} {y:.6f} {z:.6f}" for x, y in PLAN]
    for node in range(NODES):
        layer, cell = divmod(node, NCPL)
        xc, yc = CENTERS[cell]
        zc = (LEVELS[layer] + LEVELS[layer + 1]) / 2.0
        top = [layer * len(PLAN) + c + 1 for c in CORNERS[cell]]
        bot = [(layer + 1) * len(PLAN) + c + 1 for c in CORNERS[cell]]
        lines.append(
            f"{node + 1} {xc:.6f} {yc:.6f} {zc:.6f} {layer + 1} 8 "
            + " ".join(str(i) for i in top + bot)
        )
    path.write_text("\n".join(lines) + "\n")


def _write_cln(path):
    """A CLN with a 3-node trunk, a 1-node tributary touching it, and a 5-node lake.

    Nodes 1-3 are ``Creek`` running along the bottom row; node 4 is ``Wet217``,
    joined to node 3, so the two features touch and a component-based grouping
    would merge them. Nodes 5-9 are ``Pond``, mutually connected, so its mean
    degree of 4 reads as a 2-D patch rather than a chain.
    """

    # Nodes 1-3 are a chain (Creek), node 4 hangs off node 3 (Wet217), and
    # nodes 5-9 are mutually connected (Pond) -- mean degree 4, a 2-D patch.
    pond = [5, 6, 7, 8, 9]
    adjacency = {1: [2], 2: [1, 3], 3: [2, 4], 4: [3]}
    for node in pond:
        adjacency[node] = [n for n in pond if n != node]
    iac = [1 + len(adjacency[n]) for n in range(1, 10)]
    ja = []
    for node in range(1, 10):
        ja += [node] + adjacency[node]

    rows = ["         0         9         0         0         0         0         9         2",
            f"      {len(ja)}   NJA_CLN",
            "INTERNAL  1  (FREE)  1  IAC_CLN",
            " ".join(str(v) for v in iac),
            "INTERNAL  1  (FREE)  1  JA_CLN",
            " ".join(str(v) for v in ja)]
    # IFNO IFTYP IFDIR FLENG FELEV FANGLE IFLIN ICCWADI  <label>
    nodes = [(1, 1, 1, 10.0, 9.0, -999.0, -4, 0, "Creek"),
             (2, 1, 1, 10.0, 8.0, -999.0, -4, 0, "Creek"),
             (3, 1, 1, 10.0, 7.0, -999.0, -4, 0, "Creek"),
             (4, 1, 1, 10.0, 8.5, -999.0, -4, 0, "Wet217"),
             (5, 2, 1, 10.0, 6.0, -999.0, -4, 0, "Pond"),
             (6, 2, 1, 10.0, 6.0, -999.0, -4, 0, "Pond"),
             (7, 2, 1, 10.0, 6.0, -999.0, -4, 0, "Pond"),
             (8, 2, 1, 10.0, 6.0, -999.0, -4, 0, "Pond"),
             (9, 2, 1, 10.0, 6.0, -999.0, -4, 0, "Pond")]
    for row in nodes:
        rows.append("  ".join(f"{v:g}" if isinstance(v, float) else str(v) for v in row))
    # IFNO IGWNOD IFCON FSKIN FLENGW FANISO ICGWADI  <label>
    connections = [(1, 1, 3, 2.5, 10.0, 1.0, 0, "Creek"),
                   (2, 2, 3, 2.5, 10.0, 1.0, 0, "Creek"),
                   (3, 3, 3, 2.5, 10.0, 1.0, 0, "Creek"),
                   (4, 6, 3, 2.5, 10.0, 1.0, 0, "Wet217"),
                   (5, 4, 3, 0.5, 10.0, 1.0, 0, "Pond"),
                   (6, 5, 3, 4.5, 10.0, 1.0, 0, "Pond"),
                   (7, 7, 3, 2.5, 10.0, 1.0, 0, "Pond"),
                   (8, 8, 3, 2.5, 10.0, 1.0, 0, "Pond"),
                   (9, 9, 3, 2.5, 10.0, 1.0, 0, "Pond")]
    for row in connections:
        rows.append("  ".join(f"{v:g}" if isinstance(v, float) else str(v) for v in row))
    rows += ["1  6.5  1.0e-6   Creek type", "2  7.0  2.0e-6   Pond type"]
    rows += ["CONSTANT 1  (FREE)  1  IBOUND", "CONSTANT  8.0  (FREE)  -1  Initial Heads"]
    path.write_text("\n".join(rows) + "\n")


@pytest.fixture(scope="module")
def workspace(tmp_path_factory):
    """A small USG model carrying a CLN, written to a temporary directory."""

    ws = tmp_path_factory.mktemp("usg_export")
    iac, ja = _connectivity()
    model = flopy.mfusg.MfUsg(modelname="tiny", model_ws=str(ws), structured=False)
    flopy.mfusg.MfUsgDisU(
        model, nodes=NODES, nlay=NLAY, njag=int(iac.sum()), ivsd=-1, nper=4,
        itmuni=4, lenuni=1, idsymrd=0, nodelay=NCPL,
        top=[np.full(NCPL, 10.0), np.full(NCPL, 5.0)],
        bot=[np.full(NCPL, 5.0), np.full(NCPL, 0.0)],
        area=100.0, iac=iac, ja=ja, cl12=5.0, fahl=50.0,
        perlen=[30.0] * 4, nstp=[1] * 4, tsmult=[1.0] * 4, steady=[False] * 4,
    )
    flopy.mfusg.MfUsgBas(model, ibound=1, strt=8.0, hnoflo=999.0)
    flopy.mfusg.MfUsgLpf(model, laytyp=4, layvka=0, chani=1.0, hk=10.0, vka=1.0, ss=1e-5, sy=0.2)
    model.write_input()

    (ws / "tiny.drn").write_text(
        "# drn\n  3  0\n  3  0  Stress Period 1\n"
        "  1  7.5  100.0\n  2  7.4  120.0\n  3  7.3  140.0\n"
        " -1  0  Stress Period 2\n -1  0  Stress Period 3\n -1  0  Stress Period 4\n"
    )
    # Recharge repeats with a period of 2, so the cycle detector has something
    # to find and the export must write 2 rasters rather than 4.
    (ws / "tiny.rch").write_text(
        "# rch\n  1  0\n  1\nCONSTANT  0.001  (FREE)  -1  Stress Period 1\n"
        "  1\nCONSTANT  0.002  (FREE)  -1  Stress Period 2\n"
        "  1\nCONSTANT  0.001  (FREE)  -1  Stress Period 3\n"
        "  1\nCONSTANT  0.002  (FREE)  -1  Stress Period 4\n"
    )
    # CLN wells: precipitation minus evaporation over the pond.
    (ws / "tiny.wel").write_text(
        "# wel on CLNs\n  5  0  5\n"
        "  0  0  5  Stress Period 1\n" + "".join(f"  {n}  -1.0\n" for n in range(5, 10)) +
        "  0  0  5  Stress Period 2\n" + "".join(f"  {n}  2.0\n" for n in range(5, 10)) +
        "  0  0  5  Stress Period 3\n" + "".join(f"  {n}  -1.0\n" for n in range(5, 10)) +
        "  0  0  5  Stress Period 4\n" + "".join(f"  {n}  2.0\n" for n in range(5, 10))
    )
    (ws / "tiny.oc").write_text("HEAD SAVE UNIT 30\nPeriod 1 Step 1\nSAVE HEAD\n")
    _write_cln(ws / "tiny.cln")
    _write_gsf(ws / "tiny.gsf")
    (ws / "tiny.nam").write_text(
        "LIST 7 tiny.lst\nBAS6 1 tiny.bas\nDISU 10 tiny.disu\nLPF 11 tiny.lpf\n"
        "OC 22 tiny.oc\nDRN 13 tiny.drn\nRCH 12 tiny.rch\nWEL 14 tiny.wel\n"
        "CLN 71 tiny.cln\n"
    )
    return ws


@pytest.fixture(scope="module")
def usg(workspace):
    return read_usg(workspace / "tiny.nam", gsf=workspace / "tiny.gsf", crs="EPSG:2927")


@pytest.fixture(scope="module")
def exported(usg, tmp_path_factory):
    """The model exported once, with a mouth-first reference centerline."""

    out = tmp_path_factory.mktemp("export_out") / "gis"
    reference = out.parent / "creek.gpkg"
    # Digitized from the DOWNSTREAM end, which is what real hydrography does.
    gpd.GeoDataFrame(
        {"geometry": [LineString([(25.0, 5.0), (15.0, 5.0), (5.0, 5.0)])]},
        crs="EPSG:2927",
    ).to_file(reference, driver="GPKG")
    manifest = usg.export_gis(out, stream_lines={"Creek": reference}, overwrite=True)
    return out, manifest


# --- the parsing fix the names depend on ----------------------------------


def test_a_label_with_a_digit_is_not_read_as_a_number():
    """``free_row`` stops at the first non-numeric token; ``free_floats`` does not.

    ``Wet217`` contributes a phantom ``217.0`` to a regex scan, giving that row
    one more number than its neighbours -- the same shape of bug that once read
    an ``IPRN`` flag as an array multiplier.
    """

    row = "  4  1  1  10  8.5  -999  -4  0  Wet217"
    assert len(free_floats(row)) == 9  # the phantom
    numbers, label = free_row(row)
    assert len(numbers) == 8
    assert label == "Wet217"


def test_free_row_on_an_all_numeric_line():
    numbers, label = free_row("  1  2.5  -3.0e2  ")
    assert numbers == [1.0, 2.5, -300.0]
    assert label is None


# --- what the CLN reader now keeps ----------------------------------------


def test_features_are_split_by_name_not_by_component(usg):
    """The tributary touches the trunk, so components would merge them."""

    labels = sorted(f.label for f in usg.cln.features)
    assert labels == ["Creek", "Pond", "Wet217"]


def test_fskin_and_the_conduit_table_survive_the_read(usg):
    """``FSKIN`` is the only calibrated number a CLN carries; it used to be dropped."""

    by_name = {f.label: f for f in usg.cln.features}
    assert np.allclose(by_name["Creek"].fskin, 2.5)
    assert np.allclose(sorted(by_name["Pond"].fskin), [0.5, 2.5, 2.5, 2.5, 4.5])
    assert np.allclose(usg.cln.radii, [6.5, 7.0])
    assert np.allclose(usg.cln.conductivities, [1e-6, 2e-6])


# --- the export ------------------------------------------------------------


def test_export_writes_the_expected_files(exported):
    out, manifest = exported
    for name in ("boundaries.gpkg", "streams.gpkg", "lakes.gpkg", "cell_values.gpkg",
                 "periods.csv", "manifest.json", "README.md"):
        assert (out / name).is_file(), name
    assert manifest.crs == "EPSG:2927"
    assert "boundaries.gpkg:drn" in manifest.layers
    assert "boundaries.gpkg:drn_lines" in manifest.layers


def test_conductance_is_recoverable_from_the_per_foot_twin(exported):
    """The invariant that makes a line BC survive a change of grid.

    ``cond`` itself cannot be copied onto a finer mesh without inflating the
    total; ``cond_per_ft`` times the new length can.
    """

    out, _ = exported
    frame = gpd.read_file(out / "boundaries.gpkg", layer="drn_lines")
    rebuilt = float((frame["conductance_per_ft"] * frame["length_ft"]).sum())
    assert rebuilt == pytest.approx(float(frame["conductance"].sum()), rel=1e-9)
    assert rebuilt == pytest.approx(360.0, rel=1e-9)


def test_segments_tile_the_feature(exported):
    """Consecutive segments meet exactly, so lengths sum to the whole line."""

    out, _ = exported
    frame = gpd.read_file(out / "boundaries.gpkg", layer="drn_lines")
    for _, chain in frame.groupby(["layer", "chain"]):
        merged = chain.geometry.union_all()
        assert float(chain["length_ft"].sum()) == pytest.approx(merged.length, abs=1e-6)


def test_station_increases_downstream_on_a_mouth_first_centerline(exported):
    """A reference line digitized mouth-first must not invert the profile.

    Both real Ten Trails creeks are digitized this way; projecting onto them raw
    gives a station that decreases downstream, so the bed appears to climb.
    """

    out, _ = exported
    nodes = gpd.read_file(out / "streams.gpkg", layer="stream_nodes")
    creek = nodes[nodes["stream_id"] == "Creek"].sort_values("station")
    assert list(creek["order"]) == sorted(creek["order"])
    assert (np.diff(creek["bed_elev"].to_numpy()) <= 0).all()


def test_streams_carry_sfr_column_names(exported):
    out, _ = exported
    streams = gpd.read_file(out / "streams.gpkg", layer="streams")
    assert set(streams.columns) >= {"stream_id", "rwid", "rhk", "rgrd", "rbth", "man"}
    creek = streams[streams["stream_id"] == "Creek"].iloc[0]
    assert creek["rwid"] == pytest.approx(13.0)  # 2 x FRAD
    assert creek["rhk"] == pytest.approx(2.5)  # FSKIN
    assert creek.geometry.has_z


def test_lakes_carry_lak_column_names_and_the_bedleak_assumption(exported):
    out, manifest = exported
    lakes = gpd.read_file(out / "lakes.gpkg", layer="lakes")
    assert set(lakes.columns) >= {"lake_id", "strt", "lake_bottom", "bedleak", "fskin_median"}
    pond = lakes[lakes["lake_id"] == "Pond"].iloc[0]
    assert pond["bedleak"] == pytest.approx(2.5)  # median of 0.5, 2.5, 2.5, 2.5, 4.5
    assert any("bed thickness" in note for note in manifest.notes)


def test_bed_leakance_spread_is_reported_not_hidden(exported):
    """One ``bed_leakance`` per lake is a real loss when FSKIN varies 9x."""

    _, manifest = exported
    assert any("fskin varies more than 2x" in note for note in manifest.notes)


def test_a_repeated_period_is_written_once(exported):
    """Recharge alternates between two arrays over four periods."""

    import pandas as pd

    out, manifest = exported
    periods = pd.read_csv(out / "periods.csv")
    assert list(periods["recharge_array"]) == [1, 2, 1, 2]
    assert manifest.rasters["arrays/recharge_NN.tif"]["files"] == 2
    assert (out / "arrays" / "recharge_01.tif").is_file()
    assert not (out / "arrays" / "recharge_03.tif").exists()


def test_lake_forcing_is_a_rate_over_the_lake_area(exported):
    """CLN wells become LAK rainfall/evaporation, which is what mf.lak takes."""

    import pandas as pd

    out, _ = exported
    forcing = pd.read_csv(out / "lake_forcing.csv")
    pond = forcing[forcing["lake_id"] == "Pond"].sort_values("period")
    # Five cells of 100 ft2; -5.0 ft3/d over 500 ft2 is -0.01 ft/d.
    assert pond["rate"].tolist() == pytest.approx([-0.01, 0.02, -0.01, 0.02])
    assert pond["evaporation"].tolist() == pytest.approx([0.01, 0.0, 0.01, 0.0])
    assert pond["rainfall"].tolist() == pytest.approx([0.0, 0.02, 0.0, 0.02])


def test_cell_values_are_exact_where_the_rasters_are_not(exported):
    """The points carry every array without resampling."""

    out, _ = exported
    points = gpd.read_file(out / "cell_values.gpkg", layer="cell_values")
    assert len(points) == NCPL
    assert {"top", "botm_1", "botm_2", "k_1", "recharge_01"} <= set(points.columns)
    assert points["top"].tolist() == pytest.approx([10.0] * NCPL)


def test_export_refuses_a_non_empty_directory(usg, tmp_path):
    (tmp_path / "already").mkdir()
    (tmp_path / "already" / "something.txt").write_text("x")
    with pytest.raises(FileExistsError, match="overwrite=True"):
        usg.export_gis(tmp_path / "already")


def test_export_without_a_grid_is_refused(workspace, tmp_path):
    """A .gsf is where a USG model's coordinates live; without one there is
    nothing to write, and the refusal must say so rather than writing empties."""

    import shutil

    bare = tmp_path / "no_gsf"
    bare.mkdir()
    for item in workspace.iterdir():
        if item.is_file() and item.suffix != ".gsf":
            shutil.copy(item, bare / item.name)
    model = read_usg(bare / "tiny.nam")
    assert model.grid is None
    with pytest.raises(ValueError, match="needs the grid"):
        model.export_gis(tmp_path / "nope")
