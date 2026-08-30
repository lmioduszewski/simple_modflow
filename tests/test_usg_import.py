"""The MODFLOW-USG importer, exercised on a model small enough to read by eye.

Every fixture here is SYNTHESIZED. The importer was developed against a
proprietary 140 MB model that CI will never have, so the tests build their own
two-layer, four-cell USG model instead -- small enough that the expected answer
can be written down, and complete enough to exercise the paths that actually
broke during development: the ``OPEN/CLOSE`` multiplier, ``ITMP < 0`` reuse,
node-to-cellid conversion, and layered-ness detection.
"""

from __future__ import annotations

import numpy as np
import pytest

flopy = pytest.importorskip("flopy")

from myflopy.modflow.usg import read_usg  # noqa: E402
from myflopy.modflow.usg._io import ArrayCursor, NameFile  # noqa: E402
from myflopy.modflow.usg.packages import layered_report  # noqa: E402

# A 2 x 2 plan mesh of 10 ft squares, two layers thick.
NCPL, NLAY, NODES = 4, 2, 8
PLAN = [(0, 0), (10, 0), (20, 0), (0, 10), (10, 10), (20, 10), (0, 20), (10, 20), (20, 20)]
CORNERS = [[0, 1, 4, 3], [1, 2, 5, 4], [3, 4, 7, 6], [4, 5, 8, 7]]
CENTERS = [(5, 5), (15, 5), (5, 15), (15, 15)]
LEVELS = [10.0, 5.0, 0.0]  # layer 1 spans 10->5, layer 2 spans 5->0

# Plan neighbours: 0-1, 0-2, 1-3, 2-3.
PLAN_NEIGHBOURS = {0: [1, 2], 1: [0, 3], 2: [0, 3], 3: [1, 2]}


def _connectivity():
    """Return ``(iac, ja)`` for the stacked 2 x 2 x 2 mesh, 1-based."""

    iac, ja = [], []
    for node in range(NODES):
        layer, cell = divmod(node, NCPL)
        neighbours = [layer * NCPL + n for n in PLAN_NEIGHBOURS[cell]]
        neighbours += [node + NCPL] if layer == 0 else [node - NCPL]
        iac.append(1 + len(neighbours))
        ja += [node + 1] + [n + 1 for n in neighbours]
    return np.array(iac), np.array(ja)


def _write_gsf(path):
    """Write the .gsf that carries the geometry a DISU does not."""

    lines = ["UNSTRUCTURED", f"{NODES} 1 1 1", str(len(PLAN) * len(LEVELS))]
    for z in LEVELS:
        lines += [f"{x:.6f} {y:.6f} {z:.6f}" for x, y in PLAN]
    for node in range(NODES):
        layer, cell = divmod(node, NCPL)
        xc, yc = CENTERS[cell]
        zc = (LEVELS[layer] + LEVELS[layer + 1]) / 2.0
        top = [layer * len(PLAN) + c + 1 for c in CORNERS[cell]]
        bot = [(layer + 1) * len(PLAN) + c + 1 for c in CORNERS[cell]]
        ids = " ".join(str(i) for i in top + bot)
        lines.append(f"{node + 1} {xc:.6f} {yc:.6f} {zc:.6f} {layer + 1} 8 {ids}")
    path.write_text("\n".join(lines) + "\n")


@pytest.fixture(scope="module")
def usg_workspace(tmp_path_factory):
    """Build a complete little MODFLOW-USG model and return its name file."""

    ws = tmp_path_factory.mktemp("usg")
    iac, ja = _connectivity()

    model = flopy.mfusg.MfUsg(modelname="tiny", model_ws=str(ws), structured=False)
    flopy.mfusg.MfUsgDisU(
        model,
        nodes=NODES,
        nlay=NLAY,
        njag=int(iac.sum()),
        ivsd=-1,
        nper=2,
        itmuni=4,
        lenuni=1,
        idsymrd=0,
        nodelay=NCPL,
        top=[np.full(NCPL, 10.0), np.full(NCPL, 5.0)],
        bot=[np.full(NCPL, 5.0), np.full(NCPL, 0.0)],
        area=100.0,
        iac=iac,
        ja=ja,
        cl12=5.0,
        fahl=50.0,
        perlen=[30.0, 30.0],
        nstp=[1, 1],
        tsmult=[1.0, 1.0],
        steady=[False, False],
    )
    flopy.mfusg.MfUsgBas(model, ibound=1, strt=8.0, hnoflo=999.0)
    flopy.mfusg.MfUsgLpf(model, laytyp=4, layvka=0, chani=1.0, hk=10.0, vka=1.0, ss=1e-5, sy=0.2)
    flopy.mfusg.MfUsgSms(model, hclose=0.01, hiclose=0.001, mxiter=50, iter1=100, nonlinmeth=1)
    model.write_input()

    # Hand-written list/array packages: these carry the reuse conventions the
    # importer has to honour, which FloPy's writers do not expose.
    (ws / "tiny.drn").write_text(
        "# drn\n"
        "  2  0\n"
        "  2  0  Stress Period 1\n"
        "  1  7.5  100.0\n"
        "  2  7.5  120.0\n"
        " -1  0  Stress Period 2 (reuse)\n"
    )
    (ws / "tiny.chd").write_text(
        "# chd\n"
        "  2\n"
        "  2  Stress Period 1\n"
        "  3  9.0  9.0\n"
        "  4  9.0  8.5\n"
        "  2  Stress Period 2\n"
        "  3  9.0  9.0\n"
        "  4  8.5  8.0\n"
    )
    # The recharge array comes through OPEN/CLOSE with a multiplier of 2.0, which
    # is the record shape that once silently negated an entire array.
    (ws / "rech.dat").write_text("0.001 0.002 0.003 0.004\n")
    (ws / "tiny.rch").write_text(
        "# rch\n"
        "  1  0\n"
        "  1\n"
        "OPEN/CLOSE rech.dat 2.0 (FREE) -1  Stress Period 1\n"
        " -1\n"
    )
    (ws / "tiny.oc").write_text("HEAD SAVE UNIT 30\nPeriod 1 Step 1\nSAVE HEAD\n")

    nam = ws / "tiny.nam"
    nam.write_text(
        "LIST 7 tiny.lst\n"
        "BAS6 1 tiny.bas\n"
        "DISU 10 tiny.disu\n"
        "LPF 11 tiny.lpf\n"
        "SMS 19 tiny.sms\n"
        "OC 22 tiny.oc\n"
        "DRN 13 tiny.drn\n"
        "CHD 18 tiny.chd\n"
        "RCH 12 tiny.rch\n"
    )
    _write_gsf(ws / "tiny.gsf")
    return nam


@pytest.fixture(scope="module")
def usg(usg_workspace):
    """The synthetic model, read."""

    return read_usg(usg_workspace, gsf=usg_workspace.parent / "tiny.gsf", crs="EPSG:2927")


# --- the grid claim, checked rather than assumed --------------------------


def test_disu_is_recognized_as_layered(usg):
    """The DISU is layered, so DISV is a lossless target."""

    ok, why = layered_report(usg.disu)
    assert ok, why
    assert usg.nlay == NLAY and usg.ncpl == NCPL and usg.nodes == NODES


def test_non_layered_disu_is_refused():
    """A DISU whose connections are not layered must not convert silently."""

    from myflopy.modflow.usg.packages import DisuData

    iac, ja = _connectivity()
    ja = ja.copy()
    ja[1] = 8  # node 1 -> node 8: a diagonal jump no DISV can express
    broken = DisuData(
        nodes=NODES, nlay=NLAY, njag=int(iac.sum()), ivsd=-1, nper=2, itmuni=4, lenuni=1,
        nodelay=np.full(NLAY, NCPL), top=np.zeros(NODES), bot=np.zeros(NODES),
        area=np.zeros(NODES), iac=iac, ja=ja,
        perlen=np.ones(2), nstp=np.ones(2, int), tsmult=np.ones(2),
        steady=np.zeros(2, bool),
    )
    ok, why = layered_report(broken)
    assert not ok
    assert "vertical" in why


def test_grid_geometry_matches_the_disu_areas(usg):
    """The .gsf polygons reproduce the cell areas the DISU states independently."""

    assert usg.grid is not None
    areas = usg.grid.gdf_vorPolys.area.to_numpy()
    assert np.allclose(areas, 100.0)
    assert usg.grid.gdf_vorPolys.is_valid.all()


# --- node numbering -------------------------------------------------------


def test_node_to_cellid_round_trips(usg):
    """Every groundwater node maps to exactly one (layer, cell), and back."""

    nodes = np.arange(1, NODES + 1)
    cellid = usg.to_cellid(nodes)
    assert cellid.shape == (NODES, 2)
    assert (cellid[:, 0] * NCPL + cellid[:, 1] + 1 == nodes).all()
    assert len(set(map(tuple, cellid))) == NODES


def test_cln_nodes_map_to_no_cell(usg):
    """A node above the groundwater grid belongs to CLN and has no MF6 cell."""

    assert (usg.to_cellid(np.array([NODES + 1, NODES + 9])) == -1).all()


# --- the readers that broke ----------------------------------------------


def test_open_close_multiplier_is_applied_by_token_not_by_number(usg):
    """The OPEN/CLOSE multiplier is field 2, even though field 1 is a filename.

    Indexing the multiplier by position among the *numbers* on the line reads
    ``IPRN`` instead, which silently negates the whole array.
    """

    assert usg.rch is not None
    array = usg.rch.rech[0]
    assert np.allclose(array, np.array([0.001, 0.002, 0.003, 0.004]) * 2.0)
    assert (array > 0).all()


def test_itmp_negative_means_reuse_not_empty(usg):
    """A period that reused the previous one is absent, never present-and-empty.

    MODFLOW 6 reads an absent period as "carry on" and an EMPTY period as
    "delete every record", so the distinction is the whole boundary condition.
    """

    drn = usg.boundaries["DRN"]
    assert drn.n_defined == 1 and drn.reused == (1,)
    assert 1 not in drn.periods

    assert usg.rch.reused == (1,)
    assert 1 not in usg.rch.rech


def test_list_records_are_read_with_their_values(usg):
    """DRN and CHD records keep node, elevation/head and conductance."""

    drn = usg.boundaries["DRN"].periods[0]
    assert drn.shape == (2, 3)
    assert np.allclose(drn[:, 0], [1, 2])
    assert np.allclose(drn[:, 2], [100.0, 120.0])

    chd = usg.boundaries["CHD"]
    assert chd.n_defined == 2
    assert np.allclose(chd.periods[1][:, 1], [9.0, 8.5])


def test_bas_and_lpf_are_read_without_flopys_loader(usg):
    """BAS/LPF come from myflopy's own reader, which FloPy's cannot always do."""

    assert usg.idomain.shape == (NLAY, NCPL)
    assert (usg.idomain == 1).all()
    assert np.allclose(usg.strt, 8.0)
    assert np.allclose(usg.k, 10.0)
    assert np.allclose(usg.k33, 1.0)
    assert (usg.icelltype == 1).all()  # LAYTYP 4 is convertible


def test_geometry_arrays_have_disv_shapes(usg):
    """top is (ncpl,), botm/idomain are (nlay, ncpl) -- what MF6 DISV wants."""

    assert usg.top.shape == (NCPL,)
    assert usg.botm.shape == (NLAY, NCPL)
    assert np.allclose(usg.top, 10.0)
    assert np.allclose(usg.botm[0], 5.0) and np.allclose(usg.botm[1], 0.0)
    assert (usg.thickness == 5.0).all()


def test_uppermost_active_resolves_the_ets_target(usg):
    """NETSOP=3 needs the highest active cell per column; IBOUND is static."""

    assert usg.uppermost_active.shape == (NCPL,)
    assert (usg.uppermost_active == 0).all()
    assert usg.has_active_column.all()


def test_time_discretization(usg):
    """Periods, lengths and units survive."""

    assert usg.nper == 2
    assert usg.time_units == "days"
    assert np.allclose(usg.disu.perlen, 30.0)


# --- the conversion -------------------------------------------------------


def test_to_mf6_builds_a_runnable_spec(usg):
    """The converted spec carries the packages the USG model actually had."""

    sim = usg.to_mf6("tiny", start_date_time="2020-01-01")
    names = [p.name for p in sim.models[0].packages]
    assert names[:4] == ["disv", "npf", "ic", "sto"]
    # named `rch`, not `rcha`: MODFLOW 6 has ONE recharge package and READASARRAYS
    # is an option inside its file -- FloPy's two classes are a FloPy artifact.
    assert {"chd", "drn", "rch", "oc"} <= set(names)
    assert [p.name for p in sim.packages] == ["tdis", "ims"]


def test_conversion_writes_a_local_origin(usg):
    """The mesh is written near zero with its true corner declared.

    MODFLOW 6 builds DISV conductances from raw vertex coordinates and loses the
    precision to do so on projected coordinates, returning a NaN budget while
    still reporting normal termination.
    """

    sim = usg.to_mf6("tiny", start_date_time="2020-01-01")
    disv = next(p for p in sim.models[0].packages if p.name == "disv")
    assert disv.options["xorigin"] == 0.0  # this synthetic mesh already starts at 0
    xs = [v[1] for v in disv.options["vertices"]]
    assert min(xs) == 0.0


def test_transient_flags_are_not_lost(usg):
    """Every transient period is flagged; an empty period map is never passed.

    Passing ``steady_state={}`` makes FloPy emit empty period blocks, which drops
    the TRANSIENT keyword and makes MODFLOW 6 treat the run as steady state.
    """

    sim = usg.to_mf6("tiny", start_date_time="2020-01-01")
    sto = next(p for p in sim.models[0].packages if p.name == "sto")
    assert sto.options["transient"] == {0: True, 1: True}
    assert sto.options.get("steady_state") is None


def test_solver_defaults_to_complex(usg):
    """A converted USG model is nonlinear; MODERATE kills MODFLOW 6 on this class."""

    sim = usg.to_mf6("tiny", start_date_time="2020-01-01")
    ims = next(p for p in sim.packages if p.name == "ims")
    assert ims.options["complexity"] == "COMPLEX"


def test_report_names_what_it_did(usg):
    """The report accounts for the conversion in words a modeller can check."""

    text = usg.report()
    assert "DISU  -> DISV" in text
    assert "CONVERTED" in text
    assert "2 stress periods" in text


def test_conversion_without_a_grid_is_refused(usg_workspace):
    """A DISU has no coordinates, so converting without the .gsf must not proceed."""

    model = read_usg(usg_workspace, gsf=None, read_boundaries=False)
    if model.grid is None:
        with pytest.raises(ValueError, match="grid geometry"):
            model.to_mf6("tiny")


# --- the low-level reader -------------------------------------------------


def test_name_file_resolves_paths_case_insensitively(usg_workspace):
    """These files are written on Windows and read on Linux."""

    name_file = NameFile.read(usg_workspace)
    assert "DISU" in name_file.package_types
    assert name_file.package("DISU").path.is_file()


def test_array_control_records(tmp_path):
    """CONSTANT / INTERNAL / OPEN/CLOSE, including each one's multiplier field."""

    (tmp_path / "vals.dat").write_text("1 2 3 4\n")
    target = tmp_path / "arrays.txt"
    target.write_text(
        "CONSTANT 7.5\n"
        "INTERNAL 3.0 (FREE) -1  a label\n"
        "1 2 3 4\n"
        "OPEN/CLOSE vals.dat 10.0 (FREE) -1  another label\n"
    )
    cursor = ArrayCursor(target)
    assert np.allclose(cursor.read_array(4), 7.5)
    assert np.allclose(cursor.read_array(4), [3.0, 6.0, 9.0, 12.0])
    assert np.allclose(cursor.read_array(4), [10.0, 20.0, 30.0, 40.0])


def test_a_bad_control_record_says_so(tmp_path):
    """A file that is not where the reader thinks fails by name, not by silence."""

    target = tmp_path / "bad.txt"
    target.write_text("NOT_A_CONTROL_WORD 1 2 3\n")
    with pytest.raises(ValueError, match="array control record"):
        ArrayCursor(target).read_array(4)


def test_sms_under_relaxation_is_carried_into_ims(usg):
    """SMS's delta-bar-delta tuning maps field-for-field onto MODFLOW 6's IMS.

    This is the tuning that made the original model converge. Keeping only the
    tolerances and discarding it leaves MODFLOW 6 taking Newton steps of tens of
    thousands of feet on a model whose heads span a few hundred.
    """

    sms = usg.sms
    assert sms.is_newton
    assert sms.theta is not None and sms.numtrack is not None

    sim = usg.to_mf6("tiny", start_date_time="2020-01-01")
    ims = next(p for p in sim.packages if p.name == "ims").options
    assert ims["under_relaxation"] == "DBD"
    assert ims["under_relaxation_theta"] == sms.theta
    assert ims["under_relaxation_kappa"] == sms.akappa
    assert ims["under_relaxation_gamma"] == sms.gamma
    assert ims["backtracking_number"] == sms.numtrack
    assert ims["backtracking_tolerance"] == sms.btol


def test_sms_inner_tolerance_is_not_copied(usg):
    """SMS's HICLOSE is not MODFLOW 6's INNER_DVCLOSE; the preset owns the inner solve."""

    sim = usg.to_mf6("tiny", start_date_time="2020-01-01")
    ims = next(p for p in sim.packages if p.name == "ims").options
    assert "inner_dvclose" not in ims
    assert ims["outer_dvclose"] == usg.sms.hclose
    assert ims["outer_maximum"] == usg.sms.mxiter


def test_the_grid_gets_the_models_layer_elevations(usg):
    """A `.gsf` carries geometry only; the model carries the elevations.

    Publishing them as `grid.gdf_topbtm` is what lets the layer-elevation hover,
    the mounding colorscale and the surface-aware SFR/LAK builders work on an
    imported grid -- all of which read that one frame.
    """

    frame = usg.grid.gdf_topbtm
    assert frame is not None, "read_usg should attach layer elevations"
    # geometry, then 0 = model top and 1..nlay = layer bottoms
    assert list(frame.columns) == ["geometry", *range(NLAY + 1)]
    assert np.allclose(frame[0].to_numpy(), usg.top)
    for layer in range(NLAY):
        assert np.allclose(frame[layer + 1].to_numpy(), usg.botm[layer])
    assert frame.crs == usg.grid.crs
