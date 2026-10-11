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
import shapely

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


#: The hand-written list/array packages. These carry the reuse conventions the importer
#: has to honour, which FloPy's writers do not expose. A variant model overrides one by
#: name -- the below-cell-bottom tests build the same mesh with offending boundaries on
#: it rather than a second mesh of their own.
_DEFAULT_PACKAGES = {
    "DRN": (13, "tiny.drn", (
        "# drn\n"
        "  2  0\n"
        "  2  0  Stress Period 1\n"
        "  1  7.5  100.0\n"
        "  2  7.5  120.0\n"
        " -1  0  Stress Period 2 (reuse)\n"
    )),
    "CHD": (18, "tiny.chd", (
        "# chd\n"
        "  2\n"
        "  2  Stress Period 1\n"
        "  3  9.0  9.0\n"
        "  4  9.0  8.5\n"
        "  2  Stress Period 2\n"
        "  3  9.0  9.0\n"
        "  4  8.5  8.0\n"
    )),
    # The recharge array comes through OPEN/CLOSE with a multiplier of 2.0, which
    # is the record shape that once silently negated an entire array.
    "RCH": (12, "tiny.rch", (
        "# rch\n"
        "  1  0\n"
        "  1\n"
        "OPEN/CLOSE rech.dat 2.0 (FREE) -1  Stress Period 1\n"
        " -1\n"
    )),
}


def _write_usg_model(ws, *, laytyp=4, packages=None):
    """Write a complete little MODFLOW-USG model into ``ws`` and return its name file.

    ``laytyp`` becomes MODFLOW 6's ``icelltype``, which is the whole severity split for
    a boundary elevation below its cell bottom: 4 is convertible, 0 confined.
    """

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
    flopy.mfusg.MfUsgLpf(
        model, laytyp=laytyp, layvka=0, chani=1.0, hk=10.0, vka=1.0, ss=1e-5, sy=0.2
    )
    flopy.mfusg.MfUsgSms(model, hclose=0.01, hiclose=0.001, mxiter=50, iter1=100, nonlinmeth=1)
    model.write_input()

    written = dict(_DEFAULT_PACKAGES)
    written.update(packages or {})
    (ws / "rech.dat").write_text("0.001 0.002 0.003 0.004\n")
    for _, filename, text in written.values():
        (ws / filename).write_text(text)
    (ws / "tiny.oc").write_text("HEAD SAVE UNIT 30\nPeriod 1 Step 1\nSAVE HEAD\n")

    nam = ws / "tiny.nam"
    nam.write_text(
        "LIST 7 tiny.lst\n"
        "BAS6 1 tiny.bas\n"
        "DISU 10 tiny.disu\n"
        "LPF 11 tiny.lpf\n"
        "SMS 19 tiny.sms\n"
        "OC 22 tiny.oc\n"
        + "".join(
            f"{ftype} {unit} {filename}\n"
            for ftype, (unit, filename, _) in written.items()
        )
    )
    _write_gsf(ws / "tiny.gsf")
    return nam


@pytest.fixture(scope="module")
def usg_workspace(tmp_path_factory):
    """The baseline synthetic model: every boundary inside the cell it sits in."""

    return _write_usg_model(tmp_path_factory.mktemp("usg"))


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


def test_the_converted_npf_saves_what_particle_tracking_reads(usg):
    """PRT reads specific discharge and saturation from the flow budget; without
    both, every particle-tracking run on a converted model failed (ledger 217)."""

    sim = usg.to_mf6("tiny", start_date_time="2020-01-01")
    npf = next(p for p in sim.models[0].packages if p.name == "npf")
    assert npf.options["save_specific_discharge"] is True
    assert npf.options["save_saturation"] is True


# --- boundaries below the cell bottom (ledger 291) -------------------------

# A USG RIV record is `node, stage, cond, rbot`. The two rows are the two branches
# of the repair rule. Node 1 is the case a stage-only check cannot see: its BED is
# below the layer-0 bottom of 5.0 while its water surface is comfortably above it,
# which is what an ordinary river looks like -- and its stage must not be touched.
# Node 2 has BOTH below the bottom, so its stage has to move, and the declared 2 ft
# of depth has to survive the move.
_BELOW_BOTTOM = {
    "RIV": (14, "tiny.riv", (
        "# riv\n"
        "  2  0\n"
        "  2  0  Stress Period 1\n"
        "  1  7.0  50.0  2.0\n"
        "  2  3.0  60.0  1.0\n"
        " -1  0  Stress Period 2 (reuse)\n"
    )),
    # Layer 2 is bottomed at exactly 0.0, where a relative nudge degenerates to
    # the absolute one.
    "GHB": (15, "tiny.ghb", (
        "# ghb\n"
        "  1  0\n"
        "  1  0  Stress Period 1\n"
        "  5  -1.0  70.0\n"
        " -1  0  Stress Period 2 (reuse)\n"
    )),
    "DRN": (13, "tiny.drn", (
        "# drn\n"
        "  2  0\n"
        "  2  0  Stress Period 1\n"
        "  1  7.5  100.0\n"
        "  2  2.5  120.0\n"
        " -1  0  Stress Period 2 (reuse)\n"
    )),
}

#: Where the offending records sit, so a test can state the expected answer.
LAYER0_BOTTOM, LAYER1_BOTTOM, MODEL_TOP = 5.0, 0.0, 10.0


@pytest.fixture(scope="module")
def below_bottom(tmp_path_factory):
    """Factory for the same mesh with boundaries hanging below their own cells.

    ``laytyp`` is the whole point: it becomes MODFLOW 6's ``icelltype``, and that is
    what decides whether MODFLOW 6 refuses the model or runs it and quietly changes
    what the boundary does.
    """

    cache = {}

    def build(laytyp=4):
        if laytyp not in cache:
            ws = tmp_path_factory.mktemp(f"usg_below_{laytyp}")
            nam = _write_usg_model(ws, laytyp=laytyp, packages=_BELOW_BOTTOM)
            cache[laytyp] = read_usg(nam, gsf=ws / "tiny.gsf", crs="EPSG:2927")
        return cache[laytyp]

    return build


def _records(sim, name):
    """The MODFLOW 6 stress-period records the converted spec carries for a package."""

    package = next(p for p in sim.models[0].packages if p.name == name)
    return package.options["stress_period_data"][0]


def _finding(model, code):
    """One finding by code, or None."""

    return next((f for f in model.validate() if f.code == code), None)


# --- the field MODFLOW 6 actually checks -----------------------------------


def test_a_river_bed_below_its_cell_is_seen_even_when_the_stage_is_not(below_bottom):
    """The ordinary river: bed below the cell bottom, water surface above it.

    MODFLOW 6 measures ``rbot`` against the cell bottom. Inspecting ``stage`` instead
    -- which is the FIRST value in the record, so the one a positional read lands on
    -- reports nothing at all for this record, and the run dies on a field validate()
    never looked at.
    """

    model = below_bottom(laytyp=4)
    riv = model.boundaries["RIV"].periods[0]
    assert riv[0, 1] > LAYER0_BOTTOM > riv[0, 3], "node 1 must offend on rbot ALONE"

    finding = _finding(model, "riv_rbot_below_cell_bottom")
    assert finding is not None, "a river bed below its cell bottom must be named"
    assert finding.blocks_run
    assert finding.count == 2
    assert "rbot" in finding.detail


def test_fix_for_mf6_lifts_the_river_BED_not_its_stage(below_bottom):
    """``rbot`` is what MODFLOW 6 compares; lifting ``stage`` leaves the run refused.

    Measured on a real MODFLOW 6 run at ``icelltype = 1``, cell bottom 48.0: stage
    lifted with rbot left at 35.0 still fails with ``RIV BOUNDARY (1) RIVER BOTTOM
    IS LESS THAN CELL BOTTOM``.
    """

    model = below_bottom(laytyp=4)
    sim = model.to_mf6("tiny", start_date_time="2020-01-01", fix_for_mf6=True)
    stages = np.array([float(r[1]) for r in _records(sim, "riv")])
    rbots = np.array([float(r[3]) for r in _records(sim, "riv")])

    assert (rbots > LAYER0_BOTTOM).all(), "every river bed must end up inside its cell"
    # Strictly above, never equal: FloPy writes botm fixed-point and a boundary value
    # in E format, so exactly-equal can round BELOW in the file MODFLOW 6 reads.
    assert (rbots - LAYER0_BOTTOM).min() > 0
    assert (stages > rbots).all(), "MODFLOW 6 enforces stage >= rbot as well"


def test_a_river_lifted_off_the_bottom_keeps_the_depth_it_was_drawn_with(below_bottom):
    """Node 2, whose stage was below the bottom too, so the stage HAD to move.

    Clamping both numbers to the floor collapses the declared depth to zero: measured
    on a one-cell bench at cell bottom 48.0, both-clamped q = -381.82 with 0.000 ft of
    depth against pair-lifted -109.09.
    """

    model = below_bottom(laytyp=4)
    declared = model.boundaries["RIV"].periods[0]
    both_below = next(i for i, row in enumerate(declared) if row[1] < LAYER0_BOTTOM)

    sim = model.to_mf6("tiny", start_date_time="2020-01-01", fix_for_mf6=True)
    record = _records(sim, "riv")[both_below]

    assert float(record[3]) > float(declared[both_below, 3]), "the bed must have moved"
    assert float(record[1]) > float(declared[both_below, 1]), "and carried its stage"
    depth = float(record[1]) - float(record[3])
    assert depth == pytest.approx(float(declared[both_below, 1] - declared[both_below, 3]))


def test_a_river_whose_water_surface_is_legal_keeps_it_exactly(below_bottom):
    """Node 1: bed below the cell bottom, water surface inside the cell.

    A stage already inside its cell is the water surface the river imposes, and moving
    it moves the answer. Lifting it by the bed's delta -- which is what ledger 290's
    rule did, generalised from a case where the stage was ALSO below -- measured on a
    one-cell bench (bottom 48, top 60): stage 55 / rbot 40 became stage 63, three feet
    ABOVE land surface, and q went from 0.00 to +436.37. With a deeper bed the sign
    reverses: -272.73 to +1254.55. Here it would put the stage above this model's own
    top of 10.0.
    """

    model = below_bottom(laytyp=4)
    declared = model.boundaries["RIV"].periods[0]
    legal = next(i for i, row in enumerate(declared) if row[1] > LAYER0_BOTTOM > row[3])

    sim = model.to_mf6("tiny", start_date_time="2020-01-01", fix_for_mf6=True)
    record = _records(sim, "riv")[legal]

    assert float(record[1]) == float(declared[legal, 1]), "a legal stage must not move"
    assert float(record[3]) > LAYER0_BOTTOM, "while its bed still comes inside the cell"
    assert float(record[1]) <= MODEL_TOP, "and no stage may end up above land surface"


def test_an_inverted_river_is_named_rather_than_pair_lifted(tmp_path):
    """A stage below its own bed is refused by MODFLOW 6 at any icelltype.

    Carrying the stage up with the bed preserves the inversion rather than repairing
    it, and which of the two numbers is wrong is not something a converter can know.
    """

    inverted = dict(_BELOW_BOTTOM)
    inverted["RIV"] = (14, "tiny.riv", (
        "# riv\n"
        "  1  0\n"
        "  1  0  Stress Period 1\n"
        "  1  6.0  50.0  8.0\n"
        " -1  0  Stress Period 2 (reuse)\n"
    ))
    nam = _write_usg_model(tmp_path, laytyp=4, packages=inverted)
    model = read_usg(nam, gsf=tmp_path / "tiny.gsf", crs="EPSG:2927")

    finding = _finding(model, "riv_stage_below_river_bottom")
    assert finding is not None and finding.blocks_run
    assert "will not guess" in finding.fix


def test_a_short_river_record_does_not_silently_empty_its_period(tmp_path):
    """A malformed record must not be floored past its own end -- nor vanish.

    A RIV line carrying three numbers instead of four gives a values block with no
    `rbot` column at all. Indexing it raises; skipping the period instead deletes
    every record in it, and in MODFLOW 6 an ABSENT period means "reuse the previous
    one" while an EMPTY one means "delete every record" -- so a quiet skip here is a
    boundary condition silently switched off.
    """

    short = dict(_BELOW_BOTTOM)
    short["RIV"] = (14, "tiny.riv", (
        "# riv\n"
        "  1  0\n"
        "  1  0  Stress Period 1\n"
        "  1  7.0  50.0\n"
        " -1  0  Stress Period 2 (reuse)\n"
    ))
    nam = _write_usg_model(tmp_path, laytyp=4, packages=short)
    model = read_usg(nam, gsf=tmp_path / "tiny.gsf", crs="EPSG:2927")

    sim = model.to_mf6("tiny", start_date_time="2020-01-01", fix_for_mf6=True)
    records = _records(sim, "riv")
    assert len(records) == 1
    assert float(records[0][1]) == pytest.approx(7.0), "left exactly as it was read"
    # The record is still malformed, and FloPy refuses it at build time rather than
    # here -- with a message that names neither the RIV line nor its period. The
    # conversion's job is only not to make it disappear first.


# --- the DRN exemption, which was false in both directions -----------------


def test_a_drain_below_its_cell_bottom_stops_the_run_and_validate_says_so(below_bottom):
    """It was exempted on the premise that such a drain is simply inactive.

    Measured: at ``icelltype != 0`` MODFLOW 6 refuses by name -- ``DRN BOUNDARY (1)
    ELEVATION ( 45.000) IS LESS THAN CELL BOTTOM`` -- so validate(), whose documented
    job is to name what MODFLOW 6 will reject, was silent about a model that cannot
    start.
    """

    model = below_bottom(laytyp=4)
    finding = _finding(model, "drn_elev_below_cell_bottom")
    assert finding is not None, "a drain below its cell bottom must be named"
    assert finding.blocks_run
    assert finding.count == 1
    assert any(f.blocks_run for f in model.validate())


def test_a_confined_model_is_told_it_will_run_and_over_discharge(below_bottom):
    """The other half of the same false premise: at ``icelltype = 0`` it is not inactive.

    MODFLOW 6 accepts the record, terminates normally, and drains against a head below
    the cell -- measured 2.0x the discharge of the same drain placed inside its cell.
    So the finding is reported but must NOT block the run.
    """

    model = below_bottom(laytyp=0)
    finding = _finding(model, "drn_elev_below_cell_bottom")
    assert finding is not None
    assert not finding.blocks_run and finding.severity == "allow"
    assert "CONFINED" in finding.detail
    assert not [f for f in model.validate() if f.blocks_run]


def test_fix_for_mf6_lifts_a_drain_invert_too(below_bottom):
    """Otherwise validate() names a blocker whose stated remedy does nothing."""

    model = below_bottom(laytyp=4)
    sim = model.to_mf6("tiny", start_date_time="2020-01-01", fix_for_mf6=True)
    elevations = np.array([float(r[1]) for r in _records(sim, "drn")])

    assert (elevations >= LAYER0_BOTTOM).all()
    assert (elevations - LAYER0_BOTTOM).min() > 0
    # The drain that was already inside its cell is left exactly where it was drawn.
    assert elevations.max() == pytest.approx(7.5)


def test_a_bottom_of_exactly_zero_still_lands_strictly_above_it(below_bottom):
    """Where the relative nudge degenerates, the absolute one has to carry it."""

    model = below_bottom(laytyp=4)
    sim = model.to_mf6("tiny", start_date_time="2020-01-01", fix_for_mf6=True)
    heads = np.array([float(r[1]) for r in _records(sim, "ghb")])
    assert (heads > LAYER1_BOTTOM).all()


# --- one rule, not two copies of it ----------------------------------------


def test_the_importer_and_the_geopackage_sources_floor_by_the_same_rule():
    """Both layers used to carry their own copy, and the copies disagreed.

    CLAUDE.md's rule is to grep BOTH the package-first layer and the legacy engine
    before adding a capability. The USG importer's copy clamped the wrong RIV field
    for as long as it existed.
    """

    import myflopy.geopackage as geopackage
    from myflopy.modflow.usg import convert

    assert convert.river_pair_lift is geopackage.river_pair_lift
    assert convert.bottom_clearance is geopackage.bottom_clearance


# --- what MODFLOW 6 itself says --------------------------------------------


def _verdict(model, workspace, *, include, fix_for_mf6):
    """Run the converted model under MODFLOW 6 and return ``(success, listing text)``."""

    sim = model.to_mf6(
        "tiny", start_date_time="2020-01-01", include=include, fix_for_mf6=fix_for_mf6
    )
    success, report = sim.build(workspace).execute()
    return success, "\n".join(report).upper()


@pytest.mark.slow
def test_mf6_refuses_the_model_until_the_river_BED_is_inside_its_cell(
    below_bottom, tmp_path
):
    """The end of the argument: MODFLOW 6's own verdict on both spellings.

    This is the test that fails if the clamp goes back to lifting ``stage`` -- which
    it did, and which left ``rbot`` exactly where MODFLOW 6 was about to look.
    """

    model = below_bottom(laytyp=4)
    include = ("chd", "riv")

    success, listing = _verdict(model, tmp_path / "raw", include=include, fix_for_mf6=False)
    assert not success
    assert "RIVER BOTTOM" in listing and "LESS THAN CELL BOTTOM" in listing

    success, listing = _verdict(model, tmp_path / "fixed", include=include, fix_for_mf6=True)
    assert success, listing[-2000:]


@pytest.mark.slow
def test_mf6_refuses_the_model_until_the_drain_invert_is_inside_its_cell(
    below_bottom, tmp_path
):
    """The drain validate() used to exempt, put to MODFLOW 6 itself.

    Each package is read in turn and MODFLOW 6 stops at the FIRST one it objects to,
    which is the argument for naming every case up front instead of discovering them
    one failed run at a time.
    """

    model = below_bottom(laytyp=4)
    include = ("chd", "drn")

    success, listing = _verdict(model, tmp_path / "raw", include=include, fix_for_mf6=False)
    assert not success
    assert "DRN BOUNDARY" in listing and "LESS THAN CELL BOTTOM" in listing

    success, listing = _verdict(model, tmp_path / "fixed", include=include, fix_for_mf6=True)
    assert success, listing[-2000:]


@pytest.mark.slow
def test_a_confined_model_runs_with_the_very_same_records(below_bottom, tmp_path):
    """The claim behind the non-blocking severity, checked rather than asserted."""

    model = below_bottom(laytyp=0)
    success, _ = _verdict(
        model, tmp_path / "confined", include=("chd", "drn", "ghb", "riv"), fix_for_mf6=False
    )
    assert success, "at icelltype 0 MODFLOW 6 accepts every one of these records"


def test_a_model_whose_boundaries_are_inside_their_cells_says_nothing(usg):
    """The other direction: a check that fires on everything reports nothing either.

    The baseline model's DRN sits at 7.5 in a cell bottomed at 5.0 and its CHD at 8.0
    to 9.0, so every finding here would be a false one -- and a false `reject` sends a
    modeller editing boundary elevations that were right.
    """

    assert usg.validate() == []
    assert "MODFLOW 6 DIFFERS" not in usg.report()


def test_a_confined_model_is_warned_about_at_conversion_time(below_bottom, caplog):
    """The `allow` case is the one nothing else says out loud.

    MODFLOW 6 runs the model and the boundary does something other than what was drawn;
    leaving that to `validate()`/`report()`, which a caller has to think to ask for, is
    ledger 288's defect -- a degradation reported through a channel nobody hears.
    """

    model = below_bottom(laytyp=0)
    with caplog.at_level("WARNING", logger="myflopy.modflow.usg.convert"):
        model.to_mf6("tiny", start_date_time="2020-01-01")

    said = [r.getMessage() for r in caplog.records if "change what a boundary does" in r.getMessage()]
    assert said, "a confined model with below-bottom boundaries must not convert in silence"
    assert any("DRN" in message for message in said)
    assert not [r for r in caplog.records if "will reject this model" in r.getMessage()]


# --- a boundary switched off stays off (ledger 367) -------------------
#
# MODFLOW 6 reads an ABSENT period as "reuse the previous one" and an EMPTY period
# as "delete every record". A USG period with ITMP = 0 switched the package off, so
# it has to become the second, never the first.

_SWITCHED_OFF = {
    "DRN": (13, "tiny.drn", (
        "# drn\n"
        "  2  0\n"
        "  2  0  Stress Period 1\n"
        "  1  7.5  100.0\n"
        "  2  7.5  120.0\n"
        "  0  0  Stress Period 2 (switched off)\n"
    )),
}


def _period_keys(sim, name):
    package = next(p for p in sim.models[0].packages if p.name == name)
    return package.options["stress_period_data"]


def test_an_itmp_zero_period_is_written_empty_not_left_out(tmp_path):
    nam = _write_usg_model(tmp_path, packages=_SWITCHED_OFF)
    model = read_usg(nam, gsf=tmp_path / "tiny.gsf", crs="EPSG:2927")
    data = _period_keys(model.to_mf6("tiny", start_date_time="2020-01-01"), "drn")

    assert sorted(data) == [0, 1]
    assert len(data[0]) == 2
    assert data[1] == [], "ITMP=0 switches the drains off; an absent period would reuse them"
    assert "1 with no groundwater record (switched off)" in model.report()


def test_a_reused_period_is_still_left_out(usg):
    """The other half: ITMP < 0 is reuse, and reuse is an ABSENT period."""

    data = _period_keys(usg.to_mf6("tiny", start_date_time="2020-01-01"), "drn")
    assert sorted(data) == [0]


def test_a_period_whose_records_are_all_cln_is_written_empty(tmp_path):
    """Every record on a CLN node: the GROUNDWATER drains are off that period."""

    packages = {"DRN": (13, "tiny.drn", (
        "# drn\n"
        "  2  0\n"
        "  2  0  Stress Period 1\n"
        "  1  7.5  100.0\n"
        "  2  7.5  120.0\n"
        "  1  0  Stress Period 2 (CLN only)\n"
        f"  {NODES + 1}  7.5  100.0\n"
    ))}
    nam = _write_usg_model(tmp_path, packages=packages)
    model = read_usg(nam, gsf=tmp_path / "tiny.gsf", crs="EPSG:2927")
    assert model.boundaries["DRN"].n_cln_records == 1
    data = _period_keys(model.to_mf6("tiny", start_date_time="2020-01-01"), "drn")
    assert sorted(data) == [0, 1] and data[1] == []


def test_a_period_fix_for_mf6_emptied_is_written_empty(tmp_path):
    """Every CHD record of period 2 sits below its cell bottom and is omitted.

    Leaving the period out instead re-applied period 1's constant heads through it.
    """

    packages = {"CHD": (18, "tiny.chd", (
        "# chd\n"
        "  2\n"
        "  2  Stress Period 1\n"
        "  3  9.0  9.0\n"
        "  4  9.0  8.5\n"
        "  2  Stress Period 2\n"
        "  3  4.0  4.0\n"
        "  4  4.0  4.0\n"
    ))}
    nam = _write_usg_model(tmp_path, packages=packages)
    model = read_usg(nam, gsf=tmp_path / "tiny.gsf", crs="EPSG:2927")
    data = _period_keys(
        model.to_mf6("tiny", start_date_time="2020-01-01", fix_for_mf6=True), "chd"
    )
    assert sorted(data) == [0, 1] and len(data[0]) == 2 and data[1] == []


def test_a_package_with_no_groundwater_record_anywhere_is_not_written(tmp_path):
    """Empty periods are kept only beside a period that HAS records."""

    packages = {"DRN": (13, "tiny.drn", (
        "# drn\n"
        "  1  0\n"
        "  0  0  Stress Period 1\n"
        "  0  0  Stress Period 2\n"
    ))}
    nam = _write_usg_model(tmp_path, packages=packages)
    model = read_usg(nam, gsf=tmp_path / "tiny.gsf", crs="EPSG:2927")
    sim = model.to_mf6("tiny", start_date_time="2020-01-01")
    assert "drn" not in [p.name for p in sim.models[0].packages]


@pytest.mark.slow
def test_mf6_moves_no_drain_water_in_a_switched_off_period(tmp_path):
    """MODFLOW 6's own budget: the drains discharge in period 1 and not in period 2."""

    (tmp_path / "usg").mkdir()
    nam = _write_usg_model(tmp_path / "usg", packages=_SWITCHED_OFF)
    model = read_usg(nam, gsf=tmp_path / "usg" / "tiny.gsf", crs="EPSG:2927")
    sim = model.to_mf6("tiny", start_date_time="2020-01-01", include=("chd", "drn"))
    workspace = tmp_path / "mf6"
    success, report = sim.build(workspace).execute()
    assert success, "\n".join(report)[-2000:]

    budget = flopy.utils.CellBudgetFile(str(workspace / "tiny.cbc"))
    flows = {
        int(kper): float(budget.get_data(text="DRN", kstpkper=(kstp, kper))[0]["q"].sum())
        for kstp, kper in budget.get_kstpkper()
    }
    assert flows[0] < -1.0, "the drains discharge while they are on"
    assert flows[1] == 0.0, "ITMP=0 switched them off; reusing period 1 drained ~-64 ft3/d"


@pytest.mark.parametrize("barriers", [0, 1])
def test_read_hfb_is_n_by_3_even_with_no_barrier(tmp_path, barriers):
    """An HFB6 file declaring zero barriers reads as (0, 3), like a model with none."""

    lines = f"  0  0  {barriers}\n" + "".join("  1  2  0.001\n" for _ in range(barriers)) + "  0\n"
    nam = _write_usg_model(tmp_path, packages={"HFB6": (16, "tiny.hfb", lines)})
    model = read_usg(nam, gsf=tmp_path / "tiny.gsf", crs="EPSG:2927")
    assert model.hfb.shape == (barriers, 3)


# --- one CRS per grid (ledger 368) ------------------------------------


def test_to_mf6_crs_relabels_every_frame_the_grid_caches(usg_workspace):
    """``to_mf6(crs=)`` is a RELABEL, and the grid then reports ONE CRS.

    Assigning only ``grid.crs`` left ``gdf_vorPolys``/``gdf_topbtm`` (built by
    ``read_usg``) and the lat/lon frames on the old label, and readers split between
    the two: the GIS export read the frames, the heads/run readers the attribute.
    """

    model = read_usg(usg_workspace, gsf=usg_workspace.parent / "tiny.gsf", crs="EPSG:2927")
    grid = model.grid
    topbtm = grid.gdf_topbtm
    xy = grid.gdf_vorPolys.geometry.get_coordinates().to_numpy().copy()
    old_latlon = grid.gdf_latlon.total_bounds.copy()

    model.to_mf6("tiny", start_date_time="2020-01-01", crs="EPSG:2926")

    assert grid.crs == "EPSG:2926"
    assert grid.gdf_vorPolys.crs == "EPSG:2926"
    assert grid.gdf_topbtm.crs == "EPSG:2926"
    assert grid.gdf_topbtm is topbtm, "relabelled in place, so a held reference agrees"
    assert np.array_equal(grid.gdf_vorPolys.geometry.get_coordinates().to_numpy(), xy), (
        "a relabel, never a reprojection"
    )
    expected = grid.gdf_vorPolys.to_crs("EPSG:4326").total_bounds
    assert np.allclose(grid.gdf_latlon.total_bounds, expected)
    assert not np.allclose(grid.gdf_latlon.total_bounds, old_latlon), (
        "the lat/lon frame is a reprojection, so it must be rebuilt from the new label"
    )
    lon, lat = grid.points_to_latlon(np.array([10.0]), np.array([10.0]))
    assert grid.gdf_latlon.union_all().buffer(1e-9).contains(
        shapely.Point(lon[0], lat[0])
    ), "points and cells are placed from the same label"


def test_a_grid_pickled_while_crs_was_an_attribute_still_reads_its_crs(usg):
    """Old pickles carry ``crs`` in ``__dict__``, which the property would shadow."""

    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus

    state = {k: v for k, v in usg.grid.__getstate__().items() if k != "_crs"}
    state["crs"] = "EPSG:2927"
    old = VoronoiGridPlus.__new__(VoronoiGridPlus)
    old.__setstate__(state)
    assert old.crs == usg.grid.crs == "EPSG:2927"
    assert "crs" not in old.__dict__


# --- review of U1 (ledger 370) ----------------------------------------


def test_a_period_whose_records_were_not_recognised_is_refused(tmp_path):
    """ITMP = 2 with records written as integers: not read as records at all.

    The reader recognises a record by a REAL value after its node, so ``1 7 100``
    was taken for the next period's header and the period read EMPTY -- which,
    since ledger 367, converts as the drains switched off, in silence.
    """

    packages = {"DRN": (13, "tiny.drn", (
        "# drn\n"
        "  2  0\n"
        "  2  0  Stress Period 1\n"
        "  1  7.5  100.0\n"
        "  2  7.5  120.0\n"
        "  2  0  Stress Period 2\n"
        "  1  7  100\n"
        "  2  7  120\n"
    ))}
    nam = _write_usg_model(tmp_path, packages=packages)
    with pytest.raises(ValueError, match=r"stress period 2 declares ITMP = 2 record\(s\) but 0"):
        read_usg(nam, gsf=tmp_path / "tiny.gsf", crs="EPSG:2927")


def test_a_period_short_of_its_itmp_is_refused(tmp_path):
    """One real record and one integer one: one read of the two declared."""

    packages = {"DRN": (13, "tiny.drn", (
        "# drn\n"
        "  2  0\n"
        "  2  0  Stress Period 1\n"
        "  1  7.5  100.0\n"
        "  2  7  120\n"
        "  -1  0  Stress Period 2\n"
    ))}
    nam = _write_usg_model(tmp_path, packages=packages)
    with pytest.raises(ValueError, match=r"stress period 1 declares ITMP = 2 record\(s\) but 1"):
        read_usg(nam, gsf=tmp_path / "tiny.gsf", crs="EPSG:2927")


def test_a_grid_crs_of_none_is_refused_and_changes_nothing(usg_workspace):
    """``grid.crs = None`` raised from geopandas HALFWAY through the relabel."""

    model = read_usg(usg_workspace, gsf=usg_workspace.parent / "tiny.gsf", crs="EPSG:2927")
    grid = model.grid
    grid.gdf_topbtm
    centroid = grid.grid_centroid
    for bad in (None, "EPSG:not-a-code"):
        with pytest.raises(ValueError, match="grid.crs|cannot be None"):
            grid.crs = bad
        assert grid.crs == "EPSG:2927"
        assert grid.gdf_vorPolys.crs == "EPSG:2927"
        assert grid.gdf_topbtm.crs == "EPSG:2927"
        assert grid.grid_centroid.equals(centroid)


def test_a_numpy_integer_crs_is_an_epsg_code(usg_workspace):
    """``numpy.int64`` is not an ``int``: it was stored as itself, not ``EPSG:N``."""

    model = read_usg(usg_workspace, gsf=usg_workspace.parent / "tiny.gsf", crs="EPSG:2927")
    model.grid.crs = np.int64(2926)
    assert type(model.grid.crs) is str and model.grid.crs == "EPSG:2926"
    assert model.grid.gdf_vorPolys.crs == "EPSG:2926"


# --- include= names packages, and only packages (ledger 429) ------------


def _package_names(sim):
    return [p.name for p in sim.models[0].packages]


def test_a_bare_string_include_is_one_package_name(usg):
    """``include="chd"`` was iterated letter by letter and kept no package at all."""

    names = _package_names(usg.to_mf6("tiny", start_date_time="2020-01-01", include="CHD"))
    assert names == ["disv", "npf", "ic", "sto", "chd", "oc"]


@pytest.mark.parametrize("include", [("chd", "dnr"), "rech", ("disv",)])
def test_an_unknown_include_name_is_refused_naming_the_valid_ones(below_bottom, include, caplog):
    """A misspelt ``"dnr"`` dropped the drains in silence; it is refused before anything.

    ``below_bottom`` logs a finding on every conversion, so an empty log shows the
    refusal came first.
    """

    with caplog.at_level("WARNING", logger="myflopy.modflow.usg.convert"):
        with pytest.raises(ValueError, match=r"include= chooses among 'chd', 'drn'.*'hfb'") as raised:
            below_bottom().to_mf6("tiny", start_date_time="2020-01-01", include=include)
    bad = include if isinstance(include, str) else include[-1]
    assert repr(bad) in str(raised.value)
    assert "always written" in str(raised.value)
    assert not caplog.records, "refused before any conversion work or finding is logged"


def test_an_include_item_that_is_not_a_string_is_refused(usg):
    with pytest.raises(TypeError, match=r"include=\('chd', 3\)\): 3 is not a package name"):
        usg.to_mf6("tiny", start_date_time="2020-01-01", include=("chd", 3))


def test_include_names_a_package_the_model_lacks_without_complaint(usg):
    """``"ghb"`` is a valid name; the synthetic model just has no GHB to write."""

    names = _package_names(
        usg.to_mf6("tiny", start_date_time="2020-01-01", include=("ghb", "drn"))
    )
    assert names == ["disv", "npf", "ic", "sto", "drn", "oc"]
