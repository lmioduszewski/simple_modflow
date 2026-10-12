"""Build-side cellid checks read a DIS ``(layer, row, col)`` and a DISU ``(node,)`` as
their own grids (ledger 469).

``ModelSpec.build``'s range check, removed-cell check, observation check and LAK
connection check read a DISV package only, so on a DIS or DISU model nothing was
compared: a cellid past the grid or on a removed cell built clean, a DISV pair on a DIS
stopped inside FloPy's writer (``IndexError``), and a DIS triple on a DISV -- or a pair on
a DISU -- built clean and FloPy wrote only its leading parts, which MF6 RAN
(``docs/reproducers/build_cellids_dis_disu.py``). Every check now reads the built
discretization through ``CellNumbering`` (ledger 419).
"""

from __future__ import annotations

import numpy as np
import pytest
from flopy.utils.gridutil import get_disu_kwargs, get_disv_kwargs

import myflopy as mf

NLAY, NROW, NCOL = 2, 3, 4


def _discretization(kind, idomain=None):
    if kind == "dis":
        return mf.dis(nlay=NLAY, nrow=NROW, ncol=NCOL, delr=10.0, delc=10.0, top=10.0,
                      botm=[5.0, 0.0], idomain=idomain)
    if kind == "disv":
        return mf.disv(**get_disv_kwargs(NLAY, NROW, NCOL, 10.0, 10.0, 10.0, [5.0, 0.0]),
                       idomain=idomain)
    kwargs = get_disu_kwargs(NLAY, NROW, NCOL, np.full(NCOL, 10.0), np.full(NROW, 10.0), 10.0,
                             [5.0, 0.0], return_vertices=True)
    return mf.disu(**kwargs, idomain=idomain)


def _build(tmp_path, kind, *packages, idomain=None):
    model = mf.gwf("m", packages=[
        _discretization(kind, idomain), mf.npf(k=1.0), mf.ic(strt=8.5), mf.oc(), *packages
    ])
    return mf.simulation(model, name="s").build(tmp_path / "run")


def _chd(*cellids):
    return mf.chd(stress_period_data={0: [[cellid, 9.0] for cellid in cellids]})


@pytest.mark.parametrize(
    ("kind", "cellid", "form"),
    [("dis", (0, 5), r"\(layer, row, col\)"), ("disv", (0, 1, 2), r"\(layer, cell\)"),
     ("disu", (3, 5), r"\(node,\)"), ("dis", (0,), r"\(layer, row, col\)")],
)
def test_a_cellid_of_another_form_is_refused_before_flopy(tmp_path, kind, cellid, form):
    with pytest.raises(ValueError, match=rf"not this {kind.upper()}'s {form}: \({cellid[0]}"):
        _build(tmp_path, kind, _chd(cellid))


@pytest.mark.parametrize(
    ("kind", "cellid", "extent"),
    [("dis", (0, 0, 9), "nlay=2, nrow=3, ncol=4"), ("dis", (2, 0, 0), "nlay=2, nrow=3, ncol=4"),
     ("dis", (0, 3, 0), "nlay=2, nrow=3, ncol=4"), ("disu", (24,), "nodes=24"),
     ("disv", (0, 12), "nlay=2, ncpl=12")],
)
def test_a_cellid_off_the_grid_is_refused(tmp_path, kind, cellid, extent):
    with pytest.raises(ValueError, match=rf"grid does not have \({extent}\): \({cellid[0]}"):
        _build(tmp_path, kind, _chd(cellid))


@pytest.mark.parametrize(
    ("kind", "cellid", "shape", "removed"),
    [("dis", (0, 0, 1), (NLAY, NROW, NCOL), (0, 0, 1)), ("disu", (1,), (NLAY * NROW * NCOL,), (1,))],
)
def test_a_list_boundary_on_a_removed_cell_is_refused(tmp_path, kind, cellid, shape, removed):
    idomain = np.ones(shape, dtype=int)
    idomain[removed] = 0
    with pytest.raises(ValueError, match=rf"the built {kind.upper()} has REMOVED .*: \({cellid[0]}") as error:
        _build(tmp_path, kind, _chd(cellid), idomain=idomain)
    assert f"addresses cells the {kind.upper()} it was built with has removed" in str(error.value)


def test_an_observation_on_a_removed_dis_cell_is_refused(tmp_path):
    idomain = np.ones((NLAY, NROW, NCOL), dtype=int)
    idomain[0, 1, 1] = 0
    chd = mf.chd(stress_period_data={0: [[(1, 0, 0), 9.0]]},
                 observations={"c.csv": [("h", "head", (0, 1, 1))]})
    with pytest.raises(ValueError, match=r"name 1 cell\(s\) the built DIS has REMOVED .*'h' at \(0, 1, 1\)"):
        _build(tmp_path, "dis", chd, idomain=idomain)


def _lak(connection):
    return mf.lak.flopy(nlakes=1, packagedata=[[0, 9.0, 1]], perioddata={},
                        connectiondata=[[0, 0, *connection]])


def test_a_dis_lake_face_outside_its_cell_is_refused(tmp_path):
    """MF6 refuses it ('belev ... MUST BE >= cell bottom'): measured on this grid."""

    face = ((0, 1, 1), "HORIZONTAL", 0.1, 3.0, 8.0, 5.0, 10.0)  # the cell spans 5-10
    with pytest.raises(ValueError, match=r"has lake connections the DIS it was built with contradicts") as error:
        _build(tmp_path, "dis", _chd((1, 0, 0)), _lak(face))
    message = str(error.value)
    assert "reaching outside their cell in the built DIS" in message and "(0, 1, 1)" in message
    assert "correct each one against this DIS" in message


def test_a_dis_model_with_list_boundaries_sfr_and_lak_builds_and_runs(tmp_path):
    """The checks refuse nothing MF6 runs: every package on DIS cellids, in range."""

    sfr = mf.sfr.flopy(
        nreaches=2,
        packagedata=[[0, (0, 2, 1), 10.0, 1.0, 0.001, 8.0, 1.0, 1.0, 0.03, 1, 1.0, 0],
                     [1, (0, 2, 2), 10.0, 1.0, 0.001, 7.5, 1.0, 1.0, 0.03, 1, 1.0, 0]],
        connectiondata=[[0, -1], [1, 0]], perioddata={0: [[0, "INFLOW", 5.0]]},
    )
    lak = _lak(((0, 0, 2), "VERTICAL", 0.1, 5.0, 5.0, 0.0, 0.0))
    run = _build(
        tmp_path, "dis",
        mf.chd(stress_period_data={0: [[(1, 0, 0), 9.0], [(1, 2, 3), 8.0]]}),
        mf.drn(stress_period_data={0: [[(0, 1, 2), 8.0, 1.0]]}),
        mf.wel(stress_period_data={0: [[(1, 1, 1), -0.1]]}),
        sfr, lak,
    )
    success, report = run.execute()
    assert success, report[-10:]
    gwf = run.simulation.get_model("m")
    assert {p.package_type for p in gwf.packagelist} >= {"dis", "chd", "drn", "wel", "sfr", "lak"}


@pytest.mark.parametrize(
    ("kind", "cellid", "top", "bottom"),
    [("dis", (1, 2, 3), 5.0, 0.0), ("disv", (1, 11), 5.0, 0.0), ("disu", (23,), 5.0, 0.0),
     ("disu", (5,), 10.0, 5.0)],
)
def test_surface_water_validation_reads_the_grids_own_cellid(tmp_path, kind, cellid, top, bottom):
    """``validate_surface_water_configuration`` refused every cellid but a ``(layer,
    cell)`` pair as malformed, and read a cell's elevations by that pair. The stage
    warnings name the connected cell's own top and bottom, read through the grid."""

    from myflopy.modflow.mf6.surface_water_validation import validate_surface_water_configuration

    gwf = _build(tmp_path, kind, _chd((0, 0, 0) if kind == "dis" else (0, 0) if kind == "disv" else (0,)))
    gwf = gwf.simulation.get_model("m")
    codes = {}
    for stage in (50.0, -50.0):
        report = validate_surface_water_configuration(
            gwf, nlakes=1, lak_packagedata=[[0, 9.0, 1]],
            lak_connectiondata=[[0, 0, cellid, "VERTICAL", 0.1, 0.0, 0.0, 0.0, 0.0]],
            lak_perioddata={0: [[0, "STAGE", stage]]},
        )
        assert report.errors == [], report.errors
        codes.update({issue.code: issue.context for issue in report.warnings})
    assert codes["lak_stage_above_top"]["max_top"] == top
    assert codes["lak_stage_below_bottom"]["min_bottom"] == bottom


def test_surface_water_validation_refuses_another_grids_form(tmp_path):
    from myflopy.modflow.mf6.surface_water_validation import validate_surface_water_configuration

    gwf = _build(tmp_path, "dis", _chd((0, 0, 0))).simulation.get_model("m")
    report = validate_surface_water_configuration(
        gwf, nlakes=1, lak_packagedata=[[0, 9.0, 1]],
        lak_connectiondata=[[0, 0, (0, 5), "VERTICAL", 0.1, 0.0, 0.0, 0.0, 0.0]],
    )
    [issue] = report.errors
    assert issue.code == "lak_cellid_shape" and "(layer, row, col)" in issue.message
    report = validate_surface_water_configuration(
        gwf, nlakes=1, lak_packagedata=[[0, 9.0, 1]],
        lak_connectiondata=[[0, 0, (0, 0, 4), "VERTICAL", 0.1, 0.0, 0.0, 0.0, 0.0]],
    )
    assert [issue.code for issue in report.errors] == ["lak_cell_range"]


# -- the lazy areal/UZF builders on a DIS -------------------------------------------


def _removed_01():
    idomain = np.ones((NLAY, NROW, NCOL), dtype=int)
    idomain[0, 0, 1] = 0
    return idomain


def _records(run, kind):
    gwf = run.simulation.get_model("m")
    [package] = [p for p in gwf.packagelist if p.package_type == kind]
    if kind == "uzf":
        return [tuple(row[1]) for row in package.packagedata.get_data()]
    if kind in ("rcha", "evta"):
        return np.asarray(package.irch.get_data(0) if kind == "rcha" else package.ievt.get_data(0))
    return [tuple(row[0]) for row in package.stress_period_data.get_data(0)]


@pytest.mark.parametrize(
    ("label", "package", "kind"),
    [("rch", lambda: mf.rch(nper=1, recharge=1e-3), "rch"),
     ("evt", lambda: mf.evt(nper=1, rate=1e-3, depth=2.0, surface=10.0), "evt"),
     ("uzf", lambda: mf.uzf(nper=1, vks=1.0, thts=0.3, thtr=0.05, thti=0.1, eps=4.0, finf=1e-3), "uzf")],
)
def test_a_lazy_areal_or_uzf_package_resolves_on_a_dis(tmp_path, label, package, kind):
    """Declared without ``context=``, it resolved its cells from the built DIS's
    (nlay, nrow, ncol) idomain and refused, naming a reshape the caller cannot make.
    Each column's top active cell, written as ``(layer, row, col)``: column (0, 1) is
    removed in layer 0, so it is (1, 0, 1)."""

    run = _build(tmp_path, "dis", _chd((1, 2, 3)), package(), idomain=_removed_01())
    cells = _records(run, kind)
    assert cells[:3] == [(0, 0, 0), (1, 0, 0), (1, 0, 1)] if kind == "uzf" else cells[:3] == [
        (0, 0, 0), (1, 0, 1), (0, 0, 2)
    ]
    assert all(len(cell) == 3 for cell in cells)
    success, report = run.execute(silent=True)
    assert success, report[-10:]


def test_a_lazy_array_recharge_resolves_on_a_dis(tmp_path):
    run = _build(tmp_path, "dis", _chd((1, 2, 3)), mf.rch.array(recharge=1e-3), idomain=_removed_01())
    irch = _records(run, "rcha")
    assert irch.reshape(-1).tolist() == [0, 1] + [0] * 10
    assert run.execute(silent=True)[0]


def test_a_lazy_package_on_a_disu_is_refused_naming_the_records(tmp_path):
    with pytest.raises(ValueError, match=r"DISU, whose nodes have no columns .*\(node,\) cellids"):
        _build(tmp_path, "disu", _chd((12,)), mf.rch(nper=1, recharge=1e-3))


def test_a_lazy_gpkg_on_a_dis_is_refused_naming_the_records(tmp_path):
    import geopandas as gpd
    from shapely.geometry import Point

    path = tmp_path / "w.gpkg"
    gpd.GeoDataFrame({"name": ["w"]}, geometry=[Point(5.0, 5.0)], crs="EPSG:2927").to_file(path)
    with pytest.raises(ValueError, match=r"discretization is a DIS, whose cellids are \(layer, row, col\)"):
        _build(tmp_path, "dis", _chd((1, 2, 3)), mf.wel.gpkg(path, nper=1, rate=-1.0, layer_field=None))


def test_an_array_index_on_a_removed_dis_cell_says_what_mf6_drops(tmp_path):
    """``irch=0`` on a DIS whose (0, 0, 1) is removed: MF6 applies 1.1 of 1.2 ft3/d
    (measured), and the warning names the cell and that rate."""

    with pytest.warns(UserWarning, match=r"built DIS has REMOVED .*cell \(0, 0, 1\).* 0\.1 of 1\.2"):
        _build(tmp_path, "dis", _chd((1, 2, 3)), mf.rch.array(recharge=1e-3, irch=0),
               idomain=_removed_01())
