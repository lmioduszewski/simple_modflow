"""The dependent-variable reader behind ``model.hds``/``.conc``/``.temp`` (ledger 598).

Four defects of ``DependentVariableFile``, found in the review of ``model.hds.remap``
(ledger 595) and each measured before it was fixed:

- ``get_kstpkper()`` named the record BEFORE the one meant, and ``nper``/``numstp`` were
  one short: ``__init__`` stored ``kstpkper`` zero-based over FloPy's one-based list, and
  FloPy's ``get_kstpkper`` subtracts one from it again.
- The reader opened ``<model>.<suffix>`` whatever the OC declared, so a run declaring
  ``HEAD FILEOUT other.hds`` raised ``FileNotFoundError``; and ``output_files`` had no
  ``temperature`` kind at all.
- ``model.hds`` built a new reader on every access and closed none (one ``map()`` built
  five); ``model.kstpkper`` was cached once and never invalidated.
- A reader held across a re-run read the NEW file's bytes through its OLD index.

Every test here fails on the code before the fix, bar the ones naming the mutant they kill.
"""

from __future__ import annotations

import gc
import warnings
from pathlib import Path
from types import SimpleNamespace

import flopy
import numpy as np
import pytest

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.headsplus import HeadsPlus
from myflopy.modflow.mf6.output_files import _declared_field_output, _declared_output_path


def _grid():
    return irregular_voronoi_grid(nrow=4, ncol=4, cell_size=100.0)


def _simulation(vor, *, nper, nstp=1, head_file="m.hds", west=50.0):
    """A one-layer transient model whose recharge grows each period, so every record differs."""

    props = vor.get_disv_gridprops()
    edges = vor.gdf_vorPolys.geometry.bounds
    west_cells = np.flatnonzero(edges.minx <= edges.minx.min() + 1e-6)
    east_cells = np.flatnonzero(edges.maxx >= edges.maxx.max() - 1e-6)
    model = mf.gwf("m", grid=vor, packages=[
        mf.disv(nlay=1, ncpl=props["ncpl"], nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=100.0, botm=0.0),
        mf.npf(k=10.0), mf.ic(strt=50.0), mf.sto(transient={0: True}, sy=0.2, ss=1e-5),
        mf.chd(stress_period_data={0: [[(0, int(c)), west] for c in west_cells]
                                       + [[(0, int(c)), 40.0] for c in east_cells]}),
        mf.rch(stress_period_data={
            per: [[(0, c), 0.001 * (per + 1)] for c in range(props["ncpl"])]
            for per in range(nper)
        }),
        mf.oc(head_filerecord=head_file, saverecord=[("HEAD", "ALL")]),
    ])
    return mf.simulation(model, name="s", packages=[
        mf.tdis(nper=nper, perioddata=[(10.0, nstp, 1.0)] * nper), mf.ims(models=["m"])])


def _read(path):
    """The last record of a head file, read by FloPy alone."""

    with flopy.utils.HeadFile(str(path)) as raw:
        return raw.get_data()


def _execute(run):
    success, report = run.execute()
    assert success, "\n".join(report[-20:])
    return run


def _run(root, simulation, name="r", project=None):
    """Build and execute ``simulation`` as run ``name``, reusing its directory if it exists."""

    with warnings.catch_warnings():
        # A built grid cannot be frozen into the run's snapshot; that is not under test.
        warnings.filterwarnings("ignore", "Run 'r' built but its frozen snapshot", UserWarning)
        project = project or mf.Project(root / "p", name="p")
        return project, _execute(project.prepare_run(name, simulation, overwrite=True))


@pytest.fixture(scope="module")
def stepped(tmp_path_factory):
    """Three periods of two steps, the heads saved every step to ``other.hds``."""

    vor = _grid()
    _project, run = _run(tmp_path_factory.mktemp("stepped"), _simulation(
        vor, nper=3, nstp=2, head_file="other.hds"))
    return SimpleNamespace(run=run, vor=vor, workspace=Path(run.workspace))


# -- 1. the clock ---------------------------------------------------------------------------


@pytest.mark.slow
def test_get_kstpkper_nper_and_numstp_name_the_records_the_file_saved(stepped):
    hds = stepped.run.model().hds
    raw = flopy.utils.HeadFile(str(hds.output_path))
    try:
        truth = [(int(kstp), int(kper)) for kstp, kper in raw.get_kstpkper()]
        assert truth == [(0, 0), (1, 0), (0, 1), (1, 1), (0, 2), (1, 2)]
        assert not np.array_equal(raw.get_data(kstpkper=truth[-1]), raw.get_data(kstpkper=truth[-2]))
        assert hds.kstpkper == truth
        assert hds.get_kstpkper() == truth
        assert hds.get_kstpkper() == truth, "a second call drifted"
        assert (hds.nper, hds.numstp) == (3, 2)
        np.testing.assert_array_equal(
            hds.get_data(kstpkper=hds.get_kstpkper()[-1]), raw.get_data(kstpkper=truth[-1])
        )
    finally:
        raw.close()


@pytest.mark.slow
def test_a_one_period_file_has_one_period(tmp_path):
    _project, run = _run(tmp_path, _simulation(_grid(), nper=1))
    hds = run.model().hds
    assert hds.get_kstpkper() == [(0, 0)]
    assert (hds.nper, hds.numstp) == (1, 1)


# -- 2. the file the OC declares --------------------------------------------------------------


@pytest.mark.slow
def test_the_reader_opens_the_head_file_the_output_control_declares(stepped):
    assert not (stepped.workspace / "m.hds").exists()
    for model in (stepped.run.model(), mf.load_mf6_run(stepped.workspace, verbosity_level=0)):
        assert Path(model.hds.output_path) == stepped.workspace / "other.hds"
        assert model.kstpkper[-1] == (1, 2)


@pytest.mark.slow
def test_a_readable_oc_declaring_no_head_file_is_refused_not_read_from_a_stale_hds(tmp_path):
    """MF6 refuses SAVE HEAD without HEAD FILEOUT, so an ``m.hds`` there is another run's."""

    _project, run = _run(tmp_path, _simulation(_grid(), nper=2))
    workspace = Path(run.workspace)
    assert (workspace / "m.hds").is_file()
    (workspace / "m.oc").write_text("BEGIN options\nEND options\n")
    with pytest.raises(ValueError, match=(
            r"model.hds: model 'm''s output control 'm.oc' declares no HEAD FILEOUT, so the "
            r"run saved no head")):
        mf.load_mf6_run(workspace, verbosity_level=0).hds


def test_temperature_is_a_kind_the_output_files_know(tmp_path):
    """``output_files`` had no ``temperature`` kind: both of its fallbacks raised KeyError."""

    assert _declared_output_path(tmp_path, "e", "temperature") == tmp_path / "e.ucn"
    assert _declared_field_output(tmp_path, "e", "temperature", caller="x") == tmp_path / "e.ucn"
    (tmp_path / "e.nam").write_text("BEGIN packages\n  OC6  e.oc  oc\nEND packages\n")
    (tmp_path / "e.oc").write_text("BEGIN options\n  BUDGET FILEOUT e.cbc\nEND options\n")
    with pytest.raises(ValueError, match=(
            r"model.temp: model 'e''s output control 'e.oc' declares no TEMPERATURE FILEOUT"
            r".*mf.oc\(temperature_filerecord=\"e.ucn\"")):
        _declared_field_output(tmp_path, "e", "temperature", caller="model.temp")
    (tmp_path / "e.oc").write_text("BEGIN options\n  TEMPERATURE FILEOUT e.tem\nEND options\n")
    assert _declared_field_output(tmp_path, "e", "temperature", caller="x") == tmp_path / "e.tem"


@pytest.mark.slow
def test_model_temp_reads_the_temperature_file_its_oc_declares(tmp_path):
    from myflopy.modflow.mf6.simulation.accessors import get_temp
    from tests.test_gwt_gwe_results import _coupled_run

    run = _coupled_run(tmp_path, "gwe", field_file="trans.tem")
    model = run.model("trans")
    assert Path(model.temp.output_path).name == "trans.tem"
    assert np.isfinite(model.temp.array()).any()
    # The kept reader answers for the unit it was built with, and no other (mutant: drop
    # the unit from the cache key).
    assert model.temp.unit == "°C"
    assert get_temp(model, unit="K").unit == "K"
    assert model.temp.unit == "°C"


@pytest.mark.slow
def test_a_duck_typed_model_s_profile_reads_the_declared_head_file_too(stepped):
    """``InterpolatedSurface.hds`` builds its own reader for a model with no
    ``_field_reader``; it opened ``<name>.hds``."""

    from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface

    duck = SimpleNamespace(model_output_folder_path=stepped.workspace, name="m")
    surface = SimpleNamespace(_hds=None, model=duck, vor=stepped.vor)
    reader = InterpolatedSurface.hds.fget(surface)
    assert Path(reader.output_path) == stepped.workspace / "other.hds"


# -- 3. one reader per model while its file is unchanged, and every handle closed ----------


@pytest.mark.slow
def test_model_hds_is_one_reader_while_its_file_is_unchanged(stepped):
    from myflopy.modflow.mf6 import headsplus

    model = stepped.run.model()
    model.hds.close()       # a reader closed by hand is not handed back (mutant: drop the check)
    hds = model.hds
    assert not hds.file.closed
    assert model.hds is hds
    built = []
    original = headsplus.DependentVariableFile.__init__

    def counting(self, *args, **kwargs):
        built.append(self)
        original(self, *args, **kwargs)

    headsplus.DependentVariableFile.__init__ = counting
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model.hds.map()
            model.plot.map(type="hds")
    finally:
        headsplus.DependentVariableFile.__init__ = original
    assert built == [], f"a map built {len(built)} more reader(s)"


@pytest.mark.slow
def test_a_reader_closes_its_file_when_it_is_collected(stepped):
    """HeadsPlus refers to itself (``.hds``), so only the cycle collector frees it."""

    reader = HeadsPlus(hds_path=stepped.workspace / "other.hds", vor=stepped.vor)
    handle = reader.file
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ResourceWarning)
        del reader
        gc.collect()
    assert handle.closed
    assert not [w for w in caught if issubclass(w.category, ResourceWarning)]


@pytest.mark.slow
def test_a_file_that_cannot_be_indexed_leaves_no_handle_open(stepped, tmp_path):
    """FloPy opens the file, then indexes it: a file cut inside a record header -- one
    still being written -- raises in between. (An EMPTY file leaks inside FloPy's own
    precision probe, ``get_headfile_precision``, before the reader opens anything.)"""

    data = (stepped.workspace / "other.hds").read_bytes()
    (tmp_path / "m.hds").write_bytes(data[: len(data) // 6 + 20])   # record 1 + 20 bytes
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ResourceWarning)
        with pytest.raises(EOFError):
            HeadsPlus(hds_path=tmp_path / "m.hds", vor=stepped.vor)
        gc.collect()
    assert not [w for w in caught if issubclass(w.category, ResourceWarning)]


@pytest.mark.slow
def test_the_per_cell_lists_are_built_on_first_use(stepped):
    """They were a reader's cost (0.76-0.82 s on a 20,323-cell grid, built cell by cell;
    ``docs/reproducers/field_reader_per_model.py``) and nothing reads them."""

    reader = HeadsPlus(hds_path=stepped.workspace / "other.hds", vor=stepped.vor)
    names = ("vor_list", "cell_list", "area_list", "x_list", "y_list")
    assert not set(names) & set(vars(reader))
    polygons = stepped.vor.gdf_vorPolys.geometry
    assert reader.cell_list == list(range(len(polygons)))
    np.testing.assert_allclose(reader.area_list, polygons.area)
    np.testing.assert_allclose(reader.x_list, polygons.centroid.x)
    np.testing.assert_allclose(reader.y_list, polygons.centroid.y)


@pytest.mark.slow
def test_a_rerun_into_the_same_directory_gets_a_new_reader_and_clock(tmp_path):
    """MF6 rewrites the head file in place -- the same inode, a new size and mtime.

    A reader HELD across it refuses to read (its index no longer describes the bytes),
    while one whose file was not touched keeps reading: the refusal is the rewrite, not
    the age.
    """

    vor = _grid()
    project, run = _run(tmp_path, _simulation(vor, nper=2))
    model = run.model()
    held = model.hds
    assert model.kstpkper == [(0, 0), (0, 1)]
    inode = (Path(run.workspace) / "m.hds").stat().st_ino

    _run(tmp_path, _simulation(vor, nper=3, west=60.0), project=project)
    assert (Path(run.workspace) / "m.hds").stat().st_ino == inode, "MF6 no longer rewrites in place"
    assert model.hds is not held
    assert model.kstpkper == [(0, 0), (0, 1), (0, 2)]
    np.testing.assert_array_equal(model.hds.get_data(), _read(Path(run.workspace) / "m.hds"))
    for read in (lambda: held.get_data(), lambda: held.get_ts((0, 0)), lambda: held.get()):
        with pytest.raises(RuntimeError, match="m.hds changed after this reader indexed it"):
            read()

    # The OC now names another file: the kept reader is of an UNCHANGED m.hds, and is
    # replaced all the same -- the stamp compared is the NEW file's.
    kept = model.hds
    _run(tmp_path, _simulation(vor, nper=3, head_file="m2.hds"), project=project)
    assert Path(model.hds.output_path).name == "m2.hds"
    np.testing.assert_array_equal(kept.get_data(), _read(Path(run.workspace) / "m.hds"))


@pytest.mark.slow
def test_a_new_grid_gets_a_new_reader(stepped):
    """Mutant: drop the grid from the cache key -- the reader's maps would draw the old one."""

    model = stepped.run.model()
    first = model.hds
    model.vor = _grid()
    try:
        assert model.hds is not first
        assert model.hds.vor is model.vor
    finally:
        model.vor = stepped.vor
