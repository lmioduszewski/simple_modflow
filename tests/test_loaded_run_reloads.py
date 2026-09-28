"""A reopened run loads its packages ONCE, and keeps what you change.

``LoadedMf6Run`` recorded the package NAMES FloPy reports (``CHD_0``) and checked the
TYPE asked for (``CHD``) against them. On any workspace whose package names are not
their types -- FloPy's defaults, a name file without pnames, a myflopy ``name=``
override -- every read reloaded the whole simulation, and so did probing a package
that is not there. A reload replaces the in-memory simulation, so an edit made
before one was silently dropped from the files ``run_simulation`` wrote: measured,
a WEL rate of -999 written back as -50.
"""

from __future__ import annotations

import flopy
import numpy as np
import pytest

import myflopy as mf
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.project.run_model import LoadedMf6Run


@pytest.fixture
def flopy_named_workspace(tmp_path):
    """A 3x3 DISV model whose boundary packages carry FloPy's default names."""

    n, d = 3, 10.0
    verts = np.array([[i * d, j * d] for j in range(n + 1) for i in range(n + 1)], dtype=float)
    iverts = [
        [j * (n + 1) + i, (j + 1) * (n + 1) + i, (j + 1) * (n + 1) + i + 1, j * (n + 1) + i + 1]
        for j in range(n)
        for i in range(n)
    ]
    xcyc = np.array([[(i + 0.5) * d, (j + 0.5) * d] for j in range(n) for i in range(n)])
    vor = VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)

    sim = flopy.mf6.MFSimulation(sim_name="named", sim_ws=str(tmp_path), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim)
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="named", save_flows=True)
    flopy.mf6.ModflowGwfdisv(gwf, nlay=1, top=10.0, botm=[0.0], **vor.get_disv_gridprops())
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0), 6.0], [(0, 8), 4.0]])
    flopy.mf6.ModflowGwfwel(gwf, stress_period_data=[[(0, 4), -50.0]])
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord="named.hds", budget_filerecord="named.cbc",
                           saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])
    sim.write_simulation(silent=True)
    success, _ = sim.run_simulation(silent=True)
    assert success
    assert "CHD_0" in gwf.get_package_list()  # the premise: name is not type
    return tmp_path


@pytest.fixture
def counted(monkeypatch):
    """Count full or partial FloPy loads."""

    calls = []
    original = LoadedMf6Run._load_simulation

    def counting(self, **kwargs):
        calls.append(kwargs.get("reason"))
        return original(self, **kwargs)

    monkeypatch.setattr(LoadedMf6Run, "_load_simulation", counting)
    return calls


def test_repeated_reads_load_once(flopy_named_workspace, counted):
    model = mf.load_mf6_run(flopy_named_workspace)
    for _ in range(3):
        assert len(model.chd.inputs.head.get()) == 2
    assert len(counted) <= 2, counted  # core, then one full load -- never one per read


def test_probing_an_absent_package_raises_without_reloading(flopy_named_workspace, counted):
    model = mf.load_mf6_run(flopy_named_workspace)
    model.package("chd")
    loads = len(counted)
    with pytest.raises(AttributeError):
        model.package("sfr")
    assert len(counted) == loads


def test_an_edit_survives_other_reads_and_reaches_the_files(flopy_named_workspace):
    model = mf.load_mf6_run(flopy_named_workspace)
    model.package("wel").stress_period_data.set_data({0: [((0, 4), -999.0)]})
    model.package("chd")
    with pytest.raises(AttributeError):
        model.package("sfr")
    model.run_simulation()

    written = flopy.mf6.MFSimulation.load(sim_ws=str(flopy_named_workspace), verbosity_level=0)
    rates = written.get_model("named").get_package("wel").stress_period_data.get_data()[0]["q"]
    assert rates.tolist() == [-999.0]


def test_an_edit_to_a_core_package_reaches_the_files_too(flopy_named_workspace):
    """OC was served from a PARTIAL load, which the full load writing needs dropped."""

    model = mf.load_mf6_run(flopy_named_workspace)
    model.package("oc").saverecord.set_data({0: [("HEAD", "LAST"), ("BUDGET", "LAST")]})
    model.run_simulation()
    written = " ".join((flopy_named_workspace / "named.oc").read_text().upper().split())
    assert "SAVE HEAD LAST" in written


def test_a_reopened_run_is_quiet_at_verbosity_zero(flopy_named_workspace, capsys):
    """The loader floored its verbosity at 1, so the default 0 printed anyway."""

    model = mf.load_mf6_run(flopy_named_workspace)
    capsys.readouterr()
    model.chd.inputs.head.get()
    assert capsys.readouterr().out == ""


# --- the model's own options survive a reload (ledger 286) --------------------


@pytest.mark.parametrize(
    "declared, expected_block",
    [
        ({}, []),
        ({"newtonoptions": "NEWTON"}, ["NEWTON"]),
        ({"newtonoptions": ["under_relaxation"]}, ["NEWTON  UNDER_RELAXATION"]),
        ({"print_flows": True, "save_flows": True}, ["PRINT_FLOWS", "SAVE_FLOWS"]),
    ],
)
def test_model_options_round_trip_through_a_reload(tmp_path, declared, expected_block):
    """Rebuilding specs from a finished workspace reads each PACKAGE off the loaded
    model -- and the model's own options block belongs to no package, so it was dropped
    in silence.

    Measured on a real 18,107-cell model that lost `NEWTON`: 90,408 cells converted to
    DRY straight from the initial conditions, period-1 recharge fell from 2,854,736 to
    1,645 ft3/d, every head-dependent boundary went to zero, and the solver diverged at
    a 200% budget discrepancy. Restoring the one option reproduced the original run's
    iteration trace digit for digit (ledger 286).

    Both directions are asserted -- a model that declares nothing must come back
    declaring nothing. `newtonoptions=""` and `[]` write NOTHING at all (measured on all
    four spellings), so a reader that normalized to the empty form would drop exactly
    the option it exists to carry, while every "did it round-trip" test still passed.
    """

    grid = _three_by_three_grid()
    packages = [
        mf.disv(nlay=1, top=10.0, botm=[0.0], **grid.get_disv_gridprops()),
        mf.ic(strt=5.0), mf.npf(k=1.0), mf.oc(head_filerecord="m.hds"),
    ]
    context = mf.ModelContext(grid=grid, domain=np.ones((1, int(grid.ncpl)), dtype=int))

    def write(workspace, **options):
        model = mf.gwf("m", context=context, packages=packages, **options)
        sim = mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
        sim.build_flopy(workspace).simulation.write_simulation(silent=True)
        body = (workspace / "m.nam").read_text().split("BEGIN options")[1].split("END options")[0]
        return [line.strip() for line in body.splitlines() if line.strip()]

    assert write(tmp_path / "a", **declared) == expected_block

    carried = mf.load_mf6_run(str(tmp_path / "a")).model_options
    assert write(tmp_path / "b", **carried) == expected_block


def test_model_options_keeps_a_filename_option_rather_than_dropping_it(tmp_path):
    """`_options_rows` skips every `*filerecord` because an output path is noise in a
    DIFF. It is CONTENT in a rebuild, so this reader keeps them -- dropping them here
    would be the same silent loss in miniature, one tier down."""

    grid = _three_by_three_grid()
    context = mf.ModelContext(grid=grid, domain=np.ones((1, int(grid.ncpl)), dtype=int))
    model = mf.gwf(
        "m", context=context, list="custom.lst",
        packages=[mf.disv(nlay=1, top=10.0, botm=[0.0], **grid.get_disv_gridprops()),
                  mf.ic(strt=5.0), mf.npf(k=1.0), mf.oc()],
    )
    sim = mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    sim.build_flopy(tmp_path).simulation.write_simulation(silent=True)

    assert mf.load_mf6_run(str(tmp_path)).model_options.get("list") == "custom.lst"


def _three_by_three_grid():
    """A 3x3 DISV grid, built the way this module's other fixture builds one."""

    n, d = 3, 10.0
    verts = np.array([[i * d, j * d] for j in range(n + 1) for i in range(n + 1)], dtype=float)
    iverts = [
        [j * (n + 1) + i, (j + 1) * (n + 1) + i, (j + 1) * (n + 1) + i + 1, j * (n + 1) + i + 1]
        for j in range(n)
        for i in range(n)
    ]
    xcyc = np.array([[(i + 0.5) * d, (j + 0.5) * d] for j in range(n) for i in range(n)])
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)
