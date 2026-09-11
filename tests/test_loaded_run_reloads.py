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
        assert len(model.packages.chd.inputs.head.get()) == 2
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
