"""A model's package NAMES and TYPES, shown both -- live and reopened alike.

A live model listed FloPy's names (``CHD_0``); a reopened run listed file-suffix
TYPES (``CHD``) and missed every GWT/GWE/PRT package, array recharge and any
custom filename. So the same model said different things live and reopened, and
code that meant "is there a CHD?" asked a list of names (ledger 228). Now
``package_names`` are the packages' own names everywhere, ``package_types`` what
they are, once each -- array recharge counting as RCH, since MF6 has one recharge
package -- and ``packages.summary()`` shows both.
"""

from __future__ import annotations

import shutil

import flopy
import numpy as np
import pytest

import myflopy as mf
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.utils.validators import package_name_validator
from myflopy.project.run_model import LoadedMf6Run

NAMES = ["CHD_EAST", "CHD_WEST", "DISV", "IC", "NPF", "OC", "RCH_ARRAY", "RCH_LIST"]
TYPES = ["CHD", "DISV", "IC", "NPF", "OC", "RCH"]


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    n, d = 3, 10.0
    verts = np.array([[i * d, j * d] for j in range(n + 1) for i in range(n + 1)], dtype=float)
    iverts = [
        [j * (n + 1) + i, (j + 1) * (n + 1) + i, (j + 1) * (n + 1) + i + 1, j * (n + 1) + i + 1]
        for j in range(n)
        for i in range(n)
    ]
    xcyc = np.array([[(i + 0.5) * d, (j + 0.5) * d] for j in range(n) for i in range(n)])
    vor = VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)
    context = mf.ModelContext(grid=vor, domain=np.ones((1, n * n), dtype=int))
    model = mf.gwf(
        "flow",
        context=context,
        packages=[
            mf.disv(nlay=1, top=10.0, botm=[0.0], **vor.get_disv_gridprops()),
            mf.ic(strt=5.0),
            mf.npf(k=1.0),
            mf.chd(stress_period_data={0: [[(0, 0), 6.0]]}, name="chd_west"),
            mf.chd(stress_period_data={0: [[(0, 8), 4.0]]}, name="chd_east"),
            mf.rch(stress_period_data={0: [[(0, 4), 1e-4]]}, name="rch_list"),
            mf.rch.array(recharge={0: np.full(9, 5e-5)}, name="rch_array", context=context),
            mf.oc(),
        ],
    )
    run = mf.Project(tmp_path_factory.mktemp("names"), name="names").prepare_run(
        "r", mf.simulation(model, name="names")
    )
    success, report = run.execute(silent=True)
    assert success, report
    return run


def test_a_live_model_lists_names_and_types(run):
    live = run.model("flow")
    assert live.package_names == NAMES
    assert live.package_types == TYPES


def test_a_reopened_run_says_the_same_without_loading(run, monkeypatch):
    loads = []
    monkeypatch.setattr(LoadedMf6Run, "_load_simulation", lambda self, **kw: loads.append(kw))
    reopened = mf.load_mf6_run(run.workspace)
    assert (reopened.package_names, reopened.package_types) == (NAMES, TYPES)
    assert loads == []


def test_a_name_file_without_names_gets_the_names_flopy_gives(run, tmp_path):
    copy = tmp_path / "nameless"
    shutil.copytree(run.workspace, copy)
    namefile = copy / "flow.nam"
    lines, inside = [], False
    for line in namefile.read_text().splitlines():
        upper = line.strip().upper()
        if upper.startswith("BEGIN PACKAGES"):
            inside = True
        elif upper.startswith("END PACKAGES"):
            inside = False
        elif inside and line.strip():
            line = "  " + "  ".join(line.split()[:2])
        lines.append(line)
    namefile.write_text("\n".join(lines) + "\n")

    flopys = flopy.mf6.MFSimulation.load(sim_ws=str(copy), verbosity_level=0).get_model("flow")
    reopened = mf.load_mf6_run(copy)
    assert reopened.package_names == sorted(n.upper() for n in flopys.get_package_list())
    assert reopened.package_types == TYPES


def test_a_transport_model_lists_its_own_packages(tmp_path):
    """GWT/GWE/PRT packages were invisible: the suffixes knew only GWF's ten."""

    (tmp_path / "mfsim.nam").write_text(
        "BEGIN models\n  gwf6  flow.nam  flow\n  gwt6  conc.nam  conc\nEND models\n"
    )
    (tmp_path / "flow.nam").write_text("BEGIN packages\n  DISV6  flow.disv  disv\nEND packages\n")
    (tmp_path / "conc.nam").write_text(
        "BEGIN packages\n  DISV6 conc.disv disv\n  IC6 conc.ic ic\n  ADV6 conc.adv adv\n"
        "  DSP6 conc.dsp dsp\n  MST6 conc.mst mst\n  SSM6 conc.ssm ssm\n  CNC6 conc.cnc cnc\n"
        "  OC6 conc.oc oc\nEND packages\n"
    )
    reopened = mf.load_mf6_run(tmp_path, model_name="conc")
    assert reopened.package_types == ["ADV", "CNC", "DISV", "DSP", "IC", "MST", "OC", "SSM"]


def test_the_summary_shows_every_package_with_its_type(run):
    frame = run.model("flow").packages.summary()
    assert list(frame["package"]) == NAMES
    assert dict(zip(frame["package"], frame["type"], strict=True))["RCH_ARRAY"] == "RCH"
    assert list(frame["type"]).count("CHD") == 2


def test_a_package_type_is_valid_whatever_the_package_is_named(run):
    """The validator compared 'chd' against names like 'chd_west' and refused it."""

    assert package_name_validator(run.model("flow"), "chd") == "chd"
