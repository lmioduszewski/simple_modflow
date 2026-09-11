"""Particle tracking on a grid whose origin is not zero.

Every USG import is one: ``to_mf6`` writes the mesh near zero and declares its
true corner, because MF6 loses precision on raw State Plane vertices. MF6 then
reads and writes particle coordinates in the model's OWN frame, while FloPy's cell
centres and every myflopy map are in the WORLD frame -- 282 miles apart for a
(1,340,000, 650,000) ft origin. Release points were written in the world frame
(ledger 218), so MF6 released particles off the grid.
"""

from __future__ import annotations

import pytest

import myflopy as mf
from myflopy.modflow.mf6.prt import PRTReleasePoints
from myflopy.modflow.usg import read_usg
from tests.test_usg_import import (  # noqa: F401 -- usg_workspace is a fixture
    CENTERS,
    CORNERS,
    LEVELS,
    NCPL,
    NODES,
    PLAN,
    usg_workspace,
)

X0, Y0 = 1_340_000.0, 650_000.0


def _write_offset_gsf(path):
    """The synthetic model's geometry, moved to State Plane."""

    lines = ["UNSTRUCTURED", f"{NODES} 1 1 1", str(len(PLAN) * len(LEVELS))]
    for z in LEVELS:
        lines += [f"{x + X0:.6f} {y + Y0:.6f} {z:.6f}" for x, y in PLAN]
    for node in range(NODES):
        layer, cell = divmod(node, NCPL)
        xc, yc = CENTERS[cell]
        zc = (LEVELS[layer] + LEVELS[layer + 1]) / 2.0
        top = [layer * len(PLAN) + c + 1 for c in CORNERS[cell]]
        bot = [(layer + 1) * len(PLAN) + c + 1 for c in CORNERS[cell]]
        ids = " ".join(str(i) for i in top + bot)
        lines.append(f"{node + 1} {xc + X0:.6f} {yc + Y0:.6f} {zc:.6f} {layer + 1} 8 {ids}")
    path.write_text("\n".join(lines) + "\n")


@pytest.fixture(scope="module")
def offset_run(usg_workspace, tmp_path_factory):
    """The synthetic USG model at State Plane, converted and prepared (not yet run)."""

    gsf = usg_workspace.parent / "tiny_offset.gsf"
    _write_offset_gsf(gsf)
    sim = read_usg(usg_workspace, gsf=gsf, crs="EPSG:2927").to_mf6(
        "tiny", start_date_time="2020-01-01"
    )
    project = mf.Project(tmp_path_factory.mktemp("offset"), name="offset")
    return project.prepare_run("r1", sim)


def test_the_converted_grid_really_is_offset(offset_run):
    grid = offset_run.model("tiny").gwf.modelgrid
    assert (grid.xoffset, grid.yoffset) == (X0, Y0)


def test_release_points_from_cells_are_in_the_models_own_frame(offset_run):
    points = PRTReleasePoints.from_cells(offset_run.model("tiny"), [0, 3])
    assert [(row[2], row[3]) for row in points.packagedata] == [(5.0, 5.0), (15.0, 15.0)]


def test_release_points_from_world_points_are_in_the_models_own_frame(offset_run):
    points = PRTReleasePoints.from_points(offset_run.model("tiny"), [(X0 + 15.0, Y0 + 5.0)])
    (row,) = points.packagedata
    assert row[1] == (0, 1)
    assert (row[2], row[3]) == pytest.approx((15.0, 5.0))
