"""`VoronoiGridPlus.from_gsf`: the plan view of a MODFLOW-USG grid.

A `.gsf` is the only file in a USG model that knows where the cells ARE -- DISU
stores connectivity and areas but no coordinates. It describes 3-D cells (eight
vertices per hexahedron, the bottom four under the top four) and repeats the
whole grid for every layer, so the reader has two jobs: collapse each cell to
its distinct plan corners, and take one layer.

Getting either wrong is quiet rather than loud. A double-traced ring still
fills on a map -- it just reports double the area, fails `is_valid`, and loses
hover, because maplibre's hit-testing gives up on a self-overlapping polygon.
That is the failure this file is here to keep from coming back, so the geometry
assertions are on validity and exact area, not on the vertex list alone.
"""

from __future__ import annotations

import numpy as np
import pytest

import myflopy as mf
from myflopy.modflow.mf6.grid.helpers import get_griddata_from_gsf, signed_area

#: One cell of the synthetic grid, in projected feet.
CELL = 100.0
#: Cells per side; the grid is NCELL x NCELL per layer.
NCELL = 3
#: Lower-left corner, in EPSG:2927 territory so the lat/lon reprojection is sane.
X0, Y0 = 1200000.0, 600000.0


def write_gsf(path, *, nlay=2, header="UNSTRUCTURED GWF", solid=True, comment=False):
    """Write a synthetic .gsf.

    ``solid=True`` writes hexahedra (8 vertices per cell, one vertex table per
    elevation level) the way a real USG .gsf does; ``solid=False`` writes flat
    4-vertex cells, which some writers emit and which must survive the same path
    untouched.
    """

    levels = [100.0 - 50.0 * i for i in range(nlay + 1)] if solid else [0.0]
    vid, verts = {}, []
    for level, z in enumerate(levels):
        for row in range(NCELL + 1):
            for col in range(NCELL + 1):
                vid[(level, row, col)] = len(verts) + 1  # .gsf vertex ids are 1-based
                verts.append((X0 + col * CELL, Y0 + row * CELL, z))

    nodes = []
    for lay in range(1, nlay + 1):
        for row in range(NCELL):
            for col in range(NCELL):
                corners = [(row, col), (row, col + 1), (row + 1, col + 1), (row + 1, col)]
                if solid:
                    cell = [vid[(lay - 1, r, c)] for r, c in corners]
                    cell += [vid[(lay, r, c)] for r, c in corners]
                else:
                    cell = [vid[(0, r, c)] for r, c in corners]
                nodes.append(
                    f"{len(nodes) + 1} {X0 + (col + 0.5) * CELL} {Y0 + (row + 0.5) * CELL} "
                    f"{-50.0 * lay} {lay} {len(cell)} " + " ".join(str(v) for v in cell)
                )

    lines = ["# written by the test suite"] if comment else []
    lines += [header, str(len(nodes)), str(len(verts))]
    lines += [f"{x} {y} {z}" for x, y, z in verts]
    lines += nodes
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture
def gsf(tmp_path):
    return write_gsf(tmp_path / "grid.gsf")


def test_reads_one_layer_not_the_whole_stack(gsf):
    verts, iverts, xcyc = get_griddata_from_gsf(gsf)

    # 2 layers x 9 cells in the file; a plan view is 9.
    assert len(iverts) == NCELL * NCELL
    assert xcyc.shape == (NCELL * NCELL, 2)
    # 4x4 plan positions, NOT the 3 elevation levels' worth in the file.
    assert verts.shape == ((NCELL + 1) ** 2, 2)


def test_hexahedra_collapse_to_plan_corners(gsf):
    _, iverts, _ = get_griddata_from_gsf(gsf)

    assert [len(ring) for ring in iverts] == [4] * (NCELL * NCELL)
    assert all(len(set(ring)) == len(ring) for ring in iverts)


def test_rings_are_clockwise_like_the_disu_reader(gsf):
    verts, iverts, _ = get_griddata_from_gsf(gsf)

    assert all(signed_area(verts[ring]) < 0 for ring in iverts)


def test_polygons_are_valid_with_true_areas(gsf):
    vor = mf.VoronoiGridPlus.from_gsf(gsf)
    polys = vor.gdf_vorPolys.geometry

    assert vor.ncpl == NCELL * NCELL
    assert bool(polys.is_valid.all())
    # The failure this guards: a double-traced ring reports 2 * CELL**2.
    assert np.allclose(polys.area.to_numpy(), CELL**2)


def test_map_hover_carries_the_cell_identity(gsf):
    """The symptom that started this: a map that draws but does not hover.

    Plotly hover is customdata + a template, so the payload is checkable without
    a browser -- and a degenerate grid never gets this far, because maplibre
    drops the hit-test rather than the fill.
    """

    trace = mf.VoronoiGridPlus.from_gsf(gsf).plot.map().fig.data[0]

    assert "Cell No." in trace.hovertemplate
    assert len(trace.customdata) == NCELL * NCELL
    assert trace.customdata[0][0] == 0
    assert trace.customdata[0][1] == pytest.approx(CELL**2)


def test_layers_share_plan_geometry_but_not_centers(gsf):
    first, iverts_0, xcyc_0 = get_griddata_from_gsf(gsf, layer=0)
    second, iverts_1, xcyc_1 = get_griddata_from_gsf(gsf, layer=1)

    assert np.array_equal(first, second)
    assert iverts_0 == iverts_1
    assert np.array_equal(xcyc_0, xcyc_1)  # same plan centers, different node ids


def test_flat_gsf_needs_no_collapsing(tmp_path):
    flat = write_gsf(tmp_path / "flat.gsf", nlay=1, solid=False)
    _, iverts, _ = get_griddata_from_gsf(flat)

    assert [len(ring) for ring in iverts] == [4] * (NCELL * NCELL)


@pytest.mark.parametrize("header", ["UNSTRUCTURED", "UNSTRUCTURED GWF", "unstructured gwf"])
def test_header_forms_all_load(tmp_path, header):
    """flopy's own reader rejects the two-word header; ours must not."""

    path = write_gsf(tmp_path / f"{header.replace(' ', '_')}.gsf", header=header)

    assert len(get_griddata_from_gsf(path)[1]) == NCELL * NCELL


def test_comments_are_skipped(tmp_path):
    path = write_gsf(tmp_path / "commented.gsf", comment=True)

    assert len(get_griddata_from_gsf(path)[1]) == NCELL * NCELL


def test_a_file_that_is_not_a_gsf_says_so(tmp_path):
    path = tmp_path / "notagrid.gsf"
    path.write_text("BEGIN OPTIONS\nEND OPTIONS\n")

    with pytest.raises(ValueError, match="not a grid specification file"):
        get_griddata_from_gsf(path)


def test_truncated_file_is_named_not_silently_short(tmp_path):
    path = write_gsf(tmp_path / "short.gsf")
    lines = path.read_text().splitlines()
    path.write_text("\n".join(lines[:-3]) + "\n")

    with pytest.raises(ValueError, match="provides"):
        get_griddata_from_gsf(path)


def test_a_cell_that_miscounts_its_vertices_is_named(tmp_path):
    path = write_gsf(tmp_path / "miscount.gsf")
    lines = path.read_text().splitlines()
    parts = lines[-1].split()
    parts[5] = str(int(parts[5]) + 1)          # declare one more than provided
    lines[-1] = " ".join(parts)
    path.write_text("\n".join(lines) + "\n")

    with pytest.raises(ValueError, match="declares"):
        get_griddata_from_gsf(path)


def test_layer_out_of_range_names_the_range(gsf):
    with pytest.raises(ValueError, match="layer must be between 0 and 1"):
        get_griddata_from_gsf(gsf, layer=2)
