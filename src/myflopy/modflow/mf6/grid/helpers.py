from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import shapely as shp
from flopy.mf6 import MFSimulation, ModflowGwf, ModflowGwfdisu

from myflopy._logging import get_logger

logger = get_logger(__name__)


def flatten(values):
    """Flatten one level of nesting: a sequence of sequences into a single list."""

    return [item for sublist in values for item in sublist]


def signed_area(coords):
    """Shoelace formula; >0 => CCW, <0 => CW."""
    x = coords[:, 0]
    y = coords[:, 1]
    return 0.5 * np.sum(x * np.roll(y, -1) - y * np.roll(x, -1))


def rows_truncate_at_first_missing(rec, treat_nan=True):
    """
    Truncate rows of a record array at the first occurrence of a missing value.
    """
    names = rec.dtype.names

    def is_missing(value):
        """True if ``value`` is ``None`` (or NaN when ``treat_nan``)."""

        return (value is None) or (treat_nan and isinstance(value, float) and np.isnan(value))

    out = []
    for row_record in rec:
        row = []
        for name in names:
            value = row_record[name]
            value = value.item() if hasattr(value, "item") else value
            if is_missing(value):
                break
            row.append(value)
        out.append(row)
    return out


def get_griddata_from_disu(disu_path: Path):
    """
    Read a DISU file and extract vertices, cell vertex indices, and centroids.
    """
    sim = MFSimulation(
        sim_name="dummy_sim",
        sim_ws=str(disu_path.parent),
    )
    gwf = ModflowGwf(
        sim,
        modelname="dummy_model",
    )
    disu = ModflowGwfdisu(
        gwf,
        filename=disu_path.name,
    )
    disu.load(strict=True)
    verts: np.recarray = pd.DataFrame(gwf.disu.vertices.get_data()[["xv", "yv"]]).to_numpy()
    cell2d: np.recarray = gwf.disu.cell2d.get_data()
    xcyc = pd.DataFrame(cell2d[["xc", "yc"]]).to_numpy()
    iverts = rows_truncate_at_first_missing(cell2d[list(cell2d.dtype.names[4:])])

    clockwise_iverts = []
    for idxs in iverts:
        vertices = verts[idxs]
        if signed_area(vertices) > 0:
            idxs = idxs[::-1]
        clockwise_iverts.append(idxs)

    return verts, clockwise_iverts, xcyc


def densify_poly(
    polygon: shp.Polygon = None,
    distance_between: int | float = None,
) -> shp.Polygon:
    """
    Add evenly spaced points along a polygon exterior.
    """
    exterior: shp.geometry.polygon.LinearRing = polygon.exterior
    total_length = exterior.length
    current_distance = 0.0
    new_points = []

    while current_distance < total_length:
        point = exterior.interpolate(current_distance)
        new_points.append(point)
        current_distance += distance_between

    new_points.append(exterior.interpolate(total_length))

    return shp.Polygon(new_points)


def get_griddata_from_gsf(gsf_path: Path, layer: int = 0, decimals: int = 6):
    """Read a MODFLOW-USG ``.gsf`` grid specification file and return plan arrays.

    A ``.gsf`` is the companion file that carries the geometry a MODFLOW-USG
    ``DISU`` does not: DISU stores connectivity and geometric *measures* (areas,
    connection lengths, face areas) but no coordinates at all, so the ``.gsf`` is
    the only place a USG model's cell outlines live.

    It describes **3-D** cells -- typically eight vertices per hexahedron, the
    bottom four sitting directly under the top four -- and lists **every layer**
    separately. Feeding those records to a plan-view grid as they stand produces
    polygons that trace their own outline twice (invalid, and double the true
    area) and one stacked copy of the grid per layer. So this collapses each cell
    to its distinct ``(x, y)`` corners, keeping file order, and returns a single
    layer.

    Vertices are then re-indexed into a compact array of unique ``(x, y)``
    positions -- a .gsf repeats each plan vertex once per elevation level, so the
    raw table is several times larger than the plan grid needs.

    Rings come back clockwise, matching :func:`get_griddata_from_disu`.

    Parameters
    ----------
    gsf_path
        Path to the ``.gsf`` file.
    layer
        Zero-based layer to take the plan view from, as everywhere else in
        myflopy (``0`` is the model's first layer). Every layer of a USG grid
        shares one plan geometry, so this matters only for the cell centers and
        for which node numbers the cells correspond to.
    decimals
        Rounding used when deciding that two vertices are the same ``(x, y)``.
        The default is far finer than any projected coordinate needs and exists
        for files whose top and bottom vertices are written with different digits.

    Returns
    -------
    tuple
        ``(verts, iverts, xcyc)`` -- the same three arrays
        :func:`get_griddata_from_disu` returns, ready for
        :class:`~myflopy.modflow.mf6.grid.voronoi.VoronoiGridPlus`.
    """

    gsf_path = Path(gsf_path)
    lines = [
        line.strip()
        for line in gsf_path.read_text(errors="replace").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not lines:
        raise ValueError(f"{gsf_path.name} is empty")

    header = lines[0].upper().split()
    if not header or header[0] != "UNSTRUCTURED":
        raise ValueError(
            f"{gsf_path.name} is not a grid specification file: expected a first "
            f"line of 'UNSTRUCTURED' (optionally 'UNSTRUCTURED GWF'), got {lines[0]!r}"
        )

    nnodes = int(lines[1].split()[0])
    nverts = int(lines[2].split()[0])
    vert_lines = lines[3 : 3 + nverts]
    node_lines = lines[3 + nverts : 3 + nverts + nnodes]
    if len(vert_lines) != nverts or len(node_lines) != nnodes:
        raise ValueError(
            f"{gsf_path.name} declares {nverts} vertices and {nnodes} nodes but "
            f"provides {len(vert_lines)} and {len(node_lines)}"
        )

    # x, y only: the z of a .gsf vertex distinguishes the top and bottom copies
    # of one plan position, which is exactly what the plan view collapses.
    xy = np.array([[float(v) for v in line.split()[:2]] for line in vert_lines])

    records = []
    for row, line in enumerate(node_lines):
        parts = line.split()
        declared = int(parts[5])
        provided = len(parts) - 6
        if declared != provided:
            raise ValueError(
                f"{gsf_path.name}: node {row + 1} declares {declared} vertices "
                f"but provides {provided}"
            )
        records.append(
            (
                int(parts[4]),                                   # layer, 1-based
                float(parts[1]),                                 # xc
                float(parts[2]),                                 # yc
                [int(v) - 1 for v in parts[6 : 6 + declared]],   # vertex indices
            )
        )

    layers = sorted({record[0] for record in records})
    if not 0 <= layer < len(layers):
        raise ValueError(
            f"layer must be between 0 and {len(layers) - 1} "
            f"({gsf_path.name} has {len(layers)} layers); got {layer}"
        )
    wanted = layers[layer]

    # One pass, one dict: `compact` maps a rounded (x, y) to its index in
    # `positions`, which IS the vertex array being built. Re-deriving the array
    # inside the loop instead would make this quadratic -- 9405 cells against a
    # vertex table that grows to the same order of magnitude.
    compact: dict[tuple[float, float], int] = {}
    positions: list[tuple[float, float]] = []
    iverts: list[list[int]] = []
    centers: list[list[float]] = []
    collapsed = 0
    for lay, xc, yc, cell_verts in records:
        if lay != wanted:
            continue
        ring: list[int] = []
        seen: set[tuple[float, float]] = set()
        for index in cell_verts:
            key = (round(float(xy[index][0]), decimals), round(float(xy[index][1]), decimals))
            if key in seen:
                continue
            seen.add(key)
            if key not in compact:
                compact[key] = len(positions)
                positions.append(key)
            ring.append(compact[key])
        collapsed += len(cell_verts) - len(ring)
        if signed_area(np.array([positions[i] for i in ring], dtype=float)) > 0:
            ring = ring[::-1]
        iverts.append(ring)
        centers.append([xc, yc])

    verts = np.array(positions, dtype=float)
    logger.debug(
        "%s: layer %d -> %d cells, %d vertices (%d duplicate 3-D vertices collapsed)",
        gsf_path.name,
        layer,
        len(iverts),
        len(verts),
        collapsed,
    )
    return verts, iverts, np.array(centers, dtype=float)
