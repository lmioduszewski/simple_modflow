from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import shapely as shp
from flopy.mf6 import MFSimulation, ModflowGwf, ModflowGwfdisu


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
