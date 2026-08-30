"""Read a MODFLOW-USG ``CLN`` (Connected Linear Network) package.

MODFLOW 6 has no CLN, so nothing here converts. What it does is make the
omission *legible*: a CLN network is read, segmented into its physically
distinct features, and reported with the groundwater cells each one touches --
so a user can see exactly what a CLN-deferred conversion left behind, and can
later rebuild those features as ``LAK`` / ``SFR`` on whatever grid they end up
with.

The segmentation is by graph shape, because that is what actually distinguishes
the two things CLN is used for. A stream is a chain: interior nodes have two
neighbours. A water body discretized as a 2-D CLN mesh has five or six. Taking
connected components and looking at mean degree separates them without needing
to be told which is which.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from myflopy._logging import get_logger
from myflopy.modflow.usg._io import ArrayCursor, NameFile, free_floats, free_ints

logger = get_logger(__name__)

__all__ = ["ClnData", "ClnFeature", "read_cln"]

#: Mean node degree at or above which a component is treated as a 2-D water
#: body rather than a stream chain. A chain's interior degree is 2; a triangular
#: or quadrilateral mesh runs 4-6.
_MESH_DEGREE = 3.5


@dataclass(slots=True)
class ClnFeature:
    """One connected component of a CLN network."""

    index: int
    kind: str
    nodes: np.ndarray
    gwf_nodes: np.ndarray
    layers: np.ndarray
    cells: np.ndarray
    elevations: np.ndarray
    lengths: np.ndarray
    mean_degree: float

    @property
    def n_nodes(self) -> int:
        """Number of CLN nodes in this feature."""

        return int(self.nodes.size)

    @property
    def is_flat(self) -> bool:
        """True when every node shares one elevation -- the signature of a lake."""

        return bool(np.ptp(self.elevations) < 1e-6)

    def describe(self) -> str:
        """One line naming what this feature is and where it sits."""

        elevation = (
            f"{self.elevations[0]:.2f} ft (flat)"
            if self.is_flat
            else f"{self.elevations.min():.2f}-{self.elevations.max():.2f} ft"
        )
        return (
            f"{self.kind:<10} {self.n_nodes:>4} nodes, "
            f"{len(np.unique(self.cells)):>4} layer-{int(self.layers[0]) + 1} cells, "
            f"elevation {elevation}"
        )


@dataclass(slots=True)
class ClnData:
    """A parsed CLN network and the features it decomposes into."""

    ncln: int
    nclnnds: int
    nclngwc: int
    nconduityp: int
    iac: np.ndarray
    ja: np.ndarray
    gwf_nodes: np.ndarray
    lengths: np.ndarray
    elevations: np.ndarray
    features: tuple[ClnFeature, ...] = ()

    @property
    def waterbodies(self) -> tuple[ClnFeature, ...]:
        """Features whose connectivity says 2-D mesh."""

        return tuple(f for f in self.features if f.kind == "waterbody")

    @property
    def streams(self) -> tuple[ClnFeature, ...]:
        """Features whose connectivity says linear chain."""

        return tuple(f for f in self.features if f.kind == "stream")

    def cells_touched(self, ncpl: int) -> np.ndarray:
        """Every distinct plan cell (0-based) the network connects to."""

        return np.unique((self.gwf_nodes - 1) % ncpl)


def read_cln(name_file: NameFile, *, ncpl: int) -> ClnData | None:
    """Read the ``CLN`` package, or return ``None`` when the model has none.

    Parameters
    ----------
    name_file
        The parsed name file.
    ncpl
        Nodes per layer of the groundwater grid, used to split a CLN's
        groundwater node numbers into layer and plan cell.
    """

    entry = name_file.package("CLN")
    if entry is None:
        return None

    cursor = ArrayCursor(entry.path, name_file)
    header = free_ints(cursor.next_line())
    ncln, nclnnds = header[0], header[1]
    nclngwc = header[6] if len(header) > 6 else nclnnds
    nconduityp = header[7] if len(header) > 7 else 0

    if nclnnds <= 0:
        logger.warning("CLN declares %d nodes; nothing to read", nclnnds)
        return None

    njacln = free_ints(cursor.next_line())[0]
    iac = cursor.read_array(nclnnds, int)
    ja = cursor.read_array(njacln, int)

    # Node properties: IFNO IFTYP IFDIR FLENG FELEV FANGLE IFLIN ICCWADI
    lengths = np.zeros(nclnnds)
    elevations = np.zeros(nclnnds)
    for i in range(nclnnds):
        values = free_floats(cursor.next_line())
        lengths[i] = values[3]
        elevations[i] = values[4]

    # CLN-GWF connections: IFNO IGWNOD IFCON FSKIN FLENG FANISO ICGWADI
    gwf_nodes = np.zeros(nclngwc, dtype=np.int64)
    for i in range(nclngwc):
        values = free_floats(cursor.next_line())
        gwf_nodes[i] = int(values[1])

    features = _segment(iac, ja, gwf_nodes, elevations, lengths, ncpl=ncpl)
    logger.debug(
        "CLN: %d nodes, %d features (%d waterbody, %d stream)",
        nclnnds,
        len(features),
        sum(f.kind == "waterbody" for f in features),
        sum(f.kind == "stream" for f in features),
    )
    return ClnData(
        ncln=ncln,
        nclnnds=nclnnds,
        nclngwc=nclngwc,
        nconduityp=nconduityp,
        iac=iac,
        ja=ja,
        gwf_nodes=gwf_nodes,
        lengths=lengths,
        elevations=elevations,
        features=features,
    )


def _segment(
    iac: np.ndarray,
    ja: np.ndarray,
    gwf_nodes: np.ndarray,
    elevations: np.ndarray,
    lengths: np.ndarray,
    *,
    ncpl: int,
) -> tuple[ClnFeature, ...]:
    """Split the CLN graph into connected components and classify each."""

    n = iac.size
    pointer = np.concatenate([[0], np.cumsum(iac)])
    adjacency: dict[int, list[int]] = {i: [] for i in range(n)}
    for node in range(n):
        for neighbour in ja[pointer[node] + 1 : pointer[node + 1]]:
            index = abs(int(neighbour)) - 1
            if 0 <= index < n and index != node:
                adjacency[node].append(index)

    degree = iac - 1
    seen: set[int] = set()
    components: list[list[int]] = []
    for start in range(n):
        if start in seen:
            continue
        stack, component = [start], []
        seen.add(start)
        while stack:
            node = stack.pop()
            component.append(node)
            for neighbour in adjacency[node]:
                if neighbour not in seen:
                    seen.add(neighbour)
                    stack.append(neighbour)
        components.append(sorted(component))

    components.sort(key=len, reverse=True)
    features: list[ClnFeature] = []
    for index, component in enumerate(components):
        members = np.asarray(component)
        mean_degree = float(degree[members].mean())
        connected = gwf_nodes[members] if members.max() < gwf_nodes.size else np.array([], int)
        features.append(
            ClnFeature(
                index=index,
                kind="waterbody" if mean_degree >= _MESH_DEGREE else "stream",
                nodes=members + 1,
                gwf_nodes=connected,
                layers=(connected - 1) // ncpl,
                cells=(connected - 1) % ncpl,
                elevations=elevations[members],
                lengths=lengths[members],
                mean_degree=mean_degree,
            )
        )
    return tuple(features)


def cln_polygons(cln: ClnData, grid, *, crs: str | None = None):
    """Dissolve each CLN feature into one polygon on the groundwater grid.

    This is the shape worth keeping: it survives a change of grid, which the
    node numbers do not. Rebuilding these features as ``LAK``/``SFR`` on a new
    mesh starts from these polygons.

    Parameters
    ----------
    cln
        A parsed network.
    grid
        A ``VoronoiGridPlus`` (anything exposing ``gdf_vorPolys``).
    crs
        Coordinate reference system for the result; defaults to the grid's.

    Returns
    -------
    geopandas.GeoDataFrame
        One row per feature, with ``feature``, ``kind``, ``n_nodes``,
        ``elev_min``, ``elev_max`` and the dissolved geometry.
    """

    import geopandas as gpd

    cells = grid.gdf_vorPolys
    rows = []
    for feature in cln.features:
        if feature.cells.size == 0:
            continue
        merged = cells.iloc[np.unique(feature.cells)].union_all()
        rows.append(
            {
                "feature": feature.index,
                "kind": feature.kind,
                "n_nodes": feature.n_nodes,
                "n_cells": int(np.unique(feature.cells).size),
                "elev_min": float(feature.elevations.min()),
                "elev_max": float(feature.elevations.max()),
                "geometry": merged,
            }
        )
    return gpd.GeoDataFrame(rows, crs=crs or getattr(cells, "crs", None))


def write_cln_report(cln: ClnData, path: str | Path) -> Path:
    """Write the CLN inventory to a text file, for the record."""

    path = Path(path)
    lines = [f"CLN network: {cln.nclnnds} nodes, {len(cln.features)} features", ""]
    lines += [f"  [{f.index}] {f.describe()}" for f in cln.features]
    path.write_text("\n".join(lines) + "\n")
    return path
