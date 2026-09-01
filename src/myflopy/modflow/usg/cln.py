"""Read a MODFLOW-USG ``CLN`` (Connected Linear Network) package.

MODFLOW 6 has no CLN, so nothing here converts. What it does is make the
omission *legible*: a CLN network is read, segmented into its physically
distinct features, and reported with the groundwater cells each one touches --
so a user can see exactly what a CLN-deferred conversion left behind, and can
later rebuild those features as ``LAK`` / ``SFR`` on whatever grid they end up
with.

Features are split by **name** when the file labels its nodes, and by graph
shape otherwise. The names are authoritative and the graph is not: on the Ten
Trails network, connected components merge ``CrispCreek`` with its tributary
``Wlnd217`` into one feature, because they meet -- so the graph finds six
features where the file names seven. Where names are absent the shape still
separates the two things CLN is used for: a stream is a chain, whose interior
nodes have two neighbours, while a water body discretized as a 2-D mesh has
five or six.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from myflopy._logging import get_logger
from myflopy.modflow.usg._io import ArrayCursor, NameFile, free_ints, free_row

logger = get_logger(__name__)

__all__ = ["ClnData", "ClnFeature", "read_cln"]

#: Mean node degree at or above which a component is treated as a 2-D water
#: body rather than a stream chain. A chain's interior degree is 2; a triangular
#: or quadrilateral mesh runs 4-6.
_MESH_DEGREE = 3.5


@dataclass(slots=True)
class ClnFeature:
    """One feature of a CLN network -- a named waterbody or stream."""

    index: int
    kind: str
    nodes: np.ndarray
    gwf_nodes: np.ndarray
    layers: np.ndarray
    cells: np.ndarray
    elevations: np.ndarray
    lengths: np.ndarray
    mean_degree: float
    name: str | None = None
    fskin: np.ndarray | None = None
    conduit_types: np.ndarray | None = None

    @property
    def label(self) -> str:
        """The feature's name, or a positional stand-in when the file names none."""

        return self.name or f"{self.kind}_{self.index}"

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
    #: ``FSKIN`` per CLN-GWF connection -- the bed hydraulic conductivity the
    #: connection conductance is built from, and the only calibrated quantity
    #: the network carries. MODFLOW 6's ``bedleak`` is a *leakance* (K over a
    #: bed thickness, 1/T); this is a K, so converting needs a thickness.
    fskin: np.ndarray | None = None
    #: ``FLENGW`` per connection -- the length the connection acts over.
    conn_lengths: np.ndarray | None = None
    #: ``IFTYP`` per node: which row of the conduit table the node uses.
    conduit_types: np.ndarray | None = None
    #: ``FRAD`` per conduit type. A stream's SFR width is twice this.
    radii: np.ndarray | None = None
    #: ``CONDUITK`` per conduit type -- in-pipe conductivity. MODFLOW 6 has no
    #: counterpart: LAK is well mixed and SFR routes by Manning's equation.
    conductivities: np.ndarray | None = None
    #: The trailing label on each node's row, when the file writes one.
    names: np.ndarray | None = None

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

    # Node properties: IFNO IFTYP IFDIR FLENG FELEV FANGLE IFLIN ICCWADI [label]
    lengths = np.zeros(nclnnds)
    elevations = np.zeros(nclnnds)
    conduit_types = np.zeros(nclnnds, dtype=np.int64)
    labels: list[str | None] = []
    for i in range(nclnnds):
        values, label = free_row(cursor.next_line())
        conduit_types[i] = int(values[1])
        lengths[i] = values[3]
        elevations[i] = values[4]
        labels.append(label)

    # CLN-GWF connections: IFNO IGWNOD IFCON FSKIN FLENGW FANISO ICGWADI [label]
    gwf_nodes = np.zeros(nclngwc, dtype=np.int64)
    fskin = np.zeros(nclngwc)
    conn_lengths = np.zeros(nclngwc)
    for i in range(nclngwc):
        values, _ = free_row(cursor.next_line())
        gwf_nodes[i] = int(values[1])
        fskin[i] = values[3] if len(values) > 3 else np.nan
        conn_lengths[i] = values[4] if len(values) > 4 else np.nan

    radii, conductivities = _read_conduit_types(cursor, nconduityp)
    names = np.array(labels, dtype=object) if any(x is not None for x in labels) else None
    features = _segment(
        iac,
        ja,
        gwf_nodes,
        elevations,
        lengths,
        ncpl=ncpl,
        names=names,
        fskin=fskin,
        conduit_types=conduit_types,
    )
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
        fskin=fskin,
        conn_lengths=conn_lengths,
        conduit_types=conduit_types,
        radii=radii,
        conductivities=conductivities,
        names=names,
    )


def _read_conduit_types(cursor: ArrayCursor, nconduityp: int):
    """Read the circular-conduit table: ``ICONDUITYP FRAD CONDUITK`` per type.

    Returns ``(None, None)`` when the network declares no conduit types, or when
    the table is not where it should be -- the geometry it carries is a CLN
    idealization with no MODFLOW 6 counterpart, so failing to find it must not
    cost the caller the rest of the network.
    """

    if nconduityp <= 0:
        return None, None
    radii = np.full(nconduityp, np.nan)
    conductivities = np.full(nconduityp, np.nan)
    for i in range(nconduityp):
        try:
            values, _ = free_row(cursor.next_line())
        except (StopIteration, ValueError, OSError) as error:
            logger.debug(
                "CLN conduit table stops after %d of %d types; radii and in-pipe K "
                "will be missing (neither converts to MF6): %s",
                i,
                nconduityp,
                error,
            )
            break
        if len(values) >= 3:
            radii[i] = values[1]
            conductivities[i] = values[2]
    return radii, conductivities


def _segment(
    iac: np.ndarray,
    ja: np.ndarray,
    gwf_nodes: np.ndarray,
    elevations: np.ndarray,
    lengths: np.ndarray,
    *,
    ncpl: int,
    names: np.ndarray | None = None,
    fskin: np.ndarray | None = None,
    conduit_types: np.ndarray | None = None,
) -> tuple[ClnFeature, ...]:
    """Split the CLN network into features and classify each.

    Groups by the file's own node labels when it writes them, and by connected
    component otherwise. The distinction matters: components merge features that
    touch, so a tributary is swallowed by the stream it joins.
    """

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

    if names is not None:
        # One group per distinct label, ordered largest first to match the
        # component ordering a nameless file would produce.
        groups = [(str(label), np.flatnonzero(names == label)) for label in dict.fromkeys(names)]
        groups.sort(key=lambda item: item[1].size, reverse=True)
        logger.debug(
            "CLN: %d named feature(s) vs %d connected component(s)", len(groups), len(components)
        )
    else:
        groups = [(None, np.asarray(component)) for component in components]

    features: list[ClnFeature] = []
    for index, (label, members) in enumerate(groups):
        if members.size == 0:
            continue
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
                name=label,
                fskin=None if fskin is None or members.max() >= fskin.size else fskin[members],
                conduit_types=None if conduit_types is None else conduit_types[members],
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
