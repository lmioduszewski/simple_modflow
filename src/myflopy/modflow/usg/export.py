"""Extract a MODFLOW-USG model's inputs as grid-independent GIS files.

A converted model is tied to the grid it was converted on. That is the wrong
unit of reuse for the usual reason a USG model gets imported at all: the point
is to rebuild it on a better mesh, and the mesh will change again. So this
writes the model's *content* -- geometry, properties and forcings -- in forms
that outlive any particular discretization, and leaves building packages to the
caller.

Three ideas run through it.

**A boundary condition here is a line, not a patch.** Measured on Ten Trails,
every list BC is a chain of cells one cell wide: mean in-set neighbour counts of
1.87 (DRN), 1.94 (GHB) and 2.52 (CHD), against 4-6 for anything areal. So they
are written as line segments, one per source cell, which partition the feature
exactly and re-resolve onto any mesh.

**Conductance is the field that breaks.** Split a drain cell in two and the
elevation is unchanged while the conductance is not. Every conductance column is
therefore written beside a ``*_per_ft`` twin, so a new discretization recovers
it as ``value x length in the cell`` and the feature total is preserved. What
the source conductance is *not* is a leakance: measured over 205 DRN records,
``corr(cond, cell area) ~ 0`` and GHB is one constant per family regardless of a
249,000-fold spread in cell area. These are calibrated numbers, so the feature
total is the only thing worth conserving.

**A repeated period is written once.** The Ten Trails run is 72 periods and a
strict 12-month cycle, so the arrays are written per cycle-month with an index
table mapping every period onto one, rather than 72 near-identical rasters.

Rasters are the convenience form and lose something: a Voronoi mesh whose cells
span 2.8 to 697,010 ft2 cannot be resampled onto a regular grid without absorbing
the smallest cells. ``cell_values.gpkg`` carries every array again as points at
the cell centres, which is exact -- use it when the raster's resolution matters.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

from myflopy._logging import get_logger

if TYPE_CHECKING:  # pragma: no cover - typing only
    from myflopy.modflow.usg.model import UsgModel

logger = get_logger(__name__)

__all__ = ["ExportManifest", "export_gis"]

#: Value columns of each list BC, in the order the reader stores them.
_BC_VALUES = {
    "CHD": ("shead", "ehead"),
    "DRN": ("elev", "cond"),
    "GHB": ("bhead", "cond"),
    "RIV": ("stage", "cond", "rbot"),
    "WEL": ("q",),
}

#: Columns that are extensive -- they scale with the size of the piece they sit
#: on, so a change of mesh must redistribute rather than copy them.
_EXTENSIVE = {"cond", "q"}

#: USG's field name -> the name the myflopy package registry already uses for it
#: (``_PACKAGE_EXPLORER_SPECS[pkg].gpkg_defaults``). Writing the registry's names
#: is what lets ``mf.drn.gpkg(path, layer="drn", context=ctx, nper=nper)`` read the
#: export with no field arguments at all -- the whole point of the export is that
#: its output is ordinary input to the canonical package-first API.
_CANONICAL = {
    "elev": "elevation",
    "cond": "conductance",
    "bhead": "head",
    "shead": "head",
    "stage": "stage",
    "rbot": "rbot",
    "q": "q",
}

#: Bed thickness assumed when turning a USG ``FSKIN`` (a hydraulic conductivity,
#: L/T) into a MODFLOW 6 ``bedleak`` (a leakance, 1/T). CLN carries no bed
#: thickness -- the conduit has a skin, not a bed -- so this is a choice, and it
#: is recorded in the manifest so it can be changed with one multiplication.
_ASSUMED_BED_THICKNESS = 1.0


@dataclass(slots=True)
class ExportManifest:
    """What an export wrote, and what each piece means."""

    directory: Path
    crs: str
    nper: int
    layers: dict[str, dict[str, Any]] = field(default_factory=dict)
    rasters: dict[str, dict[str, Any]] = field(default_factory=dict)
    tables: dict[str, dict[str, Any]] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    def to_json(self) -> str:
        """Return the manifest as indented JSON."""

        return json.dumps(
            {
                "crs": self.crs,
                "nper": self.nper,
                "vector_layers": self.layers,
                "rasters": self.rasters,
                "tables": self.tables,
                "notes": self.notes,
            },
            indent=2,
            default=str,
        )

    def describe(self) -> str:
        """Return a plain-text summary of the export."""

        lines = [f"USG export -> {self.directory}", f"  crs {self.crs}, {self.nper} stress periods", ""]
        for title, entries in (
            ("vector layers", self.layers),
            ("rasters", self.rasters),
            ("tables", self.tables),
        ):
            if not entries:
                continue
            lines.append(f"  {title}:")
            for name, info in entries.items():
                count = info.get("features", info.get("rows", info.get("files", "")))
                lines.append(f"    {name:<28} {count!s:>8}  {info.get('what', '')}")
            lines.append("")
        if self.notes:
            lines.append("  notes:")
            lines += [f"    - {note}" for note in self.notes]
        return "\n".join(lines)


def export_gis(
    model: UsgModel,
    directory: str | Path,
    *,
    crs: str | None = None,
    stream_lines: dict[str, str | Path] | None = None,
    clip_to: str | Path | None = None,
    resolution: float | None = None,
    bed_thickness: float = _ASSUMED_BED_THICKNESS,
    overwrite: bool = False,
) -> ExportManifest:
    """Write a USG model's inputs as grid-independent vector and raster files.

    Parameters
    ----------
    model
        A model read by :func:`~myflopy.modflow.usg.reader.read_usg`, with a grid.
    directory
        Where to write. Created if absent.
    crs
        Coordinate reference system for everything written. Defaults to the
        grid's own. A ``.gsf`` carries no CRS, so check that the grid's is right
        before relying on it -- the wrong one is silent and puts the export
        nowhere near the model it is meant to rebuild.
    stream_lines
        ``{stream name: path}`` mapping a CLN stream onto a centerline you would
        rather use. When given, each node's ``station`` is measured along *your*
        line, so interpolating onto your own reaches needs no reprojection. The
        CLN's own centerline is written either way.
    clip_to
        A polygon file; supplied ``stream_lines`` are clipped to it first.
    resolution
        Raster cell size, in the grid's units. Defaults to a quarter of the
        median cell width, which resolves a typical cell and absorbs the
        smallest ones -- see the module docstring.
    bed_thickness
        Thickness assumed when deriving ``bedleak`` from ``FSKIN``.
    overwrite
        Replace existing files rather than refusing.

    Returns
    -------
    ExportManifest
    """

    if model.grid is None:
        raise ValueError(
            "export_gis() needs the grid; read the model with gsf= so it has geometry."
        )
    directory = Path(directory)
    if directory.exists() and any(directory.iterdir()) and not overwrite:
        raise FileExistsError(f"{directory} is not empty; pass overwrite=True to replace it.")
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "arrays").mkdir(exist_ok=True)

    target_crs = crs or _grid_crs(model)
    manifest = ExportManifest(directory=directory, crs=str(target_crs), nper=model.nper)

    _export_boundaries(model, directory, target_crs, manifest)
    _export_cln(model, directory, target_crs, manifest, stream_lines, clip_to, bed_thickness)
    _export_arrays(model, directory, target_crs, manifest, resolution)
    _export_tables(model, directory, manifest)
    (directory / "manifest.json").write_text(manifest.to_json())
    (directory / "README.md").write_text(_readme(model, manifest))

    logger.info("wrote USG export to %s (%s)", directory, target_crs)
    return manifest


# -- geometry helpers ------------------------------------------------------


def _grid_crs(model: UsgModel) -> str:
    """The grid's CRS as a string, or a refusal naming what is missing."""

    crs = getattr(model.grid.gdf_vorPolys, "crs", None)
    if crs is None:
        raise ValueError(
            "the grid has no CRS and none was given. A .gsf carries no coordinate "
            "system, so pass crs= to export_gis() (or to read_usg())."
        )
    return str(crs)


def _adjacency(cells: list[int], polys) -> dict[int, set[int]]:
    """Which of ``cells`` touch which, as an adjacency map."""

    subset = polys.iloc[cells]
    touch = subset.sindex.query(subset.geometry, predicate="touches")
    adjacent: dict[int, set[int]] = {c: set() for c in cells}
    for a, b in zip(*touch, strict=False):
        adjacent[cells[a]].add(cells[b])
    return adjacent


def _chains(cells: list[int], polys) -> list[list[int]]:
    """Order a set of touching cells into open chains, one per component."""

    adjacent = _adjacency(cells, polys)
    members = set(cells)
    chains, seen = [], set()
    for start in cells:
        if start in seen:
            continue
        component, stack = [], [start]
        seen.add(start)
        while stack:
            node = stack.pop()
            component.append(node)
            for other in adjacent[node]:
                if other not in seen:
                    seen.add(other)
                    stack.append(other)
        ends = [c for c in component if len(adjacent[c] & members) <= 1] or [component[0]]
        walk, current, used = [ends[0]], ends[0], {ends[0]}
        while True:
            nxt = [o for o in adjacent[current] if o in members and o not in used]
            if not nxt:
                break
            current = nxt[0]
            used.add(current)
            walk.append(current)
        # A branching feature leaves cells unvisited; keep them as their own run
        # rather than dropping them silently.
        chains.append(walk)
        leftover = [c for c in component if c not in used]
        if leftover:
            logger.debug("chain from cell %d left %d branch cell(s) over", start, len(leftover))
            chains.extend(_chains(leftover, polys) if len(leftover) > 1 else [leftover])
    return chains


def _segment_line(points: list[Any], index: int):
    """The piece of a centroid polyline belonging to vertex ``index``.

    Runs from the midpoint of the previous span to the midpoint of the next, so
    consecutive segments tile the line exactly and their lengths sum to it.
    """

    from shapely.geometry import LineString, Point

    here = points[index]
    coords = []
    if index > 0:
        before = points[index - 1]
        coords.append(Point((before.x + here.x) / 2, (before.y + here.y) / 2))
    coords.append(here)
    if index < len(points) - 1:
        after = points[index + 1]
        coords.append(Point((here.x + after.x) / 2, (here.y + after.y) / 2))
    if len(coords) == 1:  # a lone cell -- give it a degenerate but valid line
        coords = [here, here]
    return LineString([(p.x, p.y) for p in coords])


# -- list boundary conditions ----------------------------------------------


def _export_boundaries(model: UsgModel, directory: Path, crs: str, manifest: ExportManifest) -> None:
    """Write DRN/GHB/CHD/RIV as canonical point inputs, plus reference lines and the mesh.

    The **point** layer is the one to build from, and it is named for the package
    (``drn``, ``ghb``, ``chd``) with the column names the package registry already
    uses, so it is read by ``mf.drn.gpkg(path, layer="drn", context=ctx, nper=nper)``
    with no field arguments.

    Build from the ``*_lines`` layer. A line covers the same continuous run of cells
    the original model occupied, which a point layer cannot: points resolve to
    exactly the cells that contain them, so on a finer mesh the boundary comes back
    as a dotted line instead of a feature. Conductance is split along the line with
    :class:`~myflopy.geopackage.Spread`, so the total survives the change of mesh.

    The ``*`` point layers carry the same values one per source cell. They resolve to
    one cell each, so an unwrapped conductance is already exact there -- useful as a
    check, and as the fallback if a boundary must not spread at all.

    Why the split is needed at all:
    ``GeoPackageSource._boundary_data`` writes a feature's values to *every* cell the
    feature intersects, which is right for an elevation and wrong for a conductance:
    a line spanning two cells doubles it. Measured on this model's DRN, through
    ``mf.drn.gpkg`` on the *same* grid the model came from -- 205 segments became 424
    records and 66,007.43 ft2/d became 133,631.25, a factor of 2.02. As points the
    same call returns 205 records and 66,007.43, exactly.

    The ``*_lines`` layers carry the same values plus ``conductance_per_ft`` and are
    written for reading and for a future length-splitting resolver; they are not the
    build path today.
    """

    import geopandas as gpd

    path = directory / "boundaries.gpkg"
    polys = model.grid.gdf_vorPolys
    written = 0

    for ftype in ("DRN", "GHB", "CHD", "RIV"):
        records = model.boundaries.get(ftype)
        if records is None or not records.periods:
            continue
        points, lines = _bc_frames(model, ftype, crs)
        if points is None or points.empty:
            continue
        package = ftype.lower()
        value_columns = [
            c for c in points.columns
            if c not in ("geometry", "name", "layer", "old_cell", "length_ft")
        ]
        points.to_file(path, layer=package, driver="GPKG")
        manifest.layers[f"boundaries.gpkg:{package}"] = {
            "features": int(len(points)),
            "what": f"{ftype} as one point per source cell -- exact conductance, discontinuous",
            "consumed_by": f"mf.{package}.gpkg(path, layer='{package}', context=ctx, nper=nper)",
            "columns": value_columns,
            "layer_base": 1,
            "geometry": "Point",
        }
        lines.to_file(path, layer=f"{package}_lines", driver="GPKG")
        manifest.layers[f"boundaries.gpkg:{package}_lines"] = {
            "features": int(len(lines)),
            "what": f"{ftype} as lines -- BUILD FROM THIS, continuous cell coverage",
            "consumed_by": (
                f"mf.{package}.gpkg(path, layer='{package}_lines', context=ctx, nper=nper, "
                f"conductance=mf.Spread('conductance'))"
            ),
            "columns": [c for c in lines.columns if c != "geometry"],
            "geometry": "LineString",
        }
        written += 1

    cells = gpd.GeoDataFrame(
        {
            "cell": np.arange(model.ncpl),
            "area_ft2": polys.geometry.area.to_numpy(),
            "geometry": polys.geometry.to_numpy(),
        },
        crs=polys.crs,
    ).to_crs(crs)
    cells.to_file(path, layer="source_cells", driver="GPKG")
    manifest.layers["boundaries.gpkg:source_cells"] = {
        "features": int(len(cells)),
        "what": "the old Voronoi mesh, so every value above is auditable",
        "geometry": "Polygon",
    }

    active = np.asarray(model.idomain).max(axis=0) > 0
    domain = gpd.GeoDataFrame(
        {"what": ["active domain"], "geometry": [polys.iloc[np.flatnonzero(active)].union_all()]},
        crs=polys.crs,
    ).to_crs(crs)
    domain.to_file(path, layer="domain", driver="GPKG")
    manifest.layers["boundaries.gpkg:domain"] = {
        "features": 1,
        "what": "outline of cells active in any layer",
        "geometry": "Polygon",
    }
    logger.debug("boundaries.gpkg: %d BC package(s) + source_cells + domain", written)


def _bc_frames(model: UsgModel, ftype: str, crs: str):
    """One list BC as ``(points, lines)``, carrying the package registry's field names.

    Both frames hold one row per source cell and the same values; they differ only
    in geometry and in the ``*_per_ft`` twins, which are meaningless on a point.
    """

    import geopandas as gpd

    frame = model.boundary_frame(ftype)
    if frame.empty:
        return None, None
    polys = model.grid.gdf_vorPolys
    columns = [c for c in _BC_VALUES.get(ftype, ()) if c in frame.columns]
    periods = sorted(frame["period"].unique())
    cyclic = _cycle_index({int(p): frame[frame["period"] == p][columns[0]].to_numpy()
                           for p in periods}) if len(periods) > 1 else None

    package = ftype.lower()
    # CHD's second USG field is EHEAD, the end-of-period ramp; MF6 has one HEAD, so
    # only SHEAD becomes `head` and EHEAD rides along under its own name (see report()).
    canonical = {c: _CANONICAL.get(c, c) for c in columns}
    if ftype == "CHD":
        canonical["ehead"] = "ehead_usg_only"

    point_rows, line_rows = [], []
    first = frame[frame["period"] == periods[0]]
    for layer, block in first.groupby("layer"):
        values = block.set_index("cell")
        cells = list(np.unique(block["cell"].to_numpy()))
        for chain_id, walk in enumerate(_chains(cells, polys)):
            points = [polys.geometry.iloc[c].centroid for c in walk]
            for index, cell in enumerate(walk):
                line = _segment_line(points, index)
                shared = {
                    "name": f"{package}_L{int(layer) + 1}_c{int(cell)}",
                    "layer": int(layer) + 1,
                    "old_cell": int(cell),
                    "length_ft": float(line.length),
                }
                for column in columns:
                    shared[canonical[column]] = float(values.loc[cell, column])
                # The elevation-like field, expressed as a signed offset from the OLD
                # model top. Absolute elevations are tied to the old layering and land
                # below their new cell bottom when the layering changes; the offset is
                # grid-independent, and `mf.CellSurfaceOffset("cell_top", offset=...)`
                # reads it directly.
                first = canonical[columns[0]]
                shared[f"{first}_below_top"] = shared[first] - float(model.top[cell])
                point_rows.append({**shared, "geometry": points[index]})
                extra = {
                    f"{canonical[c]}_per_ft": (
                        shared[canonical[c]] / line.length if line.length else np.nan
                    )
                    for c in columns if c in _EXTENSIVE
                }
                line_rows.append({**shared, **extra, "chain": chain_id,
                                  "position": index, "geometry": line})

    out_points = gpd.GeoDataFrame(point_rows, geometry="geometry", crs=polys.crs)
    out_lines = gpd.GeoDataFrame(line_rows, geometry="geometry", crs=polys.crs)
    if cyclic is not None and len(cyclic["distinct"]) > 1:
        out_points = _attach_cycle(out_points, frame, columns[0], cyclic,
                                   name=canonical[columns[0]])
        out_lines = _attach_cycle(out_lines, frame, columns[0], cyclic,
                                  name=canonical[columns[0]])
    return out_points.to_crs(crs), out_lines.to_crs(crs)


def _cycle_index(store: dict[int, np.ndarray]) -> dict[str, Any]:
    """Map each period onto a distinct array, and report the cycle it forms.

    A transient model built from a repeating climate writes the same arrays over
    and over -- Ten Trails is 72 periods and 12 months. Writing one copy per
    distinct array with an index is both smaller and clearer than 72 files whose
    equality the reader has to discover.
    """

    keys = sorted(store)
    distinct: list[np.ndarray] = []
    index: list[int] = []
    for key in keys:
        array = np.asarray(store[key], float).ravel()
        for position, seen in enumerate(distinct):
            if array.shape == seen.shape and np.allclose(array, seen, equal_nan=True):
                index.append(position)
                break
        else:
            distinct.append(array)
            index.append(len(distinct) - 1)
    period_of_cycle = None
    for length in range(1, len(index) // 2 + 1):
        if all(index[i] == index[i % length] for i in range(len(index))):
            period_of_cycle = length
            break
    return {"periods": keys, "index": index, "distinct": distinct, "cycle": period_of_cycle}


def _attach_cycle(out, frame, column: str, cyclic: dict[str, Any], *, name: str | None = None):
    """Add one column per cycle position for a BC whose values move in time.

    Keyed on ``(layer, cell)``, not on the cell alone: a boundary is routinely
    applied to the same plan cell in several layers -- Ten Trails' CHD covers 23
    cells in each of layers 1-3 -- and keying on the cell would silently pick one
    of the three, or raise where a scalar was expected.
    """

    cycle = cyclic["cycle"] or len(cyclic["periods"])
    for position, period in enumerate(cyclic["periods"][:cycle], start=1):
        block = frame[frame["period"] == period]
        if column not in block.columns:
            continue
        lookup = {
            (int(layer) + 1, int(cell)): float(value)
            for layer, cell, value in zip(
                block["layer"], block["cell"], block[column], strict=False
            )
        }
        out[f"{name or column}_m{position:02d}"] = [
            lookup.get((int(layer), int(cell)), np.nan)
            for layer, cell in zip(out["layer"], out["old_cell"], strict=False)
        ]
    return out


# -- CLN: streams and lakes ------------------------------------------------


def _export_cln(
    model: UsgModel,
    directory: Path,
    crs: str,
    manifest: ExportManifest,
    stream_lines: dict[str, str | Path] | None,
    clip_to: str | Path | None,
    bed_thickness: float,
) -> None:
    """Write the CLN network as one file for streams and one for lakes.

    Shaped for ``mf.sfr`` and ``mf.lak`` rather than for MODFLOW-USG: column
    names are the builders' parameter names, so the files can be handed straight
    to them. What CLN cannot supply -- Manning's ``n``, a streambed thickness --
    is written as a documented default rather than omitted, so a missing column
    never becomes a silent zero.
    """

    if model.cln is None:
        manifest.notes.append("no CLN package: no streams or lakes were written")
        return

    streams = [f for f in model.cln.features if f.kind == "stream"]
    lakes = [f for f in model.cln.features if f.kind == "waterbody"]
    if streams:
        _write_streams(model, streams, directory, crs, manifest, stream_lines, clip_to)
    if lakes:
        _write_lakes(model, lakes, directory, crs, manifest, bed_thickness)


def _walk_stream(model: UsgModel, feature) -> np.ndarray:
    """Order one stream's CLN nodes from its high end to its low end."""

    cln = model.cln
    pointer = np.concatenate([[0], np.cumsum(cln.iac)])
    members = set(int(n) - 1 for n in feature.nodes)
    adjacent = {i: [] for i in members}
    for node in members:
        for neighbour in cln.ja[pointer[node] + 1 : pointer[node + 1]]:
            index = abs(int(neighbour)) - 1
            if index in members and index != node:
                adjacent[node].append(index)
    ends = [i for i in members if len(adjacent[i]) == 1] or [next(iter(members))]
    start = max(ends, key=lambda i: cln.elevations[i])
    walk, current, used = [start], start, {start}
    while True:
        nxt = [o for o in adjacent[current] if o not in used]
        if not nxt:
            break
        current = nxt[0]
        used.add(current)
        walk.append(current)
    if len(walk) < len(members):
        logger.debug(
            "%s: walked %d of %d nodes; the rest branch off the main stem",
            feature.label,
            len(walk),
            len(members),
        )
    return np.asarray(walk)


def _write_streams(
    model: UsgModel,
    features,
    directory: Path,
    crs: str,
    manifest: ExportManifest,
    stream_lines: dict[str, str | Path] | None,
    clip_to: str | Path | None,
) -> None:
    """Write ``streams.gpkg``: one line per stream plus a node point per reach."""

    import geopandas as gpd
    from shapely.geometry import LineString, Point

    cln, polys = model.cln, model.grid.gdf_vorPolys
    path = directory / "streams.gpkg"
    reference = _reference_lines(stream_lines, clip_to, polys.crs)
    active = np.asarray(model.idomain).max(axis=0) > 0

    line_rows, node_rows, short = [], [], []
    for feature in features:
        walk = _walk_stream(model, feature)
        cells = (cln.gwf_nodes[walk] - 1) % model.ncpl
        points = [polys.geometry.iloc[c].centroid for c in cells]
        bed = cln.elevations[walk]
        own_station = np.r_[
            0.0, np.cumsum([points[i].distance(points[i + 1]) for i in range(len(points) - 1)])
        ]
        if len(points) < 2:
            # One cell is a point, not a centerline. Its node still goes to
            # stream_nodes -- the bed elevation and FSKIN are real -- but there
            # is no line to write, and inventing one would put a reach somewhere
            # nobody chose.
            logger.warning(
                "%s occupies a single cell, so it has no centerline; it is in "
                "stream_nodes but not in the streams layer",
                feature.label,
            )
            short.append(feature.label)
        theirs = reference.get(feature.label)
        station = (
            _oriented_station(theirs, points, feature.label)
            if theirs is not None
            else own_station
        )
        radius = _radius_of(cln, feature)
        skin = feature.fskin if feature.fskin is not None else np.full(walk.size, np.nan)

        drop = float(bed[0] - bed[-1])
        if len(points) >= 2:
            line_rows.append(
                {
                    "stream_id": feature.label,
                    "n_nodes": int(walk.size),
                    "rwid": 2.0 * radius if np.isfinite(radius) else np.nan,
                    "rhk": float(np.nanmedian(skin)),
                    "rgrd": drop / own_station[-1] if own_station[-1] else np.nan,
                    "rbth": 1.0,
                    "man": 0.03,
                    "rtp_upstream": float(bed[0]),
                    "rtp_downstream": float(bed[-1]),
                    "length_ft": float(own_station[-1]),
                    "reference_line": (
                        str(stream_lines.get(feature.label, "")) if stream_lines else ""
                    ),
                    "geometry": LineString(
                        [(p.x, p.y, float(z)) for p, z in zip(points, bed, strict=False)]
                    ),
                }
            )
        for order, (node, point) in enumerate(zip(walk, points, strict=False)):
            node_rows.append(
                {
                    "stream_id": feature.label,
                    "order": order,
                    "station": float(station[order]),
                    "station_cln": float(own_station[order]),
                    "bed_elev": float(cln.elevations[node]),
                    "fskin": float(skin[order]) if order < skin.size else np.nan,
                    "fleng": float(cln.lengths[node]),
                    "old_cell": int(cells[order]),
                    "old_layer": int((cln.gwf_nodes[walk][order] - 1) // model.ncpl) + 1,
                    "cell_active": bool(active[cells[order]]),
                    "geometry": Point(point.x, point.y),
                }
            )

    lines = gpd.GeoDataFrame(line_rows, geometry="geometry", crs=polys.crs).to_crs(crs)
    nodes = gpd.GeoDataFrame(node_rows, geometry="geometry", crs=polys.crs).to_crs(crs)
    lines.to_file(path, layer="streams", driver="GPKG")
    nodes.to_file(path, layer="stream_nodes", driver="GPKG")

    manifest.layers["streams.gpkg:streams"] = {
        "features": int(len(lines)),
        "what": "one LineString Z per stream; Z is the CLN bed elevation",
        "columns": [c for c in lines.columns if c != "geometry"],
        "for": "mf.sfr(streams=..., stream_id='stream_id', width='rwid', "
        "streambed_k='rhk', gradient='rgrd', streambed_thickness='rbth', roughness='man')",
        "geometry": "LineStringZ",
    }
    manifest.layers["streams.gpkg:stream_nodes"] = {
        "features": int(len(nodes)),
        "what": "one point per old cell the stream crosses; interpolate on 'station'",
        "columns": [c for c in nodes.columns if c != "geometry"],
        "geometry": "Point",
    }
    manifest.notes.append(
        "SFR rbth (streambed thickness) and man (Manning's n) are NOT in CLN -- "
        "CLN routes with a pipe conductivity instead. They are written as 1.0 ft "
        "and 0.03; change them."
    )
    _note_bed_at_top(model, features, manifest)
    _note_inactive_cells(features, active, manifest)
    if short:
        manifest.notes.append(
            "single-cell stream(s) with no centerline to draw: "
            + ", ".join(short)
            + ". Their nodes are in stream_nodes; the streams layer has no row for them."
        )


def _oriented_station(line, points, label: str) -> np.ndarray:
    """Distance of each node along ``line``, measured from its upstream end.

    A centerline drawn independently has no reason to run the way the stream
    flows -- both Ten Trails creeks are digitized mouth-first, so projecting onto
    them raw gives a station that *decreases* downstream and a profile that reads
    as climbing. The nodes are already ordered upstream-first, so the sign of the
    correlation between order and raw station says which end is which.
    """

    raw = np.array([line.project(point) for point in points])
    order = np.arange(raw.size)
    if raw.size > 1 and raw.std() > 0 and np.corrcoef(order, raw)[0, 1] < 0:
        raw = line.length - raw
        logger.debug("%s: reference centerline runs mouth-first; station flipped", label)
    scrambled = int((np.diff(raw) < 0).sum())
    if scrambled:
        logger.warning(
            "%s: %d of %d nodes project out of order onto the supplied centerline -- "
            "they sit where it meanders back on itself, so their station is ambiguous. "
            "station_cln is unaffected.",
            label,
            scrambled,
            raw.size - 1,
        )
    return raw


def _reference_lines(stream_lines, clip_to, crs):
    """Load and clip caller-supplied centerlines, keyed by stream name."""

    if not stream_lines:
        return {}
    import geopandas as gpd
    from shapely.geometry import MultiLineString
    from shapely.ops import linemerge

    boundary = None
    if clip_to is not None:
        boundary = gpd.read_file(clip_to)
        boundary = boundary.to_crs(crs) if boundary.crs is not None else boundary.set_crs(crs)

    out = {}
    for name, source in stream_lines.items():
        frame = gpd.read_file(source)
        frame = frame.set_crs(crs, allow_override=True) if frame.crs is None else frame.to_crs(crs)
        if boundary is not None:
            frame = gpd.clip(frame, boundary)
        if frame.empty:
            logger.warning("reference line for %s is empty after clipping; ignoring it", name)
            continue
        geometry = frame.union_all()
        merged = linemerge(geometry) if isinstance(geometry, MultiLineString) else geometry
        if isinstance(merged, MultiLineString):
            merged = max(merged.geoms, key=lambda part: part.length)
            logger.debug("reference line for %s is in pieces; using the longest", name)
        out[name] = merged
    return out


def _radius_of(cln, feature) -> float:
    """The conduit radius for a feature, or NaN when the table is absent."""

    if cln.radii is None or feature.conduit_types is None or feature.conduit_types.size == 0:
        return float("nan")
    kinds = np.unique(feature.conduit_types)
    index = int(kinds[0]) - 1
    if not 0 <= index < cln.radii.size:
        return float("nan")
    return float(cln.radii[index])


def _write_lakes(
    model: UsgModel,
    features,
    directory: Path,
    crs: str,
    manifest: ExportManifest,
    bed_thickness: float,
) -> None:
    """Write ``lakes.gpkg``: a polygon per lake, its bed nodes, and its forcing.

    The CLN "wells" that carry precipitation minus evaporation over each lake
    become a rate over the lake's area, which is what ``mf.lak`` wants for
    ``rainfall``/``evaporation`` -- the one part of a CLN lake that transfers
    without an assumption.
    """

    import geopandas as gpd
    from shapely.geometry import Point

    cln, polys = model.cln, model.grid.gdf_vorPolys
    path = directory / "lakes.gpkg"
    forcing = _lake_forcing(model, features)
    active = np.asarray(model.idomain).max(axis=0) > 0

    lake_rows, node_rows = [], []
    for feature in features:
        cells = np.unique(feature.cells)
        footprint = polys.iloc[cells].union_all()
        skin = feature.fskin if feature.fskin is not None else np.array([np.nan])
        radius = _radius_of(cln, feature)
        lake_rows.append(
            {
                "lake_id": feature.label,
                "n_nodes": feature.n_nodes,
                "n_old_cells": int(cells.size),
                "area_ft2": float(footprint.area),
                "strt": float(np.median(feature.elevations)),
                "lake_bottom": float(np.median(feature.elevations)),
                "bottom_min": float(feature.elevations.min()),
                "bottom_max": float(feature.elevations.max()),
                "flat_bottom": bool(feature.is_flat),
                "bedleak": float(np.nanmedian(skin)) / bed_thickness,
                "fskin_min": float(np.nanmin(skin)),
                "fskin_median": float(np.nanmedian(skin)),
                "fskin_max": float(np.nanmax(skin)),
                "conduit_radius": radius,
                "old_layer": int(np.unique(feature.layers)[0]) + 1,
                "geometry": footprint,
            }
        )
        centroids = polys.geometry.iloc[feature.cells]
        for order, node in enumerate(feature.nodes):
            point = centroids.iloc[order]
            node_rows.append(
                {
                    "lake_id": feature.label,
                    "node": int(node),
                    "bed_elev": float(feature.elevations[order]),
                    "fskin": float(skin[order]) if order < skin.size else np.nan,
                    "bedleak": (float(skin[order]) / bed_thickness) if order < skin.size else np.nan,
                    "fleng": float(feature.lengths[order]),
                    "old_cell": int(feature.cells[order]),
                    "cell_active": bool(active[feature.cells[order]]),
                    "geometry": Point(point.centroid.x, point.centroid.y),
                }
            )

    _write_lakebed_raster(model, features, directory, crs, manifest)

    lakes = gpd.GeoDataFrame(lake_rows, geometry="geometry", crs=polys.crs).to_crs(crs)
    nodes = gpd.GeoDataFrame(node_rows, geometry="geometry", crs=polys.crs).to_crs(crs)
    lakes.to_file(path, layer="lakes", driver="GPKG")
    nodes.to_file(path, layer="lake_nodes", driver="GPKG")
    if forcing is not None:
        forcing.to_csv(directory / "lake_forcing.csv", index=False)
        manifest.tables["lake_forcing.csv"] = {
            "rows": int(len(forcing)),
            "what": "per lake per period: rainfall/evaporation in L/T over the lake area",
            "for": "mf.lak(rainfall=..., evaporation=...)",
        }

    manifest.layers["lakes.gpkg:lakes"] = {
        "features": int(len(lakes)),
        "what": "one polygon per CLN waterbody",
        "columns": [c for c in lakes.columns if c != "geometry"],
        "for": "mf.lak(lakes=..., lake_id_field='lake_id', starting_stage='strt', "
        "lake_bottom=Path('arrays/lakebed_bottom.tif'), bed_leakance='bedleak'). "
        "Use the RASTER for lake_bottom, not the per-lake column: the column is one "
        "number for a bed that varies, and mf.lak then rejects the cells it misses.",
        "geometry": "Polygon",
    }
    manifest.layers["lakes.gpkg:lake_nodes"] = {
        "features": int(len(nodes)),
        "what": "per-node bed elevation and bed K; the per-lake columns are medians of these",
        "columns": [c for c in nodes.columns if c != "geometry"],
        "geometry": "Point",
    }
    manifest.notes.append(
        f"bedleak = fskin / {bed_thickness:g} ft assumed bed thickness. CLN has no bed "
        "thickness (a conduit has a skin, not a bed), so this is a choice: rescale "
        "bedleak by your own thickness. fskin is written unmodified beside it."
    )
    _note_bed_at_top(model, features, manifest)
    _note_inactive_cells(features, active, manifest)
    spread = [row for row in lake_rows if row["fskin_max"] > 2.0 * max(row["fskin_min"], 1e-12)]
    if spread:
        manifest.notes.append(
            "fskin varies more than 2x within: "
            + ", ".join(f"{row['lake_id']} ({row['fskin_min']:.3g}-{row['fskin_max']:.3g})" for row in spread)
            + ". mf.lak takes one bed_leakance per lake, so the median is written; "
            "lake_nodes carries the variation."
        )


def _note_bed_at_top(model, features, manifest) -> None:
    """Say plainly when a CLN feature's bed IS the model top, because it usually is.

    Measured on Ten Trails every feature -- lakes as well as streams -- sits on the
    ground surface to within 0.00003 ft. Two consequences the caller must decide
    about rather than discover: the lakes carry **no bathymetry**, so a LAK built
    from these bottoms has no depth; and a bottom exactly at its cell's top is a
    knife edge that float noise tips either way, which is how ``mf.lak`` comes to
    reject a lake whose numbers look right.
    """

    coincident, above = [], 0
    for feature in features:
        difference = model.top[feature.cells] - feature.elevations
        if np.nanmax(np.abs(difference)) < 1e-3:
            coincident.append(feature.label)
            above += int((difference < -1e-6).sum())
    if not coincident:
        return
    manifest.notes.append(
        "the CLN bed IS the model top for "
        + ", ".join(coincident)
        + " (within 0.001 length units), so these features carry NO bathymetry -- a lake "
        f"built from them has zero depth. {above} node(s) round to just ABOVE their cell "
        "top, which mf.lak rejects as a bottom outside the cell; give lake_bottom your own "
        "depth below top rather than the exported elevation."
    )


def _note_inactive_cells(features, active: np.ndarray, manifest) -> None:
    """Report CLN nodes attached to cells that are inactive in every layer.

    MODFLOW-USG tolerates the connection; MODFLOW 6 refuses to place a reach or a
    lake connection there, so this is a rejection waiting to happen on whatever
    grid the caller builds next.
    """

    offenders = {
        feature.label: sorted(int(c) for c in np.unique(feature.cells) if not active[c])
        for feature in features
    }
    offenders = {label: cells for label, cells in offenders.items() if cells}
    if not offenders:
        return
    manifest.notes.append(
        "CLN node(s) attached to cells inactive in EVERY layer: "
        + "; ".join(f"{label} ({len(cells)}: {cells[:6]})" for label, cells in offenders.items())
        + ". USG allowed it, MODFLOW 6 will not -- see the cell_active column."
    )


def _write_lakebed_raster(model, features, directory: Path, crs: str, manifest) -> None:
    """Write the per-cell lakebed elevation as one raster covering every lake.

    ``mf.lak``'s ``lake_bottom=`` samples a raster at cell centres, which is the
    only route that survives a change of grid -- a per-lake scalar cannot, and on
    this network it is actively wrong: Keevie Lake's bed spans 505.90 to 516.77 ft,
    so one number for the lake puts most of its cells' bottoms outside the cells
    they sit in.
    """

    polys = model.grid.gdf_vorPolys
    bed = np.full(model.ncpl, np.nan)
    filled = 0
    for feature in features:
        bed[feature.cells] = feature.elevations
        # A dissolved footprint touches cells the CLN itself never had, and
        # LAKBuilder resolves lake cells from the polygon -- so a raster covering
        # only the CLN's own cells leaves the builder asking for a bottom that is
        # not there. Extend each lake to every cell its footprint reaches, taking
        # the nearest node's bed.
        known = np.unique(feature.cells)
        footprint = polys.iloc[known].union_all()
        # Through the spatial index: an elementwise `.intersects` against a large
        # dissolved polygon is O(ncpl) full comparisons and takes minutes here.
        touched = np.asarray(polys.sindex.query(footprint, predicate="intersects"))
        missing = touched[~np.isfinite(bed[touched])]
        if missing.size == 0:
            continue
        centres = polys.geometry.iloc[known].centroid
        known_xy = np.column_stack([centres.x.to_numpy(), centres.y.to_numpy()])
        gaps = polys.geometry.iloc[missing].centroid
        gap_xy = np.column_stack([gaps.x.to_numpy(), gaps.y.to_numpy()])
        nearest = np.argmin(
            np.hypot(
                gap_xy[:, 0, None] - known_xy[None, :, 0],
                gap_xy[:, 1, None] - known_xy[None, :, 1],
            ),
            axis=1,
        )
        bed[missing] = bed[known[nearest]]
        filled += int(missing.size)

    path = directory / "arrays" / "lakebed_bottom.tif"
    _write_raster(model, bed, path, crs, _resolution_for(model))
    manifest.rasters["arrays/lakebed_bottom.tif"] = {
        "files": 1,
        "what": "per-cell lakebed elevation, NaN away from the lakes"
        + (f"; {filled} fringe cell(s) filled from the nearest node" if filled else ""),
        "for": "mf.lak(lake_bottom='arrays/lakebed_bottom.tif')",
    }


def _resolution_for(model) -> float:
    """The default raster cell size for a grid: a quarter of the median cell width."""

    return float(np.sqrt(np.median(model.grid.gdf_vorPolys.geometry.area)) / 4.0)


def _lake_forcing(model: UsgModel, features):
    """Per-lake precipitation-minus-evaporation, as a rate over the lake area."""

    import pandas as pd

    wells = model.boundaries.get("WEL")
    if wells is None or not wells.cln_periods:
        return None
    polys = model.grid.gdf_vorPolys
    rows = []
    for feature in features:
        area = float(polys.iloc[np.unique(feature.cells)].union_all().area)
        nodes = np.asarray(feature.nodes, dtype=np.int64)
        for period, block in sorted(wells.cln_periods.items()):
            if block.size == 0:
                continue
            match = np.isin(block[:, 0].astype(np.int64), nodes)
            total = float(block[match, 1].sum())
            rate = total / area if area else np.nan
            rows.append(
                {
                    "lake_id": feature.label,
                    "period": int(period),
                    "net_q": total,
                    "rate": rate,
                    "rainfall": max(rate, 0.0),
                    "evaporation": max(-rate, 0.0),
                }
            )
    return pd.DataFrame(rows) if rows else None


# -- arrays: rasters, points, and the period index -------------------------


def _array_store(model: UsgModel) -> dict[str, Any]:
    """Every per-cell array the model carries, keyed by name.

    Static properties come back as one entry per layer; the transient forcings
    come back as a :func:`_cycle_index`, so a repeated period is stored once.
    """

    store: dict[str, Any] = {}
    for name in ("top",):
        store[name] = {"static": {name: np.asarray(getattr(model, name), float)}}
    for name in ("botm", "k", "k33", "ss", "sy", "strt"):
        values = np.atleast_2d(np.asarray(getattr(model, name), float))
        store[name] = {
            "static": {f"{name}_{i + 1}": values[i] for i in range(values.shape[0])}
        }
    if model.rch is not None and model.rch.rech:
        store["recharge"] = {"cyclic": _cycle_index(model.rch.rech)}
    if model.ets is not None:
        for field_name, source in (("et_rate", model.ets.rate), ("et_surface", model.ets.surf),
                                   ("et_depth", model.ets.depth)):
            if isinstance(source, dict) and source:
                store[field_name] = {"cyclic": _cycle_index(source)}
    return store


def _export_arrays(
    model: UsgModel, directory: Path, crs: str, manifest: ExportManifest, resolution: float | None
) -> None:
    """Write every per-cell array as a raster, and again as exact points."""

    import geopandas as gpd

    polys = model.grid.gdf_vorPolys
    if resolution is None:
        resolution = _resolution_for(model)
    store = _array_store(model)
    columns: dict[str, np.ndarray] = {}

    for name, entry in store.items():
        if "static" in entry:
            for label, values in entry["static"].items():
                columns[label] = values
                _write_raster(model, values, directory / "arrays" / f"{label}.tif", crs, resolution)
                manifest.rasters[f"arrays/{label}.tif"] = {
                    "files": 1,
                    "what": f"{name} (static)",
                }
        else:
            cyclic = entry["cyclic"]
            for position, values in enumerate(cyclic["distinct"], start=1):
                label = f"{name}_{position:02d}"
                columns[label] = values
                _write_raster(model, values, directory / "arrays" / f"{label}.tif", crs, resolution)
            manifest.rasters[f"arrays/{name}_NN.tif"] = {
                "files": len(cyclic["distinct"]),
                "what": f"{name}: {len(cyclic['periods'])} periods -> "
                f"{len(cyclic['distinct'])} distinct arrays"
                + (f", a {cyclic['cycle']}-period cycle" if cyclic["cycle"] else ""),
                "index": "periods.csv",
            }

    points = gpd.GeoDataFrame(
        {"cell": np.arange(model.ncpl), **{k: v for k, v in columns.items() if v.size == model.ncpl},
         "geometry": polys.geometry.centroid.to_numpy()},
        crs=polys.crs,
    ).to_crs(crs)
    points.to_file(directory / "cell_values.gpkg", layer="cell_values", driver="GPKG")
    manifest.layers["cell_values.gpkg:cell_values"] = {
        "features": int(len(points)),
        "what": "every array again at the old cell centres -- exact, unlike the rasters",
        "columns": [c for c in points.columns if c != "geometry"],
        "geometry": "Point",
    }
    worst = _raster_fidelity(directory / "arrays" / "top.tif", points, "top")
    manifest.notes.append(
        f"rasters are {resolution:g} units per pixel; the old mesh spans "
        f"{polys.geometry.area.min():.1f}-{polys.geometry.area.max():,.0f} square units, so cells "
        f"smaller than about {resolution ** 2:,.0f} square units are absorbed by their neighbours. "
        + (
            f"Measured on top.tif that costs at most {worst['max']:.2f} and typically "
            f"{worst['rms']:.3f} units of elevation. "
            if worst
            else ""
        )
        + "cell_values.gpkg carries every array at the cell centres and is exact -- use it "
        "where the difference matters."
    )


def _raster_fidelity(path: Path, points, column: str) -> dict[str, float] | None:
    """What a raster lost against the exact per-cell values, in the values' units.

    Sampling the raster back at the cell centres is the honest check: it is the
    same operation a consumer performs, so the number is what they would suffer
    rather than a bound derived from cell sizes.
    """

    if column not in points.columns or not path.exists():
        return None
    try:
        import rasterio

        with rasterio.open(path) as handle:
            sampled = np.array(
                [v[0] for v in handle.sample(zip(points.geometry.x, points.geometry.y, strict=False))],
                dtype=float,
            )
            sampled[sampled == handle.nodata] = np.nan
    except (OSError, ValueError, ImportError) as error:
        # A fidelity figure is a nicety; failing to read back the raster we just
        # wrote must not cost the caller the export.
        logger.debug("could not measure raster fidelity for %s: %s", path.name, error)
        return None
    difference = sampled - points[column].to_numpy(dtype=float)
    finite = np.isfinite(difference)
    if not finite.any():
        return None
    return {
        "max": float(np.nanmax(np.abs(difference))),
        "rms": float(np.sqrt(np.nanmean(difference[finite] ** 2))),
    }


def _write_raster(model: UsgModel, values: np.ndarray, path: Path, crs: str, resolution: float) -> None:
    """Burn one per-cell array onto a regular raster."""

    import geopandas as gpd
    import rasterio
    from rasterio.features import rasterize
    from rasterio.transform import from_origin

    polys = model.grid.gdf_vorPolys
    frame = gpd.GeoDataFrame(geometry=polys.geometry, crs=polys.crs).to_crs(crs)
    if values.size != len(frame):
        logger.debug("skipping raster %s: %d values for %d cells", path.name, values.size, len(frame))
        return
    minx, miny, maxx, maxy = frame.total_bounds
    width = max(int(np.ceil((maxx - minx) / resolution)), 1)
    height = max(int(np.ceil((maxy - miny) / resolution)), 1)
    transform = from_origin(minx, maxy, resolution, resolution)
    nodata = -9999.0
    burned = rasterize(
        ((geom, float(v)) for geom, v in zip(frame.geometry, values, strict=False) if np.isfinite(v)),
        out_shape=(height, width),
        transform=transform,
        fill=nodata,
        dtype="float32",
    )
    with rasterio.open(
        path, "w", driver="GTiff", height=height, width=width, count=1,
        dtype="float32", crs=crs, transform=transform, nodata=nodata, compress="deflate",
    ) as handle:
        handle.write(burned, 1)


def _export_tables(model: UsgModel, directory: Path, manifest: ExportManifest) -> None:
    """Write the period table: timing, and which array each package uses when."""

    import pandas as pd

    disu = model.disu
    nper = model.nper
    rows = {
        "period": np.arange(nper),
        "perlen": np.asarray(disu.perlen, float)[:nper],
        "nstp": np.asarray(disu.nstp, int)[:nper],
        "tsmult": np.asarray(disu.tsmult, float)[:nper],
        "steady": np.asarray(disu.steady)[:nper],
    }
    rows["day_start"] = np.r_[0.0, np.cumsum(rows["perlen"])[:-1]]
    store = _array_store(model)
    for name, entry in store.items():
        if "cyclic" not in entry:
            continue
        cyclic = entry["cyclic"]
        lookup = dict(zip(cyclic["periods"], cyclic["index"], strict=False))
        rows[f"{name}_array"] = [
            (lookup[p] + 1) if p in lookup else -1 for p in range(nper)
        ]
    table = pd.DataFrame(rows)
    table.to_csv(directory / "periods.csv", index=False)
    manifest.tables["periods.csv"] = {
        "rows": int(len(table)),
        "what": "period timing, and the 1-based array file each package uses",
        "columns": list(table.columns),
    }


def _readme(model: UsgModel, manifest: ExportManifest) -> str:
    """The prose that ships with an export: what each file is, and how to use it."""

    lines = [
        "# MODFLOW-USG export",
        "",
        f"Written from `{model.name_file.path.name}` in **{manifest.crs}**, "
        f"{manifest.nper} stress periods.",
        "",
        "Everything here is grid-independent: no file refers to the old mesh except",
        "`boundaries.gpkg:source_cells`, which is kept so every value can be traced back.",
        "",
        "## Files",
        "",
        manifest.describe(),
        "",
        "## Building packages from this",
        "",
        "```python",
        "import myflopy as mf, geopandas as gpd",
        "",
        "# Lakes -- the columns are already mf.lak's parameter names",
        'lakes = gpd.read_file("lakes.gpkg", layer="lakes")',
        "lak = mf.lak(context=ctx, nper=nper, lakes=lakes, lake_id_field='lake_id',",
        "             starting_stage='strt', bed_leakance='bedleak',",
        '             lake_bottom=Path("arrays/lakebed_bottom.tif"))   # the RASTER, not the column',
        "#   the per-lake `lake_bottom` column is one number for a bed that varies;",
        "#   the raster is per cell and covers the whole footprint, fringe included.",
        "",
        "# Streams -- likewise for mf.sfr",
        'streams = gpd.read_file("streams.gpkg", layer="streams")',
        "sfr = mf.sfr(context=ctx, nper=nper, streams=streams, stream_id='stream_id',",
        "             width='rwid', streambed_k='rhk', gradient='rgrd',",
        "             streambed_thickness='rbth', roughness='man')",
        "",
        "# A line BC on a new grid: conductance comes from the per-foot twin, so the",
        "# feature total survives a change of discretization.",
        'drn = gpd.read_file("boundaries.gpkg", layer="drn_segments")',
        "#   cond_new = drn.cond_per_ft * (length of the segment inside each new cell)",
        "```",
        "",
        "## Interpolating a stream profile onto finer reaches",
        "",
        "```python",
        'nodes = gpd.read_file("streams.gpkg", layer="stream_nodes")',
        'one = nodes[nodes.stream_id == "RockCreek"].sort_values("station")',
        "bed = np.interp(new_station, one.station, one.bed_elev)",
        "```",
        "",
        "`station` is measured along your own centerline when one was supplied to",
        "`export_gis(stream_lines=...)`, and along the CLN's otherwise; `station_cln`",
        "is always the CLN's.",
        "",
        "## What did not come across",
        "",
    ]
    lines += [f"- {note}" for note in manifest.notes]
    return "\n".join(lines) + "\n"
