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
exactly and re-resolve onto any mesh. **Neighbour means a shared face**, the
adjacency MODFLOW 6 itself connects: on a structured or quadtree mesh a cell pair
meeting at one corner is no more joined than two cells a mile apart, and a chain
drawn through a corner claims length the feature does not occupy. A cell left alone
on its chain -- nothing shares a face with it, or the open walk could not reach it
through a branch -- is counted and reported rather than left to be discovered.

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
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

from myflopy._logging import get_logger
from myflopy._warnings import stacklevel_outside_package
from myflopy.modflow.mf6.package_registry import _PACKAGE_EXPLORER_SPECS, agg_for_kind

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
        Pixel size of EVERY raster written (the per-cell arrays and the lakebed),
        in the units of the output ``crs`` -- the raster's own units, since the
        cells are reprojected before the raster is laid out. Defaults to a
        quarter of the median cell width, measured on the cells AS REPROJECTED
        into ``crs``, which resolves a typical cell and absorbs the smallest
        ones -- see the module docstring. A value that is not a finite number
        above zero is refused (ledger 369).
    bed_thickness
        Thickness assumed when deriving ``bedleak`` from ``FSKIN``.
    overwrite
        Replace existing files rather than refusing.

    Returns
    -------
    ExportManifest

    Raises
    ------
    ValueError
        No grid; no CRS given or on the grid; a ``resolution`` that is not finite
        and above zero.
    TypeError
        A ``resolution`` that is not a number.
    FileExistsError
        ``directory`` is not empty and ``overwrite`` is false.
    """

    if model.grid is None:
        raise ValueError(
            "export_gis() needs the grid; read the model with gsf= so it has geometry."
        )
    target_crs = crs or _grid_crs(model)
    # ONE pixel size for every raster, in the OUTPUT CRS's units (ledger 369),
    # resolved before anything is written so a refusal leaves no directory behind.
    pixel = _pixel_size(model, target_crs, resolution)

    directory = Path(directory)
    if directory.exists() and any(directory.iterdir()) and not overwrite:
        raise FileExistsError(f"{directory} is not empty; pass overwrite=True to replace it.")
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "arrays").mkdir(exist_ok=True)

    manifest = ExportManifest(directory=directory, crs=str(target_crs), nper=model.nper)

    _export_boundaries(model, directory, target_crs, manifest)
    _export_cln(
        model, directory, target_crs, manifest, stream_lines, clip_to, bed_thickness, pixel
    )
    _export_arrays(model, directory, target_crs, manifest, pixel)
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


def _adjacency(cells: list[int], grid) -> dict[int, set[int]]:
    """Which of ``cells`` share a FACE with which, as an adjacency map.

    Adjacency here is MODFLOW 6's: a shared **edge**, not a shared corner. Two
    cells meeting at a single vertex -- every diagonal step of a structured or
    quadtree mesh -- are not neighbours, because nothing flows across a corner
    and no one-cell-wide feature runs through one. ``grid.shared_face`` is the
    test -- the rule the HFB views already resolve a barrier's own face against,
    so "adjacent" means one thing across the package. It agrees pair-for-pair
    with the older ``grid.adjacent_cells_idx`` on every mesh measured (a 16-cell
    square mesh; a 441-cell irregular Voronoi, 2,494 pairs, 0 differences) and
    additionally handles a face interrupted into a MultiLineString and a face
    tolerance, which a bare ``isinstance(intersection, LineString)`` does not.

    The spatial index only proposes candidate pairs; the face test decides.

    Parameters
    ----------
    cells
        Zero-based plan cell indices, the subset to relate.
    grid
        The model's grid. Needs ``gdf_vorPolys`` *and* ``shared_face`` -- a
        stand-in exposing only the polygon frame is not enough here.
    """

    polys = grid.gdf_vorPolys
    subset = polys.iloc[cells]
    candidates = subset.sindex.query(subset.geometry, predicate="intersects")
    adjacent: dict[int, set[int]] = {c: set() for c in cells}
    for a, b in zip(*candidates, strict=False):
        if a == b:
            continue
        left, right = cells[a], cells[b]
        if grid.shared_face(left, right) is not None:
            adjacent[left].add(right)
    return adjacent


def _chains(cells: list[int], grid) -> list[list[int]]:
    """Order a set of face-sharing cells into open chains, one per component.

    Two cells meeting only at a corner are not connected, so each is its own
    chain. So is a branch arm the open walk below cannot extend through -- a
    chain has two ends and a branch point has three or more. Either way a cell
    alone on its chain gets the degenerate segment :func:`_segment_line` gives a
    lone cell, and :func:`_export_boundaries` counts and reports them.
    """

    adjacent = _adjacency(cells, grid)
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
            chains.extend(_chains(leftover, grid) if len(leftover) > 1 else [leftover])
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

    A source cell is every ``(layer, cell)`` ANY period gives a record, with the
    values of the first period that does, and a ``periods`` column naming the model
    periods a record on it is in force (see :func:`_bc_frames`). The model periods in
    which NO record of the package is in force -- ``ITMP = 0``, all records on CLN
    nodes, or before the first definition -- go into both manifest entries as
    ``switched_off_periods``, with a note, because a switch-off has no row of its own
    (ledger 430). Both entries count, as ``features_in_some_periods``, the
    features in force in only some of the periods the package is -- a cell added or
    dropped between periods -- with a note too, since a layer read with ``nper=``
    applies every feature in every period.

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

    A chain runs between cells that share a **face**. A cell alone on its chain --
    nothing else in the feature shares a face with it, or it is a branch arm the open
    walk could not extend through -- keeps the degenerate segment a lone cell has
    always had: ``length_ft`` 0, and a NaN ``*_per_ft`` twin where the package has
    one. The count goes into the manifest as ``cells_alone_on_a_chain``, with a note
    and a warning pointing at the point layer.
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
            if c not in ("geometry", "name", "layer", "old_cell", "length_ft",
                         "first_period", "periods")
        ]
        # A switch-off has no row to show it, so the manifest says it (ledger 430).
        switched_off = _switched_off(model, ftype)
        if switched_off:
            manifest.notes.append(
                f"{ftype}: no groundwater record is in force in {len(switched_off)} of "
                f"{model.nper} period(s) (zero-based: {_period_ranges(switched_off)}) -- "
                "ITMP = 0, every record on a CLN node, or before the first period that "
                "defines records. Each feature's 'periods' column leaves them out; a "
                "layer read with nper= applies its features in every period."
            )
        # A feature the package's records name in only SOME of its periods is the same
        # trap without a switch-off: read with nper= it acts in all of them, and only
        # its own 'periods' column said otherwise (Review of U6, ledger 430).
        everywhere = _period_ranges(
            [p for p in range(int(model.nper)) if p not in set(switched_off)]
        )
        partial = int((points["periods"] != everywhere).sum())
        # Several records on one cell in one period are ONE feature (ledger 467).
        combined = int(points.attrs.get("records_combined", 0))
        disagreeing = int(points.attrs.get("cells_disagreeing", 0))
        if combined:
            manifest.notes.append(
                f"{ftype}: {combined} record(s) shared a cell with another record of the "
                "same period and were combined into its feature -- conductance and rate "
                "summed, an elevation or head the mean"
                + (
                    f"; on {disagreeing} cell(s) the copies' elevations differ, so the "
                    "feature is not the boundary MF6 applied there"
                    if disagreeing else ""
                )
                + "."
            )
        if partial:
            manifest.notes.append(
                f"{ftype}: {partial} of {len(points)} feature(s) are in force in only some "
                f"of the periods the package is ({everywhere or 'none'}) -- their cell is "
                "added or dropped between periods. A layer read with nper= applies them in "
                "every period; each one's 'periods' column names where it acts."
            )
        points.to_file(path, layer=package, driver="GPKG")
        manifest.layers[f"boundaries.gpkg:{package}"] = {
            "features": int(len(points)),
            "what": f"{ftype} as one point per source cell -- exact conductance, discontinuous",
            "consumed_by": f"mf.{package}.gpkg(path, layer='{package}', context=ctx, nper=nper)",
            "columns": value_columns,
            "layer_base": 1,
            "period_base": 0,
            "switched_off_periods": switched_off,
            "features_in_some_periods": partial,
            "records_combined": combined,
            "geometry": "Point",
        }
        lines.to_file(path, layer=f"{package}_lines", driver="GPKG")
        # Count the CHAINS of one cell, not the zero lengths: the two agree on any
        # frame this writes, but the chain is what the sentence below is about, and
        # a cell alone on its chain is not necessarily one with no face neighbour --
        # a branch arm the open walk could not extend through lands here too.
        alone = int((lines.groupby(["layer", "chain"]).size() == 1).sum())
        per_ft = [c for c in lines.columns if c.endswith("_per_ft")]
        # CHD carries no extensive column, so it has no twin to go NaN.
        twins = f" and their {' and '.join(per_ft)} is NaN" if per_ft else ""
        if alone:
            logger.warning(
                "%s: %d of %d source cell(s) are alone on their chain -- no face neighbour "
                "in the feature, or a branch the open walk could not extend through -- so "
                "their segment has no length%s; take those cells from the '%s' point layer",
                package, alone, len(lines), twins, package,
            )
            manifest.notes.append(
                f"{ftype}: {alone} of {len(lines)} source cell(s) are alone on their chain -- "
                "either no other cell of the feature shares a face with them (a shared corner "
                "is not a MODFLOW connection), or they are a branch the open walk could not "
                f"extend through. Their line segment has zero length{twins}; build those cells "
                f"from the '{package}' point layer, which is exact for one cell."
            )
        manifest.layers[f"boundaries.gpkg:{package}_lines"] = {
            "features": int(len(lines)),
            "cells_alone_on_a_chain": alone,
            "switched_off_periods": switched_off,
            "features_in_some_periods": partial,
            "records_combined": combined,
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


def _value_kind(ftype: str, column: str) -> str:
    """What a USG list-BC value column IS -- ``"extensive"``, ``"elevation"`` -- per the registry.

    Read off the package registry's :class:`FieldSpec` for the MF6 package, by the USG
    name (``elev``, ``cond``, ``bhead``, ``stage``, ``rbot``, ``q``), else by its
    canonical one (CHD's ``shead`` and ``ehead`` are both its ``head``).
    """

    inputs = _PACKAGE_EXPLORER_SPECS[ftype.lower()].inputs
    spec = inputs.get(column) or inputs.get(_CANONICAL.get(column, column))
    if spec is None and ftype == "CHD":
        spec = inputs["head"]
    return spec.kind


def _combine_records(frame, ftype: str, columns: list[str]):
    """One row per ``(period, layer, cell)``: several records on one cell combined by KIND.

    MF6 applies every record a period gives one cell -- two drains, two wells -- and
    the export writes one feature per cell, so they are combined by what each value
    IS (:func:`~myflopy.modflow.mf6.package_registry.agg_for_kind`, ledger 297): a
    conductance or a rate ADDS, an elevation or a head is the MEAN. ``_bc_frames``
    read such a cell's values as a Series and stopped on ``float()`` (ledger 467).

    Returns
    -------
    tuple
        ``(frame, disagreeing)``: the combined frame (its ``attrs`` kept), and the
        ``(period, layer, cell, column)`` of every non-extensive value whose copies on
        one cell DIFFER, where the mean is not the boundary MF6 applied.
    """

    keys = ["period", "layer", "cell"]
    if not frame.duplicated(keys).any():
        return frame, []
    rules = {column: agg_for_kind(_value_kind(ftype, column)) for column in columns}
    grouped = frame.groupby(keys, sort=False)
    disagreeing = [
        (int(period), int(layer), int(cell), column)
        for column, rule in rules.items()
        if rule != "sum"
        for (period, layer, cell), spread in (grouped[column].max() - grouped[column].min()).items()
        if spread > 0
    ]
    others = {column: "first" for column in frame.columns if column not in keys and column not in rules}
    combined = grouped.agg({**rules, **others}).reset_index()[list(frame.columns)]
    combined = type(frame)(combined, geometry="geometry", crs=frame.crs) if "geometry" in frame else combined
    combined.attrs = dict(frame.attrs)
    return combined, disagreeing


def _bc_frames(model: UsgModel, ftype: str, crs: str):
    """One list BC as ``(points, lines)``, carrying the package registry's field names.

    Both frames hold one row per source cell and the same values; they differ only
    in geometry and in the ``*_per_ft`` twins, which are meaningless on a point.

    The rows are every ``(layer, cell)`` ANY period gives a record, not the first
    defined period's alone (ledger 430): a cell a later period adds was absent
    from the export. Each row's values are those of the first period that names it
    (``first_period``), and ``periods`` lists every model period in which a record
    on it is in force, reuse resolved -- so a period the package was switched off in
    (``ITMP = 0``, or every record on a CLN node) is missing from every row's list.
    When the first value column moves in time, one ``<column>_mNN`` column per model
    period of its cycle carries it, NaN where the cell has no record in force.

    Several records one period gives one cell are combined first, by kind
    (:func:`_combine_records`): MF6 applies each, and a feature is one per cell. Both
    frames' ``attrs`` say how many records were combined (``records_combined``) and on
    how many cells the copies' elevations differ (``cells_disagreeing``), which also
    warns.

    The chains come from face adjacency (see :func:`_adjacency`), so a corner-only
    contact starts a new chain rather than extending one through the corner.
    """

    import geopandas as gpd

    frame = model.boundary_frame(ftype)
    if frame.empty:
        return None, None
    polys = model.grid.gdf_vorPolys
    columns = [c for c in _BC_VALUES.get(ftype, ()) if c in frame.columns]
    rows_before = len(frame)
    frame, disagreeing = _combine_records(frame, ftype, columns)
    if disagreeing:
        cells = sorted({(layer + 1, cell) for _period, layer, cell, _column in disagreeing})
        message = (
            f"{ftype}: {len(cells)} cell(s) hold several records in one period whose "
            f"{', '.join(sorted({column for *_, column in disagreeing}))} values differ -- "
            f"(layer, cell) {cells[:5]}{' ...' if len(cells) > 5 else ''}. MF6 applies each "
            "record; the export writes ONE feature per cell, its conductance or rate the "
            "sum and its elevation the MEAN, which is not the same boundary there."
        )
        logger.warning(message)
        warnings.warn(message, UserWarning, stacklevel=stacklevel_outside_package())
    acting = model.boundaries[ftype].in_force(model.nper)

    package = ftype.lower()
    # CHD's second USG field is EHEAD, the end-of-period ramp; MF6 has one HEAD, so
    # only SHEAD becomes `head` and EHEAD rides along under its own name (see report()).
    canonical = {c: _CANONICAL.get(c, c) for c in columns}
    if ftype == "CHD":
        canonical["ehead"] = "ehead_usg_only"

    # Each (layer, cell) any period names, with the first defined period naming it.
    first_period = (
        frame.groupby(["layer", "cell"], sort=True)["period"].min().astype(int).to_dict()
    )
    defining = {
        int(period): set(zip(block["layer"].astype(int), block["cell"].astype(int), strict=True))
        for period, block in frame.groupby("period")
    }

    blocks = {
        (int(period), int(layer)): block.set_index("cell")
        for (period, layer), block in frame.groupby(["period", "layer"])
    }

    point_rows, line_rows = [], []
    for layer in sorted({layer for layer, _ in first_period}):
        cells = [cell for (lay, cell) in first_period if lay == layer]
        for chain_id, walk in enumerate(_chains(cells, model.grid)):
            points = [polys.geometry.iloc[c].centroid for c in walk]
            for index, cell in enumerate(walk):
                key = (int(layer), int(cell))
                source = first_period[key]
                values = blocks[(source, int(layer))]
                line = _segment_line(points, index)
                shared = {
                    "name": f"{package}_L{int(layer) + 1}_c{int(cell)}",
                    "layer": int(layer) + 1,
                    "old_cell": int(cell),
                    "length_ft": float(line.length),
                    "first_period": int(source),
                    "periods": _period_ranges(
                        [p for p, act in enumerate(acting)
                         if act is not None and key in defining.get(act, ())]
                    ),
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
    out_points = _attach_cycle(out_points, frame, columns[0], acting, name=canonical[columns[0]])
    out_lines = _attach_cycle(out_lines, frame, columns[0], acting, name=canonical[columns[0]])
    out_points, out_lines = out_points.to_crs(crs), out_lines.to_crs(crs)
    combined = {
        "records_combined": rows_before - len(frame),
        "cells_disagreeing": len({(layer, cell) for _p, layer, cell, _c in disagreeing}),
    }
    out_points.attrs.update(combined)
    out_lines.attrs.update(combined)
    return out_points, out_lines


def _period_ranges(periods: list[int]) -> str:
    """Zero-based periods as compact text: ``[0, 1, 2, 5]`` -> ``"0-2,5"``; none -> ``""``."""

    runs: list[str] = []
    start = previous = None
    for period in periods:
        if previous is not None and period == previous + 1:
            previous = period
            continue
        if start is not None:
            runs.append(f"{start}-{previous}" if previous != start else f"{start}")
        start = previous = period
    if start is not None:
        runs.append(f"{start}-{previous}" if previous != start else f"{start}")
    return ",".join(runs)


def _switched_off(model: UsgModel, ftype: str) -> list[int]:
    """Model periods in which no groundwater record of ``ftype`` is in force."""

    records = model.boundaries[ftype]
    return [
        period for period, act in enumerate(records.in_force(model.nper))
        if act is None or len(records.periods[act]) == 0
    ]


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


def _attach_cycle(out, frame, column: str, acting: list[int | None], *, name: str | None = None):
    """Add one column per cycle position for a BC whose values move in time.

    Each MODFLOW period ``p`` reads the records in force there (``acting[p]``, reuse
    resolved), so a column is one model period: ``<name>_m01`` is period 0, and a
    period the package was switched off in is NaN on every row -- it used to be
    skipped, along with every reused period, so ``_m02`` meant "the second period
    that DEFINED records" (ledger 430). When the periods repeat in a cycle only
    its first cycle is written; when every period agrees nothing is.

    Keyed on ``(layer, cell)``, not on the cell alone: a boundary is routinely
    applied to the same plan cell in several layers -- Ten Trails' CHD covers 23
    cells in each of layers 1-3 -- and keying on the cell would silently pick one
    of the three, or raise where a scalar was expected.
    """

    lookups: dict[int, dict[tuple[int, int], float]] = {}
    for period, block in frame.groupby("period"):
        lookups[int(period)] = {
            (int(layer) + 1, int(cell)): float(value)
            for layer, cell, value in zip(
                block["layer"], block["cell"], block[column], strict=False
            )
        }
    keys = list(zip(out["layer"].astype(int), out["old_cell"].astype(int), strict=False))
    store = {
        period: np.array(
            [lookups.get(act, {}).get(key, np.nan) if act is not None else np.nan
             for key in keys],
            float,
        )
        for period, act in enumerate(acting)
    }
    if len(store) < 2:
        return out
    cyclic = _cycle_index(store)
    if len(cyclic["distinct"]) < 2:
        return out
    cycle = cyclic["cycle"] or len(cyclic["periods"])
    for position, period in enumerate(cyclic["periods"][:cycle], start=1):
        out[f"{name or column}_m{position:02d}"] = store[period]
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
    resolution: float,
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
        _write_lakes(model, lakes, directory, crs, manifest, bed_thickness, resolution)


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
    resolution: float,
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

    _write_lakebed_raster(model, features, directory, crs, manifest, resolution)

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
    knife edge that float noise tips either way. ``mf.lak`` used to reject such a
    lake outright; through the raster it now caps the ones above at their top
    and warns, so the lake builds -- still with no depth.
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
        "top, which mf.lak caps at the top with a warning; give lake_bottom your own "
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


def _write_lakebed_raster(
    model, features, directory: Path, crs: str, manifest, resolution: float
) -> None:
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
    _write_raster(model, bed, path, crs, resolution)
    manifest.rasters["arrays/lakebed_bottom.tif"] = {
        "files": 1,
        "what": "per-cell lakebed elevation, NaN away from the lakes"
        + (f"; {filled} fringe cell(s) filled from the nearest node" if filled else ""),
        "for": "mf.lak(lake_bottom='arrays/lakebed_bottom.tif')",
    }


def _cells_in(model, crs):
    """The grid's cell polygons, reprojected into the output ``crs``."""

    polys = model.grid.gdf_vorPolys
    return polys.geometry.to_crs(crs)


def _resolution_for(model, crs) -> float:
    """The default raster cell size: a quarter of the median cell width, IN ``crs``.

    Measured on the cells as reprojected, because that is the frame the raster is
    laid out in: taken in the grid's units it was 3.28x too coarse for a feet grid
    written to a metre CRS (ledger 369). In a geographic CRS the width is in
    degrees, which is what that raster's pixels are measured in, so geopandas'
    warning about areas in degrees is not news here.
    """

    import warnings

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Geometry is in a geographic CRS")
        areas = _cells_in(model, crs).area
    return float(np.sqrt(np.median(areas)) / 4.0)


def _pixel_size(model, crs, resolution: float | None) -> float:
    """The one pixel size every raster of an export uses, in ``crs`` units."""

    if resolution is None:
        return _resolution_for(model, crs)
    if isinstance(resolution, bool) or not isinstance(resolution, (int, float, np.number)):
        raise TypeError(
            f"resolution={resolution!r}: give the pixel size as a number, in the units "
            f"of the output crs ({crs})."
        )
    if not np.isfinite(resolution) or resolution <= 0:
        raise ValueError(
            f"resolution={resolution!r}: the pixel size must be a finite number above "
            f"zero, in the units of the output crs ({crs})."
        )
    return float(resolution)


def _unit_name(crs) -> str:
    """The output CRS's length unit, for a sentence; ``"units"`` when it has none."""

    import pyproj

    try:
        axes = pyproj.CRS.from_user_input(crs).axis_info
    except pyproj.exceptions.CRSError as error:
        logger.debug("no unit name for %s: %s", crs, error)
        return "units"
    return axes[0].unit_name if axes else "units"


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
    model: UsgModel, directory: Path, crs: str, manifest: ExportManifest, resolution: float
) -> None:
    """Write every per-cell array as a raster, and again as exact points.

    ``resolution`` is the pixel size :func:`_pixel_size` resolved, in ``crs`` units.
    """

    import geopandas as gpd

    polys = model.grid.gdf_vorPolys
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
    # Every number in this sentence is in the OUTPUT crs's units, the pixel size's
    # own: the areas were the grid's before, beside a pixel size in the crs's.
    import warnings

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Geometry is in a geographic CRS")
        areas = _cells_in(model, crs).area
    unit = _unit_name(crs)
    manifest.notes.append(
        f"rasters are {resolution:g} {unit} per pixel; the old mesh spans "
        f"{areas.min():.1f}-{areas.max():,.0f} square {unit}, so cells "
        f"smaller than about {resolution ** 2:,.0f} square {unit} are absorbed by their neighbours. "
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
        'drn = gpd.read_file("boundaries.gpkg", layer="drn_lines")',
        "#   cond_new = drn.conductance_per_ft * (length of the segment in each new cell)",
        "#   A segment of length 0 is a cell alone on its chain -- nothing shares a face",
        "#   with it, or the walk could not reach it through a branch. Its per-foot twin",
        "#   is NaN, so take that cell from the point layer. The manifest counts them.",
        "#   A feature is every cell ANY period gives a record: `periods` lists the",
        "#   zero-based periods one is in force on it, `first_period` the period its",
        "#   values come from, and a period missing from every row is a switch-off.",
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
