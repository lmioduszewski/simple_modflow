"""Zone definitions for ``cal.parameterize(..., style="zone")``.

A *zone* is a per-cell integer label; ``style="zone"`` gives one adjustable
multiplier per distinct label. This module is the front door for saying where
those labels come from -- a raster, a polygon layer, or an array you already
have -- and, just as importantly, for normalizing them into the shape pyEMU
wants.

**Why normalization is not optional.** pyEMU's zone contract is asymmetric and
its failure modes are opaque. Measured on this stack (2026-07-30):

============ ================================= ==============================
family       accepted                          rejected
============ ================================= ==============================
array        ``(ncpl,)``  ``(ncpl,1)``          ``(nlay,ncpl)`` ->
             ``(1,ncpl)``                       ``write_array_tpl() error``
list         ``(nlay,ncpl)``, ``nlay >= 2``     per-cell shapes ->
                                                ``IndexError: index 1 is out
                                                of bounds for axis 1``
============ ================================= ==============================

So the same zone definition cannot be handed to both families, and passing the
wrong one raises from inside pyEMU naming nothing the caller wrote. Worse, on a
**single-layer** DISV model a list target cannot be zoned by pyEMU at all: the
correct ``(1, ncpl)`` shape hits ``checker2`` (``pst_from.py:2121-2134``), which
assumes a ``(1, n)`` array on a vertex grid is an idomain array and reshapes it
to ``(n, 1)``. :func:`resolve_zone_array` pads to ``(2, ncpl)`` to sidestep that;
the padding row addresses no real cell, so nothing indexes it.

That table is a vertex (DISV) grid's. On a structured (DIS) grid pyEMU reads an
array file as ``(nrow, ncol)`` and a list record's zone at its ``(layer, row, col)``,
so :func:`resolve_zone_array` returns ``(nrow, ncol)`` and ``(nlay, nrow, ncol)``
there, with no padding: the ``checker2`` reshape is a vertex-grid branch (ledger
420).

Zone ids are integers. ``0`` and negatives mean "not parameterized" for array
targets (pyEMU skips them); for list targets pyEMU makes a real parameter for
zone 0, so :func:`resolve_zone_array` warns rather than letting the two families
disagree silently.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from myflopy.modflow.mf6.grid.geometry import centroid_zones, model_cell_polygons, where_they_are

__all__ = ["ZoneSpec", "resolve_zone_array"]


@dataclass(frozen=True)
class ZoneSpec:
    """Per-cell integer zone labels, plus where they came from.

    Build one with :meth:`from_raster`, :meth:`from_polygons` or
    :meth:`from_array` and pass it as ``zones=``::

        cal.parameterize("k", style="zone", zones=ZoneSpec.from_raster("hsu.tif"))

    A bare numpy array still works as ``zones=`` -- this class adds the sources
    that need the grid to resolve, and a place to record provenance.

    Parameters
    ----------
    values : array_like of int
        One integer zone label per cell. Anything ``numpy.asarray`` accepts;
        it is stored as an ``int`` array, so a float label is truncated toward
        zero (``1.7`` becomes ``1``). An empty array, or one holding NaN or
        inf, raises ``ValueError``. A shape of ``(ncpl,)``, ``(1, ncpl)`` or
        ``(ncpl, 1)`` is one label per cell; ``(nlay, ncpl)`` is accepted only
        when every layer carries the same labels. The shape is checked against
        the grid when the parameter is built (``resolve_zone_array``), not
        here. For array targets a label ``<= 0`` is left unparameterized; for
        list targets zone ``0`` becomes a real parameter and a warning says so.
    source : str, optional
        Human-readable description of where the labels came from, quoted in
        error and warning messages. Default ``"array"``; the ``from_*``
        constructors fill it in.
    """

    #: One integer label per cell, length ``ncpl`` (2-D input is flattened per
    #: layer; see :func:`resolve_zone_array`).
    values: np.ndarray
    #: Human-readable description of where the labels came from, for errors.
    source: str = "array"

    def __post_init__(self):
        values = np.asarray(self.values)
        if values.size == 0:
            raise ValueError(f"Zone definition from {self.source} is empty.")
        if not np.all(np.isfinite(values.astype(float))):
            raise ValueError(
                f"Zone definition from {self.source} contains NaN/inf. Zones are "
                "integer labels; give unmapped cells an explicit id (e.g. 0)."
            )
        object.__setattr__(self, "values", values.astype(int))

    @property
    def ids(self) -> list[int]:
        """The distinct zone labels, sorted."""

        return sorted(int(value) for value in np.unique(self.values))

    @classmethod
    def from_array(cls, values, *, source: str = "array") -> ZoneSpec:
        """Zones from an array you already have (per-cell or per-layer).

        Parameters
        ----------
        values : array_like of int
            One integer label per cell, or one row per layer carrying the same
            labels; see the class ``values`` entry for the shapes, the int cast
            and what is refused.
        source : str, optional
            Description quoted in error and warning messages. Default
            ``"array"``.
        """

        return cls(values=np.asarray(values), source=source)

    @classmethod
    def from_raster(cls, path, model, *, nodata_zone: int = 0) -> ZoneSpec:
        """Zones from a CATEGORICAL raster, by majority vote per cell.

        Each cell takes the label held by the most raster pixels inside its
        polygon. Deliberately not the area-weighted MEAN that
        :func:`~myflopy.modflow.mf6.grid.surfaces.get_raster_vals_at_centroids`
        uses by default: averaging category ids is meaningless -- zones 1 and 3
        would average to 2, a zone the cell does not touch and may not exist.

        A cell too small to hold a pixel CENTRE takes the label under its
        centroid, as the surface sampler does; a cell whose centroid is off the
        raster, or on its no-data, gets ``nodata_zone``.

        Parameters
        ----------
        path : str or pathlib.Path
            Categorical raster opened with ``rasterio``; band 1 is read and its
            values cast to integers. Pixels equal to the file's declared nodata
            value never vote, and neither does a NaN or infinite pixel, whatever
            the file declares (a float raster's nodata is often NaN); a centroid
            on one counts as on nodata.
        model : SimulationBase
            The myflopy model the calibration is built on (``cal.model``). Its
            ``model.vor.gdf_vorPolys`` cell polygons are what is sampled,
            reprojected to the raster's CRS when both declare one; when either
            lacks a CRS a ``UserWarning`` is issued and the polygons are used
            as they are. The result has one label per cell of that grid.
        nodata_zone : int, optional
            Label for a cell that holds no valid pixel centre AND whose centroid
            is off the raster or on its nodata (a cell holding no pixel centre
            otherwise takes the label under its centroid). Default ``0``, which
            array targets leave unparameterized.
        """

        return cls(
            values=_majority_sample(Path(path), model, nodata_zone=nodata_zone),
            source=f"raster {Path(path).name}",
        )

    @classmethod
    def from_polygons(cls, polygons, model, *, column: str | None = None,
                      nodata_zone: int = 0) -> ZoneSpec:
        """Zones from a polygon layer: each cell takes the polygon it falls in.

        ``column`` names the integer label field; without it, polygons are
        numbered ``1..n`` in row order. Cells in no polygon get ``nodata_zone``.

        A cell falls in the zone holding its CENTROID, and a zone is CLOSED -- a
        centroid exactly on a zone's edge is in it (:func:`centroid_zones`). It
        used to be ``within``'s strict inside, so a boundary drawn along a row of
        centroids left that row in no zone: pyEMU gives a zone ``<= 0`` no
        multiplier, and the row went uncalibrated in silence (ledger 243). A
        cell in several zones takes the FIRST ROW, which is why overlapping zone
        polygons are a modelling mistake rather than something to resolve here.

        A row whose LABEL is missing -- a NULL attribute in the file -- is not a
        zone: its cells fall to whichever other row covers them, or to
        ``nodata_zone``. Casting that label would make it
        ``-9223372036854775808``, a number pyEMU calibrates like any other.
        A layer reaching NO cell warns, with both extents and both CRSs, since
        array targets skip ``nodata_zone`` entirely and would be parameterizing
        nothing.

        Parameters
        ----------
        polygons : geopandas.GeoDataFrame or str or pathlib.Path
            The zone polygons. A GeoDataFrame is used as given; anything else is
            opened with ``geopandas.read_file`` (its first layer). Reprojected
            to the grid's CRS when both declare one and they differ; with
            either CRS missing the coordinates are compared as they are.
        model : SimulationBase
            The myflopy model the calibration is built on (``cal.model``). Its
            ``model.vor.gdf_vorPolys`` cell polygons supply the centroids, and
            the result has one label per cell of that grid.
        column : str or None, optional
            Column holding each polygon's integer label; its values are cast
            with ``astype(int)``, so a float label is truncated. A row whose
            label is missing (NULL) is not a zone. Default ``None`` numbers the
            polygons ``1..n`` in row order. A name absent from the layer raises
            ``KeyError`` listing the columns it has.
        nodata_zone : int, optional
            Label for a cell whose centroid falls in no labelled polygon.
            Default ``0``, which array targets leave unparameterized.
        """

        import geopandas as gpd

        frame = polygons if isinstance(polygons, gpd.GeoDataFrame) else gpd.read_file(polygons)
        source = f"polygons ({column or 'row order'})"
        if column is not None and column not in frame.columns:
            raise KeyError(
                f"Zone column {column!r} is not in the polygon layer; it has "
                f"{sorted(frame.columns)}."
            )
        if column is not None:
            given = frame[column]
            labels, unlabelled = given.to_numpy(), given.isna().to_numpy()
        else:
            labels = np.arange(1, len(frame) + 1)
            unlabelled = np.zeros(len(frame), dtype=bool)

        cells = model_cell_polygons(model)
        source_crs = frame.crs
        if frame.crs is not None and cells.crs is not None and frame.crs != cells.crs:
            # Reprojection preserves row order, so `labels` still lines up.
            frame = frame.to_crs(cells.crs)

        values = np.full(len(cells), int(nodata_zone), dtype=int)
        # A row whose LABEL is missing is not a zone. It keeps its POSITION, so
        # `labels` still lines up, but is handed no geometry, which the shared rule
        # never matches: casting the missing label instead gives every one of its
        # cells -9223372036854775808, a label pyEMU calibrates (ledger 243).
        geometries = [
            None if missing else geometry
            for geometry, missing in zip(frame.geometry, unlabelled, strict=True)
        ]
        rows = centroid_zones(cells, geometries)
        matched = rows >= 0
        values[matched] = np.asarray(labels)[rows[matched]].astype(int)
        if not matched.any():
            warnings.warn(
                f"No cell centroid falls in any zone of {source}, so every cell is "
                f"zone {int(nodata_zone)} and array targets would be parameterizing "
                f"nothing. "
                + where_they_are(frame, cells, noun="zone polygons",
                                 source_crs=source_crs, grid_crs=cells.crs),
                UserWarning,
                stacklevel=2,
            )
        return cls(values=values, source=source)


def _majority_sample(raster_path: Path, model, *, nodata_zone: int) -> np.ndarray:
    """Per-cell majority label of a categorical raster.

    Mirrors ``grid/surfaces.py``'s ``_area_weighted_sample`` -- burn cell ids
    onto the raster grid, then reduce per cell -- but reduces by MODE instead of
    mean, which is the only meaningful reduction for categories. ``rasterize``
    burns a pixel only where its CENTRE falls in a cell, so a cell smaller than
    the raster's resolution holds none; it takes the label under its centroid,
    as the surface sampler's own fallback does. It used to take ``nodata_zone``,
    which froze 24% of a refined grid at zone 0.
    """

    import rasterio
    from rasterio.errors import WindowError
    from rasterio.features import rasterize
    from rasterio.windows import Window, from_bounds

    cells = model_cell_polygons(model)
    ncpl = len(cells)
    values = np.full(ncpl, int(nodata_zone), dtype=int)
    covered = np.zeros(ncpl, dtype=bool)

    with rasterio.open(raster_path) as src:
        polygons = cells.to_crs(src.crs) if (cells.crs and src.crs) else cells
        if not (cells.crs and src.crs):
            warnings.warn(
                f"CRS missing (grid={cells.crs}, raster={src.crs}); sampling "
                f"'{raster_path}' without reprojection.",
                stacklevel=3,
            )
        minx, miny, maxx, maxy = polygons.total_bounds
        window = None
        if np.all(np.isfinite([minx, miny, maxx, maxy])) and maxx > minx and maxy > miny:
            try:
                window = from_bounds(minx, miny, maxx, maxy, transform=src.transform)
                window = window.round_offsets().round_lengths()
                window = window.intersection(Window(0, 0, src.width, src.height))
            except WindowError:
                window = None
        if window is not None and window.width >= 1 and window.height >= 1:
            band = src.read(1, window=window)
            labels = rasterize(
                (
                    (geom, index + 1)
                    for index, geom in enumerate(polygons.geometry.values)
                    if geom is not None and not geom.is_empty
                ),
                out_shape=band.shape,
                transform=src.window_transform(window),
                fill=0,
                dtype="int32",
            )
            flat_cells = labels.ravel()
            flat_band = band.ravel()
            valid = (flat_cells > 0) & _has_a_label(flat_band, src.nodata)
            if valid.any():
                winners, modes = _cell_modes(flat_cells[valid] - 1,
                                             flat_band[valid].astype(np.int64))
                covered[winners] = True
                values[winners] = modes
        missing = np.flatnonzero(~covered)
        if missing.size:
            centroids = cells.geometry.iloc[missing].centroid
            if cells.crs and src.crs:
                centroids = centroids.to_crs(src.crs)
            values[missing] = _centroid_labels(src, centroids, nodata_zone=nodata_zone)
    return values


def _has_a_label(band: np.ndarray, nodata) -> np.ndarray:
    """Which pixels hold a label: a finite value that is not the file's nodata.

    A NaN never equals itself, so ``band != nodata`` kept every pixel of a float
    raster whose nodata is NaN, and the cast to integer made each such pixel the label
    ``-9223372036854775808``: ``from_raster`` died with ``OverflowError`` where the
    pixels outvoted a cell, and a cell whose centroid sat on one took that label in
    silence. A non-finite value is no label whatever the file declares.
    """

    band = np.asarray(band)
    keep = np.isfinite(band) if np.issubdtype(band.dtype, np.floating) else np.ones(band.shape, bool)
    if nodata is not None and not (isinstance(nodata, float) and np.isnan(nodata)):
        keep &= band != nodata
    return keep


def _cell_modes(cell_ids: np.ndarray, zone_ids: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Each cell's most frequent zone label (the lowest label on a tie).

    Counts the distinct ``(cell, label)`` PAIRS, so memory follows the pixels read.
    It counted into a ``ncells x (label span)`` table, which an undeclared nodata of
    65,535 beside labels 1 and 2 made 463 MB at 441 cells and would make ~9.5 GB at
    18,107.

    Parameters
    ----------
    cell_ids : numpy.ndarray of int
        The zero-based cell each pixel lies in.
    zone_ids : numpy.ndarray of int
        The pixel's label, the same length.

    Returns
    -------
    tuple of (numpy.ndarray, numpy.ndarray)
        The cells holding any pixel, and each one's label.
    """

    labels, label_index = np.unique(zone_ids, return_inverse=True)
    pairs, counts = np.unique(
        cell_ids.astype(np.int64) * len(labels) + label_index.reshape(-1), return_counts=True
    )
    cells, label_of = pairs // len(labels), pairs % len(labels)
    # Per cell, the highest count first and, among equal counts, the lowest label.
    order = np.lexsort((label_of, -counts, cells))
    first = np.ones(order.size, dtype=bool)
    first[1:] = cells[order][1:] != cells[order][:-1]
    chosen = order[first]
    return cells[chosen], labels[label_of[chosen]]


def _centroid_labels(src, points, *, nodata_zone: int) -> np.ndarray:
    """The label under each point; ``nodata_zone`` off the raster or on its no-data.

    Bounds are checked here rather than left to ``src.sample``, which hands back
    its fill value -- 0 when the file declares no nodata -- for a point off the
    raster: a real-looking label for a cell the raster never covered.
    """

    xs, ys = points.x.to_numpy(), points.y.to_numpy()
    left, bottom, right, top = src.bounds
    inside = np.flatnonzero((xs >= left) & (xs < right) & (ys > bottom) & (ys <= top))
    labels = np.full(len(xs), int(nodata_zone), dtype=int)
    if inside.size:
        sampled = np.array(
            [value[0] for value in src.sample(zip(xs[inside], ys[inside], strict=True))]
        )
        keep = _has_a_label(sampled, src.nodata)
        labels[inside[keep]] = sampled[keep].astype(int)
    return labels


def resolve_zone_array(zones: Any, *, family: str, model) -> np.ndarray:
    """Normalize a zone definition into the shape pyEMU wants for ``family``.

    The result is per-cell ``(ncpl,)`` for array targets and ``(nlay, ncpl)`` for
    list targets on a vertex grid. See the module docstring for why the two differ
    and why a single-layer list target is padded. On a structured (DIS) grid pyEMU
    reads an array file as ``(nrow, ncol)`` and a list record's zone at its
    ``(layer, row, col)``, so the result is ``(nrow, ncol)`` and
    ``(nlay, nrow, ncol)`` there (ledger 420). On DISU an array file is one node
    per line and a list record's cellid is its ``node``, so the result is
    ``(nodes, 1)`` and ``(nodes,)`` (ledger 452).

    Parameters
    ----------
    zones : ZoneSpec or array_like
        The zone labels: a :class:`ZoneSpec`, or any integer array-like
        (:meth:`ZoneSpec.from_array`) of one label per cell of a layer -- flat,
        ``(nrow, ncol)`` on DIS, or one copy per layer, which must agree.
    family : str
        The target's family: ``"array"`` or ``"list"``; anything else is read as
        ``"list"``.
    model : SimulationBase
        The calibrated model; its grid view's ``vor.ncpl`` (else the FloPy grid's
        ``ncpl``) and its FloPy ``gwf.modelgrid`` (``nlay``, ``grid_type``, and
        ``nrow``/``ncol`` on DIS) are read.

    Returns
    -------
    numpy.ndarray
        Integer labels in the shape above.

    Raises
    ------
    ValueError
        When the labels are not one per cell, or their per-layer copies differ.
    """

    spec = zones if isinstance(zones, ZoneSpec) else ZoneSpec.from_array(zones)
    values = spec.values
    modelgrid = model.gwf.modelgrid
    vor = getattr(model, "vor", None)
    ncpl = int(vor.ncpl) if vor is not None else int(modelgrid.ncpl)
    nlay = int(modelgrid.nlay)

    flat = values.reshape(-1)
    if values.ndim >= 2 and flat.size != ncpl and flat.size % ncpl == 0:
        # Per-layer input: zones are a property of the map, so require the
        # layers to agree rather than silently taking the first.
        layers = flat.reshape(-1, ncpl)
        if not np.all(layers == layers[0]):
            raise ValueError(
                f"Zone definition from {spec.source} differs between layers. "
                "Zones are per-cell labels applied to every layer the target "
                "covers; use `layers=` to restrict which layers are "
                "parameterized instead."
            )
        flat = layers[0]
    elif flat.size != ncpl:
        raise ValueError(
            f"Zone definition from {spec.source} has {flat.size} values, but the "
            f"grid has {ncpl} cells. Zones are one integer label per cell."
        )

    per_cell = np.asarray(flat, dtype=int).reshape(ncpl)
    structured = getattr(modelgrid, "grid_type", None) == "structured"
    if getattr(modelgrid, "grid_type", None) == "unstructured":
        # DISU: pyEMU reads an array file one node per line as `(nodes, 1)` and keys
        # a list record's zone by its one `node` column, so the zones are those two
        # shapes; given a flat array for a file it raised `passed shape (24, 1) !=
        # zone_array.shape (24,)` (ledger 452).
        return per_cell.reshape(ncpl, 1) if family == "array" else per_cell
    if family == "array":
        return per_cell.reshape(int(modelgrid.nrow), int(modelgrid.ncol)) if structured else per_cell

    # List targets are indexed (layer, cell), so the array must be 2-D -- and at
    # least 2 rows deep even on a one-layer model, because pyEMU reshapes any
    # (1, n) array on a vertex grid assuming it is an idomain array
    # (pst_from.py checker2). The padding row addresses no real cell, so nothing
    # ever indexes it.
    if structured:
        stacked = np.repeat(
            per_cell.reshape(1, int(modelgrid.nrow), int(modelgrid.ncol)), nlay, axis=0
        )
    else:
        stacked = np.repeat(per_cell[None, :], max(nlay, 2), axis=0)
    if 0 in set(per_cell.tolist()):
        warnings.warn(
            f"Zone id 0 in the definition from {spec.source} becomes a real, "
            "adjustable parameter for list targets, while array targets skip it. "
            "Use a positive id for every zone you want calibrated, and keep 0 "
            "for 'leave alone'.",
            UserWarning,
            stacklevel=3,
        )
    return stacked
