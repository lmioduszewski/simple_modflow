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

__all__ = ["ZoneSpec", "resolve_zone_array"]


@dataclass(frozen=True)
class ZoneSpec:
    """Per-cell integer zone labels, plus where they came from.

    Build one with :meth:`from_raster`, :meth:`from_polygons` or
    :meth:`from_array` and pass it as ``zones=``::

        cal.parameterize("k", style="zone", zones=ZoneSpec.from_raster("hsu.tif"))

    A bare numpy array still works as ``zones=`` -- this class adds the sources
    that need the grid to resolve, and a place to record provenance.
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
        """Zones from an array you already have (per-cell or per-layer)."""

        return cls(values=np.asarray(values), source=source)

    @classmethod
    def from_raster(cls, path, model, *, nodata_zone: int = 0) -> ZoneSpec:
        """Zones from a CATEGORICAL raster, by majority vote per cell.

        Each cell takes the label held by the most raster pixels inside its
        polygon. Deliberately not the area-weighted MEAN that
        :func:`~myflopy.modflow.mf6.grid.surfaces.get_raster_vals_at_centroids`
        uses by default: averaging category ids is meaningless -- zones 1 and 3
        would average to 2, a zone the cell does not touch and may not exist.

        Cells covering no pixel (smaller than the raster resolution, or outside
        coverage) get ``nodata_zone``.
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
        numbered ``1..n`` in row order. Cells in no polygon get ``nodata_zone``;
        cells in several take the first match, which is why overlapping zone
        polygons are a modelling mistake rather than something to resolve here.
        """

        import geopandas as gpd

        frame = polygons if isinstance(polygons, gpd.GeoDataFrame) else gpd.read_file(polygons)
        source = f"polygons ({column or 'row order'})"
        if column is not None and column not in frame.columns:
            raise KeyError(
                f"Zone column {column!r} is not in the polygon layer; it has "
                f"{sorted(frame.columns)}."
            )
        labels = (
            frame[column].to_numpy() if column is not None
            else np.arange(1, len(frame) + 1)
        )

        cells = model.vor.gdf_vorPolys
        if frame.crs is not None and cells.crs is not None and frame.crs != cells.crs:
            # Reprojection preserves row order, so `labels` still lines up.
            frame = frame.to_crs(cells.crs)

        values = np.full(len(cells), int(nodata_zone), dtype=int)
        centroids = gpd.GeoDataFrame(geometry=cells.geometry.centroid, crs=cells.crs)
        joined = gpd.sjoin(centroids, frame.assign(_zone=labels)[["_zone", frame.geometry.name]],
                           how="left", predicate="within")
        # A centroid inside two overlapping polygons produces two rows; keep the
        # first so the result stays one label per cell.
        joined = joined[~joined.index.duplicated(keep="first")]
        matched = joined["_zone"].notna().to_numpy()
        values[matched] = joined.loc[matched, "_zone"].to_numpy().astype(int)
        return cls(values=values, source=source)


def _majority_sample(raster_path: Path, model, *, nodata_zone: int) -> np.ndarray:
    """Per-cell majority label of a categorical raster.

    Mirrors ``grid/surfaces.py``'s ``_area_weighted_sample`` -- burn cell ids
    onto the raster grid, then reduce per cell -- but reduces by MODE instead of
    mean, which is the only meaningful reduction for categories.
    """

    import rasterio
    from rasterio.errors import WindowError
    from rasterio.features import rasterize
    from rasterio.windows import Window, from_bounds

    cells = model.vor.gdf_vorPolys
    ncpl = len(cells)
    values = np.full(ncpl, int(nodata_zone), dtype=int)

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
        if window is None or window.width < 1 or window.height < 1:
            return values

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
        valid = flat_cells > 0
        if src.nodata is not None:
            valid &= flat_band != src.nodata
        if not valid.any():
            return values

        cell_ids = flat_cells[valid] - 1
        zone_ids = flat_band[valid].astype(np.int64)
        # Mode per cell in one pass: count every (cell, zone) pair, then take
        # each cell's most frequent zone. Offsetting by the minimum keeps
        # negative category ids usable as bincount indices.
        offset = zone_ids.min()
        span = int(zone_ids.max() - offset) + 1
        counts = np.bincount(cell_ids * span + (zone_ids - offset),
                             minlength=ncpl * span).reshape(ncpl, span)
        covered = counts.sum(axis=1) > 0
        values[covered] = counts[covered].argmax(axis=1) + offset
    return values


def resolve_zone_array(zones: Any, *, family: str, model) -> np.ndarray:
    """Normalize a zone definition into the shape pyEMU wants for ``family``.

    ``zones`` may be a :class:`ZoneSpec` or any array-like. The result is
    per-cell ``(ncpl,)`` for array targets and ``(nlay, ncpl)`` for list
    targets. See the module docstring for why the two differ and why a
    single-layer list target is padded.
    """

    spec = zones if isinstance(zones, ZoneSpec) else ZoneSpec.from_array(zones)
    values = spec.values
    ncpl = int(model.vor.ncpl)
    nlay = int(model.gwf.modelgrid.nlay)

    flat = values.reshape(-1)
    if values.ndim >= 2 and values.shape[-1] == ncpl and values.shape[0] not in (1, ncpl):
        # Per-layer input: zones are a property of the map, so require the
        # layers to agree rather than silently taking the first.
        if not np.all(values == values[0]):
            raise ValueError(
                f"Zone definition from {spec.source} differs between layers. "
                "Zones are per-cell labels applied to every layer the target "
                "covers; use `layers=` to restrict which layers are "
                "parameterized instead."
            )
        flat = values[0]
    elif flat.size != ncpl:
        raise ValueError(
            f"Zone definition from {spec.source} has {flat.size} values, but the "
            f"grid has {ncpl} cells. Zones are one integer label per cell."
        )

    per_cell = np.asarray(flat, dtype=int).reshape(ncpl)
    if family == "array":
        return per_cell

    # List targets are indexed (layer, cell), so the array must be 2-D -- and at
    # least 2 rows deep even on a one-layer model, because pyEMU reshapes any
    # (1, n) array on a vertex grid assuming it is an idomain array
    # (pst_from.py checker2). The padding row addresses no real cell, so nothing
    # ever indexes it.
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
