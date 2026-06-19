"""Layer top/bottom surfaces for MODFLOW 6 discretization.

A :class:`Surface` describes one elevation surface -- a raster, a flat constant,
an interpolation from drawn elevation contours (via GRASS), or interpolation
from scattered points. :class:`LayerSurfaces` stacks them into ordered model
layers, samples them onto a grid's cells, and reconciles overlaps into
DISV-ready top/bottom elevations.

The grid does not own this logic. ``LayerSurfaces.attach(grid)`` writes the
result back to ``grid.gdf_topbtm`` so existing consumers (choropleth mapping,
layer-elevation lookups, ``mf.disv``) keep working unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

import geopandas as gpd
import numpy as np

from myflopy.modflow.mf6.grid.surfaces import get_raster_vals_at_centroids
from myflopy.modflow.mf6.grid.geometry import reconcile_surfaces


@dataclass(frozen=True)
class Surface:
    """One elevation surface, sampled onto a grid's cells when needed.

    Build one with a constructor rather than the raw fields:
    :meth:`raster`, :meth:`flat`, :meth:`from_contours`, or :meth:`from_points`.
    """

    kind: str
    value: float | None = None
    path: Path | None = None
    # contour-interpolation parameters
    contours: Path | None = None
    z: str = "Elev"
    region_raster: Path | None = None
    resolution: float = 4
    out: Path | None = None
    epsg: str = "2927"
    grass_kwargs: dict = field(default_factory=dict)
    # point-interpolation parameters
    points: tuple[Any, Any, Any] | None = None
    method: str = "linear"

    @classmethod
    def raster(cls, path: Path | str) -> Surface:
        """A surface sampled from an existing raster."""

        return cls(kind="raster", path=Path(path))

    @classmethod
    def flat(cls, value: float) -> Surface:
        """A flat, constant-elevation surface (replaces the magic integer)."""

        return cls(kind="flat", value=float(value))

    @classmethod
    def from_contours(
        cls,
        contours: Path | str,
        *,
        z: str = "Elev",
        region_raster: Path | str | None = None,
        resolution: float = 4,
        out: Path | str | None = None,
        epsg: str = "2927",
        **grass_kwargs: Any,
    ) -> Surface:
        """A surface interpolated from drawn elevation contours via GRASS.

        The interpolated raster is written to ``out`` (default: alongside the
        contours with a ``.interp.tif`` suffix) the first time the surface is
        resolved, then reused -- interpolate once, sample many times.
        """

        contours = Path(contours)
        out = Path(out) if out is not None else contours.with_suffix(".interp.tif")
        return cls(
            kind="contours",
            contours=contours,
            z=z,
            region_raster=region_raster,
            resolution=resolution,
            out=out,
            epsg=epsg,
            grass_kwargs=grass_kwargs,
        )

    @classmethod
    def from_points(cls, xs, ys, zs, *, method: str = "linear") -> Surface:
        """A surface interpolated from scattered (x, y, z) points (SciPy)."""

        return cls(
            kind="points",
            points=(np.asarray(xs), np.asarray(ys), np.asarray(zs)),
            method=method,
        )

    def resolve_source(self) -> Path | float:
        """Return a raster path or flat constant for centroid sampling.

        Interpolates contours to a cached raster on first use.
        """

        if self.kind == "flat":
            return self.value
        if self.kind == "raster":
            return self.path
        if self.kind == "contours":
            if not self.out.exists():
                from myflopy.modflow.utils.contour_interp import (
                    interpolate_contours_to_raster,
                )

                interpolate_contours_to_raster(
                    self.contours,
                    out=self.out,
                    z_field=self.z,
                    region_raster=self.region_raster,
                    resolution=self.resolution,
                    epsg=self.epsg,
                    **self.grass_kwargs,
                )
            return self.out
        raise TypeError(f"Surface kind '{self.kind}' has no raster/constant source.")

    def values(self, vor) -> np.ndarray:
        """Return per-cell elevation values sampled on ``vor``."""

        if self.kind == "points":
            from scipy.interpolate import griddata

            cx = np.asarray(vor.centroids[0])
            cy = np.asarray(vor.centroids[1])
            xs, ys, zs = self.points
            return griddata((xs, ys), zs, (cx, cy), method=self.method)
        gdf = get_raster_vals_at_centroids(vor, [self.resolve_source()], [0])
        return gdf[0].to_numpy()


class LayerSurfaces:
    """An ordered stack of layer surfaces sampled onto a grid's cells.

    The first surface is the model top; each subsequent surface is the bottom of
    a layer. ``sample`` returns a centroid GeoDataFrame with one column per
    surface (matching ``gdf_topbtm``); ``top_botm`` returns DISV-ready arrays.
    """

    def __init__(self, surfaces: Sequence[Surface], labels: Sequence | None = None):
        self.surfaces = list(surfaces)
        self.labels = (
            list(labels) if labels is not None else list(range(len(self.surfaces)))
        )
        if len(self.labels) != len(self.surfaces):
            raise ValueError("labels length must match surfaces length")

    def sample(
        self,
        vor,
        *,
        reconcile: bool = True,
        min_sep: float = 0.1,
        trigger_sep: float = 1,
        which: str = "bottom",
    ) -> gpd.GeoDataFrame:
        """Sample every surface onto ``vor`` cells and reconcile overlaps.

        Raster/flat/contour surfaces are sampled in one pass; point surfaces are
        interpolated separately. With ``reconcile=True`` (default), a flat
        surface stays flat where it fits and is lowered where it would intrude
        on the surface above -- the "flat where possible, fit between" behavior.
        """

        sources: list[Path | float] = []
        source_labels: list = []
        point_items: list = []
        for label, surface in zip(self.labels, self.surfaces):
            if surface.kind == "points":
                point_items.append((label, surface))
            else:
                sources.append(surface.resolve_source())
                source_labels.append(label)

        if sources:
            gdf = get_raster_vals_at_centroids(vor, sources, source_labels)
        else:
            gdf = gpd.GeoDataFrame(
                geometry=vor.gdf_vorPolys.centroid.values, crs=vor.crs
            )
        for label, surface in point_items:
            gdf[label] = surface.values(vor)

        gdf = gdf[["geometry"] + self.labels]
        if reconcile:
            df = reconcile_surfaces(
                vor, df=gdf, min_sep=min_sep, trigger_sep=trigger_sep, which=which
            )
            gdf = gpd.GeoDataFrame(df, geometry=gdf.geometry.values, crs=gdf.crs)
        return gdf

    def attach(self, vor, **kwargs) -> gpd.GeoDataFrame:
        """Sample and store the result on ``vor.gdf_topbtm`` for downstream use."""

        gdf = self.sample(vor, **kwargs)
        vor.gdf_topbtm = gdf
        return gdf

    @staticmethod
    def _split_top_botm(gdf: gpd.GeoDataFrame) -> tuple[np.ndarray, np.ndarray]:
        """Split a sampled centroid frame into ``(top, botm)`` arrays."""

        columns = [c for c in gdf.columns if c != "geometry"]
        elevations = gdf[columns].to_numpy().T  # (n_surfaces, n_cells)
        return elevations[0], elevations[1:]

    def top_botm(self, vor, **kwargs) -> tuple[np.ndarray, np.ndarray]:
        """Return ``(top, botm)`` elevation arrays ready for ``mf.disv``."""

        return self._split_top_botm(self.sample(vor, **kwargs))

    def to_disv(
        self,
        vor,
        *,
        nlay: int | None = None,
        idomain: Any | None = None,
        name: str = "disv",
        attach: bool = True,
        **sample_kwargs,
    ):
        """Build a ready-to-use ``mf.disv`` package spec from ``vor`` + surfaces.

        Samples the surfaces onto ``vor`` (reconciling overlaps by default),
        takes the cell geometry from ``vor.get_disv_gridprops()``, and returns a
        :class:`~myflopy.specs.PackageSpec` for DISV -- the single call that
        turns a grid and a layer stack into discretization. With ``attach=True``
        (default) the sampled elevations are also written to ``vor.gdf_topbtm``
        so choropleth/mapping code keeps working.
        """

        from myflopy.package_api import disv  # lazy: avoid an import cycle

        gdf = self.sample(vor, **sample_kwargs)
        if attach:
            vor.gdf_topbtm = gdf
        top, botm = self._split_top_botm(gdf)
        props = vor.get_disv_gridprops()
        return disv(
            nlay=len(self.surfaces) - 1 if nlay is None else nlay,
            ncpl=props["ncpl"],
            nvert=len(props["vertices"]),
            vertices=props["vertices"],
            cell2d=props["cell2d"],
            top=top,
            botm=botm,
            idomain=idomain,
            name=name,
        )


__all__ = ["Surface", "LayerSurfaces"]
