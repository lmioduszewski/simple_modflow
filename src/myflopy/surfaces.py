"""Layer top/bottom surfaces for MODFLOW 6 discretization.

A :class:`Surface` describes one elevation surface -- a raster, a flat constant,
a surface offset below the one above it, an interpolation from drawn elevation
contours (via GRASS), or interpolation from scattered points. :class:`LayerSurfaces`
stacks them into ordered model layers, samples them onto a grid's cells, and
reconciles overlaps into DISV-ready top/bottom elevations. Thin layers can
optionally be pinched out (``idomain = -1``, vertical pass-through).

The grid does not own this logic. ``LayerSurfaces.attach(grid)`` writes the
result back to ``grid.gdf_topbtm`` so existing consumers (choropleth mapping,
layer-elevation lookups, ``mf.disv``) keep working unchanged.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import geopandas as gpd
import numpy as np

from myflopy.modflow.mf6.grid.geometry import reconcile_surfaces
from myflopy.modflow.mf6.grid.surfaces import get_raster_vals_at_centroids

# Surfaces defined relative to the surface immediately above them.
_RELATIVE_KINDS = ("offset_below", "constant_thickness")

# Pinch-out policies for thin cells (thickness < minimum_thickness):
#   passthrough -> idomain = -1 (cell removed, vertical flow passes through)
#   inactive    -> idomain = 0  (cell removed, blocks vertical flow)
#   floor       -> stay active  (reconcile enforces minimum spacing)
_PINCH_POLICIES = ("passthrough", "inactive", "floor")


def _as_per_layer(value, nlay: int, name: str) -> list:
    """Broadcast a scalar to one value per layer, or validate a per-layer sequence."""

    if np.isscalar(value):
        return [value] * nlay
    seq = list(value)
    if len(seq) != nlay:
        raise ValueError(
            f"{name} has {len(seq)} entries but there are {nlay} layers."
        )
    return seq


# Length-unit handling. Elevations are converted by a single scale factor (same
# datum, different unit), so only a per-unit "value in feet" table is needed.
_LENGTH_TO_FEET = {
    "feet": 1.0, "ft": 1.0, "foot": 1.0, "us-ft": 1.0, "usft": 1.0,
    "meters": 3.280839895013123, "meter": 3.280839895013123,
    "metres": 3.280839895013123, "m": 3.280839895013123,
    "centimeters": 0.03280839895013123, "cm": 0.03280839895013123,
}
# Model length unit -> MF6 LENGTH_UNITS token.
_MF6_LENGTH_UNITS = {"feet": "FEET", "meters": "METERS", "centimeters": "CENTIMETERS"}


def _length_factor(from_units: str, to_units: str) -> float:
    """Multiplier converting an elevation from ``from_units`` to ``to_units``."""

    try:
        return _LENGTH_TO_FEET[from_units.lower()] / _LENGTH_TO_FEET[to_units.lower()]
    except KeyError as exc:
        raise ValueError(
            f"Unknown length unit {exc.args[0]!r}; known: {sorted(_LENGTH_TO_FEET)}."
        ) from None


def _cells_in_zone(vor, zone) -> np.ndarray:
    """Boolean per-cell mask: True where the cell centroid lies inside ``zone``.

    ``zone`` is a shapely geometry (in the grid CRS) or a path to a vector file.
    """

    if isinstance(zone, (str, Path)):
        gdf = gpd.read_file(zone)
        vor_crs = getattr(vor, "crs", None)
        if vor_crs is not None and gdf.crs is not None:
            gdf = gdf.to_crs(vor_crs)
        geom = gdf.union_all()
    else:
        geom = zone
    return vor.gdf_vorPolys.geometry.centroid.within(geom).to_numpy()


@dataclass(frozen=True)
class Surface:
    """One elevation surface, sampled onto a grid's cells when needed.

    Build one with a constructor rather than the raw fields: :meth:`raster`,
    :meth:`flat`, :meth:`offset_below`, :meth:`constant_thickness`,
    :meth:`from_contours`, or :meth:`from_points`.
    """

    kind: str
    value: float | None = None
    path: Path | None = None
    fill: str | None = None
    units: str | None = None  # source units, if different from the model units
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
    # direct per-cell elevation array (e.g. an existing model's top/botm)
    array: Any = None
    # surface-algebra / isopach parameters
    operands: tuple = ()
    zone: Any = None

    @classmethod
    def raster(
        cls, path: Path | str, *, fill: str | None = None, units: str | None = None
    ) -> Surface:
        """A surface sampled from an existing raster.

        ``fill="propagate"`` carries the surface above down into cells where the
        raster has no data (nodata or outside coverage), so the layer pinches
        against the surface above instead of leaving a ``NaN`` elevation.
        ``units`` (e.g. ``"meters"``) converts the sampled elevations to the
        model length units when sampled through a stack that declares them.
        """

        return cls(kind="raster", path=Path(path), fill=fill, units=units)

    @classmethod
    def flat(cls, value: float, *, units: str | None = None) -> Surface:
        """A flat, constant-elevation surface (replaces the magic integer)."""

        return cls(kind="flat", value=float(value), units=units)

    @classmethod
    def offset_below(cls, distance: float) -> Surface:
        """A surface a fixed vertical distance below the surface above it."""

        return cls(kind="offset_below", value=float(distance))

    @classmethod
    def constant_thickness(cls, thickness: float) -> Surface:
        """A surface giving the layer above it a fixed thickness.

        Numerically identical to :meth:`offset_below`, but named for the common
        "this layer is N units thick" intent.
        """

        return cls(kind="constant_thickness", value=float(thickness))

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
        units: str | None = None,
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
            units=units,
            grass_kwargs=grass_kwargs,
        )

    @classmethod
    def from_points(
        cls, xs, ys, zs, *, method: str = "linear", units: str | None = None
    ) -> Surface:
        """A surface interpolated from scattered (x, y, z) points (SciPy)."""

        return cls(
            kind="points",
            points=(np.asarray(xs), np.asarray(ys), np.asarray(zs)),
            method=method,
            units=units,
        )

    @classmethod
    def from_array(cls, values, *, units: str | None = None) -> Surface:
        """A surface whose per-cell elevations are supplied directly.

        ``values`` must have one entry per grid cell (``len == ncpl``); use this
        to reuse an existing model's top/botm arrays on the *same* grid. To move
        an existing model's surfaces onto a *different* grid, interpolate instead
        (:meth:`from_points`, or :meth:`LayerStack.from_modflow`).
        """

        return cls(kind="array", array=np.asarray(values, dtype=float), units=units)

    @staticmethod
    def _coerce(x) -> Surface:
        """Coerce a Surface / path / number operand into a Surface."""
        if isinstance(x, Surface):
            return x
        if isinstance(x, (int, float)):
            return Surface.flat(float(x))
        if isinstance(x, (str, Path)):
            return Surface.raster(x)
        raise TypeError(f"Expected a Surface, path, or number; got {type(x).__name__}.")

    @staticmethod
    def _reject_lone_envelope(name: str, fluent: str, surfaces: tuple) -> None:
        """Refuse a one-surface envelope -- nearly always a misbound instance call.

        ``minimum``/``maximum`` are classmethods taking ``*surfaces``, so
        ``a.maximum(b)`` binds ``a`` to ``cls`` and arrives here as the single
        operand ``b``: a plausible-looking Surface that quietly ignores ``a``,
        and downstream a wrong contact elevation nobody sees until the heads
        look odd. The variadic signature is what hides it -- ``a.shift(-25)``
        raises for the missing ``distance``, but ``*surfaces`` swallows the
        arity error. An envelope of one surface is a no-op even when it is
        deliberate, so it costs nothing to make it the tripwire.
        """

        if len(surfaces) >= 2:
            return
        raise TypeError(
            f"Surface.{name}() takes two or more surfaces, got {len(surfaces)}. "
            f"If you wrote `a.{name}(b)`: {name} is a classmethod, so `a` binds to "
            f"`cls` and is dropped -- the result would be `b` alone. Write "
            f"`Surface.{name}(a, b)`, or the instance-bound `a.{fluent}(b)`."
        )

    @classmethod
    def minimum(cls, *surfaces) -> Surface:
        """Per-cell minimum (lower envelope) of two or more surfaces.

        A classmethod: call it as ``Surface.minimum(a, b)``, or reach for the
        instance-bound :meth:`capped_at`. One surface raises :class:`TypeError`,
        because that is what a misbound ``a.minimum(b)`` looks like from here.
        """
        cls._reject_lone_envelope("minimum", "capped_at", surfaces)
        return cls(kind="min", operands=tuple(cls._coerce(s) for s in surfaces))

    @classmethod
    def maximum(cls, *surfaces) -> Surface:
        """Per-cell maximum (upper envelope) of two or more surfaces.

        A classmethod: call it as ``Surface.maximum(a, b)``, or reach for the
        instance-bound :meth:`floored_at`. One surface raises :class:`TypeError`,
        because that is what a misbound ``a.maximum(b)`` looks like from here.
        """
        cls._reject_lone_envelope("maximum", "floored_at", surfaces)
        return cls(kind="max", operands=tuple(cls._coerce(s) for s in surfaces))

    @classmethod
    def clamp(cls, surface, *, lower=None, upper=None) -> Surface:
        """Constrain ``surface`` to stay within ``[lower, upper]`` (per cell).

        ``lower``/``upper`` may be surfaces, paths, or constants; either may be
        a relative surface (e.g. ``upper=Surface.offset_below(5)`` to keep this
        surface at least 5 below the one above). At least one is required: with
        neither, the result is ``surface`` unchanged, which is also exactly what
        a misbound ``a.clamp(b)`` produces, so both raise :class:`TypeError`.
        """
        if lower is None and upper is None:
            raise TypeError(
                "Surface.clamp() needs lower= or upper=; with neither it would "
                "return `surface` unchanged. If you wrote `a.clamp(b)`: clamp is a "
                "classmethod, so `a` binds to `cls` and is dropped. Write "
                "`Surface.clamp(a, lower=..., upper=...)`, or the instance-bound "
                "`a.between(lower=..., upper=...)`."
            )
        return cls(
            kind="clamp",
            operands=(
                cls._coerce(surface),
                None if lower is None else cls._coerce(lower),
                None if upper is None else cls._coerce(upper),
            ),
        )

    @classmethod
    def where(cls, zone, inside, outside) -> Surface:
        """Use ``inside`` for cells whose centroid is in ``zone``, else ``outside``.

        ``zone`` is a shapely geometry or a path to a vector file.
        """
        return cls(
            kind="where", zone=zone, operands=(cls._coerce(inside), cls._coerce(outside))
        )

    @classmethod
    def isopach(cls, source) -> Surface:
        """A layer bottom defined by a *thickness* map subtracted from above."""
        return cls(kind="isopach", operands=(cls._coerce(source),))

    @classmethod
    def shift(cls, source, distance: float) -> Surface:
        """``source`` displaced vertically by ``distance`` (negative = downward)."""
        return cls(kind="shift", value=float(distance), operands=(cls._coerce(source),))

    # -- fluent surface algebra ------------------------------------------- #
    # Read base-surface-first, left to right: ``Raster(bedrock).capped_at(ground - 5)``.
    # These are thin wrappers over the constructors above (same engine), so
    # ``Min``/``Max``/``Clamp``/``Where`` keep working unchanged.
    def capped_at(self, other) -> Surface:
        """This surface, but never *above* ``other`` (per cell) -- a ceiling.

        The fluent form of :meth:`minimum` (lower envelope); e.g. a layer top
        ``capped_at`` an eroding valley floor."""
        return Surface.minimum(self, other)

    def floored_at(self, other) -> Surface:
        """This surface, but never *below* ``other`` (per cell) -- a floor.

        The fluent form of :meth:`maximum` (upper envelope)."""
        return Surface.maximum(self, other)

    def between(self, *, lower=None, upper=None) -> Surface:
        """This surface constrained to ``[lower, upper]`` per cell (fluent :meth:`clamp`)."""
        return Surface.clamp(self, lower=lower, upper=upper)

    def below(self, distance: float) -> Surface:
        """A surface ``distance`` units *below* this one (same as ``self - distance``)."""
        return Surface.shift(self, -float(distance))

    def above(self, distance: float) -> Surface:
        """A surface ``distance`` units *above* this one (same as ``self + distance``)."""
        return Surface.shift(self, float(distance))

    def within(self, zone, *, outside) -> Surface:
        """Use this surface for cells inside ``zone``, ``outside`` elsewhere.

        The fluent form of :meth:`where`; e.g.
        ``Flat(15).within(channel, outside=Flat(2))``."""
        return Surface.where(zone, self, outside)

    def thickness(self) -> Surface:
        """A layer bottom set by this *thickness* map below the surface above
        (fluent :meth:`isopach`)."""
        return Surface.isopach(self)

    def __sub__(self, distance) -> Surface:
        """``surface - distance``: a new surface shifted down by a constant elevation."""

        if not isinstance(distance, (int, float, np.integer, np.floating)):
            return NotImplemented
        return Surface.shift(self, -float(distance))

    def __add__(self, distance) -> Surface:
        """``surface + distance``: a new surface shifted up by a constant elevation."""

        if not isinstance(distance, (int, float, np.integer, np.floating)):
            return NotImplemented
        return Surface.shift(self, float(distance))

    __radd__ = __add__

    @property
    def is_relative(self) -> bool:
        """True if this surface is defined relative to the surface above it."""

        return self.kind in _RELATIVE_KINDS

    @property
    def is_derived(self) -> bool:
        """True if this surface resolves to a cached, regenerable raster."""

        return self.kind == "contours"

    def _derived_raster(self):
        """Return the :class:`DerivedRaster` cache handle for a derived surface."""

        from myflopy.modflow.utils.derived_raster import DerivedRaster

        def produce() -> None:
            """Regenerate the derived raster by interpolating this surface's contours."""

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

        sources = [self.contours]
        if self.region_raster is not None:
            sources.append(self.region_raster)
        params = {
            "z": self.z,
            "resolution": self.resolution,
            "epsg": self.epsg,
            "grass_kwargs": {k: str(v) for k, v in self.grass_kwargs.items()},
        }
        return DerivedRaster(self.out, sources, params, produce)

    def cache_status(self) -> str | None:
        """``"missing"``/``"fresh"``/``"stale"`` for derived surfaces, else ``None``."""

        return self._derived_raster().status() if self.is_derived else None

    def resolve_source(self, *, refresh: bool = False, warn: bool = True) -> Path | float:
        """Return a raster path or flat constant for centroid sampling.

        Derived surfaces (contours) interpolate to a cached raster on first use;
        a stale cache is reused with a warning unless ``refresh=True``.
        """

        if self.kind == "flat":
            return self.value
        if self.kind == "raster":
            return self.path
        if self.kind == "contours":
            return self._derived_raster().ensure(refresh=refresh, warn=warn)
        raise TypeError(f"Surface kind '{self.kind}' has no raster/constant source.")

    def values(
        self,
        vor,
        previous: np.ndarray | None = None,
        *,
        method: str = "area",
        refresh: bool = False,
    ) -> np.ndarray:
        """Return per-cell elevation values sampled on ``vor``.

        ``previous`` is the surface immediately above (or ``None`` for the model
        top). Relative surfaces (:meth:`offset_below`, :meth:`constant_thickness`)
        and raster ``fill="propagate"`` require it. ``method`` controls raster
        sampling: ``"area"`` (area-weighted, default) or ``"centroid"``.
        ``refresh`` rebuilds a derived (contour) surface's cached raster.
        """

        def _op(operand):
            """Sample one operand surface on ``vor`` as a float array (same args as the parent)."""

            return np.asarray(
                operand.values(vor, previous, method=method, refresh=refresh),
                dtype=float,
            )

        if self.kind in ("min", "max"):
            stacked = np.vstack([_op(op) for op in self.operands])
            return np.nanmin(stacked, axis=0) if self.kind == "min" else np.nanmax(stacked, axis=0)

        if self.kind == "shift":
            return _op(self.operands[0]) + float(self.value)

        if self.kind == "clamp":
            base, lower, upper = self.operands
            vals = _op(base)
            if lower is not None:
                vals = np.maximum(vals, _op(lower))
            if upper is not None:
                vals = np.minimum(vals, _op(upper))
            return vals

        if self.kind == "where":
            inside, outside = self.operands
            return np.where(_cells_in_zone(vor, self.zone), _op(inside), _op(outside))

        if self.kind == "isopach":
            if previous is None:
                raise ValueError("Surface.isopach cannot be the model top.")
            return previous - _op(self.operands[0])

        if self.is_relative:
            if previous is None:
                raise ValueError(
                    f"Surface.{self.kind} cannot be the model top -- there is no "
                    "surface above to offset from."
                )
            return previous - float(self.value)

        if self.kind == "flat":
            return np.full(len(np.asarray(vor.centroids[0])), float(self.value))

        if self.kind == "array":
            vals = np.asarray(self.array, dtype=float).ravel()
            n = len(np.asarray(vor.centroids[0]))
            if vals.size != n:
                raise ValueError(
                    f"from_array surface has {vals.size} values but the grid has {n} "
                    "cells; use Surface.from_points to interpolate onto another grid."
                )
            return vals

        if self.kind == "points":
            from scipy.interpolate import griddata

            cx = np.asarray(vor.centroids[0])
            cy = np.asarray(vor.centroids[1])
            xs, ys, zs = self.points
            return griddata((xs, ys), zs, (cx, cy), method=self.method)

        gdf = get_raster_vals_at_centroids(
            vor, [self.resolve_source(refresh=refresh)], [0], method=method
        )
        vals = gdf[0].to_numpy(dtype=float)
        if self.fill == "propagate":
            mask = np.isnan(vals)
            if mask.any():
                if previous is None:
                    raise ValueError(
                        "fill='propagate' cannot be used on the model top -- there "
                        "is no surface above to inherit from."
                    )
                vals = vals.copy()
                vals[mask] = previous[mask]
        elif self.fill is not None:
            raise ValueError(
                f"Unknown fill mode {self.fill!r}; expected None or 'propagate'."
            )
        return vals


class LayerSurfaces:
    """An ordered stack of layer surfaces sampled onto a grid's cells.

    The first surface is the model top; each subsequent surface is the bottom of
    a layer. ``sample`` returns a centroid GeoDataFrame with one column per
    surface (matching ``gdf_topbtm``); ``top_botm`` returns DISV-ready arrays.
    """

    def __init__(self, surfaces: Sequence[Surface], labels: Sequence | None = None):
        """Wrap an ordered surface stack (top first) with matching column ``labels``.

        ``labels`` default to integer positions; their count must match
        ``surfaces``.
        """

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
        method: str = "area",
        length_units: str | None = None,
        refresh: bool = False,
    ) -> gpd.GeoDataFrame:
        """Sample every surface onto ``vor`` cells and reconcile overlaps.

        Surfaces are resolved top-down so that each one can depend on the surface
        above it (relative surfaces and ``fill="propagate"``). With
        ``reconcile=True`` (default), a flat surface stays flat where it fits and
        is lowered where it would intrude on the surface above -- the "flat where
        possible, fit between" behavior.

        Parameters
        ----------
        vor : VoronoiGridPlus
            The grid to sample onto.
        reconcile : bool, default True
            Enforce layer ordering where surfaces cross or crowd. ``False``
            returns the raw sampled surfaces, crossings and all.
        min_sep : float, default 0.1
            The separation reconcile enforces where it acts.
        trigger_sep : float, default 1
            How close counts as a conflict. Surfaces nearer than this are
            separated even if they never actually cross, so the default treats
            "within 1 length unit" as touching.
        which : {'bottom', 'top'}, default 'bottom'
            Which side gives way. ``'bottom'`` lowers the deeper surface to
            ``upper - min_sep`` and cascades downward, leaving the model top
            untouched. ``'top'`` raises the shallower one to ``lower + min_sep``
            and propagates upward, so it CAN MOVE THE MODEL TOP -- with
            ``top=100, a=105`` it returns ``top=105.1``. Prefer ``'bottom'``
            unless the deeper contact is the one you trust.
        method : {'area', 'centroid'}, default 'area'
            Raster sampling: area-weighted over the cell, or a single centroid
            lookup.
        length_units : str, optional
            The model length unit. A surface declaring different ``units`` is
            converted to it; relative surfaces cannot be, and raise.
        refresh : bool, default False
            Rebuild any derived (contour) caches before sampling.

        Returns
        -------
        geopandas.GeoDataFrame
            Cell centroids with one column per surface, in stack order.
        """

        columns: dict = {}
        previous: np.ndarray | None = None
        for label, surface in zip(self.labels, self.surfaces, strict=False):
            vals = np.asarray(
                surface.values(vor, previous=previous, method=method, refresh=refresh),
                dtype=float,
            )
            if length_units is not None and surface.units:
                if surface.is_relative:
                    raise ValueError(
                        "units= is only supported on absolute surfaces, not "
                        "relative (offset_below / constant_thickness)."
                    )
                vals = vals * _length_factor(surface.units, length_units)
            columns[label] = vals
            previous = vals

        gdf = gpd.GeoDataFrame(
            columns, geometry=vor.gdf_vorPolys.centroid.values, crs=vor.crs
        )
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

    @staticmethod
    def _thickness(gdf: gpd.GeoDataFrame) -> np.ndarray:
        """Per-layer thickness ``(nlay, n_cells)`` from a sampled centroid frame."""

        columns = [c for c in gdf.columns if c != "geometry"]
        elev = gdf[columns].to_numpy().T  # (n_surfaces, n_cells)
        return elev[:-1] - elev[1:]

    @staticmethod
    def _idomain_from_thickness(
        thickness: np.ndarray, minimum_thickness=1.0, pinch="passthrough"
    ) -> np.ndarray:
        """Per-layer idomain from thickness, applying each layer's pinch policy.

        ``minimum_thickness`` and ``pinch`` may each be a scalar (applied to
        every layer) or a per-layer sequence. Policies: ``"passthrough"`` ->
        ``-1`` (vertical pass-through), ``"inactive"`` -> ``0`` (blocks flow),
        ``"floor"`` -> stay active (reconcile enforces minimum spacing).
        """

        nlay = thickness.shape[0]
        min_thk = _as_per_layer(minimum_thickness, nlay, "minimum_thickness")
        policy = _as_per_layer(pinch, nlay, "pinch")
        idomain = np.ones(thickness.shape, dtype=int)
        for i in range(nlay):
            if policy[i] not in _PINCH_POLICIES:
                raise ValueError(
                    f"Unknown pinch policy {policy[i]!r}; expected one of "
                    f"{_PINCH_POLICIES}."
                )
            thin = thickness[i] < float(min_thk[i])
            if policy[i] == "passthrough":
                idomain[i, thin] = -1
            elif policy[i] == "inactive":
                idomain[i, thin] = 0
            # "floor": leave active (idomain stays 1)
        return idomain

    @staticmethod
    def _validate_pinch_invariant(
        minimum_thickness, pinch, nlay: int, sample_kwargs: dict
    ) -> None:
        """Guard ``min_sep < minimum_thickness`` for each *pinching* layer.

        Reconcile spaces layers ``min_sep`` apart, so if a pinching layer's
        threshold is <= ``min_sep`` nothing would ever pinch -- a silent no-op.
        ``"floor"`` layers do not pinch and are exempt.
        """

        if not sample_kwargs.get("reconcile", True):
            return
        min_sep = sample_kwargs.get("min_sep", 0.1)
        min_thk = _as_per_layer(minimum_thickness, nlay, "minimum_thickness")
        policy = _as_per_layer(pinch, nlay, "pinch")
        for i in range(nlay):
            if policy[i] in ("passthrough", "inactive") and not (
                min_sep < float(min_thk[i])
            ):
                raise ValueError(
                    f"layer {i}: min_sep ({min_sep}) must be < its minimum_thickness "
                    f"({min_thk[i]}) for pinch policy {policy[i]!r} to take effect; "
                    "reconcile would otherwise inflate thin cells above the threshold."
                )

    def top_botm(self, vor, **kwargs) -> tuple[np.ndarray, np.ndarray]:
        """Return ``(top, botm)`` elevation arrays ready for ``mf.disv``."""

        return self._split_top_botm(self.sample(vor, **kwargs))

    def top_botm_idomain(
        self,
        vor,
        *,
        minimum_thickness=1.0,
        pinch="passthrough",
        **sample_kwargs,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return ``(top, botm, idomain)`` with thin cells pinched out.

        ``minimum_thickness`` and ``pinch`` may be scalars (every layer) or
        per-layer sequences. Default policy ``"passthrough"`` gives thin cells
        ``idomain = -1`` (MODFLOW 6 vertical pass-through).
        """

        gdf = self.sample(vor, **sample_kwargs)
        top, botm = self._split_top_botm(gdf)
        thickness = self._thickness(gdf)
        self._validate_pinch_invariant(
            minimum_thickness, pinch, thickness.shape[0], sample_kwargs
        )
        idomain = self._idomain_from_thickness(thickness, minimum_thickness, pinch)
        return top, botm, idomain

    def thickness_report(
        self, vor, *, minimum_thickness: float = 1.0, **sample_kwargs
    ) -> str:
        """Per-layer thickness summary and thin-cell counts (for diagnostics)."""

        thickness = self._thickness(self.sample(vor, **sample_kwargs))
        nlay, ncell = thickness.shape
        min_thk = _as_per_layer(minimum_thickness, nlay, "minimum_thickness")
        lines = [f"LayerSurfaces: {nlay} layers, {ncell} cells"]
        for i in range(nlay):
            t = thickness[i]
            thin = int((t < float(min_thk[i])).sum())
            lines.append(
                f"  layer {i}: thickness min={t.min():.2f} mean={t.mean():.2f} "
                f"max={t.max():.2f}  thin={thin} (<{min_thk[i]})"
            )
        return "\n".join(lines)

    def to_disv(
        self,
        vor,
        *,
        nlay: int | None = None,
        idomain: Any | None = None,
        pinch_out: bool = False,
        minimum_thickness=1.0,
        pinch="passthrough",
        length_units: str | None = None,
        name: str = "disv",
        attach: bool = True,
        **sample_kwargs,
    ):
        """Build a ready-to-use ``mf.disv`` package spec from ``vor`` + surfaces.

        Samples the surfaces onto ``vor`` (reconciling overlaps by default),
        takes the cell geometry from ``vor.get_disv_gridprops()``, and returns a
        :class:`~myflopy.specs.PackageSpec` for DISV. With ``pinch_out=True``,
        cells thinner than ``minimum_thickness`` get ``idomain = -1`` (vertical
        pass-through). With ``attach=True`` (default) the sampled elevations are
        also written to ``vor.gdf_topbtm`` so choropleth/mapping code keeps
        working.
        """

        from myflopy.package_api import disv  # lazy: avoid an import cycle

        gdf = self.sample(vor, length_units=length_units, **sample_kwargs)
        if attach:
            vor.gdf_topbtm = gdf
        top, botm = self._split_top_botm(gdf)
        if pinch_out:
            if idomain is not None:
                raise ValueError(
                    "Pass either pinch_out=True or an explicit idomain, not both."
                )
            thickness = self._thickness(gdf)
            self._validate_pinch_invariant(
                minimum_thickness, pinch, thickness.shape[0], sample_kwargs
            )
            idomain = self._idomain_from_thickness(
                thickness, minimum_thickness, pinch
            )
        props = vor.get_disv_gridprops()
        disv_options = {}
        if length_units is not None:
            disv_options["length_units"] = _MF6_LENGTH_UNITS.get(
                length_units.lower(), length_units.upper()
            )
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
            **disv_options,
        )


__all__ = ["Surface", "LayerSurfaces"]
