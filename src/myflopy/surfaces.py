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
from dataclasses import MISSING, dataclass, field, fields
from pathlib import Path
from typing import Any

import geopandas as gpd
import numpy as np

from myflopy._defaults import DEFAULT_EPSG
from myflopy.modflow.mf6.grid.geometry import (
    centroid_zones,
    geometry_first,
    reconcile_surfaces,
)
from myflopy.modflow.mf6.grid.surfaces import get_raster_vals_at_centroids

# Surfaces defined relative to the surface immediately above them.
_RELATIVE_KINDS = ("offset_below", "constant_thickness")

# Kinds that need `previous` but are NOT a fixed offset from it, so `values`
# handles them before the `is_relative` branch.
_PREVIOUS_KINDS = ("isopach", "toward")


def _depends_on_previous(surface) -> bool:
    """True if sampling ``surface`` consults the surface above it, at any depth.

    Composites carry it upward: ``Flat(150).capped_at(Isopach(map))`` is a `min`
    node, which looks absolute, over an operand that is not. `toward` needs this
    because a target measured from `previous` re-resolves against each cut it
    makes, drifting instead of dividing.
    """

    if surface.kind in _RELATIVE_KINDS or surface.kind in _PREVIOUS_KINDS:
        return True
    return any(
        _depends_on_previous(op) for op in surface.operands if op is not None
    )

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

    ``zone`` is a shapely geometry (in the grid CRS) or a path to a vector file,
    whose rows are unioned into one zone. The zone is CLOSED -- a centroid exactly
    on its edge is inside it, through the shared
    :func:`~myflopy.modflow.mf6.grid.geometry.centroid_zones` rule. It used to be
    ``within``'s strict inside, so a zone whose edge ran along a row of centroids
    reached none of that row: ``box(0, 0, 100, 37.5)`` on a 4x4 grid of 25 ft cells
    selected cells 0-3 and not 4-7 (ledger 243).
    """

    if isinstance(zone, (str, Path)):
        gdf = gpd.read_file(zone)
        vor_crs = getattr(vor, "crs", None)
        if vor_crs is not None and gdf.crs is not None:
            gdf = gdf.to_crs(vor_crs)
        geom = gdf.union_all()
    else:
        geom = zone
    return centroid_zones(vor.gdf_vorPolys, [geom]) >= 0


#: Surface kinds whose content is per-cell arrays, so they cannot ride in a spec
#: (ledger 193 refused inlining arrays; ledger 269 names them here).
_UNSERIALIZABLE_SURFACES = frozenset({"array", "points"})


@dataclass(frozen=True)
class Surface:
    """One elevation surface, sampled onto a grid's cells when needed.

    Build one with a constructor rather than the raw fields: :meth:`raster`,
    :meth:`flat`, :meth:`offset_below`, :meth:`constant_thickness`,
    :meth:`from_contours`, or :meth:`from_points`.
    """

    kind: str
    #: Token :func:`myflopy.specs._json_value` writes this object under, and
    #: :func:`myflopy.specs._spec_value` revives it from. A `$`-prefixed key rather
    #: than a ``kind`` tag because ``Surface`` HAS a ``kind`` field of its own
    #: ("flat", "raster", "shift", ...) and the source-spec dispatcher keys on
    #: ``kind`` -- tagging one as the other makes a surface come back a plain dict.
    _SPEC_TOKEN = "$surface"

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-ready payload, or raise for a surface that holds arrays.

        ``from_array`` and ``from_points`` are REFUSED by name. Their content is
        per-cell data, and inlining arrays into a spec is what ledger 193 refused --
        ~3.3 MB per model at 9,405 cells, and it floods `Project.validate`. They stay
        usable; they just cannot be written down, which is what a run's frozen
        snapshot is for (ledger 269).
        """

        if self.kind in _UNSERIALIZABLE_SURFACES:
            raise ValueError(
                f"A '{self.kind}' surface holds per-cell arrays and cannot be written "
                "into a spec. Declare the layer from a raster, contours or a constant "
                "so the stack can be rebuilt on another mesh, or keep the built arrays "
                "and let the run's frozen snapshot carry them."
            )
        payload: dict[str, Any] = {}
        for item in fields(self):
            value = getattr(self, item.name)
            # `default_factory` fields have `default is MISSING`, so comparing against
            # `item.default` alone wrote an empty `grass_kwargs: {}` into every payload.
            default = (
                item.default_factory()
                if item.default is MISSING and item.default_factory is not MISSING
                else item.default
            )
            if item.name != "kind" and type(value) is type(default) and value == default:
                continue
            if item.name == "operands":
                payload[item.name] = [operand.to_dict() for operand in value]
            elif isinstance(value, Path):
                payload[item.name] = value.as_posix()
            else:
                payload[item.name] = value
        return payload

    def digest_payload(self) -> dict[str, Any]:
        """The declaration with file CONTENT standing in for file LOCATION.

        Same rule as :meth:`myflopy.sources.BoundarySource.digest_payload`, and it
        exists for the same measured reason one carrier over (ledger 293). A surface
        reaches ``PackageSpec.inputs`` through a lazy ``mf.disv`` or ``mf.evt``, where
        the file it names is read at BUILD time -- so digesting the path meant
        ``save(bundle_sources=True)`` moved the package's digest on a byte-identical
        archive, and a different raster at the same path digested identically.

        A path is replaced only where :attr:`source_shas` has a hash for it, so a
        surface whose file could not be read still digests by path rather than
        colliding with every other unreadable one.

        ``out`` is dropped outright rather than hashed. It is the GRASS interpolation
        CACHE -- a file this surface WRITES -- and its content is fully determined by
        ``contours`` plus the build parameters, both of which are digested here. So it
        carries nothing the digest needs, while bundling repoints it like any other
        path and would move the digest for a file nobody edited.
        """

        payload = self.to_dict()
        payload.pop("out", None)
        for name in self.source_shas:
            payload.pop(name, None)
        return payload

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Surface:
        """Rebuild a surface from :meth:`to_dict`'s payload."""

        known = {item.name for item in fields(cls)}
        unknown = set(data) - known
        if unknown:
            raise ValueError(f"Unknown Surface fields: {sorted(unknown)}")
        values = dict(data)
        for name in ("path", "contours", "region_raster", "region_vector", "out"):
            if values.get(name) is not None:
                values[name] = Path(values[name])
        if "operands" in values:
            values["operands"] = tuple(cls.from_dict(item) for item in values["operands"])
        return cls(**values)

    value: float | None = None
    path: Path | None = None
    fill: str | None = None
    #: Overrides the raster's declared no-data value; see `Surface.raster`.
    nodata: float | None = None
    units: str | None = None  # source units, if different from the model units
    # contour-interpolation parameters
    contours: Path | None = None
    z: str = "Elev"
    region_raster: Path | None = None
    region_vector: Path | None = None
    clip: bool = True
    resolution: float = 4
    out: Path | None = None
    epsg: str = str(DEFAULT_EPSG)
    grass_kwargs: dict = field(default_factory=dict)
    #: Show GRASS's own per-step progress while interpolating. Deliberately NOT
    #: part of `grass_kwargs`: those feed the derived-raster cache signature, so
    #: routing a display flag through them would invalidate the cache and force a
    #: silent re-interpolation the first time anyone asked to watch.
    progress: bool = False
    # point-interpolation parameters
    points: tuple[Any, Any, Any] | None = None
    method: str = "linear"
    # direct per-cell elevation array (e.g. an existing model's top/botm)
    array: Any = None
    # surface-algebra / isopach parameters
    operands: tuple = ()
    zone: Any = None
    #: ``{field name: 16-char content hash}`` for the files this surface reads, filled
    #: in when the surface becomes part of a LAZY package (``mf.disv`` on a gridless
    #: stack, ``mf.evt(surface=...)``). Empty everywhere else: a bare
    #: ``Surface.raster(path)`` deliberately opens nothing, and hashing a multi-hundred-
    #: megabyte DEM at construction would make the cheapest call in the layer API the
    #: most expensive. It is the DIGEST's stand-in for the path (:meth:`digest_payload`)
    #: -- see ledger 293's second pass.
    source_shas: dict = field(default_factory=dict)

    @classmethod
    def raster(
        cls,
        path: Path | str,
        *,
        fill: str | None = None,
        units: str | None = None,
        nodata: float | None = None,
    ) -> Surface:
        """A surface sampled from an existing raster.

        Parameters
        ----------
        path : Path or str
            The raster to sample.
        fill : {'propagate'}, optional
            ``"propagate"`` carries the surface above down into cells where the
            raster has no data (nodata or outside coverage), so the layer
            pinches against the surface above instead of leaving a ``NaN``
            elevation.
        units : str, optional
            e.g. ``"meters"``. Converts the sampled elevations to the model
            length units when sampled through a stack that declares them.
        nodata : float, optional
            Treat this value as no-data, overriding whatever the file declares.

            For rasters clipped to a boundary but written WITHOUT a no-data
            value in the header -- a very common GIS export -- where the fill
            (usually ``0``) is then read as a real elevation. Symptoms: a
            surface that is rectangular rather than domain-shaped, ground at
            elevation zero outside the boundary, and every layer beneath it
            collapsing to ``min_sep`` once the contacts get capped to it.

            Masking happens BEFORE the per-cell average, so a cell straddling
            the edge is the mean of its valid pixels rather than a blend with
            the sentinel.

            ``mf.Raster(dem, nodata=0)`` is the fix without rewriting the file;
            declaring it properly in the raster (``gdalwarp -dstnodata``) is the
            fix that helps every other tool too.

        Returns
        -------
        Surface

        Examples
        --------
        >>> mf.Raster("ground.tif")                       # header nodata is honoured
        >>> mf.Raster("clipped.tif", nodata=0)            # 0 is fill, not ground
        >>> mf.Raster("bedrock.tif", fill="propagate")    # gaps inherit from above
        """

        return cls(
            kind="raster", path=Path(path), fill=fill, units=units, nodata=nodata
        )

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
        region_vector: Path | str | None = None,
        resolution: float = 4,
        out: Path | str | None = None,
        epsg: str = str(DEFAULT_EPSG),
        units: str | None = None,
        clip: bool = True,
        progress: bool = False,
        **grass_kwargs: Any,
    ) -> Surface:
        """A surface interpolated from drawn elevation contours via GRASS.

        The interpolated raster is written to ``out`` (default: alongside the
        contours with a ``.interp.tif`` suffix) the first time the surface is
        resolved, then reused -- interpolate once, sample many times.

        **When does the interpolation run?** Not here. Declaring the surface is
        free; the first call that needs its values -- or merely its extent, which
        includes ``LayerStack.draft_grid`` and so a gridless ``stack.plot`` --
        resolves it. After that the GeoTIFF is reused across builds, sessions and
        processes. :meth:`cache_status` reports ``fresh``/stale/absent, and
        ``refresh=True`` on ``build``/``qc``/``to_disv`` forces a rebuild.

        Parameters
        ----------
        contours : Path or str
            Vector file of labelled elevation contour lines. Stored as a ``Path``.
        z : str, default "Elev"
            The attribute holding each contour's elevation.
        region_raster : Path or str, optional
            Raster whose extent becomes the GRASS region. **Exactly one of this
            or ``region_vector`` is required** -- without a region GRASS has no
            extent to work in. Stored as a ``Path``.

            Note it uses the WHOLE raster. A DEM covering more ground than your
            model makes every interpolation proportionally slower, and the extra
            area is the worst case for ``r.surf.contour`` because it has no
            contours to work from out there. Prefer ``region_vector``.
        region_vector : Path or str, optional
            Polygon whose extent becomes the GRASS region -- usually your model
            domain, which is normally what you want interpolated and no more.
            Stored as a ``Path``: a ``str`` was kept as given while
            :meth:`from_dict` revives a ``Path``, so a saved-and-reloaded contour
            surface compared ``!=`` the original with an identical payload (found
            2026-09-25).
        resolution : float, default 4
            Output cell size, in the CRS's units.
        out : Path or str, optional
            Where the interpolated GeoTIFF lands. Defaults to the contours path
            with a ``.interp.tif`` suffix. Stored as a ``Path``.
        epsg : str, default "2927"
            EPSG code for the GRASS location.
        units : str, optional
            Length unit of the elevations, converted to the model's if they differ.
        clip : bool, default True
            Null out everything beyond the region. With ``region_vector`` that
            means the POLYGON, not merely its bounding box -- a domain-shaped
            result, with the bbox corners left as nodata. Already what you want
            in nearly every case; ``False`` keeps the full rectangle, which is
            occasionally useful when a neighbouring model needs the overlap.
        progress : bool, default False
            Print GRASS's own per-step percentages while it works.

            GRASS modules are subprocesses that write to the inherited file
            descriptors, so **in a terminal you already see them** and this flag
            changes nothing. In a **Jupyter** notebook they go to the kernel's
            console instead of the cell, and `contextlib.redirect_stdout` cannot
            catch them because the writing happens below Python. This captures
            them at the fd level and re-emits them through Python, so they land
            in the cell.

            Off by default: library code logs rather than prints, and the flag is
            the human asking. It is NOT part of the cache signature, so turning it
            on does not re-interpolate anything.
        **grass_kwargs
            Forwarded to :class:`~myflopy.modflow.utils.contour_interp
            .ContourSurfaceInterpolator`, whose remaining keywords are
            ``grass_bin`` (the GRASS launcher), ``grassdata`` (the GRASS database
            directory, default ``~/grassdata``) and ``location`` (default
            ``"myflopy_interp"``). These DO feed the cache signature, so changing
            one invalidates the cached raster. A name it does not take is accepted
            HERE and refused only when the surface is first resolved, with a
            ``TypeError`` from that class (measured with ``extra=(1, 2)``; found
            2026-09-25, not changed).

        Returns
        -------
        Surface
            A ``kind="contours"`` surface. Nothing is read or interpolated yet; see
            "When does the interpolation run?" above.

        Raises
        ------
        ValueError
            Neither or both of ``region_raster`` and ``region_vector``.

        Examples
        --------
        >>> import myflopy as mf
        >>> ground = mf.Contours("ground.gpkg", z="Elev", region_vector="domain.gpkg")
        >>> ground.kind, ground.region_vector.as_posix(), ground.out.as_posix()
        ('contours', 'domain.gpkg', 'ground.interp.tif')
        >>> ground == mf.Surface.from_dict(ground.to_dict())   # what a reload builds
        True
        >>> watch = mf.Contours("clay.gpkg", region_vector="domain.gpkg", progress=True)
        """

        # Checked here, not in `_set_region`: that runs after GRASS has started a
        # session, so the message arrives late and buried. This is a typo you can
        # catch while writing the line.
        if (region_raster is None) == (region_vector is None):
            both = region_raster is not None
            raise ValueError(
                "Contours need exactly one region, and got "
                + ("both" if both else "neither")
                + ". Pass region_vector=<your domain polygon> (usually what you "
                "want) or region_raster=<a raster whose extent to use>. GRASS has "
                "no extent to interpolate in without one."
            )

        contours = Path(contours)
        out = Path(out) if out is not None else contours.with_suffix(".interp.tif")
        return cls(
            kind="contours",
            contours=contours,
            z=z,
            # `Path`, as `contours` and `out` are and as `from_dict` revives them: a
            # `str` kept as given made a round-tripped surface compare unequal.
            region_raster=Path(region_raster) if region_raster is not None else None,
            region_vector=Path(region_vector) if region_vector is not None else None,
            clip=clip,
            resolution=resolution,
            out=out,
            epsg=epsg,
            units=units,
            grass_kwargs=grass_kwargs,
            progress=progress,
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

        ``zone`` is a shapely geometry or a path to a vector file. A centroid
        exactly ON the zone's edge counts as inside it (ledger 243).
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

    @classmethod
    def toward(cls, target, fraction: float) -> Surface:
        """A surface ``fraction`` of the way from the one above down to ``target``.

        ``previous + fraction * (target - previous)``, per cell. This is the one
        primitive vertical subdivision needs: ``fraction=1/3`` of the way from a
        unit's top to its base is a contact that no amount of ``shift``/envelope
        composition can express, because the interval is per-cell and unknown
        until the surfaces are sampled.

        It exists INSTEAD of general ``Surface - Surface`` and ``Surface *
        scalar``. Those are deliberately absent (see
        ``test_subtracting_a_surface_from_a_surface_is_unsupported``) and adding
        them would widen the elevation-vs-thickness ambiguity that already makes
        a thickness map easy to confuse with a contact. The recurrence form also
        dodges the fact that :meth:`values` carries only ONE ``previous`` slot:
        each cut resolves against the contact immediately above it, so a caller
        subdividing a unit must renormalise cumulative fractions into per-step
        ones. :meth:`myflopy.layers.LayerStack.add`'s ``split=`` does that for
        you, and is the spelling to reach for.

        ``target`` must be ABSOLUTE. A relative target (``offset_below`` /
        ``constant_thickness``) re-resolves against each sub-contact rather than
        the unit top, which yields plausible, wrong, monotonically drifting
        elevations -- ``toward(constant_thickness(60), f)`` at thirds gives
        ``[80, 50, -10]`` where the answer is ``[80, 60, 40]``.

        Raises
        ------
        TypeError
            If ``target`` is a relative surface, or ``fraction`` is outside
            ``(0, 1]``.
        """

        target = cls._coerce(target)
        if _depends_on_previous(target):
            raise TypeError(
                f"Surface.toward() needs a target that does not depend on the "
                f"surface above it; {target.kind!r} does. Such a target measures "
                f"from each cut rather than from the unit top, so the contacts "
                f"drift downward instead of dividing the unit. Give the unit's "
                f"actual base surface -- e.g. Raster('sand_base.tif') -- or, for "
                f"a unit declared by thickness, split it with "
                f"`.add(name, thickness=T, split=N)`, which needs no target."
            )
        fraction = float(fraction)
        if not 0.0 < fraction <= 1.0:
            raise TypeError(
                f"Surface.toward() fraction must be in (0, 1], got {fraction}. "
                f"It is the share of the REMAINING interval to descend, so 1.0 "
                f"lands exactly on the target and 0 would be a zero-thickness "
                f"layer."
            )
        return cls(kind="toward", value=fraction, operands=(target,))

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
                region_vector=self.region_vector,
                clip=self.clip,
                resolution=self.resolution,
                epsg=self.epsg,
                progress=self.progress,
                **self.grass_kwargs,
            )

        sources = [self.contours]
        # BOTH region kinds belong here, not just the raster. `sources` is
        # content-tracked while `params` holds only values, so while
        # `region_vector` rode `grass_kwargs` its path was recorded and its
        # CONTENT was not -- editing the domain polygon left the cache "fresh"
        # and returned a raster interpolated over the old extent.
        for region in (self.region_raster, self.region_vector):
            if region is not None:
                sources.append(region)
        params = {
            "z": self.z,
            "resolution": self.resolution,
            "epsg": self.epsg,
            # Geometry-affecting, so it belongs here. `progress` deliberately
            # does NOT: it changes what you see, not what gets written.
            "clip": self.clip,
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

        if self.kind == "toward":
            # Handled here rather than via `is_relative`, matching `isopach`
            # above: both NEED `previous`, but neither is a fixed offset from
            # it, so `_RELATIVE_KINDS` (which `values` reads as
            # `previous - value`) would be exactly wrong for them.
            if previous is None:
                raise ValueError(
                    "Surface.toward cannot be the model top -- there is no "
                    "surface above to descend from."
                )
            return previous + float(self.value) * (_op(self.operands[0]) - previous)

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
            vor, [self.resolve_source(refresh=refresh)], [0], method=method,
            nodata=self.nodata,
            # A raster covering part of the model is normal here: `fill="propagate"`
            # inherits the rest from the surface above, and the pinch-out rules turn
            # what is left into idomain. Warning per cell would be noise.
            gaps="ignore",
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


def sample_bed(surface, grid, ground, *, length_units: str | None = None) -> np.ndarray:
    """Sample a lakebed or excavation floor onto every cell of ``grid``.

    The one bed sampler behind ``mf.lak(lake_bottom=)`` and
    ``LayerStack.excavate(bottom=)``, so a stack cut to a bed and a lake declared on
    it read the SAME float at every cell (the excavation design, section 3.5) --
    which is what lets ``mf.lak`` compare its bed with an excavation's record
    exactly. It is always AREA-weighted, whatever a stack's ``build(method=)`` says:
    that is how the lake samples, and a centroid sample of the same raster was
    measured up to 1.027 ft away from it.

    Parameters
    ----------
    surface : Surface, str or pathlib.Path
        The bed: any :class:`Surface` (``mf.Raster``, ``mf.Contours``, ``mf.Flat``,
        algebra), or a raster file path, read as ``Surface.raster(path)``.
    grid : VoronoiGridPlus
        The grid to sample onto.
    ground : array_like of float, shape (ncpl,)
        The model top -- the surface ABOVE the bed. A surface that consults the one
        above it (``offset_below``, ``fill="propagate"``) is resolved against it; an
        absolute one ignores it.
    length_units : str, optional
        The model's length unit. A surface declaring ``units=`` is converted to it
        (``"meters"`` -> ``"feet"`` multiplies by 3.2808...). ``None`` leaves the
        values as sampled.

    Returns
    -------
    numpy.ndarray of float, shape (ncpl,)
        One elevation per cell; ``NaN`` where the source has no data.

    Raises
    ------
    ValueError
        ``units=`` on a relative surface (``offset_below`` / ``constant_thickness``),
        which a unit factor cannot convert; an unknown length unit; any refusal
        :meth:`Surface.values` makes (a ``from_array`` of the wrong length, ...).

    Examples
    --------
    >>> import numpy as np
    >>> from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
    >>> from myflopy.surfaces import Surface, sample_bed
    >>> vor = irregular_voronoi_grid(nrow=2, ncol=2, cell_size=100.0)
    >>> sample_bed(Surface.flat(30.0, units="meters"), vor, np.full(vor.ncpl, 200.0),
    ...            length_units="feet").round(2).tolist()
    [98.43, 98.43, 98.43, 98.43]
    >>> sample_bed(Surface.offset_below(5.0), vor, np.full(vor.ncpl, 200.0)).tolist()
    [195.0, 195.0, 195.0, 195.0]
    """

    if isinstance(surface, (str, Path)):
        surface = Surface.raster(surface)
    previous = np.asarray(ground, dtype=float)
    values = np.asarray(surface.values(grid, previous=previous), dtype=float)
    if surface.units and length_units:
        if surface.is_relative:
            raise ValueError(
                "units= is only supported on absolute surfaces, not relative "
                "(offset_below / constant_thickness)."
            )
        values = values * _length_factor(surface.units, length_units)
    return values


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
            # Re-apply the ordering: `reconcile_surfaces` returns a plain frame
            # and `gpd.GeoDataFrame(df, geometry=...)` APPENDS geometry at the
            # end, undoing the `[["geometry"] + labels]` three lines above. That
            # made a keyword about physics move the column every positional
            # reader depends on -- see `geometry_first`.
            gdf = geometry_first(
                gpd.GeoDataFrame(df, geometry=gdf.geometry.values, crs=gdf.crs)
            )
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
        """Guard the fallback separation against each *pinching* layer's threshold.

        A ``"floor"`` layer has reconcile enforce its own ``min_thickness``, so it is
        exempt -- that is the whole mechanism by which ``min_thickness`` sets geometry.

        A ``"passthrough"`` or ``"inactive"`` layer is different: it must be allowed to
        come out thinner than its threshold, because that is the signal it pinches out.
        Reconcile gives it only the fallback ``min_sep``, so if the fallback is not
        below the threshold nothing can ever be thin enough to pinch and the policy is
        a silent no-op.
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
                    f"layer {i}: the fallback min_sep ({min_sep}) must be below its "
                    f"min_thickness ({min_thk[i]}) for pinch={policy[i]!r} to ever "
                    "fire -- reconcile would otherwise inflate every thin cell past "
                    "the threshold. Either lower min_sep, raise this layer's "
                    "min_thickness, or use pinch='floor' if you want the layer kept "
                    "at min_thickness rather than pinched out."
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
