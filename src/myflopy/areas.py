"""Areas of interest: named places to look at, shared by a model and its project.

``model.aoi`` / ``project.aoi`` hold :class:`Area` objects by name, so
``model.plot.map(extent="site")`` -- and ``select="site"``, and
``features=model.aoi["site"]`` -- mean the same place on every map, mosaic and
array call (ledger 310).

**Defined on the model, kept by the project** (the user's design): a model that
came from a project -- ``prepare_run``, ``reopen_run``, ``project.runs[...]``,
``mf.load_run`` of a run inside one -- shares that project's LIVE library while the
``Project`` object exists, so ``model.add_aoi(...)`` IS ``project.add_aoi(...)``
and ``project.save()`` keeps it. With no live project, a model reads the project's
SAVED areas once, and anything added there lasts for the session only; a model in
no project at all keeps its areas in memory, like ``model.regions``.

This module is the meeting point, in import layer 0, because the model (layer 8)
cannot import the workspace (layer 10) that owns the ``Project``: a live project
REGISTERS its library here under its root, and a model looks it up by the root its
run recorded.

An area is a GEOMETRY plus the CRS it is in, independent of any grid -- a
``region`` (``model.add_region``) is a set of CELLS on one grid, a different thing,
and a name defined as both is refused where it is used.
"""

from __future__ import annotations

import json
import weakref
from dataclasses import dataclass, replace
from numbers import Real
from pathlib import Path
from typing import Any

import geopandas as gpd
import numpy as np
import shapely as shp
from shapely.geometry.base import BaseGeometry

from myflopy._logging import get_logger

logger = get_logger(__name__)

__all__ = ["Area", "live_areas", "register_project_areas", "saved_areas"]


def _crs_string(crs) -> str | None:
    """A CRS as the string a document stores, or None."""

    if crs is None:
        return None
    to_string = getattr(crs, "to_string", None)
    return str(to_string()) if callable(to_string) else str(crs)


@dataclass(frozen=True)
class Area:
    """One area of interest: a geometry and the CRS it is in.

    ``crs=None`` means "the grid's frame" -- the rule ``features=`` and ``extent=``
    use for numbers and bare geometries -- so a box typed in model coordinates
    lands where it was typed on any model, including one left at the origin on
    purpose. It is never filled in with a default.

    Build one with :meth:`coerce` (what ``add_aoi`` does): four numbers, a run of
    ``(x, y)`` points, any shapely geometry, a GeoDataFrame, or a path to a vector
    file. The WHOLE geometry is kept -- ``extent=`` uses its bounds, ``select=`` the
    shape itself.
    """

    geometry: BaseGeometry
    crs: str | None = None
    note: str | None = None
    name: str | None = None
    source: str | None = None

    @classmethod
    def coerce(cls, value, *, crs=None, layer=None, note=None, name=None) -> Area:
        """``value`` as an :class:`Area` -- see the class docstring for the spellings.

        A file or frame carrying its own CRS keeps it; a ``crs=`` that disagrees
        with it raises rather than silently relabelling the data.
        """

        if isinstance(value, Area):
            if crs is not None and _crs_string(crs) != value.crs:
                raise ValueError(
                    f"this area is already in {value.crs}; crs={crs!r} would relabel it. "
                    "Reproject its geometry instead."
                )
            updates = {k: v for k, v in (("note", note), ("name", name)) if v is not None}
            return replace(value, **updates)
        source = None
        if isinstance(value, (str, Path)):
            source = str(Path(value))
            value = gpd.read_file(value, layer=layer)
        if isinstance(value, (gpd.GeoSeries, gpd.GeoDataFrame)):
            carried = _crs_string(value.crs)
            if carried is not None and crs is not None:
                wanted = gpd.GeoSeries([], crs=crs).crs
                if not value.crs.equals(wanted):
                    raise ValueError(
                        f"the geometry already declares CRS {carried}; drop crs= or "
                        f"reproject it instead of relabelling it as {crs!r}."
                    )
            if len(value) == 0:
                raise ValueError("an area needs a geometry; this frame is empty.")
            geometry = shp.union_all(np.asarray(value.geometry))
            return cls(geometry, carried or _crs_string(crs), note, name, source)
        if hasattr(value, "geom_type"):
            return cls(value, _crs_string(crs), note, name)
        array = np.asarray(value, dtype=object) if not isinstance(value, dict) else None
        if array is not None and array.size and all(
            isinstance(item, Real) and not isinstance(item, bool) for item in array.reshape(-1)
        ):
            numbers = np.asarray(value, dtype=float)
            if numbers.shape == (4,):
                xmin, ymin, xmax, ymax = numbers.tolist()
                if not (xmin < xmax and ymin < ymax):
                    raise ValueError(
                        f"{tuple(numbers.tolist())} is not (xmin, ymin, xmax, ymax) with "
                        "xmin < xmax and ymin < ymax -- shapely's order, not FloPy's."
                    )
                return cls(shp.box(xmin, ymin, xmax, ymax), _crs_string(crs), note, name)
            if numbers.ndim == 2 and numbers.shape[1] == 2:
                box = shp.box(*shp.MultiPoint(numbers).bounds)
                return cls(box, _crs_string(crs), note, name)
        raise TypeError(
            f"{value!r} is not an area: pass (xmin, ymin, xmax, ymax), a list of (x, y) "
            "points, a shapely geometry, a GeoDataFrame, or a path to a vector file."
        )

    def to_geoseries(self, crs=None) -> gpd.GeoSeries:
        """This area as a one-row GeoSeries, reprojected to ``crs`` when it has one."""

        series = gpd.GeoSeries([self.geometry], crs=self.crs)
        if crs is not None and self.crs is not None:
            series = series.to_crs(crs)
        return series

    def bounds(self, crs=None) -> tuple[float, float, float, float]:
        """``(xmin, ymin, xmax, ymax)``, in ``crs`` when both it and the area have one."""

        return tuple(float(v) for v in self.to_geoseries(crs).total_bounds)

    def to_dict(self) -> dict[str, Any]:
        """A JSON-safe document: WKT plus CRS -- plain Python types only.

        Inline rather than a path to a file, like an eager package: the geometry is
        read once, when the area is defined, so a saved project is complete without
        the file it came from (``source`` stays as provenance).
        """

        # Full precision: shapely rounds WKT to 6 decimals by default, which moves
        # a lon/lat vertex ~0.1 m and a State Plane one not at all -- a silent
        # CRS-dependent loss is not worth the bytes.
        wkt = shp.to_wkt(self.geometry, rounding_precision=-1)
        doc: dict[str, Any] = {"kind": "Area", "wkt": wkt, "crs": self.crs}
        for key in ("note", "source"):
            value = getattr(self, key)
            if value is not None:
                doc[key] = str(value)
        return doc

    @classmethod
    def from_dict(cls, doc: dict[str, Any], *, name: str | None = None) -> Area:
        """The inverse of :meth:`to_dict`."""

        return cls(
            shp.from_wkt(doc["wkt"]),
            doc.get("crs"),
            doc.get("note"),
            name,
            doc.get("source"),
        )

    def __repr__(self) -> str:
        xmin, ymin, xmax, ymax = self.geometry.bounds
        label = f"{self.name!r}, " if self.name else ""
        return (
            f"Area({label}{self.geometry.geom_type}, "
            f"bounds=({xmin:.6g}, {ymin:.6g}, {xmax:.6g}, {ymax:.6g}), crs={self.crs!r})"
        )


# -- the meeting point between a model and a live project ---------------------------

#: root -> the LIVE project's area library. Weak: a project that is gone must not be
#: kept alive by a model that once looked at it.
_LIVE: weakref.WeakValueDictionary[str, Any] = weakref.WeakValueDictionary()


def _key(root) -> str:
    return str(Path(root).expanduser().resolve())


def register_project_areas(root, library) -> None:
    """Record ``library`` as the live area library of the project at ``root``."""

    _LIVE[_key(root)] = library


def live_areas(root):
    """The live project's area library at ``root``, or None."""

    return None if root is None else _LIVE.get(_key(root))


def saved_areas(root) -> dict[str, Area] | None:
    """The areas a SAVED project at ``root`` holds, or None when there is none there.

    Read straight from the project's JSON (``project.json`` names the spec file,
    which carries an ``aoi`` block), so nothing else in the project is loaded --
    ``Project.load`` would unpickle grids and packages to answer a question about
    boxes.
    """

    if root is None:
        return None
    manifest = Path(root) / "project.json"
    if not manifest.is_file():
        return None
    try:
        spec = json.loads(manifest.read_text())["project"]["spec"]
        document = json.loads((Path(root) / spec).read_text())
    except (OSError, KeyError, TypeError, ValueError):
        # OSError: a spec file gone missing; KeyError/TypeError: a manifest of an
        # older shape; ValueError: unparseable JSON. Each means no saved areas to
        # offer -- the model keeps a session library instead.
        logger.debug("no saved areas readable under %s", root, exc_info=True)
        return None
    return {
        name: Area.from_dict(doc, name=name) for name, doc in (document.get("aoi") or {}).items()
    }
