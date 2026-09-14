"""Custom geometry drawn over a map -- wells, pond footprints, a boundary.

``map(features=...)`` takes a path to any vector file GeoPandas reads, a
``GeoDataFrame``/``GeoSeries``, a shapely geometry, a :class:`Feature` carrying
its own styling, or a list or dict of those. This module is the one place that
turns any of those spellings into *resolved layers* -- a WGS84 frame plus a
style -- and renders them, as Plotly ``Scattermap`` traces for the interactive
map and onto an ``Axes`` for the static one. One resolver, two renderers, so the
two backends of one map cannot disagree about what a feature layer is.
"""

from __future__ import annotations

import html
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import shapely as shp

from myflopy._logging import get_logger
from myflopy.modflow.utils.datatypes.hover import HoverStyle, format_number
from myflopy.viz import category_colors, mpl_color

logger = get_logger(__name__)

__all__ = [
    "Feature",
    "ResolvedFeature",
    "draw_features_mpl",
    "feature_traces",
    "resolve_features",
]

#: Attribute columns tried, in order, when a :class:`Feature` names no ``label``.
#: ``ExploName`` is the field the project's own location files carry.
DEFAULT_LABEL_COLUMNS = ("name", "Name", "NAME", "label", "Label", "id", "ID", "ExploName")

#: How many attribute columns a hover shows before it is truncated. A wells
#: table can carry forty columns; a tooltip cannot.
MAX_HOVER_COLUMNS = 12

@dataclass(frozen=True, eq=False)
class Feature:
    """One layer of custom geometry to draw over a map, with its styling.

    Compared by identity (``eq=False``): the source may be a ``GeoDataFrame``,
    whose ``==`` is elementwise and whose truth value raises.

    Pass it as ``map(features=Feature(...))`` or as an item of the list or dict
    handed to ``features=``. A bare path, GeoDataFrame or geometry is shorthand
    for ``Feature(source)`` with every default below.

    Parameters
    ----------
    source
        A path to a vector file (GeoPackage, shapefile, GeoJSON -- anything
        GeoPandas reads), a ``GeoDataFrame`` or ``GeoSeries``, a shapely
        geometry, or a list of shapely geometries.
    name : str, optional
        Legend label. Defaults to the file stem (plus the ``layer`` name when
        one is given), or ``"features"`` for in-memory geometry.
    layer : str or int, optional
        Which layer to read from a multi-layer file such as a GeoPackage.
        ``None`` reads the first. Refused for an in-memory source.
    color : str, optional
        A color both backends accept: a hex string, ``rgb(...)``/``rgba(...)``,
        or a CSS color name. Validated at construction. Defaults to the house
        categorical palette keyed by ``name``, so a layer keeps its color
        across figures.
    width : float, default 2.0
        Line width for lines and polygon outlines, in pixels.
    size : float, default 9.0
        Marker diameter for points, in pixels.
    fill : bool, default False
        Fill polygons (at ``opacity``) instead of outlining them.
    opacity : float, default 0.35
        Fill opacity for filled polygons, 0 to 1.
    label : str, optional
        Attribute column whose value heads each feature's hover. Defaults to
        the first of ``name``/``label``/``id``/``ExploName`` present, else the
        layer name.
    hover : bool or str or sequence of str, default True
        Which attribute columns the hover lists under the label: ``True`` for
        all of them (at most 12), ``False`` for none, or the column names.
    crs
        Declares the coordinate system of a source that carries none. A source
        that already has a CRS keeps it; a different ``crs=`` then raises.
        Defaults to the grid's CRS.
    show_legend : bool, default True
        Give the layer a legend entry.
    """

    source: Any
    name: str | None = None
    layer: str | int | None = None
    color: str | None = None
    width: float = 2.0
    size: float = 9.0
    fill: bool = False
    opacity: float = 0.35
    label: str | None = None
    hover: bool | str | Sequence[str] = True
    crs: Any = None
    show_legend: bool = True

    def __post_init__(self):
        if self.source is None:
            raise ValueError("Feature needs a source: a vector file path, a GeoDataFrame or a geometry.")
        if self.layer is not None and not isinstance(self.source, (str, Path)):
            raise ValueError(
                "layer= selects a layer inside a vector FILE, but the source is a "
                f"{type(self.source).__name__}. Drop layer=, or pass the file path."
            )
        for attribute in ("width", "size"):
            value = getattr(self, attribute)
            if not np.isfinite(value) or value < 0:
                raise ValueError(f"Feature {attribute} must be zero or greater, not {value!r}.")
        if not 0.0 <= float(self.opacity) <= 1.0:
            raise ValueError(f"Feature opacity must be between 0 and 1, not {self.opacity!r}.")
        if self.color is not None:
            try:
                mpl_color(self.color)
            except ValueError as error:
                raise ValueError(
                    f"Feature color {self.color!r} is not a color both backends accept; "
                    "use a hex string, rgb(...)/rgba(...), or a CSS color name."
                ) from error
        if not isinstance(self.hover, (bool, str)) and not (
            isinstance(self.hover, Sequence) and all(isinstance(c, str) for c in self.hover)
        ):
            raise TypeError("Feature hover must be True, False, a column name, or a sequence of column names.")


@dataclass(frozen=True, eq=False)
class ResolvedFeature:
    """A :class:`Feature` read, reprojected to WGS84 and given its final look."""

    name: str
    frame: gpd.GeoDataFrame
    color: str
    spec: Feature
    hover_columns: tuple[str, ...]
    label_column: str | None

    @property
    def show_legend(self) -> bool:
        """Whether this layer earns a legend entry."""

        return bool(self.spec.show_legend)


# ---------------------------------------------------------------------------
# resolving
# ---------------------------------------------------------------------------


def _is_geometry(value) -> bool:
    return isinstance(value, shp.geometry.base.BaseGeometry)


def _as_spec(item, name: str | None = None) -> Feature:
    """Coerce one ``features=`` item to a :class:`Feature`.

    A dict key (``name``) wins over the spec's own name: it is the outer, explicit
    label the caller wrote next to the layer.
    """

    if isinstance(item, Feature):
        return replace(item, name=name) if name is not None else item
    if isinstance(item, (str, Path, gpd.GeoDataFrame, gpd.GeoSeries)) or _is_geometry(item):
        return Feature(item, name=name)
    if isinstance(item, Sequence) and not isinstance(item, str) and item and all(_is_geometry(g) for g in item):
        return Feature(list(item), name=name)
    raise TypeError(
        "features= accepts a vector file path, a GeoDataFrame, a shapely geometry, "
        f"a Feature, or a list/dict of those; got {type(item).__name__}."
    )


def _specs(features) -> list[Feature]:
    """Flatten any accepted ``features=`` spelling into a list of specs."""

    if features is None:
        return []
    if isinstance(features, Mapping):
        return [_as_spec(item, name=str(key)) for key, item in features.items()]
    if isinstance(features, (list, tuple)):
        if features and all(_is_geometry(g) for g in features):
            return [Feature(list(features))]
        return [_as_spec(item) for item in features]
    return [_as_spec(features)]


def combine_features(features, extra):
    """``features`` with one more layer appended, whatever spelling it uses.

    Used by the retired location-marker option, which is folded into whatever
    the caller also passed as ``features=`` rather than replacing it.
    """

    if features is None:
        return extra
    if isinstance(features, Mapping):
        key = getattr(extra, "name", None) or "locations"
        candidate, n = key, 1
        while candidate in features:
            n += 1
            candidate = f"{key} {n}"
        return {**features, candidate: extra}
    if isinstance(features, (list, tuple)) and not (
        features and all(_is_geometry(g) for g in features)
    ):
        return [*features, extra]
    return [features, extra]


def _read_vector(path: Path, layer) -> gpd.GeoDataFrame:
    """Read one vector file, naming which layer of a multi-layer file to use.

    ``layer=None`` is GDAL's default -- the first layer. A wrong name or index
    raises listing the layers actually present, because GDAL's own error names
    none of them. Written here rather than imported from
    ``mf6.grid.selection``: that package's ``__init__`` imports the choropleth,
    which imports this module, so the shared reader is one import away from a
    cycle that only shows in a fresh interpreter.
    """

    if layer is None:
        return gpd.read_file(path)
    lister = getattr(gpd, "list_layers", None)
    available = list(lister(path)["name"]) if lister is not None else None
    if available is not None:
        if isinstance(layer, int):
            if not 0 <= layer < len(available):
                raise ValueError(
                    f"{path.name} has {len(available)} layer(s), so index {layer} is out "
                    f"of range; layers: {available}"
                )
        elif layer not in available:
            raise ValueError(f"{path.name} has no layer {layer!r}; layers: {available}")
    return gpd.read_file(path, layer=layer)


def _read_source(spec: Feature) -> gpd.GeoDataFrame:
    source = spec.source
    if isinstance(source, (str, Path)):
        path = Path(source)
        if not path.exists():
            raise FileNotFoundError(f"features= names a vector file that does not exist: {path}")
        return _read_vector(path, spec.layer)
    if isinstance(source, gpd.GeoDataFrame):
        return source.copy()
    if isinstance(source, gpd.GeoSeries):
        return gpd.GeoDataFrame(geometry=source.copy())
    if _is_geometry(source):
        return gpd.GeoDataFrame(geometry=[source])
    return gpd.GeoDataFrame(geometry=list(source))


def _default_name(spec: Feature) -> str:
    if spec.name:
        return str(spec.name)
    if isinstance(spec.source, (str, Path)):
        stem = Path(spec.source).stem
        return f"{stem}:{spec.layer}" if spec.layer is not None else stem
    return "features"


def _unique(names: list[str]) -> list[str]:
    """Suffix duplicate legend names so two layers never share one entry.

    The suffixed candidate is re-checked: ``["a", "a", "a 2"]`` must not hand
    two layers the same ``"a 2"`` (one legend toggle would hide both).
    """

    used: set[str] = set()
    out = []
    for name in names:
        candidate, n = name, 1
        while candidate in used:
            n += 1
            candidate = f"{name} {n}"
        used.add(candidate)
        out.append(candidate)
    return out


def _to_wgs84(frame: gpd.GeoDataFrame, spec: Feature, grid_crs, name: str) -> gpd.GeoDataFrame:
    if frame.crs is None:
        crs = spec.crs if spec.crs is not None else grid_crs
        if crs is None:
            raise ValueError(
                f"features layer {name!r} has no CRS and the grid declares none; "
                "pass Feature(..., crs=...) so it can be placed on the map."
            )
        frame = frame.set_crs(crs)
    elif spec.crs is not None and not frame.crs.equals(gpd.GeoDataFrame(geometry=[], crs=spec.crs).crs):
        raise ValueError(
            f"features layer {name!r} already declares CRS {frame.crs.to_string()}; "
            f"drop crs= or reproject the source instead of overriding it with {spec.crs!r}."
        )
    return frame.to_crs(epsg=4326)


def _geometry_columns(frame: gpd.GeoDataFrame) -> set[str]:
    """Every geometry-typed column, active or not (``set_geometry`` leaves the old one)."""

    return set(frame.select_dtypes("geometry").columns) | {frame.geometry.name}


def _hover_columns(frame: gpd.GeoDataFrame, spec: Feature, label_column: str | None) -> tuple[str, ...]:
    attributes = [c for c in frame.columns if c not in _geometry_columns(frame)]
    if spec.hover is False:
        return ()
    if spec.hover is True:
        columns = [c for c in attributes if c != label_column]
        if len(columns) > MAX_HOVER_COLUMNS:
            logger.info(
                "features layer hover lists the first %d of %d attribute columns; "
                "pass Feature(hover=[...]) to choose them",
                MAX_HOVER_COLUMNS, len(columns),
            )
        return tuple(columns[:MAX_HOVER_COLUMNS])
    wanted = [spec.hover] if isinstance(spec.hover, str) else list(spec.hover)
    missing = [c for c in wanted if c not in attributes]
    if missing:
        raise KeyError(
            f"features hover names column(s) {missing} that the layer does not have; "
            f"available: {attributes}"
        )
    return tuple(wanted)


def resolve_features(features, *, crs=None) -> list[ResolvedFeature]:
    """Turn any accepted ``features=`` spelling into resolved, WGS84 layers.

    ``crs`` is the grid's coordinate system, assumed for a source that carries
    none. Raises -- rather than warning and drawing nothing -- for a missing
    file, an unknown layer, an unplaceable frame or an unknown hover column,
    because a map that silently omits the wells it was asked for looks like an
    answer. An EMPTY source is the one thing that only warns: there is nothing
    wrong with the request, just nothing in it.
    """

    specs = _specs(features)
    if not specs:
        return []
    names = _unique([_default_name(spec) for spec in specs])
    colors = category_colors([n for spec, n in zip(specs, names, strict=False) if spec.color is None])
    resolved = []
    for spec, name in zip(specs, names, strict=False):
        frame = _read_source(spec)
        # shapely's vectorised predicates, not GeoSeries.notna(): the latter
        # warns on every call for a frame holding one empty geometry.
        geometries = frame.geometry.to_numpy()
        frame = frame[~shp.is_missing(geometries) & ~shp.is_empty(geometries)]
        if frame.empty:
            logger.warning("features layer %r holds no geometry; nothing is drawn for it", name)
            continue
        frame = _to_wgs84(frame, spec, crs, name)
        label_column = spec.label
        if label_column is not None and label_column not in frame.columns:
            raise KeyError(
                f"features label column {label_column!r} is not in layer {name!r}; "
                f"available: {[c for c in frame.columns if c != frame.geometry.name]}"
            )
        if label_column is None:
            geometry_columns = _geometry_columns(frame)
            label_column = next(
                (c for c in DEFAULT_LABEL_COLUMNS if c in frame.columns and c not in geometry_columns),
                None,
            )
        resolved.append(
            ResolvedFeature(
                name=name,
                frame=frame,
                color=spec.color if spec.color is not None else colors[name],
                spec=spec,
                hover_columns=_hover_columns(frame, spec, label_column),
                label_column=label_column,
            )
        )
    return resolved


# ---------------------------------------------------------------------------
# rendering: plotly
# ---------------------------------------------------------------------------


def _is_missing(value) -> bool:
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def _hover_value(value) -> str:
    """One attribute value as hover text: numbers formatted, text escaped."""

    if _is_missing(value):
        return ""
    if isinstance(value, (bool, np.bool_)):
        return str(bool(value))
    if isinstance(value, (int, float, np.integer, np.floating)):
        return format_number(value)
    return html.escape(str(value))


def _hover_text(layer: ResolvedFeature, row: pd.Series) -> str:
    head = layer.name
    if layer.label_column is not None:
        value = row[layer.label_column]
        head = layer.name if _is_missing(value) else html.escape(str(value))
    lines = [f"<b>{head}</b>"]
    for column in layer.hover_columns:
        lines.append(f"{html.escape(str(column))}: {_hover_value(row[column])}")
    return "<br>".join(lines)


def _parts(geometry):
    """The single-part geometries inside ``geometry`` (itself, if single)."""

    if geometry.geom_type in ("MultiPoint", "MultiLineString", "MultiPolygon", "GeometryCollection"):
        parts = []
        for part in shp.get_parts(geometry):
            parts.extend(_parts(part))
        return parts
    return [geometry]


def _rings(polygon):
    yield polygon.exterior
    yield from polygon.interiors


def _split(layer: ResolvedFeature):
    """``(points, lines, polygons)``: per kind, ``(lon, lat, text)`` lists.

    Lines and rings are separated with ``None`` so one trace draws the whole
    layer; a per-row trace would give a 500-well file 500 legend entries.
    """

    points: dict[str, list] = {"lon": [], "lat": [], "text": []}
    lines: dict[str, list] = {"lon": [], "lat": [], "text": []}
    polygons: dict[str, list] = {"lon": [], "lat": [], "text": []}
    exteriors: dict[str, list] = {"lon": [], "lat": []}

    def _append(target, xs, ys, text):
        if target["lon"]:
            target["lon"].append(None)
            target["lat"].append(None)
            target["text"].append(None)
        target["lon"].extend(float(x) for x in xs)
        target["lat"].extend(float(y) for y in ys)
        target["text"].extend([text] * len(xs))

    # `row.geometry` on a plain row is the column LABELLED "geometry", which
    # after `set_geometry("buffer")` is a stale, un-reprojected column. Read the
    # active geometry alongside the row instead.
    for geometry, (_, row) in zip(layer.frame.geometry, layer.frame.iterrows(), strict=True):
        text = _hover_text(layer, row)
        for part in _parts(geometry):
            if part.is_empty:
                continue
            kind = part.geom_type
            if kind == "Point":
                points["lon"].append(float(part.x))
                points["lat"].append(float(part.y))
                points["text"].append(text)
            elif kind in ("LineString", "LinearRing"):
                xs, ys = part.xy
                _append(lines, xs, ys, text)
            elif kind == "Polygon":
                for ring in _rings(part):
                    xs, ys = ring.xy
                    _append(polygons, xs, ys, text)
                xs, ys = part.exterior.xy
                if exteriors["lon"]:
                    exteriors["lon"].append(None)
                    exteriors["lat"].append(None)
                exteriors["lon"].extend(float(x) for x in xs)
                exteriors["lat"].extend(float(y) for y in ys)
            else:  # pragma: no cover - shapely has no other single-part kinds
                logger.warning("features layer %r: skipping a %s geometry", layer.name, kind)
    return points, lines, polygons, exteriors


def _rgba(color: str, opacity: float) -> str:
    """``color`` as an ``rgba(...)`` string at ``opacity`` (hex or rgb in)."""

    r, g, b = (round(c * 255) for c in mpl_color(color)[:3])
    return f"rgba({r},{g},{b},{float(opacity):.3f})"


def feature_traces(resolved: Sequence[ResolvedFeature]) -> list[go.Scattermap]:
    """Plotly ``Scattermap`` traces for resolved layers, one legend entry each.

    A layer holding several geometry kinds becomes several traces sharing one
    ``legendgroup`` (and one legend entry), so toggling the entry hides all of
    them. Filled polygons add a line-less fill trace beneath their outline.
    """

    hoverlabel = HoverStyle().to_hoverlabel()
    traces = []
    for layer in resolved:
        spec = layer.spec
        points, lines, polygons, exteriors = _split(layer)
        common = {
            "name": layer.name,
            "legendgroup": layer.name,
            "hovertemplate": "%{text}<extra></extra>",
            "hoverlabel": hoverlabel,
        }
        first = True

        def _legend(first):
            return {"showlegend": bool(first and layer.show_legend)}

        if spec.fill and exteriors["lon"]:
            traces.append(
                go.Scattermap(
                    mode="lines",
                    lon=exteriors["lon"],
                    lat=exteriors["lat"],
                    fill="toself",
                    fillcolor=_rgba(layer.color, spec.opacity),
                    line={"width": 0, "color": layer.color},
                    name=layer.name,
                    legendgroup=layer.name,
                    hoverinfo="skip",
                    showlegend=False,
                )
            )
        if polygons["lon"]:
            traces.append(
                go.Scattermap(
                    mode="lines",
                    lon=polygons["lon"],
                    lat=polygons["lat"],
                    text=polygons["text"],
                    line={"color": layer.color, "width": spec.width},
                    **common,
                    **_legend(first),
                )
            )
            first = False
        if lines["lon"]:
            traces.append(
                go.Scattermap(
                    mode="lines",
                    lon=lines["lon"],
                    lat=lines["lat"],
                    text=lines["text"],
                    line={"color": layer.color, "width": spec.width},
                    **common,
                    **_legend(first),
                )
            )
            first = False
        if points["lon"]:
            traces.append(
                go.Scattermap(
                    mode="markers",
                    lon=points["lon"],
                    lat=points["lat"],
                    text=points["text"],
                    marker={"size": spec.size, "color": layer.color},
                    **common,
                    **_legend(first),
                )
            )
            first = False
    return traces


# ---------------------------------------------------------------------------
# rendering: matplotlib
# ---------------------------------------------------------------------------


def draw_features_mpl(resolved: Sequence[ResolvedFeature], ax, *, crs) -> None:
    """Draw resolved layers onto a matplotlib ``Axes`` in the grid's ``crs``.

    The static map is drawn in grid coordinates (``vor.gdf_vorPolys``), so the
    WGS84 frames are projected back onto it. Adds a legend when any layer asks
    for one.
    """

    from matplotlib.lines import Line2D

    handles = []
    for layer in resolved:
        spec = layer.spec
        frame = layer.frame.to_crs(crs) if crs is not None else layer.frame
        color = mpl_color(layer.color)
        kinds = frame.geometry.geom_type
        points = frame[kinds.isin(["Point", "MultiPoint"])]
        lines = frame[kinds.isin(["LineString", "MultiLineString", "LinearRing"])]
        polygons = frame[kinds.isin(["Polygon", "MultiPolygon"])]
        collections = frame[kinds == "GeometryCollection"]
        while not collections.empty:
            # Explode until no collection remains, so a collection inside a
            # collection is drawn -- the Plotly renderer recurses the same way.
            exploded = collections.explode(index_parts=False)
            kinds = exploded.geometry.geom_type
            points = pd.concat([points, exploded[kinds.isin(["Point", "MultiPoint"])]])
            lines = pd.concat([lines, exploded[kinds.isin(["LineString", "MultiLineString", "LinearRing"])]])
            polygons = pd.concat([polygons, exploded[kinds.isin(["Polygon", "MultiPolygon"])]])
            collections = exploded[kinds == "GeometryCollection"]
        if not polygons.empty:
            polygons.plot(
                ax=ax,
                facecolor=(*color[:3], spec.opacity) if spec.fill else "none",
                edgecolor=color,
                linewidth=spec.width,
                zorder=4,
            )
        if not lines.empty:
            lines.plot(ax=ax, color=color, linewidth=spec.width, zorder=4)
        if not points.empty:
            points.plot(ax=ax, color=color, markersize=spec.size ** 2 / 2.0, zorder=5)
        if layer.show_legend:
            if not points.empty and polygons.empty and lines.empty:
                handles.append(Line2D([], [], marker="o", linestyle="", color=color,
                                      markersize=spec.size * 0.7, label=layer.name))
            else:
                handles.append(Line2D([], [], color=color, linewidth=max(spec.width, 1.0),
                                      label=layer.name))
    if handles:
        ax.legend(handles=handles, loc="best", frameon=True)
