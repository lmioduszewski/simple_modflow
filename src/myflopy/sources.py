"""Serializable source specifications for project-backed model recipes."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, ClassVar


def _path_text(path: Path | str) -> str:
    """Serialize a path as a forward-slash string for portable JSON/YAML output."""

    return str(Path(path).as_posix())


@dataclass(frozen=True, slots=True)
class DataSourceSpec:
    """Durable, serializable reference to a file that feeds a grid or package.

    Source specs are how a :class:`~myflopy.specs.SimulationSpec` records *where*
    its input data lives without inlining the data itself, so a recipe can be
    written to JSON/YAML, version-controlled, and resolved again later. Each
    subclass adds the fields a given file type needs (a table name, a CRS, a
    raster band, ...); this base captures the parts they all share -- the file
    ``path``, whether it is ``external`` to the project, and free-form
    ``metadata``. ``to_dict``/``from_dict`` round-trip a spec through its
    ``kind`` tag, dispatching to the right subclass on the way back.

    You rarely build these by hand: the project library helpers
    (``project.add_package``/``add_grid``) and the GIS-aware package forms
    (``mf.ghb.gpkg`` etc.) construct the appropriate subclass for you. Reach for
    them directly only when authoring a serialized recipe.

    Parameters
    ----------
    path
        Filesystem path to the source file. Stored as a :class:`~pathlib.Path`.
    external
        ``True`` marks the file as living outside the project tree (it is
        referenced in place rather than copied into the run workspace).
    metadata
        Arbitrary JSON-serializable annotations carried alongside the reference.

    See Also
    --------
    TableSource, ShapeSource, GeoPackageSourceSpec, RasterSource, LiteralSource
    """

    path: Path | str
    external: bool = False
    metadata: dict[str, Any] = field(default_factory=dict)

    kind: ClassVar[str] = "DataSourceSpec"

    def __post_init__(self) -> None:
        """Coerce ``path`` to a :class:`~pathlib.Path` (on this frozen dataclass)."""

        object.__setattr__(self, "path", Path(self.path))

    def to_dict(self) -> dict[str, Any]:
        """Serialize to a ``kind``-tagged dict, omitting default ``external``/``metadata``."""

        payload: dict[str, Any] = {
            "kind": self.kind,
            "path": _path_text(self.path),
        }
        if self.external:
            payload["external"] = True
        if self.metadata:
            payload["metadata"] = self.metadata
        return payload

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> DataSourceSpec:
        """Rebuild a source spec, dispatching on its ``kind`` tag to the right subclass."""

        kind = str(data.get("kind"))
        try:
            source_type = _SOURCE_TYPES[kind]
        except KeyError as error:
            raise ValueError(f"Unknown data source kind: {kind!r}") from error
        return source_type._from_dict(data)

    @classmethod
    def _from_dict(cls, data: dict[str, Any]) -> DataSourceSpec:
        """Construct this concrete class from an already-dispatched payload dict."""

        return cls(
            path=data["path"],
            external=bool(data.get("external", False)),
            metadata=dict(data.get("metadata", {})),
        )


@dataclass(frozen=True, slots=True)
class TableSource(DataSourceSpec):
    """Reference to a tabular file (CSV/Excel/Parquet) that supplies package records.

    Use this for stress-period data, well lists, observation tables, and similar
    row-oriented inputs. ``table`` names the sheet/table within a multi-table
    workbook (an Excel sheet name, for example) and may be omitted for a single
    flat file such as a CSV.

    Parameters
    ----------
    path
        Path to the ``.csv`` / ``.xlsx`` / ``.parquet`` file.
    table
        Optional sheet or table name inside a multi-table workbook.

    Examples
    --------
    >>> TableSource("wells.csv")
    >>> TableSource("inputs.xlsx", table="ghb_stress")
    """

    table: str | None = None
    kind: ClassVar[str] = "TableSource"

    def to_dict(self) -> dict[str, Any]:
        """Serialize like the base spec, adding ``table`` when set."""

        payload = DataSourceSpec.to_dict(self)
        if self.table is not None:
            payload["table"] = self.table
        return payload

    @classmethod
    def _from_dict(cls, data: dict[str, Any]) -> TableSource:
        """Rebuild a :class:`TableSource` from a payload dict."""

        return cls(
            path=data["path"],
            table=data.get("table"),
            external=bool(data.get("external", False)),
            metadata=dict(data.get("metadata", {})),
        )


@dataclass(frozen=True, slots=True)
class ShapeSource(DataSourceSpec):
    """Reference to a vector file (shapefile/GeoJSON) used for grid or BC geometry.

    Use this for single-layer vector files -- a boundary polygon for grid
    generation, refinement lines, or features that drive a list boundary
    condition. For a multi-layer GeoPackage, use :class:`GeoPackageSourceSpec`
    instead (it can name a layer and carry a field mapping).

    Pass everything after ``path`` **by keyword**: the inherited fields come
    first in positional order -- ``(path, external, metadata, crs)`` -- so a
    second positional argument sets ``external``, not ``crs``.

    Parameters
    ----------
    path : Path or str
        Path to the ``.shp`` / ``.geojson`` file. Relative paths are anchored at
        resolve time to the caller's ``project_root`` unless ``external`` is set.
    crs : str or int, optional
        A CRS to **assume when the file declares none** -- a fallback, not an
        override, and not a reprojection target. A file carrying its own CRS
        keeps it and this is ignored; grid sources are reprojected to the CRS on
        the :class:`~myflopy.specs.GridSpec`. Either spelling works:
        ``"EPSG:2927"`` or ``2927``.
    external : bool, default False
        ``True`` marks the file as outside the project tree: referenced in
        place, and its relative path is not anchored to ``project_root``.
    metadata : dict, optional
        Free-form JSON-serializable annotations. ``metadata["crs"]`` doubles as
        a last-resort CRS fallback.

    Examples
    --------
    >>> ShapeSource("domain_boundary.shp")
    >>> ShapeSource("streams.geojson", crs="EPSG:2927")   # .geojson has no CRS block

    See Also
    --------
    GeoPackageSourceSpec : Multi-layer GeoPackage, with a layer name and field map.
    """

    crs: str | int | None = None
    kind: ClassVar[str] = "ShapeSource"

    def to_dict(self) -> dict[str, Any]:
        """Serialize like the base spec, adding ``crs`` when set."""

        payload = DataSourceSpec.to_dict(self)
        if self.crs is not None:
            payload["crs"] = self.crs
        return payload

    @classmethod
    def _from_dict(cls, data: dict[str, Any]) -> ShapeSource:
        """Rebuild a :class:`ShapeSource` from a payload dict."""

        return cls(
            path=data["path"],
            crs=data.get("crs"),
            external=bool(data.get("external", False)),
            metadata=dict(data.get("metadata", {})),
        )


@dataclass(frozen=True, slots=True)
class GeoPackageSourceSpec(DataSourceSpec):
    """Serializable reference to one layer inside a GeoPackage (``.gpkg``).

    A ``.gpkg`` can hold many layers, so a bare path is ambiguous in a way a
    shapefile path is not. This spec pins **which** layer, **which** of its
    columns carry the values a reader needs, and optionally a subset of its
    rows -- all as plain data, so a recipe can be written to JSON/YAML and
    resolved again later.

    Use it wherever a source may be a multi-layer GeoPackage or needs a column
    mapping; use :class:`ShapeSource` for a single-layer vector file that needs
    neither. Today the only reader is **grid resolution** -- the ``boundary`` /
    ``refinement`` / ``breaklines`` / ``points`` sources of
    :meth:`~myflopy.specs.GridSpec.voronoi`. The runtime
    :class:`~myflopy.geopackage.GeoPackageSource` behind ``mf.ghb.gpkg(...)``
    is a separate class and does not read any of these fields; it names its
    columns with ``name_field`` / ``layer_field`` / ``period_field``.

    .. warning::
       Pass everything after ``path`` **by keyword**. This is a dataclass
       extending :class:`DataSourceSpec`, so the inherited fields come first in
       positional order -- ``(path, external, metadata, layer, query, crs,
       fields)``. ``GeoPackageSourceSpec("x.gpkg", "areas")`` therefore sets
       ``external="areas"``, which is truthy, rather than naming a layer.

    Parameters
    ----------
    path : Path or str
        Path to the ``.gpkg`` file. Relative paths are anchored at resolve time
        to the ``project_root`` given to
        :meth:`~myflopy.specs.GridSpec.resolve`, unless ``external`` is set.
        Coerced to a :class:`~pathlib.Path`.
    layer : str, optional
        Which layer to read. ``None`` reads the file's first/default layer,
        which is fine for a single-layer GeoPackage and a silent source of
        surprise in a multi-layer one -- name it explicitly when the file has
        more than one. The same file may be referenced by several specs with
        different ``layer`` values, which is how one GeoPackage supplies both
        refinement polygons and breakline centerlines.
    query : str, optional
        Row filter, applied **after** the layer is read. This is a
        :meth:`pandas.DataFrame.query` expression, **not** SQL: write
        ``"active == 1"``, not ``"active = 1"`` (the latter raises
        ``ValueError: cannot assign without a target object``). Column names
        are bare identifiers; quote strings, e.g. ``"kind == 'stream'"``.
        Because it filters after reading, it subsets features but does not
        reduce IO.
    crs : str or int, optional
        A CRS to **assume when the file declares none** -- a fallback, not an
        override. A GeoPackage that carries its own CRS keeps it, and passing
        something different here does nothing rather than reinterpreting the
        coordinates. Accepts either spelling: ``"EPSG:2927"`` or ``2927``.

        This is not the reprojection target either. Grid sources are reprojected
        to the CRS declared on the :class:`~myflopy.specs.GridSpec`, so a layer
        in EPSG:4326 inside a spec built with ``crs=2927`` is converted for you
        and ``crs=`` on this spec is not what does it. ``metadata["crs"]`` is
        consulted as a further fallback.
    fields : dict of str to str, optional
        Mapping of **logical key -> the column in this layer that holds it**.
        The KEY is the name myflopy looks up; the VALUE is your column. So
        ``{"area": "max_area"}`` reads each feature's target cell area from a
        column called ``max_area``.

        Grid sources understand exactly three keys, all optional:

        ``"area"``
            Target cell area for that feature, in CRS units squared. Falls back
            to ``refinement_max_area`` / ``breakline_max_area``; one of the two
            is required or resolution raises naming both.
        ``"label"``
            Region name used in diagnostics and in
            ``tri._prepared_regions``. Defaults to ``refinement_<i>`` /
            ``breakline_<i>``.
        ``"priority"``
            Who claims the overlap where regions cross. Defaults to
            ``refinement_priority`` (0) or ``breakline_priority`` (1).

        An unrecognised key is **ignored silently**, so a typo degrades to the
        global default rather than raising -- check the spelling here first when
        a per-feature value seems not to apply. A row whose mapped column is
        null or missing also falls back.
    external : bool, default False
        ``True`` marks the file as living outside the project tree: it is
        referenced in place and its relative path is **not** anchored to
        ``project_root``.
    metadata : dict, optional
        Free-form JSON-serializable annotations carried with the reference.
        ``metadata["crs"]`` doubles as a last-resort CRS fallback.

    Notes
    -----
    Serialization writes ``kind`` as ``"GeoPackageSource"`` -- the tag, not the
    class name -- and omits ``layer`` / ``query`` / ``crs`` / ``fields`` when
    unset, so ``to_dict``/``from_dict`` round-trip a minimal payload.

    Frozen and slotted: build a modified copy with
    :func:`dataclasses.replace`, not by assignment.

    Examples
    --------
    Refinement polygons carrying their own cell sizes:

    >>> GeoPackageSourceSpec("grid_refinement.gpkg", layer="areas",
    ...                      fields={"area": "max_area", "label": "name"})

    Two layers of one file feeding two different grid channels:

    >>> refine = GeoPackageSourceSpec("refine.gpkg", layer="areas",
    ...                               fields={"area": "max_area"})
    >>> creeks = GeoPackageSourceSpec("refine.gpkg", layer="centerlines")

    A subset of one layer, filtered with pandas syntax:

    >>> GeoPackageSourceSpec("bcs.gpkg", layer="boundaries",
    ...                      query="kind == 'ghb' and active == 1")

    A layer written without a CRS, told what it is:

    >>> GeoPackageSourceSpec("no_crs.gpkg", crs=2927)

    See Also
    --------
    ShapeSource : Single-layer vector file; no layer name, no field map.
    DataSourceSpec : The shared ``path``/``external``/``metadata`` base.
    myflopy.specs.GridSpec.voronoi : Where these sources are consumed.
    """

    layer: str | None = None
    query: str | None = None
    crs: str | int | None = None
    fields: dict[str, str] = field(default_factory=dict)
    kind: ClassVar[str] = "GeoPackageSource"

    def to_dict(self) -> dict[str, Any]:
        """Serialize like the base spec, adding ``layer``/``query``/``crs``/``fields`` when set."""

        payload = DataSourceSpec.to_dict(self)
        if self.layer is not None:
            payload["layer"] = self.layer
        if self.query is not None:
            payload["query"] = self.query
        if self.crs is not None:
            payload["crs"] = self.crs
        if self.fields:
            payload["fields"] = self.fields
        return payload

    @classmethod
    def _from_dict(cls, data: dict[str, Any]) -> GeoPackageSourceSpec:
        """Rebuild a :class:`GeoPackageSourceSpec` from a payload dict."""

        return cls(
            path=data["path"],
            layer=data.get("layer"),
            query=data.get("query"),
            crs=data.get("crs"),
            fields=dict(data.get("fields", {})),
            external=bool(data.get("external", False)),
            metadata=dict(data.get("metadata", {})),
        )


@dataclass(frozen=True, slots=True)
class RasterSource(DataSourceSpec):
    """Reference to a raster (GeoTIFF/IMG) sampled onto the grid for a surface or array.

    Use this for elevation surfaces (top/botm), recharge grids, or any gridded
    package input that lives in a raster. The raster is sampled at cell
    locations using ``method`` and, optionally, assigned to a named target via
    ``map_to``.

    Parameters
    ----------
    path
        Path to the raster file (``.tif`` / ``.img`` ...).
    band
        1-based band index to read (default: the first band).
    map_to
        Optional name of the package input/surface this raster populates.
    method
        Sampling/resampling method (e.g. ``"nearest"``, ``"linear"``). ``None``
        uses the caller's default.
    crs
        Optional CRS override; defaults to the raster's stored CRS.

    Examples
    --------
    >>> RasterSource("ground_surface.tif", map_to="top")
    >>> RasterSource("recharge_mm_yr.tif", band=1, method="linear")
    """

    band: int | None = None
    map_to: str | None = None
    method: str | None = None
    crs: str | int | None = None
    kind: ClassVar[str] = "RasterSource"

    def to_dict(self) -> dict[str, Any]:
        """Serialize like the base spec, adding any set ``band``/``map_to``/``method``/``crs``."""

        payload = DataSourceSpec.to_dict(self)
        for key in ("band", "map_to", "method", "crs"):
            value = getattr(self, key)
            if value is not None:
                payload[key] = value
        return payload

    @classmethod
    def _from_dict(cls, data: dict[str, Any]) -> RasterSource:
        """Rebuild a :class:`RasterSource` from a payload dict."""

        return cls(
            path=data["path"],
            band=data.get("band"),
            map_to=data.get("map_to"),
            method=data.get("method"),
            crs=data.get("crs"),
            external=bool(data.get("external", False)),
            metadata=dict(data.get("metadata", {})),
        )


@dataclass(frozen=True, slots=True)
class LiteralSource:
    """Inline literal value for a package input that has no external file.

    The escape hatch among the source specs: when an input is a plain scalar,
    list, or small mapping that you simply want to embed in a serialized recipe
    (rather than point at a file), wrap it in ``LiteralSource``. The ``value`` is
    stored as-is, so it must be JSON-serializable to round-trip. Unlike the
    file-backed specs it has no ``path`` and is not a :class:`DataSourceSpec`
    subclass.

    Parameters
    ----------
    value
        The literal payload (must be JSON-serializable to survive ``to_dict``).
    metadata
        Arbitrary JSON-serializable annotations.

    Examples
    --------
    >>> LiteralSource(1.0e-4)                 # a constant K value
    >>> LiteralSource({0: True})              # steady-state flags per period
    """

    value: Any
    metadata: dict[str, Any] = field(default_factory=dict)

    kind: ClassVar[str] = "LiteralSource"

    def to_dict(self) -> dict[str, Any]:
        """Serialize the inline value to a ``kind``-tagged dict (omitting empty metadata)."""

        payload = {
            "kind": self.kind,
            "value": self.value,
        }
        if self.metadata:
            payload["metadata"] = self.metadata
        return payload

    @classmethod
    def _from_dict(cls, data: dict[str, Any]) -> LiteralSource:
        """Rebuild a :class:`LiteralSource` from a payload dict."""

        return cls(
            value=data.get("value"),
            metadata=dict(data.get("metadata", {})),
        )


_SOURCE_TYPES: dict[str, type[Any]] = {
    DataSourceSpec.kind: DataSourceSpec,
    TableSource.kind: TableSource,
    ShapeSource.kind: ShapeSource,
    GeoPackageSourceSpec.kind: GeoPackageSourceSpec,
    RasterSource.kind: RasterSource,
    LiteralSource.kind: LiteralSource,
}


def source_from_dict(data: dict[str, Any]) -> DataSourceSpec | LiteralSource:
    """Recreate a source specification from its serialized representation."""

    kind = str(data.get("kind"))
    try:
        source_type = _SOURCE_TYPES[kind]
    except KeyError as error:
        raise ValueError(f"Unknown source kind: {kind!r}") from error
    return source_type._from_dict(data)


__all__ = [
    "DataSourceSpec",
    "GeoPackageSourceSpec",
    "LiteralSource",
    "RasterSource",
    "ShapeSource",
    "TableSource",
    "source_from_dict",
]
