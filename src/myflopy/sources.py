"""Serializable source specifications for project-backed model recipes."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, ClassVar


def _path_text(path: Path | str) -> str:
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
        object.__setattr__(self, "path", Path(self.path))

    def to_dict(self) -> dict[str, Any]:
        payload = {
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
        kind = data.get("kind")
        try:
            source_type = _SOURCE_TYPES[kind]
        except KeyError as error:
            raise ValueError(f"Unknown data source kind: {kind!r}") from error
        return source_type._from_dict(data)

    @classmethod
    def _from_dict(cls, data: dict[str, Any]) -> DataSourceSpec:
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
        payload = DataSourceSpec.to_dict(self)
        if self.table is not None:
            payload["table"] = self.table
        return payload

    @classmethod
    def _from_dict(cls, data: dict[str, Any]) -> TableSource:
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

    Parameters
    ----------
    path
        Path to the ``.shp`` / ``.geojson`` file.
    crs
        Optional coordinate reference system override (e.g. ``"EPSG:2927"``).
        When omitted, the file's own CRS is used.

    Examples
    --------
    >>> ShapeSource("domain_boundary.shp")
    >>> ShapeSource("streams.geojson", crs="EPSG:2927")
    """

    crs: str | None = None
    kind: ClassVar[str] = "ShapeSource"

    def to_dict(self) -> dict[str, Any]:
        payload = DataSourceSpec.to_dict(self)
        if self.crs is not None:
            payload["crs"] = self.crs
        return payload

    @classmethod
    def _from_dict(cls, data: dict[str, Any]) -> ShapeSource:
        return cls(
            path=data["path"],
            crs=data.get("crs"),
            external=bool(data.get("external", False)),
            metadata=dict(data.get("metadata", {})),
        )


@dataclass(frozen=True, slots=True)
class GeoPackageSourceSpec(DataSourceSpec):
    """Serializable reference to one layer/table inside a GeoPackage (``.gpkg``).

    The durable, recipe-friendly counterpart to the runtime
    :class:`~myflopy.geopackage.GeoPackageSource`. It pins which ``layer`` to
    read, an optional attribute ``query`` to subset features, a ``crs``
    override, and a ``fields`` mapping that renames source columns to the names
    a package builder expects (e.g. ``{"stage": "BHEAD", "cond": "COND"}``).
    This is what the ``mf.ghb.gpkg(...)`` / ``mf.drn.gpkg(...)`` forms record so
    a GIS-driven boundary condition can be rebuilt from the saved spec.

    Parameters
    ----------
    path
        Path to the ``.gpkg`` file.
    layer
        Layer/table name within the GeoPackage. ``None`` uses the first/default
        layer.
    query
        Optional attribute filter (an OGR/SQL ``WHERE`` expression) applied when
        reading features.
    crs
        Optional CRS override; defaults to the layer's stored CRS.
    fields
        Mapping of source column name -> target field name expected downstream.

    Examples
    --------
    >>> GeoPackageSourceSpec("bcs.gpkg", layer="ghb_cells",
    ...                      fields={"head": "bhead", "k": "cond"})
    """

    layer: str | None = None
    query: str | None = None
    crs: str | None = None
    fields: dict[str, str] = field(default_factory=dict)
    kind: ClassVar[str] = "GeoPackageSource"

    def to_dict(self) -> dict[str, Any]:
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
    crs: str | None = None
    kind: ClassVar[str] = "RasterSource"

    def to_dict(self) -> dict[str, Any]:
        payload = DataSourceSpec.to_dict(self)
        for key in ("band", "map_to", "method", "crs"):
            value = getattr(self, key)
            if value is not None:
                payload[key] = value
        return payload

    @classmethod
    def _from_dict(cls, data: dict[str, Any]) -> RasterSource:
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
        payload = {
            "kind": self.kind,
            "value": self.value,
        }
        if self.metadata:
            payload["metadata"] = self.metadata
        return payload

    @classmethod
    def _from_dict(cls, data: dict[str, Any]) -> LiteralSource:
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

    kind = data.get("kind")
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
