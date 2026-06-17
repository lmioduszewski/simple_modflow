"""Serializable source specifications for project-backed model recipes."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, ClassVar


def _path_text(path: Path | str) -> str:
    return str(Path(path).as_posix())


@dataclass(frozen=True, slots=True)
class DataSourceSpec:
    """Base class for durable references to package/grid input data."""

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
    """Tabular source used to build package records."""

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
    """Vector file source used for grid or package construction."""

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
    """Serializable reference to one layer/table inside a GeoPackage."""

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
    """Raster source used for gridded package inputs or surfaces."""

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
    """Named literal value used when a package input is intentionally inline."""

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
