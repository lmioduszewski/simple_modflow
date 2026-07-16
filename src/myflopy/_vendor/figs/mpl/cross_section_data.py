from __future__ import annotations

"""
Data-loading and normalization helpers for Matplotlib cross sections.

This module is responsible for turning a variety of input sources into a
consistent dataframe shape that `figs.mpl.plot_cross_section(...)` can render.
Keeping this logic separate from plotting makes it easier to reuse the loading
path in notebooks, preprocessing steps, and future adapters without importing
the full plotting surface.
"""

from pathlib import Path
from typing import Any

import geopandas as gpd
import pandas as pd
import shapely as shp

GIS_EXTENSIONS = {".gpkg", ".geojson", ".json", ".shp"}
EXCEL_EXTENSIONS = {".xls", ".xlsx", ".xlsm", ".xlsb", ".ods"}
TABULAR_EXTENSIONS = EXCEL_EXTENSIONS | {".csv", ".tsv"}


def _normalize_tabular_frame(data: pd.DataFrame, x: str, y: str) -> pd.DataFrame:
    normalized = data.copy()
    if x in normalized.columns and y in normalized.columns:
        return normalized

    if x == "x" and y == "y" and len(normalized.columns) >= 2:
        normalized["x"] = normalized.iloc[:, 0]
        normalized["y"] = normalized.iloc[:, 1]
        return normalized

    raise ValueError(
        f"Cross-section data must include '{x}' and '{y}' columns or provide at least two columns."
    )


def _as_linestring(geometry) -> shp.LineString:
    if isinstance(geometry, shp.LineString):
        return geometry

    if isinstance(geometry, shp.MultiLineString):
        merged = shp.line_merge(geometry)
        if isinstance(merged, shp.LineString):
            return merged
        coords = []
        for line in geometry.geoms:
            line_coords = list(line.coords)
            if coords and coords[-1] == line_coords[0]:
                coords.extend(line_coords[1:])
            else:
                coords.extend(line_coords)
        return shp.LineString(coords)

    raise TypeError(f"Expected a LineString-like geometry, got {type(geometry)}")


def read_section_line(
    source: Path | str | shp.LineString | shp.MultiLineString | gpd.GeoDataFrame,
    layer: str | int | None = None,
) -> shp.LineString:
    """
    Normalize a cross-section line from a GIS source or Shapely geometry.

    Parameters
    ----------
    source : Path, str, LineString, MultiLineString, or GeoDataFrame
        Cross-section source geometry. Supported inputs are:

        - a file path to a GIS dataset such as `.gpkg`, `.shp`, or `.geojson`
        - a `geopandas.GeoDataFrame`
        - a Shapely `LineString`
        - a Shapely `MultiLineString`

        When a GIS file or `GeoDataFrame` contains multiple line features, the
        geometries are merged into a single logical section line.
    layer : str or int, optional
        GIS layer name or layer index to read when `source` is a file path.
        This is especially useful for multi-layer GeoPackages.

    Returns
    -------
    shapely.LineString
        Normalized section line geometry.

    Raises
    ------
    TypeError
        If `source` is not a supported line-like input.

    Notes
    -----
    This function is intentionally geometry-focused. If you want a dataframe
    ready for plotting, use `line_to_profile_frame(...)` or
    `load_cross_section_data(...)` instead.
    """
    if isinstance(source, (shp.LineString, shp.MultiLineString)):
        return _as_linestring(source)

    if isinstance(source, gpd.GeoDataFrame):
        geometry = source.geometry.union_all()
        return _as_linestring(geometry)

    if isinstance(source, (Path, str)):
        read_kwargs = {}
        if layer is not None:
            read_kwargs["layer"] = layer
        gdf = gpd.read_file(Path(source), **read_kwargs)
        geometry = gdf.geometry.union_all()
        return _as_linestring(geometry)

    raise TypeError(f"Unsupported section source type: {type(source)}")


def line_to_profile_frame(
    source: Path | str | shp.LineString | shp.MultiLineString | gpd.GeoDataFrame,
    simplify_tolerance: float | None = None,
    layer: str | int | None = None,
) -> pd.DataFrame:
    """
    Convert a section line into a profile dataframe.

    Parameters
    ----------
    source : Path, str, LineString, MultiLineString, or GeoDataFrame
        Cross-section line geometry. This accepts the same source types as
        `read_section_line(...)`.
    simplify_tolerance : float, optional
        Optional Shapely simplification tolerance applied before extracting
        coordinates. This is useful when a GIS-exported line contains more
        vertices than are needed for a printed section.
    layer : str or int, optional
        GIS layer name or layer index to read when `source` is a file path.

    Returns
    -------
    pandas.DataFrame
        Dataframe with these columns:

        - `x`: horizontal coordinate from the source geometry
        - `y`: vertical coordinate from the source geometry
        - `distance`: cumulative distance along the line
        - `elevation`: alias for `y`, provided for convenience

    Notes
    -----
    This helper is intended for geometry-derived sections, such as a QGIS
    section line exported to GeoPackage. It does not sample terrain or raster
    values; it only converts existing line coordinates into a plottable profile.
    """
    line = read_section_line(source, layer=layer)
    if simplify_tolerance is not None:
        line = _as_linestring(line.simplify(simplify_tolerance))

    coords = list(line.coords)
    rows = []
    distance = 0.0
    for i, (x, y) in enumerate(coords):
        if i > 0:
            px, py = coords[i - 1]
            distance += ((x - px) ** 2 + (y - py) ** 2) ** 0.5
        rows.append(
            {
                "x": x,
                "y": y,
                "distance": distance,
                "elevation": y,
            }
        )
    return pd.DataFrame(rows)


def load_cross_section_data(
    source: pd.DataFrame | Path | str | shp.LineString | shp.MultiLineString | gpd.GeoDataFrame,
    x: str = "x",
    y: str = "y",
    layer: str | int | None = None,
    simplify_tolerance: float | None = None,
    sheet_name: str | int = 0,
    read_kwargs: dict[str, Any] | None = None,
) -> pd.DataFrame:
    """
    Load cross-section data from a dataframe, GIS source, or tabular file.

    Parameters
    ----------
    source : DataFrame, Path, str, LineString, MultiLineString, or GeoDataFrame
        Cross-section source. Supported forms are:

        - `pandas.DataFrame`
        - path to a GIS file such as `.gpkg`, `.shp`, or `.geojson`
        - path to an Excel file such as `.xlsx` or `.xls`
        - path to a delimited text file such as `.csv` or `.tsv`
        - Shapely `LineString` or `MultiLineString`
        - `geopandas.GeoDataFrame`

    x : str, default "x"
        Name of the x column to use for tabular sources. If a dataframe or file
        does not contain this column and `x="x"` with `y="y"`, the first two
        columns are assumed to be x and y respectively.
    y : str, default "y"
        Name of the y column to use for tabular sources. If a dataframe or file
        does not contain this column and `x="x"` with `y="y"`, the first two
        columns are assumed to be x and y respectively.
    layer : str or int, optional
        GIS layer name or index to read when `source` points to a GIS dataset.
    simplify_tolerance : float, optional
        Optional simplification tolerance applied to GIS-derived geometry before
        conversion to a dataframe.
    sheet_name : str or int, default 0
        Excel sheet name or zero-based sheet index when reading Excel files.
    read_kwargs : dict, optional
        Additional keyword arguments passed to the underlying pandas reader for
        tabular file sources.

    Returns
    -------
    pandas.DataFrame
        Normalized dataframe ready for `plot_cross_section(...)`.

        For geometry-derived sources, the returned dataframe includes `x`, `y`,
        `distance`, and `elevation`.

        For tabular sources, the returned dataframe preserves the original
        columns and guarantees that the requested x and y columns are available.

    Raises
    ------
    TypeError
        If `source` is not one of the supported input types.
    ValueError
        If a tabular source does not provide the requested x/y columns and does
        not have at least two columns available for fallback assignment.
    """
    read_kwargs = read_kwargs or {}

    if isinstance(source, (shp.LineString, shp.MultiLineString, gpd.GeoDataFrame)):
        return line_to_profile_frame(
            source,
            simplify_tolerance=simplify_tolerance,
            layer=layer,
        )

    if isinstance(source, pd.DataFrame):
        return _normalize_tabular_frame(source, x=x, y=y)

    if isinstance(source, (Path, str)):
        source_path = Path(source)
        suffix = source_path.suffix.lower()
        if suffix in GIS_EXTENSIONS:
            return line_to_profile_frame(
                source_path,
                simplify_tolerance=simplify_tolerance,
                layer=layer,
            )
        if suffix in EXCEL_EXTENSIONS:
            frame = pd.read_excel(source_path, sheet_name=sheet_name, **read_kwargs)
            return _normalize_tabular_frame(frame, x=x, y=y)
        if suffix in TABULAR_EXTENSIONS:
            if suffix == ".csv":
                frame = pd.read_csv(source_path, **read_kwargs)
            elif suffix == ".tsv":
                frame = pd.read_csv(source_path, sep="\t", **read_kwargs)
            else:
                frame = pd.read_table(source_path, **read_kwargs)
            return _normalize_tabular_frame(frame, x=x, y=y)

    raise TypeError(f"Unsupported cross-section source type: {type(source)}")


__all__ = [
    "line_to_profile_frame",
    "load_cross_section_data",
    "read_section_line",
]
