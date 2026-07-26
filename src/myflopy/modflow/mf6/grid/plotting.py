"""Plotting helpers for Voronoi grids and grid-derived cross sections."""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING

import geopandas as gpd
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import shapely as shp
from flopy.discretization.vertexgrid import VertexGrid
from flopy.plot.crosssection import PlotCrossSection

from myflopy import viz as f
from myflopy.modflow.utils.datatypes.choros import Choro
from myflopy.modflow.utils.datatypes.readers import read_shp_gpkg

if TYPE_CHECKING:
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
    from myflopy.modflow.mf6.simulation.base import SimulationBase


def mapit(vor, crs: str = 'EPSG:2927'):
    """
    Explore the Voronoi polygons in an interactive map.
    """
    poly = vor.get_voronoi_polygons()
    return gpd.GeoDataFrame(geometry=poly, crs=crs).explore()


def build_choropleth(
    vor,
    model: SimulationBase = None,
    kstpkper: tuple = None,
    per: int = None,
    layer: int = 0,
    type: str = 'hds',
    custom_hover: dict = None,
    custom_zs: list = None,
    zmin: float | int = None,
    zmax: float | int = None,
    zoom: int = 13,
    show_layer_elevs: bool = False,
    show_mounding: bool = False,
    hover_heads: bool = True,
    hover_ks: bool = False,
    locs: Path = None,
    colorscale: str | list | tuple = None,
    logscale: bool = False,
    hover_spec=None,
    **choro_kwargs,
) -> Choro:
    """
    Create a Choro wrapper for Voronoi plotting.

    ``show_layer_elevs`` defaults to ``False`` here and to ``True`` on
    :class:`Choro` on purpose: this is the grid-only front door (``model=None``),
    where a ``vor`` carrying no layer elevations makes Choro's own default raise
    ``AttributeError`` as soon as the hover is built. Do not "simplify" this into
    a bare passthrough.

    ``colorscale``/``logscale``/``hover_spec`` are named because callers reach
    for them constantly; anything else in ``**choro_kwargs`` rides through to the
    ``go.Choroplethmap`` trace (``zmid``, ``colorbar``, ``reversescale``, ...).
    Those are validated LATE, by Plotly at ``plot()`` time, not here -- ``title=``
    in particular is a matplotlib-only argument and raises there.
    """
    return Choro(
        vor=vor,
        model=model,
        kstpkper=kstpkper,
        per=per,
        layer=layer,
        type=type,
        custom_hover=custom_hover,
        custom_zs=custom_zs,
        zmin=zmin,
        zmax=zmax,
        zoom=zoom,
        show_layer_elevs=show_layer_elevs,
        show_mounding=show_mounding,
        hover_heads=hover_heads,
        hover_ks=hover_ks,
        locs=locs,
        colorscale=colorscale,
        logscale=logscale,
        hover_spec=hover_spec,
        **choro_kwargs,
    )


def get_dash_selector(vor):
    """Return the Dash selector widget from the default choropleth."""
    return build_choropleth(vor).dash_selector()


def show(vor):
    """Display the default Voronoi choropleth plot."""
    return build_choropleth(vor).plot()


def map_nodes(vor) -> go.Figure:
    """Build an interactive map showing cell polygons keyed by node ID."""
    latlonselect = vor.gdf_vorPolys.to_crs('EPSG:4326')
    latlonselect['cellidx'] = latlonselect.index.astype(str)
    geojsonselect = json.loads(latlonselect['geometry'].to_json())
    centroid_grid = vor.grid_centroid
    fig_sel = go.Figure(go.Choroplethmap(
        geojson=geojsonselect,
        locations=latlonselect['cellidx'].to_list(),
        featureidkey='id',
        z=latlonselect['cellidx'].to_list(),
        colorscale='earth'
    ))

    fig_sel.update_layout(
        margin={"r": 0, "t": 20, "l": 0, "b": 0},
        map_style="carto-positron",
        map_zoom=15,
        map_center={"lat": centroid_grid.y, "lon": centroid_grid.x},
    )

    return fig_sel


def show_selected_cells(vor, cell_list: list = None, **kwargs):
    """
    Show selected cells on the Voronoi choropleth.
    """
    choro = build_choropleth(vor, **kwargs).choropleth
    choro.data[0].selectedpoints = tuple(cell_list)
    return go.Figure(choro).show(renderer='browser')


def show_overlapping_geometry(vor, shp_gpkg):
    """
    Show cells overlapping the provided geometry.
    """
    cells = vor.get_vor_cells_as_series(shp_gpkg).to_list()
    return show_selected_cells(vor, cells)


def plot3d(vor, z=None) -> go.Figure:
    """Build a simple 3D view of Voronoi edges or centroids."""
    if hasattr(vor, 'x_vor_regions') and hasattr(vor, 'y_vor_regions'):
        x_lines = vor.x_vor_regions
        y_lines = vor.y_vor_regions
    else:
        x_lines = []
        y_lines = []
        for xs, ys in zip(vor.x_coords_by_node, vor.y_coords_by_node, strict=False):
            x_lines.extend(list(xs) + [None])
            y_lines.extend(list(ys) + [None])

    if z is None:
        z = [0 for _ in x_lines]
    fig3d = go.Figure(
        data=go.Scatter3d(
            x=x_lines,
            y=y_lines,
            z=z,
            opacity=0.5,
            mode='lines',
            line_color='black',
        ),
        layout={'height': 1000},
    )
    if all(hasattr(vor, attr) for attr in ('x_vor', 'y_vor', 'i', 'j', 'k')):
        fig3d.add_trace(
            go.Mesh3d(
                x=vor.x_vor,
                y=vor.y_vor,
                z=[0 for _ in range(len(vor.x_vor))],
                i=vor.i,
                j=vor.j,
                k=vor.k,
                colorscale='Viridis',
                intensity=vor.x_vor,
            )
        )
    else:
        cx, cy = vor.centroids
        fig3d.add_trace(
            go.Scatter3d(
                x=cx,
                y=cy,
                z=[0 for _ in cx],
                mode='markers',
                marker=dict(size=4, color='royalblue'),
                name='centroids',
            )
        )
    fig3d.update_layout(
        title='Voronoi Diagram',
        scene=dict(
            xaxis=dict(title='X'),
            yaxis=dict(title='Y'),
            zaxis=dict(title=''),
        ),
    )
    return fig3d


def plot2d(vor) -> go.Figure:
    """Build a 2D Plotly view of the Voronoi cell edges."""
    fig2d = go.Figure(layout=vor.scatt_layout)
    for cell in range(len(vor.x_coords_by_node)):
        fig2d.add_scattergl(
            x=vor.x_coords_by_node[cell],
            y=vor.y_coords_by_node[cell],
            opacity=1,
            mode='lines',
            line_color='black',
            line_width=1,
        )

    vor.fig2d = fig2d
    return fig2d


def plottri(vor):
    """
    Plot the triangulated mesh generated by Triangle.
    """
    df_tricells = np.asarray(vor.tri.get_cell2d())
    ilist = df_tricells[:, 4].astype(int).tolist()
    jlist = df_tricells[:, 5].astype(int).tolist()
    klist = df_tricells[:, 6].astype(int).tolist()
    xtriverts = np.asarray(vor.tri.verts)[:, 0].tolist()
    ytriverts = np.asarray(vor.tri.verts)[:, 1].tolist()

    xtricells = []
    ytricells = []
    for i, j, k in zip(ilist, jlist, klist, strict=False):
        xtricells.extend([xtriverts[i], xtriverts[j], xtriverts[k], None])
        ytricells.extend([ytriverts[i], ytriverts[j], ytriverts[k], None])

    trifig2d = go.Figure(
        go.Scattergl(
            x=xtricells,
            y=ytricells,
            mode='lines',
            line_color='black',
            line_width=1,
        ),
        layout=vor.scatt_layout,
    )
    return trifig2d.show(config=vor.config)


def build_grid_section(vor: VoronoiGridPlus, line: shp.LineString | Path):
    """
    Build a cross-section helper for the Voronoi grid.
    """
    return GridSection(vor=vor, line=line)


def _as_linestring(geometry) -> shp.LineString:
    """Coerce a ``LineString`` or ``MultiLineString`` to a single merged ``LineString`` (raises otherwise)."""

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

    raise ValueError(f'line arg must resolve to a LineString, not {type(geometry)}')


class GridSection:
    """
    Represents a section of a grid and provides tools for creating and plotting
    cross-sections.
    """

    def __init__(self, vor, line: shp.LineString | shp.MultiLineString | Path):
        """Build a cross-section of grid ``vor`` along ``line`` (a geometry or vector file)."""

        self.vor = vor
        props = vor.get_disv_gridprops()
        self.grid = VertexGrid(
            vertices=props['vertices'],
            top=vor.gdf_topbtm[0].values,
            botm=vor.gdf_topbtm.loc[:, 1:].values.T,
            cell2d=props['cell2d'],
            lenuni='feet',
            ncpl=props['ncpl'],
            crs=vor.crs,
            nlay=vor.nlay,
        )

        if isinstance(line, Path):
            geometry = read_shp_gpkg(line).union_all()
            self.coords = _as_linestring(geometry).coords
        elif isinstance(line, (shp.LineString, shp.MultiLineString)):
            self.coords = _as_linestring(line).coords
        else:
            raise ValueError(f'line arg must be a Path or LineString, not {type(line)}')

        self.xy = np.array([xy for xy in self.coords])

    @property
    def polys(self):
        """Return FloPy cross-section polygons for the selected line."""
        polys = PlotCrossSection(
            modelgrid=self.grid,
            line={'line': self.xy},
        ).polygons
        return polys

    @property
    def poly_coords(self):
        """Return raw vertex arrays for each cross-section polygon."""
        poly_coords = []
        for poly in self.polys.values():
            verts = poly[0].get_xy()
            poly_coords.append(verts)
        return poly_coords

    def to_frame(self) -> pd.DataFrame:
        """
        Return section polygon outlines as a long-form DataFrame.
        """
        rows = []
        for i, verts in enumerate(self.poly_coords):
            series_name = f"polygon_{i}"
            for distance, elevation in verts:
                rows.append(
                    {
                        "distance": float(distance),
                        "elevation": float(elevation),
                        "series": series_name,
                        "polygon_id": i,
                    }
                )
        return pd.DataFrame(rows)

    def plot_mpl(self, **kwargs):
        """
        Plot the grid cross-section with the figs matplotlib cross-section helper.
        """
        from myflopy.viz import plot_cross_section

        data = self.to_frame()
        kwargs.setdefault("show_legend", False)
        return plot_cross_section(
            data=data,
            x="distance",
            y="elevation",
            series_col="series",
            **kwargs,
        )

    @property
    def figure(self):
        """Build a Plotly figure of the current cross section."""
        fig = f.Fig()
        for verts in self.poly_coords:
            xs, ys = verts[:, 0], verts[:, 1]
            fig.add_trace(
                go.Scatter(
                    x=xs,
                    y=ys,
                    fill='toself',
                    mode='lines',
                    line=dict(color='black'),
                    name='Polygon',
                )
            )
        return fig

    def plot(self):
        """Display the current cross-section figure."""
        self.figure.show()
