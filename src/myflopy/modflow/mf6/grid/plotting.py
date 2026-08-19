"""Plotting helpers for Voronoi grids and grid-derived cross sections."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import shapely as shp
from flopy.discretization.vertexgrid import VertexGrid
from flopy.plot.crosssection import PlotCrossSection

from myflopy import viz as f
from myflopy.viz import Picture
from myflopy.modflow.utils.datatypes.choros import Choro
from myflopy.modflow.utils.datatypes.readers import read_shp_gpkg

if TYPE_CHECKING:
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
    from myflopy.modflow.mf6.simulation.base import SimulationBase


def _choropleth_factory(
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
    show_layer_elevs: bool | None = None,
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

    ``show_layer_elevs=None`` (the default) means **decide from the grid**: the
    layer-elevation hover needs ``vor.gdf_topbtm``, and a grid without it makes
    Choro's own ``True`` default raise ``AttributeError`` as soon as the hover is
    built. Pass ``True``/``False`` to force it.

    That resolution used to live at the CALL SITES -- 25 of them repeated
    ``kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(model))``
    around a hardcoded ``False`` here, and the two that forgot silently lost five
    hover rows. It belongs in one place, next to the attribute it depends on. Do
    not "simplify" this back into a bare passthrough.

    ``colorscale``/``logscale``/``hover_spec`` are named because callers reach
    for them constantly; anything else in ``**choro_kwargs`` rides through to the
    ``go.Choroplethmap`` trace (``zmid``, ``colorbar``, ``reversescale``, ...).
    Those are validated LATE, by Plotly at ``plot()`` time, not here -- ``title=``
    in particular is a matplotlib-only argument and raises there.
    """
    if show_layer_elevs is None:
        show_layer_elevs = getattr(vor, "gdf_topbtm", None) is not None
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


def _grid_section_factory(vor: VoronoiGridPlus, line: shp.LineString | Path):
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


class GridSection(Picture):
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
    def fig(self):
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



class GridMesh(Picture):
    """The bare mesh -- cell edges, no values, no basemap (plan 8.4).

    The picture ``map()`` cannot give you. A :class:`Choro` colours cells against
    a web basemap and so hard-requires a CRS; this draws the geometry in the
    grid's own coordinates and needs none, which is what makes it the right view
    for a grid you are still refining, before there is a model or a projection to
    speak of.

    Absorbs the old ``vor.plot2d()`` (Plotly) and the inherited FloPy
    ``VoronoiGrid.plot()`` (Matplotlib, still reachable as :meth:`plot_mpl`).
    """

    def __init__(self, vor):
        """Bind a mesh view to grid ``vor``."""

        self.vor = vor

    @property
    def fig(self):
        """Cell edges as one Plotly trace per cell, in the grid's own coordinates."""

        fig = f.Fig(layout=getattr(self.vor, "scatt_layout", None))
        for cell in range(len(self.vor.x_coords_by_node)):
            fig.add_scattergl(
                x=self.vor.x_coords_by_node[cell],
                y=self.vor.y_coords_by_node[cell],
                opacity=1,
                mode='lines',
                line_color='black',
                line_width=1,
                showlegend=False,
            )
        return fig

    def plot_mpl(self, ax=None, plot_title: bool = True, **kwargs):
        """Render with Matplotlib, via FloPy's own patch-collection renderer.

        This IS ``VoronoiGrid.plot`` -- called unbound, because ``vor.plot`` is
        now the namespace this object came from. Nothing about the drawing
        changed; only the spelling did.
        """

        from flopy.utils.voronoi import VoronoiGrid

        return VoronoiGrid.plot(self.vor, ax=ax, plot_title=plot_title, **kwargs)


class GridPlots:
    """The plotting verbs for a Voronoi grid -- ``vor.plot.map()`` (plan 8.4).

    Three verbs, because a bare grid can answer three questions: what are the
    values over it (:meth:`map`), what does it look like in section
    (:meth:`section`), and what does the mesh itself look like (:meth:`grid`).
    Everything they return is a :class:`~myflopy.viz.Picture`.

    Calling the namespace is shorthand for the mesh: ``vor.plot()`` is
    ``vor.plot.grid()``. That is deliberate -- ``vor.plot`` used to BE FloPy's
    ``VoronoiGrid.plot()``, and "draw the grid" is what people already reach for
    that name to do.

    Replaces eleven aliases. Five were genuinely distinct pictures and survive as
    these verbs or as options on them; the rest were dead, duplicated, or not
    pictures at all -- see the ledger.

    Lives here rather than in :mod:`myflopy.plot` for a layering reason:
    ``voronoi.py`` sits BELOW that module in the import graph, so binding this
    from there would point upward and force a deferred import.
    """

    def __init__(self, vor):
        """Bind the plotting verbs to grid ``vor``."""

        self.vor = vor

    def __repr__(self):
        """Name the verbs, since tab-completion is how this gets discovered."""

        return f"GridPlots({getattr(self.vor, 'ncpl', '?')} cells: map, section, grid)"

    def __call__(self, **kwargs) -> GridMesh:
        """``vor.plot()`` -> the mesh. See :meth:`grid`."""

        return self.grid(**kwargs)

    def map(self, values=None, *, select=None, **kwargs) -> Choro:
        """A plan-view map of this grid, on a basemap.

        ``values`` is any per-cell array; with none, cells are keyed by node id
        -- which is what the old ``map_nodes()`` drew. ``select`` highlights a
        subset, either as cell indices or as a vector file/geometry to intersect,
        replacing ``show_selected_cells()`` and ``show_overlapping_geometry()``.

        Requires a CRS, because the basemap does. Use :meth:`grid` for a grid
        that does not have one yet.
        """

        if values is not None:
            kwargs["custom_zs"] = list(values)
        picture = _choropleth_factory(self.vor, **kwargs)
        if select is not None:
            if isinstance(select, (str, Path)) or hasattr(select, "geom_type"):
                select = self.vor.get_vor_cells_as_series(select).to_list()
            picture.fig.data[0].selectedpoints = tuple(select)
        return picture

    def section(self, line) -> GridSection:
        """A vertical slice of the grid geometry along ``line``.

        Layers and cell edges, no results -- for a section through a model's
        results use ``model.plot.section(...)``.
        """

        return _grid_section_factory(self.vor, line=line)

    def grid(self, **kwargs) -> GridMesh:
        """The bare mesh: cell edges, no values, no basemap, no CRS needed."""

        return GridMesh(self.vor, **kwargs)
