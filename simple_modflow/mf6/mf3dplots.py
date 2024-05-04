import plotly.graph_objs as go
import pandas as pd
import numpy as np
from figs import Fig
from simple_modflow.mf6.voronoiplus import VoronoiGridPlus


class SurfacePlot(Fig):

    def __init__(self, vor: VoronoiGridPlus=None):
        super().__init__()
        self.vor = vor

    def _get_gridded_xy_from_vor(self, spacing=10):

        if self.vor is None:
            return ValueError('No valid voronoi grid')
        else:
            vor = self.vor

        """Get Voronoi Grid Centroid x,y-coords"""
        centroids_xy = list(zip(vor.centroids_x, vor.centroids_y))

        """"Get gridded xs and ys of vornoi grid based on provided spacing"""
        x_surface, y_surface = vor.generate_grid_coordinates(spacing)
        xy = list(zip(x_surface, y_surface))

        unqiue_x = list(set(x_surface))
        unqiue_y = list(set(y_surface))
        xcol = unqiue_x.count(unqiue_x[0])
        yrow = unqiue_y.count(unqiue_y[0])


        """Make DataFrame of x,y,z coordinates, sort by y,
        and drop any N/A values"""
        df_xyz_gridded = pd.DataFrame(xyz_heads, columns=["x", "y", "head", "botm"])
        df_xyz_gridded.sort_values(by=["y", "x"], inplace=True, ascending=False)
        df_xyz_gridded.reset_index(drop=True, inplace=True)


    def _interpolate_z(self):

        """Interpolate z at gridded xy-data"""
        kstpkper = (14, 0)  # stress period to get data from
        gridded_heads = griddata(
            points=centroids_xy, values=heads.all_heads.loc[kstpkper, :]['elev'].to_list(), xi=xy, method="cubic"
        )
        gridded_model_botm = griddata(
            points=centroids_xy, values=vor.gdf_topbtm["bottom"], xi=xy, method="cubic"
        )
        xyz_heads = list(zip(x_surface, y_surface, gridded_heads, gridded_model_botm))
    def plot(self):

        """Reshape x,y,z arrays for Surface plot"""
        z_reshaped = df_xyz_gridded["z"].to_numpy().reshape(xcol, yrow)
        x_reshaped = df_xyz_gridded["x"].to_numpy().reshape(xcol, yrow)
        y_reshaped = df_xyz_gridded["y"].to_numpy().reshape(xcol, yrow)

        z_min = df_xyz_gridded["head"].min()
        z_max = df_xyz_gridded["head"].max()

        surface_scene = go.layout.Scene(zaxis_range=[400, 800])

        fig_surf = go.Figure(
            go.Surface(
                z=botm_reshaped,
                x=x_reshaped,
                y=y_reshaped,
                opacity=0.8,
                scene="scene",
                contours={
                    "z": {
                        "show": True,
                        "start": z_min,
                        "end": z_max,
                        "size": 2,
                        "highlight": True,
                        "width": 2,
                    }
                },
            )
        )
        fig_surf.update_scenes(patch=surface_scene)
