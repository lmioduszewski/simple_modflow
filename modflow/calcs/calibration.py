from __future__ import annotations
from typing import TYPE_CHECKING

from simple_modflow import SimulationBase

if TYPE_CHECKING:
    import plotly.graph_objects as go
    import pandas as pd

import numpy as np
import figs as f
from pathlib import Path
from pandas import IndexSlice as idxx
import pandas as pd


def calculate_calibration_statistics(observed, simulated):
    """
    Calculate common calibration statistics.

    Parameters:
    observed (np.array): Array of observed values.
    simulated (np.array): Array of simulated values.

    Returns:
    dict: Dictionary of calibration statistics.
    """
    if len(observed) != len(simulated):
        raise ValueError(f"Observed and simulated arrays must have the same length."
                         f"simulated length: {len(simulated)}, observed length: {len(observed)}")

    # Mean Error (ME)
    me = np.mean(simulated - observed)

    # Mean Absolute Error (MAE)
    mae = np.mean(np.abs(simulated - observed))

    # Root Mean Square Error (RMSE)
    rmse = np.sqrt(np.mean((simulated - observed) ** 2))

    # Normalized Root Mean Square Error (NRMSE)
    observed_range = np.max(observed) - np.min(observed)
    nrmse = rmse / observed_range

    # Nash-Sutcliffe Efficiency (NSE)
    nse = 1 - (np.sum((simulated - observed) ** 2) / np.sum((observed - np.mean(observed)) ** 2))

    # Coefficient of Determination (R^2)
    r_squared = np.corrcoef(observed, simulated)[0, 1] ** 2

    stats = {
        'ME': me,
        'MAE': mae,
        'RMSE': rmse,
        'NRMSE': nrmse,
        'NSE': nse,
        'R^2': r_squared
    }

    return stats


class CalibrationPlot(f.Fig):

    def __init__(
            self,
            observed: pd.Series | list | Path = None,
            simulated: pd.Series | list = None,
            model: SimulationBase = None,
            loc_cell_dict: dict = None,
            loc_shp_gpkg: Path = None,
            cal_range: list | tuple = None,
            exploration_name_field: str = 'ExploName',
            crs=None,
            obs_layers: int | list[int] = None
    ):
        """
        Initializes an object to plot observed and simulated calibration data.

        :param observed: Observed data, provided as a pandas Series, a list, or a Path
            object pointing to the data.
        :param simulated: Simulated model data, provided as a pandas Series or a list.
        :param model: Simulation model object adhering to the SimulationBase class.
        :param loc_cell_dict: Dictionary linking geospatial locations to cell mappings
            or identifiers.
        :param loc_shp_gpkg: Path to shapefile or GeoPackage used for geospatial
            location definitions and data extraction.
        :param cal_range: List or tuple specifying the calibration range, typically as
            time steps or indices for calibration analysis.
        :param exploration_name_field: Field name in the geospatial file (e.g.,
            shapefile or GeoPackage) representing exploration or location names.
        :param crs: Coordinate reference system used for geospatial data processing
            and mapping.
        :param obs_layers: Integer or list of integers representing observation layers
        """
        super().__init__()
        self._observed = None
        self._simulated = None
        self._model = None
        self._crs = crs
        self._cal_range = cal_range

        self._loc_cell_dict = loc_cell_dict
        self._loc_shp_gpkg = None
        self._exploration_name_field = exploration_name_field
        self._obs_layers = obs_layers

        self.model = model
        self.loc_shp_gpkg = loc_shp_gpkg
        self.simulated = simulated
        self.observed = observed

        if self._loc_shp_gpkg is not None:
            sim = self.model.hds.get_obs_heads(
                locs=self._loc_shp_gpkg,
                loc_name_field=self._exploration_name_field,
                long_format=True
            )
            if self._simulated is not None:
                print('replacing simulated data with obs heads based on'
                      'the provided shapefile or geopackage')
            self._simulated = sim

        if self.observed is not None and self.simulated is not None:
            self.add_calib_scatter()
            self.add_stats()
        self.set_layout()

        """if self._observed is not None and self._simulated is not None:
            self.add_calib_annotation(self._observed, self._simulated)"""

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, value):
        if value is not None:
            assert isinstance(value, SimulationBase), 'model must be a SimulationBase object'
        self._model = value

    @property
    def loc_shp_gpkg(self):
        return self._loc_shp_gpkg

    @loc_shp_gpkg.setter
    def loc_shp_gpkg(self, value):
        if value is not None:
            assert isinstance(value, Path), 'loc_shp_gpkg must be a Path object'
        self._loc_shp_gpkg = value

    @property
    def simulated(self):
        if self._simulated is None and self.loc_shp_gpkg and self.model:
            print('generating simulated data using model object and spatial data provided')
            simulated = self.model.hds.get_obs_heads(locs=self.loc_shp_gpkg, long_format=True)
            self._simulated = simulated
        return self._simulated

    @simulated.setter
    def simulated(self, value):
        if self.loc_shp_gpkg is not None:
            print('path to shapefile or geopackage provided, '
                  'so provided simulated data is ignored, and'
                  'will generate observed data using model'
                  'object and spatial data provided')
        else:
            self._simulated = value

    @property
    def observed(self):
        return self._observed

    @observed.setter
    def observed(self, value):
        if isinstance(value, Path):
            try:
                obs_data = pd.read_excel(value)
                if self._obs_layers is None:
                    self._obs_layers = 0  # default to layer 1
                col1 = obs_data.columns[0]
                obs_data = obs_data[sorted(obs_data.columns)]
                obs_data['layer'] = self._obs_layers
                obs_data = obs_data.set_index([col1, 'layer'])
                obs_data = obs_data.melt(ignore_index=False, value_name='elev', var_name='locs')
                obs_data = obs_data.reset_index().set_index(['locs', 'layer', col1])
                value = obs_data
            except ValueError:
                print('observed data path not readable, must be Excel file')
        self._observed = value

    def add_stats(
            self,
            observed: list = None,
            simulated: list = None,
    ):
        """
        Computes calibration statistics between observed and simulated data, formats
        them into a string for annotation, and adds the annotation to a plotly figure.

        The function optionally accepts 'observed' and 'simulated' data as input
        parameters. If not provided, default attributes of the class, namely
        'self.observed' and 'self.simulated', are used. The calibration statistics
        are computed using the 'calculate_calibration_statistics' utility and are
        annotated at the upper-left corner of the plot.

        :param observed: A list of observed values. If not provided, defaults to the
                         attribute 'self.observed'. Can also accept a pandas DataFrame,
                         in which case, the 'elev' column is extracted as a list.
        :param simulated: A list of simulated values. If not provided, defaults to the
                          attribute 'self.simulated'. Can also accept a pandas DataFrame,
                          in which case, the 'elev' column is extracted as a list.
        :return: A string representation of the calibration statistics formatted
                 with line breaks. Returns None when there is a mismatch in the length
                 or data type of the 'observed' and 'simulated' inputs.
        """

        observed = self.observed if observed is None else observed
        simulated = self.simulated if simulated is None else simulated

        if all(isinstance(x, pd.DataFrame) for x in [observed, simulated]):
            print('observed and simulated data are provided as DataFrames')
            obs = observed['elev'].to_list()
            sim = simulated['elev'].to_list()
            w = pd.concat([pd.Series(sim), pd.Series(obs)], axis=1).dropna()
            simulated = w.loc[:, 0]
            observed = w.loc[:, 1]

        stats = calculate_calibration_statistics(observed, simulated)
        calib_stats = [f'{k}: {v:.3f}' for k, v in stats.items()]

        # Create a single string with line breaks
        stats_text = "<br>".join(calib_stats)

        # Add annotation in the upper-left area
        self.add_annotation(
            x=0.02,  # small offset from the left border
            y=0.98,  # near top
            xref="paper",  # 'paper' means figure reference (0=left, 1=right)
            yref="paper",  # 'paper' means figure reference (0=bottom, 1=top)
            text=stats_text,
            showarrow=False,
            align="left",  # text alignment
            bordercolor="black",
            borderwidth=1,
            borderpad=4,
            bgcolor="white",
            opacity=0.8
        )
        return stats_text

    def add_calib_scatter(self):
        """
        Adds scatter plot elements for calibration visualization in the plot. This method
        handles two scenarios: calibration using simple datasets (e.g., lists or
        pandas Series) and calibration using observation data linked to shapefiles
        or geopackages. If calibration data is provided with a shapefile or geopackage,
        the observed and calibrated data provided to class must be in DataFrame format.

        :param self: The instance of the class where this method is defined.
        :return: None
        """
        cal_range = self._cal_range if self._cal_range is not None else None

        # if observed and simulated data are simple lists or Series, add all data as one Scattergl trace
        if all([isinstance(x, pd.Series | list) for x in [self._observed, self.simulated]]):
            self.add_scattergl(x=self._observed, y=self._simulated, mode='markers', name='data')
            data_min = self._observed.min()
            data_max = self._observed.max()

        # if a shapefile or geopackage was provided, add each observation
        # location as a separate Scattergl trace'
        if self.loc_shp_gpkg is not None:
            data_min = None
            data_max = None
            print('using observation data based on provided'
                  'shapefile or geopackage for calibration plot')
            assert isinstance(self.simulated, pd.DataFrame)
            loc_names = self.simulated.index.get_level_values(level=0).unique().to_list()
            layers = self.simulated.index.get_level_values(level=1).unique().to_list()

            for loc_name in loc_names:
                for layer in layers:

                    sim_data = self.simulated.loc[idxx[loc_name, layer, :], 'elev']
                    obs_data = self.observed.loc[idxx[loc_name, layer, :], 'elev']

                    # get max and min of obs data to draw cal line
                    if data_min is None:
                        data_min = obs_data.min()
                    if data_max is None:
                        data_max = obs_data.max()
                    if obs_data.min() < data_min:
                        data_min = obs_data.min()
                    if obs_data.max() > data_max:
                        data_max = obs_data.max()

                    self.add_scattergl(
                        x=sim_data,
                        y=obs_data,
                        mode='markers',
                        name=f'{loc_name} {layer}')

        if cal_range is None:
            cal_range = [data_min, data_max]

        # add 1:1 calibration line
        self.add_scattergl(
            x=cal_range,
            y=cal_range,
            mode='lines',
            name='1:1 line'
        )

    def set_layout(self):
        """
        Updates the layout configuration for a plot with specific titles for axes, legend,
        and the plot itself.
        :return: None
        """
        self.update_layout(
            title='Calibration',
            xaxis_title='Simulated',
            yaxis_title='Observed'
        )


if __name__ == '__main__':
    fig = CalibrationPlot()
