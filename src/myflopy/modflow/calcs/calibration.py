from __future__ import annotations
from typing import TYPE_CHECKING

from myflopy.modflow.mf6.simulation.base import SimulationBase

if TYPE_CHECKING:
    import plotly.graph_objects as go
    import pandas as pd

import numpy as np
import figs as f
from pathlib import Path
from pandas import IndexSlice as idxx
import pandas as pd
import itertools
import traceback


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

    observed = np.array(observed, dtype=float).flatten()
    simulated = np.array(simulated, dtype=float).flatten()

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


def _require_compare_columns(frame: pd.DataFrame, required: set[str], *, caller: str) -> pd.DataFrame:
    """Validate that a compare/residual frame includes the required columns."""

    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{caller} requires columns: {', '.join(sorted(missing))}")
    return frame.copy()


def _sort_compare_by_time(frame: pd.DataFrame) -> pd.DataFrame:
    """Sort a compare-style DataFrame by time/per when possible."""

    if "time" not in frame.columns:
        return frame.reset_index(drop=True)
    sort_values = pd.to_numeric(frame["time"], errors="coerce")
    if sort_values.notna().all():
        return frame.assign(_sort=sort_values).sort_values("_sort").drop(columns="_sort").reset_index(drop=True)
    return frame.sort_values("time").reset_index(drop=True)


def _multiindex_head_frames(compare: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build observed/simulated MultiIndex head frames from compare output."""

    frame = _require_compare_columns(
        compare,
        {"name", "layer", "head_target", "sim_head"},
        caller="CalibrationPlot.from_compare(type='heads')",
    )
    time_values = pd.to_numeric(frame["per"] if "per" in frame.columns else frame["time"], errors="coerce")
    if time_values.isna().any():
        raise ValueError(
            "CalibrationPlot.from_compare(type='heads') requires numeric per/time values."
        )
    tuples = list(
        zip(
            frame["name"].astype(str),
            pd.to_numeric(frame["layer"], errors="coerce").fillna(0).astype(int),
            [(0, int(value)) for value in time_values],
        )
    )
    index = pd.MultiIndex.from_tuples(tuples, names=["locs", "layer", "kstpkper"])
    observed = pd.DataFrame({"elev": frame["head_target"].to_numpy()}, index=index)
    simulated = pd.DataFrame({"elev": frame["sim_head"].to_numpy()}, index=index)
    return observed, simulated


def _infer_compare_value_columns(
    frame: pd.DataFrame,
    *,
    target_column: str | None = None,
    simulated_column: str | None = None,
    caller: str,
) -> tuple[str, str]:
    """Resolve target/simulated value columns for compare-style frames."""

    if target_column is not None and simulated_column is not None:
        _require_compare_columns(
            frame,
            {target_column, simulated_column},
            caller=caller,
        )
        return target_column, simulated_column

    known_pairs = [
        ("head_target", "sim_head"),
        ("stage_target", "sim_stage"),
        ("flow_target", "sim_flow"),
        ("target_value", "sim_value"),
        ("value_target", "sim_value"),
    ]
    for candidate_target, candidate_simulated in known_pairs:
        if {candidate_target, candidate_simulated}.issubset(frame.columns):
            return candidate_target, candidate_simulated

    raise ValueError(
        f"{caller} could not infer target/simulated value columns. "
        "Pass target_column= and simulated_column= explicitly."
    )


class CalibrationPlot(f.Fig):

    @classmethod
    def from_compare(
        cls,
        compare: pd.DataFrame,
        *,
        type: str = "calibration",
        target_column: str | None = None,
        simulated_column: str | None = None,
        title: str | None = None,
        yaxis_title: str | None = None,
    ):
        """Build a calibration plot directly from a compare/residual DataFrame."""

        resolved_target, resolved_simulated = _infer_compare_value_columns(
            compare,
            target_column=target_column,
            simulated_column=simulated_column,
            caller="CalibrationPlot.from_compare(...)",
        )
        frame = compare.copy()

        if type == "heads":
            if {resolved_target, resolved_simulated} != {"head_target", "sim_head"}:
                frame = frame.rename(columns={resolved_target: "head_target", resolved_simulated: "sim_head"})
            observed, simulated = _multiindex_head_frames(frame)
            return cls(observed=observed, simulated=simulated, type=type)

        if type == "obs_vs_sim":
            return cls.from_obs_vs_sim(
                frame,
                target_column=resolved_target,
                simulated_column=resolved_simulated,
                title=title or "Observed vs simulated",
                yaxis_title=yaxis_title or "Simulated",
            )

        if type == "residuals_by_period":
            return cls.from_residuals_by_period(frame)

        return cls.from_obs_vs_sim(
            frame,
            target_column=resolved_target,
            simulated_column=resolved_simulated,
            title=title or "Calibration",
            yaxis_title=yaxis_title or "Simulated",
            add_stats=True,
        )

    @classmethod
    def from_targets(cls, targets, *, model=None, type: str = "calibration"):
        """Build a calibration plot from ``HeadTargets`` or a bound target helper."""

        if hasattr(targets, "targets") and model is None:
            compare = targets.compare()
        else:
            if model is None:
                raise ValueError(
                    "CalibrationPlot.from_targets(...) requires model= when given raw HeadTargets."
                )
            compare = targets.compare(model)
        return cls.from_compare(compare, type=type)

    @classmethod
    def from_obs_vs_sim(
        cls,
        compare: pd.DataFrame,
        *,
        baseline_compare: pd.DataFrame | None = None,
        target_column: str = "head_target",
        simulated_column: str = "sim_head",
        title: str = "Observed vs simulated heads",
        xaxis_title: str = "Observed",
        yaxis_title: str = "Simulated",
        add_stats: bool = False,
    ):
        """Build a target-vs-simulated cross plot from compare-style DataFrames."""

        frame = _require_compare_columns(
            compare,
            {target_column, simulated_column},
            caller="CalibrationPlot.from_obs_vs_sim(...)",
        ).dropna(subset=[target_column, simulated_column])
        fig = cls()
        fig.add_scattergl(
            x=frame[target_column],
            y=frame[simulated_column],
            mode="markers",
            name="Model",
        )
        values = [frame[target_column], frame[simulated_column]]
        if baseline_compare is not None:
            baseline = _require_compare_columns(
                baseline_compare,
                {target_column, simulated_column},
                caller="CalibrationPlot.from_obs_vs_sim(..., baseline_compare=...)",
            ).dropna(subset=[target_column, simulated_column])
            fig.add_scattergl(
                x=baseline[target_column],
                y=baseline[simulated_column],
                mode="markers",
                name="Baseline",
            )
            values.extend([baseline[target_column], baseline[simulated_column]])

        merged_values = pd.concat(values, axis=0).dropna()
        if not merged_values.empty:
            lower = float(merged_values.min())
            upper = float(merged_values.max())
            fig.add_scattergl(
                x=[lower, upper],
                y=[lower, upper],
                mode="lines",
                name="1:1 line",
            )
        fig.update_layout(
            title=title,
            xaxis_title=xaxis_title,
            yaxis_title=yaxis_title,
        )
        if add_stats and not frame.empty:
            fig.add_stats(observed=frame[target_column], simulated=frame[simulated_column])
        return fig

    @classmethod
    def from_timeseries(
        cls,
        compare: pd.DataFrame,
        *,
        name: str | None = None,
        baseline_compare: pd.DataFrame | None = None,
        target_column: str = "head_target",
        simulated_column: str = "sim_head",
        target_label: str = "Target",
        simulated_label: str = "Model",
        baseline_label: str = "Baseline",
        yaxis_title: str = "Head",
        title: str | None = None,
    ):
        """Build a time-series plot of targets and simulated heads for one target."""

        frame = _require_compare_columns(
            compare,
            {"name", "time", target_column, simulated_column},
            caller="CalibrationPlot.from_timeseries(...)",
        )
        names = frame["name"].astype(str)
        if name is None:
            unique_names = sorted(names.unique().tolist())
            if len(unique_names) != 1:
                raise ValueError(
                    "CalibrationPlot.from_timeseries(...) requires name= when more than one target location is present."
                )
            name = unique_names[0]
        key = str(name).strip().lower()
        frame = frame.loc[names.str.lower() == key].copy()
        if frame.empty:
            raise ValueError(f"No compare rows found for observation name {name!r}.")
        frame = _sort_compare_by_time(frame)

        fig = cls()
        fig.add_scattergl(
            x=frame["time"],
            y=frame[target_column],
            mode="lines+markers",
            name=target_label,
        )
        fig.add_scattergl(
            x=frame["time"],
            y=frame[simulated_column],
            mode="lines+markers",
            name=simulated_label,
        )
        if baseline_compare is not None:
            baseline = _require_compare_columns(
                baseline_compare,
                {"name", "time", simulated_column},
                caller="CalibrationPlot.from_timeseries(..., baseline_compare=...)",
            )
            baseline_names = baseline["name"].astype(str)
            baseline = baseline.loc[baseline_names.str.lower() == key].copy()
            baseline = _sort_compare_by_time(baseline)
            fig.add_scattergl(
                x=baseline["time"],
                y=baseline[simulated_column],
                mode="lines+markers",
                name=baseline_label,
            )
        fig.update_layout(
            title=title or str(frame["name"].iloc[0]),
            xaxis_title="Time / period",
            yaxis_title=yaxis_title,
        )
        return fig

    @classmethod
    def from_residuals_by_period(
        cls,
        compare: pd.DataFrame,
        *,
        baseline_compare: pd.DataFrame | None = None,
        title: str = "Residual MAE by period",
        yaxis_title: str = "Mean absolute error",
    ):
        """Build a by-period residual summary plot from compare-style data."""

        frame = _require_compare_columns(
            compare,
            {"time", "abs_residual"},
            caller="CalibrationPlot.from_residuals_by_period(...)",
        )
        current_stats = (
            frame.groupby("time")
            .agg(mae=("abs_residual", "mean"))
            .reset_index()
        )
        current_stats = _sort_compare_by_time(current_stats.rename(columns={"mae": "value"}))

        fig = cls()
        fig.add_scattergl(
            x=current_stats["time"],
            y=current_stats["value"],
            mode="lines+markers",
            name="Model MAE",
        )

        if baseline_compare is not None:
            baseline = _require_compare_columns(
                baseline_compare,
                {"time", "abs_residual"},
                caller="CalibrationPlot.from_residuals_by_period(..., baseline_compare=...)",
            )
            baseline_stats = (
                baseline.groupby("time")
                .agg(mae=("abs_residual", "mean"))
                .reset_index()
            )
            baseline_stats = _sort_compare_by_time(baseline_stats.rename(columns={"mae": "value"}))
            fig.add_scattergl(
                x=baseline_stats["time"],
                y=baseline_stats["value"],
                mode="lines+markers",
                name="Baseline MAE",
            )
        fig.update_layout(
            title=title,
            xaxis_title="Time / period",
            yaxis_title=yaxis_title,
        )
        return fig

    def __init__(
            self,
            observed: pd.Series | list | Path = None,
            simulated: pd.Series | list = None,
            model: SimulationBase = None,
            loc_cell_dict: dict = None,
            loc_shp_gpkg: Path = None,
            lak_obs_dict: dict = None,
            cal_range: list | tuple = None,
            exploration_name_field: str = 'ExploName',
            crs=None,
            obs_layers: int | list[int] = None,
            verbose: bool = False,
            type: str = 'calibration'
    ):
        """
        Initializes an object to plot observed and simulated calibration data.

        :param observed: Observed data, provided as a pandas Series, a list, or a Path
            object pointing to the data.
        :param simulated: Simulated model data, provided as a pandas Series or a list.
            if not provided, simulated data is generated using the model object and
            and spatial data provided in loc_shp_gpkg.
        :param model: Simulation model object adhering to the SimulationBase class.
        :param loc_cell_dict: Dictionary linking geospatial locations to cell mappings
            or identifiers.
        :param loc_shp_gpkg: Path to shapefile or GeoPackage used for geospatial
            location of point observations of heads (i.e. wells)
        :param lak_obs_dict: Dictionary linking lake names to lake number ids.
        :param cal_range: List or tuple specifying the calibration range, typically as
            time steps or indices for calibration analysis.
        :param exploration_name_field: Field name in the geospatial file (e.g.,
            shapefile or GeoPackage) representing exploration or location names.
        :param crs: Coordinate reference system used for geospatial data processing
            and mapping.
        :param obs_layers: Integer or list of integers representing observation layers
        :param type: can be 'calibration' for cross plot or 'heads' to show actual heads as a plot.
        """
        super().__init__()
        self._observed = None
        self._simulated = None
        self._model = None
        self._crs = crs
        self._cal_range = cal_range

        self._loc_cell_dict = loc_cell_dict
        self._loc_shp_gpkg = None
        self._lak_obs_dict = None
        self._exploration_name_field = exploration_name_field
        self._obs_layers = obs_layers
        self._verbose = verbose

        self.model = model
        self.lak_obs_dict = lak_obs_dict
        self.loc_shp_gpkg = loc_shp_gpkg
        self.simulated = simulated
        self.observed = observed

        if self.observed is not None and self.simulated is not None:
            if type == 'calibration':
                self.add_calib_scatter()
                self.add_stats()
                self.set_layout()
            if type == 'heads':
                self.add_heads()


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
        """gets simulated data from model object if no simulated data is provided"""
        if self._simulated is None and self.loc_shp_gpkg and self.model:

            if self._verbose:
                print('generating simulated data using model object and spatial data provided')
            simulated = self.model.hds.get_obs_heads(locs=self.loc_shp_gpkg, long_format=True)

            # add lake stage data to simulated data
            if self.lak_obs_dict is not None:
                lak_obs = {}
                # get lake stage data for each lake
                for k, v in self.lak_obs_dict.items():
                    lak_obs[k] = self.model.outputs.lak.stage.get()[:, v].tolist()
                print(len(lak_obs['Deep Lake']))
                # create a multi-index for lake stage data
                lake_names = list(lak_obs.keys())
                layers = [0]  # set layer to 0 for all lake observations
                kstpkper = self.model.gwf.output.head().get_kstpkper()
                assert len(lak_obs[lake_names[0]]) == len(kstpkper), \
                    f'length of lake stage data does not match number of stress periods in model'
                lak_idx = pd.MultiIndex.from_product(
                    [lake_names, layers, kstpkper], names=['locs', 'layer', 'kstpkper'])
                lake_data = list(itertools.chain.from_iterable(lak_obs.values()))
                # create a DataFrame with lake stage data and append to simulated data
                lak_df = pd.DataFrame(lake_data, index=lak_idx, columns=['elev'])
                simulated = pd.concat([simulated, lak_df], axis=0)

            simulated = simulated[sorted(simulated.columns)]
            self._simulated = simulated

        return self._simulated

    @simulated.setter
    def simulated(self, value):
        if self.loc_shp_gpkg is not None:
            if self._verbose:
                print('path to shapefile or geopackage provided, '
                      'so provided simulated data is ignored, and'
                      'will generate observed data using model'
                      'object and spatial data provided')
        else:
            self._simulated = value

    @property
    def lak_obs_dict(self):
        return self._lak_obs_dict

    @lak_obs_dict.setter
    def lak_obs_dict(self, value):
        if value is not None:
            assert isinstance(value, dict), 'lak_obs_dict must be a dictionary'
            for lak_id in value.values():
                assert lak_id in range(self.model.outputs.lak.stage.nlakes), \
                    'lak_id not found in lak package'
        self._lak_obs_dict = value

    @property
    def observed(self):
        return self._observed

    @observed.setter
    def observed(self, value):
        if isinstance(value, Path):
            obs_data = pd.read_excel(value)
            assert len(obs_data) == self.model.nper, \
                'length of observed data does not match number of stress periods in model'
            if self._obs_layers is None:
                self._obs_layers = 0  # default to layer 1
            col1 = obs_data.columns[0]  # should be stress periods
            if len(obs_data) != len(self.model.kstpkper):
                print(f'model kstpkper {len(self.model.kstpkper)} does not match length of observed data,'
                      f'{len(obs_data)}, trimming observed data to match model kstpkper')
                obs_data = obs_data.iloc[:len(self.model.kstpkper), :]
            obs_data['kstpkper'] = self.model.kstpkper
            # drop stress period column in favor of model kstpkper to match simulated data
            obs_data.drop(columns=col1, inplace=True)
            obs_data = obs_data[sorted(obs_data.columns)]
            obs_data['layer'] = self._obs_layers
            obs_data = obs_data.set_index(['kstpkper', 'layer'])
            obs_data = obs_data.melt(ignore_index=False, value_name='elev', var_name='locs')
            obs_data = obs_data.reset_index().set_index(['locs', 'layer', 'kstpkper'])
            value = obs_data
            """except ValueError as e:
                print('observed data path not readable, must be Excel file')
                traceback.print_exc()"""
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
            sim_obs = pd.concat([simulated, observed], axis=1).dropna()
            sim_obs.columns = ['simulated', 'observed']
            # drop nan rows so the calib stats are calculated correctly
            sim_obs.dropna()
            simulated = sim_obs.loc[:, 'simulated'].to_numpy()
            observed = sim_obs.loc[:, 'observed'].to_numpy()

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
            if self._verbose:
                print('using observation data based on provided'
                      'shapefile or geopackage for calibration plot')
            assert isinstance(self.simulated, pd.DataFrame)
            loc_names = self.simulated.index.get_level_values(level=0).unique().to_list()
            layers = self.simulated.index.get_level_values(level=1).unique().to_list()

            for loc_name in loc_names:
                for layer in layers:

                    sim_data = pd.Series(self.simulated.loc[idxx[loc_name, layer, :], 'elev'].to_list())
                    obs_data = pd.Series(self.observed.loc[idxx[loc_name, layer, :], 'elev'].to_list())
                    sim_obs = pd.concat([sim_data, obs_data], axis=1)
                    sim_obs.columns = ['simulated', 'observed']
                    sim_obs = sim_obs.dropna()

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
                        x=sim_obs.loc[:, 'simulated'],
                        y=sim_obs.loc[:, 'observed'],
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

    def add_heads(self):

        locs = list(self.observed.index.get_level_values(0).unique())
        layers = list(self.observed.index.get_level_values(1).unique())

        for loc in locs:
            for layer in layers:

                # get observed data for location and layer
                loc_obs = self.observed.loc[idxx[loc, layer, :], :].dropna().reset_index().drop(
                    ['locs', 'layer'], axis=1)
                loc_obs.kstpkper = loc_obs.kstpkper.apply(lambda x: x[1])
                loc_obs.set_index('kstpkper', inplace=True)
                if len(loc_obs) == 0:
                    continue

                self.add_scattergl(
                    name=f'{loc} observed',
                    x=loc_obs.index,
                    y=loc_obs.elev,
                    line_color='blue'
                )

                # get simulated data for location and layer
                sim_obs = self.simulated.loc[idxx[loc, layer, :], :].dropna().reset_index().drop(
                    ['locs', 'layer'], axis=1)
                sim_obs.kstpkper = sim_obs.kstpkper.apply(lambda x: x[1])
                sim_obs.set_index('kstpkper', inplace=True)

                self.add_scattergl(
                    name=f'{loc} simulated',
                    x=sim_obs.index,
                    y=sim_obs.elev,
                    line_color='red'
                )



if __name__ == '__main__':
    fig = CalibrationPlot()
