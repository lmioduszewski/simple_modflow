from pathlib import Path
import plotly.graph_objs as go
import pandas as pd
import scipy
import numpy as np

data = pd.read_excel(Path.home().joinpath('Python', 'data', 'EB-48W.xlsx'))


def recovery_extrapolate(x: list, y: list, time_off: int | float, y0=None, step=1):
    """
    Function to extrapolate the drawdown curve of an aquifer test based on recovery data
    :param x: a list of the time data (x-values of plot)
    :param y: a list of the water level data for the full test including recovery. Should be elevation data.
    :param time_off: the time value where the pump was shut off and recovery started
    :param y0: optional parameter to define the starting water level. If not provided, the function will assume that the lowest
    :param step: time interval for which to interpolate values, defaults to 1
    :return:
    """
    if y0 is None:
        y0 = min(y)
    interpolated_data = scipy.interpolate.interp1d(x, y)
    x_extrapolated = [time for time in range(int(time_off), int(time_off) * 2, step)]  # need to extend past two times test length if recovery extends that long
    x_since_off = np.array(x_extrapolated) - int(time_off)
    sd_recovery = interpolated_data(x_extrapolated) - y0  # why y0?
    sd_since_off = interpolated_data(x_since_off) - y0
    sd_extrapolated = sd_recovery + sd_since_off
    y_extrapolated = sd_extrapolated + y0
    return x_extrapolated, y_extrapolated, sd_extrapolated


x, y, sd = recovery_extrapolate(
    x=data.loc[:, 'min'],
    y=data.loc[:, 'dd'],
    time_off=4324.97,
)
fig = go.Figure()
fig.add_scattergl(
    x=x,
    y=y
)
fig.add_scattergl(
    x=data.loc[:, 'min'],
    y=data.loc[:, 'dd'],
)
fig.update_layout(
    xaxis_type='log',
    yaxis_type='log'
)
fig.show(renderer='browser')
