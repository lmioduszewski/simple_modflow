import zipfile
from pathlib import Path

import pandas as pd
import pandas.core.indexes.datetimes
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from myflopy.modflow.mf6.paths import *
from myflopy.viz import Fig, Template

from myflopy._logging import get_logger

logger = get_logger(__name__)

colors = ['rgb(31, 119, 180)', 'rgb(255, 127, 14)', 'rgb(44, 160, 44)', 'rgb(214, 39, 40)',
          'rgb(148, 103, 189)', 'rgb(140, 86, 75)', 'rgb(227, 119, 194)', 'rgb(127, 127, 127)',
          'rgb(188, 189, 34)', 'rgb(23, 190, 207)']
dashes = ['solid', 'dash', 'dashdot', 'dot', 'longdash', 'longdashdot']


class WaterLevelPlot(Fig):
    def __init__(self):
        super().__init__()

        self._df_excel_dict = None
        self._trace_colors_dict = {}
        self._colors_already_assigned = []
        self._dash_dict = {}
        self._precip_fig = None

    def read_excel_files_in_dir(self, data_path: Path = None) -> dict:
        """
        Iterate through a directory (specified as data_path) and import *.xlsx Excel found there.
        Returns a dictionary where keys are the filename stems and the values are the DataFrames for
        each Excel file.
        :param data_path: path to directory of Excel files
        :return: Dict of DataFrames
        """
        self._df_excel_dict = {}
        df_excel_dict = self._df_excel_dict  # empty dict to store excel data
        """Interate through the data directory. Will only read .xlsx files"""
        for filename in data_path.iterdir():
            if filename.suffix == '.xlsx':
                try:  # try reading file
                    thisdf = pd.read_excel(filename)  # import excel data
                except (OSError, ValueError, zipfile.BadZipFile):
                    # A .xlsx that is unreadable (OSError), not actually a zip
                    # container (BadZipFile -- which subclasses Exception, not
                    # OSError), or malformed inside (ValueError). This loop
                    # scans a whole directory, so skipping is the job.
                    logger.warning('%s is not a readable .xlsx file; skipping it', filename)
                    continue  # skip to next file in for loop
                """create a dictionary of Pandas dataframes. Each item in the
                dict is an imported excel file"""
                df_excel_dict[filename.stem] = thisdf  # store df in dict
            else:
                logger.debug('%s is not a .xlsx file; skipping it', filename)
        return df_excel_dict

    @staticmethod
    def _resample_timeseries_df(dframe: pd.DataFrame = None, time_step='D', time_col_name=None, by='mean'):
        if time_col_name is None:
            time_col_name = dframe.columns[0]
        dframe = dframe.set_index(time_col_name)
        try:
            resampled_df = dframe.resample(time_step)
        except (ValueError, TypeError) as error:
            # ValueError: `time_step` is not a pandas offset alias. TypeError:
            # the index is not datetime-like. This used to `return print(...)`,
            # which printed the message and handed the caller None -- so the
            # failure surfaced later as an AttributeError on None, far from the
            # bad argument that caused it.
            raise ValueError(
                f"cannot resample by {time_step!r}: check that it is a valid "
                f"pandas resample rule and that {time_col_name!r} is a "
                f"datetime-like column."
            ) from error
        valid_bys = ['mean', 'sum', 'max', 'min']
        if by in valid_bys:
            return getattr(resampled_df, by)()
        else:
            raise ValueError(f'`by` must be one of {valid_bys}, not {by!r}.')

    def _get_colors_for_traces(self, names=None, color_list=colors) -> dict:
        """
        Helper method to get color names for each trace being added to a fig.
        Returns a dictionary of names and colors. Can be used to sync colors between traces
        with certain names.
        :param names: iterable with names of traces
        :param color_list: list of eligible colors. Defaults to default plotly color list.
        :return: Dict where keys are names and values are CSS colors
        """
        if names is None:
            return self._trace_colors_dict
        trace_colors_dict = self._trace_colors_dict
        for idx, name in enumerate(names):
            if name in trace_colors_dict.keys():
                continue
            color = color_list[idx % len(color_list)]
            trace_colors_dict[name] = color
        return trace_colors_dict

    def _get_dashes_for_traces(self, names=None, dash_list=dashes):
        """Helper method to assign dash types to a list of list. Returns a dict
        of name keys and dash type values."""
        if not names:
            return self._dash_dict
        dash_dict = self._dash_dict
        for idx, name in enumerate(names):
            dash_type = dash_list[idx % len(dash_list)]
            dash_dict[name] = dash_type
        return dash_dict

    def plot_from_excel_dict(
            self,
            resample_time_step: str = None,
            time_col_name: str = None,
            by: str = 'mean',
            group_legend: bool = True,
            vary_dash_by_df=False
    ):
        if self._df_excel_dict is None:
            raise RuntimeError(
                'no Excel data has been read yet; call read_excel_files_in_dir() '
                'before plotting.'
            )
        else:
            df_excel_dict = self._df_excel_dict
        if vary_dash_by_df:
            dash_dict = self._get_dashes_for_traces(names=df_excel_dict.keys())
        for filename, df in df_excel_dict.items():
            if time_col_name:
                time_col_check = isinstance(pd.Index(df.loc[:, time_col_name]), pd.DatetimeIndex)
            else:
                time_col_check = isinstance(pd.Index(df.iloc[:, 0]), pd.DatetimeIndex)
            if time_col_check:
                if resample_time_step:
                    df = self._resample_timeseries_df(
                        df, time_step=resample_time_step,
                        time_col_name=time_col_name,
                        by=by
                    )
                elif time_col_name:
                    df = df.set_index(time_col_name)
                else:
                    df = df.set_index(df.columns[0])
            else:
                df = df.set_index(df.columns[0])

            if group_legend:
                legendgroup = filename
            else:
                legendgroup = None

            names = df.columns  # obs names for this file
            column_nums = range(len(names))  # num of obs in this file
            trace_colors_dict = self._get_colors_for_traces(names)
            """for each dataframe add the x and y values for each column"""
            dash_type = dash_dict[filename] if vary_dash_by_df else 'solid'
            for column_num in column_nums:
                x = df.index  # all times in the first column
                y = df.iloc[:, column_num]  # y vals for this obs
                name = names[column_num]
                """add data to figure"""
                self.add_scattergl(
                    x=x, y=y, name=name,
                    mode='lines', marker_size=3,
                    legendgroup=legendgroup,
                    legendgrouptitle_text=legendgroup,
                    line_color=trace_colors_dict[name],
                    line_dash=dash_type
                )
        return self.show()

    def add_precip(
            self,
            precip_path: Path = None,
            datetime_column_idx: int = 0,
            precip_column_idx: int = 1,
            y_title: str = None,
            x_title: str = None,
            precip_units: str = 'Inches',
            legend_title: str = 'Precipitation',
            trace_color: str = 'blue',
            plot_title: str = None
    ) -> go.Figure:
        """
        Adds a precipitation subplot and returns a new figure with the same data as the old figure.
        :param precip_path: Path of precipitation file
        :param datetime_column_idx: Integer index of column with datetime data (defaults to 0).
        :param precip_column_idx: Integer index of column with precip data (defaults to 1).
        :param y_title:
        :param x_title:
        :param precip_units: str: precip units
        :param legend_title: str: legend title
        :param trace_color: color of precip Bar trace
        :param plot_title: title of plot
        :return: go.Figure
        """
        precip_df = pd.read_excel(precip_path)
        precip_df = precip_df.set_index(precip_df.columns[datetime_column_idx])
        self._precip_fig = make_subplots(
            rows=2, cols=1, column_widths=[1],
            row_heights=[0.7, 0.3], vertical_spacing=0.06,
            shared_xaxes=True, x_title=x_title,
            y_title=y_title,
        )
        precip_fig = self._precip_fig
        precip_fig.add_traces(
            list(self.data),
            rows=[1 for i in range(len(self.data))],
            cols=[1 for i in range(len(self.data))]
        )
        precip_fig.add_trace(
            go.Bar(
                x=precip_df.index,
                y=precip_df.iloc[:, precip_column_idx - 1],
                name=precip_units,
                legendgroup=legend_title,
                legendgrouptitle_text=legend_title,
                legendgrouptitle_font_size=20,
                marker_color=trace_color
            ), row=2, col=1
        )
        precip_fig.update_layout(Fig().layout)
        precip_fig.update_yaxes(title_text=f'{legend_title} ({precip_units})', row=2, col=1)
        precip_fig.update_yaxes(title_text='Groundwater Elevation (feet)', row=1, col=1)
        precip_fig.update_xaxes(Template()._xaxis_template)
        precip_fig.update_layout(
            title=go.layout.Title(text=plot_title, font_size=30, x=0.5, xanchor='center'))
        return precip_fig
