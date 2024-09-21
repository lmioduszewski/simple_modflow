import pandas as pd
import numpy as np
import pandas.core.indexes.datetimes
import plotly.graph_objects as go
import plotly.subplots
from plotly.subplots import make_subplots
from pathlib import Path
from simple_modflow.modflow.mf6.paths import *
from figs import Fig, Template
from plotly.colors import DEFAULT_PLOTLY_COLORS as colorsbo

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
                except:  # if error
                    print(f'{filename.stem} is not a valid .xlsx file! Skipping...')
                    continue  # skip to next file in for loop
                """create a dictionary of Pandas dataframes. Each item in the 
                dict is an imported excel file"""
                df_excel_dict[filename.stem] = thisdf  # store df in dict
            else:
                print(f'{filename.stem} is not a valid .xlsx file! Skipping...')
        return df_excel_dict

    @staticmethod
    def _resample_timeseries_df(dframe: pd.DataFrame = None, time_step='D', time_col_name=None, by='mean'):
        if time_col_name is None:
            time_col_name = dframe.columns[0]
        dframe = dframe.set_index(time_col_name)
        try:
            resampled_df = dframe.resample(time_step)
        except:
            return print('Error: Check that a valid resample rule (time_step) is set, and a valid datetime-like'
                         'time_col_name is set.')
        valid_bys = ['mean', 'sum', 'max', 'min']
        if by in valid_bys:
            return getattr(resampled_df, by)()
        else:
            return print(f'ValueError: "by" attribute must be one of {valid_bys}')

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
            return print('No Excel dict defined. Run read_excel() first.')
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


class ChoroplethPlot(Fig):

    def __init__(self, vor=None, zoom=13):
        super().__init__()

        if vor:
            mapbox_center = {"lat": vor.grid_centroid.y, "lon": vor.grid_centroid.x}
        else:
            mapbox_center = None
        self.update_layout(
            margin={"r": 0, "t": 20, "l": 0, "b": 0},
            mapbox_style="carto-positron",
            mapbox_zoom=zoom,
            mapbox_center=mapbox_center
        )


if __name__ == "__main__":
    date_range = pd.date_range('2021.10.1', '2022.04.01')
    start_day = pd.to_datetime('2021.10.01')
    end_day = pd.to_datetime('2022.04.01')
    excel_path = Path().home().joinpath('mf6', 'CDF_Mounding', 'xls', 'wls')
    precip_path = Path().home().joinpath('mf6', 'CDF_Mounding', 'xls', 'Precip_daily2.xlsx')


    def fig_D_8():
        plt = WaterLevelPlot()
        df = plt.read_excel_files_in_dir(excel_path)
        top_Qvt_dict = pd.read_excel(path_wellSurvey).set_index('Unnamed: 0', drop=True).loc[
            'Top of Qvr/Qvt Elevation ( Feet )'].to_dict()
        actual_df = df[list(df.keys())[0]].copy()
        for well, qvt_contact in top_Qvt_dict.items():
            try:
                actual_df.loc[actual_df[well] <= qvt_contact, well] = qvt_contact
            except:
                continue
        df[list(df.keys())[0]] = actual_df
        plt.plot_from_excel_dict(resample_time_step='D', by='mean', group_legend=True, vary_dash_by_df=True)
        precip_plt = plt.add_precip(precip_path,
                                    plot_title='Figure D-8 | Simulated vs. Actual Water Levels - 2021 - 2022 Wet Season')
        precip_plt.show(config={'scrollZoom': True})


    def fig_D_9():

        excel_path = Path().home().joinpath('mf6', 'CDF_Mounding', 'xls', 'ts_precip')
        sep_plt = WaterLevelPlot()
        sep_df = sep_plt.read_excel_files_in_dir(excel_path)
        fig_num = 9
        for filename, df in sep_df.items():
            sep_plt_per_file = WaterLevelPlot()
            df = df.set_index(df.columns[0])
            names = df.columns  # obs names for this file
            column_nums = range(len(names))  # num of obs in this file
            trace_colors_dict = sep_plt_per_file._get_colors_for_traces(names)
            """for each dataframe add the x and y values for each column"""
            for column_num in column_nums:
                x = df.index  # all times in the first column
                y = df.iloc[:, column_num]  # y vals for this obs
                name = names[column_num]
                """add data to figure"""
                sep_plt_per_file.add_scattergl(
                    x=x, y=y, name=name,
                    mode='lines', marker_size=3,
                    line_color=trace_colors_dict[name],
                )
            sep_plt_per_file.update_yaxes(title_text='Separation from 3.75-foot Trigger Elevation (feet)',
                                          range=[-20, 2])
            sep_plt_title = go.layout.Title(
                text=f'Figure D-{fig_num} | {filename[4:]}',
                font_size=30, x=0.5,
                xanchor='center')
            sep_plt_per_file.update_layout(
                title=sep_plt_title,
            )
            fig_num += 1
            for trace in sep_plt_per_file.data:
                trace.x = date_range
            sep_plt_per_file.show()
            sep_path = Path().home().joinpath('mf6', 'CDF_Mounding', 'html_figs', f'{filename}.jpg')
            #plotly.io.write_image(sep_plt_per_file, sep_path, format='jpg', scale=6, width=1000, height=750)


    fig_D_8()
