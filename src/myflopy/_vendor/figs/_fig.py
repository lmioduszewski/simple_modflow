from __future__ import annotations

from pathlib import Path
import pandas as pd
import shapely as shp
from plotly import graph_objects as go
from plotly import express as px
from plotly.subplots import make_subplots
import geopandas as gpd
import plotly
import pickle
from . import layout_template as lt
import numpy as np
from ._layout_presets import LayoutPresets
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import MultipleLocator

data_dir = Path.cwd().joinpath('sample_data')


def _validate_dataframe(value, *, name: str):
    if not isinstance(value, pd.DataFrame):
        raise TypeError(f"{name} must be a pandas DataFrame, got {type(value)}")
    if value.empty:
        raise ValueError(f"{name} must not be empty")
    return value


def _coerce_name_collection(values, *, name: str):
    if values is None:
        return None
    if isinstance(values, str):
        return {values}
    if isinstance(values, (list, tuple, set)):
        return set(values)
    raise TypeError(f"{name} must be a string or a collection of strings")


class Template:
    """Base Template for figure objects in the Fig class"""

    def __init__(self):

        # FIG TEMPLATE
        self._config = {
            'scrollZoom': True,
            'displaylogo': False
        }
        self.default_trace_colors = plotly.colors.DEFAULT_PLOTLY_COLORS
        self._trace_colors_dict = {}
        self.layout = lt.layout
        self.layout_template = lt.template
        self.modebar = lt.modebar
        self.scattergl_template = lt.scattergl_template

    def _get_colors_for_traces(self, names=None, color_list=None) -> dict:
        """
        Helper method to get color names for each trace being added to a fig.
        Returns a dictionary of names and colors. Can be used to sync colors between traces
        with certain names.
        :param names: iterable with names of traces
        :param color_list: list of eligible colors. Defaults to default plotly color list.
        :return: Dict where keys are names and values are CSS colors
        """
        if color_list is None:
            color_list = self.default_trace_colors
        if names is None:
            return self._trace_colors_dict
        trace_colors_dict = self._trace_colors_dict
        for idx, name in enumerate(names):
            if name in trace_colors_dict.keys():
                continue
            color = color_list[idx % len(color_list)]
            trace_colors_dict[name] = color
        return trace_colors_dict

    @staticmethod
    def create_hover(name_dict: dict = None):
        """
        Function to return a hover template for a Plotly figure

        :args:
        name_dict [dict]: dict where the keys are the names of each hover data label,
        such as 'Proj #'. The values associated with each dict key are a list of data for that key,
        such as a list of project numbers, one for each data point in the figure trace

        :return: customdata, hovertemplate
        """
        names = list(name_dict.keys())
        lists = list(name_dict.values())
        list_len = len(lists[0])
        custom_data = []
        hover_template_list = []

        for i in range(list_len):
            this_list = [param[i] for param in lists]
            custom_data.append(this_list)

        for i, name in enumerate(names):
            if type(custom_data[0][i]) is float:
                hover_template_list.append(
                    f'<b>{name}: </b>%{{customdata[{i}]:.2f}}<br>'
                )
            else:
                hover_template_list.append(
                    f'<b>{name}: </b>%{{customdata[{i}]}}<br>'
                )
        hover_template_list.append('<extra></extra>')
        hover_template = ''.join(hover_template_list)

        return custom_data, hover_template

    @staticmethod
    def add_to_hover_dict(existing_hover_dict: dict, dict_to_add: dict):
        """add key,value pairs from a dict to another dict"""
        new_hover_dict = existing_hover_dict
        for hover_name, hover_data in dict_to_add.items():
            new_hover_dict[hover_name] = hover_data
        return new_hover_dict

    @staticmethod
    def add_aesi_logo(fig: go.Figure):
        """adds the AESI logo to the center bottom of the figure"""
        fig.add_layout_image(
            source=f'data:image/png;base64,{lt.aesi_logo}',
            layer='above',
            xref="paper", yref="paper",
            x=0.5, y=-0.15,
            sizex=0.1, sizey=0.1,
            xanchor="center", yanchor="bottom"
        )


template = Template()

# Keep this allowlist in sync with `PlotlyExpressProxy` in `src/figs/_fig.pyi`.
# The runtime proxy and the IDE stub intentionally mirror one another so users
# get the same Plotly Express surface in both execution and autocomplete.
_PX_METHODS = {
    "area",
    "bar",
    "bar_polar",
    "box",
    "choropleth",
    "choropleth_map",
    "choropleth_mapbox",
    "density_contour",
    "density_heatmap",
    "density_map",
    "density_mapbox",
    "ecdf",
    "funnel",
    "funnel_area",
    "histogram",
    "icicle",
    "imshow",
    "line",
    "line_3d",
    "line_geo",
    "line_map",
    "line_mapbox",
    "line_polar",
    "line_ternary",
    "parallel_categories",
    "parallel_coordinates",
    "pie",
    "scatter",
    "scatter_3d",
    "scatter_geo",
    "scatter_map",
    "scatter_mapbox",
    "scatter_matrix",
    "scatter_polar",
    "scatter_ternary",
    "strip",
    "sunburst",
    "timeline",
    "treemap",
    "violin",
}


class PlotlyExpressProxy:
    """Expose Plotly Express constructors while merging results into a Fig."""

    def __init__(self, fig: "Fig" = None, fig_cls: type["Fig"] = None):
        self._fig = fig
        self._fig_cls = fig_cls if fig_cls is not None else type(fig)

    def _get_target_fig(self) -> "Fig":
        if self._fig is not None:
            return self._fig
        if self._fig_cls is None:
            raise RuntimeError("PlotlyExpressProxy requires a Fig instance or Fig class")
        return self._fig_cls()

    def __getattr__(self, name: str):
        if name not in _PX_METHODS:
            raise AttributeError(f"figs.plotly.Fig.px has no Plotly Express method {name!r}")

        px_method = getattr(px, name, None)
        if px_method is None or not callable(px_method):
            raise AttributeError(f"plotly.express has no callable attribute {name!r}")

        def _wrapped(*args, **kwargs):
            fig = self._get_target_fig()
            return fig._add_px_figure(name, *args, **kwargs)

        _wrapped.__name__ = name
        _wrapped.__doc__ = getattr(px_method, "__doc__", None)
        return _wrapped

    def __dir__(self):
        return sorted(_PX_METHODS)


class _PlotlyExpressAccessor:
    """Descriptor exposing Plotly Express helpers on both Fig and Fig instances."""

    def __get__(self, instance: Fig | None, owner: type[Fig]) -> PlotlyExpressProxy:
        if instance is None:
            return PlotlyExpressProxy(fig_cls=owner)
        return instance._px_proxy


class Fig(go.Figure):
    """Use this to instantiate figures."""

    px = _PlotlyExpressAccessor()

    def __init__(
            self,
            add_logo=False,
            subplot: go.Figure = None,
            preset: str = None,
            *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._config = template._config
        self._post_scripts = []
        self._subplot = subplot
        self.update_layout(template=template.layout_template)
        if preset is not None:
            # sets a layout preset if provided
            self.update_layout(LayoutPresets.set(preset))
        self._titles = lt.titles
        self._trace_specs = {
            'precip_color': 'blue',
            'precip_width': 1,
            'water_levels_width': 1.5
        }
        self._px_proxy: PlotlyExpressProxy = PlotlyExpressProxy(fig=self)
        if self._subplot is not None:
            self.data = self._subplot.data
            self.layout = self._subplot.layout
            self._grid_ref = self._subplot._grid_ref
        if add_logo:
            template.add_aesi_logo(self)

    @property
    def trace_specs(self):
        return self._trace_specs

    def subplot(self, show_precip=False, show_map=False, show_flow=False,
                num_rows=None, num_cols=None, row_heights=None, col_widths=None, *args, **kwargs) -> "Fig":
        subplot = Subplot(show_precip=show_precip, show_map=show_map,
                          show_flow=show_flow, num_rows=num_rows,
                          row_heights=row_heights, col_widths=col_widths,
                          num_cols=num_cols, *args, **kwargs)._subplot

        return Fig(subplot=subplot)

    @property
    def titles(self):
        return self._titles

    @staticmethod
    def _sanitize_px_layout(layout_dict: dict, *, for_subplot: bool) -> dict:
        sanitized = dict(layout_dict)
        sanitized.pop("template", None)
        if for_subplot:
            sanitized.pop("annotations", None)
            for key in list(sanitized):
                if key.startswith("xaxis") or key.startswith("yaxis"):
                    sanitized.pop(key, None)
        return sanitized

    def _add_plotly_figure(
            self,
            figure: go.Figure,
            *,
            row: int = None,
            col: int = None,
            secondary_y: bool = None,
            include_layout: bool = True,
    ) -> "Fig":
        trace_target = {}
        if row is not None:
            trace_target["row"] = row
        if col is not None:
            trace_target["col"] = col
        if secondary_y is not None:
            trace_target["secondary_y"] = secondary_y

        for trace in figure.data:
            if trace_target:
                self.add_trace(trace, **trace_target)
            else:
                self.add_trace(trace)

        if include_layout:
            layout_json = self._sanitize_px_layout(
                figure.layout.to_plotly_json(),
                for_subplot=bool(trace_target),
            )
            if layout_json:
                self.update_layout(layout_json)

        return self

    def _add_px_figure(self, method: str, *args, **kwargs) -> "Fig":
        row = kwargs.pop("row", None)
        col = kwargs.pop("col", None)
        secondary_y = kwargs.pop("secondary_y", None)

        if method not in _PX_METHODS:
            raise AttributeError(f"figs.plotly.Fig.px has no Plotly Express method {method!r}")

        px_method = getattr(px, method, None)
        if px_method is None or not callable(px_method):
            raise AttributeError(f"plotly.express has no callable attribute {method!r}")

        px_figure = px_method(*args, **kwargs)
        return self._add_plotly_figure(
            px_figure,
            row=row,
            col=col,
            secondary_y=secondary_y,
        )

    def add_df(self, df: pd.DataFrame, method: str = "scatter", *args, **kwargs) -> "Fig":
        """
        Build a Plotly Express figure from a DataFrame and merge it into this Fig.

        Example:
            fig.add_df(df, method="line", x="Date", y="WaterLevel", color="Well")
        """
        df = _validate_dataframe(df, name="df")
        if not isinstance(method, str):
            raise TypeError(f"method must be a string, got {type(method)}")
        return self._add_px_figure(method, df, *args, **kwargs)

    def add_scattergl(self, **kwargs):
        # Overrides add_scattergl and applies the defined scattergl template first then adds all other args
        trace = template.scattergl_template.update(**kwargs)
        self.add_trace(trace)

    def add_to_subplot(self, subplot_fig: go.Figure, row: int, col: int):
        """Add go.Figure traces to a subplot figure

        Args:
            subplot_fig (make_subplots): subplot figure to add to
            row (int): row of subplot to add to
            col (int): column of subplot to add to
        """
        num_traces = range(len(self.data))
        subplot_fig.add_traces(
            list(self.data),
            rows=[row for i in num_traces],
            cols=[col for i in num_traces]
        )

    def add_water_levels(
            self,
            excel_paths: Path | list = None,
            df: pd.DataFrame = None,
            secondary_ys=None,
            row=1,
            col=1,
            **kwargs
    ):
        """
        Add water level data to the main water level subplot. Can provide excel paths or a single DataFrame. The excel
        paths take presidence over the df if provided. Default trace mode is a line, but you can specify marker or
        line by including 'HD' and 'DD' in the name of the excel worksheet. 'HD' = 'hand data' = markers. 'DD' =
        'data logger data' = lines.
        :param col: column in main subplot
        :param row: row of main subplot
        :param excel_paths: list of excel paths with water level data to plot
        :param df: alternatively can provide a pandas DataFrame of the water level data
        :param secondary_ys: list of column names to plot on secondary y axis:
        :param kwargs: keyword arguments passed to go.Scattergl constructor
        """
        secondary_ys = _coerce_name_collection(secondary_ys, name="secondary_ys")
        if excel_paths is not None:
            from ._datatypes import ExcelDateData

            excel_data = ExcelDateData(excel_paths=excel_paths)
            data = excel_data.dfs
            sheet_names = [v['sheet_name'] for v in excel_data.excel_dict.values()]
            line_marker = []
            for name in sheet_names:
                if 'HD' in str(name):
                    line_marker.append('markers')
                elif 'DD' in str(name):
                    line_marker.append('lines')
                else:
                    line_marker.append('lines')
        elif df is not None:
            df = _validate_dataframe(df, name="df")
            line_marker = ['lines']
            if pd.api.types.is_datetime64_any_dtype(df.index):
                data = [df.copy()]
            else:
                if len(df.columns) < 2:
                    raise ValueError("df must have at least two columns when using the first column as the x-axis")
                working_df = df.copy()
                working_df.set_index(working_df.columns[0], inplace=True)
                data = [working_df]
        else:
            raise ValueError('No excel_paths or df data provided')

        trace_names = []
        for df in data:
            names = df.columns.values.tolist()
            trace_names += names
        trace_colors = template._get_colors_for_traces(trace_names)
        for i, df in enumerate(data):
            line_marker_mode = line_marker[i]
            for column in df.columns:
                if secondary_ys and column in secondary_ys:
                    secondary_y = True
                else:
                    secondary_y = False
                self.add_trace(
                    go.Scattergl(
                        arg=template.scattergl_template,
                        x=df.index,
                        y=df[column],
                        name=column,
                        line_width=self.trace_specs['water_levels_width'],
                        line_color=trace_colors[column],
                        marker_color=trace_colors[column],
                        mode=line_marker_mode,
                        **kwargs
                    ),
                    row=row,
                    col=col,
                    secondary_y=secondary_y
                ),
        self.update_yaxes(
            title_text="Elevation (ft)",
            showticklabels=True,
            automargin=True,
        )

    def add_precip(self, df: pd.DataFrame = None, row=2, col=1, cols_to_plot=None, type='bars', **kwargs):
        df = _validate_dataframe(df, name="df")
        if len(df.columns) < 2:
            raise ValueError("df must have at least two columns for precipitation plotting")
        if type not in {'bars', 'lines'}:
            raise ValueError("type must be either 'bars' or 'lines'")
        columns = None
        if cols_to_plot is None:
            columns = df.columns[1:]
        else:
            columns = df.columns[cols_to_plot:cols_to_plot + 1]
        if len(columns) > 1:
            self.update_layout(barmode='group')
        for loc in columns:
            if type == 'bars':
                self.add_trace(
                    go.Bar(
                        x=list(df.iloc[:, 0]),
                        y=list(df.loc[:, loc]),
                        marker_color=self._trace_specs['precip_color'],
                        marker_line_color=self._trace_specs['precip_color'],
                        name=loc,
                        **kwargs
                    ),
                    row=row,
                    col=col
                )
            if type == 'lines':
                self.add_trace(
                    go.Scattergl(
                        x=list(df.iloc[:, 0]),
                        y=list(df.loc[:, loc]),
                        name=loc,
                        mode='lines',
                        line_color=self._trace_specs['precip_color'],
                        line_width=self._trace_specs['precip_width'],
                        **kwargs
                    ),
                    row=row,
                    col=col
                )
        self.update_yaxes(
            title_text="Rainfall (in)",
            showticklabels=True,
            automargin=True,
            row=row,
            col=col
        )

    def add_map(self, locs: Path, loc_name_field='ExploName', row=1, col=2, map_zoom=13):
        gdf_locs = gpd.read_file(locs)
        map_center = shp.MultiPoint(gdf_locs.geometry).centroid
        map_center = dict(
            lat=map_center.y,
            lon=map_center.x
        )
        self.add_trace(
            go.Scattermap(
                lat=gdf_locs.geometry.y.to_list(),
                lon=gdf_locs.geometry.x.to_list(),
                mode='markers+text',
                textposition="top right",
                textfont=dict(size=12, color='black'),
                text=gdf_locs.loc[:, loc_name_field].to_list(),
                hovertemplate=gdf_locs.loc[:, loc_name_field] +
                              '<br>Lat: %{lat:.4f}<br>' +
                              'Lon: %{lon:.4f}<extra></extra>',
                showlegend=False,
                marker_color='blue'
            ),
            row=row,
            col=col
        )
        self.update_layout(
            xaxis=dict(domain=[0.0, 0.7]),
            xaxis2=dict(domain=[0.0, 0.7]),
            map=dict(
                style='open-street-map',
                center=map_center,
                zoom=map_zoom,
                domain=dict(x=[0.7, 1.0])
            )
        )

    def add_post_script(self, script: str) -> "Fig":
        """Attach a JavaScript snippet to run after this figure renders.

        Plotly runs each snippet inside the ``Plotly.newPlot(...).then(...)``
        callback of ``show()`` / ``write_html()``; use the ``{plot_id}``
        placeholder to reach the plot div. Snippets are ignored on inline
        (notebook mimebundle) display, which has no post-render hook -- the
        figure's static state still shows, only the JS behavior is skipped.
        """

        if script:
            scripts = getattr(self, "_post_scripts", None)
            if scripts is None:
                scripts = []
                self._post_scripts = scripts
            scripts.append(str(script))
        return self

    def _combined_post_script(self):
        scripts = getattr(self, "_post_scripts", None)
        return "\n".join(scripts) if scripts else None

    def show(self, renderer='browser', config=None, *args, **kwargs):
        if config is None:
            config = template._config
        post_script = self._combined_post_script()
        if post_script is not None and 'post_script' not in kwargs:
            kwargs['post_script'] = post_script
        super().show(renderer=renderer, config=config, *args, **kwargs)

    def write_html(self, file=None, config=None, *args, **kwargs):
        if config is None:
            config = self._config
        post_script = self._combined_post_script()
        if post_script is not None and 'post_script' not in kwargs:
            kwargs['post_script'] = post_script
        return super().write_html(file, config=config, *args, **kwargs)

    def _repr_mimebundle_(self, *args, **kwargs):
        """Inline display (Jupyter, last-expression) carries this Fig's config.

        Plotly's mimetype bundle has a ``config`` key the JupyterLab/notebook
        renderer applies via ``Plotly.newPlot``; the base figure doesn't put a
        figure-specific config there, so scroll-zoom (and ``displaylogo``) only
        engaged on ``show()`` / ``write_html()``. Merge ``self._config`` into the
        bundle so the interaction defaults work on inline display too.
        """
        bundle = super()._repr_mimebundle_(*args, **kwargs) or {}
        key = "application/vnd.plotly.v1+json"
        spec = bundle.get(key)
        if isinstance(spec, dict):
            merged = dict(spec)
            merged["config"] = {**merged.get("config", {}), **self._config}
            bundle = {**bundle, key: merged}
        return bundle


class Subplot:
    """
    Class for making subplots for water level plots.
    """

    def __init__(
            self,
            show_precip=False,
            show_map=False,
            show_flow=False,
            num_rows: int = None,
            num_cols: int = None,
            row_heights: list = None,
            col_widths: list = None,
            *args,
            **kwargs
    ) -> go.Figure:

        self._subplot_titles = []
        self._titles = lt.titles
        self._trace_specs = {
            'precip_color': 'blue',
            'precip_width': 1,
            'water_levels_width': 1.5
        }
        self._num_rows = None
        self._num_cols = None
        self._col_widths = None
        self._row_heights = None
        self._specs = None
        self._show_precip = show_precip
        self._show_map = show_map
        self._show_flow = show_flow

        self.num_rows = num_rows
        self.num_cols = num_cols
        if row_heights is not None:
            self.row_heights = row_heights
        if col_widths is not None:
            self.col_widths = col_widths

        self._subplot = self.make(
            rows=self.num_rows,
            cols=self.num_cols,
            specs=self.specs,
            column_widths=self.col_widths,
            row_heights=self.row_heights,
            *args, **kwargs
        )
        self.apply_template_layout()

    @property
    def num_rows(self):
        """Defines row heights based on the number of rows. The number of rows is determined by what the figure
        will show - water levels, precipitation, and flows. Up to three rows if all are shown."""
        if self._num_rows is None:
            if self._show_precip and self._show_flow:
                rows = 3
                self._row_heights = [0.7, 0.1, 0.2]
                titles = ['water levels', 'rain', 'flow']
            elif self._show_precip or self._show_flow:
                rows = 2
                self._row_heights = [0.75, 0.25]
                titles = ['water levels']
                if self._show_precip:
                    titles.append('rain')
                else:
                    titles.append('flow')
            else:
                rows = 1
                self._row_heights = [1]
                titles = ['water levels']
            self._subplot_titles = [self._titles[title] for title in titles if title in self._titles]
            self._num_rows = rows
        return self._num_rows

    @num_rows.setter
    def num_rows(self, num_rows):
        if num_rows is None:
            self._num_rows = None
            return
        if not isinstance(num_rows, int):
            raise TypeError(f'num_rows should be an integer, got {type(num_rows)}')
        if num_rows < 1:
            raise ValueError(f"num_rows must be greater than or equal to 1, got {num_rows}")
        self._num_rows = num_rows

    @property
    def num_cols(self):
        """Defines column widths based the number of columns. If show_map is True, there will be two columns.
        Otherwise there will be one column."""
        if self._num_cols is None:
            if self._show_map:
                cols = 2
                self._col_widths = [0.7, 0.3]
                self._subplot_titles.insert(1, self._titles['map'])
            else:
                cols = 1
                self._col_widths = [1]
            self._num_cols = cols
        return self._num_cols

    @num_cols.setter
    def num_cols(self, num_cols):
        if num_cols is None:
            self._num_cols = None
            return
        if not isinstance(num_cols, int):
            raise TypeError(f'num_cols should be an integer, got {type(num_cols)}')
        if num_cols < 1:
            raise ValueError(f"num_cols must be greater than or equal to 1, got {num_cols}")
        self._num_cols = num_cols

    @property
    def row_heights(self):
        return self._row_heights

    @row_heights.setter
    def row_heights(self, row_heights):
        if not isinstance(row_heights, list):
            raise TypeError(f'row_heights should be a list of heights, got {type(row_heights)}')
        if len(row_heights) != self.num_rows:
            raise ValueError(
                f'row_heights should be a list of length {self.num_rows}, got {len(row_heights)}'
            )
        if not np.isclose(np.array(row_heights).sum(), 1.0):
            raise ValueError('all row heights must sum to 1.0')
        self._row_heights = row_heights

    @property
    def col_widths(self):
        return self._col_widths

    @col_widths.setter
    def col_widths(self, col_widths):
        if not isinstance(col_widths, list):
            raise TypeError(f'col_widths should be a list of widths, got {type(col_widths)}')
        if len(col_widths) != self.num_cols:
            raise ValueError(
                f'col_widths should be a list of length {self.num_cols}, got {len(col_widths)}'
            )
        if not np.isclose(np.array(col_widths).sum(), 1.0):
            raise ValueError('all column widths must sum to 1.0')
        self._col_widths = col_widths

    @property
    def specs(self):

        specs = [[{} for col in range(self.num_cols)] for row in range(self.num_rows)]
        if self._show_map and self._show_precip:
            specs[0][1].update({'rowspan': 2, 'type': 'scattermap'})
        elif self._show_map and not self._show_precip:
            specs[0][1].update({'rowspan': 1, 'type': 'scattermap'})
        specs[0][0].update({'secondary_y': True})
        self._specs = specs

        return self._specs

    def make(self, rows, cols, specs, column_widths, row_heights,
             horizontal_spacing=0.0,
             vertical_spacing=0.05,
             shared_xaxes=True,
             *args, **kwargs
             ):
        return make_subplots(
            rows=rows,
            cols=cols,
            specs=specs,
            subplot_titles=self._subplot_titles,
            column_widths=column_widths,
            row_heights=row_heights,
            horizontal_spacing=horizontal_spacing,
            vertical_spacing=vertical_spacing,
            shared_xaxes=shared_xaxes,
            *args, **kwargs
        )

    def apply_template_layout(self):
        self._subplot.update_layout(**lt.layout.to_plotly_json())

def snsfig(
        data = None,
        x = None,
        y = None,
        hue = None,
        title: str = None,
        xlabel: str = None,
        ylabel: str = None,
):

    sns.set_context("paper")
    sns.set_theme(style="whitegrid")
    major_color = "#C0C0C0"
    minor_color = "#E0E0E0"
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Calibri"],
        "font.size": 9,
        "figure.figsize": (6.5, 4),
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        'axes.linewidth': 1,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        'grid.linewidth': 0.5,
    })

    fig, ax = plt.subplots(constrained_layout=True)

    # --- Titles and labels ---
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # Force tick marks on
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    ax.tick_params(axis="both", which="both", bottom=True, left=True, top=False, right=False)

    ax.tick_params(axis='both', which='major', labelsize=8, length=4, width=0.6, direction='out', color=major_color)
    ax.minorticks_on()
    ax.tick_params(axis='both', which='minor', length=2, width=0.4, direction='out', color=minor_color, )
    ax.grid(True, which='major', axis='both', linewidth=0.6, color=major_color)
    ax.grid(True, which='minor', axis='both', linewidth=0.4, linestyle=':', alpha=1, color=minor_color)

    """ax.xaxis.set_major_locator(mdates.YearLocator(1, 10))
    ax.xaxis.set_major_formatter(
        mdates.DateFormatter('%Y-%m')
    )
    ax.xaxis.set_minor_locator(mdates.MonthLocator())
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))"""

    sns.lineplot(data=data, x=x, y=y, hue=hue, ax=ax, )

    ax.legend(
        title="",
        frameon=True,
        loc="best"
    )
    return fig
