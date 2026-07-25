from __future__ import annotations

import copy
from typing import TYPE_CHECKING

from myflopy.viz import mpl_axes

if TYPE_CHECKING:
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor
    from myflopy.modflow.mf6.simulation.base import SimulationBase

import base64
import mimetypes
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import rasterio
import shapely as shp
from pandas import IndexSlice as idxx
from PIL import Image
from rasterio.warp import transform_bounds

from myflopy.modflow.mf6.contour_plotting import contour_line_segments_latlon
from myflopy.modflow.utils.animations import Animation
from myflopy.viz import Fig

# Map the Plotly colorscale names this class understands to the nearest
# matplotlib colormap, so plot_mpl() honors a Choro's configured colorscale.
_PLOTLY_TO_MPL_CMAP = {
    "earth": "gist_earth", "viridis": "viridis", "cividis": "cividis",
    "blues": "Blues", "reds": "Reds", "greens": "Greens", "greys": "Greys",
    "hot": "hot", "jet": "jet", "rainbow": "rainbow", "electric": "plasma",
    "portland": "turbo", "rdbu": "RdBu_r", "bluered": "coolwarm",
    "picnic": "coolwarm", "ylgnbu": "YlGnBu", "ylorrd": "YlOrRd",
}


def _content_aware_hover(name_dict: dict[str, list]):
    """Build hover metadata with readable significant digits for numeric values."""

    names = list(name_dict)
    values = [list(name_dict[name]) for name in names]
    custom_data = [list(row) for row in zip(*values, strict=False)]
    template = []
    for index, name in enumerate(names):
        non_null = next(
            (
                value
                for value in values[index]
                if value is not None and not (isinstance(value, float) and np.isnan(value))
            ),
            None,
        )
        if isinstance(non_null, (float, np.floating)):
            value_template = f"%{{customdata[{index}]:.3g}}"
        else:
            value_template = f"%{{customdata[{index}]}}"
        template.append(f"<b>{name}: </b>{value_template}<br>")
    template.append("<extra></extra>")
    return custom_data, "".join(template)


def _as_cell_vector(values, *, ncpl: int, label: str):
    """Flatten DIS/DISV per-cell arrays to one value per cell."""

    arr = np.asarray(values)
    arr = np.squeeze(arr)
    if arr.ndim == 0:
        return np.full(int(ncpl), arr.item())
    if arr.size == int(ncpl):
        return arr.reshape(-1)
    raise ValueError(f"{label} must contain one value per cell; got shape={arr.shape}, ncpl={ncpl}.")


def _initial_map_zoom(bounds, *, padding: float = 0.05, width: int = 1000, height: int = 700) -> float:
    """Estimate a Plotly map zoom that initially fits unconstrained WGS84 bounds."""

    west, south, east, north = map(float, bounds)
    south = max(south, -85.051129)
    north = min(north, 85.051129)
    scale = 1.0 + (2.0 * float(padding))

    lon_fraction = max(abs(east - west) * scale / 360.0, 1.0e-12)

    def mercator_y(latitude):
        """Web Mercator vertical fraction (0 at north pole, 1 at south) for a latitude."""

        radians = np.radians(latitude)
        return (1.0 - np.log(np.tan(radians) + (1.0 / np.cos(radians))) / np.pi) / 2.0

    lat_fraction = max(abs(mercator_y(north) - mercator_y(south)) * scale, 1.0e-12)
    lon_zoom = np.log2(float(width) / 512.0 / lon_fraction)
    lat_zoom = np.log2(float(height) / 512.0 / lat_fraction)
    return float(np.clip(min(lon_zoom, lat_zoom), 0.0, 20.0))


# Choro ``type`` -> the model's dependent-variable reader attribute. These three
# types share one code path: read the reader's ``(kstpkper, layer, cell)`` table,
# mask MODFLOW's 1e30 sentinel, and label the layer-hover table by value name.
_DEPVAR_READER_ATTR = {"hds": "hds", "conc": "conc", "temp": "temp"}


class Choro:

    def __init__(
            self,
            model: SimulationBase = None,
            vor: Vor = None,
            kstpkper: tuple = None,
            per: int = None,
            per_timestep: int | str = "last",
            layer: int = 0,
            type: str = 'hds',
            custom_hover: dict = None,
            custom_zs: list = None,
            zmin: float | int = None,
            zmax: float | int = None,
            zoom: int = 13,
            fit_bounds: bool = True,
            bounds_padding: float = 0.05,
            show_layer_elevs: bool = True,
            show_mounding: bool = False,
            hover_heads: bool = True,
            hover_ks: bool = False,
            locs: Path = None,
            rch_scale: float = None,
            bgs: bool = False,
            hillshade_path: Path = None,
            colorscale: str = None,
            logscale: bool = False,
            contours: bool | str = False,
            contour_values: list | np.ndarray | pd.Series = None,
            contour_levels: int | float | list[float] = 10,
            contour_color: str = "black",
            contour_width: float = 1.5,
            contour_name: str = None,
            contour_clip: bool = True,
            contour_resolution: int = 150,
            contour_method: str = "linear",
            animation_kstpkpers=None,
            hover_spec=None,
            hover=None,
            hover_layers=None,
            hover_surfaces=None,
            hover_fields=None,
            **kwargs

    ):
        """
        Class defining the basic choropleth plots generated from a modflow model.
        :param model: a modflow simulation object. SimulationBase is the base class for all modflow simulations.
        :param vor: a VoronoiGridPlus object.
        :param kstpkper: tuple of stress period and time step to plot
        :param per: stress period to plot. Defaults to the final saved timestep
        in that stress period.
        :param per_timestep: saved timestep to use with ``per``. Accepts
        ``"last"``, ``"first"``, a zero-based saved-output index, or an exact
        zero-based MODFLOW timestep number.
        :param layer: what layer to plot, a zero-index. 0 equals layer 1.
        :param type: type of choropleth to plot - options are 'hds', 'ks', 'input_rch', 'output_rch'
        :param custom_hover: custom dictionary of hover labels to use for choropleth.
        Must be same length as no. of cells in model.
        :param custom_zs: custom list of z values to use for choropleth. Must be same length as no. of cells in model.
        :param zmin: minimum z value to use for choropleth colorscale.
        :param zmax: maximum z value to use for choropleth colorscale.
        :param zoom: zoom level for choropleth map. Default is 13.
        :param fit_bounds: calculate an unconstrained initial view fitted to the model polygons.
        :param bounds_padding: fractional padding added around fitted model bounds.
        :param show_layer_elevs: Default is True. To show elevations of all layers on hover
        :param show_mounding: if True, colorscale will be mounding over given Layer
        :param hover_heads: boolean to show heads on hover.
        :param hover_ks: boolean to show Kh on hover.
        :param locs: specify the path of a shapefile or geopackage with location points to show on choropleth
        :param rch_scale: value to scale z-values in choropleth. For example to convert units in recharge to another L/T
        :param bgs: if True, will take precedence, and will plot water levels in given Layer relative to top of model
        :param hillshade_path: path to hillshade raster.
        :param colorscale: colorscale to use for choropleth. Default is 'earth'.
        """

        self._vor = None
        self.model = model
        self.vor = model.vor if model is not None else vor
        # Resolve which dependent-variable reader (heads/conc/temp) this map reads;
        # None for non-field types (ks/rch/...), which fall back to heads for timing.
        self._depvar_attr = _DEPVAR_READER_ATTR.get(type)
        if model is not None:
            self.kstpkper = (self._depvar_reader or self.model.hds).kstpkper[0] if kstpkper is None else kstpkper
        else:
            self.kstpkper = None
        self._per = None
        self.per_timestep = per_timestep
        self.per = per
        self._layer = layer
        self.type = type
        self._custom_hover = custom_hover
        self._custom_zs = custom_zs
        self._zmin = zmin
        self._locs = None
        self._overlays = []
        self._zmax = zmax
        self.zoom = zoom
        self.fit_bounds = fit_bounds
        self.bounds_padding = float(bounds_padding)
        if self.bounds_padding < 0:
            raise ValueError("bounds_padding must be zero or greater.")
        self.show_layer_elevs = show_layer_elevs
        self.show_mounding = show_mounding
        self.hover_heads = hover_heads
        self.hover_ks = hover_ks
        self.locs = locs
        self.rch_scale = rch_scale
        self.bgs = bgs
        self._show_mounding_above_ground = False
        self._colorscale = None
        self.logscale = logscale
        self.contours = contours
        self.contour_values = contour_values
        self.contour_levels = contour_levels
        self.contour_color = contour_color
        self.contour_width = contour_width
        self.contour_name = contour_name
        self.contour_clip = contour_clip
        self.contour_resolution = contour_resolution
        self.contour_method = contour_method
        self.animation_kstpkpers = (
            list(self.model.kstpkper)
            if animation_kstpkpers is None and self.model is not None
            else list(animation_kstpkpers or [])
        )
        self._contour_segments = []
        self.hover_spec = hover_spec
        self._hover_override = hover
        self._hover_layers = hover_layers
        self._hover_surfaces = hover_surfaces
        self._hover_fields = hover_fields
        self.kwargs = kwargs

        self.fig = Fig()
        self.vor_list = self.vor.gdf_vorPolys.geometry.to_list()
        self.cell_list = [i for i in range(len(self.vor_list))]
        self.area_list = [cell.area for cell in self.vor_list]
        self.x_list = [cell.centroid.xy[0][0] for cell in self.vor_list]
        self.y_list = [cell.centroid.xy[1][0] for cell in self.vor_list]

        self._hover_dict = self.hover_dict_default.copy()
        self._all_heads = None
        self._all_ks = None
        self._hover_dict = self.hover_dict_default

        self.hillshade_path = hillshade_path
        self.colorscale = colorscale

    @property
    def per(self):
        """The stress period being plotted (``None`` means the ``kstpkper`` given directly)."""

        return self._per

    @per.setter
    def per(self, per):
        """Select the ``kstpkper`` for stress period ``per`` per the ``per_timestep`` rule.

        Validates ``per`` against the model's available periods, then picks the
        saved timestep in that period: ``"last"``/``"first"``, an exact zero-based
        MODFLOW timestep, or a positional index into that period's saved steps.
        """

        if per is not None:
            matches = [tuple(value) for value in self.model.kstpkper if int(value[1]) == int(per)]
            if not matches:
                available = sorted({int(value[1]) for value in self.model.kstpkper})
                raise ValueError(f"Stress period {per} is unavailable. Available periods: {available}")

            selector = self.per_timestep
            if selector == "last":
                selected = matches[-1]
            elif selector == "first":
                selected = matches[0]
            elif isinstance(selector, int):
                exact = [value for value in matches if int(value[0]) == selector]
                if exact:
                    selected = exact[0]
                elif -len(matches) <= selector < len(matches):
                    selected = matches[selector]
                else:
                    raise ValueError(
                        f"per_timestep={selector} is unavailable for stress period {per}. "
                        f"Available kstpkper values: {matches}"
                    )
            else:
                raise ValueError("per_timestep must be 'first', 'last', or an integer.")
            self.kstpkper = selected
        self._per = per

    @property
    def _depvar_reader(self):
        """The dependent-variable reader (heads/conc/temp) for this map's ``type``, or None."""

        depvar_attr = getattr(self, "_depvar_attr", None)
        if depvar_attr is None or self.model is None:
            return None
        return getattr(self.model, depvar_attr)

    @property
    def _value_column(self):
        """Stored value column of the active reader (``'elev'`` for heads)."""

        reader = self._depvar_reader
        return reader.store_column if reader is not None else "elev"

    @property
    def _value_name(self):
        """Public value name of the active reader (``head``/``conc``/``temp``)."""

        reader = self._depvar_reader
        return reader.value_name if reader is not None else "head"

    @property
    def all_heads(self):
        """All-layer, all-time value table for the map's field (loaded and cached).

        Despite the name, this resolves the map's dependent-variable table --
        heads, concentration, or temperature -- per :attr:`_depvar_reader`.
        """

        if self._all_heads is None:
            reader = self._depvar_reader or self.model.hds
            self._all_heads = reader.all_values
        return self._all_heads

    @property
    def all_ks(self):
        """Per-layer horizontal K array from the model's NPF package (cached)."""

        if self._all_ks is None:
            self._all_ks = self.model.gwf.npf.k.data
        return self._all_ks

    def _cell_vector(self, values, label: str):
        """Coerce ``values`` to a one-per-cell array sized to this grid's ``ncpl``."""

        return _as_cell_vector(values, ncpl=self.vor.ncpl, label=label)

    def _top_vector(self):
        """The model-top elevation as a one-value-per-cell array."""

        return self._cell_vector(self.model.gwf.modelgrid.top, "model top")

    def _bottom_vector(self, layer: int):
        """The bottom elevation of ``layer`` (zero-based) as a per-cell array."""

        return self._cell_vector(self.model.gwf.modelgrid.botm[layer], f"layer {layer + 1} bottom")

    def _contour_vector(self):
        """Return values used for optional contour overlays."""

        if self.contour_values is not None:
            return self._cell_vector(self.contour_values, "contour values")
        if self.contours is True:
            return self._cell_vector(self.zs, "choropleth values")
        contour_key = str(self.contours).lower()
        if contour_key in {"top", "model_top", "model top"}:
            return self._top_vector()
        if contour_key in {"bottom", "botm", "layer_bottom", "layer bottom"}:
            return self._bottom_vector(self.layer)
        if contour_key in {"heads", "head", "hds"}:
            return self._cell_vector(
                self.all_heads.loc[idxx[self.kstpkper, self.layer], self._value_column].reset_index(drop=True),
                f"layer {self.layer + 1} {self._value_name}",
            )
        raise ValueError(
            "contours must be False, True, 'top', 'bottom', 'heads', or use contour_values."
        )

    def _active_cell_mask(self):
        """Return active cells for the selected layer, defaulting to all cells."""

        ncpl = int(self.vor.ncpl)
        if self.model is None:
            return np.ones(ncpl, dtype=bool)
        idomain = getattr(getattr(self.model, "gwf", None).modelgrid, "idomain", None)
        if idomain is None:
            dis = getattr(self.model.gwf, "disv", None) or getattr(self.model.gwf, "disu", None) or getattr(self.model.gwf, "dis", None)
            idomain_data = getattr(dis, "idomain", None)
            idomain = getattr(idomain_data, "array", None)
        if idomain is None:
            return np.ones(ncpl, dtype=bool)
        arr = np.asarray(idomain).squeeze()
        if arr.ndim == 0:
            return np.full(ncpl, bool(arr.item()))
        if arr.ndim == 1:
            if arr.size == ncpl:
                return arr.astype(int) != 0
            if arr.size == int(self.nlay) * ncpl:
                return arr.reshape(int(self.nlay), ncpl)[int(self.layer)].astype(int) != 0
        if arr.ndim >= 2:
            if arr.shape[0] == int(self.nlay) and int(np.prod(arr.shape[1:])) == ncpl:
                return arr.reshape(int(self.nlay), ncpl)[int(self.layer)].astype(int) != 0
            if arr.size == int(self.nlay) * ncpl:
                return arr.reshape(int(self.nlay), ncpl)[int(self.layer)].astype(int) != 0
        return np.ones(ncpl, dtype=bool)

    def _contour_clip_geometry(self):
        """Return active model-domain geometry for contour clipping."""

        if not self.contour_clip:
            return None
        mask = self._active_cell_mask()
        if mask.size != int(self.vor.ncpl) or not mask.any():
            return None
        active = self.vor.gdf_vorPolys.loc[mask]
        if active.empty:
            return None
        return active.union_all()

    def _series_minus_cell_vector(self, series, values, label: str):
        """``series`` minus a per-cell vector of ``values`` (e.g. head above a datum)."""

        vector = self._cell_vector(values, label)
        return pd.Series(pd.to_numeric(series, errors="coerce").to_numpy() - vector, index=series.index)

    @property
    def model(self):
        """The bound MODFLOW simulation (``SimulationBase``), or ``None`` for a grid-only map."""

        return self._model

    @model.setter
    def model(self, model):
        """Bind the source simulation (accepts ``None`` for a grid-only map)."""

        # assert isinstance(model, SimulationBase), 'model must be an instance of SimulationBase'
        self._model = model

    @property
    def hover_dict_default(self):
        """The hover_dict default"""
        hover_dict_default = {
            'Cell No.': self.cell_list,
            'Area': self.area_list,
            'x': self.x_list,
            'y': self.y_list,
        }
        return hover_dict_default

    @property
    def vor(self):
        """The :class:`VoronoiGridPlus` grid whose cells this map draws."""

        return self._vor

    @vor.setter
    def vor(self, vor):
        """Bind the Voronoi grid whose polygons define the map's cells."""

        # assert isinstance(vor, Vor), 'vor must be an instance of VoronoiGridPlus'
        self._vor = vor

    @property
    def kstpkper(self):
        """The ``(timestep, stress_period)`` currently selected for plotting."""

        return self._kstpkper

    @kstpkper.setter
    def kstpkper(self, kstpkper):
        """Set the plotted timestep, asserting the tuple exists in the head file."""

        if kstpkper is not None:
            assert isinstance(kstpkper, tuple), 'kstpkper must be an instance of tuple'
            reader = self._depvar_reader or self.model.hds
            assert kstpkper in reader.kstpkper, f'kstpkper {kstpkper} invalid, not listed in output file'
        self._kstpkper = kstpkper

    @property
    def nlay(self):
        """Number of layers in the bound model's grid."""

        return self.model.gwf.modelgrid.nlay

    @property
    def hover_dict(self):
        """The legacy flat ``{label: per-cell list}`` hover payload for this map.

        Adds time-step/period, per-layer heads or Kh, recharge, layer-elevation,
        mounding, below-ground, and custom-hover entries according to the plot
        ``type`` and the ``show_*`` / ``hover_*`` flags. Consumed by
        :func:`_content_aware_hover`; the newer sectioned hover uses ``hover_spec``.
        """

        if self.model is not None:
            if self.kstpkper is not None:
                kstp, kper = self.kstpkper
                self._hover_dict['Time Step'] = [kstp] * self.vor.ncpl
                self._hover_dict['Stress Period'] = [kper] * self.vor.ncpl
            if self._depvar_attr is not None or self.hover_heads is True:
                label = self._value_name.capitalize()
                for lyr in range(self.nlay):
                    lyr_heads = self.all_heads.loc[idxx[self.kstpkper, lyr], self._value_column].to_list()
                    self._hover_dict[f'Layer {lyr + 1} {label}'] = lyr_heads

            if self.type == 'ks' or self.hover_ks is True:
                for lyr in range(self.nlay):
                    lyr_ks = self.all_ks[lyr].tolist()
                    self._hover_dict[f'Layer {lyr + 1} Kh'] = lyr_ks

            if self.type == 'rch':
                self._hover_dict['Recharge'] = self.output_rch_zs

        if self.show_layer_elevs:
            layer_nums = self.vor.gdf_topbtm.columns[2:].to_list()
            if self.model is not None:
                botms = [
                    self._bottom_vector(lyr)
                    for lyr in range(self.model.gwf.modelgrid.nlay)
                ]
                top = self._top_vector()
            else:
                botms = self.vor.gdf_topbtm.iloc[:, 2:].to_numpy().reshape(-1, self.vor.nlay).transpose()
                top = self.vor.gdf_topbtm.iloc[:, 1].to_numpy().reshape(-1, 1).transpose()[
                    0]  # TODO why do i have to add [0]
            layer_nums = list(range(len(botms)))
            self._hover_dict.update(
                {'Top of Model': np.round(top, 2).tolist()})
            self._hover_dict.update(
                {
                    f'Layer {lyr + 1} Bottom': np.round(botm, 2).tolist() for lyr, botm in enumerate(botms)
                }
            )
        if self.show_mounding:
            if self.layer == -1:
                self._show_mounding_above_ground = True
                self.layer = 0
            if self._show_mounding_above_ground is True:
                layer_bottom = self._top_vector()
            else:
                layer_bottom = self._bottom_vector(self.layer)
            z_hd = self.all_heads.loc[idxx[self.kstpkper, self.layer], self._value_column].reset_index(drop=True)
            zs = pd.Series(pd.to_numeric(z_hd, errors="coerce").to_numpy() - layer_bottom, index=z_hd.index)
            # remove negative mounding values
            zs = zs.mask(zs < 0, 0)
            self._hover_dict.update(
                {
                    f'Layer {self.layer + 1} Mounding': zs
                }
            )
        else:
            self._hover_dict.update(
                {
                    'zs': self.zs
                }
            )
        if self._custom_hover:
            for key in self._custom_hover.keys():
                self._hover_dict.update(
                    {
                        f'{key}': self._custom_hover[key]
                    }
                )
        if self.bgs:
            z_hd = self.all_heads.loc[idxx[self.kstpkper, self.layer], self._value_column].reset_index(drop=True)
            zs = self._series_minus_cell_vector(z_hd, self.model.gwf.modelgrid.top, "model top")
            self._hover_dict.update(
                {
                    f'Layer {self.layer + 1} Below Ground': zs
                }
            )

        return self._hover_dict

    @property
    def layer(self):
        """The zero-based layer being plotted (0 == layer 1)."""

        return self._layer

    @layer.setter
    def layer(self, layer):
        """Set the plotted layer, asserting it is within the model's layer range."""

        assert int(layer) in list(range(self.model.gwf.modelgrid.nlay)), \
            'layer must be an integer, from 0 up to 1 less than the number of model layers'
        self._layer = layer

    @property
    def custom_zs(self):
        """Caller-supplied per-cell color values that override the computed ``zs``, if set."""

        return self._custom_zs

    @custom_zs.setter
    def custom_zs(self, custom_zs):
        """Set explicit per-cell color values (a list with one entry per grid cell)."""

        if not isinstance(custom_zs, list):
            raise ValueError('custom_zs must be an instance of list')
        assert len(custom_zs) == self.vor.ncpl, 'customs zs must be provided for every cell'
        self._custom_zs = custom_zs

    @staticmethod
    def _logscaled(values):
        """``log10`` of ``values``, with non-positive entries blanked to NaN.

        A zero or negative value has no logarithm; letting it through as -inf
        would drag the shared color range to negative infinity and blank the
        whole map. NaN draws as a gap, which is what "no log-scale value here"
        looks like.
        """

        array = np.asarray(values, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            scaled = np.log10(array)
        return np.where(np.isfinite(scaled), scaled, np.nan).tolist()

    @property
    def zs(self):
        """Defines the z values of the choropleth plot, which will be represented by a varying colorscale"""
        if self.custom_zs is not None:
            # ``logscale`` used to be dropped on this branch, so every custom map
            # that offered it silently drew a linear scale (found while adding the
            # PRT travel-time map, plan 6.3B).
            return self._logscaled(self.custom_zs) if self.logscale else self.custom_zs

        if self._depvar_attr is not None and self.model is not None:
            if self.show_mounding is True:
                # if layer is specified as -1, show mounding above ground
                if self.layer == -1:
                    self._show_mounding_above_ground = True
                    self.layer = 0
                if self._show_mounding_above_ground is True:
                    layer_bottom = self._top_vector()
                else:
                    layer_bottom = self._bottom_vector(self.layer)
                z_hd = self.all_heads.loc[idxx[self.kstpkper, self.layer], self._value_column].reset_index(drop=True)
                z_hd.loc[z_hd == 1e+30] = np.nan  # make modflow empty elevations NaN
                zs = pd.Series(pd.to_numeric(z_hd, errors="coerce").to_numpy() - layer_bottom, index=z_hd.index)
                # remove negative mounding values
                zs = zs.mask(zs < 0, 0)
            elif self.bgs is True:
                z_hd = self.all_heads.loc[idxx[self.kstpkper, self.layer], self._value_column].reset_index(drop=True)
                zs = self._series_minus_cell_vector(z_hd, self.model.gwf.modelgrid.top, "model top")
            else:
                zs = self.all_heads.loc[idxx[self.kstpkper, self.layer], self._value_column].reset_index(drop=True)
                zs.loc[zs == 1e+30] = np.nan  # make modflow empty elevations NaN

        elif self.type == 'ks' and self.model is not None:
            zs = self._cell_vector(self.all_ks[self.layer], f"layer {self.layer + 1} Kh").tolist()

        elif self.type == 'output_rch' and self.model is not None:
            zs = self.output_rch_zs

        elif self.type == 'input_rch' and self.model is not None:
            zs = self.input_rch_zs

        else:
            zs = self.vor.gdf_vorPolys.index.to_list()

        if self.logscale:
            zs = np.log10(zs)

        return zs

    @property
    def output_rch_zs(self):
        """gets model output recharge values for the specified kstpkper"""
        df = self.model.bud('rch').df
        df = df[~df.index.duplicated()].reset_index()  # remove duplicates...TODO why are there duplicates sometimes
        areas = pd.Series(self.model.vor.area_list, name='area')
        areas.index.name = 'node'
        m = pd.merge(df, areas, on='node')
        full_idx = range(self.model.vor.ncpl)
        m['rch_per_period'] = (m['q'] / m['area'])  # calc L/T recharge
        m = m[m.loc[:, 'kstpkper'] == self.kstpkper]
        m.loc[:, 'node'] = m['node'] - 1  # adjust for zero-based index
        m = m.set_index('node')
        m = m.reindex(full_idx, fill_value=0)
        zs = m['rch_per_period'].to_list()
        if self.rch_scale is not None:
            zs = (np.array(zs) * self.rch_scale).tolist()
        return zs

    @property
    def input_rch_zs(self):
        """gets model input recharge values for the specified kstpkper"""
        d = pd.DataFrame(self.model.gwf.rch.stress_period_data.data[37])
        d['cell'] = d['cellid'].apply(lambda x: x[1])
        d = d.set_index('cell')
        d = d[~d.index.duplicated()]  # remove duplicates if they exist
        full_idx = range(self.model.vor.ncpl)
        d = d.reindex(full_idx, fill_value=0)  # fill in missing cell indices
        zs = d['recharge'].to_list()
        if self.rch_scale is not None:
            zs = (np.array(zs) * self.rch_scale).tolist()
        return zs

    @property
    def colorscale(self):
        """
        Gets the current colorscale used by an instance. If no colorscale
        has been explicitly set, it defaults to 'earth'.

        :return: The current colorscale.
        :rtype: str
        """
        if self._colorscale is None:
            self._colorscale = 'earth'
        return self._colorscale

    @colorscale.setter
    def colorscale(self, colorscale):
        """Set the colorscale: a named Plotly scale (validated) or explicit ``[[pos, color], ...]`` stops.

        Named strings are matched case-insensitively against the supported set and
        fall back to ``'earth'`` with a warning if unknown. A list/tuple of stops
        (e.g. the gaining/losing blue-white-red scale) passes straight through to Plotly.
        """

        valid_colorscales = [
            'Blackbody', 'Bluered', 'Blues', 'Cividis', 'Earth', 'Electric',
            'Greens', 'Greys', 'Hot', 'Jet', 'Picnic', 'Portland', 'Rainbow',
            'RdBu', 'Reds', 'Viridis', 'YlGnBu', 'YlOrRd']
        if isinstance(colorscale, str):
            if colorscale.lower() in [scale.lower() for scale in valid_colorscales]:
                self._colorscale = colorscale
            else:
                print(f'colorscale {colorscale} not recognized, using default: "earth".\n'
                      f'colorscale options are: {valid_colorscales}')
                self._colorscale = 'earth'
        elif isinstance(colorscale, (list, tuple)):
            # explicit [[position, color], ...] stops (e.g. the gaining/losing
            # blue-white-red scale) pass through to Plotly untouched
            self._colorscale = [list(stop) for stop in colorscale]

    @property
    def locs(self):
        """The location-marker GeoDataFrame (WGS84) overlaid on the map, or ``None``."""

        return self._locs

    @locs.setter
    def locs(self, locs):
        """Set overlay locations from a shapefile/GeoPackage ``Path`` or a GeoDataFrame.

        The features are read (if a path) and reprojected to EPSG:4326; an
        unreadable or unrecognized value leaves ``locs`` unset with a printed warning.
        """

        if locs is not None:
            if isinstance(locs, Path):
                try:
                    gdf = gpd.read_file(locs)
                    gdf = gdf.to_crs(epsg=4326)  # convert to lat/lon
                    locations = gdf
                except ValueError:
                    print(f'Unable to read {locs}')
            elif isinstance(locs, gpd.GeoDataFrame):
                locs = locs.to_crs(epsg=4326)
                locations = locs
            else:
                print('could not determine locs; locs should be path or geodataframe')
                locations = None
            self._locs = locations

    def update_layout(self):
        """Apply the default map style, center, and (fitted or fixed) zoom to ``self.fig``."""

        # Set up default choropleth map styles
        if self.vor:
            map_center = {"lat": self.vor.grid_centroid.y, "lon": self.vor.grid_centroid.x}
        else:
            map_center = None
        map_layout = {
            "style": "carto-voyager",
            "center": map_center,
        }
        if self.fit_bounds and self.vor is not None:
            map_layout["zoom"] = _initial_map_zoom(
                self.vor.gdf_latlon.total_bounds,
                padding=self.bounds_padding,
            )
        else:
            map_layout["zoom"] = self.zoom
        self.fig.update_layout(
            margin={"r": 0, "t": 20, "l": 0, "b": 0},
            map=map_layout,
            uirevision="lock",
        )

    @property
    def latlon_bounds(self):
        """``(west, south, east, north)`` WGS84 grid bounds, or ``None``.

        Consumed by :func:`myflopy.viz.mosaic` to fit -- and optionally to
        synchronize -- the view of each map subplot to the model data.
        """

        if self.vor is None:
            return None
        west, south, east, north = self.vor.gdf_latlon.total_bounds
        return (float(west), float(south), float(east), float(north))

    def map_view(self, *, bounds=None):
        """Return the ``{style, center, zoom}`` layout for this map's subplot.

        Mirrors :meth:`update_layout`'s framing so a choropleth composed into a
        subplot grid zooms to the data just like the standalone map does.
        ``bounds`` -- ``(west, south, east, north)`` in WGS84 -- overrides this
        map's own extent so several panels can share one synchronized view; it
        defaults to :attr:`latlon_bounds`.
        """

        view = {"style": "carto-voyager"}
        extent = self.latlon_bounds if bounds is None else tuple(float(b) for b in bounds)
        if extent is None:
            return view
        west, south, east, north = extent
        view["center"] = {"lat": (south + north) / 2.0, "lon": (west + east) / 2.0}
        if self.fit_bounds:
            view["zoom"] = _initial_map_zoom(extent, padding=self.bounds_padding)
        else:
            view["zoom"] = self.zoom
        return view

    def add_choropleth(self):
        """creates a choropleth map based on the provided params and adds to the fig"""
        custom_data, hover_template = _content_aware_hover(self.hover_dict)
        self.update_layout()
        self.fig.add_trace(self.get_choropleth())

    def _contour_traces(self):
        """Build -- without adding -- the contour polyline traces for this map."""

        if not self.contours and self.contour_values is None:
            return []
        values = self._contour_vector()
        segments = contour_line_segments_latlon(
            self.vor,
            values,
            levels=self.contour_levels,
            label=self.contour_name or "contour",
            clip_geometry=self._contour_clip_geometry(),
            resolution=self.contour_resolution,
            method=self.contour_method,
        )
        self._contour_segments = segments
        traces = []
        for segment in segments:
            hover_text = f"{self.contour_name or 'Contour'}: {segment['level']:.6g}"
            traces.append(
                go.Scattermap(
                    mode="lines",
                    lon=segment["lon"].tolist(),
                    lat=segment["lat"].tolist(),
                    line={"color": self.contour_color, "width": self.contour_width},
                    name=self.contour_name or "Contour",
                    text=[hover_text] * len(segment["lon"]),
                    hovertemplate="%{text}<extra></extra>",
                    showlegend=False,
                )
            )
        return traces

    def add_contours(self):
        """Add optional contour lines to the choropleth map."""

        for trace in self._contour_traces():
            self.fig.add_trace(trace)

    def add_overlay(self, *traces):
        """Register extra map traces (pathlines, features, ...) drawn over the cells.

        Overlays ride with the map wherever it goes: :attr:`choropleth` draws them
        onto the figure, and :func:`myflopy.viz.mosaic` copies them into the panel
        alongside the cell trace -- which is what keeps a composed small-multiple
        from quietly showing only half the picture.
        """

        self._overlays.extend(traces)
        return self

    def overlay_traces(self):
        """Every non-cell trace this map carries: contours, locations, overlays.

        Built fresh and detached from :attr:`fig`, so a composer can place them in
        its own subplot cell without disturbing this map.
        """

        return [*self._contour_traces(), *self._locs_traces(), *self._overlays]

    def _build_hover_context(self):
        """Assemble a :class:`HoverContext` from this map's data for ``hover_spec``."""

        from myflopy.modflow.utils.datatypes.hover import HoverContext

        ncpl = int(self.vor.ncpl)
        period = step = date = None
        if self.kstpkper is not None:
            step, period = self.kstpkper
        if self.model is not None and getattr(self.model, "per_dates", None) is not None:
            try:
                dates = self.model.per_dates
                if period is not None and period < len(dates):
                    date = dates[period].strftime("%b %Y")
            except Exception:
                date = None

        payload = dict(self._custom_hover) if self._custom_hover else {}
        layer_fields: dict[str, list] = {}
        top = botm = None
        if self.model is not None and (self._depvar_attr is not None or self.hover_heads):
            heads = []
            for lyr in range(self.nlay):
                lyr_heads = self.all_heads.loc[idxx[self.kstpkper, lyr], self._value_column].to_list()
                heads.append([np.nan if h == 1e30 else h for h in lyr_heads])
            layer_fields[self._value_name] = heads
            try:
                topbtm = self.vor.gdf_topbtm
                top = topbtm.iloc[:, 1].to_list()
                botm = [topbtm.iloc[:, 2 + i].to_list() for i in range(self.nlay)]
            except Exception:
                top = botm = None

        return HoverContext(
            ncpl=ncpl,
            active_layer=self.layer,
            payload=payload,
            layer_fields=layer_fields,
            top=top,
            botm=botm,
            area=list(self.area_list),
            cells=list(self.cell_list),
            period=period,
            step=step,
            date=date,
        )

    def _resolved_hover_spec(self):
        """Merge the base spec with any call-site sugar (hover/hover_layers/...).

        ``hover=`` replaces the base spec outright; ``hover_layers`` /
        ``hover_surfaces`` / ``hover_fields`` tweak it. Sugar that does not apply
        (e.g. ``hover_layers`` on a single-value-per-cell field) is a no-op --
        the layer table simply finds no per-layer field and renders nothing.
        """

        from dataclasses import replace

        spec = self._hover_override if self._hover_override is not None else self.hover_spec
        if spec is None:
            return None
        if self._hover_layers is not None:
            spec = replace(spec, layers=self._hover_layers)
        if self._hover_surfaces is not None:
            spec = replace(spec, surfaces=self._hover_surfaces)
        if self._hover_fields:
            spec = spec.with_fields(*self._hover_fields)
        return spec

    def get_choropleth(self):
        """Build and return the ``go.Choroplethmap`` trace for the current selection.

        Uses the sectioned ``hover_spec`` (rendered to customdata + template +
        hoverlabel) when one resolves, otherwise the legacy content-aware hover, and
        colors cells by ``zs`` with the configured colorscale and z-range.
        """

        extra = {}
        resolved_spec = self._resolved_hover_spec()
        if resolved_spec is not None:
            context = self._build_hover_context()
            custom_data, hover_template, hoverlabel = resolved_spec.render(context)
            extra["hoverlabel"] = hoverlabel
        else:
            custom_data, hover_template = _content_aware_hover(self.hover_dict)
        choropleth = go.Choroplethmap(
            geojson=self.vor.latlon,
            featureidkey="id",
            locations=self.vor.gdf_latlon.index.to_list(),
            z=self.zs,
            hovertemplate=hover_template,
            customdata=custom_data,
            colorscale=self.colorscale,
            zmax=self._zmax,
            zmin=self._zmin,
            **extra,
            **self.kwargs,
        )
        return choropleth

    def _locs_traces(self, name_field='ExploName'):
        """Build -- without adding -- the location-marker traces for this map.

        ``name_field`` names the attribute column used for the trace label (falling
        back to the row index); points become markers and polygons become outlines.
        Geometry types that are neither (lines, collections) are skipped rather
        than raising, so a composer can always ask a map for its overlays.
        """

        if self.locs is None:
            return []
        traces = []
        for idx, row in self.locs.iterrows():
            try:
                name = row[name_field]
            except:
                name = idx
            geom = row.geometry
            if isinstance(geom, shp.Polygon):
                coords = geom.exterior.xy
                mode = 'lines'
            elif isinstance(geom, shp.Point):
                coords = geom.xy
                mode = 'markers'
            else:
                continue
            traces.append(
                go.Scattermap(
                    mode=mode,
                    lat=coords[1].tolist(),
                    lon=coords[0].tolist(),
                    name=name,
                    marker_color='black',
                    showlegend=False,
                )
            )
        return traces

    def add_locs(self, name_field='ExploName'):
        """Overlay each location feature as a labeled point/line trace on ``self.fig``."""

        for trace in self._locs_traces(name_field):
            self.fig.add_trace(trace)

    @property
    def choropleth(self):
        """The fully assembled figure: cells + contours + location markers + overlays + hillshade."""

        self.add_choropleth()
        self.add_contours()
        if self.locs is not None:
            self.add_locs()
        for trace in self._overlays:
            self.fig.add_trace(trace)
        if self.hillshade_path is not None:
            self.add_hillshade(self.hillshade_path)
        return self.fig

    def add_hillshade(
            self,
            tif_path: Path = None,
    ):
        """Render a hillshade GeoTIFF as a PNG image layer beneath the map traces.

        Reuses an already-exported PNG next to the model output when present,
        otherwise converts ``tif_path`` first.
        """

        hillshade = ChoroplethHillshadeBackground(tif_path, self.model)
        if hillshade.png_path.exists():
            print('hillshade png already exists, using existing file')
        else:
            hillshade.save_geotiff_as_png()
        hillshade.update_fig_layout(self.fig)

    @property
    def ani(self):
        """get animation frames for a choropleth plot"""

        frames = []
        zmin = np.inf
        zmax = -np.inf

        base_data = None
        periods = getattr(self, "animation_kstpkpers", None)
        if periods is None:
            periods = list(self.model.kstpkper)
        for per in periods:
            print(f'reading kstpkper {per}', end='\r')
            self.kstpkper = per
            choropleth = self.get_choropleth()
            frame_values = np.asarray(choropleth.z, dtype=float)
            finite_values = frame_values[np.isfinite(frame_values)]
            if finite_values.size:
                zmin = min(zmin, float(finite_values.min()))
                zmax = max(zmax, float(finite_values.max()))
            frames.append(
                go.Frame(
                    data=[choropleth],
                    name=str(per),
                    baseframe=str(periods[0]),
                )
            )

        if frames:
            requested_zmin = getattr(self, "_zmin", None)
            requested_zmax = getattr(self, "_zmax", None)
            shared_zmin = requested_zmin if requested_zmin is not None else zmin
            shared_zmax = requested_zmax if requested_zmax is not None else zmax
            base_data = copy.deepcopy(frames[0].data[0])
            coloraxis = {
                "colorscale": base_data.colorscale,
                "cmin": shared_zmin,
                "cmax": shared_zmax,
                "cauto": False,
                "colorbar": base_data.colorbar.to_plotly_json(),
            }
            for trace in [base_data, *(frame.data[0] for frame in frames)]:
                trace.coloraxis = "coloraxis"
                trace.colorscale = None
                trace.zmin = None
                trace.zmax = None
        else:
            base_data = None
            coloraxis = None

        self.fig = Fig(
            data=[base_data],
            frames=frames,
            layout=go.Layout(
                coloraxis=coloraxis,
                updatemenus=[
                    {
                        "type": "buttons",
                        "buttons": [
                            {"label": "Play", "method": "animate", "args": [None]},
                            {
                                "label": "Pause",
                                "method": "animate",
                                "args": [[None], {"mode": "immediate", "frame": {"duration": 0, "redraw": False}}],
                            },
                        ],
                    }
                ],
                sliders=Animation(self.model, periods=periods).sliders,
            ),
        )
        self.update_layout()

        return self.fig

    def plot(self):
        """Return the Plotly figure for notebook display or explicit export."""

        return self.choropleth

    def plot_mpl(
        self,
        *,
        ax=None,
        cmap: str | None = None,
        vmin: float | None = None,
        vmax: float | None = None,
        colorbar: bool = True,
        title: str | None = None,
        outline_regions: tuple[str, ...] = (),
        edgecolor: str = "none",
        **plot_kwargs,
    ):
        """Render this choropleth with matplotlib -- the static counterpart to
        :meth:`plot` (Plotly), mirroring ``GridSection.plot_mpl`` / ``.plot``.

        Colors the Voronoi cells by the same ``zs`` the interactive map uses
        (``custom_zs`` when provided, otherwise the resolved heads / Kh /
        recharge for this ``type``), so a single ``Choro`` gives both backends.

        Parameters
        ----------
        ax
            Existing matplotlib ``Axes`` to draw into (e.g. one panel of a
            multi-layer mosaic). A new figure is created when omitted.
        cmap
            Matplotlib colormap. Defaults to the nearest equivalent of this
            Choro's Plotly ``colorscale``.
        vmin, vmax
            Color limits; default to the Choro's ``zmin`` / ``zmax`` when set.
        colorbar
            Draw the colorbar legend (default ``True``).
        outline_regions
            Names of model regions to outline in black for context (requires a
            parent ``model``), e.g. ``("all_streams", "all_lakes")``.
        edgecolor
            Cell edge color (default ``"none"`` for a clean fill).

        Returns
        -------
        matplotlib.figure.Figure
            The figure the choropleth was drawn on.
        """


        values = np.asarray(self.zs, dtype=float)
        gdf = self.vor.gdf_vorPolys.copy()
        gdf["_choro"] = values

        if ax is None:
            fig, ax = mpl_axes(figsize=(7, 6))
        else:
            fig = ax.figure

        if cmap is None:
            scale = self.colorscale
            if isinstance(scale, (list, tuple)):
                # explicit color stops -> equivalent matplotlib colormap
                from matplotlib.colors import LinearSegmentedColormap

                cmap = LinearSegmentedColormap.from_list(
                    "choro_custom", [(float(pos), color) for pos, color in scale]
                )
            else:
                cmap = _PLOTLY_TO_MPL_CMAP.get(str(scale).lower(), "gist_earth")
        vmin = self._zmin if vmin is None else vmin
        vmax = self._zmax if vmax is None else vmax
        gdf.plot(column="_choro", ax=ax, cmap=cmap, legend=colorbar,
                 vmin=vmin, vmax=vmax, edgecolor=edgecolor, **plot_kwargs)

        for name in outline_regions:
            if self.model is None:
                continue
            cells = list(self.model.get_region_cells(name))
            if cells:
                gdf.iloc[cells].boundary.plot(ax=ax, color="black", linewidth=0.5)

        ax.set_aspect("equal")
        ax.set_axis_off()
        if title:
            ax.set_title(title)
        return fig

    def dash_selector(self):
        """Launch an interactive Dash app for box/lasso-selecting cells on the map.

        The selection callback returns the picked cell indices; used to interactively
        gather a set of cell numbers (runs an external Jupyter-mode server on port 8050).
        """

        # Dash is the `viz` extra, not a core dependency — import only when the
        # interactive selector is actually launched.
        import dash
        import dash_bootstrap_components as dbc
        from dash import Input, Output, dcc, html

        self.add_choropleth()

        app = dash.Dash()
        app.layout = html.Div(
            [
                dbc.Row(
                    dbc.Col(
                        [
                            dcc.Graph(
                                figure=self.fig,
                                className="flex-grow-1",
                                style={"height": "95vh"},
                                id="fig",
                            )
                        ],
                        class_name="h-100 d-flex flex-column",
                        style={"height": "95vh"},
                    ),
                    style={"height": "95vh"},
                ),
                dcc.Store(id="selected"),
                html.Div(
                    id='cell_print'
                )
            ],
            style={"height": "95vh"},
        )

        @app.callback(
            Output(component_id="selected", component_property="data"),
            Output(component_id='cell_print', component_property="children"),
            Input(component_id="fig", component_property="selectedData"),
            prevent_initial_callbacks=True,
        )
        def on_select(selectedData):
            """Dash callback: map the box/lasso selection to a list of cell indices."""

            if not selectedData:
                return None, None
            selected_cells = []
            for cell in selectedData["points"]:
                selected_cells.append(int(cell["location"]))

            return selected_cells, str(selected_cells)

        app.run(debug=True, port=8050, jupyter_mode='external', use_reloader=False)


class ChoroplethHillshadeBackground:

    def __init__(
            self,
            tif_path: Path = None,
            model: SimulationBase = None,
            png_path: Path = None,
    ):
        """Wrap a hillshade GeoTIFF for use as a Plotly map background image.

        ``png_path`` defaults to ``hillshade.png`` in the model's output folder (or
        the working directory when no ``model`` is given); ``bounds`` are read lazily.
        """

        self.tif_path = tif_path
        self.model = model
        if png_path is None:
            if model is not None:
                self.png_path = model.model_output_folder_path / 'hillshade.png'
            else:
                self.png_path = Path('hillshade.png')
        self._bounds = None

    @property
    def bounds(self):
        """
        Retrieves the geographical bounds of the raster file in EPSG:4326
        coordinate reference system. It reads the bounds from
        the source file and transforms them into the desired CRS.

        :return: A tuple containing the bounds in the format
                 (left, bottom, right, top).
        :rtype: tuple[float, float, float, float]
        """
        if self._bounds is None:
            with rasterio.open(self.tif_path) as src:
                left, bottom, right, top = transform_bounds(
                    src.crs, "EPSG:4326",
                    src.bounds.left, src.bounds.bottom,
                    src.bounds.right, src.bounds.top
                )
            self._bounds = (left, bottom, right, top)
        return self._bounds

    def save_geotiff_as_png(self):
        """
           Read a GeoTIFF and save as PNG suitable for Plotly image layer.
           - Handles 1, 2, 3, 4 bands.
           Returns: ((left, bottom, right, top) in EPSG:4326)
           """
        with rasterio.open(self.tif_path) as src:
            data = src.read()  # (bands, rows, cols)

        bands = data.shape[0]

        # Normalize to 0–255 uint8 if needed
        def norm255(a):
            """Linearly rescale a band to 0-255 ``uint8`` (all-zero if degenerate/NaN)."""

            a = a.astype(np.float32)
            amin, amax = np.nanmin(a), np.nanmax(a)
            if not np.isfinite(amin) or not np.isfinite(amax) or amax == amin:
                # fall back to zeros if degenerate
                return np.zeros_like(a, dtype=np.uint8)
            return np.clip(255 * (a - amin) / (amax - amin), 0, 255).astype(np.uint8)

        if bands == 1:
            # (H, W) grayscale
            arr = norm255(data[0])
            Image.fromarray(arr, mode="L").save(self.png_path)
        elif bands == 2:
            # Treat band2 as alpha if that makes sense; otherwise just ignore it
            rgb = np.stack([norm255(data[0])] * 3, axis=-1)  # fake RGB from band1
            alpha = norm255(data[1])
            rgba = np.dstack([rgb, alpha])
            Image.fromarray(rgba, mode="RGBA").save(self.png_path)
        elif bands == 3:
            # RGB
            rgb = np.moveaxis(data[:3], 0, 2)
            if rgb.dtype != np.uint8:
                rgb = np.dstack([norm255(data[0]), norm255(data[1]), norm255(data[2])])
            Image.fromarray(rgb, mode="RGB").save(self.png_path)
        else:
            # 4+ bands -> RGBA using first 4
            r, g, b, a = (data[0], data[1], data[2], data[3])
            rgba = np.dstack([*(norm255(x) if x.dtype != np.uint8 else x for x in (r, g, b, a))])
            Image.fromarray(rgba, mode="RGBA").save(self.png_path)

        return

    def file_to_data_uri(self):
        """
        Converts the file specified by the instance attribute `png_path` to a Data URI format.

        This method reads the content of a file, encodes its bytes in base64, deduces
        its MIME type, and constructs a Data URI string in the format:
        `data:[MIME-type];base64,[base64-string]`. By default, if the MIME type is not
        detected, it assumes `image/png`.

        :return: A string representing the file as a Data URI in base64-encoded format.
        :rtype: str
        """
        mime = mimetypes.guess_type(self.png_path)[0] or "image/png"
        with open(self.png_path, "rb") as f:
            b64 = base64.b64encode(f.read()).decode("ascii")
        return f"data:{mime};base64,{b64}"

    def update_fig_layout(self, fig: Fig, zoom: int = 12):
        """
        Updates the layout of a choropleth map figure to show a png
        hillshade background.

        :param fig: The choropleth map figure to update, which should be an instance
            of Choro.
        :param zoom: The zoom level of the map. Defaults to 12.
        :type zoom: int
        :return: None
        """
        left, bottom, right, top = self.bounds
        coords = [[left, top], [right, top], [right, bottom], [left, bottom]]

        fig.update_layout(
            map=dict(
                style="white-bg",  # minimalist canvas
                layers=[
                    dict(
                        sourcetype="image",
                        source=self.file_to_data_uri(),  # local path or URL
                        coordinates=coords,
                        below="traces",  # rendered beneath all data traces
                        opacity=1.0,
                        visible=True,
                        name='hillshade'
                    )
                ],
                center={"lat": (top + bottom) / 2, "lon": (left + right) / 2},
                zoom=zoom
            ))
        fig.update_traces(marker=dict(opacity=0.6), selector=dict(type="choroplethmap"))
