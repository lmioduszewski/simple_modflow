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

import numpy as np
import pandas as pd
import plotly.graph_objs as go
import rasterio
from pandas import IndexSlice as idxx
from PIL import Image
from rasterio.errors import RasterioError
from rasterio.warp import transform_bounds

from myflopy._deprecation import warn_deprecated
from myflopy._logging import get_logger
from myflopy._warnings import stacklevel_outside_package
from myflopy.modflow.mf6.contour_plotting import contour_line_segments_latlon
from myflopy.modflow.mf6.headsplus import _is_sentinel
from myflopy.modflow.mf6.hydraulics import flow_model, model_transmissivity
from myflopy.modflow.utils.animations import Animation
from myflopy.modflow.utils.datatypes.extent import Extent, resolve_extent
from myflopy.modflow.utils.datatypes.features import (
    Feature,
    combine_features,
    draw_features_mpl,
    feature_traces,
    resolve_features,
)
from myflopy.modflow.utils.datatypes.hover import thousands_formatter
from myflopy.viz import (
    COLORBAR_MARGIN,
    HIGHLIGHT_DIM_OPACITY,
    HIGHLIGHT_WIDTH,
    PALETTE,
    PLOTLY_DEFAULT_SIZE,
    Fig,
    Picture,
    log_decade_colorbar,
    map_canvas,
    mpl_color,
    mpl_colormap_for,
    plotly_colorscale_stops,
)

__all__ = [
    "Choro",
    "ChoroplethHillshadeBackground",
    "mpl_color",
    "mpl_colormap_for",
    "plotly_colorscale_stops",
]

logger = get_logger(__name__)

#: What "ask the grid for its layer surfaces" can raise.
#:
#: `gdf_topbtm` is frequently LAZY: building it samples rasters, so the call
#: reaches rasterio, pyproj and geopandas. Two MRO facts drive the membership --
#: `RasterioError` does NOT subclass OSError (only `RasterioIOError` does, by
#: multiple inheritance), and pyproj's `CRSError` is a RuntimeError.
_SURFACES_UNAVAILABLE = (
    AttributeError, IndexError, KeyError, TypeError, ValueError, RuntimeError,
    OSError, RasterioError,
)

# `mpl_colormap_for` / `plotly_colorscale_stops` / `mpl_color` live in
# `myflopy.viz` (layer 0) so the noun-tier composers can import them at module
# level; they stay reachable here because this is where callers learned them.


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
# mask MODFLOW's +/-1e30 sentinels (``Choro._layer_field``), and label the
# layer-hover table by value name.
_DEPVAR_READER_ATTR = {"hds": "hds", "conc": "conc", "temp": "temp"}

#: Every ``type`` a map draws, and whether it needs the model's RESULTS. Checked at
#: the call (ledger 305): an unknown value used to fall through to the cell-index
#: branch, so ``type='trans'`` or a typo drew node numbers 0..ncpl-1 under a
#: colorbar, with no error anywhere. ``'custom'`` is the package-grammar path,
#: which always supplies its own values.
MAP_TYPES = {
    "hds": True,
    "conc": True,
    "temp": True,
    "rch": True,
    "output_rch": True,
    "ks": False,
    "trans": True,
    "trans_max": False,
    "custom": False,
}

#: The transmissivity types (ledger 308): ``'trans'`` is SATURATED at the map's
#: time, so it needs results; ``'trans_max'`` is K times the full thickness.
TRANSMISSIVITY_TYPES = ("trans", "trans_max")


def check_map_type(value) -> str:
    """``value`` as a known map ``type``, lower-cased, or raise naming the set.

    Lower-cased because the default-hover lookup always was and the reader
    lookup was not, so ``type='HDS'`` drew node numbers under a Head hover.
    """

    key = str(value).lower()
    if key == "input_rch":
        raise ValueError(
            "type='input_rch' is retired: it read stress period 37 (zero-based) "
            "whatever period was asked for, and raised KeyError on a model with "
            "fewer. Map "
            "declared recharge with its package noun, "
            "model.packages.rch.inputs.recharge.map(per=...)."
        )
    if key not in MAP_TYPES:
        raise ValueError(
            f"type={value!r} is not a map type; expected one of {sorted(MAP_TYPES)}. "
            "Pass values= to draw any other per-cell array."
        )
    return key



def _resolve_cells(select, *, vor, model=None):
    """Turn any accepted ``select=`` spelling into ``(cells, label)``.

    One resolver for every scope, because the alternative is each caller
    re-deriving the same rules and diverging -- which is what happened to the
    geometry branch this replaces.

    Accepts cell indices, a boolean mask of length ``ncpl``, a registered model
    region name, an area-of-interest name or :class:`~myflopy.areas.Area`
    (ledger 310), a path to a vector file, or a shapely/GeoPandas geometry. A
    name that is BOTH a region and an area is refused rather than guessed, and a
    name that is neither -- and no file -- raises naming what is defined.
    Returns a flat, sorted, de-duplicated ``list[int]`` and a legend label.

    Three bugs are deliberately NOT inherited from the code this supersedes:
    ``get_vor_cells_as_series(...).to_list()`` is a list *of lists*, one per
    feature, so it must be flattened; a geometry matching nothing comes back as
    an empty ``GeoDataFrame`` rather than a Series; and a bare ``str`` path hits
    a helper that RETURNS a ``ValueError`` instead of raising it, so paths are
    coerced to ``Path`` first.
    """

    if select is None:
        return [], None

    label = "selection"

    # A registered model region ("all_streams"), an area of interest, else a
    # vector file path.
    if isinstance(select, str):
        cells = None
        if model is not None:
            try:
                cells = model.get_region_cells(select)
            except (KeyError, ValueError, AttributeError):
                cells = None
        areas = getattr(model, "aoi", None) if model is not None else None
        in_areas = areas is not None and select in areas
        if cells is not None and in_areas:
            raise ValueError(
                f"select={select!r} names BOTH a region (a set of cells) and an area "
                "of interest (a geometry), which may not select the same cells. Pass "
                f"model.aoi[{select!r}] for the area, or rename one of them."
            )
        if cells is not None:
            return sorted({int(c) for c in cells}), select
        if in_areas:
            select = areas[select]
        elif not Path(select).exists():
            registry = getattr(model, "regions", None)
            known_regions = sorted(registry.list()) if hasattr(registry, "list") else []
            known_areas = sorted(areas) if areas is not None else []
            raise ValueError(
                f"select={select!r} is not a region, an area of interest or a file. "
                f"Regions: {known_regions or 'none'}; areas: {known_areas or 'none'}."
            )
        else:
            select = Path(select)

    if hasattr(select, "to_geoseries"):
        label = getattr(select, "name", None) or label
        select = vor.get_vor_cells_as_series(select.to_geoseries())

    if isinstance(select, Path):
        label = select.stem
        select = vor.get_vor_cells_as_series(select)

    elif hasattr(select, "geom_type") or hasattr(select, "geometry"):
        select = vor.get_vor_cells_as_series(select)

    # `get_vor_cells_as_series` yields one entry PER FEATURE, each a list of
    # cells; an empty match yields an empty frame with no such entries.
    if hasattr(select, "to_list"):
        select = select.to_list()

    flat = []
    for item in np.atleast_1d(np.asarray(select, dtype=object)).ravel():
        if isinstance(item, (list, tuple, set, np.ndarray, pd.Series)):
            flat.extend(int(c) for c in np.asarray(item).ravel())
        elif item is not None:
            flat.append(item)

    values = np.asarray(flat)
    if values.size and values.dtype == bool:
        values = np.flatnonzero(values)

    cells = sorted({int(c) for c in np.asarray(values, dtype=float).astype(int)}) if values.size else []

    ncpl = len(vor.gdf_vorPolys)
    bad = [c for c in cells if not 0 <= c < ncpl]
    if bad:
        raise ValueError(
            f"select= names cell(s) outside the grid (0..{ncpl - 1}): {bad[:5]}"
            + (f" and {len(bad) - 5} more" if len(bad) > 5 else "")
        )
    if not cells:
        logger.warning("select= matched no cells; nothing is highlighted")
    return cells, label


class Choro(Picture):

    #: Whether `.fig` has already assembled its traces. A CLASS attribute, not
    #: just an `__init__` assignment: `Choro` is constructed via
    #: `object.__new__` by several test doubles, and a `.fig` that raises
    #: AttributeError on those is a contract that only half-holds.
    _assembled = False

    #: The model a map of an UNRUN model reads its inputs and names from. A class
    #: default for the same `object.__new__` reason as `_assembled`.
    _input_model = None


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
            features=None,
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
            select=None,
            select_style: str = "outline",
            select_color: str = None,
            input_model=None,
            extent=None,
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
        :param type: what to draw -- one of :data:`MAP_TYPES` ('hds', 'conc', 'temp', 'rch', 'ks',
        'custom'); anything else raises at the call.
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
        :param features: custom geometry drawn over the cells -- a vector file path, a
        GeoDataFrame, a shapely geometry, a :class:`~myflopy.plot.Feature`, or a
        list/dict of those (see :func:`myflopy.plot.map`).
        :param rch_scale: value to scale z-values in choropleth. For example to convert units in recharge to another L/T
        :param bgs: if True, will take precedence, and will plot water levels in given Layer relative to top of model
        :param hillshade_path: path to hillshade raster.
        :param colorscale: any Plotly colorscale name (``_r`` reverses), or explicit
        ``[[position, color], ...]`` stops. Default is 'earth'.
        :param logscale: colour ``log10(value)``; the colorbar is labelled in real
        units at whole decades, and ``zmin``/``zmax``/``zmid`` are given in DATA units.
        :param input_model: the model to read INPUTS and names from (K, idomain,
        region names) when ``model`` is None because it has no results yet.
        ``model`` wins when both are given.
        :param extent: the window the map opens on -- a box, points, a geometry, a
        GeoDataFrame, an area or its name (:mod:`.extent`). Sets the VIEW only.
        """

        self._vor = None
        self.model = model
        # Two different questions, which one attribute used to answer: `model`
        # is "whose RESULTS", and `_input_model` is "whose inputs and names".
        # A model that has not been run has the second and not the first, and
        # binding neither made `type='ks'` draw node numbers and every region
        # name raise as a missing file (ledger 305).
        self._input_model = input_model
        self.vor = model.vor if model is not None else vor
        type = check_map_type(type)
        if (
            MAP_TYPES[type]
            and model is None
            and input_model is not None
            and custom_zs is None
        ):
            instead = (
                "type='trans_max' (K x FULL thickness) needs no heads"
                if type == "trans"
                else "map an input instead (type='ks', or a package noun such as "
                "model.packages.npf.k.map())"
            )
            raise ValueError(
                f"type={type!r} maps a RESULT and this model has none yet -- it has "
                f"not been run, or its output files are unreadable. Run it first; "
                f"{instead}; or pass values=."
            )
        if type in TRANSMISSIVITY_TYPES:
            # Refused HERE, naming the flow model, rather than as an AttributeError
            # on `.npf` at render: a transport model has no K to multiply.
            flow_model(self._inputs)
            # Its hover shows T, K and thickness; the heads rows would be noise.
            hover_heads = False
        # Resolve which dependent-variable reader (heads/conc/temp) this map reads;
        # None for non-field types (ks/rch/...), which fall back to heads for timing.
        self._depvar_attr = _DEPVAR_READER_ATTR.get(type)
        if model is not None:
            self.kstpkper = self._timing_reader.kstpkper[0] if kstpkper is None else kstpkper
        else:
            self.kstpkper = None
        self._per = None
        self.per_timestep = per_timestep
        self.per = per
        self._layer = self._checked_layer(layer, show_mounding=show_mounding)
        self.type = type
        self._custom_hover = custom_hover
        self._custom_zs = custom_zs
        self._zmin = zmin
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
        self.rch_scale = rch_scale
        self.bgs = bgs
        self._show_mounding_above_ground = False
        if show_mounding and self.layer == -1:
            # `layer=-1` means "mound above ground". Resolved here rather than on
            # first use, so the label and the numbers cannot disagree about the
            # datum depending on which is asked for first.
            self._show_mounding_above_ground = True
            self.layer = 0
        self._colorscale = None
        self.logscale = bool(logscale)
        if self.logscale:
            # The limits are in DATA units, converted to log10 for the trace,
            # so a zero or negative one has no place on the scale. Refused here
            # rather than drawn as NaN limits, which Plotly renders as an
            # autoscaled map that quietly ignores what was asked for.
            for label, value in (("zmin", zmin), ("zmax", zmax), ("zmid", kwargs.get("zmid"))):
                if value is not None and not (np.isfinite(float(value)) and float(value) > 0):
                    raise ValueError(
                        f"{label}={value!r} cannot be drawn on a log scale: with "
                        "logscale=True the color limits are in data units and must "
                        "be finite and positive. (A signed field's symmetric limits "
                        "come from the noun's own policy -- a log scale cannot draw "
                        "a signed field; map its magnitude instead.)"
                    )
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
        legacy_locations = kwargs.pop("locs", None)
        if legacy_locations is not None:
            # Behavior deprecation, not an alias: the retired spelling reaches
            # here through every scope's **kwargs tail and is drawn through
            # the features engine, so it keeps working while nothing an editor
            # reads still offers it.
            warn_deprecated(
                "myflopy.plot.map(locs=)", "features=", since="0.2.0",
                stacklevel=stacklevel_outside_package() + 1,
            )
            features = combine_features(features, Feature(legacy_locations, name="locations"))
        self.kwargs = kwargs
        # Resolved HERE, like `select=` below, so a missing file or an unknown
        # layer raises at the call rather than from a notebook's display hook.
        self._features = resolve_features(features, crs=getattr(self.vor, "crs", None))

        if select_style not in ("outline", "dim", "both"):
            raise ValueError(
                f"select_style must be 'outline', 'dim' or 'both', not {select_style!r}"
            )
        self.select_style = select_style
        self.select_color = select_color
        # Resolved HERE, not lazily, so a bad `select=` raises at the call rather
        # than from a notebook's display hook cells later.
        self._select_cells, self._select_label = _resolve_cells(
            select, vor=self.vor, model=self._inputs
        )

        # Resolved HERE, like `select=`, so a bad box raises at the call.
        self._extent = self.resolve_extent(extent)
        if self._extent is not None and not self.fit_bounds:
            raise ValueError(
                "extent= is a window to FIT the view to, and fit_bounds=False asks for "
                "no fit; drop one of them (zoom= only acts without extent=)."
            )

        self._fig = Fig()
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

        if per is not None and self.model is None:
            # A model-less choropleth -- an INPUT map -- carries `per` only as a
            # label: which stress period's records were selected is already
            # baked into the values, and there is no results file to pick a
            # timestep from. Resolving a kstpkper here dereferenced a model that
            # is legitimately absent, so every input map on a model whose
            # results were unreadable died with
            # `AttributeError: 'NoneType' object has no attribute 'kstpkper'`.
            self._per = per
            return
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
    def _timing_reader(self):
        """The reader this map borrows its time axis and layer table from.

        Non-field maps -- ``type='custom'`` (every package-grammar map), ``ks``,
        ``rch``, ... -- have no dependent variable of their own but still need a
        ``kstpkper``. They used to borrow ``model.hds`` unconditionally, which
        the kind guard makes RAISE on a GWT/GWE model: that one fallback is why
        no package map, group diff, or budget map worked on a transport model at
        all. Borrow the model's OWN dependent variable instead.
        """

        reader = self._depvar_reader
        if reader is not None or self.model is None:
            return reader
        # getattr for the benefit of duck-typed stand-ins that are not a
        # SimulationBase; a real model always has it.
        return getattr(self.model, "_field_reader", None) or self.model.hds

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
            reader = self._timing_reader
            self._all_heads = reader.all_values
        return self._all_heads

    def _layer_values(self, layer: int) -> pd.Series:
        """The stored field on ``layer`` at :attr:`kstpkper`, one float per cell, sentinels kept."""

        return pd.to_numeric(
            self.all_heads.loc[idxx[self.kstpkper, layer], self._value_column],
            errors="coerce",
        ).reset_index(drop=True)

    def _layer_field(self, layer: int) -> pd.Series:
        """:meth:`_layer_values` with MF6's no-value sentinels as NaN: what every view draws.

        +1e30 is an inactive cell and -1e30 a DRY one, and the rule is
        ``headsplus._is_sentinel``, the one the readers use. The views here each
        tested their own sign, or none: a dry cell mapped and hovered as -1e30,
        set the colour range and made the contour levels impossible, so no contour
        drew at all; an inactive one hovered as 1e30 (ledger 315).
        """

        values = self._layer_values(layer)
        return values.mask(_is_sentinel(values))

    @property
    def all_ks(self):
        """Per-layer horizontal K array from the model's NPF package (cached)."""

        if self._all_ks is None:
            self._all_ks = self._inputs.gwf.npf.k.data
        return self._all_ks

    def _cell_vector(self, values, label: str):
        """Coerce ``values`` to a one-per-cell array sized to this grid's ``ncpl``."""

        return _as_cell_vector(values, ncpl=self.vor.ncpl, label=label)

    def _top_vector(self):
        """The model-top elevation as a one-value-per-cell array."""

        return self._cell_vector(self._inputs.gwf.modelgrid.top, "model top")

    def _bottom_vector(self, layer: int):
        """The bottom elevation of ``layer`` (zero-based) as a per-cell array."""

        return self._cell_vector(self._inputs.gwf.modelgrid.botm[layer], f"layer {layer + 1} bottom")

    def _contour_vector(self):
        """Return values used for optional contour overlays."""

        if self.contour_values is not None:
            return self._cell_vector(self.contour_values, "contour values")
        if self.contours is True:
            # DATA units, not the log10 the trace carries: a contour labelled
            # "-2" on a log K map is a reading nobody asked for.
            return self._cell_vector(self.linear_zs, "choropleth values")
        contour_key = str(self.contours).lower()
        if contour_key in {"top", "model_top", "model top"}:
            return self._top_vector()
        if contour_key in {"bottom", "botm", "layer_bottom", "layer bottom"}:
            return self._bottom_vector(self.layer)
        if contour_key in {"heads", "head", "hds"}:
            return self._cell_vector(
                self._layer_field(self.layer), f"layer {self.layer + 1} {self._value_name}"
            )
        raise ValueError(
            "contours must be False, True, 'top', 'bottom', 'heads', or use contour_values."
        )

    def _active_cell_mask(self):
        """Return active cells for the selected layer, defaulting to all cells."""

        ncpl = int(self.vor.ncpl)
        source = self._inputs
        if source is None:
            return np.ones(ncpl, dtype=bool)
        idomain = getattr(getattr(source, "gwf", None).modelgrid, "idomain", None)
        if idomain is None:
            dis = getattr(source.gwf, "disv", None) or getattr(source.gwf, "disu", None) or getattr(source.gwf, "dis", None)
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
            reader = self._timing_reader
            assert kstpkper in reader.kstpkper, f'kstpkper {kstpkper} invalid, not listed in output file'
        self._kstpkper = kstpkper

    @property
    def nlay(self):
        """Number of layers in the bound model's grid."""

        return self._inputs.gwf.modelgrid.nlay

    @property
    def number_style(self) -> str:
        """How this map writes its numbers on a colorbar: ``"thousands"`` for a
        transmissivity map (112,000, ledger 313), else ``"sig"``. Read by the
        composers, which pool the panels' colorbars into one."""

        return "thousands" if getattr(self, "type", None) in TRANSMISSIVITY_TYPES else "sig"

    @property
    def _transmissivity(self):
        """This map's :class:`~myflopy.modflow.mf6.hydraulics.Transmissivity`, cached.

        Saturated for ``type='trans'`` -- at THIS map's ``kstpkper``, read from the
        head file -- and the maximum for ``'trans_max'``. One computation for the
        free map and the npf nouns, which call this map (ledger 308).
        """

        # Keyed by the TIME it was read at: `ani` moves `kstpkper` frame by frame
        # and `per` can be reassigned, and an unkeyed cache drew every frame with
        # the first time's field (ledger 311).
        key = (self.type, self.kstpkper if self.type == "trans" else None)
        cached = getattr(self, "_transmissivity_cache", None)
        if cached is None or cached[0] != key:
            heads = None
            if self.type == "trans":
                reader = self.model.hds
                heads = np.vstack([
                    reader.array(layer=layer, kstpkper=self.kstpkper, masked=False)
                    for layer in range(int(self.nlay))
                ])
            cached = (key, model_transmissivity(self._inputs, heads=heads))
            self._transmissivity_cache = cached
        return cached[1]

    @property
    def _inputs(self):
        """The model this map reads INPUTS and names from, or None for a bare grid.

        The results model when there is one, else the input-only model a map of
        an unrun model is given. `getattr` because several test doubles build a
        `Choro` through `object.__new__` and set `model` alone.
        """

        model = getattr(self, "_model", None)
        return model if model is not None else self._input_model

    def _checked_layer(self, layer, *, show_mounding: bool) -> int:
        """``layer`` as a plain int, refused at the call if it cannot be drawn.

        `__init__` wrote it straight past the asserting setter, so ``layer=4`` on a
        4-layer model and ``layer=None`` both constructed a picture that died at
        ``.fig``, and ``layer=-1`` silently drew the BOTTOM layer through numpy's
        negative indexing (ledger 305). ``-1`` keeps its one meaning, "mound above
        the ground surface", and only beside ``show_mounding``.
        """

        if isinstance(layer, bool) or not isinstance(layer, (int, np.integer)):
            raise TypeError(f"layer must be a zero-based integer, not {layer!r}")
        layer = int(layer)
        if layer == -1 and show_mounding:
            return layer
        if self._inputs is None:
            # A bare grid carries no layer count to check against; the layer is a
            # label there, and every value was accepted before.
            return layer
        nlay = int(self.nlay)
        if not 0 <= layer < nlay:
            raise ValueError(
                f"layer={layer} is outside this model's layers 0..{nlay - 1} "
                "(zero-based; layer=-1 means 'above ground' only with show_mounding=True)."
            )
        return layer

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
                    self._hover_dict[f'Layer {lyr + 1} {label}'] = self._layer_field(lyr).to_list()

            if self.type == 'rch':
                self._hover_dict['Recharge'] = self.output_rch_zs

        # K is an INPUT, so it reads from the input model -- a model that has not
        # been run still has one to show.
        if self._inputs is not None and (self.type == 'ks' or self.hover_ks is True):
            for lyr in range(self.nlay):
                lyr_ks = self.all_ks[lyr].tolist()
                self._hover_dict[f'Layer {lyr + 1} Kh'] = lyr_ks

        if self.show_layer_elevs:
            # NOTE: this used to open by reading `self.vor.gdf_topbtm.columns`,
            # which crashed with `'NoneType' has no attribute 'columns'` on any
            # grid without a layer frame -- a `.gsf`-built one, for instance.
            # The value was assigned to `layer_nums`, immediately overwritten
            # below, and never read: the model branch takes its elevations from
            # the MODEL, which necessarily has them or it could not have run.
            if self._inputs is not None:
                botms = [
                    self._bottom_vector(lyr)
                    for lyr in range(self.nlay)
                ]
                top = self._top_vector()
            else:
                botms = self.vor.gdf_topbtm.iloc[:, 2:].to_numpy().reshape(-1, self.vor.nlay).transpose()
                top = self.vor.gdf_topbtm.iloc[:, 1].to_numpy().reshape(-1, 1).transpose()[
                    0]  # TODO why do i have to add [0]
            self._hover_dict.update(
                {'Top of Model': np.round(top, 2).tolist()})
            self._hover_dict.update(
                {
                    f'Layer {lyr + 1} Bottom': np.round(botm, 2).tolist() for lyr, botm in enumerate(botms)
                }
            )
        if self.show_mounding:
            self._hover_dict.update(
                {
                    self._mounding_label: self._mounding_series()
                }
            )
        else:
            self._hover_dict.update(
                {
                    'zs': self.linear_zs
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
            zs = self._series_minus_cell_vector(
                self._layer_field(self.layer), self.model.gwf.modelgrid.top, "model top"
            )
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

        assert int(layer) in list(range(self.nlay)), \
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
    def linear_zs(self):
        """The mapped values in DATA units, one per cell, before any log transform.

        What the hover, the contours and the static backend's ``LogNorm`` read.
        :attr:`zs` is this, log10-transformed when ``logscale`` is on, because
        Plotly has no log color axis and the trace must carry the transformed
        numbers itself.
        """

        if self.custom_zs is not None:
            return self.custom_zs

        if self._depvar_attr is not None and self.model is not None:
            if self.show_mounding is True:
                # if layer is specified as -1, show mounding above ground
                if self.layer == -1:
                    self._show_mounding_above_ground = True
                    self.layer = 0
                zs = self._mounding_series()
            elif self.bgs is True:
                zs = self._series_minus_cell_vector(
                    self._layer_field(self.layer), self.model.gwf.modelgrid.top, "model top"
                )
            else:
                zs = self._layer_field(self.layer)

        elif self.type == 'ks' and self._inputs is not None:
            zs = self._cell_vector(self.all_ks[self.layer], f"layer {self.layer + 1} Kh").tolist()

        elif self.type in TRANSMISSIVITY_TYPES and self._inputs is not None:
            zs = self._transmissivity.transmissivity[self.layer].tolist()

        elif self.type in ('rch', 'output_rch') and self.model is not None:
            # `'rch'` was documented, carried a Recharge hover row, and had no
            # branch here -- it drew node numbers (ledger 305).
            zs = self.output_rch_zs

        else:
            zs = self.vor.gdf_vorPolys.index.to_list()

        return zs

    @property
    def zs(self):
        """The z values the cell TRACE carries: :attr:`linear_zs`, log10 when ``logscale``.

        `_logscaled`, not a bare `np.log10`: the head branch reads a series that
        can arrive object-dtype, and the raw ufunc then raises "loop of ufunc
        does not support argument 0 of type float which has no callable log10
        method" -- so `model.hds.map(logscale=True)` failed outright while the
        `custom_zs` branch beside it, which already used the helper, worked. The
        helper also blanks non-positive values to NaN instead of letting -inf
        drag the shared colour range to negative infinity and blank the whole map.
        """

        values = self.linear_zs
        return self._logscaled(values) if self.logscale else values

    def _trace_color_kwargs(self):
        """``(color_kwargs, trace_kwargs)`` for the cell trace, honouring ``logscale``.

        The caller's ``zmin``/``zmax``/``zmid`` are in DATA units. Linear, they
        pass straight through. Log-scaled, they are converted to log10 to match
        the trace's values, and the colorbar gets ticks at whole decades labelled
        in real units -- under any ``colorbar`` keys the caller supplied, so a
        title or a unit suffix wins. Without this a log map read ``-3 ... 2`` for
        a field running 0.001 to 100 (ledger 74).
        """

        lo, hi = self._zmin, self._zmax
        trace_kwargs = dict(self.kwargs)
        thousands = self.number_style == "thousands"
        if not self.logscale:
            color_kwargs = {"zmin": lo, "zmax": hi}
            if thousands:
                # d3's "," groups thousands and keeps a small value's decimals, so
                # the bar reads 112,000 like the hover (ledger 313). Under any
                # `colorbar` keys the caller gave.
                caller = trace_kwargs.pop("colorbar", None)
                caller = {} if caller is None else (
                    caller if isinstance(caller, dict) else caller.to_plotly_json()
                )
                color_kwargs["colorbar"] = {"tickformat": ",", **caller}
            return color_kwargs, trace_kwargs
        lo = None if lo is None else float(np.log10(lo))
        hi = None if hi is None else float(np.log10(hi))
        if trace_kwargs.get("zmid") is not None:
            trace_kwargs["zmid"] = float(np.log10(trace_kwargs["zmid"]))
        color_kwargs = {"zmin": lo, "zmax": hi}
        values = np.asarray(self.zs, dtype=float)
        finite = values[np.isfinite(values)]
        low = lo if lo is not None else (float(finite.min()) if finite.size else None)
        high = hi if hi is not None else (float(finite.max()) if finite.size else None)
        caller = trace_kwargs.pop("colorbar", None)
        if low is not None and high is not None:
            colorbar = log_decade_colorbar(low, high, thousands=thousands)
            if caller is not None:
                caller = caller if isinstance(caller, dict) else caller.to_plotly_json()
                colorbar = {**colorbar, **caller}
            color_kwargs["colorbar"] = colorbar
        elif caller is not None:
            color_kwargs["colorbar"] = caller
        return color_kwargs, trace_kwargs

    @property
    def output_rch_zs(self):
        """The simulated recharge RATE (flux / cell area) per cell at this map's time.

        ``bud('rch').df`` is ZERO-based over the whole model -- ``node`` is
        ``layer * ncpl + cell`` -- so a record maps to its column by ``node % ncpl``
        and is summed there. This used to subtract one more, which drew every
        cell's rate on the PREVIOUS cell (map[i] == truth[i+1], the last cell 0),
        and matched on the whole node, which dropped recharge applied below
        layer 1 (ledger 311).
        """

        frame = self.model.bud('rch').df.reset_index()
        key = tuple(int(v) for v in self.kstpkper)
        frame = frame[frame['kstpkper'].map(lambda k: tuple(int(v) for v in k) == key)]
        ncpl = int(self.vor.ncpl)
        cells = frame['node'].to_numpy(dtype=int) % ncpl
        flux = np.bincount(cells, weights=frame['q'].to_numpy(dtype=float), minlength=ncpl)
        zs = flux / np.asarray(self.vor.area_list, dtype=float)
        if self.rch_scale is not None:
            zs = zs * self.rch_scale
        return zs.tolist()

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
        """Set the colorscale: any named Plotly scale, or explicit ``[[pos, color], ...]`` stops.

        Every name Plotly knows is accepted, case-insensitively, and an ``_r``
        suffix reverses it (``"Viridis_r"``); an unknown name falls back to
        ``'earth'`` with a warning that lists the options. Until 2026-09-13 this
        checked against a hand-written list of eighteen names, so
        ``colorscale="Turbo"`` warned and drew ``earth``. A list/tuple of stops
        (e.g. the gaining/losing blue-white-red scale) passes straight through to Plotly.
        """

        if isinstance(colorscale, str):
            if plotly_colorscale_stops(colorscale) is not None:
                self._colorscale = colorscale
            else:
                import plotly.colors as pc

                logger.warning(
                    'colorscale %r is not a Plotly colorscale name, using "earth" '
                    'instead; the names are %s (append "_r" to reverse one)',
                    colorscale, ", ".join(pc.named_colorscales()),
                )
                self._colorscale = 'earth'
        elif isinstance(colorscale, (list, tuple)):
            # explicit [[position, color], ...] stops (e.g. the gaining/losing
            # blue-white-red scale) pass through to Plotly untouched
            self._colorscale = [list(stop) for stop in colorscale]

    @property
    def features(self):
        """The resolved ``features=`` layers drawn over the cells (WGS84), possibly empty."""

        return list(self._features)

    #: The margins a standalone map is drawn with.
    _MAP_MARGIN = {"r": 0, "t": 20, "l": 0, "b": 0}

    #: The view :meth:`update_layout` fitted, so an export can tell it from one the
    #: caller set on ``.fig`` (a class default for `object.__new__` doubles).
    _fitted_view = None

    def _shows_colorbar(self) -> bool:
        """Whether plotly will draw (and so make room for) this map's colorbar."""

        return self.kwargs.get("showscale", True) is not False if hasattr(self, "kwargs") else True

    def update_layout(self):
        """Apply the default map style, center, and (fitted or fixed) zoom to ``self._fig``."""

        self._fig.update_layout(margin=self._MAP_MARGIN, uirevision="lock")
        view = self.map_view(
            size=map_canvas(self._fig.layout, colorbar=self._shows_colorbar())
        )
        if "center" not in view:
            view["zoom"] = self.zoom
        self._fig.update_layout(map=view)
        self._fitted_view = {key: view.get(key) for key in ("center", "zoom")}

    def export_map_view(self, size):
        """The ``map`` layout to write this picture at ``size`` (a PNG's pixels).

        Read by :meth:`~myflopy.viz.Picture.save`, which swaps it in for the
        export and puts the screen view back after. None -- write the figure as
        it stands -- when there is nothing to refit, or when the caller has MOVED
        the view on ``.fig`` since it was fitted: that view is theirs (ledger 311).
        """

        if self.vor is None or not self.fit_bounds:
            return None
        current = self.fig.layout.map  # assembles (and so fits) the figure first
        if self._fitted_view is None:
            return None
        fitted = self._fitted_view
        center = current.center.to_plotly_json() if current.center is not None else None
        moved = (
            fitted.get("zoom") is None
            or current.zoom is None
            or abs(float(current.zoom) - float(fitted["zoom"])) > 1e-9
            or center is None
            or any(abs(float(center.get(k, 0.0)) - float(fitted["center"][k])) > 1e-12
                   for k in ("lat", "lon"))
        )
        if moved:
            return None
        view = self.map_view(
            size=map_canvas(self._fig.layout, size=size, colorbar=self._shows_colorbar())
        )
        return {key: view[key] for key in ("center", "zoom") if key in view}

    #: The resolved ``extent=`` (a class default for `object.__new__` doubles).
    _extent = None

    @property
    def extent(self) -> Extent | None:
        """The window this map opens on (``extent=``), or None for the whole grid."""

        return self._extent

    def resolve_extent(self, value) -> Extent | None:
        """``value`` resolved against THIS map's grid and names, without applying it.

        What a composer calls to frame a panel it did not build; the map's own
        ``extent=`` goes through the same call.
        """

        vor = self.vor
        # Only a NAME needs the areas, and reading them can touch the project's
        # saved spec -- so every other extent (and None) never asks.
        names = (
            getattr(self._inputs, "aoi", None)
            if isinstance(value, str) and self._inputs is not None
            else None
        )
        return resolve_extent(
            value,
            crs=getattr(vor, "crs", None),
            bounds=None if vor is None else vor.gdf_vorPolys.total_bounds,
            names=names,
        )

    @property
    def latlon_bounds(self):
        """``(west, south, east, north)`` WGS84 bounds to frame: the ``extent`` when
        one was given, else the grid's; ``None`` without a grid.

        Consumed by :func:`myflopy.viz.mosaic` to fit -- and optionally to
        synchronize -- the view of each map subplot to the model data.
        """

        if self._extent is not None:
            return self._extent.wgs84
        if self.vor is None:
            return None
        west, south, east, north = self.vor.gdf_latlon.total_bounds
        return (float(west), float(south), float(east), float(north))

    def map_view(self, *, bounds=None, size=None):
        """Return the ``{style, center, zoom}`` layout for this map's subplot.

        The one framing rule, shared by the standalone map (:meth:`update_layout`),
        the mosaic subplots and the animation frames. ``bounds`` -- ``(west,
        south, east, north)`` in WGS84 -- overrides this map's own extent so
        several panels can share one synchronized view; it defaults to
        :attr:`latlon_bounds`. ``size`` is the ``(width, height)`` map canvas in
        pixels the zoom is fitted for; it defaults to a standalone map at
        plotly's default figure size (:func:`~myflopy.viz.map_canvas`).

        Centred on the extent's MIDPOINT: the standalone map used the grid's
        area centroid while the subplots used the midpoint, so an irregular
        domain was framed two ways, and only the midpoint keeps a fitted box
        inside the view.
        """

        view = {"style": "carto-voyager"}
        extent = self.latlon_bounds if bounds is None else tuple(float(b) for b in bounds)
        if extent is None:
            return view
        west, south, east, north = extent
        view["center"] = {"lat": (south + north) / 2.0, "lon": (west + east) / 2.0}
        if self.fit_bounds:
            if size is None:
                # A standalone map at plotly's default size, colorbar included.
                size = (
                    PLOTLY_DEFAULT_SIZE[0] - self._MAP_MARGIN["l"]
                    - max(self._MAP_MARGIN["r"], COLORBAR_MARGIN),
                    PLOTLY_DEFAULT_SIZE[1] - self._MAP_MARGIN["t"] - self._MAP_MARGIN["b"],
                )
            view["zoom"] = _initial_map_zoom(
                extent, padding=self.bounds_padding, width=size[0], height=size[1]
            )
        else:
            view["zoom"] = self.zoom
        return view

    def add_choropleth(self):
        """creates a choropleth map based on the provided params and adds to the fig"""
        custom_data, hover_template = _content_aware_hover(self.hover_dict)
        self.update_layout()
        self._fig.add_trace(self.get_choropleth())

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
            self._fig.add_trace(trace)

    def add_overlay(self, *traces):
        """Register extra map traces (pathlines, features, ...) drawn over the cells.

        Overlays ride with the map wherever it goes: :attr:`choropleth` draws them
        onto the figure, and :func:`myflopy.viz.mosaic` copies them into the panel
        alongside the cell trace -- which is what keeps a composed small-multiple
        from quietly showing only half the picture.
        """

        self._overlays.extend(traces)
        return self

    def _highlight_traces(self):
        """Build -- without adding -- the ``select=`` highlight outline.

        One ``go.Scattermap`` tracing the dissolved boundary of the selected
        cells, so the field underneath keeps its full opacity. That is the whole
        argument for an outline over Plotly's own ``selectedpoints``: selection
        styling on a choropleth exposes *only* ``opacity``, so "highlight" there
        can only mean dimming everything else -- which costs a measured 6.9x of
        readable contrast across the cells you did not select.

        ``mode="lines"`` is load-bearing. Plotly treats a scatter-like trace as
        selectable only when it has markers or text, so a lines-only overlay is
        immune to the user's own box/lasso gesture. Do not add markers.
        """

        if not self._select_cells or self.select_style == "dim":
            return []

        selected = self.vor.gdf_latlon.iloc[self._select_cells]
        merged = selected.geometry.union_all()
        parts = list(merged.geoms) if merged.geom_type.startswith("Multi") else [merged]

        lon, lat = [], []
        for part in parts:
            # Interiors too: a selection wrapped around an inactive island draws
            # wrong from the exterior alone.
            for ring in (part.exterior, *part.interiors):
                x, y = ring.xy
                lon.extend([*x, None])
                lat.extend([*y, None])

        return [
            go.Scattermap(
                mode="lines",
                lon=lon,
                lat=lat,
                line={"color": self.select_color or PALETTE.highlight, "width": HIGHLIGHT_WIDTH},
                name=self._select_label or "selection",
                hoverinfo="skip",
                showlegend=True,
            )
        ]

    def overlay_traces(self):
        """Every non-cell trace this map carries: contours, features, highlight, overlays.

        Built fresh and detached from :attr:`fig`, so a composer can place them in
        its own subplot cell without disturbing this map.
        """

        return [
            *self._contour_traces(),
            *self._feature_traces(),
            *self._highlight_traces(),
            *self._overlays,
        ]

    @property
    def _mounding_label(self):
        """What the mounding is measured FROM, said in the hover's own words.

        "Layer 1 Mounding" was ambiguous and, above ground, wrong: `layer=-1`
        measures from the model top but sets `self.layer = 0`, so the label
        claimed layer 1's bottom as the datum when it was not.
        """

        if self._show_mounding_above_ground:
            return "Mounding above ground surface"
        return f"Mounding above layer {self.layer + 1} bottom"

    def _mounding_series(self):
        """Head above the mounding datum, clipped at zero, one value per cell.

        Shared by the legacy content-aware hover, the sectioned `hover_spec`
        hover and :attr:`linear_zs`: three renderers of one number, and computing
        it twice is how they drift. The map DID compute it again, blanking only
        +1e30, so an inactive cell drew blank and hovered as 1e30 (ledger 315).

        An INACTIVE cell (+1e30) is blank: it has no head. A DRY cell (-1e30)
        reads 0: it has no water above the datum, like every cell whose head
        sits below it. That was the user's call (ledger 315).
        """

        datum = self._top_vector() if self._show_mounding_above_ground else self._bottom_vector(self.layer)
        heads = self._layer_values(self.layer)
        no_value = _is_sentinel(heads)
        values = pd.Series(heads.to_numpy() - datum, index=heads.index)
        values = values.mask(values < 0, 0).mask(no_value)
        return values.mask(no_value & (heads.to_numpy() < 0), 0.0)

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
            except (AttributeError, ValueError):
                # AttributeError: `per_dates` holds something without
                # `.strftime` (a plain int period label). ValueError: a date
                # that pandas cannot format. The hover loses its date footer,
                # which is strictly better than losing the map.
                logger.debug("no date label for period %s", period)
                date = None

        payload = dict(self._custom_hover) if self._custom_hover else {}
        if self.show_mounding:
            # `hover_spec` renders from this payload and never reads
            # `hover_dict`, so a mounding map showed no mounding in its hover --
            # on a model map, which is the only place mounding is meaningful,
            # the spec path is always the one taken.
            payload[self._mounding_label] = np.round(
                self._mounding_series().to_numpy(), 2
            ).tolist()
        layer_fields: dict[str, list] = {}
        if self.type in TRANSMISSIVITY_TYPES and self._inputs is not None:
            result = self._transmissivity
            layer_fields["T"] = result.transmissivity.tolist()
            payload.setdefault("K", result.k[self.layer].tolist())
            payload.setdefault("thickness", result.thickness[self.layer].tolist())
        top = botm = None
        reads_field = self.model is not None and (
            self._depvar_attr is not None or self.hover_heads
        )
        if reads_field:
            layer_fields[self._value_name] = [
                self._layer_field(lyr).tolist() for lyr in range(self.nlay)
            ]
        # The layer elevations are the GRID's, not the mapped field's, so they
        # are fetched whenever the resolved spec asks for them -- not only when
        # this map happens to read heads. Gating them on the heads branch is
        # what made `show_layer_elevs=True` INERT on every input and static-array
        # map (`npf.k`, `drn.inputs.cond`, ...): each passes `hover_heads=False`,
        # so the flag reached the spec, the spec asked for surfaces, and the
        # context had no elevations to hand it. Reported 2026-09-14 against
        # `packages.npf.k.mosaic(show_layer_elevs=True)`.
        #
        # Still GATED rather than unconditional: `gdf_topbtm` is lazy, and
        # building it samples rasters, so a map that does not ask must not pay.
        if reads_field or self._wants_surfaces():
            try:
                topbtm = self.vor.gdf_topbtm
                # The layer count comes from the FRAME when there is no model to
                # ask -- `nlay` dereferences `model.gwf`, and a bare grid that
                # carries surfaces can still show them. Columns are
                # `[cell, top, botm0, botm1, ...]`.
                nlay = self.nlay if self._inputs is not None else topbtm.shape[1] - 2
                top = topbtm.iloc[:, 1].to_list()
                botm = [topbtm.iloc[:, 2 + i].to_list() for i in range(nlay)]
            except _SURFACES_UNAVAILABLE:
                # `gdf_topbtm` is often LAZY -- building it samples rasters
                # through rasterio/pyproj/geopandas, so the reachable set runs
                # well past IndexError: RasterioError (which is NOT an OSError
                # subclass -- only RasterioIOError is), and pyproj's CRSError,
                # which is a RuntimeError. Without top/botm the hover simply
                # omits the layer elevations.
                logger.debug("no layer surfaces for the hover", exc_info=True)
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

    def _wants_surfaces(self) -> bool:
        """Whether the resolved hover spec asks for the grid's layer elevations.

        Read BEFORE the context is built, because it decides whether the context
        pays for `gdf_topbtm` at all. Safe to call from `_build_hover_context`:
        `_resolved_hover_spec` only merges dataclasses and never builds a context.
        """

        return bool(getattr(self._resolved_hover_spec(), "surfaces", False))

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
        if self.show_mounding:
            # The payload carries the numbers, but a spec only renders fields it
            # names -- so without this a `show_mounding=True` map showed none.
            spec = spec.with_fields(self._mounding_label)
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
        # Set at CONSTRUCTION rather than on the finished figure, so the dim
        # survives `mosaic`, `plot.animate` and `.ani` -- all of which rebuild
        # the cell trace from this method and would drop a post-render value.
        if self._select_cells and self.select_style in ("dim", "both"):
            extra["selectedpoints"] = tuple(self._select_cells)
            extra["unselected"] = {"marker": {"opacity": HIGHLIGHT_DIM_OPACITY}}
        color_kwargs, trace_kwargs = self._trace_color_kwargs()
        choropleth = go.Choroplethmap(
            geojson=self.vor.latlon,
            featureidkey="id",
            locations=self.vor.gdf_latlon.index.to_list(),
            z=self.zs,
            hovertemplate=hover_template,
            customdata=custom_data,
            colorscale=self.colorscale,
            **color_kwargs,
            **extra,
            **trace_kwargs,
        )
        return choropleth

    def _feature_traces(self):
        """Build -- without adding -- the ``features=`` overlay traces for this map.

        One legend entry per feature layer; points are markers, lines and
        polygon rings are polylines, filled polygons add a fill beneath their
        outline. Built by :func:`~myflopy.modflow.utils.datatypes.features
        .feature_traces` from the layers resolved at construction, so a composer
        can always ask a map for its overlays.
        """

        return feature_traces(self._features)

    @property
    def fig(self):
        """The assembled figure: cells + contours + feature layers + overlays + hillshade.

        **Assembled once.** `add_choropleth()` and friends append traces
        unconditionally, so before plan 8.1 -- when this was the `.choropleth`
        property and `.plot()` returned it -- touching the figure twice silently
        drew every trace twice. The `Picture` contract requires idempotence, so
        the assembly is cached and `choro.fig.update_layout(...)` followed by
        `choro.show()` operates on one figure, not two.
        """

        if not self._assembled:
            self.add_choropleth()
            self.add_contours()
            for trace in self._feature_traces():
                self._fig.add_trace(trace)
            # Above the cells and contours, below caller-supplied overlays, so a
            # PRT pathline still draws on top of the highlight.
            for trace in self._highlight_traces():
                self._fig.add_trace(trace)
            for trace in self._overlays:
                self._fig.add_trace(trace)
            if self.hillshade_path is not None:
                self.add_hillshade(self.hillshade_path)
            self._assembled = True
        return self._fig

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
            logger.info('reusing the existing hillshade png at %s', hillshade.png_path)
        else:
            hillshade.save_geotiff_as_png()
        hillshade.update_fig_layout(self._fig)

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
            logger.info('building animation frame for kstpkper %s', per)
            self.kstpkper = per
            choropleth = self.get_choropleth()
            frame_values = np.asarray(choropleth.z, dtype=float)
            finite_values = frame_values[np.isfinite(frame_values)]
            if finite_values.size:
                zmin = min(zmin, float(finite_values.min()))
                zmax = max(zmax, float(finite_values.max()))
            # The cell trace PLUS everything drawn over it, rebuilt per period
            # so contours follow the field. Until 2026-09-13 a frame carried
            # only the cells, so `.ani` silently dropped every overlay --
            # features, contours, the select= highlight -- that the static
            # picture and `plot.animate` both kept.
            frames.append(
                go.Frame(
                    data=[choropleth, *self.overlay_traces()],
                    name=str(per),
                    baseframe=str(periods[0]),
                )
            )

        if frames:
            # Trace-space limits: log10 of the data-unit request on a log map.
            requested = self._trace_color_kwargs()[0]
            requested_zmin = requested.get("zmin")
            requested_zmax = requested.get("zmax")
            shared_zmin = requested_zmin if requested_zmin is not None else zmin
            shared_zmax = requested_zmax if requested_zmax is not None else zmax
            base_data = copy.deepcopy(frames[0].data[0])
            colorbar = base_data.colorbar.to_plotly_json()
            if self.logscale and np.isfinite([shared_zmin, shared_zmax]).all():
                # The first frame's decades span only ITS range; the shared axis
                # spans every frame, so the labels are rebuilt over that range,
                # still under any colorbar keys the caller gave.
                caller = self.kwargs.get("colorbar")
                caller = {} if caller is None else (
                    caller if isinstance(caller, dict) else caller.to_plotly_json()
                )
                colorbar = {**colorbar, **log_decade_colorbar(shared_zmin, shared_zmax), **caller}
            coloraxis = {
                "colorscale": base_data.colorscale,
                "cmin": shared_zmin,
                "cmax": shared_zmax,
                "cauto": False,
                "colorbar": colorbar,
            }
            for trace in [base_data, *(frame.data[0] for frame in frames)]:
                trace.coloraxis = "coloraxis"
                trace.colorscale = None
                trace.zmin = None
                trace.zmax = None
        else:
            base_data = None
            coloraxis = None

        self._fig = Fig(
            data=[base_data, *copy.deepcopy(list(frames[0].data[1:]))] if frames else [],
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
        # The animation figure IS this Choro's figure from here on. Without this,
        # `_assembled` stays False and the next `.fig` access -- which is what
        # `.show()`, `.html()` and `.save()` all go through -- appends the static
        # choropleth, contours, features and overlays ON TOP of the animation.
        self._assembled = True
        self.update_layout()

        return self._fig


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
        norm=None,
        **plot_kwargs,
    ):
        """Render this choropleth with matplotlib -- the static counterpart to
        :meth:`plot` (Plotly), mirroring ``GridSection.plot_mpl`` / ``.plot``.

        Colors the Voronoi cells by the same values the interactive map uses
        (``custom_zs`` when provided, otherwise the resolved heads / Kh /
        recharge for this ``type``), so a single ``Choro`` gives both backends.
        A ``logscale`` map colours the DATA-unit values through a
        ``matplotlib.colors.LogNorm``, so its colorbar is a real log axis in
        real units rather than a linear bar of log10 exponents; ``vmin``/``vmax``
        are then in data units too. Feature layers (``features=``) and the
        ``select=`` outline are drawn; caller-registered ``add_overlay`` traces
        are Plotly objects and are not.

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
        norm
            A ready ``matplotlib.colors.Normalize`` shared with other panels --
            what the static ``mosaic``/``animate`` pass so every panel draws
            through this one renderer on one scale. Replaces ``vmin``/``vmax``
            (and the automatic ``LogNorm``) when given.

        Returns
        -------
        matplotlib.figure.Figure
            The figure the choropleth was drawn on.
        """


        values = np.asarray(self.linear_zs, dtype=float)
        gdf = self.vor.gdf_vorPolys.copy()

        if ax is None:
            fig, ax = mpl_axes(figsize=(7, 6))
        else:
            fig = ax.figure

        if cmap is None:
            cmap = mpl_colormap_for(self.colorscale)
        vmin = self._zmin if vmin is None else vmin
        vmax = self._zmax if vmax is None else vmax
        style = {"cmap": cmap, "legend": colorbar, "edgecolor": edgecolor}
        if norm is not None:
            # A SHARED scale from a composer. Non-positive values still have no
            # logarithm on a log panel, so they stay gaps.
            if self.logscale:
                values = np.where(values > 0, values, np.nan)
            style["norm"] = norm
        elif self.logscale:
            from matplotlib.colors import LogNorm

            for label, value in (("vmin", vmin), ("vmax", vmax)):
                if value is not None and not (np.isfinite(float(value)) and float(value) > 0):
                    raise ValueError(
                        f"{label}={value!r} cannot be drawn on a log scale: with "
                        "logscale=True the color limits are in data units and must "
                        "be finite and positive."
                    )
            # Non-positive values have no logarithm and draw as gaps, exactly as
            # the interactive map blanks them.
            values = np.where(values > 0, values, np.nan)
            positive = values[np.isfinite(values)]
            if positive.size:
                if vmin is None:
                    vmin = float(positive.min())
                if vmax is None:
                    vmax = float(positive.max())
                style["norm"] = LogNorm(vmin=vmin, vmax=vmax)
        else:
            style.update(vmin=vmin, vmax=vmax)
        gdf["_choro"] = values
        gdf.plot(column="_choro", ax=ax, **style, **plot_kwargs)
        if colorbar and self.number_style == "thousands" and len(fig.axes) > 1:
            # The hover's own number, on the static bar: 112,000 (ledger 313).
            from matplotlib.ticker import NullFormatter

            cbar_axis = fig.axes[-1]
            cbar_axis.yaxis.set_major_formatter(thousands_formatter())
            if self.logscale:
                cbar_axis.yaxis.set_minor_formatter(NullFormatter())
        if self.logscale and colorbar and len(fig.axes) > 1:
            self._apply_caller_colorbar_ticks(fig.axes[-1])

        # `outline_regions` and `select=` draw the same picture, so they go
        # through the same resolver -- otherwise the two backends of one map
        # object disagree about what a highlight is, which is what happened
        # before ledger 160 and is what sent `canonical_02` out of the
        # interactive map into `plot_mpl` to get cell outlines at all.
        for name in outline_regions:
            cells, _ = _resolve_cells(name, vor=self.vor, model=self._inputs)
            if cells:
                gdf.iloc[cells].boundary.plot(
                    ax=ax, color=PALETTE.mpl_highlight, linewidth=0.5
                )

        if self._select_cells and self.select_style in ("outline", "both"):
            gdf.iloc[self._select_cells].boundary.plot(
                ax=ax,
                color=self.select_color or PALETTE.mpl_highlight,
                linewidth=HIGHLIGHT_WIDTH / 2,
            )

        if self._features:
            draw_features_mpl(self._features, ax, crs=gdf.crs)

        ax.set_aspect("equal")
        ax.set_axis_off()
        if self._extent is not None:
            # LAST, so anything drawn after the cells -- features far off the grid,
            # a caller's own overlay -- cannot autoscale it away. Only with an
            # extent: a default static map keeps autoscaling to what it draws.
            xmin, ymin, xmax, ymax = self._extent.padded(self.bounds_padding)
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
        if title:
            ax.set_title(title)
        return fig

    def _apply_caller_colorbar_ticks(self, axis) -> None:
        """Put caller-supplied log ``colorbar`` ticks on the static colorbar too.

        A ``colorbar={"tickvals": ..., "ticktext": ...}`` given to a log map is
        in the TRACE's log10 positions (a PEST change map labels its decades
        ``0.1x / 1x / 10x``). The static colorbar is a real log axis in data
        units, so the positions are raised back to data units before they are
        set. Without caller ticks matplotlib's own decade labels stand.
        """

        caller = self.kwargs.get("colorbar")
        if caller is None:
            return
        caller = caller if isinstance(caller, dict) else caller.to_plotly_json()
        tickvals, ticktext = caller.get("tickvals"), caller.get("ticktext")
        if not tickvals or not ticktext:
            return
        axis.set_yticks([10.0 ** float(tick) for tick in tickvals])
        axis.set_yticklabels([str(text) for text in ticktext])
        axis.minorticks_off()

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
                                figure=self._fig,
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

    def update_fig_layout(self, fig: Fig):
        """
        Updates the layout of a choropleth map figure to show a png
        hillshade background.

        Adds the raster as an image layer and nothing else: it used to set the
        view too -- zoom 12, centred on the RASTER -- after the map had been
        fitted, so a hillshade map ignored ``fit_bounds``/``zoom`` and opened on
        the raster rather than the grid (ledger 307).

        :param fig: The choropleth map figure to update, which should be an instance
            of Choro.
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
            ))
        fig.update_traces(marker=dict(opacity=0.6), selector=dict(type="choroplethmap"))
