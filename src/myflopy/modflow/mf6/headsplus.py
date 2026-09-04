"""Dependent-variable output readers built on FloPy's ``HeadFile``.

MODFLOW 6 writes heads (GWF), concentration (GWT), and temperature (GWE) as the
*same* binary layout -- FloPy's ``HeadFile`` reads all three, differing only by
its ``text=`` tag. So the file-agnostic machinery (tidy frames, the ``(kstpkper,
layer, cell)`` table, the unified grammar verbs, choropleth maps) lives once in
:class:`DependentVariableFile`; :class:`HeadsPlus`, :class:`ConcResults`, and
:class:`TempResults` are thin subclasses that set the file suffix / ``text`` tag /
value label and (for heads) add the head-only extras (observations, mounding,
legacy choropleth helpers).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:

    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor

from pathlib import Path

import flopy.utils.binaryfile as bf
import numpy as np
import pandas as pd

from myflopy._optional import require
from myflopy.modflow.mf6.heads_observations import (
    get_obs_cells as _get_obs_cells,
)
from myflopy.modflow.mf6.heads_observations import (
    get_obs_heads as _get_obs_heads,
)
from myflopy.modflow.mf6.heads_observations import (
    sort_dict_by_keys as _sort_dict_by_keys,
)
from myflopy.modflow.mf6.package_plotting import (
    SpatialView,
    _apply_backend,
    refuse_noun_parameters,
    resolve_noun_hover,
)
from myflopy.modflow.utils.datatypes.datalists import convert_nested_to_int
from myflopy.modflow.utils.datatypes.hover import HoverSpec, conc_hover, head_hover, temp_hover

idxx = pd.IndexSlice  # for easy index slicing in a MultiIndex DataFrame


def _as_layer_cell_heads(data, *, nlay: int, ncpl: int):
    """Return dependent-variable data as a consistent ``(nlay, ncpl)`` array."""

    values = np.asarray(data)
    values = np.squeeze(values)

    if values.ndim == 1:
        if nlay == 1 and values.size == ncpl:
            return values.reshape(1, ncpl)
        if values.size == nlay * ncpl:
            return values.reshape(nlay, ncpl)

    if values.ndim == 2:
        if values.shape == (nlay, ncpl):
            return values
        if values.shape == (ncpl, nlay):
            return values.T
        if nlay == 1 and values.size == ncpl:
            return values.reshape(1, ncpl)
        if values.size == nlay * ncpl:
            return values.reshape(nlay, ncpl)

    if values.ndim >= 3:
        if values.shape[0] == nlay and int(np.prod(values.shape[1:])) == ncpl:
            return values.reshape(nlay, ncpl)
        if values.shape[-1] == nlay and int(np.prod(values.shape[:-1])) == ncpl:
            return np.moveaxis(values, -1, 0).reshape(nlay, ncpl)

    if values.size == nlay * ncpl:
        return values.reshape(nlay, ncpl)

    raise ValueError(
        "Could not reshape values to layer/cell form. "
        f"Got shape={values.shape}, expected nlay={nlay}, ncpl={ncpl}."
    )


class DependentVariableFile(SpatialView, bf.HeadFile):
    """A per-cell dependent-variable output reader (heads / concentration / temperature).

    Also a full :class:`SpatialView` leaf in the unified grammar --
    ``.get()/.summary()`` (tables), ``.map()/.plot()/.section()`` (panels),
    ``.mosaic()/.animate()`` (composers). ``SpatialView`` is first in the MRO
    deliberately so the grammar's ``plot`` shadows flopy's legacy
    ``LayerFile.plot``. Subclasses set the class attributes below.
    """

    #: public value label + tidy-frame column used by the unified grammar
    value_name = "value"
    #: internal column name in the ``(kstpkper, layer, cell)`` table
    store_column = "value"
    #: output-file suffix appended to ``<model name>`` when no path is given
    _output_suffix = ".bin"
    #: FloPy ``HeadFile`` binary ``text=`` selector for this field
    _binary_text = "head"
    #: label used by :meth:`summary`
    _summary_label = "value"
    #: ``Choro`` type key so the map reads this field from the model
    _choro_type = "hds"

    def __init__(
            self,
            path: Path = None,
            model=None,
            vor: Vor = None,
    ):
        """Parameters
        ----------
        path
            Path to a binary MF6 output file. If omitted, ``model`` is used with
            this class's ``_output_suffix``.
        model
            Optional parent model used to infer file paths and geometry helpers.
        vor
            Optional Voronoi/grid helper, defaulting to ``model.vor``.
        """
        from myflopy.modflow.mf6.simulation.base import SimulationBase
        if path is None:
            if model is None:
                raise ValueError("Must provide an output file or model")
            assert isinstance(model, SimulationBase), 'no valid model provided'
            self.output_path = model.model_output_folder_path / f'{model.name}{self._output_suffix}'
        else:
            self.output_path = path

        super().__init__(filename=self.output_path, text=self._binary_text)

        if model is None:
            self.model = None
        elif isinstance(model, SimulationBase):
            self.model = model
        else:
            raise ValueError("model must be an instance of SimulationBase")

        self.kstpkper = convert_nested_to_int(self.get_kstpkper())
        self.vor = self.model.vor if vor is None else vor
        self._value_table = None
        self.nper = pd.DataFrame(self.get_kstpkper()).iloc[:, 1].max() + 1
        self.numstp = pd.DataFrame(self.get_kstpkper()).iloc[:, 0].max() + 1
        self.vor_list = self.vor.gdf_vorPolys.geometry.to_list()
        self.cell_list = [i for i in range(len(self.vor_list))]
        self.area_list = [cell.area for cell in self.vor_list]
        self.x_list = [cell.centroid.xy[0][0] for cell in self.vor_list]
        self.y_list = [cell.centroid.xy[1][0] for cell in self.vor_list]
        self.crs = self.vor.crs

    @property
    def all_values(self):
        """The full ``kstpkper/layer/cell`` value table, built and cached on first access."""

        if self._value_table is None:
            self._value_table = self._build_value_table()
        return self._value_table

    def _default_hover(self, *, layers: str, surfaces: bool):
        """The default hover spec for this field's :meth:`map` (overridden per kind)."""

        return HoverSpec(
            primary=self.value_name,
            title=self.value_name,
            layers=layers,
            surfaces=surfaces,
            footer=("period", "date"),
        )

    # -- unified grammar: data verbs ---------------------------------------
    def get(
        self,
        *,
        per: int | list[int] | None = None,
        layer: int | list[int] | None = None,
        cells: int | list[int] | None = None,
    ) -> pd.DataFrame:
        """Return period-end values as a tidy ``per``/``layer``/``cell``/``<value>`` frame.

        Each stress period is reduced to its last saved timestep; flopy dry/no-data
        sentinels (``|value| >= 1e29``) become ``NaN``. All selectors are zero-based.
        """

        value = self.value_name
        frame = self.all_values.reset_index()
        frame["kstp"] = [int(key[0]) for key in frame["kstpkper"]]
        frame["per"] = [int(key[1]) for key in frame["kstpkper"]]
        period_end = frame.groupby("per")["kstp"].transform("max")
        frame = frame[frame["kstp"] == period_end].copy()
        frame = frame.rename(columns={self.store_column: value})
        frame[value] = pd.to_numeric(frame[value], errors="coerce")
        frame.loc[frame[value].abs() >= 1e29, value] = np.nan
        for column, selector in (("per", per), ("layer", layer), ("cells", cells)):
            if selector is None:
                continue
            values = (
                [int(selector)]
                if isinstance(selector, (int, np.integer))
                else [int(item) for item in selector]
            )
            frame = frame[frame["cell" if column == "cells" else column].isin(values)]
        return frame[["per", "layer", "cell", value]].reset_index(drop=True)

    def summary(self) -> pd.DataFrame:
        """Return a one-row summary of the saved field."""

        value = self.value_name
        frame = self.get()
        data = frame[value].to_numpy(dtype=float)
        finite = data[np.isfinite(data)]
        return pd.DataFrame(
            [
                {
                    "label": self._summary_label,
                    "records": int(len(frame)),
                    "periods": int(frame["per"].nunique()),
                    "layers": int(frame["layer"].nunique()),
                    "cells": int(frame["cell"].nunique()),
                    "min": float(finite.min()) if finite.size else float("nan"),
                    "max": float(finite.max()) if finite.size else float("nan"),
                    "mean": float(finite.mean()) if finite.size else float("nan"),
                }
            ]
        )

    # -- unified grammar: dimension + series hooks --------------------------
    def _spatial_layers(self) -> list[int]:
        """Every grid layer -- the field is saved for all layers."""

        return list(range(int(self.nlay)))

    def _spatial_periods(self) -> list[int]:
        """Stress periods present in the saved output."""

        return sorted({int(key[1]) for key in self.kstpkper})

    def _series_default_agg(self) -> str:
        """Collapse cells within a plotted line by mean (averaging the field, not summing)."""

        return "mean"  # heads/conc/temp average over cells; summing them is meaningless

    def _sections(
        self,
        model=None,
        *,
        line=None,
        cells: int | list[int] | None = None,
        per: int | None = None,
        layer: int | list[int] = 0,
        **kwargs,
    ):
        """Return this model's :class:`XSection` for the grammar's ``xs`` verbs."""

        del model  # single-model surface; SpatialView validates the selector
        if self.model is None:
            raise ValueError(f"{type(self).__name__}.section() requires a parent model.")
        from myflopy.modflow.utils.datatypes.xsections import XSection

        return {
            self.model.name: XSection(
                model=self.model, line=line, cells=cells, per=per, layer=layer, **kwargs
            )
        }

    def long(self, *, values: str | None = None) -> pd.Series:
        """Return the field as a long series indexed by ``kstpkper/layer/cell``."""

        values = values or self.store_column
        if values not in self.all_values.columns:
            raise KeyError(f"Value column {values!r} was not found.")
        series = pd.to_numeric(self.all_values[values], errors="coerce")
        series.name = values
        return series

    def wide(
        self,
        *,
        index: list[str] | tuple[str, ...] = ("layer", "cell"),
        values: str | None = None,
        agg: str = "first",
    ) -> pd.DataFrame:
        """Pivot the field to one row per layer/cell and one column per ``kstpkper``."""

        values = values or self.store_column
        if values not in self.all_values.columns:
            raise KeyError(f"Value column {values!r} was not found.")
        frame = self.all_values.reset_index()
        missing_index = [column for column in index if column not in frame.columns]
        if missing_index:
            raise KeyError(f"Wide index columns were not found: {missing_index}")
        wide = frame.pivot_table(
            index=list(index),
            columns="kstpkper",
            values=values,
            aggfunc=agg,
        )
        wide.columns = [f"kstpkper_{kstp}_{kper}" for kstp, kper in wide.columns]
        return wide.reset_index()

    def array(
        self,
        *,
        layer: int = 0,
        per: int | None = -1,
        kstpkper: tuple[int, int] | None = None,
        masked: bool = True,
    ) -> np.ndarray:
        """Return one layer's field at one time as a 1-D ``(ncpl,)`` array.

        The quick spatial accessor behind the choropleth maps: it picks a saved
        time, picks a ``layer``, reshapes to one value per Voronoi cell, and (by
        default) converts MODFLOW's dry/no-flow sentinel (``1e30``) to ``NaN`` so
        the result drops straight into a plot or a per-cell calculation.

        Parameters
        ----------
        layer
            Zero-based model layer (0 == layer 1).
        per
            Stress period; the *last* saved timestep in that period is returned.
            ``-1`` (default) or ``None`` returns the last saved time overall.
            Ignored when ``kstpkper`` is given.
        kstpkper
            Exact ``(kstp, kper)`` to read, taking precedence over ``per``.
        masked
            Replace MODFLOW's ``1e30`` dry/no-flow sentinel with ``NaN``
            (default ``True``).

        Returns
        -------
        numpy.ndarray
            A ``(ncpl,)`` float array, one value per Voronoi cell.
        """

        # Use the normalized kstpkper list the class trusts (the value table reads
        # with these); flopy's raw get_kstpkper() is offset for some files.
        keys = [tuple(int(v) for v in key) for key in self.kstpkper]
        if kstpkper is not None:
            key = tuple(int(v) for v in kstpkper)
            if key not in keys:
                raise ValueError(f"kstpkper {key} is unavailable. Available: {keys}")
        elif per in (None, -1):
            key = keys[-1]
        else:
            matches = [k for k in keys if k[1] == int(per)]
            if not matches:
                available = sorted({k[1] for k in keys})
                raise ValueError(f"Stress period {per} is unavailable. Available periods: {available}")
            key = matches[-1]

        ncpl = int(self.vor.ncpl)
        field = _as_layer_cell_heads(self.get_data(kstpkper=key), nlay=self.nlay, ncpl=ncpl)
        values = np.asarray(field[int(layer)], dtype=float)
        if masked:
            values[np.abs(values) > 1.0e29] = np.nan
        return values

    def to_xugrid(self, *, layers=None, times=None, name: str | None = None, masked: bool = True):
        """Export the simulated field across layers and time as an xugrid object.

        Stacks the saved field into a single ``(time, layer, cell)``
        :class:`xugrid.UgridDataArray` on this model's Voronoi mesh -- the
        convenience behind ``model.to_xugrid()``. Unlike :meth:`array` (one
        layer at one time), this returns the whole history at once, ready for
        xarray-style slicing (``uda.isel(time=-1, layer=0)``), native unstructured
        plotting (``.ugrid.plot()``), and UGRID-NetCDF sharing
        (``.ugrid.to_netcdf(...)``). Uses the same normalized ``kstpkper`` and
        dry-cell masking as :meth:`array`, and the grid topology from
        :meth:`~myflopy.modflow.mf6.grid.voronoi.VoronoiGridPlus.ugrid2d`.

        Parameters
        ----------
        layers
            Zero-based layers to include (default: all ``nlay`` layers).
        times
            ``(kstp, kper)`` keys to include (default: every saved time). Order is
            preserved; an unavailable key raises ``ValueError``.
        name
            Variable name for the field (default: this reader's ``value_name``).
        masked
            Replace MODFLOW's ``1e30`` dry/no-flow sentinel with ``NaN`` (default).

        Returns
        -------
        xugrid.UgridDataArray
            Dims ``("time", "layer", <face_dim>)``, with ``kstp``/``kper`` as
            coordinates on the ``time`` axis and the zero-based layer index on
            ``layer``.

        Raises
        ------
        ImportError
            If the optional ``xugrid`` / ``xarray`` packages are not installed.
        ValueError
            If a requested ``times`` key is not a saved output time.

        Notes
        -----
        **Analysis vs. plotting shape.** The result keeps its full
        ``(time, layer, cell)`` shape -- that is what you want for differencing,
        aggregation, and NetCDF export. But ``.ugrid.plot()`` draws *one value per
        cell*, so reduce to a single field first with ``.isel(time=-1, layer=0)``
        (or ``.sel(...)``). xugrid does **not** facet over extra dimensions, so
        ``.ugrid.plot(col="layer")`` will not work; plot several layers by looping
        and passing ``ax=`` (see Examples). The same applies to NetCDF if you only
        want one slice -- though the full 3-D array also writes fine.

        Examples
        --------
        >>> uda = model.to_xugrid()                  # (time, layer, cell)
        >>> uda.isel(time=-1, layer=0).ugrid.plot()  # one map: last time, top layer
        >>> uda.isel(time=-1).ugrid.to_netcdf("heads.nc")   # share all layers

        Plot every layer as its own subplot:

        >>> import matplotlib.pyplot as plt
        >>> final = uda.isel(time=-1)                       # (layer, cell)
        >>> fig, axes = plt.subplots(1, final.sizes["layer"])
        >>> for k, ax in zip(final["layer"].values, axes):
        ...     final.sel(layer=k).ugrid.plot(ax=ax)
        """

        xu = require("xugrid", feature="to_xugrid() unstructured-grid export")
        import xarray as xr  # guaranteed present: xarray is a xugrid dependency

        name = name or self.value_name
        ncpl = int(self.vor.ncpl)
        available = [tuple(int(v) for v in key) for key in self.kstpkper]
        if times is None:
            keys = available
        else:
            keys = [tuple(int(v) for v in key) for key in times]
            missing = [key for key in keys if key not in available]
            if missing:
                raise ValueError(
                    f"to_xugrid(): kstpkper {missing} unavailable. "
                    f"Available: {available}"
                )
        lays = list(range(self.nlay)) if layers is None else [int(layer) for layer in layers]

        data = np.empty((len(keys), len(lays), ncpl), dtype=float)
        for ti, key in enumerate(keys):
            field = _as_layer_cell_heads(self.get_data(kstpkper=key), nlay=self.nlay, ncpl=ncpl)
            for li, layer in enumerate(lays):
                data[ti, li, :] = np.asarray(field[layer], dtype=float)
        if masked:
            data[np.abs(data) > 1.0e29] = np.nan

        grid = self.vor.ugrid2d()
        face_dim = grid.face_dimension
        array = xr.DataArray(
            data,
            dims=("time", "layer", face_dim),
            coords={
                "time": np.arange(len(keys)),
                "kstp": ("time", [key[0] for key in keys]),
                "kper": ("time", [key[1] for key in keys]),
                "layer": lays,
            },
            name=name,
        )
        return xu.UgridDataArray(array, grid)

    def map(
        self,
        *,
        # -- which numbers ---------------------------------------------------
        per: int | None = None,
        kstpkper: tuple[int, int] | None = None,
        per_timestep: str | int = "last",
        layer: int = 0,
        # -- colour ----------------------------------------------------------
        zmin: float | None = None,
        zmax: float | None = None,
        colorscale: str | list | tuple | None = None,
        logscale: bool = False,
        # -- contours --------------------------------------------------------
        contours: bool | str = False,
        contour_values=None,
        contour_levels: int | float | list = 10,
        contour_color: str = "black",
        contour_width: float = 1.5,
        contour_name: str | None = None,
        contour_clip: bool = True,
        contour_resolution: int = 150,
        contour_method: str = "linear",
        # -- highlighting ----------------------------------------------------
        select=None,
        select_style: str = "outline",
        select_color: str | None = None,
        # -- overlays --------------------------------------------------------
        locs=None,
        hillshade_path=None,
        bgs: bool = False,
        # -- framing ---------------------------------------------------------
        fit_bounds: bool = True,
        bounds_padding: float = 0.05,
        # -- hover -----------------------------------------------------------
        hover=None,
        hover_layers: str = "active+strip",
        hover_surfaces: bool = False,
        show_layer_elevs: bool | None = None,
        show_mounding: bool = False,
        # -- renderer --------------------------------------------------------
        backend: str = "plotly",
        **trace_kwargs,
    ):
        """This field as a choropleth, optionally with contours and overlays.

        Every parameter is named rather than forwarded through ``**kwargs``:
        PyCharm and Pylance read the ``def`` line and never run the module, so a
        parameter that arrives through a tail is one an editor can never offer
        (plan 8.8). ``show_mounding`` is the reason it matters here -- it does not
        decorate the map, it rewrites every value on it.

        ``hover`` accepts a :class:`~myflopy.modflow.utils.datatypes.hover.HoverSpec`
        for full control of the sectioned, styled hover; otherwise this field's
        default spec is used. ``hover_layers`` (``"active"`` | ``"active+strip"`` |
        ``"all"``) sets how the vertical profile shows, and ``hover_surfaces=True``
        adds the model-top/layer-bottom column (merged on ``"all"``).

        Parameters
        ----------
        backend : {'plotly', 'mpl'}, default 'plotly'
            Renderer. ``'mpl'`` returns a Matplotlib figure instead of a Picture.
        per : int, optional
            Stress period to read (0-based). Mutually exclusive with ``kstpkper``;
            with neither, the model's first output time is used.
        kstpkper : tuple of (int, int), optional
            Exact ``(timestep, period)`` to read, as MODFLOW reports it. Use
            ``model.kstpkper`` to list what is available.
        per_timestep : {'last', 'first'} or int, default 'last'
            Which timestep WITHIN ``per`` to read, when a period has several. Ignored
            when ``kstpkper`` names the timestep outright.
        layer : int, default 0
            Zero-based layer. Layer 0 is the top.
        zmin, zmax : float, optional
            Fixed color-scale limits. Set both to hold the scale steady across
            frames or panels; ``zmin >= zmax`` raises rather than rendering one flat
            color. With neither, the range comes from the data.
        colorscale : str or list of (float, str), optional
            A Plotly colorscale name, or explicit stops. **Pass stops for a diverging
            scale** -- names round-trip through a plotly-to-matplotlib table that
            maps ``'rdbu'`` to the REVERSED colormap, so a named diverging scale
            renders mirrored between backends (ledger 69/70).
        logscale : bool, default False
            Color on a log scale. Non-positive values are masked.
        contours : bool or str, default False
            Overlay contour lines. ``True`` contours the mapped values; a string
            names a different field to contour instead.
        contour_values : sequence of float, optional
            Contour a supplied array rather than the mapped one.
        contour_levels : int or float or list of float, default 10
            A count of levels, a fixed interval, or explicit level values.
        contour_color : str, default 'black'
            Line colour for the contour trace.
        contour_width : float, default 1.5
            Line width for the contour trace.
        contour_name : str, optional
            Legend name for the contour trace.
        contour_clip : bool, default True
            Clip contours to the active domain instead of the full grid extent.
        contour_resolution : int, default 150
            Grid size used to build the contours, per axis. Higher is smoother and
            slower. **Only acts under ``contour_method="cubic"``**, which is the one
            that interpolates onto a square grid before contouring; the default
            linear method triangulates the cell centres directly, so there is no
            grid for this to size and passing it changes nothing (measured: 338
            contour points at 40, 150 and 400 alike).
        contour_method : {'linear', 'cubic'}, default 'linear'
            How the scattered cell values become contours. ``'linear'``
            triangulates the cell centres and contours the triangulation --
            fast, exact at the centres, and faceted. ``'cubic'`` interpolates onto
            a ``contour_resolution``-square grid first (Clough-Tocher) and contours
            that -- smoother, slower, and able to overshoot between cells.
            ``'tri'``/``'tricontour'`` and ``'clough'``/``'clough_tocher'``/
            ``'cloughtocher'`` are accepted as aliases. Anything else raises naming
            both; ``'nearest'`` in particular was documented here for a while and
            has never been implemented.
        select : sequence of int, ndarray, str, Path or geometry, optional
            Cells to highlight. Cell indices, a boolean mask of length ``ncpl``, the
            name of a registered model region (``"all_streams"``), a path to a
            vector file, or a shapely/GeoPandas geometry to intersect. Highlights
            nothing (with a warning) when the selection is empty.
        select_style : {'outline', 'dim', 'both'}, default 'outline'
            How the highlight is drawn. ``'outline'`` traces the dissolved boundary
            of the selection and leaves every cell at full opacity. ``'dim'`` fades
            the *unselected* cells to 20% instead, which suits a bare grid where the
            selection is the subject, but on a field it costs a measured 6.9x of
            readable contrast everywhere you did not select -- and one box/lasso
            gesture in the browser overwrites it. ``'both'`` draws each.
        select_color : str, optional
            Highlight colour, defaulting to ``viz.PALETTE.highlight``. Worth setting
            when the default red collides with a red-blue diverging colorscale.
        locs : Path or GeoDataFrame, optional
            Point locations to mark -- wells, observations, samples. A path is read
            as a vector file.
        hillshade_path : Path, optional
            A hillshade GeoTIFF to draw beneath the cells for topographic context.
        bgs : bool, default False
            Draw the basemap beneath a semi-transparent cell layer.
        fit_bounds : bool, default True
            Fit the initial view to the grid extent rather than using ``zoom``.
        bounds_padding : float, default 0.05
            Fractional padding around the fitted bounds.
        hover : HoverSpec, optional
            Replace the sectioned hover outright. See
            :mod:`myflopy.modflow.utils.datatypes.hover`.
        hover_layers : {'active', 'active+strip', 'all', 'none'}, optional
            How the per-layer profile renders in the hover.
        hover_surfaces : bool, optional
            Add the model-top / layer-bottom table to the sectioned hover.
        show_layer_elevs : bool, optional
            Add model-top and per-layer-bottom rows to the hover. Defaults to
            whether the grid actually carries layer elevations (``vor.gdf_topbtm``),
            because forcing it on a grid without them raises.
        show_mounding : bool, default False
            Add head-above-initial (mounding) to the hover.

        Returns
        -------
        Choro or matplotlib.figure.Figure
            With ``backend='plotly'`` (the default), a
            :class:`~myflopy.viz.Picture`: it renders itself in Jupyter, and answers
            ``.fig``, ``.show()``, ``.save(path)`` and ``.html(path)``. It also
            carries ``.plot_mpl()`` for a static rendering and ``.ani`` for the
            animation over periods.

            With ``backend='mpl'``, a bare Matplotlib ``Figure`` -- not a Picture, so
            use ``.savefig(path)`` and ``.axes[0]`` rather than the picture verbs.

        Raises
        ------
        TypeError
            If given ``values=`` or ``type=`` -- this noun draws its own field, so
            an override would repaint the cells while the hover, title and
            colorscale went on describing the real one -- or one of the legacy
            hover arguments, which a noun's own ``hover_spec`` supersedes.

        See Also
        --------
        section : the same field as a vertical slice.
        array : the values behind the picture, as a plain array.
        myflopy.plot.map : the free verb, which DOES take ``values=``.

        Examples
        --------
        >>> model.hds.map()                                  # last time, layer 0
        >>> model.hds.map(per=3, layer=1)
        >>> model.hds.map(kstpkper=(0, 5), contours=True, contour_levels=10)
        >>> model.hds.map(zmin=100, zmax=125)                # a scale held steady
        >>> model.hds.map(show_mounding=True)                # head above initial
        >>> model.hds.map(select="all_streams", select_style="both")
        >>> model.hds.map(backend="mpl").savefig("heads.png")
        """

        if self.model is None:
            raise ValueError(f"{type(self).__name__}.map() requires a parent model.")
        refuse_noun_parameters(f"model.{self._choro_type}", self._choro_type, trace_kwargs)
        hover = resolve_noun_hover(hover, trace_kwargs)
        if hover is None:
            hover = self._default_hover(layers=hover_layers, surfaces=hover_surfaces)
        choro = self.model.plot.map(
            per=per,
            kstpkper=kstpkper,
            per_timestep=per_timestep,
            layer=layer,
            type=self._choro_type,
            zmin=zmin,
            zmax=zmax,
            colorscale=colorscale,
            logscale=logscale,
            contours=contours,
            contour_values=contour_values,
            contour_levels=contour_levels,
            contour_color=contour_color,
            contour_width=contour_width,
            contour_name=contour_name,
            contour_clip=contour_clip,
            contour_resolution=contour_resolution,
            contour_method=contour_method,
            select=select,
            select_style=select_style,
            select_color=select_color,
            locs=locs,
            hillshade_path=hillshade_path,
            bgs=bgs,
            fit_bounds=fit_bounds,
            bounds_padding=bounds_padding,
            show_layer_elevs=show_layer_elevs,
            show_mounding=show_mounding,
            hover_spec=hover,
            **trace_kwargs,
        )
        return _apply_backend(choro, backend)

    def _build_value_table(self):
        """Build the full ``(kstpkper, layer, cell)`` value table for this model."""

        vor_cell_list = list(self.vor.gdf_vorPolys.index)

        """set generic MultiIndex for all stress periods and all cells"""
        mdx = pd.MultiIndex.from_product(
            iterables=[
                self.kstpkper,
                list(range(self.nlay)),
                vor_cell_list
            ],
            names=['kstpkper', 'layer', 'cell']
        )
        """Set up a MultiIndex DataFrame to hold the values
            for all cells and stress periods"""
        df_values = pd.DataFrame(
            index=mdx,
            columns=[self.store_column]
        )
        """get data for each stress period"""
        for kstpkper in self.kstpkper:
            sp_values = _as_layer_cell_heads(
                self.get_data(kstpkper=kstpkper),
                nlay=self.nlay,
                ncpl=len(vor_cell_list),
            )
            """copy and paste this stress period data to the MultiIndex DataFrame"""
            for layer in range(self.nlay):
                assert self.vor.ncpl == sp_values.shape[1], (
                    'Are you using the wrong voronoi grid??? \n'
                    f'The provided vor grid has {self.vor.ncpl} cells, but there are {sp_values.shape[1]} values '
                    f'in the model'
                )
                df_values.loc[idxx[kstpkper, layer, :]] = sp_values[layer].reshape(-1, 1)
        return df_values


class HeadsPlus(DependentVariableFile):
    """Extended heads-file reader (GWF) with head-only observation + choropleth extras.

    ``model.hds``. Adds the head observation helpers, the legacy ``choropleth``
    builders, and the head-below-bottom "dry" affordances on top of the shared
    :class:`DependentVariableFile` grammar.
    """

    value_name = "head"
    store_column = "elev"
    _output_suffix = ".hds"
    _binary_text = "head"
    _summary_label = "hds"
    _choro_type = "hds"

    def __init__(
            self,
            hds_path: Path = None,
            model=None,
            vor: Vor = None,
            obs_path: Path = None
    ):
        """Parameters
        ----------
        hds_path
            Path to a binary MF6 heads file. If omitted, ``model`` is used.
        model
            Optional parent model used to infer file paths and geometry helpers.
        vor
            Optional Voronoi/grid helper, defaulting to ``model.vor``.
        obs_path
            Optional observation-point file used by observation helper methods.
        """

        super().__init__(path=hds_path, model=model, vor=vor)
        # Historical aliases: ``.hds_path`` and a ``.hds`` self-reference (no
        # second binary reader is opened for the same file).
        self.hds_path = self.output_path
        self.hds = self
        self.obs_heads_df = None
        self._obs = {}
        self._obs_heads = None
        self.obs_path = obs_path

    def _default_hover(self, *, layers: str, surfaces: bool):
        """The default heads hover spec."""

        return head_hover(layers=layers, surfaces=surfaces)

    @property
    def all_heads(self):
        """The full ``kstpkper/layer/cell`` head table (alias of ``all_values``)."""

        return self.all_values

    def get_all_heads(self):
        """Build the full head table (alias of the generic value-table builder)."""

        return self._build_value_table()

    @property
    def obs(self):
        """Mapping of registered head-observation location sets."""

        return self._obs

    @property
    def obs_path(self):
        """Filesystem path to the head observation output file, if configured."""

        return self._obs_path

    @obs_path.setter
    def obs_path(self, val):
        """Set the path to the head observation output file."""

        self._obs_path = val

    @property
    def obs_heads(self):
        """Observed/simulated heads at the observation locations (loaded and cached)."""

        if self._obs_heads is None:
            self._obs_heads = self.get_obs_heads()
        return self._obs_heads

    @staticmethod
    def sort_dict_by_keys(
            dict_to_sort: dict = None
    ):
        """Sorts the given dict by its keys and returns the sorted dict"""
        return _sort_dict_by_keys(dict_to_sort)

    def get_obs_cells(self, locs: Path, crs: str = None, loc_name_field='ExploName'):
        """Compatibility wrapper for the dedicated observation-cell helper."""

        return _get_obs_cells(self, locs=locs, crs=crs, loc_name_field=loc_name_field)

    def get_obs_heads(
            self, locs: Path | list[int] = None,
            crs: str = None,
            loc_name_field='ExploName',
            long_format=False
    ):
        """Compatibility wrapper for the dedicated observation-head helper."""

        return _get_obs_heads(
            self,
            locs=locs,
            crs=crs,
            loc_name_field=loc_name_field,
            long_format=long_format,
        )



class ConcResults(DependentVariableFile):
    """GWT concentration output reader (``model.conc``).

    The transport twin of :class:`HeadsPlus`: same grammar (``get/summary/map/
    xs/plot/mosaic/animate``), reading the ``.ucn`` concentration binary. ``unit``
    is model-dependent (mass/volume) -- the ``"mg/L"`` default is a convention,
    thread the real unit through when known.
    """

    value_name = "conc"
    store_column = "conc"
    _output_suffix = ".ucn"
    _binary_text = "concentration"
    _summary_label = "conc"
    _choro_type = "conc"
    _default_unit = "mg/L"

    def __init__(self, path: Path = None, model=None, vor: Vor = None, unit: str | None = None):
        super().__init__(path=path, model=model, vor=vor)
        self.unit = unit or self._default_unit

    def _default_hover(self, *, layers: str, surfaces: bool):
        """The default concentration hover spec (unit threaded from the reader)."""

        return conc_hover(unit=self.unit, layers=layers, surfaces=surfaces)

    @property
    def all_conc(self):
        """The full ``kstpkper/layer/cell`` concentration table (alias of ``all_values``)."""

        return self.all_values


class TempResults(DependentVariableFile):
    """GWE temperature output reader (``model.temp``).

    The energy-transport twin of :class:`HeadsPlus`: same grammar, reading the
    temperature binary. ``unit`` defaults to ``"°C"`` by convention.
    """

    value_name = "temp"
    store_column = "temp"
    _output_suffix = ".ucn"
    _binary_text = "temperature"
    _summary_label = "temp"
    _choro_type = "temp"
    _default_unit = "°C"

    def __init__(self, path: Path = None, model=None, vor: Vor = None, unit: str | None = None):
        super().__init__(path=path, model=model, vor=vor)
        self.unit = unit or self._default_unit

    def _default_hover(self, *, layers: str, surfaces: bool):
        """The default temperature hover spec (unit threaded from the reader)."""

        return temp_hover(unit=self.unit, layers=layers, surfaces=surfaces)

    @property
    def all_temp(self):
        """The full ``kstpkper/layer/cell`` temperature table (alias of ``all_values``)."""

        return self.all_values
