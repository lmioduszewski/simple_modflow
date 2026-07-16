"""Legacy head-output exploration helpers built on top of FloPy's HeadFile."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:

    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor
    from myflopy.modflow.mf6.simulation.base import SimulationBase

from pathlib import Path

import flopy.utils.binaryfile as bf
import numpy as np
import pandas as pd

from myflopy.modflow.mf6.heads_observations import (
    get_obs_cells as _get_obs_cells,
)
from myflopy.modflow.mf6.heads_observations import (
    get_obs_heads as _get_obs_heads,
)
from myflopy.modflow.mf6.heads_observations import (
    sort_dict_by_keys as _sort_dict_by_keys,
)
from myflopy.modflow.mf6.heads_plotting import (
    choropleth as _choropleth,
)
from myflopy.modflow.mf6.heads_plotting import (
    multimodel_plot_heads as _multimodel_plot_heads,
)
from myflopy.modflow.mf6.heads_plotting import (
    plot_choropleth as _plot_choropleth,
)
from myflopy.modflow.mf6.heads_plotting import (
    plot_heads as _plot_heads,
)
from myflopy.modflow.mf6.package_plotting import SpatialView, _apply_backend
from myflopy.modflow.utils.datatypes.datalists import convert_nested_to_int

idxx = pd.IndexSlice  # for easy index slicing in a MultiIndex DataFrame


def _as_layer_cell_heads(data, *, nlay: int, ncpl: int):
    """Return head data as a consistent ``(nlay, ncpl)`` array."""

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
        "Could not reshape heads to layer/cell form. "
        f"Got shape={values.shape}, expected nlay={nlay}, ncpl={ncpl}."
    )


def multimodel_plot_heads(models: list[SimulationBase], locs: int | list[int] | Path, **kwargs):
    """Compatibility wrapper for the dedicated multi-model plotting helper."""

    return _multimodel_plot_heads(models, locs, **kwargs)




class HeadsPlus(SpatialView, bf.HeadFile):
    """Extended heads-file reader with spatial/model-aware convenience methods.

    Also a full :class:`SpatialView` leaf in the unified grammar --
    ``hds.get()/.summary()`` (tables), ``hds.map()/.plot()/.xs()`` (panels),
    ``hds.mosaic()/.animate()`` (composers). ``SpatialView`` is first in the
    MRO deliberately so the grammar's ``plot`` shadows flopy's legacy
    ``LayerFile.plot``.
    """

    #: value label + tidy-frame column used by the unified grammar
    value_name = "head"

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
        from myflopy.modflow.mf6.simulation.base import SimulationBase
        if hds_path is None:
            if model is None:
                raise ValueError("Must provide heads file or model")
            assert isinstance(model, SimulationBase), 'no valid model provided'
            self.hds_path = model.model_output_folder_path / f'{model.name}.hds'
        else:
            self.hds_path = hds_path

        super().__init__(filename=self.hds_path)

        if model is None:
            self.model = None
        elif isinstance(model, SimulationBase):
            self.model = model
        else:
            raise ValueError("model must be an instance of SimulationBase")

        # Preserve the historical ``.hds`` attribute without opening a second
        # binary reader for the same file.
        self.hds = self
        self.kstpkper = convert_nested_to_int(self.get_kstpkper())
        self.vor = self.model.vor if vor is None else vor
        self.obs_heads_df = None
        self._all_heads = None
        self.nper = pd.DataFrame(self.get_kstpkper()).iloc[:, 1].max() + 1
        self.numstp = pd.DataFrame(self.get_kstpkper()).iloc[:, 0].max() + 1
        self.vor_list = self.vor.gdf_vorPolys.geometry.to_list()
        self.cell_list = [i for i in range(len(self.vor_list))]
        self.area_list = [cell.area for cell in self.vor_list]
        self.x_list = [cell.centroid.xy[0][0] for cell in self.vor_list]
        self.y_list = [cell.centroid.xy[1][0] for cell in self.vor_list]
        self._obs = {}
        self._obs_heads = None
        self.obs_path = obs_path
        self.crs = self.vor.crs

    @property
    def all_heads(self):
        """The full ``kstpkper/layer/cell`` head table, built and cached on first access."""

        if self._all_heads is None:
            self._all_heads = self.get_all_heads()
        return self._all_heads

    # -- unified grammar: data verbs ---------------------------------------
    def get(
        self,
        *,
        per: int | list[int] | None = None,
        layer: int | list[int] | None = None,
        cells: int | list[int] | None = None,
    ) -> pd.DataFrame:
        """Return period-end heads as a tidy ``per``/``layer``/``cell``/``head`` frame.

        Each stress period is reduced to its last saved timestep; flopy dry/no-data
        sentinels (``|head| >= 1e29``) become ``NaN``. All selectors are zero-based.
        """

        frame = self.all_heads.reset_index()
        frame["kstp"] = [int(key[0]) for key in frame["kstpkper"]]
        frame["per"] = [int(key[1]) for key in frame["kstpkper"]]
        period_end = frame.groupby("per")["kstp"].transform("max")
        frame = frame[frame["kstp"] == period_end].copy()
        frame = frame.rename(columns={"elev": "head"})
        frame["head"] = pd.to_numeric(frame["head"], errors="coerce")
        frame.loc[frame["head"].abs() >= 1e29, "head"] = np.nan
        for column, selector in (("per", per), ("layer", layer), ("cells", cells)):
            if selector is None:
                continue
            values = (
                [int(selector)]
                if isinstance(selector, (int, np.integer))
                else [int(value) for value in selector]
            )
            frame = frame[frame["cell" if column == "cells" else column].isin(values)]
        return frame[["per", "layer", "cell", "head"]].reset_index(drop=True)

    def summary(self) -> pd.DataFrame:
        """Return a one-row summary of the saved heads."""

        frame = self.get()
        heads = frame["head"].to_numpy(dtype=float)
        finite = heads[np.isfinite(heads)]
        return pd.DataFrame(
            [
                {
                    "label": "hds",
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
        """Every grid layer -- heads are saved for all layers."""

        return list(range(int(self.nlay)))

    def _spatial_periods(self) -> list[int]:
        """Stress periods present in the saved head output."""

        return sorted({int(key[1]) for key in self.kstpkper})

    def _series_default_agg(self) -> str:
        """Collapse cells within a plotted line by mean (averaging heads, not summing)."""

        return "mean"  # heads average over cells; summing elevations is meaningless

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
            raise ValueError("HeadsPlus.xs() requires a parent model.")
        from myflopy.modflow.utils.datatypes.xsections import XSection

        return {
            self.model.name: XSection(
                model=self.model, line=line, cells=cells, per=per, layer=layer, **kwargs
            )
        }

    def long(self, *, values: str = "elev") -> pd.Series:
        """Return heads as a long series indexed by ``kstpkper/layer/cell``."""

        if values not in self.all_heads.columns:
            raise KeyError(f"Heads value column {values!r} was not found.")
        series = pd.to_numeric(self.all_heads[values], errors="coerce")
        series.name = values
        return series

    def wide(
        self,
        *,
        index: list[str] | tuple[str, ...] = ("layer", "cell"),
        values: str = "elev",
        agg: str = "first",
    ) -> pd.DataFrame:
        """Pivot heads to one row per layer/cell and one column per ``kstpkper``."""

        if values not in self.all_heads.columns:
            raise KeyError(f"Heads value column {values!r} was not found.")
        frame = self.all_heads.reset_index()
        missing_index = [column for column in index if column not in frame.columns]
        if missing_index:
            raise KeyError(f"Heads wide index columns were not found: {missing_index}")
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
        """Return one layer's head field at one time as a 1-D ``(ncpl,)`` array.

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
            A ``(ncpl,)`` float array, one head per Voronoi cell.
        """

        # Use the normalized kstpkper list the class trusts (get_all_heads reads
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

    def to_xugrid(self, *, layers=None, times=None, name: str = "head", masked: bool = True):
        """Export simulated heads across layers and time as an xugrid object.

        Stacks the saved head field into a single ``(time, layer, cell)``
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
            Variable name for the head field (default ``"head"``).
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

        from myflopy._optional import require

        xu = require("xugrid", feature="to_xugrid() unstructured-grid export")
        import xarray as xr  # guaranteed present: xarray is a xugrid dependency

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
        per: int | None = None,
        kstpkper: tuple[int, int] | None = None,
        layer: int = 0,
        contours: bool | str = False,
        contour_levels: int | float | list[float] = 10,
        contour_resolution: int = 150,
        contour_method: str = "linear",
        backend: str = "plotly",
        hover=None,
        hover_layers: str = "active+strip",
        hover_surfaces: bool = False,
        **kwargs,
    ):
        """Return a choropleth map of heads, optionally with contour overlays.

        ``hover`` accepts a :class:`~myflopy.modflow.utils.datatypes.hover.HoverSpec`
        for full control of the sectioned, styled hover; otherwise a default heads
        spec is used. ``hover_layers`` (``"active"`` | ``"active+strip"`` |
        ``"all"``) sets how the vertical head profile shows, and
        ``hover_surfaces=True`` adds the model-top/layer-bottom column (merged with
        the layer table on ``"all"``).
        """

        if self.model is None:
            raise ValueError("HeadsPlus.map() requires a parent model.")
        if hover is None:
            from myflopy.modflow.utils.datatypes.hover import head_hover

            hover = head_hover(layers=hover_layers, surfaces=hover_surfaces)
        choro = self.model.cor(
            per=per,
            kstpkper=kstpkper,
            layer=layer,
            type="hds",
            contours=contours,
            contour_levels=contour_levels,
            contour_resolution=contour_resolution,
            contour_method=contour_method,
            hover_spec=hover,
            **kwargs,
        )
        return _apply_backend(choro, backend)

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

    def get_all_heads(self):
        """Method to get all heads for this model and store in
            a dataframe"""

        vor_cell_list = list(self.vor.gdf_vorPolys.index)

        """set generic MultiIndex for all stress periods and all cells"""
        hds_mdx = pd.MultiIndex.from_product(
            iterables=[
                self.kstpkper,
                list(range(self.nlay)),
                vor_cell_list
            ],
            names=['kstpkper', 'layer', 'cell']
        )
        """Set up a MultiIndex DataFrame to hold the heads
            for all cells and stress periods"""
        df_heads = pd.DataFrame(
            index=hds_mdx,
            columns=['elev']
        )
        """get data for each stress period"""
        for kstpkper in self.kstpkper:
            spHds = _as_layer_cell_heads(
                self.get_data(kstpkper=kstpkper),
                nlay=self.nlay,
                ncpl=len(vor_cell_list),
            )
            """copy and paste this stress period data to the MultiIndex DataFrame"""
            for layer in range(self.nlay):
                assert self.vor.ncpl == spHds.shape[1], (
                    'Are you using the wrong voronoi grid??? \n'
                    f'The provided vor grid has {self.vor.ncpl} cells, but there are {spHds.shape[1]} heads '
                    f'in the model'
                )
                df_heads.loc[idxx[kstpkper, layer, :]] = spHds[layer].reshape(-1, 1)
        return df_heads

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

    def plot_heads(
            self,
            locs: Path | int | list,
            crs: str = None,
            layer: int = 0,
            loc_name_field='ExploName',
            plot_fig: bool = True,
            return_fig: bool = False,
            show_dates: bool = False,
            show_times: bool = False,
            start_period: int = 0,
            loc_names: list = None,
            times: pd.DatetimeIndex = None
    ):
        """Compatibility wrapper for the dedicated heads plotting helper."""

        return _plot_heads(
            self,
            locs,
            crs=crs,
            layer=layer,
            loc_name_field=loc_name_field,
            plot_fig=plot_fig,
            return_fig=return_fig,
            show_dates=show_dates,
            show_times=show_times,
            start_period=start_period,
            loc_names=loc_names,
            times=times,
        )

    def plot_choropleth(self, *args, **kwargs):
        """Compatibility wrapper for the dedicated choropleth plotting helper."""

        return _plot_choropleth(self, *args, **kwargs)

    def choropleth(
            self,
            kstpkper: tuple = (0, 0),
            plot_mounding: bool = False,
            zmin=None,
            zmax=None,
            zoom=13,
            custom_hover: dict = None,
            bottom=None,
            bottom_array=None,
            all_layers: bool = False,
            layer: int = 1,
            obs: Path = None,
            obs_name: str = 'ExploName'
    ):
        """Compatibility wrapper for the dedicated choropleth builder."""

        return _choropleth(
            self,
            kstpkper=kstpkper,
            plot_mounding=plot_mounding,
            zmin=zmin,
            zmax=zmax,
            zoom=zoom,
            custom_hover=custom_hover,
            bottom=bottom,
            bottom_array=bottom_array,
            all_layers=all_layers,
            layer=layer,
            obs=obs,
            obs_name=obs_name,
        )
