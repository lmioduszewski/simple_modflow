"""Shared accessor helpers used by ``SimulationBase`` and loaded MF6 runs."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from myflopy.modflow.mf6.budget import Budget
from myflopy.modflow.mf6.headsplus import HeadsPlus as Hp
from myflopy.modflow.mf6.package_explorer import ModelPackages
from myflopy.modflow.utils.datatypes.choros import Choro
from myflopy.modflow.utils.datatypes.surface_data import ModelSurface
from myflopy.modflow.utils.datatypes.xsections import XSection
from myflopy.modflow.utils.inputs import Inputs
from myflopy.modflow.utils.outputs import LakOutputData, SFROutputData, UzfOutputData


class ModelOutputs:
    """Namespace for package-specific output helper objects.

    This keeps ``SimulationBase`` package access separate from output-helper
    access. For example, ``model.uzf`` can refer to the MF6 UZF package while
    ``model.outputs.uzf`` exposes convenience output helpers such as
    ``ifno_to_cellid``.
    """

    def __init__(self, model):
        """Bind the package-output accessor namespace (``.lak`` / ``.sfr`` / ``.uzf``) to ``model``."""

        self.model = model

    @property
    def lak(self) -> LakOutputData:
        """Return the lake output helper."""

        return get_lak_output(self.model)

    @property
    def uzf(self) -> UzfOutputData:
        """Return the UZF output helper."""

        return get_uzf_output(self.model)

    @property
    def sfr(self) -> SFROutputData:
        """Return the SFR output helper."""

        return get_sfr_output(self.model)


def get_hds(model):
    """Build and cache the ``HeadsPlus`` helper for ``model``."""
    model._hds = Hp(model=model, vor=model.vor)
    return model._hds


def get_all_heads(model):
    """Return the full heads table from the cached ``HeadsPlus`` helper."""
    return get_hds(model).all_heads


def get_surface(model):
    """Return a ``ModelSurface`` helper for the model."""
    return ModelSurface(model=model)


def build_choro(
    model,
    *,
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
    locs=None,
    rch_scale=None,
    bgs=False,
    hillshade_path: Path = None,
    colorscale: str = None,
    logscale: bool = False,
    contours: bool | str = False,
    contour_values=None,
    contour_levels: int | float | list[float] = 10,
    contour_color: str = "black",
    contour_width: float = 1.5,
    contour_name: str = None,
    contour_clip: bool = True,
    contour_resolution: int = 150,
    contour_method: str = "linear",
    **kwargs,
) -> Choro:
    """Build the standard choropleth wrapper for model result exploration.

    Parameters mirror :class:`~myflopy.modflow.utils.datatypes.choros.Choro`
    so the shared ``SimulationBase.cor(...)`` surface can stay thin.
    """
    return Choro(
        model=model,
        kstpkper=kstpkper,
        per=per,
        per_timestep=per_timestep,
        layer=layer,
        type=type,
        custom_hover=custom_hover,
        custom_zs=custom_zs,
        zmin=zmin,
        zmax=zmax,
        zoom=zoom,
        fit_bounds=fit_bounds,
        bounds_padding=bounds_padding,
        show_layer_elevs=show_layer_elevs,
        show_mounding=show_mounding,
        hover_heads=hover_heads,
        hover_ks=hover_ks,
        locs=locs,
        rch_scale=rch_scale,
        bgs=bgs,
        hillshade_path=hillshade_path,
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
        **kwargs,
    )


def build_xsection(
    model,
    *,
    per: int = None,
    kstpkper: tuple = None,
    layer: int = 0,
    cells: int | list[int] = None,
    line=None,
    x_or_y: str = None,
    spacing: int = 10,
    num_points: int = 100,
    extrapolate_beyond_section_ends: bool = False,
    interpolate: bool = False,
    use_rbf: bool = False,
    show_model_top=True,
    show_model_btm=False,
    animation_kstpkpers=None,
):
    """Build the standard cross-section helper for a model."""
    return XSection(
        model=model,
        per=per,
        kstpkper=kstpkper,
        layer=layer,
        cells=cells,
        line=line,
        x_or_y=x_or_y,
        spacing=spacing,
        num_points=num_points,
        extrapolate_beyond_section_ends=extrapolate_beyond_section_ends,
        interpolate=interpolate,
        use_rbf=use_rbf,
        show_model_top=show_model_top,
        show_model_btm=show_model_btm,
        animation_kstpkpers=animation_kstpkpers,
    )


def get_inputs(model):
    """Return an ``Inputs`` helper for inspecting model inputs."""
    return Inputs(model)


def get_lak_output(model):
    """Return the lake-package output helper."""
    return LakOutputData(model)


def get_uzf_output(model):
    """Return the UZF-package output helper."""
    return UzfOutputData(model)


def get_sfr_output(model):
    """Return the SFR-package output helper."""
    return SFROutputData(model)


def get_outputs(model):
    """Return the grouped output-helper namespace for a model."""

    return ModelOutputs(model)


def get_packages(model):
    """Return the preferred package-exploration namespace for a model."""

    return ModelPackages(model)


def get_budget(model, package: str = None):
    """Return a budget helper, optionally scoped to one package."""
    if package is None:
        return Budget(model)
    return Budget(model, package)


def get_budget_cumulative(model):
    """Return the cumulative listing budget as a DataFrame."""
    return pd.DataFrame(model.gwf.output.list().get_cumulative())


def get_budget_incremental(model):
    """Return the incremental listing budget as a DataFrame."""
    return pd.DataFrame(model.gwf.output.list().get_incremental())


def get_kstpkper(model):
    """Return and cache available ``(kstp, kper)`` combinations."""
    if model._kstpkper is None:
        model._kstpkper = get_hds(model).kstpkper
    return model._kstpkper
