from __future__ import annotations

from pathlib import Path

import pandas as pd

from simple_modflow.modflow.mf6.budget import Budget
from simple_modflow.modflow.mf6.headsplus import HeadsPlus as Hp
from simple_modflow.modflow.utils.datatypes.choros import Choro
from simple_modflow.modflow.utils.datatypes.surface_data import ModelSurface
from simple_modflow.modflow.utils.datatypes.xsections import XSection
from simple_modflow.modflow.utils.inputs import Inputs
from simple_modflow.modflow.utils.outputs import LakOutputData, SFROutputData, UzfOutputData


def get_hds(model):
    model._hds = Hp(model=model, vor=model.vor)
    return model._hds


def get_all_heads(model):
    return get_hds(model).all_heads


def get_surface(model):
    return ModelSurface(model=model)


def build_choro(
    model,
    *,
    kstpkper: tuple = None,
    per: int = None,
    layer: int = 0,
    type: str = 'hds',
    custom_hover: dict = None,
    custom_zs: list = None,
    zmin: float | int = None,
    zmax: float | int = None,
    zoom: int = 13,
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
    **kwargs,
) -> Choro:
    return Choro(
        model=model,
        kstpkper=kstpkper,
        per=per,
        layer=layer,
        type=type,
        custom_hover=custom_hover,
        custom_zs=custom_zs,
        zmin=zmin,
        zmax=zmax,
        zoom=zoom,
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
        **kwargs,
    )


def build_xsection(
    model,
    *,
    per: int = None,
    kstpkper: tuple = None,
    layer: int = 0,
    cells: int | list[int] = None,
    x_or_y: str = None,
    spacing: int = 10,
    num_points: int = 100,
    extrapolate_beyond_section_ends: bool = False,
    interpolate: bool = False,
    use_rbf: bool = False,
    show_model_top=True,
    show_model_btm=False,
):
    return XSection(
        model=model,
        per=per,
        kstpkper=kstpkper,
        layer=layer,
        cells=cells,
        x_or_y=x_or_y,
        spacing=spacing,
        num_points=num_points,
        extrapolate_beyond_section_ends=extrapolate_beyond_section_ends,
        interpolate=interpolate,
        use_rbf=use_rbf,
        show_model_top=show_model_top,
        show_model_btm=show_model_btm,
    )


def get_inputs(model):
    return Inputs(model)


def get_lak_output(model):
    return LakOutputData(model)


def get_uzf_output(model):
    return UzfOutputData(model)


def get_sfr_output(model):
    return SFROutputData(model)


def get_budget(model, package: str = None):
    if package is None:
        return Budget(model)
    return Budget(model, package)


def get_budget_cumulative(model):
    return pd.DataFrame(model.gwf.output.list().get_cumulative())


def get_budget_incremental(model):
    return pd.DataFrame(model.gwf.output.list().get_incremental())


def get_kstpkper(model):
    if model._kstpkper is None:
        model._kstpkper = get_hds(model).kstpkper
    return model._kstpkper
