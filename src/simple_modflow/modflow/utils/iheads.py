from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.simulation.base import SimulationBase

import pickle
from simple_modflow.modflow.mf6.headsplus import HeadsPlus as Hp
from pandas import IndexSlice as idxx
from pathlib import Path


def get_iheads(model: SimulationBase, kstpkper: tuple = None,
               output_folder: Path = None,):

    hds = Hp(model.model_output_folder_path / f"{model.name}.hds", vor=model.vor)
    heads = []
    if kstpkper is None:
        kstpkper = model.kstpkper[0]
    assert kstpkper in model.kstpkper, f'kstpkper not valid stress period, must be one of {model.kstpkper}'

    for lyr in range(model.modelgrid.nlay):
        lyr_hds = hds.all_heads.loc[idxx[kstpkper, lyr], 'elev'].to_list()
        heads.append(lyr_hds)

    if output_folder is None:
        output_folder = model.model_output_folder_path
    with open(output_folder / 'iheads.hds', 'wb') as file:
        pickle.dump(heads, file=file)

    return heads
