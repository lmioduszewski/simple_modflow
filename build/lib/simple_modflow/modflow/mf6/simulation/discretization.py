from __future__ import annotations

from typing import TYPE_CHECKING

import flopy

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
    from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor


class DisuGrid:

    def __init__(
            self,
            vor: "Vor",
            model: "SimulationBase",
            top=None,
            bottom=None
    ):
        self.nlay = 1
        if top is None:
            try:
                top = vor.gdf_topbtm["top"].to_list()
            except ValueError:
                print('no top in voronoi grid. cannot find top')
        if bottom is None:
            try:
                bottom = vor.gdf_topbtm["bottom"].to_list()
            except ValueError:
                print('no bottom in voronoi grid. cannot find bottom')
        grid_props = vor.get_disv_gridprops()
        self.disu = flopy.mf6.ModflowGwfdisu(
            model.gwf,
            vertices=grid_props["vertices"],
            cell2d=grid_props["cell2d"],
            length_units="FEET",
            top=top,
            bot=bottom,
            filename=f"{model.name}.disu",
            nvert=len(grid_props["vertices"]),
            nodes=len(grid_props["cell2d"]),
            nja=vor.nja,
            iac=vor.iac,
            ja=vor.ja,
            area=vor.get_cell_areas(),
            ihc=1,
            cl12=vor.cl12,
            hwva=vor.hwva,
            idomain=[1 for _ in range(vor.ncpl)],
        )


class DisvGrid:

    def __init__(
            self,
            vor: "Vor" = None,
            model: "SimulationBase" = None,
            top=None,
            bottom=None,
            nlay=1,
            idomain=None
    ):
        self.nlay = nlay
        model.nlay = nlay
        vor = model.vor if vor is None else vor
        grid_props = vor.get_disv_gridprops()
        self.disv = flopy.mf6.ModflowGwfdisv(
            model.gwf,
            length_units="FEET",
            nlay=nlay,
            ncpl=grid_props['ncpl'],
            nvert=len(grid_props["vertices"]),
            vertices=grid_props['vertices'],
            cell2d=grid_props['cell2d'],
            pname='disv',
            filename=f'{model.name}.disv',
            top=top,
            botm=bottom,
            idomain=idomain
        )


class TemporalDiscretization:

    def __init__(
            self,
            model: "SimulationBase",
            time_units: str = 'DAYS',
            per_len: int = 1,
            period_data: list = None,
            num_steps=10,
            multiplier=1.1
    ):
        nper = model.nper
        if period_data is None:
            period_data = [[per_len, num_steps, multiplier] for _ in range(nper)]
        model.num_steps = num_steps
        model.per_len = per_len
        self.tdis = flopy.mf6.modflow.mftdis.ModflowTdis(
            model.sim,
            pname="tdis",
            time_units=time_units,
            nper=nper,
            perioddata=period_data,
            filename=f"{model.name}.tdis"
        )
