from __future__ import annotations

from collections import OrderedDict
from pathlib import Path
from typing import TYPE_CHECKING

import flopy
from flopy.utils import CellBudgetFile

from myflopy.modflow.utils.datatypes.readers import read_shp_gpkg

if TYPE_CHECKING:
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
    from myflopy.modflow.mf6.simulation.base import SimulationBase


__all__ = ["PRT", "PrtMip", "PrtOc", "PrtPrp", "PrtDisv", "PrtFmi"]


class PRT:
    """Legacy FloPy-backed MODFLOW PRT model helpers.

    These classes are kept for experimental workflows that target the native
    MODFLOW PRT model rather than MP3DU. They are not part of the preferred
    public particle-tracking API because current MODFLOW PRT support is still
    constrained for many DISV grids, especially when cells exceed the practical
    vertex limit accepted by the current implementation.
    """

    def __init__(
        self,
        model: SimulationBase = None,
        name: str = None,
        mf_folder_path: Path = Path().home().joinpath("mf6"),
    ):
        self.model = model
        self.name = self.model.name if name is None else name
        self.model_output_folder_path = mf_folder_path.joinpath(f"{name}")
        self.nper = self.model.gwf.modeltime.nper
        self.bud_file = self.model.model_output_folder_path.joinpath(f"{model.name}.cbc")
        self.hds_file = self.model.model_output_folder_path.joinpath(f"{model.name}.hds")

        self.prt_sim = flopy.mf6.MFSimulation(
            sim_name="sim_prt",
            exe_name="mf6",
            version="mf6",
            sim_ws=self.model_output_folder_path,
        )
        self.tdis = flopy.mf6.modflow.ModflowTdis(
            simulation=self.prt_sim,
            time_units=self.model.gwf.modeltime.time_units,
            nper=1,
            perioddata=[(1, 1, 1)],
        )
        self.prt = flopy.mf6.modflow.ModflowPrt(
            simulation=self.prt_sim,
            modelname=f"{self.name}.prt",
            model_nam_file=f"{self.name}.prt.nam",
            version="mf6",
            exe_name="mf6",
            print_input=True,
            print_flows=False,
            save_flows=True,
        )
        self.ims = flopy.mf6.ModflowIms(
            simulation=self.prt_sim,
            pname="ims",
            complexity="COMPLEX",
        )

    @staticmethod
    def generate_tdis_from_cbc(cbc_path: str) -> list[tuple[float, int, float]]:
        cbc = CellBudgetFile(cbc_path, precision="double")
        kstpkper = cbc.get_kstpkper()
        times = cbc.get_times()

        last_time_by_kper = OrderedDict()
        prev_time = 0.0

        for (kstp, kper), totim in zip(kstpkper, times, strict=False):
            if kper not in last_time_by_kper or kstp > last_time_by_kper[kper][0]:
                dt = totim - prev_time
                last_time_by_kper[kper] = (kstp, dt)
            prev_time = totim

        return [(dt, 1, 1.0) for (kstp, dt) in last_time_by_kper.values()]


class PrtMip:
    def __init__(
        self,
        prt_model: PRT = None,
        porosity: int | list = 0.2,
        retfactor: int = 1,
        izone: int = 0,
    ):
        self.prt_mip = flopy.mf6.modflow.ModflowPrtmip(
            model=prt_model.prt,
            porosity=porosity,
            retfactor=retfactor,
            izone=izone,
            filename=f"prt-{prt_model.name}.mip",
            pname="mip",
        )


class PrtOc:
    def __init__(
        self,
        prt_model: PRT = None,
        save_record=("BUDGET", "LAST"),
        print_record=None,
    ):
        self.prt_oc = flopy.mf6.ModflowPrtoc(
            model=prt_model.prt,
            pname="oc",
            filename=f"prt-{prt_model.name}.oc",
            budget_filerecord=f"prt-{prt_model.name}.cbc",
            track_filerecord=f"prt-{prt_model.name}.trk",
            saverecord=save_record,
            printrecord=print_record,
        )


class PrtPrp:
    def __init__(
        self,
        prt_model: PRT = None,
        exit_solve_tolerance=0.00001,
        stoptime=None,
        stoptraveltime=None,
        istopzone=0,
        shp_gpkg_path: Path = None,
        vor: VoronoiGridPlus = None,
        local_z=0.5,
    ):
        particle_data = read_shp_gpkg(shp_gpkg_path).geometry
        pnts_vor = vor.get_vor_cells_as_series(particle_data).to_list()
        pnt_cells = vor.gdf_vorPolys.loc[pnts_vor].geometry
        nprt = len(pnt_cells)
        packagedata = []
        for i, pnt in enumerate(pnt_cells):
            packagedata.append([i, (0, pnts_vor[i]), pnt.centroid.x, pnt.centroid.y, local_z])

        self.prt_prp = flopy.mf6.modflow.ModflowPrtprp(
            model=prt_model.prt,
            pname="prtprp",
            filename=f"prt-{prt_model.name}.prp",
            print_input=True,
            local_z=True,
            stop_at_weak_sink=True,
            drape=True,
            nreleasepts=nprt,
            packagedata=packagedata,
        )


class PrtDisv:
    def __init__(
        self,
        gwf_model: SimulationBase = None,
        prt_model: PRT = None,
    ):
        modelgrid = gwf_model.modelgrid
        vor = gwf_model.vor
        self.prt_disv = flopy.mf6.modflow.ModflowPrtdisv(
            model=prt_model.prt,
            pname="disv",
            filename=f"prt-{prt_model.name}.disv",
            length_units="feet",
            export_array_ascii=False,
            nlay=modelgrid.nlay,
            ncpl=modelgrid.ncpl,
            nvert=modelgrid.nvert,
            top=modelgrid.top,
            botm=modelgrid.botm,
            idomain=modelgrid.idomain,
            vertices=vor.get_disv_gridprops()["vertices"],
            cell2d=modelgrid.cell2d,
        )


class PrtFmi:
    def __init__(
        self,
        prt_model: PRT = None,
        gwf_model: SimulationBase = None,
    ):
        self.prt_fmi = flopy.mf6.modflow.ModflowPrtfmi(
            save_flows=True,
            model=prt_model.prt,
            filename=f"prt-{prt_model.name}.fmi",
            pname="fmi",
            packagedata=[
                ["GWFBUDGET", prt_model.bud_file.as_posix()],
                ["GWFHEAD", prt_model.hds_file.as_posix()],
            ],
        )
