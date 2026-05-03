"""Unsaturated-zone-flow helpers for building MF6 UZF package inputs."""

from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor
    from simple_modflow.modflow.mf6.simulation.base import SimulationBase

import flopy
import numpy as np
import pandas as pd
from pathlib import Path
import pickle
from simple_modflow.modflow.mf6.boundary_support import build_cell_id, expand_periodic_cell_input
from simple_modflow.modflow.mf6.recharge import RechargeFromShp
from simple_modflow.modflow.mf6.simulation.packages import _maybe_create_package_artifact


class UZFPackageData:
    """Build UZF packagedata/perioddata and optionally attach the MF6 UZF package."""

    def __init__(
            self,
            model: SimulationBase = None,
            vor: Vor = None,
            idomain: np.ndarray = None,
            vks: list | float = None,
            thtr: list | float = None,
            thts: list | float = None,
            thti: list | float = None,
            eps: float = 4,
            finf: dict | float = None,
            pet: dict | float = None,
            extdp: dict | float = None,
            extwc: dict | float = None,
            ha: dict | float = None,
            hroot: dict | float = None,
            rootact: dict | float = None,
            boundnames: list[str] = None,
            nuzfcells: int = None,
            uzf_cells: list = None,
            aux: list[str] = None,
            add_uzf: bool = True,
            mover: bool = False,
            rch_from_shp: RechargeFromShp = None,
            register_regions: bool = False,
            region_name: str | None = None,
            region_tags: list[str] | None = None,
            overwrite_regions: bool = False,
            artifact_id: str | None = None,
            artifact_catalog=None,
            artifact_description: str | None = None,
            artifact_tags: list[str] | None = None,
            artifact_metadata: dict | None = None,
            artifact_overwrite: bool = False,
    ):
        """
        Initialize and manage the UZF (Unsaturated Zone Flow) package for the model. This
        class handles the setup of model variables required for simulating the unsaturated
        zone, including parameters such as unsaturated zone hydraulic properties, infiltration,
        evapotranspiration, and additional UZF options. It allows for dynamic handling of UZF
        cells and configuration based on model's active domain.

        :param model: SimulationBase instance representing the parent simulation model.
        :param vor: Vor instance defining the vertical discretization of grid cells.
        :param idomain: NumPy array defining the active domain of the model grid where
            UZF cells are located. Will get from model.modelgrid.idomain if not provided.
        :param vks: Vertical unsaturated zone hydraulic conductivity, provided as a value
            or a list of values.
        :param thtr: Residual water content of the unsaturated zone, provided as a value or
            a list of values.
        :param thts: Saturated water content of the unsaturated zone, provided as a value
            or a list of values.
        :param thti: Initial water content of the unsaturated zone, provided as a value or
            a list of values.
        :param eps: Brooks-Corey exponent controlling the capillary effects (default is 4).
        :param finf: Infiltration rate for the unsaturated zone, provided as a dictionary
            or a single value. Will get from rch_from_shp if not provided.
        :param pet: Potential evapotranspiration rate for the unsaturated zone, provided as
            a dictionary or a single value.
        :param extdp: Evapotranspiration extinction depth, provided as a dictionary or a
            single value.
        :param extwc: Evapotranspiration soil wetting curve, provided as a dictionary or a
            single value.
        :param ha: Surface depression storage depth, provided as a dictionary or a single
            value.
        :param hroot: Root zone thickness, provided as a dictionary or single value.
        :param rootact: Root zone activity factor, provided as a dictionary or single
            value.
        :param boundnames: List of names for UZF boundaries.
        :param nuzfcells: Total number of UZF cells; if None, it will be dynamically
            determined from the active domain.
        :param uzf_cells: List of tuples specifying the UZF cells; each tuple contains
            layer and cell indices. If None, UZF cells will be inferred from the active
            domain.
        :param aux: List of auxiliary variable names used by the UZF package.
        :param add_uzf: Boolean indicating whether to add the UZF package to the simulation
            model (default is True).
        :param mover: Boolean indicating whether to include the mover option for UZF flow
            simulation (default is False).
        :param rch_from_shp: RechargeFromShp instance defining recharge rate from a shapefile.
        """
        self.model = model
        if self.model is not None:
            self.vor = model.vor if vor is None else vor
        else:
            self.vor = vor if vor is not None else None
        self.idomain = idomain
        self.vks = vks
        self.thtr = thtr
        self.thts = thts
        self.thti = thti
        self.eps = eps
        self._finf = finf
        self.pet = pet
        self.extdp = extdp
        self.extwc = extwc
        self.ha = ha
        self.hroot = hroot
        self.rootact = rootact
        self.boundnames = boundnames
        self._uzf_vks = None
        self.aux = aux if aux else []
        self.mover = mover

        if uzf_cells is not None:
            self.uzf_cells = self._coerce_uzf_cells(uzf_cells)
        elif self.idomain is None and self.model is not None:
            try:
                idom = model.modelgrid.idomain[0]
                active_cell_indices = pd.Series(idom).loc[idom == 1].index.to_list()
                # assumes all active cells in layer 1 are uzf cells
                self.uzf_cells = self._coerce_uzf_cells(active_cell_indices)
            except ValueError:
                self.uzf_cells = None
        elif self.idomain is not None:
            active_cell_indices = pd.Series(self.idomain).loc[self.idomain == 1].index.tolist()
            # assumes all active cells in layer 1 are uzf cells
            self.uzf_cells = self._coerce_uzf_cells(active_cell_indices)
        elif self.model is not None:
            # assumes all model cells are uzf cells
            self.uzf_cells = self._coerce_uzf_cells(model.cellids)
        else:
            self.uzf_cells = None

        self.nuzfcells = len(self.uzf_cells) if nuzfcells is None else nuzfcells
        self._rch_from_shp = None
        self._rch_dict = None
        self._uzf_packagedata = None
        self._uzf_perioddata = None
        self.register_regions = register_regions
        self.region_name = region_name
        self.region_tags = [] if region_tags is None else list(region_tags)
        self.overwrite_regions = overwrite_regions
        self.artifact_id = artifact_id
        self.artifact_catalog = artifact_catalog
        self.artifact_description = artifact_description
        self.artifact_tags = artifact_tags
        self.artifact_metadata = artifact_metadata
        self.artifact_overwrite = artifact_overwrite

        self.rch_from_shp = rch_from_shp
        self.uzf = None
        self.package_artifact = None

        if self.model is None:
            add_uzf = False
        if add_uzf:
            self.add_uzf()

    @staticmethod
    def _coerce_uzf_cells(cells: list) -> list[tuple[int, int]]:
        normalized = []
        for cell in cells:
            if isinstance(cell, tuple):
                normalized.append((int(cell[0]), int(cell[1])))
            else:
                normalized.append(build_cell_id(cell, grid_type="disv", layer=0))
        return normalized

    def _coerce_cell_parameter(self, value, *, name: str):
        if isinstance(value, (float, int)):
            return [value] * self.nuzfcells
        if isinstance(value, list):
            if len(value) != self.nuzfcells:
                raise ValueError(f"{name} must have {self.nuzfcells} values, got {len(value)}")
            return value
        raise TypeError(f"{name} must be a scalar or list")

    def _period_parameter(self, value, *, default=0.0) -> dict[int, list]:
        return expand_periodic_cell_input(value, nper=self.model.nper, nitems=self.nuzfcells, default=default)

    @property
    def rch_from_shp(self):
        return self._rch_from_shp

    @rch_from_shp.setter
    def rch_from_shp(self, val):
        if val is not None:
            assert isinstance(val, RechargeFromShp), 'rch_from_shp must be a RechargeFromShp instance'
        self._rch_from_shp = val

    @property
    def rch_dict(self):

        if self.rch_from_shp is not None and self._rch_dict is None:
            self._rch_dict = self.rch_from_shp.get_rch()
        return self._rch_dict

    @property
    def finf(self):
        if self._finf is None:
            if self.rch_dict is not None:
                uzfdata = {}
                for per, rows in self.rch_dict.items():
                    by_cell = {row[0]: row[1] for row in rows}
                    uzfdata[per] = [by_cell.get(cell_id, 0.0) for cell_id in self.uzf_cells]
                self._finf = uzfdata
        return self._finf

    @property
    def uzf_vks(self):
        return self._uzf_vks

    @uzf_vks.setter
    def uzf_vks(self, val):
        assert isinstance(val, list), 'uzf_vks must be a list of length number of uzf cells'
        assert len(val) == self.nuzfcells, 'uzf_vks must be a list of length number of uzf cells'
        self._uzf_vks = val

    def get_packagedata(self):
        """Generate UZF packagedata list.

        If vks, thtr, thts, and thti are provided as floats, they are applied uniformly to all UZF cells.
        If provided as lists, they must contain values for each UZF cell.
        """
        vks = self._coerce_cell_parameter(self.vks, name="vks")
        thtr = self._coerce_cell_parameter(self.thtr, name="thtr")
        thts = self._coerce_cell_parameter(self.thts, name="thts")
        thti = self._coerce_cell_parameter(self.thti, name="thti")
        uzf_data = []
        for index, cellid in enumerate(self.uzf_cells):
            uzf_data.append([
                index,  # Feature index
                cellid,  # Cell ID
                1,  # Land flag (1 for surface cells)
                -1,  # Vertical connection index (0 = no connection)
                0.001,  # Surface depression depth
                vks[index],
                # Saturated vertical hydraulic conductivity
                thtr[index],  # Residual water content
                thts[index],  # Saturated water content
                thti[index],  # Initial water content
                self.eps,  # Brooks-Corey exponent
                # None if not self.boundnames else f'UZF_{index}'
            ])
        return uzf_data

    def get_perioddata(self):
        """Generate UZF perioddata dictionary.
         If an argument (e.g., finf, pet, extdp) is provided as a float/int, the value is used for all UZF cells.
        If an argument is provided as a dictionary, it is expected to contain stress periods as keys,
        and lists of values corresponding to UZF cells as values.
        """
        finf = self._period_parameter(self.finf, default=0.0)
        pet = self._period_parameter(self.pet, default=0.0)
        extdp = self._period_parameter(self.extdp, default=0.0)
        extwc = self._period_parameter(self.extwc, default=0.0)
        ha = self._period_parameter(self.ha, default=0.0)
        hroot = self._period_parameter(self.hroot, default=0.0)
        rootact = self._period_parameter(self.rootact, default=0.0)
        period_data = {per: [] for per in range(self.model.nper)}
        bad_inf_cells = []
        for per in period_data.keys():
            for index in range(self.nuzfcells):

                # check for bad finf entries for cells
                finf_val = finf[per][index]
                if isinstance(finf_val, float) and np.isnan(finf_val):
                    if index not in bad_inf_cells:
                        bad_inf_cells.append(index)
                    finf_val = 0.0

                period_data[per].append([
                    index,  # UZF cell index
                    finf_val,
                    pet[per][index],
                    extdp[per][index],
                    extwc[per][index],
                    ha[per][index],
                    hroot[per][index],
                    rootact[per][index],
                ])
        return period_data

    @property
    def perioddata(self):
        if self._uzf_perioddata is None:
            self._uzf_perioddata = self.get_perioddata()
        return self._uzf_perioddata

    @property
    def packagedata(self):
        if self._uzf_packagedata is None:
            self._uzf_packagedata = self.get_packagedata()
        return self._uzf_packagedata

    def add_uzf(self):
        """Add UZF package to MODFLOW 6 model."""
        simulate_et = self.pet is not None or self.extdp is not None or self.extwc is not None

        self.uzf = flopy.mf6.ModflowGwfuzf(
            self.model.gwf,
            save_flows=True,
            print_flows=False,
            pname='uzf',
            nuzfcells=self.nuzfcells,
            packagedata=self.packagedata,
            perioddata=self.perioddata,
            mover=self.mover,
            simulate_et=simulate_et,
            linear_gwet=False,
            square_gwet=False,
            simulate_gwseep=False,
            unsat_etwc=False,
            unsat_etae=False,
            budget_filerecord=f'{self.model.name}_budget.uzf',
            budgetcsv_filerecord=f'{self.model.name}_uzf_budget.csv',
            package_convergence_filerecord=f'{self.model.name}_uzf_package_convergence.csv',
            ntrailwaves=10,
            nwavesets=50,
        )
        self.model._uzf_input = self
        if self.register_regions and self.model is not None:
            region_name = self.region_name or "uzf_cells"
            self.model.add_region_from_cells(
                region_name,
                cellids=self.uzf_cells,
                category="boundary",
                package="uzf",
                tags=self.region_tags or ["uzf"],
                metadata={"nuzfcells": self.nuzfcells},
                overwrite=self.overwrite_regions,
            )
        self.package_artifact = _maybe_create_package_artifact(
            self.model,
            "uzf",
            artifact_id=self.artifact_id,
            artifact_catalog=self.artifact_catalog,
            artifact_description=self.artifact_description,
            artifact_tags=self.artifact_tags,
            artifact_metadata=self.artifact_metadata,
            artifact_overwrite=self.artifact_overwrite,
        )
        return self.uzf


if __name__ == '__main__':

    with open(Path(r"C:\Users\lukem\mf6\cumb_v5g\cumb_v5g.model"), 'rb') as file:
        model = pickle.load(file)

    uzf = UZFPackageData(model)
