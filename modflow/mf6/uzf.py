from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
    from .voronoiplus import VoronoiGridPlus as Vor

import flopy
import numpy as np
import pandas as pd
from pathlib import Path
import pickle


class UZFPackageData:
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
    ):
        """
        Initialize the UZF package class.

        If an argument (e.g., finf, pet, extdp) is provided as a float/int, the value is used for all UZF cells.
        If an argument is provided as a dictionary, it is expected to contain stress periods as keys,
        and lists of values corresponding to UZF cells as values.
        """
        print('Initializing UZF package')

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
        self.finf = finf
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
            self.uzf_cells = uzf_cells
        elif self.idomain is None and self.model is not None:
            try:
                idom = model.modelgrid.idomain[0]
                active_cell_indices = pd.Series(idom).loc[idom == 1].index.to_list()
                # assumes all active cells in layer 1 are uzf cells
                self.uzf_cells = list(zip([0 for _ in active_cell_indices], active_cell_indices))
            except ValueError:
                print('cannot get uzf cells from model idomain')
        elif self.idomain is not None:
            active_cell_indices = pd.Series(self.idomain).loc[self.idomain == 1].index.tolist()
            # assumes all active cells in layer 1 are uzf cells
            self.uzf_cells = list(zip([0 for _ in active_cell_indices], active_cell_indices))
        elif self.model is not None:
            # assumes all model cells are uzf cells
            self.uzf_cells = model.cellids
        else:
            self.uzf_cells = None

        self.nuzfcells = len(self.uzf_cells) if nuzfcells is None else nuzfcells
        self._uzf_packagedata = None
        self._uzf_perioddata = None
        self.uzf = None

        if self.model is None:
            add_uzf = False
        if add_uzf:
            print('Adding UZF package')
            self.add_uzf()

    @property
    def uzf_vks(self):
        return self._uzf_vks

    @uzf_vks.setter
    def uzf_vks(self, val):
        assert isinstance(val, list), 'uzf_vks must be a list of length number of uzf cells'
        assert len(val) == self.nuzfcells, 'uzf_vks must be a list of length number of uzf cells'
        self._uzf_vks = vks

    def get_packagedata(self):
        """Generate UZF packagedata list.

        If vks, thtr, thts, and thti are provided as floats, they are applied uniformly to all UZF cells.
        If provided as lists, they must contain values for each UZF cell.
        """
        uzf_data = []
        for index, cellid in enumerate(self.uzf_cells):
            uzf_data.append([
                index,  # Feature index
                cellid,  # Cell ID
                1,  # Land flag (1 for surface cells)
                -1,  # Vertical connection index (0 = no connection)
                0.001,  # Surface depression depth
                self.vks if isinstance(self.vks, (float, int)) else self.vks[index],
                # Saturated vertical hydraulic conductivity
                self.thtr if isinstance(self.thtr, (float, int)) else self.thtr[index],  # Residual water content
                self.thts if isinstance(self.thts, (float, int)) else self.thts[index],  # Saturated water content
                self.thti if isinstance(self.thti, (float, int)) else self.thti[index],  # Initial water content
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
        period_data = {per: [] for per in range(self.model.nper)}
        bad_inf_cells = []
        for per in period_data.keys():
            for index in range(self.nuzfcells):

                # check for bad finf entries for cells
                finf_val = self.finf if isinstance(self.finf, (float, int)) else self.finf.get(
                    per, [0] * self.nuzfcells)[index]
                if isinstance(finf_val, float) and np.isnan(finf_val):
                    if index not in bad_inf_cells:
                        bad_inf_cells.append(index)
                    finf_val = 0.0

                period_data[per].append([
                    index,  # UZF cell index
                    finf_val,
                    self.pet if isinstance(self.pet, (float, int)) else self.pet.get(per, [0] * self.nuzfcells)[
                        index],
                    self.extdp if isinstance(self.extdp, (float, int)) else
                    self.extdp.get(per, [0] * self.nuzfcells)[index],
                    self.extwc if isinstance(self.extwc, (float, int)) else
                    self.extwc.get(per, [0] * self.nuzfcells)[index],
                    self.ha if isinstance(self.ha, (float, int)) else self.ha.get(per, [0] * self.nuzfcells)[
                        index],
                    self.hroot if isinstance(self.hroot, (float, int)) else
                    self.hroot.get(per, [0] * self.nuzfcells)[index],
                    self.rootact if isinstance(self.rootact, (float, int)) else
                    self.rootact.get(per, [0] * self.nuzfcells)[index]
                ])
        if len(bad_inf_cells) > 0:
            print(f'some uzf cell data is nan. Made these cells zero finf. Check these cells....\n'
                  f'{bad_inf_cells}')
        return period_data

    @property
    def perioddata(self):
        if self._uzf_perioddata is None:
            print('Generating UZF perioddata')
            self._uzf_perioddata = self.get_perioddata()
        return self._uzf_perioddata

    @property
    def packagedata(self):
        if self._uzf_packagedata is None:
            print('Generating UZF packagedata')
            self._uzf_packagedata = self.get_packagedata()
        return self._uzf_packagedata

    def add_uzf(self):
        """Add UZF package to MODFLOW 6 model."""
        simulate_et = self.pet is not None or self.extdp is not None or self.extwc is not None
        simulate_et = False

        self.uzf = flopy.mf6.ModflowGwfuzf(
            self.model.gwf,
            save_flows=True,
            print_flows=False,
            pname='uzf',
            nuzfcells=self.nuzfcells,
            packagedata=self.packagedata,
            perioddata=self.perioddata,
            mover=self.mover,
            simulate_et=False,
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
        return self.uzf


if __name__ == '__main__':

    with open(Path(r"C:\Users\lukem\mf6\cumb_v5g\cumb_v5g.model"), 'rb') as file:
        model = pickle.load(file)

    uzf = UZFPackageData(model)
