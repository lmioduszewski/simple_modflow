### classes and methods to define and get locations in the model to use in other operations ###

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase

from pathlib import Path

import pandas as pd
from pandas import IndexSlice as idxx


class ObservationLocations:

    def __init__(
            self,
            model: SimulationBase,
            package: str = None,
            locs: Path | int | list = None,
            loc_name_field: str = 'ExploName',
            names: str | list[str] = None
    ):
        """

        :param model:
        :param locs:
        :param loc_name_field:
        :param package:
        :param names:
        """

        self.model = model
        self.vor = self.model.vor
        self._locs = None
        self.loc_name_field = loc_name_field
        self._package = None
        self._names = None

        self.package = package
        self.locs = locs
        self.names = names

        self.valid_cells = list(range(self.model.vor.ncpl))

    @property
    def package(self):
        """The MODFLOW package these locations are associated with."""

        return self._package

    @package.setter
    def package(self, package: str):
        """Set the associated package name."""

        self._package = package

    @property
    def locs(self):
        """The resolved observation cell locations."""

        return self._locs

    @locs.setter
    def locs(self, locs):
        """Set the locations, resolving a path/geometry to grid cells via validation."""

        locs = self._validate_and_return_locs(locs)
        self._locs = locs

    @property
    def names(self):
        """The observation location names."""

        return self._names

    @names.setter
    def names(self, names):
        """Set the observation location names."""

        self._names = names

    def filter_cells_by_pkg(self, pkg: str, cells):
        """Restrict ``cells`` to those active in a package's first-period budget (zero-based)."""

        filter_per = self.model.kstpkper[0]  # use first stress period as basis for filtering
        pkg_cells = self.model.bud(pkg).df.loc[idxx[:, filter_per], :].index.get_level_values(0).to_list()
        pkg_cells = [i - 1 for i in pkg_cells]

    def _validate_and_return_locs(self, locs):
        """Resolve ``locs`` (a shapefile/GeoPackage path or existing cells) to observation cell indices."""


        if isinstance(locs, Path):
            try:
                locs_dict = self.vor.get_vor_cells_as_dict(
                    locs=locs,
                    crs=self.vor.crs,
                    predicate='contains',
                    loc_name_field=self.loc_name_field
                )
                # remove dict entries where the loc was not contained in a cell (outside the grid)
                locs_dict = {key: value for key, value in locs_dict.items() if len(value) > 0}
                # create DataFrame with obs_locs as index and obs cell indices as the values
                locs_df = pd.DataFrame.from_dict(locs_dict).transpose()
                obs_locs = locs_df.index
                return obs_locs

            except ValueError:
                print('location path not readable')
                return None

        elif isinstance(locs, int):
            assert locs in self.valid_cells, 'location index out of range'
            obs_locs = [locs]
            return obs_locs

        elif isinstance(locs, list):
            assert all(loc in self.valid_cells for loc in locs), 'at least one location index not a valid cell index'
            obs_locs = locs
            return obs_locs

        else:
            return None
