from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.mfsimbase import SimulationBase

import numpy as np
from pathlib import Path
import simple_modflow.modflow.utils.validators as validators


def convert_nested_to_int(obj):
    """
    Helper method to convert nested objects to integers
    :param obj: object to convert, typically a list or tuple
    :return:
    """
    if isinstance(obj, (list, tuple)):
        return type(obj)(convert_nested_to_int(x) for x in obj)
    else:
        try:
            return int(obj)
        except:
            return obj


class MfDataList:

    def __init__(
            self,
            data=None,
            default=None,
            model: SimulationBase = None,
            name: str = None,
            package: str = None,
            cell_ids=None,
    ):

        self._default = None
        self.model = model
        self.name = name
        self._cell_ids = None
        self._data = None
        self._package = None

        self.data = data
        self.default = default
        if model is not None:
            self.package = package
        self.cell_ids = cell_ids

    @property
    def default(self):
        return self._default

    @default.setter
    def default(self, val):
        self._default = val

    @property
    def valid_package_names_3ltr(self):
        valid = validators.valid_package_names_3ltr(self.model)
        return valid

    @property
    def package(self):
        return self._package

    @package.setter
    def package(self, val):
        valid = validators.package_name_validator(self.model, val)
        self._package = valid

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, val):
        self._data = val

    @property
    def as_numpy(self):
        return np.ndarray(self.data)

    def search_and_return_indices(self, search_for: list) -> list:
        # Find first occurrence index for each search term
        indices = []
        for thing in search_for:
            if thing in self.data:
                indices.append(self.data.index(thing))
        return indices
