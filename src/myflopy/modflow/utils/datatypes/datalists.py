from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase


import numpy as np

import myflopy.modflow.utils.validators as validators


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
        """Wrap a MODFLOW list-data payload, validating ``package`` against the model when given."""

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
        """The default value used where the list has no explicit entry."""

        return self._default

    @default.setter
    def default(self, val):
        """Set the default value."""

        self._default = val

    @property
    def valid_package_names_3ltr(self):
        """The model's valid three-letter package names (for validating ``package``)."""

        valid = validators.valid_package_names_3ltr(self.model)
        return valid

    @property
    def package(self):
        """The MODFLOW package this list belongs to (validated against the model)."""

        return self._package

    @package.setter
    def package(self, val):
        """Set the owning package, validating the name against the model."""

        valid = validators.package_name_validator(self.model, val)
        self._package = valid

    @property
    def data(self):
        """The wrapped list-data payload."""

        return self._data

    @data.setter
    def data(self, val):
        """Set the list-data payload."""

        self._data = val

    @property
    def as_numpy(self):
        """The data as a NumPy array."""

        return np.ndarray(self.data)

    def search_and_return_indices(self, search_for: list) -> list:
        """The first index in the data of each item in ``search_for`` (missing items are skipped)."""

        # Find first occurrence index for each search term
        indices = []
        for thing in search_for:
            if thing in self.data:
                indices.append(self.data.index(thing))
        return indices
