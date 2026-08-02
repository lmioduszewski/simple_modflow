from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase



def valid_package_names_3ltr(model: SimulationBase) -> list | None:
    """
    Returns a list of valid package names for a modflow model.
    :param model: model object
    :return: list
    """
    if isinstance(model, SimulationBase):
        # `model.package_names` reads flopy's supported `get_package_list()`;
        # `gwf.package_name_dict` has warned since flopy 3.9. It returns the
        # same names uppercased, and every caller of this compares lowercase
        # three-letter types, so the case is normalized here.
        return [str(name).lower() for name in model.package_names]
    else:
        return None


def package_name_validator(model: SimulationBase, package: str):
    """
    Checks that a package name is valid for a given modflow model.
    :param model: modflow model object
    :param package: str of package name, 3 letters (ex. 'drn', 'ghb', 'lak')
    :return: returns package name if valid, otherwise raises error
    """
    valid_names = valid_package_names_3ltr(model)
    assert package in valid_names, \
        f'Package not valid, must be one of {valid_names}'
    return package


def valid_list_of_cell_ints(model: SimulationBase = None, cells: list = None):
    """Assert ``cells`` is a list of ints, each a valid grid cell when ``model`` is given."""

    assert all(isinstance(loc, int) for loc in cells), 'cells must be a list of ints'
    if model:
        valid_cells = model.vor.gdf_vorPolys.index.to_list()
        assert all(loc in valid_cells for loc in cells), 'Not valid cells'

    return True
