"""`mf.disv` was the one package you could not get from `mf.<name>(...)`.

Every other package in a model comes from a `mf.` helper. The discretization
came from `result.to_disv()` -- a different shape for no reason the user could
see, and the inconsistency that prompted "I don't understand your answer".

`mf.disv(layers)` delegates to that same `to_disv` rather than re-deriving the
geometry. Re-deriving would be a SECOND resolution of the same stack, and
`LayerStack.to_disv`'s own docstring records 8.5 ft of drift from exactly that.
"""

from __future__ import annotations

import numpy as np
import pytest

import myflopy as mf
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.surfaces import Surface


def _grid():
    verts = np.array(
        [[0, 0], [1, 0], [2, 0], [3, 0], [0, 1], [1, 1], [2, 1], [3, 1],
         [0, 2], [1, 2], [2, 2], [3, 2]], dtype=float,
    )
    iverts = [[0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6],
              [4, 5, 9, 8], [5, 6, 10, 9], [6, 7, 11, 10]]
    centers = np.array(
        [[0.5, 0.5], [1.5, 0.5], [2.5, 0.5], [0.5, 1.5], [1.5, 1.5], [2.5, 1.5]]
    )
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=centers)


@pytest.fixture
def built():
    grid = _grid()
    stack = (
        mf.LayerStack(grid, top=Surface.flat(50.0))
        .add("sand", thickness=10.0)
        .add("clay", thickness=10.0)
    )
    return grid, stack.build(grid, attach=True)


def _same_options(a, b) -> bool:
    if set(a.options) != set(b.options):
        return False
    for key, value in a.options.items():
        other = b.options[key]
        if isinstance(value, (np.ndarray, list)):
            if not np.array_equal(np.asarray(value), np.asarray(other), equal_nan=True):
                return False
        elif value != other:
            return False
    return True


def test_mf_disv_of_a_layer_result_equals_its_own_to_disv(built):
    _, layers = built
    assert _same_options(mf.disv(layers), layers.to_disv())


def test_the_mesh_comes_with_the_result(built):
    _, layers = built
    spec = mf.disv(layers)
    assert spec.options["nlay"] == 2
    assert spec.options["ncpl"] == 6


def test_extra_flopy_options_pass_through(built):
    _, layers = built
    assert mf.disv(layers, length_units="METERS").options["length_units"] == "METERS"


def test_the_package_can_be_renamed(built):
    _, layers = built
    assert mf.disv(layers, name="grid").name == "grid"


def test_mixing_a_layer_result_with_hand_written_geometry_is_refused(built):
    """Two descriptions of one mesh; silently preferring either is worse."""

    _, layers = built
    with pytest.raises(TypeError, match="takes the geometry FROM the layer result"):
        mf.disv(layers, nlay=99)


def test_the_hand_written_form_still_works(built):
    grid, layers = built
    props = grid.get_disv_gridprops()
    spec = mf.disv(
        nlay=2, ncpl=props["ncpl"], nvert=len(props["vertices"]),
        vertices=props["vertices"], cell2d=props["cell2d"],
        top=layers.top, botm=layers.botm,
    )
    assert spec.options["nlay"] == 2


def test_calling_with_nothing_names_what_is_missing_and_the_shortcut():
    with pytest.raises(TypeError) as excinfo:
        mf.disv()
    message = str(excinfo.value)
    assert "nlay" in message and "vertices" in message
    assert "mf.disv(stack.build(vor))" in message
