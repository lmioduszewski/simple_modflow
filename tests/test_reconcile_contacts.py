"""Reconciling layer surfaces: what moves, what is kept, and whose numbers it uses.

Five defects of one engine (ledger 377): a frame not indexed ``0..n-1`` had the
wrong rows edited; ``LayerStack`` crushed a pinching unit at or above its own
``min_thickness`` to ``min_sep`` and switched it off; ``which="top"`` left the contact
above a raised surface closer than ``min_sep``; ``adjust_top_btm_overlaps`` reconciled
with the defaults whatever ``min_sep`` it was given and wrote into the caller's frame;
and a per-layer ``min_sep`` crashed the pinch check with ``TypeError``.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import box

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.grid.geometry import adjust_top_btm_overlaps, reconcile_surfaces
from myflopy.surfaces import LayerSurfaces, Surface

Flat = mf.Flat


@pytest.fixture(scope="module")
def vor():
    return irregular_voronoi_grid(nrow=3, ncol=3, cell_size=100.0)


# --- [10] rows are POSITIONS, whatever the frame's index ------------------------- #

def test_a_frame_indexed_from_ten_is_reconciled_by_position():
    """It raised ``IndexError: positional indexers are out-of-bounds``."""

    frame = pd.DataFrame({0: [10.0, 10.0], 1: [5.0, 9.99]}, index=[10, 11])
    out = reconcile_surfaces(None, df=frame)
    assert list(out.index) == [10, 11]
    assert out.loc[10, 1] == 5.0
    assert out.loc[11, 1] == pytest.approx(9.9)


def test_a_reversed_index_edits_the_conflicting_row_and_no_other():
    """Label 1 (5.0, clear of the top) was moved to 9.90 and the conflicting label 0
    (9.99) was left where it was."""

    frame = pd.DataFrame({0: [10.0, 10.0], 1: [5.0, 9.99]}, index=[1, 0])
    out = reconcile_surfaces(None, df=frame)
    assert list(out.index) == [1, 0]
    assert out.loc[1, 1] == 5.0
    assert out.loc[0, 1] == pytest.approx(9.9)


def test_adjust_cells_by_id_reconciles_a_frame_with_cell_ids_as_labels():
    """``adjust_cells_by_id`` selects by LABEL and reconciles: both halves agree now."""

    class _Grid:
        def reconcile_surfaces(self, df, min_sep=0.1, trigger_sep=1, which="bottom"):
            return reconcile_surfaces(self, df=df, min_sep=min_sep,
                                      trigger_sep=trigger_sep, which=which)

    from myflopy.modflow.mf6.grid.geometry import adjust_cells_by_id

    frame = pd.DataFrame({0: [10.0, 10.0, 10.0], 1: [5.0, 5.0, 5.0]}, index=[7, 3, 5])
    out = adjust_cells_by_id(_Grid(), [3], adjustment=-5.0, df=frame, layer=0)
    assert out.loc[3, 1] == pytest.approx(4.9)          # its top fell to 5.0
    assert out.loc[7, 1] == 5.0 and out.loc[5, 1] == 5.0


# --- [11] a LayerStack unit at or above its own minimum is not touched ----------- #

def _clay(vor, bottom, *, min_thickness=0.5, pinch="inactive", split=None, **build):
    return (mf.LayerStack(top=Flat(10.0))
            .add("clay", bottom=Flat(bottom), min_thickness=min_thickness,
                 pinch=pinch, split=split)
            .add("below", bottom=Flat(0.0))
            .build(vor, **build))


@pytest.mark.parametrize("pinch, absent", [("inactive", 0), ("passthrough", -1)])
def test_a_pinching_unit_above_its_minimum_keeps_its_thickness_and_stays_active(vor, pinch, absent):
    """A 0.8 ft clay with min_thickness 0.5 came back 0.1 ft thick and switched off:
    the stack-wide trigger (1.0) was larger than the clay's own minimum."""

    result = _clay(vor, 9.2, pinch=pinch)
    assert np.allclose(result.thickness[0], 0.8)
    assert np.allclose(result.botm[0], 9.2), "the clay's base moved"
    assert set(np.unique(result.idomain[0])) == {1}
    # and a clay genuinely below its minimum still pinches, exactly as before
    thin = _clay(vor, 9.7, pinch=pinch)
    assert np.allclose(thin.thickness[0], 0.1)
    assert set(np.unique(thin.idomain[0])) == {absent}


def test_a_unit_exactly_at_its_minimum_is_present(vor):
    """The pinch rule is STRICTLY below the minimum, so a unit exactly at it is kept,
    and reconcile must not reach it either."""

    result = _clay(vor, 9.5)
    assert np.array_equal(result.botm[0], np.full(result.botm.shape[1], 9.5))
    assert set(np.unique(result.idomain[0])) == {1}


def test_a_split_pinching_unit_above_its_minimum_keeps_every_slice(vor):
    """Each slice is held to its share of the unit minimum: slices of 0.4 against a
    trigger of 0.5 (1.0 / 2) were crushed to [0.1, 0.7]."""

    result = _clay(vor, 9.2, split=2)
    assert np.allclose(result.thickness[:2], 0.4)
    assert set(np.unique(result.idomain[:2])) == {1}


def test_an_explicit_trigger_does_not_switch_off_a_present_pinching_unit(vor):
    """The verdict on a pinching unit is its own ``min_thickness``; ``trigger_sep``
    only decides which thin cells are crushed to ``min_sep``."""

    result = _clay(vor, 9.2, trigger_sep=2.0)
    assert np.allclose(result.thickness[0], 0.8)
    assert set(np.unique(result.idomain[0])) == {1}


def test_a_unit_below_its_minimum_but_beyond_the_trigger_is_pinched_where_it_lies(vor):
    """Below its minimum is the pinch verdict; the crush to ``min_sep`` is reconcile's
    and needs the stack trigger too. A 1.5 ft clay with min 2 is pinched in place."""

    result = _clay(vor, 8.5, min_thickness=2.0)
    assert np.allclose(result.thickness[0], 1.5)
    assert np.allclose(result.botm[0], 8.5)
    assert set(np.unique(result.idomain[0])) == {0}


def test_a_floor_unit_split_unequally_keeps_its_base_and_its_shares(vor):
    """A floor unit's minimum is divided among its slices BY SHARE. Divided equally, a
    2 ft unit at its 2 ft minimum split [0.9, 0.1] had its 0.2 ft slice held at 1.0 ft
    and its base pushed 0.8 ft down."""

    result = (mf.LayerStack(top=Flat(10.0))
              .add("sand", bottom=Flat(8.0), min_thickness=2.0, split=[0.9, 0.1])
              .add("below", bottom=Flat(0.0))
              .build(vor))
    assert np.allclose(result.botm[1], 8.0), "the unit's base moved"
    assert np.allclose(result.thickness[:2, 0], [1.8, 0.2])
    # held open BY SHARE when the unit is thinner than its minimum
    thin = (mf.LayerStack(top=Flat(10.0))
            .add("sand", bottom=Flat(9.0), min_thickness=2.0, split=[0.9, 0.1])
            .add("below", bottom=Flat(0.0))
            .build(vor))
    assert np.allclose(thin.thickness[:2, 0], [1.8, 0.2])


# --- [12] which="top" leaves no contact closer than min_sep ---------------------- #

def test_top_reconcile_rechecks_the_contact_above_a_raised_surface():
    """``[10.0, 9.5, 9.45]`` came back ``[9.6, 9.55, 9.45]``: the top 0.05 above the
    raised bottom."""

    out = reconcile_surfaces(None, df=pd.DataFrame({0: [10.0], 1: [9.5], 2: [9.45]}), which="top")
    assert out.iloc[0].tolist() == pytest.approx([9.65, 9.55, 9.45])


@pytest.mark.parametrize("which", ["bottom", "top"])
def test_every_contact_ends_at_min_sep_or_clear_of_the_trigger(which):
    """The engine's rule, checked on the RESULT for random stacks: a contact is either
    set to exactly ``min_sep`` or was left farther apart than ``trigger_sep``."""

    rng = np.random.default_rng(7)
    frame = pd.DataFrame(np.cumsum(-rng.uniform(-0.5, 1.5, size=(400, 6)), axis=1) + 100.0)
    out = reconcile_surfaces(None, df=frame, min_sep=0.1, trigger_sep=0.5, which=which)
    gaps = -np.diff(out.to_numpy(), axis=1)
    at_sep = np.isclose(gaps, 0.1)
    assert np.all(at_sep | (gaps > 0.5)), np.sort(gaps[~at_sep & (gaps <= 0.5)])[:5]


def test_a_layer_stack_reconciled_from_below_has_no_contact_closer_than_min_sep(vor):
    result = (mf.LayerStack(top=Flat(10.0))
              .add("a", bottom=Flat(9.5), pinch="floor", min_thickness=0.1)
              .add("b", bottom=Flat(9.45), pinch="floor", min_thickness=0.1)
              .add("c", bottom=Flat(0.0))
              .build(vor, reconcile="top", min_sep=0.1))
    assert np.all(result.thickness >= 0.1 - 1e-9), result.thickness[:, 0]


# --- [13] adjust_top_btm_overlaps uses the caller's numbers and frame ----------- #

class _StubGrid:
    """Three cells: cell 0 is selected and its neighbours are cells 1 and 2."""

    def __init__(self, frame):
        self.gdf_topbtm = frame

    def reconcile_surfaces(self, df=None, min_sep=0.1, trigger_sep=1, which="bottom"):
        return reconcile_surfaces(self, df=df, min_sep=min_sep, trigger_sep=trigger_sep, which=which)

    def get_vor_cells_as_series(self, shp):
        return pd.Series([[0]])

    def find_adjacent_cells(self, cell):
        return [1, 2]


def test_the_final_reconcile_uses_the_callers_min_sep():
    """``min_sep=1`` gave 0.1 ft separations: the final reconcile ran on the defaults,
    which found the 1.0 ft gaps the first one opened and closed them to 0.1."""

    grid = _StubGrid(pd.DataFrame({0: [10.0, 10.0, 10.0], 1: [9.5, 5.0, 5.0]}))
    out = adjust_top_btm_overlaps(grid, min_sep=1.0, buffer=0.5)
    assert out[1].tolist() == pytest.approx([9.0, 5.0, 5.0])
    # ... and its trigger_sep: a 0.5 ft gap is not a conflict under trigger 0.25
    loose = adjust_top_btm_overlaps(grid, min_sep=1.0, trigger_sep=0.25, buffer=0.5)
    assert loose[1].tolist() == pytest.approx([9.5, 5.0, 5.0])


def test_the_callers_frame_is_not_written():
    frame = pd.DataFrame({0: [10.0, 10.0, 10.0], 1: [5.0, 5.0, 5.0]})
    before = frame.copy()
    grid = _StubGrid(frame)
    out = adjust_top_btm_overlaps(grid, elev_df=frame, min_sep=1.0)
    pd.testing.assert_frame_equal(frame, before)
    assert out.loc[0, 1] == pytest.approx(9.0)


def test_the_grid_method_forwards_trigger_sep(vor, tmp_path):
    layers = (mf.LayerStack(top=Flat(100.0)).add("a", thickness=20.0)
              .add("b", thickness=20.0).build(vor, attach=True))
    assert layers.nlay == 2
    table = vor.gdf_topbtm.copy()
    centre = vor.gdf_vorPolys.geometry.loc[4].centroid
    gpkg = tmp_path / "one.gpkg"
    gpd.GeoDataFrame(geometry=[box(centre.x - 1, centre.y - 1, centre.x + 1, centre.y + 1)],
                     crs=vor.crs).to_file(gpkg)
    out = vor.adjust_top_btm_overlaps(shp=gpkg, buffer=0.5, min_sep=1.0, trigger_sep=1.0)
    assert out.loc[4, 1] == pytest.approx(99.0)
    # a trigger below the 0.5 ft gap leaves it: the method's own default (1) would
    # have closed it to 1.0, so this is what fails if the method stops forwarding
    loose = vor.adjust_top_btm_overlaps(shp=gpkg, buffer=0.5, min_sep=1.0, trigger_sep=0.25)
    assert loose.loc[4, 1] == pytest.approx(99.5)
    pd.testing.assert_frame_equal(vor.gdf_topbtm, table)


# --- [22] a per-layer min_sep is checked per contact ------------------------------ #

def test_a_per_layer_min_sep_builds_an_idomain(vor):
    """It raised ``TypeError: '<' not supported between instances of 'list' and 'float'``."""

    stack = LayerSurfaces([Surface.flat(100.0), Surface.flat(90.0), Surface.flat(80.0)])
    _top, _botm, idomain = stack.top_botm_idomain(vor, min_sep=[0.1, 0.1], pinch="inactive")
    assert idomain.shape == (2, vor.ncpl) and set(np.unique(idomain)) == {1}
    spec = stack.to_disv(vor, pinch_out=True, min_sep=[0.1, 0.1], attach=False)
    assert spec is not None


def test_a_per_layer_min_sep_is_compared_with_its_own_layers_minimum(vor):
    stack = LayerSurfaces([Surface.flat(100.0), Surface.flat(90.0), Surface.flat(80.0)])
    # layer 0 is a floor layer: exempt, whatever its separation
    stack.top_botm_idomain(vor, min_sep=[2.0, 0.1], minimum_thickness=1.0,
                           pinch=["floor", "inactive"])
    with pytest.raises(ValueError, match=r"layer 1: the fallback min_sep \(2.0\)"):
        stack.top_botm_idomain(vor, min_sep=[0.1, 2.0], minimum_thickness=1.0,
                               pinch=["floor", "inactive"])
