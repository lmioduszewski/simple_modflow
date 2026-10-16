"""A DRN zone is read by layer and cell, once per zone (ledger 342).

``DrnFlowTargets`` summed a zone on its budget-file route by matching MF6's
cell-budget NODE -- the model-wide user node, ``layer * ncpl + cell`` on a DISV --
against the zone's IN-LAYER cell numbers, so a zone below layer 0 read whichever
layer-0 drain shared a number (a zone whose drains held -5975.73 read -6.77), and
the PEST zone snapshot dropped ``layer``, so every forward run read layer 0 too.
Beside it: a ``GeoDataFrame`` mixing explicit-cells rows with geometry rows never
resolved the geometry rows; ``to_flopy_obs`` wrote one record per locations ROW;
the observation route's ``compare`` joined the zone table's rows, doubling every
row of a zone given as two; a dict with a ``geometry`` key lost it;
``residuals`` of a set given no values raised ``KeyError: 'per'``; and a partial
miss -- a DRN zone or a head target on no active cell -- went only to a logger
nobody hears (ledger 288).

Real MF6 runs of one small model: the 6 x 6 irregular Voronoi grid of
``test_flow_targets_period_end.py`` (``ncpl`` 36) in TWO confined layers (top 110
ft, bottoms 50 and 0), CHD at 102 ft on the west and east edges of layer 0, heads
filling from 100 ft; two 10-day periods of two steps, ``SAVE HEAD LAST`` and ``SAVE
BUDGET LAST``. Cell 3 of layer 0 is INACTIVE, so every later cell's REDUCED node
number is one less than its USER number. Drains: cells 20, 21 and 27 of LAYER 1
(elevation 60, conductance 50) and a DECOY at cell 20 of LAYER 0 (elevation 99,
conductance 5). An MF6 DRN observation package records zone ``deep`` (those three
layer-1 cells) and zone ``both`` (the decoy plus the three), so a set naming
``deep`` or ``both`` reads the OBSERVATION route and any other name the budget
file. Expected values are read from the cell-budget file with FloPy.
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path
from types import SimpleNamespace

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from flopy.utils import CellBudgetFile
from shapely.geometry import Point

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.pest.forward_run import write_named_series_targets
from tests.test_drn_zones_on_faces import attach_budget

NPER, NSTP, PERLEN = 2, 2, 10.0
ZONE = [20, 21, 27]
DECOY = 20
INACTIVE = 3
#: The two zones the model's own DRN observation package records.
OBSERVED = [
    {"name": "deep", "layer": 1, "cells": ZONE},
    {"name": "both", "layer": 0, "cells": [DECOY]},
    {"name": "both", "layer": 1, "cells": ZONE},
]
#: What ``{"every": [DECOY, 21]}`` -- cells listed with NO layer -- resolves to on this
#: model (both columns active in both layers), declared so the observation route reads
#: it (ledger 550). Cell 21 of layer 0 holds no drain: MF6 writes 3e30.
EVERY = [DECOY, 21]
EVERY_OBSERVED = [{"name": "every", "layer": k, "cells": EVERY} for k in (0, 1)]


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    """The model in the module docstring, run."""

    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = [int(c) for c in np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
            if int(c) != INACTIVE]
    idomain = np.ones((2, grid.ncpl), dtype=int)
    idomain[0, INACTIVE] = 0
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=2, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=110.0,
                botm=[50.0, 0.0], idomain=idomain),
        mf.npf(k=10.0, icelltype=0),
        mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=100.0),
        mf.chd(stress_period_data={0: [[(0, c), 102.0] for c in edge]}),
        mf.drn(stress_period_data={0: [[(1, c), 60.0, 50.0] for c in ZONE]
                                   + [[(0, DECOY), 99.0, 5.0]]},
               observations=mf.DrnFlowTargets(locations=OBSERVED + EVERY_OBSERVED).to_flopy_obs(
                   None, csv_name="drn_flow_targets.csv")),
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=NPER, perioddata=[(PERLEN, NSTP, 1.0)] * NPER)])
    with warnings.catch_warnings():
        # A built grid object cannot be frozen into the run's snapshot, and this run
        # is never reopened.
        warnings.filterwarnings("ignore", "Run 'r' built but its frozen snapshot", UserWarning)
        run = mf.Project(tmp_path_factory.mktemp("drn_layers") / "p", name="p").run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    return run.model()


def _cbc(model):
    """``[(kstpkper, {1-based node: q})]`` of the DRN record, read with FloPy."""

    cbc = CellBudgetFile(Path(model.workspace) / "m.cbc", precision="double")
    try:
        return [
            (kstpkper, {int(n): float(q) for n, q in zip(record["node"], record["q"], strict=True)})
            for kstpkper in cbc.get_kstpkper()
            for record in cbc.get_data(text="DRN", kstpkper=kstpkper)
        ]
    finally:
        cbc.close()


def _raw(model, cells):
    """Each period end's DRN flow summed over ``cells`` -- ``(layer, cell)`` pairs --
    by USER node, straight from the budget file."""

    ncpl = model.vor.ncpl
    nodes = [layer * ncpl + cell + 1 for layer, cell in cells]
    return [sum(flows.get(node, 0.0) for node in nodes) for _kstpkper, flows in _cbc(model)]


DEEP = [(1, c) for c in ZONE]


def _cells(zones, name):
    return zones.set_index("name").loc[name, "cells"]


# -- the node a budget record carries ------------------------------------------------


def test_a_budget_node_is_the_user_node_past_an_inactive_cell(toy):
    """MF6 writes a boundary's USER node, ``layer * ncpl + cell + 1``, whatever the
    idomain removes -- the assumption the budget route's (layer, cell) match rests
    on. Cell 3 of layer 0 is inactive, so a REDUCED number would be one less for
    every drain here: 20, 56, 57, 63."""

    ncpl = toy.vor.ncpl
    assert ncpl == 36
    for _kstpkper, flows in _cbc(toy):
        assert sorted(flows) == [DECOY + 1, *(ncpl + c + 1 for c in ZONE)]
    assert [toy.node_to_lni[n] for n in (DECOY, ncpl + 20, ncpl + 21, ncpl + 27)] == [
        (0, DECOY), (1, 20), (1, 21), (1, 27)]


# -- (layer, cell), on both routes ---------------------------------------------------


def test_a_zone_below_layer_0_reads_its_own_drains_on_both_routes(toy):
    """It read the layer-0 decoy that shares cell number 20 (-7.25 against -5988.33 at
    the end of period 0 here; -6.77 against -5975.73 in the reproducer)."""

    expected = _raw(toy, DEEP)
    assert abs(expected[0]) > 100 * abs(_raw(toy, [(0, DECOY)])[0])  # the decoy is small
    budget = mf.DrnFlowTargets(locations=[{"name": "basin", "layer": 1, "cells": ZONE}])
    np.testing.assert_allclose(budget.simulated_series(toy)["basin"], expected, rtol=1e-12)
    observed = mf.DrnFlowTargets(locations=[{"name": "deep", "layer": 1, "cells": ZONE}])
    series = observed.simulated_series(toy)
    # The observation route: one row per stress period, at its end (ledger 346).
    assert list(series["time"]) == list(range(NPER))
    np.testing.assert_allclose(series["deep"], expected, rtol=1e-10)


def test_a_zone_reaching_two_layers_is_one_zone_with_a_label_per_layer(toy):
    """Its rows are united as (layer, cell) pairs, and cell 20 of each layer is its
    own observation -- the labels carry ``_lay`` so they cannot collide."""

    expected = _raw(toy, [(0, DECOY), *DEEP])
    records = mf.DrnFlowTargets(locations=OBSERVED[1:]).to_flopy_obs(None)["drn_flow_targets.csv"]
    assert records == [
        ("both_lay0_c20", "DRN", (0, 20)), ("both_lay1_c20", "DRN", (1, 20)),
        ("both_lay1_c21", "DRN", (1, 21)), ("both_lay1_c27", "DRN", (1, 27))]
    renamed = [dict(row, name="both_b") for row in OBSERVED[1:]]
    np.testing.assert_allclose(
        mf.DrnFlowTargets(locations=renamed).simulated_series(toy)["both_b"], expected,
        rtol=1e-12)
    series = mf.DrnFlowTargets(locations=OBSERVED[1:]).simulated_series(toy)
    np.testing.assert_allclose(series["both"], expected, rtol=1e-10)


@pytest.mark.slow
@pytest.mark.parametrize("drawn", [False, True, "listed"], ids=["cells", "geometry", "listed"])
def test_the_pest_forward_run_reads_the_zone_in_its_layer(toy, tmp_path, drawn):
    """The build's base table and every forward run's read the zone's own layer.
    The snapshot the forward run rebuilds the zone from carried no ``layer``, so it
    read layer 0 -- a CSV snapshot for explicit cells, a GeoPackage for geometry.
    ``listed``: cells given with NO layer read every layer of their columns (ledger
    550), snapshotted as one row per layer -- the decoy included."""

    pyemu = pytest.importorskip("pyemu")
    expected = _raw(toy, DEEP)
    layers = [1]
    if drawn == "listed":
        locations = {"basin": ZONE}
        expected = _raw(toy, [(0, DECOY), *DEEP])
        layers = [0, 1]
    elif drawn:
        polys = toy.vor.gdf_vorPolys.geometry
        locations = gpd.GeoDataFrame({"name": ["basin"] * 3, "layer": [1] * 3},
                                     geometry=[polys.loc[c] for c in ZONE], crs=toy.vor.crs)
    else:
        locations = [{"name": "basin", "layer": 1, "cells": ZONE}]
    targets = mf.DrnFlowTargets(locations=locations, values=pd.DataFrame(
        {"time": range(NPER), "basin": [-5000.0] * NPER}))
    cal = toy.pest("zone", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    cal.observe(targets, prefix="drn")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build("zone.pst")
    template = Path(cal.template_workspace)
    gpkg = drawn is True
    snapshot = template / f"drn_target_locations.{'gpkg' if gpkg else 'csv'}"
    assert (gpd.read_file(snapshot) if gpkg else pd.read_csv(snapshot))["layer"].tolist() == (
        [1] * 3 if gpkg else layers)

    table = "drn_simulated_drn_flow.csv"
    forward = tmp_path / "forward.csv"
    write_named_series_targets(
        [{"kind": "drn_flow", "locations_file": snapshot,
          "output_csv": forward}],
        sim_ws=toy.workspace, model_name="m")
    for output in (template / table, forward):
        ins = pyemu.pst_utils.InstructionFile(str(template / f"{table}.ins"))
        read = ins.read_output_file(str(output))["obsval"]
        np.testing.assert_allclose(sorted(read, reverse=True), sorted(expected, reverse=True),
                                   rtol=1e-10)


# -- one zone, once ------------------------------------------------------------------


def test_a_geometry_zone_given_as_several_rows_is_written_and_counted_once(toy):
    """Every cell was written once per ROW under one label (4 records for 2 cells),
    and ``n_cells`` counted every row's copy of the zone's union."""

    polys = toy.vor.gdf_vorPolys.geometry
    zone = mf.DrnFlowTargets(locations=gpd.GeoDataFrame(
        {"name": ["g", "g"], "layer": [1, 1]}, geometry=[polys.loc[20], polys.loc[21]],
        crs=toy.vor.crs))
    assert zone.to_flopy_obs(toy)["drn_flow_targets.csv"] == [
        ("g_c20", "DRN", (1, 20)), ("g_c21", "DRN", (1, 21))]
    assert zone.summary(toy).loc[0, "n_cells"] == 2


def test_the_observation_route_compares_a_zone_given_as_two_rows_once(toy):
    """``compare`` joined the zone table's ROWS, so a zone of two rows doubled every
    comparison row. One row per simulated row, as on the budget route."""

    zone = mf.DrnFlowTargets(locations=pd.DataFrame(
        {"name": ["deep", "deep"], "layer": [1, 1], "cells": [[20], [21, 27]]}))
    compare = zone.compare(toy)
    assert len(compare) == len(zone.simulated_series(toy)) == NPER
    assert compare["cells"].iloc[0] == DEEP


def test_a_table_mixing_cells_and_geometry_resolves_the_geometry_rows(toy):
    """The geometry row came back with no cells whenever a ``cells`` column existed."""

    polys = toy.vor.gdf_vorPolys.geometry
    zones = mf.DrnFlowTargets(locations=gpd.GeoDataFrame(
        {"name": ["a", "b"], "layer": [1, 1], "cells": [[20], np.nan]},
        geometry=[None, polys.loc[27]], crs=toy.vor.crs))
    resolved = zones.zone_definitions(toy)
    assert _cells(resolved, "a") == [20] and _cells(resolved, "b") == [27]
    np.testing.assert_allclose(zones.simulated_series(toy)["b"], _raw(toy, [(1, 27)]),
                               rtol=1e-12)
    with pytest.raises(ValueError, match=r"zone\(s\) 'b' are drawn as geometry.*pass the model"):
        zones.zone_definitions(None)


def test_a_row_with_neither_cells_nor_geometry_is_refused_naming_it():
    with pytest.raises(ValueError, match=r"zone\(s\) 'b' give no cells and no geometry"):
        mf.DrnFlowTargets(locations=pd.DataFrame(
            {"name": ["a", "b"], "cells": [[1], None]})).zone_definitions(None)
    with pytest.raises(ValueError, match=r"zone\(s\) 'b' give neither cells nor a geometry"):
        mf.DrnFlowTargets(locations=gpd.GeoDataFrame(
            {"name": ["a", "b"], "cells": [[1], None]},
            geometry=[Point(0.0, 0.0), None])).zone_definitions(None)


def test_a_geometry_key_in_a_dict_or_record_dicts_is_a_geometry_zone(toy):
    """``pd.DataFrame`` kept the key as an ordinary column, which the column filter
    dropped, and the zones were refused as needing a GeoDataFrame."""

    cell27 = toy.vor.gdf_vorPolys.geometry.loc[27]
    for locations in ({"name": ["g"], "layer": [1], "geometry": [cell27]},
                      [{"name": "g", "layer": 1, "geometry": cell27}]):
        assert mf.DrnFlowTargets(locations=locations).zone_definitions(toy)["cells"].tolist() == [[27]]
    with pytest.raises(ValueError, match=r"must be a shapely geometry.*'g' \(str\).*wkt"):
        mf.DrnFlowTargets(locations={"name": ["g"], "geometry": [cell27.wkt]})


def test_a_layer_that_is_not_a_layer_number_is_refused():
    for layer in (1.5, -1, "top"):
        with pytest.raises(ValueError, match="is not a zero-based layer number"):
            mf.DrnFlowTargets(locations=[{"name": "z", "layer": layer, "cells": [1]}])


def test_two_zones_whose_labels_coincide_are_refused():
    """MF6 accepts a duplicate observation name (measured, normal termination), and
    the reader then takes one column for both zones."""

    with pytest.raises(ValueError, match=r"'a_b_c1' from zones 'a b', 'a_b'"):
        mf.DrnFlowTargets(locations=[{"name": "a b", "layer": 0, "cells": [1]},
                                     {"name": "a_b", "layer": 0, "cells": [1]}]).to_flopy_obs(None)


# -- a set given no values -----------------------------------------------------------


def test_the_residuals_of_a_set_given_no_values_are_its_comparison(toy):
    """``compare`` renamed ``per`` to ``time`` on the budget route, and ``residuals``
    then selected ``per``: ``KeyError``."""

    budget = mf.DrnFlowTargets(locations=[{"name": "basin", "layer": 1, "cells": ZONE}])
    residuals = budget.residuals(toy)
    assert residuals["per"].tolist() == residuals["time"].tolist() == list(range(NPER))
    np.testing.assert_allclose(residuals["sim_flow"], _raw(toy, DEEP), rtol=1e-12)
    assert residuals["residual"].isna().all()
    observed = mf.DrnFlowTargets(locations=[{"name": "deep", "layer": 1, "cells": ZONE}])
    assert len(observed.residuals(toy)) == len(observed.simulated_series(toy))


def test_a_budget_node_the_grid_does_not_have_is_refused():
    """A record that does not line up with the model raises rather than matching
    nothing (a stub: a real run cannot write one)."""

    model = attach_budget(
        SimpleNamespace(node_to_lni={0: (0, 0), 1: (0, 1)}, workspace=None, name="m"),
        pd.DataFrame({"kstpkper": [(0, 0)], "node": [5], "q": [-1.0]}))
    with pytest.raises(ValueError, match=r"node\(s\) \[5\].*not cells of this model's grid"):
        mf.DrnFlowTargets(locations={"z": [1]}).simulated_series(model)


@pytest.mark.parametrize(
    ("layer", "cells", "stray"),
    [(1, [20, 99], "[(1, 99)]"), (5, [20], "[(5, 20)]"), (1, [-1, 20, 21], "[(1, -1)]")],
    ids=["cell past the layer", "layer below the bottom", "negative cell"],
)
def test_a_zone_cell_the_grid_does_not_have_is_refused_on_the_budget_route(toy, layer, cells, stray):
    """Review of ledger 342 (ledger 343): an explicit pair the grid does not have
    matched no budget row, so ``[20, 99]`` in layer 1 read cell 20 alone as the zone
    and layer 5 read nothing -- silently, where MF6 refuses the same cellid on the
    observation route. A cell the grid has but ``idomain`` removes is not refused:
    it holds no drain and adds zero."""

    zone = mf.DrnFlowTargets(locations=[{"name": "basin", "layer": layer, "cells": cells}])
    pattern = rf"'basin' {re.escape(stray)}.*2 layer\(s\) of 36 cells"
    with pytest.raises(ValueError, match=pattern):
        zone.simulated_series(toy)
    with pytest.raises(ValueError, match=pattern):
        zone.compare(toy)
    kept = mf.DrnFlowTargets(locations=[{"name": "basin", "layer": 0, "cells": [INACTIVE, DECOY]}])
    np.testing.assert_allclose(kept.simulated_series(toy)["basin"], _raw(toy, [(0, DECOY)]),
                               rtol=1e-12)


# -- a drawn zone reads every layer; a period with no drain reads 0.0 ------------------


def test_a_drawn_zone_with_no_layer_reads_every_layer_of_its_columns(toy):
    """flows-1 (ledger 535): a zone drawn as geometry with no ``layer`` read
    layer 0 alone, in silence -- on the user's model a polygon over a drain stacked in
    four layers read 1.8% of its flow. It reads every layer of the columns it covers
    now, ledger 482's rule for an observation area; a ``layer`` still picks one."""

    polygon = toy.vor.gdf_vorPolys.geometry.loc[DECOY].buffer(-1.0)
    drawn = mf.DrnFlowTargets(locations=gpd.GeoDataFrame(
        {"name": ["z"]}, geometry=[polygon], crs=toy.vor.crs))
    zones = drawn.zone_definitions(toy)
    # Layer 0's cell 3 is inactive; cell 20 is active in both layers.
    assert zones["layer"].tolist() == [0, 1] and zones["cells"].tolist() == [[DECOY], [DECOY]]
    np.testing.assert_allclose(drawn.simulated_series(toy)["z"],
                               _raw(toy, [(0, DECOY), (1, DECOY)]), rtol=1e-12)
    assert drawn.to_flopy_obs(toy)["drn_flow_targets.csv"] == [
        ("z_lay0_c20", "DRN", (0, DECOY)), ("z_lay1_c20", "DRN", (1, DECOY))]
    top = mf.DrnFlowTargets(locations=gpd.GeoDataFrame(
        {"name": ["z"], "layer": [0]}, geometry=[polygon], crs=toy.vor.crs))
    np.testing.assert_allclose(top.simulated_series(toy)["z"], _raw(toy, [(0, DECOY)]),
                               rtol=1e-12)


def test_a_key_past_the_periods_is_refused_on_the_budget_route(toy):
    """Ledger 538: keys [0, 7] on two periods read NaN at 7 and ``stats``
    counted n = 1, with no message (every other file route refused it, ledger 500)."""

    targets = mf.DrnFlowTargets(locations=[{"name": "basin", "layer": 1, "cells": ZONE}],
                                values=pd.DataFrame({"time": [0, 7], "basin": [-1.0, -1.0]}))
    with pytest.raises(ValueError, match=(
            r"DrnFlowTargets\.compare\(\.\.\.\).*the periods the model's cell-budget file "
            r"holds are 0 \.\. 1; 1 target row\(s\) name none of them: 'basin' at time 7")):
        targets.stats(toy)


#: A zone over the layer-1 drains given as a cell list in LAYER 0, where cells 21 and
#: 27 hold no drain in any period.
_WRONG_LAYER = [{"name": "shallow", "layer": 0, "cells": [ZONE[1], ZONE[2]]}]
_REFUSED = (
    r"DrnFlowTargets: zone\(s\) 'shallow' \(layer 0\) hold no drain of package=None "
    r"\(every DRN package: 'drn'\) in any period the model's cell-budget file saved.*"
    r"A zone given a `layer` reads that layer alone: check its `layer` -- or give none, "
    r"to read every layer of its columns"
)


def test_a_zone_whose_cells_never_hold_a_drain_is_refused_on_the_budget_route(toy):
    """Ledger 537, as amended: such a zone came back with no rows, then read 0.0 in
    every period -- the user's drains one layer down, calibrated against a flow no
    parameter moves. It is refused naming the zone, its layer, the packages read and
    that a zone given no layer reads every layer; the same cells in their own layer
    read, and a zone holding a
    drain in only SOME periods still reads 0.0 in the others (tested on the hygiene
    toy, ``test_target_hygiene.py``)."""

    for zones in (_WRONG_LAYER, [{"name": "shallow", "layer": 0, "cells": [ZONE[1]]},
                                 {"name": "basin", "layer": 0, "cells": [DECOY]}]):
        wrong = mf.DrnFlowTargets(locations=zones)
        with pytest.raises(ValueError, match=_REFUSED):
            wrong.simulated_series(toy)
        with pytest.raises(ValueError, match=r"'shallow' \(layer 0\) hold no drain"):
            wrong.compare(toy)
    deep = mf.DrnFlowTargets(locations=[{"name": "shallow", "layer": 1, "cells": ZONE[1:]}])
    np.testing.assert_allclose(deep.simulated_series(toy)["shallow"],
                               _raw(toy, [(1, c) for c in ZONE[1:]]), rtol=1e-12)
    # A zone spanning both layers holds the decoy, so it reads: the zone, not each row.
    spanning = mf.DrnFlowTargets(locations=[{"name": "z", "layer": 0, "cells": [0, 1]},
                                            {"name": "z", "layer": 1, "cells": [ZONE[0]]}])
    np.testing.assert_allclose(spanning.simulated_series(toy)["z"],
                               _raw(toy, [(1, ZONE[0])]), rtol=1e-12)


@pytest.mark.slow
def test_the_pest_build_refuses_a_zone_whose_cells_never_hold_a_drain(toy, tmp_path):
    """The build died on a bare ``IndexError: single positional indexer is
    out-of-bounds``, then built a base table of zeros (ledger 537); it refuses
    the zone by name now, before anything is handed to PEST."""

    wrong = mf.DrnFlowTargets(locations=_WRONG_LAYER,
                              values=pd.DataFrame({"time": [0], "shallow": [-1.0]}))
    cal = toy.pest("dry", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    cal.observe(wrong)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pytest.raises(ValueError, match=_REFUSED):
            cal.build("dry.pst")


# -- a partial miss is heard ---------------------------------------------------------


def test_a_partial_miss_is_a_python_warning_for_zones_and_head_targets(toy):
    """``import myflopy`` installs a NullHandler, so the log line alone reached
    nobody (ledger 288's finding, for the ``.gpkg`` sources)."""

    polys = toy.vor.gdf_vorPolys
    _xmin, _ymin, xmax, ymax = polys.total_bounds
    far = Point(xmax + 5000.0, ymax + 5000.0)
    zones = mf.DrnFlowTargets(locations=gpd.GeoDataFrame(
        {"name": ["in", "out"], "layer": [1, 1]},
        geometry=[polys.geometry.loc[27], far.buffer(10.0)], crs=toy.vor.crs))
    with pytest.warns(UserWarning, match=r"1 of 2 DRN zones land on no active cell.*'out'"):
        zones.zone_definitions(toy)
    wells = mf.HeadTargets(
        locations=gpd.GeoDataFrame({"name": ["IN", "OUT"], "layer": [1, 1]},
                                   geometry=[polys.geometry.loc[27].centroid, far],
                                   crs=toy.vor.crs),
        values=pd.DataFrame({"time": [0], "IN": [100.0], "OUT": [100.0]}))
    with pytest.warns(UserWarning, match=r"1 of 2 head targets land on no active cell.*'OUT'"):
        wells.match_to_model(toy)


# -- a LISTED zone given no layer reads every layer (ledger 550) -------


def test_a_listed_zone_with_no_layer_reads_every_layer_of_its_columns(toy):
    """A ``cells`` list given with no ``layer`` read layer 0 alone: here the decoy,
    -6.77-scale, for a zone whose columns hold the layer-1 drains too (on the user's
    model, 245 of the south drains' 935 records, -33,300.5 ft3/d of -2,613,784.4).
    It reads every layer of its columns now, each layer holding the listed columns
    ACTIVE in it -- exactly as a zone drawn with no layer (ledger 535) -- on the
    budget route, in ``compare``, in ``summary`` and in the labels ``to_flopy_obs``
    writes. A ``layer`` still picks one, and a blank in a ``layer`` column is the
    same as no column."""

    listed = [INACTIVE, DECOY, 21, 27]
    expected = _raw(toy, [(0, DECOY), *DEEP])
    assert not np.allclose(expected, _raw(toy, [(0, DECOY)]))  # what layer 0 read
    pairs = [(0, DECOY), (0, 21), (0, 27), (1, INACTIVE), (1, DECOY), (1, 21), (1, 27)]
    for locations in ({"z": listed},
                      pd.DataFrame({"name": ["z"], "layer": [np.nan], "cells": [listed]}),
                      pd.DataFrame({"name": ["z"], "cells": [",".join(map(str, listed))]})):
        zone = mf.DrnFlowTargets(locations=locations)
        zones = zone.zone_definitions(toy)
        # Layer 0's cell 3 is inactive, so only layer 1 holds it.
        assert zones["layer"].tolist() == [0, 1]
        assert zones["cells"].tolist() == [[DECOY, 21, 27], [INACTIVE, DECOY, 21, 27]]
        np.testing.assert_allclose(zone.simulated_series(toy)["z"], expected, rtol=1e-12)
        compare = zone.compare(toy)
        np.testing.assert_allclose(compare["sim_flow"], expected, rtol=1e-12)
        assert compare["cells"].iloc[0] == pairs
        assert zone.summary(toy).loc[0, "n_cells"] == len(pairs)
        # Without the model the layers are unknown: no cell count, the listed columns
        # counted (it read 0 cells, ledger 579).
        assert zone.summary().loc[0, "n_cells"] is pd.NA
        assert zone.summary().loc[0, "n_columns"] == len(listed)
        assert [record[2] for record in zone.to_flopy_obs(toy)["drn_flow_targets.csv"]] == pairs
        assert zone.to_flopy_obs(toy)["drn_flow_targets.csv"][0][0] == "z_lay0_c20"
    one = mf.DrnFlowTargets(locations=[{"name": "z", "layer": 1, "cells": listed}])
    np.testing.assert_allclose(one.simulated_series(toy)["z"], _raw(toy, DEEP), rtol=1e-12)


def test_the_observation_route_reads_a_listed_zone_in_every_layer(toy):
    """The model declares ``every``'s cells in both layers (``EVERY_OBSERVED``), as
    ``to_flopy_obs`` writes a listed zone given no layer, so ``{"every": [20, 21]}``
    reads MF6's DRN observations -- cell 21 of layer 0 holds no drain, and its 3e30
    reads 0.0 -- and equals the budget file's drains in those cells."""

    zone = mf.DrnFlowTargets(locations={"every": EVERY})
    assert zone._observed(toy)[1] is not None  # the observation route
    expected = _raw(toy, [(k, cell) for k in (0, 1) for cell in EVERY])
    np.testing.assert_allclose(zone.simulated_series(toy)["every"], expected, rtol=1e-10)
    budget = mf.DrnFlowTargets(locations={"other": EVERY}).simulated_series(toy)
    np.testing.assert_allclose(budget["other"], expected, rtol=1e-12)


def test_a_listed_zone_with_no_layer_needs_the_model_and_columns_it_has(toy):
    """The layers are the model's, so it cannot be resolved -- or written as OBS
    records -- with ``model=None``; a listed column the grid does not have, or one
    active in NO layer, is refused naming the zone rather than dropped."""

    zone = mf.DrnFlowTargets(locations={"z": [DECOY]})
    for call in (lambda: zone.zone_definitions(None), lambda: zone.to_flopy_obs(None)):
        with pytest.raises(ValueError, match=(
                r"DrnFlowTargets zone\(s\) 'z' list `cells` with no `layer`, which reads "
                r"EVERY layer.*Pass the model \(it is None here\), or give each row its `layer`")):
            call()
    with pytest.raises(ValueError, match=(
            r"zone\(s\) 'z' \[99\] list cells this model's grid does not have: "
            r"2 layer\(s\) of 36 cells")):
        mf.DrnFlowTargets(locations={"z": [DECOY, 99]}).simulated_series(toy)
    dead = SimpleNamespace(
        name="m", node_to_lni={0: (0, 0), 1: (0, 1), 2: (1, 0), 3: (1, 1)},
        myflopy_context=SimpleNamespace(domain=np.array([[1, 0], [1, 0]])))
    with pytest.raises(ValueError, match=r"zone\(s\) 'z' \[1\] list columns active in NO layer"):
        mf.DrnFlowTargets(locations={"z": [0, 1]}).zone_definitions(dead)
    assert mf.DrnFlowTargets(locations={"z": [0]}).zone_definitions(dead)["cells"].tolist() == [
        [0], [0]]


@pytest.mark.slow
def test_a_one_column_zone_over_an_inactive_layer_runs_forward(tmp_path):
    """A zone LISTING one column with no layer, over a layer where that column is
    INACTIVE, is snapshotted as ``""`` beside ``"27"`` -- and pandas read that CSV
    column as floats, so the forward run parsed ``"27.0"`` and raised on every
    realization after a clean build (ledger 550). The snapshot's
    ``cells`` are read as their text now (``forward_run.read_target_table``), so the
    forward run writes what the base run did, and the IES review parses the same."""

    pyemu = pytest.importorskip("pyemu")
    from myflopy.modflow.mf6.pest.forward_run import read_target_table
    from myflopy.modflow.mf6.pest.ies import _parse_cell_list

    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = [int(c) for c in np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))]
    pinched = 27
    assert pinched not in edge
    idomain = np.ones((2, grid.ncpl), dtype=int)
    idomain[0, pinched] = 0
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=2, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=110.0,
                botm=[50.0, 0.0], idomain=idomain),
        mf.npf(k=10.0, icelltype=0),
        mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=100.0),
        mf.chd(stress_period_data={0: [[(0, c), 102.0] for c in edge]}),
        mf.drn(stress_period_data={0: [[(1, pinched), 60.0, 50.0]]}),
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=NPER, perioddata=[(PERLEN, NSTP, 1.0)] * NPER)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = mf.Project(tmp_path / "p", name="p").run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    built = run.model()
    expected = _raw(built, [(1, pinched)])
    assert all(value < 0 for value in expected)

    targets = mf.DrnFlowTargets(locations={"z": [pinched]}, values=pd.DataFrame(
        {"time": range(NPER), "z": [-5000.0] * NPER}))
    assert targets.zone_definitions(built)["cells"].tolist() == [[], [pinched]]
    cal = built.pest("pinched", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    cal.observe(targets, prefix="drn")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build("pinched.pst")
    template = Path(cal.template_workspace)
    snapshot = template / "drn_target_locations.csv"

    table = "drn_simulated_drn_flow.csv"
    forward = tmp_path / "forward.csv"
    write_named_series_targets(
        [{"kind": "drn_flow", "locations_file": snapshot, "output_csv": forward}],
        sim_ws=built.workspace, model_name="m")
    for output in (template / table, forward):
        ins = pyemu.pst_utils.InstructionFile(str(template / f"{table}.ins"))
        read = ins.read_output_file(str(output))["obsval"]
        np.testing.assert_allclose(sorted(read), sorted(expected), rtol=1e-10)
    # The IES review reads the same snapshot through the same reader.
    cells = read_target_table(snapshot)["cells"].tolist()
    assert cells == ["", str(pinched)]
    assert [_parse_cell_list(value) for value in cells] == [[], [pinched]]


def test_an_obs_declaring_a_layered_zones_one_layer_labels_is_named_not_ignored(toy):
    """``{"deep": ZONE}`` (no layer) reads both layers since ledger 550, so it is
    labelled ``deep_lay<k>_c<cell>``; the model's OBS declares ``deep_c20``/``_c21``/
    ``_c27`` -- the one-layer labels -- so the set reads the budget file, which it did
    without a word (ledger 579)."""

    zone = mf.DrnFlowTargets(locations={"deep": ZONE})
    with pytest.warns(UserWarning, match=(
            r"DrnFlowTargets: zone\(s\) 'deep' reach several layers, so their MF6 observations "
            r"are labelled <name>_lay<layer>_c<cell> .* the DRN package's OBS declares "
            r"DEEP_C20; DEEP_C21; DEEP_C27 -- the one-layer labels written for them before\. "
            r"They are read from the model's cell-budget file instead")):
        flows = zone.simulated_series(toy)
    expected = _raw(toy, [(k, cell) for k in (0, 1) for cell in ZONE])
    np.testing.assert_allclose(flows["deep"], expected, rtol=1e-12)
    # A zone whose labels the OBS does not declare in any form reads the budget in silence.
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        mf.DrnFlowTargets(locations={"other": ZONE}).simulated_series(toy)
        mf.DrnFlowTargets(locations=OBSERVED[:1]).simulated_series(toy)  # the declared route


@pytest.mark.slow
def test_the_forward_run_reads_zone_sets_from_the_runs_outputs_alone(toy, tmp_path, monkeypatch):
    """A realization's DRN zone sets read the run's OUTPUT files through the build's
    record of the inputs they need (ledger 582): handed a directory holding
    the budget file and the OBS output and NO input file, the forward run writes the
    base run's tables -- on the observation route (``deep``, declared) and on the
    budget route (``other``). It loaded the whole model for them before."""

    import json
    import shutil

    from myflopy.modflow.mf6.pest.run_outputs import OUTPUTS_SNAPSHOT

    observed = mf.DrnFlowTargets(locations=OBSERVED[:1], values=pd.DataFrame(
        {"time": range(NPER), "deep": [-5000.0] * NPER}))
    budget = mf.DrnFlowTargets(locations=[{"name": "other", "layer": 1, "cells": ZONE}],
                               values=pd.DataFrame({"time": range(NPER), "other": [-5000.0] * NPER}))
    cal = toy.pest("outputs", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    cal.observe(observed, prefix="obsd")
    cal.observe(budget, prefix="bud")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build("outputs.pst")
    template = Path(cal.template_workspace)
    script = (template / "forward_run.py").read_text()
    assert f"outputs={OUTPUTS_SNAPSHOT!r}" in script
    recorded = json.loads((template / OUTPUTS_SNAPSHOT).read_text())
    assert recorded["budget_file"] == "m.cbc" and recorded["ncpl"] == [36, 36]
    assert [package["name"] for package in recorded["packages"]["drn"]] == ["drn"]

    outputs_only = tmp_path / "outputs_only"
    outputs_only.mkdir()
    for name in ("m.cbc", "drn_flow_targets.csv"):
        shutil.copy(Path(toy.workspace) / name, outputs_only / name)
    sets = [{"kind": "drn_flow", "locations_file": str(template / f"{prefix}_target_locations.csv"),
             "output_csv": str(tmp_path / f"{prefix}.csv")} for prefix in ("obsd", "bud")]
    monkeypatch.chdir(template)
    write_named_series_targets(sets, sim_ws=outputs_only, model_name="m", outputs=OUTPUTS_SNAPSHOT)
    for prefix in ("obsd", "bud"):
        pd.testing.assert_frame_equal(pd.read_csv(tmp_path / f"{prefix}.csv"),
                                      pd.read_csv(template / f"{prefix}_simulated_drn_flow.csv"))
    np.testing.assert_allclose(pd.read_csv(tmp_path / "bud.csv")["other"], _raw(toy, DEEP),
                               rtol=1e-12)
    # Without the record the same directory cannot be read: it holds no model.
    from flopy.mf6.mfbase import MFDataException

    with pytest.raises((FileNotFoundError, MFDataException)):
        write_named_series_targets(sets[1:], sim_ws=outputs_only, model_name="m")
