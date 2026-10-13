"""A DRN zone's observation route on a structured (DIS) and an unstructured (DISU) model
(ledger 455), and a head target's OBS records beside it.

Both ends assumed DISV: ``DrnFlowTargets.to_flopy_obs`` wrote every record's cellid as
``(layer, cell)``, which MF6 reads as layer/row on DIS and has no form for on DISU, and
the observation route compared a declared cellid with the zone's ``(layer, cell)`` as a
raw tuple, so a correctly declared DIS ``(1, 1, 2)`` or DISU ``(6,)`` was refused as
"declared for" another cell. Both read the cellid through the model's
``CellNumbering`` now (ledger 419), so the records the library writes are the ones MF6
reads, and the route reads them: the zone's flow equals its drains' flow in the budget
file, read by the budget route.

Real MF6 runs of a 2 x 3 x 4 model, FloPy-written: CHD at 100 ft down column 0 of the
top layer, two drains at 90 ft, observed by the OBS package the library writes.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import myflopy as mf

flopy = pytest.importorskip("flopy")
from flopy.utils.gridutil import get_disu_kwargs  # noqa: E402

NLAY, NROW, NCOL = 2, 3, 4
#: The drains, as ``(layer, row, col)``.
DRAINS = [(1, 1, 2), (1, 2, 3)]


def _numbered(kind):
    """Each drain as the views number it, and as the grid's own cellid."""

    if kind == "dis":
        return [(layer, row * NCOL + col) for layer, row, col in DRAINS], list(DRAINS)
    nodes = [layer * NROW * NCOL + row * NCOL + col for layer, row, col in DRAINS]
    return [(0, node) for node in nodes], [(node,) for node in nodes]


@pytest.fixture(scope="module", params=["dis", "disu"])
def run(request, tmp_path_factory):
    kind = request.param
    pairs, cellids = _numbered(kind)
    workspace = tmp_path_factory.mktemp(f"drn_obs_{kind}")
    sim = flopy.mf6.MFSimulation("s", sim_ws=str(workspace), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=2, perioddata=[(1.0, 1, 1.0)] * 2)
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m", save_flows=True)
    if kind == "dis":
        flopy.mf6.ModflowGwfdis(gwf, nlay=NLAY, nrow=NROW, ncol=NCOL, delr=10.0, delc=10.0,
                                top=100.0, botm=[50.0, 0.0])
        chd = [(0, row, 0) for row in range(NROW)]
    else:
        flopy.mf6.ModflowGwfdisu(gwf, **get_disu_kwargs(
            NLAY, NROW, NCOL, np.full(NCOL, 10.0), np.full(NROW, 10.0), 100.0, [50.0, 0.0],
            return_vertices=True))
        chd = [(row * NCOL,) for row in range(NROW)]
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0, icelltype=0)
    flopy.mf6.ModflowGwfic(gwf, strt=100.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[cellid, 100.0] for cellid in chd])
    drn = flopy.mf6.ModflowGwfdrn(gwf, stress_period_data=[[c, 90.0, 5.0] for c in cellids])
    layer = pairs[0][0]
    records = mf.DrnFlowTargets(
        locations=[{"name": "seep", "layer": layer, "cells": [cell for _l, cell in pairs]}]
    ).to_flopy_obs(gwf, csv_name="seep.csv")
    drn.obs.initialize(filename="m.drn.obs", continuous=records)
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord="m.cbc", saverecord=[("BUDGET", "ALL")])
    sim.write_simulation(silent=True)
    success, report = sim.run_simulation(silent=True)
    assert success, "\n".join(report[-20:])
    return kind, records["seep.csv"], mf.load_mf6_run(workspace)


def test_the_records_name_each_drain_by_the_grids_own_cellid(run):
    kind, records, _model = run
    pairs, cellids = _numbered(kind)
    assert [record[2] for record in records] == cellids
    assert [record[0] for record in records] == [f"seep_c{cell}" for _l, cell in pairs]


def test_the_observation_route_reads_the_zone_its_package_declares(run):
    """Read by the observation route (its labels are declared, for the zone's own
    cells), the zone's flow is the budget route's for the same cells."""

    kind, _records, model = run
    pairs, _cellids = _numbered(kind)
    cells = [cell for _l, cell in pairs]
    observed = mf.DrnFlowTargets(locations=[{"name": "seep", "layer": pairs[0][0], "cells": cells}])
    budget = mf.DrnFlowTargets(locations=[{"name": "basin", "layer": pairs[0][0], "cells": cells}])
    flow = budget.simulated_series(model)["basin"].to_numpy()
    assert (flow < 0).all()
    np.testing.assert_allclose(observed.simulated_series(model)["seep"], flow, rtol=1e-10)


def test_the_forward_runs_outputs_only_reader_reads_what_the_loaded_model_reads(run):
    """PEST's forward run reads a zone set through ``RunOutputs`` -- the build's JSON
    record of the model plus the run's output files (ledger 582) -- and the
    build's base table through the loaded model. On a DIS model (a cellid one column
    wider) and a DISU one (one narrower) both routes read the same table either way,
    the declared cellids recorded as the model numbers them."""

    import json

    from myflopy.modflow.mf6.pest.observations import snapshot_run_outputs
    from myflopy.modflow.mf6.pest.run_outputs import RunOutputs

    kind, _records, model = run
    pairs, _cellids = _numbered(kind)
    record = json.loads(json.dumps(snapshot_run_outputs(model, ["drn"])))
    assert record["ncpl"] == ([NROW * NCOL] * NLAY if kind == "dis" else [NLAY * NROW * NCOL])
    assert [row[2] for row in record["packages"]["drn"][0]["obs"][0]["seep.csv"]] == [
        list(pair) for pair in pairs]
    outputs = RunOutputs(model.workspace, record)
    cells = [cell for _l, cell in pairs]
    for name in ("seep", "basin"):  # the observation route, then the budget route
        targets = mf.DrnFlowTargets(locations=[{"name": name, "layer": pairs[0][0], "cells": cells}])
        pd.testing.assert_frame_equal(targets.simulated_series(outputs),
                                      targets.simulated_series(model))


@pytest.mark.parametrize("kind", ["dis", "disu"])
def test_a_head_targets_obs_records_observe_its_own_cell(kind, tmp_path):
    """``HeadTargets.to_flopy_obs`` wrote ``(layer, cell)`` on every grid too (ledger
    455, review): MF6 refused it on DIS, and on DISU read it as node ``layer + 1``
    -- a target at node 18 observed MF6 node 1 (a CHD cell), under Normal termination. Written in the
    grid's own cellid, the OBS file's head is the head file's at the target's cell."""

    sim = flopy.mf6.MFSimulation("s", sim_ws=str(tmp_path), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
    if kind == "dis":
        flopy.mf6.ModflowGwfdis(gwf, nlay=NLAY, nrow=NROW, ncol=NCOL, delr=10.0, delc=10.0,
                                top=100.0, botm=[50.0, 0.0])
        chd = [[(0, row, 0), 100.0] for row in range(NROW)] + [
            [(0, row, NCOL - 1), 90.0] for row in range(NROW)]
        layer, cell, expected = 1, 6, (1, 1, 2)
    else:
        flopy.mf6.ModflowGwfdisu(gwf, **get_disu_kwargs(
            NLAY, NROW, NCOL, np.full(NCOL, 10.0), np.full(NROW, 10.0), 100.0, [50.0, 0.0],
            return_vertices=True))
        chd = [[(row * NCOL,), 100.0] for row in range(NROW)] + [
            [(row * NCOL + NCOL - 1,), 90.0] for row in range(NROW)]
        layer, cell, expected = 0, NROW * NCOL + 6, (NROW * NCOL + 6,)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0, icelltype=0)
    flopy.mf6.ModflowGwfic(gwf, strt=95.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd)
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord="m.hds", saverecord=[("HEAD", "ALL")])
    sim.write_simulation(silent=True)
    assert sim.run_simulation(silent=True)[0]

    model = mf.load_mf6_run(tmp_path)
    targets = mf.HeadTargets(
        locations=pd.DataFrame({"name": ["V"], "layer": [layer], "cell": [cell]}),
        values=pd.DataFrame({"time": [0], "V": [93.0]}))
    assert targets.to_flopy_obs(model)["head_targets.csv"] == [("V", "HEAD", expected)]
    targets.attach_flopy_obs(model)
    model.gwf.simulation.write_simulation(silent=True)
    success, report = model.gwf.simulation.run_simulation(silent=True)
    assert success, "\n".join(report[-20:])

    heads = flopy.utils.HeadFile(str(tmp_path / "m.hds")).get_data()
    truth = float(heads[expected] if kind == "dis" else heads.reshape(-1)[expected[0]])
    observed = pd.read_csv(tmp_path / "head_targets.csv")["V"].to_numpy()
    assert not np.isclose(truth, 100.0)  # node 1, which the old record observed
    np.testing.assert_allclose(observed, [truth], rtol=1e-10)
