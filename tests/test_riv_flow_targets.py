"""RIV exchange-flow targets: a creek's baseflow as a calibration target (ledger 534).

``RivFlowTargets`` is ``DrnFlowTargets``' river twin: both are thin subclasses of the
one zone implementation in ``observations/_zones.py``. What is RIV's own, and what the
shared reader now does for both families, is pinned here on real MF6 runs of one
small model shaped like the user's two-creek model (a ``riv`` main stem and a
``crisp`` tributary):

The 6 x 6 irregular Voronoi grid of ``test_flow_targets_period_end.py`` (``ncpl`` 36)
in TWO confined layers (top 110 ft, bottoms 50 and 0), CHD at 102 ft on the west and
east edges of layer 0, heads filling from 100 ft; two 10-day periods of two steps,
``SAVE BUDGET ALL`` (so a period's END step must be picked, not its steps added).
RIV package ``riv`` GAINS (stage 95 ft) at cells 14, 15, 16 of layer 0 and at cell
14 of LAYER 1 (stage 90 ft); RIV package ``crisp`` LOSES (stage 106 ft) at cells 20
and 21 of layer 0 in period 0 and at cell 20 alone in period 1. ``crisp`` declares
MF6 ``RIV`` observations for zone ``losers`` (its two cells). Truth is read straight
from the cell-budget file with FloPy, by ``paknam2`` and user node.

A second run carries ``riv`` beside a RIV package ``quiet`` declared WITHOUT
``save_flows``, and DRN packages ``drn_a`` (saved) and ``drn_b`` (not).
"""

from __future__ import annotations

import subprocess
import sys
import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from flopy.utils import CellBudgetFile

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.observations.riv import BoundRivFlowTargets
from myflopy.modflow.mf6.pest import project as pest_project
from myflopy.modflow.mf6.pest.forward_run import write_named_series_targets

NPER, NSTP, PERLEN = 2, 2, 10.0
GAIN = [14, 15, 16]
DEEP = 14
LOSE = [20, 21]
#: crisp's declared zone, in its layer: a zone given no layer reads every layer of its
#: columns (ledger 550), whose labels crisp's OBS does not declare.
LOSERS = [{"name": "losers", "layer": 0, "cells": LOSE}]


def _grid_and_edges():
    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = [int(c) for c in np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))]
    return grid, props, edge


def _disv(grid, props):
    return mf.disv(nlay=2, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                   vertices=props["vertices"], cell2d=props["cell2d"], top=110.0,
                   botm=[50.0, 0.0])


def _run(root, simulation):
    with warnings.catch_warnings():
        # A built grid object cannot be frozen into the run's snapshot, and these
        # runs are never reopened.
        warnings.filterwarnings("ignore", "Run 'r' built but its frozen snapshot", UserWarning)
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    return run.model()


@pytest.fixture(scope="module")
def rivers(tmp_path_factory):
    """The two-creek model of the module docstring, run."""

    grid, props, edge = _grid_and_edges()
    crisp_obs = mf.RivFlowTargets(locations=LOSERS).to_flopy_obs(
        None, csv_name="crisp_obs.csv")
    model = mf.gwf("m", grid=grid, packages=[
        _disv(grid, props),
        mf.npf(k=10.0, icelltype=0),
        mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=100.0),
        mf.chd(stress_period_data={0: [[(0, c), 102.0] for c in edge]}),
        mf.riv.flopy(name="riv", filename="m.riv", save_flows=True, boundnames=True,
                     stress_period_data={0: [[(0, c), 95.0, 50.0, 93.0, "rock"] for c in GAIN]
                                         + [[(1, DEEP), 90.0, 20.0, 60.0, "rock"]]}),
        mf.riv.flopy(name="crisp", filename="m.crisp", save_flows=True, boundnames=True,
                     stress_period_data={0: [[(0, c), 106.0, 50.0, 104.0, "crisp"] for c in LOSE],
                                         1: [[(0, LOSE[0]), 106.0, 50.0, 104.0, "crisp"]]},
                     observations=crisp_obs),
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "LAST"), ("BUDGET", "ALL")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=NPER, perioddata=[(PERLEN, NSTP, 1.0)] * NPER)])
    return _run(tmp_path_factory.mktemp("rivers"), simulation)


def _ends(model):
    """``{(package, user node): [q at each period's END step]}``, read with FloPy."""

    cbc = CellBudgetFile(Path(model.workspace) / "m.cbc", precision="double")
    try:
        flows: dict[tuple[str, int], list[float]] = {}
        for per in range(NPER):
            for index, header in enumerate(cbc.recordarray):
                if header["text"].decode().strip() != "RIV":
                    continue
                if (int(header["kstp"]), int(header["kper"])) != (NSTP, per + 1):
                    continue
                owner = header["paknam2"].decode().strip()
                record = cbc.get_record(index)
                for node, q in zip(record["node"], record["q"], strict=True):
                    flows.setdefault((owner, int(node) - 1), [0.0] * NPER)[per] = float(q)
        return flows
    finally:
        cbc.close()


def _truth(model, cells, packages=("RIV", "CRISP")):
    """Each period end's RIV flow summed over ``cells`` -- ``(layer, cell)`` pairs."""

    ncpl = 36
    flows = _ends(model)
    return [
        sum(flows.get((owner, layer * ncpl + cell), [0.0] * NPER)[per]
            for owner in packages for layer, cell in cells)
        for per in range(NPER)
    ]


# -- MF6's sign and the shape of the run, as the reader assumes them ------------------


def test_the_toy_gains_into_riv_and_loses_from_crisp(rivers):
    """RIV writes aquifer -> river NEGATIVE, river -> aquifer POSITIVE, and one record
    per package per saved time, each named by its package in ``paknam2``."""

    flows = _ends(rivers)
    assert all(flows[("RIV", c)][per] < 0 for c in GAIN for per in range(NPER))
    assert all(flows[("CRISP", c)][0] > 0 for c in LOSE)
    assert flows[("CRISP", LOSE[1])][1] == 0.0  # absent in period 1: no row
    assert flows[("RIV", 36 + DEEP)][0] != 0.0


# -- the budget route ------------------------------------------------------------------


def test_the_budget_route_sums_every_package_and_package_reads_one(rivers):
    """``package=None`` sums every RIV package's records in the zone's cells (ledger 352's
    cells rule); ``package=`` reads that package alone (ledger 534) -- the
    DRN family refused it there, reading a budget table that carried no package."""

    zone = {"z": [GAIN[1], LOSE[0]]}
    both = mf.RivFlowTargets(locations=zone).simulated_series(rivers)
    assert both["time"].tolist() == list(range(NPER))
    np.testing.assert_allclose(both["z"], _truth(rivers, [(0, GAIN[1]), (0, LOSE[0])]),
                               rtol=1e-12)
    for name in ("crisp", "CRISP"):
        alone = mf.RivFlowTargets(locations=zone, package=name).simulated_series(rivers)
        np.testing.assert_allclose(alone["z"], _truth(rivers, [(0, LOSE[0])], ("CRISP",)),
                                   rtol=1e-12)
    main = mf.RivFlowTargets(locations=zone, package="riv").simulated_series(rivers)
    np.testing.assert_allclose(main["z"], _truth(rivers, [(0, GAIN[1])], ("RIV",)), rtol=1e-12)
    with pytest.raises(ValueError, match=r"RivFlowTargets: package='rock' names no RIV package"):
        mf.RivFlowTargets(locations=zone, package="rock").simulated_series(rivers)


def test_a_period_is_its_end_step_not_the_sum_of_its_steps(rivers):
    """``SAVE BUDGET ALL`` writes two steps a period; the zone reads the second (336)."""

    cbc = CellBudgetFile(Path(rivers.workspace) / "m.cbc", precision="double")
    try:
        first = sum(float(q) for record in cbc.get_data(text="RIV", kstpkper=(0, 0))
                    for node, q in zip(record["node"], record["q"], strict=True)
                    if int(node) - 1 in GAIN)
    finally:
        cbc.close()
    series = mf.RivFlowTargets(locations=[{"name": "g", "layer": 0, "cells": GAIN}],
                               package="riv").simulated_series(rivers)
    end = _truth(rivers, [(0, c) for c in GAIN], ("RIV",))
    np.testing.assert_allclose(series["g"], end, rtol=1e-12)
    assert not np.isclose(series["g"].iloc[0], first + end[0])


def test_a_period_with_no_river_in_the_zone_reads_zero(rivers):
    """flows-2 (ledger 537): crisp's cell 21 is absent in period 1. The zone
    read NaN there -- an empty PEST field that failed every realization's instruction
    read -- where the observation route reads 0.0."""

    series = mf.RivFlowTargets(locations={"c21": [LOSE[1]]}).simulated_series(rivers)
    assert series["c21"].iloc[0] > 0.0
    assert series["c21"].tolist()[1:] == [0.0]


def test_a_zone_whose_cells_never_hold_a_river_of_the_packages_read_is_refused(rivers):
    """Ledger 537, as amended: a zone over cells that never hold a river of
    the packages read read 0.0 in every period (and before that, built nothing and
    died on a bare IndexError) -- a zone in the wrong layer or package, calibrated
    against a flow no parameter moves. It is refused naming the zone, its layer(s)
    and ``package=``: with a package named, the hint says ``package=None`` reads every
    RIV package."""

    with pytest.raises(ValueError, match=(
            r"RivFlowTargets: zone\(s\) 'dry' \(layers 0, 1\) hold no river boundary of "
            r"package=None \(every RIV package: 'riv', 'crisp'\) in any period the "
            r"model's cell-budget file saved.*or give none, to read every layer")):
        mf.RivFlowTargets(locations={"dry": [0, 1], "trib": LOSE}).simulated_series(rivers)
    # Not ``losers``: crisp's OBS declares that name, which would read the other route.
    wrong_package = mf.RivFlowTargets(locations={"trib": LOSE}, package="riv",
                                      values={"time": [0], "trib": [100.0]})
    for read in (wrong_package.simulated_series, wrong_package.stats):
        with pytest.raises(ValueError, match=(
                r"zone\(s\) 'trib' \(layers 0, 1\) hold no river boundary of package='riv' "
                r"\(RIV package 'riv'\).*package=None reads every RIV package\.")):
            read(rivers)
    # riv's own layer-1 boundary reads; crisp holds none in layer 1.
    deep = mf.RivFlowTargets(locations=[{"name": "d", "layer": 1, "cells": [DEEP]}])
    np.testing.assert_allclose(deep.simulated_series(rivers)["d"],
                               _truth(rivers, [(1, DEEP)], ("RIV",)), rtol=1e-12)
    with pytest.raises(ValueError, match=r"'d' \(layer 1\) hold no river boundary of package='crisp'"):
        mf.RivFlowTargets(locations=[{"name": "d", "layer": 1, "cells": [DEEP]}],
                          package="crisp").simulated_series(rivers)


def test_a_zone_reads_its_own_layer_and_one_given_no_layer_every_layer(rivers):
    """Ledger 342: a cell is a (layer, cell) pair. A zone given no layer reads every
    layer of its columns (ledger 482's area rule), whether DRAWN as geometry (ledger
    535) or LISTED as cells (ledger 550) -- both read layer 0
    alone, so the listed zone here read riv's layer-0 cell and missed its layer-1 one."""

    top = _truth(rivers, [(0, DEEP)], ("RIV",))
    deep = _truth(rivers, [(1, DEEP)], ("RIV",))
    assert np.all(np.abs(deep) > 0) and np.all(np.abs(top) > 0)
    for listed in ({"z": [DEEP]}, pd.DataFrame({"name": ["z"], "cells": [[DEEP]],
                                                  "layer": [np.nan]})):
        explicit = mf.RivFlowTargets(locations=listed).simulated_series(rivers)
        np.testing.assert_allclose(explicit["z"], np.add(top, deep), rtol=1e-12)
    zero = mf.RivFlowTargets(locations=[{"name": "z", "layer": 0, "cells": [DEEP]}])
    np.testing.assert_allclose(zero.simulated_series(rivers)["z"], top, rtol=1e-12)
    layered = mf.RivFlowTargets(locations=[{"name": "z", "layer": 1, "cells": [DEEP]}])
    np.testing.assert_allclose(layered.simulated_series(rivers)["z"], deep, rtol=1e-12)

    polygon = rivers.vor.gdf_vorPolys.geometry.loc[DEEP].buffer(-1.0)
    drawn = mf.RivFlowTargets(locations=gpd.GeoDataFrame(
        {"name": ["z"]}, geometry=[polygon], crs=rivers.vor.crs))
    zones = drawn.zone_definitions(rivers)
    assert zones["layer"].tolist() == [0, 1] and zones["cells"].tolist() == [[DEEP], [DEEP]]
    np.testing.assert_allclose(drawn.simulated_series(rivers)["z"],
                               np.add(top, deep), rtol=1e-12)
    one = mf.RivFlowTargets(locations=gpd.GeoDataFrame(
        {"name": ["z"], "layer": [1]}, geometry=[polygon], crs=rivers.vor.crs))
    np.testing.assert_allclose(one.simulated_series(rivers)["z"], deep, rtol=1e-12)


def test_either_sign_is_a_target_and_the_residual_is_simulated_minus_target(rivers):
    """No sign is refused for RIV (a losing reach is physical); DRN still refuses a
    positive discharge (ledger 350)."""

    losing = _truth(rivers, [(0, c) for c in LOSE], ("CRISP",))
    targets = mf.RivFlowTargets(locations={"losers": LOSE}, package="crisp",
                                values={"time": [0, 1], "losers": [100.0, -5.0]})
    compare = targets.compare(rivers)
    np.testing.assert_allclose(compare["residual"], np.subtract(losing, [100.0, -5.0]),
                               rtol=1e-10)
    with pytest.raises(ValueError, match=r"Enter a measured discharge Q as -Q"):
        mf.DrnFlowTargets(locations={"s": [1]}, values={"time": [0], "s": [100.0]})


def test_a_key_past_the_periods_is_refused_on_the_budget_route(rivers):
    """It read NaN and ``stats()`` dropped it in silence (ledger 538), while
    every other file route and the PEST build refused it (ledger 500)."""

    targets = mf.RivFlowTargets(locations={"g": GAIN},
                                values={"time": [0, 7], "g": [-1.0, -1.0]})
    with pytest.raises(ValueError, match=(
            r"RivFlowTargets\.compare\(\.\.\.\): a target is keyed by zero-based stress "
            r"period.*the model's cell-budget file holds are 0 \.\. 1.*'g' at time 7")):
        targets.stats(rivers)


# -- the observation route --------------------------------------------------------------


def test_the_observation_route_reads_crisps_declared_riv_observations(rivers):
    """Ledger 344's declared identity with MF6's ``RIV`` type: the zone ``losers`` is
    declared on ``crisp``, so it reads MF6's CSV -- one row per period at its end (346)
    -- and equals the budget route to MF6's printed digits; cell 21's period-1 3e30
    (no river there) reads 0.0. Without ``package=`` the two RIV packages are refused
    (ledger 352)."""

    observed = mf.RivFlowTargets(locations=LOSERS, package="crisp")
    series = observed.simulated_series(rivers)
    raw = pd.read_csv(Path(rivers.workspace) / "crisp_obs.csv")
    assert len(raw) == NPER * NSTP and series["time"].tolist() == list(range(NPER))
    np.testing.assert_allclose(series["losers"], _truth(rivers, [(0, c) for c in LOSE], ("CRISP",)),
                               rtol=1e-10)
    budget = mf.RivFlowTargets(locations={"other": LOSE}, package="crisp").simulated_series(rivers)
    np.testing.assert_allclose(series["losers"], budget["other"], rtol=1e-10)
    with pytest.raises(ValueError, match=r"2 RIV packages \('riv', 'crisp'\).*ONE package's "
                                         r"river boundary in its cell.*package=<its name>"):
        mf.RivFlowTargets(locations=LOSERS).simulated_series(rivers)
    with pytest.raises(ValueError, match=r"RivFlowTargets writes RIV observations.*kind='DRN'"):
        mf.RivFlowTargets(locations=LOSERS).to_flopy_obs(None, kind="DRN")


def test_attach_flopy_obs_goes_on_the_named_package(rivers, monkeypatch):
    """Found by TYPE then NAME (ledger 352), never the first RIV package; the shared
    model is left untouched (FloPy's OBS constructor is stood in for)."""

    calls = {}

    def _fake_obs(owner, pname, continuous, filename):
        calls.update(owner=owner, pname=pname, continuous=continuous, filename=filename)

    monkeypatch.setattr("flopy.mf6.modflow.mfutlobs.ModflowUtlobs", _fake_obs)
    mf.RivFlowTargets(locations={"z": [LOSE[0]]}, package="crisp").attach_flopy_obs(rivers)
    assert calls["owner"] is rivers.package("crisp")
    assert calls["pname"] == "riv_flow_obs" and calls["filename"] == "riv_flow_targets.obs"
    # Given no layer, the zone is cell 20 of EVERY layer (ledger 550).
    assert calls["continuous"] == {"riv_flow_targets.csv": [
        ("z_lay0_c20", "RIV", (0, LOSE[0])), ("z_lay1_c20", "RIV", (1, LOSE[0]))]}
    with pytest.raises(ValueError, match=r"2 RIV packages.*ONE of them"):
        mf.RivFlowTargets(locations={"z": [LOSE[0]]}).attach_flopy_obs(rivers)


# -- SAVE_FLOWS -------------------------------------------------------------------------


@pytest.fixture(scope="module")
def unsaved(tmp_path_factory):
    """``riv`` and ``drn_a`` saved, ``quiet`` (RIV) and ``drn_b`` declared without SAVE_FLOWS."""

    grid, props, edge = _grid_and_edges()
    model = mf.gwf("m", grid=grid, packages=[
        _disv(grid, props), mf.npf(k=10.0, icelltype=0), mf.ic(strt=100.0),
        mf.chd(stress_period_data={0: [[(0, c), 102.0] for c in edge]}),
        mf.riv.flopy(name="riv", filename="m.riv", save_flows=True,
                     stress_period_data={0: [[(0, 14), 95.0, 50.0, 93.0]]}),
        mf.riv.flopy(name="quiet", filename="m.quiet", save_flows=False,
                     stress_period_data={0: [[(0, 15), 95.0, 50.0, 93.0]]}),
        mf.drn.flopy(name="drn_a", filename="m.drn_a", save_flows=True,
                     stress_period_data={0: [[(0, 20), 95.0, 50.0]]}),
        mf.drn.flopy(name="drn_b", filename="m.drn_b", save_flows=False,
                     stress_period_data={0: [[(0, 21), 95.0, 50.0]]}),
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    return _run(tmp_path_factory.mktemp("unsaved"), simulation)


@pytest.mark.parametrize(
    ("family", "cells", "missing", "saved"),
    [(mf.RivFlowTargets, [14, 15], "quiet", "riv"), (mf.DrnFlowTargets, [20, 21], "drn_b", "drn_a")],
    ids=["riv", "drn"],
)
def test_a_package_that_saved_no_flows_is_refused_naming_it(unsaved, family, cells, missing, saved):
    """Ledger 536: a package declared without SAVE_FLOWS writes no record,
    and its boundaries were left out of the zone in silence -- the DRN zone over both
    read -339.52 where both drains removed about 679. The named saved package reads."""

    with pytest.raises(ValueError, match=(
            rf"{family.__name__}: [A-Z]+ package\(s\) '{missing}' wrote no record to the model's "
            r"cell-budget file.*save_flows=True, or name the package")):
        family(locations={"z": cells}).simulated_series(unsaved)
    with pytest.raises(ValueError, match=rf"'{missing}' wrote no record.*save_flows=True\."):
        family(locations={"z": cells}, package=missing).simulated_series(unsaved)
    series = family(locations={"z": cells}, package=saved).simulated_series(unsaved)
    assert series["z"].iloc[0] < 0.0


# -- registry and PEST --------------------------------------------------------------------


def test_model_targets_riv_flow_binds_the_bound_helper(rivers):
    targets = mf.RivFlowTargets(locations=[{"name": "g", "layer": 0, "cells": GAIN}],
                                package="riv")
    rivers.targets.riv_flow = targets
    bound = rivers.targets.riv_flow
    assert isinstance(bound, BoundRivFlowTargets) and bound.targets is targets
    np.testing.assert_allclose(bound.simulated_series()["g"],
                               _truth(rivers, [(0, c) for c in GAIN], ("RIV",)), rtol=1e-12)
    assert bound.to_flopy_obs()["riv_flow_targets.csv"][0] == ("g_c14", "RIV", (0, 14))
    del rivers.targets["riv_flow"]


def test_the_pest_build_hands_pest_the_zone_and_every_realization_rewrites_it(rivers, tmp_path):
    """``cal.observe`` wraps a ``RivFlowTargets`` in its spec (removing it from
    ``_TARGET_SPEC_TYPES`` fails here by name); the build's base table and a forward run
    write the same bytes, through ONE ``write_named_series_targets`` call for every
    named-series set (ledger 540); a zone's group and weight reach the
    control file. A package name read from a table is often ``np.str_``, and the
    forward run's entry is written by ``repr``: it is written as a plain ``str``,
    never ``np.str_('crisp')``."""

    assert any(target is mf.RivFlowTargets for target, _spec in pest_project._TARGET_SPEC_TYPES)
    gaining = mf.RivFlowTargets(
        locations=pd.DataFrame({"name": ["rock"], "layer": [0], "cells": [GAIN],
                                "group": ["baseflow"], "weight": [0.01]}),
        values={"time": [0, 1], "rock": [-900.0, -950.0]}, package="riv")
    losing = mf.RivFlowTargets(locations=LOSERS, package=np.str_("crisp"),
                               values={"time": [1], "losers": [120.0]})
    cal = rivers.pest("riv", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    cal.observe(gaining)
    cal.observe(losing, prefix="crisp_flow")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pst = cal.build("riv.pst")
    template = Path(cal.template_workspace)
    obs = pst.observation_data
    rock = obs.loc[obs["weight"] > 0]
    assert sorted(set(rock["obgnme"])) == ["baseflow", "crisp_flow"]
    assert sorted(rock["weight"].round(6).tolist()) == [0.01, 0.01, 1.0]
    forward = (template / "forward_run.py").read_text()
    assert forward.count("write_named_series_targets(sets=") == 1
    assert "'kind': 'riv_flow'" in forward and "'package': 'crisp'" in forward
    assert "np.str_" not in forward

    tables = ["riv_flow_simulated_riv_flow.csv", "crisp_flow_simulated_riv_flow.csv"]
    bases = {name: (template / name).read_text() for name in tables}
    np.testing.assert_allclose(pd.read_csv(template / tables[0])["rock"],
                               _truth(rivers, [(0, c) for c in GAIN], ("RIV",)), rtol=1e-10)
    write_named_series_targets(
        [{"kind": "riv_flow", "locations_file": template / "riv_flow_target_locations.csv",
          "output_csv": tmp_path / "a.csv", "package": "riv"},
         {"kind": "riv_flow", "locations_file": template / "crisp_flow_target_locations.csv",
          "output_csv": tmp_path / "b.csv", "package": "crisp"}],
        sim_ws=rivers.workspace, model_name="m")
    assert (tmp_path / "a.csv").read_text() == bases[tables[0]]
    assert (tmp_path / "b.csv").read_text() == bases[tables[1]]
    for name in tables:
        (template / name).unlink()
    subprocess.run([sys.executable, "forward_run.py"], cwd=template, check=True,
                   capture_output=True, text=True)
    for name in tables:
        assert (template / name).read_text() == bases[name], name


@pytest.mark.slow
def test_pestpp_ies_runs_a_riv_flow_calibration_end_to_end(rivers, tmp_path):
    """Real PESTPP-IES: two iterations of three realizations reading two RIV sets,
    each realization's riv_flow observations read back from the forward run. The
    ``losers`` zone is DRAWN with no layer, so its snapshot is one row per layer
    (ledger 535); the review reads it back as ONE zone -- it listed the
    zone's residual once per layer in ``obs_residuals()``. The ``rock`` zone is a cell
    LIST with no layer, snapshotted one row per layer the same way (ledger
    550) and read back as one zone too."""

    targets = mf.RivFlowTargets(locations={"rock": GAIN}, package="riv",
                                values={"time": [0, 1], "rock": [-900.0, -950.0]})
    cal = rivers.pest("ies", workspace=tmp_path / "ies", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    cal.observe(targets)
    polys = rivers.vor.gdf_vorPolys.geometry
    drawn = gpd.GeoDataFrame({"name": ["losers"]}, crs=rivers.vor.crs, geometry=[
        polys.loc[LOSE[0]].union(polys.loc[LOSE[1]]).buffer(-1.0)])
    cal.observe(mf.RivFlowTargets(locations=drawn, package="crisp",
                                  values={"time": [0], "losers": [150.0]}), prefix="crisp_flow")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build("ies.pst")
    result = cal.run_ies(reals=3, iterations=1, ies_lambda_mults=1.0, lambda_scale_fac=1.0)
    assert result.iterations == [0, 1]
    posterior = result.posterior._df
    columns = [c for c in posterior.columns if c.startswith("oname:riv_flow")]
    assert columns and np.isfinite(posterior[columns].to_numpy(dtype=float)).all()
    assert (posterior[columns].to_numpy(dtype=float) < 0).all()  # gaining: MF6's sign
    residuals = result.obs_residuals()
    assert sorted(zip(residuals["prefix"], residuals["location"], strict=True)) == [
        ("crisp_flow", "losers"), ("riv_flow", "rock")]
    losers = residuals.set_index("location").loc["losers"]
    assert losers["cells"] == LOSE and losers["n"] == 1
    snapshot = pd.read_csv(Path(cal.template_workspace) / "riv_flow_target_locations.csv")
    assert snapshot["layer"].tolist() == [0, 1]
    assert residuals.set_index("location").loc["rock"]["cells"] == GAIN
