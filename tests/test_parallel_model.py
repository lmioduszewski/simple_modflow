from __future__ import annotations

import shutil
from pathlib import Path
from types import SimpleNamespace

import flopy
import numpy as np
import pytest

from myflopy.modflow.mf6.canonical import canonical_partition_mask
from myflopy.modflow.mf6.parallel import (
    ParallelCompatibilityError,
    ParallelModelWorkflow,
    ParallelSplitRun,
    _mover_column_groups,
)


class _FakeModel:
    def __init__(self, name, model_type):
        self.name = name
        self.model_type = model_type


class _FakeSimulation:
    def __init__(self, models):
        self._models = {model.name: model for model in models}
        self.model_names = tuple(self._models)
        self.path = None

    def get_model(self, name):
        return self._models.get(name)

    def set_sim_path(self, path):
        self.path = Path(path)


class _FakeSplitter:
    calls = []

    def __init__(self, sim, modelname=None):
        self.sim = sim
        self.modelname = modelname
        self._modelname = modelname or sim.model_names[0]
        self._fdigits = 1

    def split_model(self, mask):
        self.calls.append(("split_model", np.asarray(mask).copy()))
        return _FakeSimulation([_FakeModel("gwf_0", "gwf6"), _FakeModel("gwf_1", "gwf6")])

    def split_multi_model(self, mask):
        self.calls.append(("split_multi_model", np.asarray(mask).copy()))
        return _FakeSimulation(
            [
                _FakeModel("gwf_0", "gwf6"),
                _FakeModel("gwf_1", "gwf6"),
                _FakeModel("gwt_0", "gwt6"),
                _FakeModel("gwt_1", "gwt6"),
            ]
        )


def _source(models):
    sim = _FakeSimulation(models)
    return SimpleNamespace(sim=sim, gwf=models[0])


def test_unified_split_dispatches_single_and_coupled_simulations(monkeypatch, tmp_path):
    import flopy.mf6.utils

    monkeypatch.setattr(flopy.mf6.utils, "Mf6Splitter", _FakeSplitter)
    _FakeSplitter.calls.clear()
    mask = np.array([0, 1])

    single = ParallelModelWorkflow(_source([_FakeModel("gwf", "gwf6")]))
    single_run = single.split_model(workspace=tmp_path / "single", mask=mask, write=False)
    assert isinstance(single_run, ParallelSplitRun)
    assert single.topology()["operation"] == "split_model"
    assert _FakeSplitter.calls[-1][0] == "split_model"

    coupled = ParallelModelWorkflow(
        _source([_FakeModel("gwf", "gwf6"), _FakeModel("gwt", "gwt6")])
    )
    coupled_run = coupled.split_model(workspace=tmp_path / "coupled", mask=mask, write=False)
    assert coupled_run.nparts == 2
    assert coupled.topology()["operation"] == "split_multi_model"
    assert _FakeSplitter.calls[-1][0] == "split_multi_model"
    assert isinstance(single.prepare(workspace=tmp_path / "prepared", mask=mask, write=False), ParallelSplitRun)


def test_unified_split_rejects_unsupported_simulation_topologies():
    workflow = ParallelModelWorkflow(
        _source([_FakeModel("gwf_a", "gwf6"), _FakeModel("gwf_b", "gwf6")])
    )
    with pytest.raises(ParallelCompatibilityError, match="exactly one foundational GWF"):
        workflow.topology()

    workflow = ParallelModelWorkflow(
        _source([_FakeModel("gwf", "gwf6"), _FakeModel("prt", "prt6")])
    )
    with pytest.raises(ParallelCompatibilityError, match="prt"):
        workflow.topology()


def test_automatic_partitioning_has_actionable_optional_dependency_error(monkeypatch, tmp_path):
    import builtins

    workflow = ParallelModelWorkflow(_source([_FakeModel("gwf", "gwf6")]))
    original_import = builtins.__import__

    def fail_pymetis(name, *args, **kwargs):
        if name == "pymetis":
            raise ImportError("missing")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fail_pymetis)
    with pytest.raises(ImportError, match="Automatic partitioning requires pymetis"):
        workflow.split_model(workspace=tmp_path / "auto", nparts=2, write=False)


def test_parallel_run_summary_environment_and_model_names(tmp_path):
    source = SimpleNamespace(
        gwf=SimpleNamespace(
            name="gwf",
            modelgrid=SimpleNamespace(idomain=np.array([[1, 1], [1, 0]])),
        )
    )
    splitter = SimpleNamespace(_fdigits=2)
    run = ParallelSplitRun(
        source,
        splitter,
        _FakeSimulation([]),
        np.array([0, 1]),
        tmp_path / "split",
    )

    summary = run.summary()
    assert summary["active_cells"].tolist() == [2, 1]
    assert summary["active_fraction"].sum() == pytest.approx(1.0)
    assert run.partition_model_name("gwf", 1) == "gwf_01"
    assert isinstance(run.environment.parallel_ready, bool)


def test_parallel_run_requires_mpi_only_for_parallel_execution(monkeypatch, tmp_path):
    calls = []

    class _RunnableSimulation(_FakeSimulation):
        def write_simulation(self):
            calls.append(("write",))

        def run_simulation(self, **kwargs):
            calls.append(("run", kwargs))
            return True, []

    source = SimpleNamespace(gwf=SimpleNamespace(name="gwf"))
    splitter = SimpleNamespace(_fdigits=1)
    run = ParallelSplitRun(
        source,
        splitter,
        _RunnableSimulation([]),
        np.array([0, 1]),
        tmp_path / "split",
    )

    assert run.run_serial(write=False)[0] is True
    assert calls[-1][1].get("processors") is None
    monkeypatch.setattr(
        ParallelSplitRun,
        "environment",
        property(lambda _self: SimpleNamespace(parallel_ready=False)),
    )
    with pytest.raises(RuntimeError, match="mpiexec"):
        run.run(processors=2, write=False)


def test_parallel_write_filters_bound_names_from_saved_node_mapping(tmp_path):
    class WritableSimulation(_FakeSimulation):
        def write_simulation(self):
            return None

    class MappingSplitter:
        _fdigits = 1
        _ncpl = 2

        def __init__(self):
            self._node_map = {
                np.int64(0): (0, np.int64(0)),
                np.int64(1): (1, np.int64(0)),
                "unconfined_supply": [(1, "unconfined_supply")],
            }
            self.saved_map = None

        def save_node_mapping(self, path):
            self.saved_map = dict(self._node_map)
            Path(path).write_text("mapping", encoding="utf-8")

    splitter = MappingSplitter()
    original = dict(splitter._node_map)
    run = ParallelSplitRun(
        SimpleNamespace(gwf=SimpleNamespace(name="gwf")),
        splitter,
        WritableSimulation([]),
        np.array([0, 1]),
        tmp_path / "split",
    )

    run.write()

    assert splitter.saved_map == {0: (0, 0), 1: (1, 0)}
    assert splitter._node_map == original
    assert (run.workspace / "node_mapping.hdf5").exists()


def test_workflow_run_uses_same_unified_split_path(monkeypatch, tmp_path):
    workflow = ParallelModelWorkflow(_source([_FakeModel("gwf", "gwf6")]))
    prepared = SimpleNamespace(run_serial=lambda **_kwargs: (True, []))
    captured = {}

    def fake_split_model(**kwargs):
        captured.update(kwargs)
        return prepared

    monkeypatch.setattr(workflow, "split_model", fake_split_model)
    assert workflow.run(workspace=tmp_path / "run", mask=np.array([0, 1])) is prepared
    assert captured["workspace"] == tmp_path / "run"


def test_unified_split_builds_real_flopy_submodels_and_exchange(tmp_path):
    workspace = tmp_path / "source"
    sim = flopy.mf6.MFSimulation(sim_name="split_demo", sim_ws=workspace, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf", save_flows=True)
    ims = flopy.mf6.ModflowIms(sim, complexity="SIMPLE")
    sim.register_ims_package(ims, [gwf.name])
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=1,
        nrow=2,
        ncol=2,
        delr=1.0,
        delc=1.0,
        top=10.0,
        botm=0.0,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[((0, 0, 0), 6.0), ((0, 1, 1), 4.0)],
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord="gwf.hds",
        budget_filerecord="gwf.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    source = SimpleNamespace(sim=sim, gwf=gwf)

    run = ParallelModelWorkflow(source).split_model(
        workspace=tmp_path / "split",
        mask=np.array([[0, 0], [1, 1]]),
        write=False,
    )

    assert run.nparts == 2
    assert set(run.simulation.model_names) == {"gwf_0", "gwf_1"}
    exchanges = run.simulation.name_file.exchanges.array
    assert len(exchanges) == 1
    assert exchanges[0]["exgtype"] == "GWF6-GWF6"

    if shutil.which("mf6") is None:
        pytest.skip("mf6 executable is required for numerical split-model validation")
    sim.write_simulation()
    source_success, _ = sim.run_simulation(silent=True, report=True)
    run.write(save_mapping=False)
    split_success, _ = run.simulation.run_simulation(silent=True, report=True)
    assert source_success and split_success
    assert run.results.heads() == pytest.approx(gwf.output.head().get_data(), abs=1.0e-5)
    assert run.compare_heads()["max_absolute_error"].iloc[0] < 1.0e-5
    assert run.validate()["contiguous"].all()


def test_unified_split_builds_real_coupled_gwf_gwt_partitions(tmp_path):
    sim = flopy.mf6.MFSimulation(sim_name="coupled_split", sim_ws=tmp_path / "source", exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf", save_flows=True)
    gwt = flopy.mf6.ModflowGwt(sim, modelname="gwt", save_flows=True)
    flow_ims = flopy.mf6.ModflowIms(sim, filename="flow.ims", complexity="SIMPLE")
    transport_ims = flopy.mf6.ModflowIms(sim, filename="transport.ims", complexity="SIMPLE")
    sim.register_ims_package(flow_ims, [gwf.name])
    sim.register_ims_package(transport_ims, [gwt.name])
    for model, constructor in (
        (gwf, flopy.mf6.ModflowGwfdis),
        (gwt, flopy.mf6.ModflowGwtdis),
    ):
        constructor(
            model,
            nlay=1,
            nrow=2,
            ncol=2,
            delr=1.0,
            delc=1.0,
            top=10.0,
            botm=0.0,
        )
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.25)
    flopy.mf6.ModflowGwtadv(gwt)
    flopy.mf6.ModflowGwfgwt(sim, exgmnamea=gwf.name, exgmnameb=gwt.name)

    run = ParallelModelWorkflow(SimpleNamespace(sim=sim, gwf=gwf)).split_model(
        workspace=tmp_path / "split",
        mask=np.array([[0, 0], [1, 1]]),
        write=False,
    )

    assert set(run.simulation.model_names) == {"gwf_0", "gwf_1", "gwt_0", "gwt_1"}
    exchange_types = run.simulation.name_file.exchanges.array["exgtype"].tolist()
    assert exchange_types.count("GWF6-GWF6") == 1
    assert exchange_types.count("GWT6-GWT6") == 1
    assert exchange_types.count("GWF6-GWT6") == 2


def test_partition_validation_rejects_disconnected_custom_mask(tmp_path):
    sim = flopy.mf6.MFSimulation(sim_name="bad_mask", sim_ws=tmp_path / "source")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf")
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=1,
        nrow=2,
        ncol=2,
        delr=1.0,
        delc=1.0,
        top=10.0,
        botm=0.0,
    )
    source = SimpleNamespace(gwf=gwf)
    run = ParallelSplitRun(
        source,
        SimpleNamespace(_fdigits=1),
        _FakeSimulation([]),
        np.array([[0, 1], [1, 0]]),
        tmp_path / "split",
    )
    with pytest.raises(ParallelCompatibilityError, match="contiguous"):
        run.validate()


@pytest.mark.canonical
def test_canonical_model_prepares_contiguous_partitions_across_representative_part_counts(
    canonical_model, tmp_path
):
    # Sweep three representative partition counts (min / mid / max of the
    # supported 2..8 range) rather than all seven -- splitting the canonical
    # model is expensive, and these endpoints exercise the boundary assertions
    # (max-fraction at low counts, small-partition floor at high counts).
    part_counts = (2, 5, 8)
    source_observations = [
        package.package_name
        for package in canonical_model.gwf.packagelist
        if package.package_type == "obs"
    ]
    for nparts in part_counts:
        mask = canonical_partition_mask(canonical_model, nparts)
        lak = canonical_model.gwf.get_package("lak")
        connectiondata = lak.connectiondata.get_data()
        for lake_id in np.unique(connectiondata["ifno"]):
            lake_cells = {
                int(cellid[-1])
                for cellid in connectiondata[connectiondata["ifno"] == lake_id]["cellid"]
            }
            assert len(set(mask[list(lake_cells)])) == 1
        for mover_group in _mover_column_groups(canonical_model.gwf):
            assert len(set(mask[mover_group])) == 1
        run = canonical_model.parallel.split_model(
            workspace=tmp_path / f"split_{nparts}",
            mask=mask,
            write=False,
        )
        validation = run.validate()
        assert run.nparts == nparts
        assert validation["contiguous"].all()
        assert validation["columns"].sum() == canonical_model.gwf.modelgrid.ncpl
        assert validation["columns"].min() >= canonical_model.gwf.modelgrid.ncpl * 0.025
        maximum_fraction = max(0.40, 1.10 / nparts)
        assert validation["columns"].max() <= canonical_model.gwf.modelgrid.ncpl * maximum_fraction
        for model_name in run.simulation.model_names:
            split_model = run.simulation.get_model(model_name)
            lak = split_model.get_package("lak")
            if lak is None:
                continue
            connectiondata = lak.connectiondata.get_data()
            packagedata = lak.packagedata.get_data()
            for lake_record in packagedata:
                lake_id = int(lake_record["ifno"])
                local = connectiondata[connectiondata["ifno"] == lake_id]
                assert int(lake_record["nlakeconn"]) == len(local)
                assert local["iconn"].tolist() == list(range(len(local)))
    assert [
        package.package_name
        for package in canonical_model.gwf.packagelist
        if package.package_type == "obs"
    ] == source_observations


@pytest.mark.canonical
def test_split_restores_source_observations_after_outputs_are_loaded(canonical_run, tmp_path):
    canonical_run.hds.map().plot()
    source_packages = list(canonical_run.gwf.packagelist)
    source_observation_children = {
        id(package): list(package.obs._packages)
        for package in source_packages
        if getattr(getattr(package, "obs", None), "_packages", None)
    }

    for nparts in (2, 3):
        run = canonical_run.parallel.split_model(
            workspace=tmp_path / f"split_after_output_load_{nparts}",
            mask=canonical_partition_mask(canonical_run, nparts),
            write=False,
        )

        assert run.nparts == nparts
        assert canonical_run.gwf.packagelist == source_packages
        assert {
            id(package): list(package.obs._packages)
            for package in canonical_run.gwf.packagelist
            if getattr(getattr(package, "obs", None), "_packages", None)
        } == source_observation_children


@pytest.mark.canonical
def test_canonical_eight_part_split_runs_with_local_lake_connections(canonical_run, tmp_path):
    if shutil.which("mf6") is None:
        pytest.skip("mf6 executable is required for split-model execution validation")

    run = canonical_run.parallel.split_model(
        workspace=tmp_path / "split_8_run",
        mask=canonical_partition_mask(canonical_run, 8),
        save_mapping=False,
    )
    success, report = run.run_serial(write=False, silent=True)

    assert success, "\n".join(report[-30:])
    assert not any("NO DATA SPECIFIED FOR LAKE" in line for line in report)
    comparison = run.compare_heads()
    assert run.results.heads().shape == canonical_run.gwf.modelgrid.shape
    assert comparison.loc[0, "max_absolute_error"] < 0.1, comparison


@pytest.mark.canonical
def test_canonical_eight_part_split_runs_with_mpi(canonical_run, tmp_path):
    if not canonical_run.parallel.environment.parallel_ready:
        pytest.skip("MPI-enabled mf6 and mpiexec are required for parallel execution validation")

    run = canonical_run.parallel.split_model(
        workspace=tmp_path / "split_8_mpi",
        mask=canonical_partition_mask(canonical_run, 8),
        save_mapping=False,
    )
    success, report = run.run(processors=run.nparts, write=False, silent=True)

    assert success, "\n".join(report[-30:])
    assert any("PARALLEL mode" in line for line in report)
    comparison = run.compare_heads()
    assert comparison.loc[0, "max_absolute_error"] < 0.1, comparison
