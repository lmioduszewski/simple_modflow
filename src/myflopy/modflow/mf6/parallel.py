"""Unified model-splitting and parallel MODFLOW 6 workflow."""

from __future__ import annotations
from myflopy.viz import mpl_axes

import json
import os
import shutil
import sys
import warnings
from contextlib import contextmanager, nullcontext
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def _find_executable(name: str) -> str | None:
    """Find an executable beside active Python first, then on PATH."""

    suffix = ".exe" if os.name == "nt" else ""
    sibling = Path(sys.executable).with_name(f"{name}{suffix}")
    if sibling.is_file():
        return str(sibling)
    return shutil.which(name)


@contextmanager
def _parallel_execution_environment(simulation):
    """Expose venv-local MPI and MF6 executables to FloPy for one run."""

    mf6 = _find_executable("mf6")
    mpiexec = _find_executable("mpiexec")
    if mf6 is None or mpiexec is None:
        raise RuntimeError(
            "Parallel MODFLOW execution requires both 'mf6' and 'mpiexec' on PATH "
            "or beside the active Python interpreter."
        )
    original_path = os.environ.get("PATH", "")
    original_exe = simulation.exe_name
    os.environ["PATH"] = os.pathsep.join((str(Path(mpiexec).parent), original_path))
    simulation.exe_name = mf6
    try:
        yield
    finally:
        simulation.exe_name = original_exe
        os.environ["PATH"] = original_path


@contextmanager
def _numeric_splitter_node_map(splitter):
    """Temporarily expose only reconstructable numeric grid-node mappings.

    FloPy may add bound-name mappings to ``_node_map`` while splitting packages.
    Those entries are useful during splitting but cannot be serialized by
    ``Mf6Splitter.save_node_mapping()``, which accepts integer grid nodes only.
    """

    original = splitter._node_map
    numeric = {}
    for key, value in original.items():
        try:
            numeric[int(key)] = (int(value[0]), int(value[1]))
        except (IndexError, TypeError, ValueError):
            continue
    expected = int(getattr(splitter, "_ncpl", len(numeric)))
    if len(numeric) != expected:
        raise ParallelCompatibilityError(
            f"Cannot save reconstruction mapping: expected {expected} numeric grid nodes, "
            f"found {len(numeric)}."
        )
    splitter._node_map = numeric
    try:
        yield
    finally:
        splitter._node_map = original


@contextmanager
def _source_without_observation_packages(model):
    """Temporarily hide source data artifacts unsupported by FloPy's splitter.

    FloPy registers utility observation packages both as children of their
    parent package and as standalone entries in the model package list. After
    output data have been loaded, transient package storage can also contain
    numeric keys with ``None`` values that FloPy's splitter treats as record
    arrays. Hiding both OBS registration paths and empty transient entries
    avoids remapping runtime artifacts while preserving the source exactly.
    """

    from flopy.mf6.data.mfdatalist import MFTransientList

    package_container = getattr(model, "_package_container", None)
    if package_container is None:
        yield
        return
    package_list = package_container.packagelist
    original_packages = list(package_list)
    observation_collections = []
    transient_storage_snapshots = []
    for package in original_packages:
        collection = getattr(package, "obs", None)
        children = getattr(collection, "_packages", None)
        if children:
            observation_collections.append((collection, list(children)))
        for block in package.blocks.values():
            for dataset in block.datasets.values():
                if not isinstance(dataset, MFTransientList):
                    continue
                storage = getattr(dataset, "_data_storage", None)
                if not storage:
                    continue
                original_storage = list(storage.items())
                current_key = getattr(dataset, "_current_key", None)
                retained = [
                    (key, value)
                    for key, value in original_storage
                    if dataset.get_data(key=key) is not None
                ]
                if len(retained) != len(original_storage):
                    transient_storage_snapshots.append(
                        (dataset, original_storage, current_key)
                    )
                    storage.clear()
                    storage.update(retained)
                if hasattr(dataset, "_current_key"):
                    dataset._current_key = current_key

    package_list[:] = [
        package for package in original_packages if package.package_type.lower() != "obs"
    ]
    for collection, _ in observation_collections:
        collection._packages = []
    try:
        yield
    finally:
        package_list[:] = original_packages
        for collection, children in observation_collections:
            collection._packages = children
        for dataset, original_storage, current_key in transient_storage_snapshots:
            dataset._data_storage.clear()
            dataset._data_storage.update(original_storage)
            if hasattr(dataset, "_current_key"):
                dataset._current_key = current_key


def _repair_split_lak_packages(simulation) -> None:
    """Repair FloPy-split LAK connection numbering and local connection counts."""

    for model_name in simulation.model_names:
        model = simulation.get_model(model_name)
        if not hasattr(model, "get_package"):
            continue
        package = model.get_package("lak")
        if package is None:
            continue
        connectiondata = package.connectiondata.get_data()
        packagedata = package.packagedata.get_data()
        if connectiondata is None or packagedata is None:
            continue

        repaired_connections = connectiondata.copy()
        lake_ids = np.asarray(repaired_connections["ifno"], dtype=int)
        counts: dict[int, int] = {}
        for lake_id in np.unique(lake_ids):
            indices = np.flatnonzero(lake_ids == int(lake_id))
            repaired_connections["iconn"][indices] = np.arange(len(indices), dtype=int)
            counts[int(lake_id)] = int(len(indices))

        repaired_packagedata = packagedata.copy()
        for index, lake_id in enumerate(np.asarray(repaired_packagedata["ifno"], dtype=int)):
            repaired_packagedata["nlakeconn"][index] = counts.get(int(lake_id), 0)

        package.connectiondata.set_data(repaired_connections)
        package.packagedata.set_data(repaired_packagedata)


def _lake_column_groups(model) -> tuple[np.ndarray, ...]:
    """Return the unique model-grid columns occupied by each physical lake."""

    package = model.get_package("lak")
    if package is None:
        return ()
    connectiondata = package.connectiondata.get_data()
    if connectiondata is None or not len(connectiondata):
        return ()

    modelgrid = model.modelgrid
    ncpl = int(modelgrid.ncpl)
    groups = []
    for lake_id in np.unique(np.asarray(connectiondata["ifno"], dtype=int)):
        records = connectiondata[np.asarray(connectiondata["ifno"], dtype=int) == lake_id]
        columns = []
        for cellid in records["cellid"]:
            cellid = tuple(int(value) for value in cellid)
            try:
                node = int(np.asarray(modelgrid.get_node([cellid])).reshape(-1)[0])
                columns.append(node % ncpl)
            except (AttributeError, IndexError, TypeError, ValueError):
                # DISV cellids are (layer, cell); this fallback also keeps the
                # helper usable with light-weight test doubles.
                columns.append(cellid[-1])
        groups.append(np.unique(np.asarray(columns, dtype=int)))
    return tuple(groups)


def _mover_column_groups(model) -> tuple[np.ndarray, ...]:
    """Return connected SFR/LAK feature columns that must remain co-located."""

    mover = model.get_package("mvr")
    if mover is None:
        return ()
    perioddata = mover.perioddata.get_data()
    if not perioddata:
        return ()

    def feature_columns(package_name, feature_id) -> set[int]:
        """The grid columns (cells) an SFR reach or LAK lake feature occupies."""

        if isinstance(package_name, bytes):
            package_name = package_name.decode()
        package = model.get_package(str(package_name))
        if package is None:
            return set()
        package_type = str(package.package_type).lower()
        if package_type == "sfr":
            data = package.packagedata.get_data()
            rows = data[np.asarray(data["ifno"], dtype=int) == int(feature_id)]
        elif package_type == "lak":
            data = package.connectiondata.get_data()
            rows = data[np.asarray(data["ifno"], dtype=int) == int(feature_id)]
        else:
            return set()
        return {int(tuple(cellid)[-1]) for cellid in rows["cellid"]}

    groups = []
    for records in perioddata.values():
        for record in records:
            columns = feature_columns(record["pname1"], record["id1"])
            columns |= feature_columns(record["pname2"], record["id2"])
            if columns:
                groups.append(columns)

    merged = []
    while groups:
        group = groups.pop()
        overlapping = [candidate for candidate in groups if group & candidate]
        while overlapping:
            for candidate in overlapping:
                group |= candidate
                groups.remove(candidate)
            overlapping = [candidate for candidate in groups if group & candidate]
        merged.append(np.asarray(sorted(group), dtype=int))
    return tuple(merged)


def _preserve_lake_partitions(model, mask) -> np.ndarray:
    """Assign physical lakes and connected MVR features to one partition."""

    preserved = np.asarray(mask, dtype=int).copy()
    flat = preserved.reshape(-1)
    if not hasattr(model, "modelgrid"):
        return preserved
    if flat.size != int(model.modelgrid.ncpl):
        raise ParallelCompatibilityError(
            f"Partition mask has {flat.size} columns, but the model grid has "
            f"{model.modelgrid.ncpl}."
        )
    groups = _lake_column_groups(model) + _mover_column_groups(model)
    neighbors = model.modelgrid.neighbors(reset=False, fast=True)
    for columns in groups:
        partitions, counts = np.unique(flat[columns], return_counts=True)
        # Keep the feature group in its majority partition. np.unique sorts
        # values, so ties resolve deterministically to the lower partition id.
        target = int(partitions[np.argmax(counts)])
        flat[columns] = target
        _connect_partition_components(model, preserved, target, neighbors=neighbors)
    return preserved


def _partition_components(model, mask, partition: int, *, neighbors=None) -> list[set[int]]:
    """The connected components (via grid neighbors) of the cells assigned to ``partition``."""

    if neighbors is None:
        neighbors = model.modelgrid.neighbors(reset=False, fast=True)
    remaining = set(int(node) for node in np.flatnonzero(mask == partition))
    components = []
    while remaining:
        visited = set()
        pending = [next(iter(remaining))]
        while pending:
            node = pending.pop()
            if node in visited:
                continue
            visited.add(node)
            pending.extend(
                int(neighbor)
                for neighbor in neighbors.get(node, [])
                if int(neighbor) in remaining and int(neighbor) not in visited
            )
        components.append(visited)
        remaining -= visited
    return components


def _connect_partition_components(model, mask, partition: int, *, neighbors=None) -> None:
    """Route the shortest grid corridor between disconnected partition pieces."""

    if neighbors is None:
        neighbors = model.modelgrid.neighbors(reset=False, fast=True)
    flat = np.asarray(mask).reshape(-1)
    while True:
        components = _partition_components(model, flat, partition, neighbors=neighbors)
        if len(components) <= 1:
            return
        largest = max(components, key=len)
        destinations = set().union(*(component for component in components if component != largest))
        pending = list(largest)
        predecessor = {node: None for node in largest}
        hit = None
        for node in pending:
            if node in destinations:
                hit = node
                break
            for neighbor in neighbors.get(node, []):
                neighbor = int(neighbor)
                if neighbor not in predecessor:
                    predecessor[neighbor] = node
                    pending.append(neighbor)
        if hit is None:
            raise ParallelCompatibilityError(
                f"Cannot connect protected partition {int(partition)}."
            )
        node = hit
        while node is not None:
            flat[node] = partition
            node = predecessor[node]


def _repair_partition_contiguity(model, mask) -> np.ndarray:
    """Move disconnected partition islands into their strongest neighbor."""

    repaired = np.asarray(mask, dtype=int).copy()
    neighbors = model.modelgrid.neighbors(reset=False, fast=True)
    for _ in range(repaired.size):
        changed = False
        for partition in np.unique(repaired):
            components = _partition_components(
                model,
                repaired,
                int(partition),
                neighbors=neighbors,
            )
            if len(components) <= 1:
                continue
            largest = max(components, key=len)
            for component in components:
                if component == largest:
                    continue
                adjacent = [
                    int(repaired[int(neighbor)])
                    for node in component
                    for neighbor in neighbors.get(node, [])
                    if int(neighbor) not in component
                    and int(repaired[int(neighbor)]) != int(partition)
                ]
                if not adjacent:
                    raise ParallelCompatibilityError(
                        f"Cannot reconnect isolated partition {int(partition)}."
                    )
                targets, counts = np.unique(adjacent, return_counts=True)
                repaired[list(component)] = int(targets[np.argmax(counts)])
                changed = True
        if not changed:
            return repaired
    raise ParallelCompatibilityError("Partition contiguity repair did not converge.")


class ParallelCompatibilityError(ValueError):
    """Raised when a simulation's topology can't be split for parallel solving.

    Domain decomposition supports a single foundational GWF model (optionally
    coupled to GWT/GWE). This error is raised when that contract is violated --
    multiple connected GWF models, or unsupported model types (e.g. PRT) -- so the
    workflow fails clearly instead of producing an invalid partition. A
    :class:`ValueError`.
    """


@dataclass(frozen=True)
class ParallelEnvironment:
    """The resolved MODFLOW 6 executables available for serial vs MPI runs.

    Captures whether the toolchain needed to run a (possibly partitioned) model is
    present: the ``mf6`` executable for serial runs and ``mpiexec`` for MPI-parallel
    runs. ``ParallelModelWorkflow`` discovers these from ``PATH`` and checks
    ``serial_ready`` / parallel readiness before launching, so a missing binary is
    reported up front.

    Attributes
    ----------
    mf6
        Path to the MODFLOW 6 executable, or ``None`` if not found.
    mpiexec
        Path to the MPI launcher, or ``None`` if not found.
    """

    mf6: str | None
    mpiexec: str | None

    @property
    def serial_ready(self) -> bool:
        """True if the MODFLOW 6 executable was found (a serial run is possible)."""

        return self.mf6 is not None

    @property
    def parallel_ready(self) -> bool:
        """True if both MODFLOW 6 and an MPI launcher were found (a parallel run is possible)."""

        return self.mf6 is not None and self.mpiexec is not None


class ParallelSplitResults:
    """Stitch a partitioned run's per-domain outputs back onto the whole grid.

    After a split simulation runs, each subdomain writes its own head/budget
    output. This helper, reached as ``split_run.results``, reads those partial
    outputs and reassembles them into arrays indexed by the *original* (unsplit)
    model grid -- so downstream code sees one coherent result regardless of how
    many partitions ran. Use ``.array(...)`` to pull a reconstructed array for a
    model/variable.

    Parameters
    ----------
    run
        The :class:`ParallelSplitRun` whose outputs are being reconstructed.
    """

    def __init__(self, run: "ParallelSplitRun"):
        """Bind the results reader to a partitioned :class:`ParallelSplitRun`."""

        self.run = run

    def array(
        self,
        *,
        model_name: str | None = None,
        output: str = "head",
        kstpkper: tuple[int, int] | None = None,
        totim: float | None = None,
    ) -> np.ndarray:
        """Return one reconstructed output array from all partitions."""

        base_name = model_name or self.run.source.gwf.name
        arrays = {}
        for partition_id in self.run.partition_ids:
            split_name = self.run.partition_model_name(base_name, partition_id)
            model = self.run.simulation.get_model(split_name)
            if model is None:
                raise KeyError(f"Split model {split_name!r} was not found.")
            output_method = getattr(model.output, output, None)
            if output_method is None:
                raise AttributeError(f"Model {split_name!r} has no {output!r} output accessor.")
            reader = output_method()
            kwargs = {}
            if kstpkper is not None:
                kwargs["kstpkper"] = kstpkper
            if totim is not None:
                kwargs["totim"] = totim
            arrays[partition_id] = reader.get_data(**kwargs)
        reconstructed = np.asarray(self.run.splitter.reconstruct_array(arrays))
        target_shape = tuple(int(value) for value in self.run.source.gwf.modelgrid.shape)
        if reconstructed.shape != target_shape and reconstructed.size == int(np.prod(target_shape)):
            reconstructed = reconstructed.reshape(target_shape)
        return reconstructed

    def heads(
        self,
        *,
        model_name: str | None = None,
        kstpkper: tuple[int, int] | None = None,
        totim: float | None = None,
    ) -> np.ndarray:
        """Return reconstructed groundwater heads."""

        return self.array(model_name=model_name, output="head", kstpkper=kstpkper, totim=totim)


class ParallelSplitRun:
    """A prepared, partitioned simulation: write it, run it, reassemble its results.

    The object returned by ``ParallelModelWorkflow.split_model(...)``. It holds the
    partitioned FloPy ``simulation``, the cell-to-subdomain ``mask`` produced by the
    splitter, and the output ``workspace``, and exposes the run lifecycle: write the
    partitioned input, execute it serially or under MPI, and -- via ``.results``
    (a :class:`ParallelSplitResults`) -- reconstruct outputs on the original grid.

    Parameters
    ----------
    source
        The original (unsplit) model/simulation.
    splitter
        The FloPy splitter that produced the partition.
    simulation
        The partitioned FloPy simulation to run.
    mask
        Per-cell subdomain assignment array.
    workspace
        Directory for the partitioned input/output files.
    """

    def __init__(self, source, splitter, simulation, mask, workspace: str | Path):
        """Hold a partitioned simulation, its cell mask, and workspace, and wire the results reader."""

        self.source = source
        self.splitter = splitter
        self.simulation = simulation
        self.mask = np.asarray(mask, dtype=int)
        self.workspace = Path(workspace)
        self.results = ParallelSplitResults(self)

    @property
    def partition_ids(self) -> tuple[int, ...]:
        """The sorted distinct subdomain ids present in the partition mask."""

        return tuple(int(value) for value in sorted(np.unique(self.mask)))

    @property
    def nparts(self) -> int:
        """The number of subdomains the model was split into."""

        return len(self.partition_ids)

    @property
    def environment(self) -> ParallelEnvironment:
        """The discovered MF6/MPI executables as a :class:`ParallelEnvironment`."""

        return ParallelEnvironment(mf6=_find_executable("mf6"), mpiexec=_find_executable("mpiexec"))

    def partition_model_name(self, base_name: str, partition_id: int) -> str:
        """The FloPy split-model name for one partition (``<base>_<zero-padded id>``)."""

        digits = int(getattr(self.splitter, "_fdigits", len(str(max(self.partition_ids)))))
        return f"{base_name}_{partition_id:0{digits}d}"

    def summary(self) -> pd.DataFrame:
        """Summarize partition sizes and active-cell workloads."""

        grid = getattr(self.source.gwf, "modelgrid", None)
        if grid is None:
            flat_mask = self.mask.reshape(-1)
            return pd.DataFrame(
                [
                    {
                        "partition": partition_id,
                        "columns": int((flat_mask == partition_id).sum()),
                        "active_cells": int((flat_mask == partition_id).sum()),
                        "active_fraction": float((flat_mask == partition_id).mean()),
                    }
                    for partition_id in self.partition_ids
                ]
            )
        idomain = getattr(grid, "idomain", None)
        if idomain is None:
            active_per_column = np.ones(self.mask.size, dtype=int)
        else:
            active = np.asarray(idomain).reshape(-1, self.mask.size) > 0
            active_per_column = active.sum(axis=0)
        flat_mask = self.mask.reshape(-1)
        rows = []
        for partition_id in self.partition_ids:
            selected = flat_mask == partition_id
            rows.append(
                {
                    "partition": partition_id,
                    "columns": int(selected.sum()),
                    "active_cells": int(active_per_column[selected].sum()),
                }
            )
        result = pd.DataFrame(rows)
        total = result["active_cells"].sum()
        result["active_fraction"] = result["active_cells"] / total if total else 0.0
        return result

    def validate(self, *, raise_on_error: bool = True) -> pd.DataFrame:
        """Validate partition contiguity and return the workload summary."""

        result = self.summary()
        modelgrid = getattr(self.source.gwf, "modelgrid", None)
        if modelgrid is None:
            result["contiguous"] = True
            return result
        neighbors = modelgrid.neighbors(reset=True, fast=False)
        flat_mask = self.mask.reshape(-1)
        contiguous = {}
        for partition_id in self.partition_ids:
            nodes = {int(node) for node in np.flatnonzero(flat_mask == partition_id)}
            if not nodes:
                contiguous[partition_id] = False
                continue
            visited = set()
            pending = [next(iter(nodes))]
            while pending:
                node = pending.pop()
                if node in visited:
                    continue
                visited.add(node)
                pending.extend(
                    int(neighbor)
                    for neighbor in neighbors.get(node, [])
                    if int(neighbor) in nodes and int(neighbor) not in visited
                )
            contiguous[partition_id] = visited == nodes
        result["contiguous"] = result["partition"].map(contiguous)
        if raise_on_error and not bool(result["contiguous"].all()):
            bad = result.loc[~result["contiguous"], "partition"].tolist()
            raise ParallelCompatibilityError(f"Partitions must be contiguous; invalid partitions: {bad}.")
        return result

    def write(self, *, save_mapping: bool = True) -> Path:
        """Write the split simulation, partition mask, and reconstruction mapping."""

        self.workspace.mkdir(parents=True, exist_ok=True)
        self.simulation.set_sim_path(self.workspace)
        self.simulation.write_simulation()
        np.save(self.workspace / "partition_mask.npy", self.mask)
        metadata = {
            "source_model": self.source.gwf.name,
            "partition_ids": list(self.partition_ids),
            "nparts": self.nparts,
            "node_mapping": None,
        }
        if save_mapping:
            mapping_path = self.workspace / "node_mapping.hdf5"
            try:
                with _numeric_splitter_node_map(self.splitter):
                    self.splitter.save_node_mapping(mapping_path)
                metadata["node_mapping"] = mapping_path.name
            except ImportError:
                warnings.warn(
                    "Saving the reusable splitter node mapping requires h5py. "
                    "Results can still be reconstructed from this in-memory run; "
                    "install h5py to reopen and reconstruct results later.",
                    RuntimeWarning,
                    stacklevel=2,
                )
            except Exception as exc:  # noqa: BLE001 - persistence is best-effort
                # flopy's HDF5 node-mapping writer can fail on some unstructured
                # grids (e.g. a DISV mapping record whose width does not match the
                # HDF5 dataset dtype). The on-disk mapping is only needed to
                # *reopen* this run later; results are still reconstructed from the
                # in-memory splitter this session, so warn and continue instead of
                # failing the whole split.
                warnings.warn(
                    f"Could not save the reusable splitter node mapping "
                    f"({type(exc).__name__}: {exc}). Results can still be "
                    f"reconstructed from this in-memory run; the on-disk mapping "
                    f"for reopening later was skipped.",
                    RuntimeWarning,
                    stacklevel=2,
                )
                try:  # drop any partial/corrupt file h5py may have created
                    mapping_path.unlink(missing_ok=True)
                except OSError:
                    pass
        (self.workspace / "parallel_split.json").write_text(
            json.dumps(metadata, indent=2),
            encoding="utf-8",
        )
        return self.workspace

    def run_serial(self, *, write: bool = True, **run_kwargs):
        """Run the split simulation serially for equivalence validation."""

        if write:
            self.write()
        run_kwargs.setdefault("silent", False)
        run_kwargs.setdefault("report", True)
        return self.simulation.run_simulation(**run_kwargs)

    def run(self, *, processors: int | None = None, write: bool = True, **run_kwargs):
        """Run serially or through MPI when ``processors`` is greater than one."""

        if processors is None or int(processors) <= 1:
            return self.run_serial(write=write, **run_kwargs)
        if not self.environment.parallel_ready:
            raise RuntimeError(
                "Parallel MODFLOW execution requires both 'mf6' and 'mpiexec' on PATH. "
                "The split simulation can still be validated with run_serial()."
            )
        if write:
            self.write()
        run_kwargs.setdefault("silent", False)
        run_kwargs.setdefault("report", True)
        with _parallel_execution_environment(self.simulation):
            return self.simulation.run_simulation(processors=int(processors), **run_kwargs)

    def compare_heads(
        self,
        *,
        kstpkper: tuple[int, int] | None = None,
        totim: float | None = None,
    ) -> pd.DataFrame:
        """Compare reconstructed split heads with the source model output."""

        reader = self.source.gwf.output.head()
        kwargs = {}
        if kstpkper is not None:
            kwargs["kstpkper"] = kstpkper
        if totim is not None:
            kwargs["totim"] = totim
        original = np.asarray(reader.get_data(**kwargs), dtype=float)
        reconstructed = np.asarray(self.results.heads(kstpkper=kstpkper, totim=totim), dtype=float)
        if reconstructed.shape != original.shape:
            if reconstructed.size != original.size:
                raise ParallelCompatibilityError(
                    f"Cannot compare source heads shaped {original.shape} with reconstructed "
                    f"heads shaped {reconstructed.shape}."
                )
            reconstructed = reconstructed.reshape(original.shape)
        difference = reconstructed - original
        finite = difference[np.isfinite(difference) & (np.abs(difference) < 1.0e29)]
        return pd.DataFrame(
            [
                {
                    "count": int(finite.size),
                    "mean_error": float(finite.mean()) if finite.size else np.nan,
                    "mean_absolute_error": float(np.abs(finite).mean()) if finite.size else np.nan,
                    "max_absolute_error": float(np.abs(finite).max()) if finite.size else np.nan,
                }
            ]
        )

    def plot_partitions(self, *, ax=None, layer: int = 0, cmap: str = "tab20"):
        """Plot the partition mask on the original model grid."""

        import matplotlib.pyplot as plt
        from flopy.plot import PlotMapView

        if ax is None:
            _, ax = mpl_axes(figsize=(10, 8))
        view = PlotMapView(model=self.source.gwf, modelgrid=self.source.gwf.modelgrid, layer=layer, ax=ax)
        image = view.plot_array(self.mask, cmap=cmap)
        view.plot_grid(linewidth=0.2, color="0.35")
        ax.figure.colorbar(image, ax=ax, label="Partition")
        ax.set_title(f"{self.nparts} model partitions")
        ax.set_aspect("equal")
        return ax


class ParallelModelWorkflow:
    """The entry point for domain-decomposing a model to run it in parallel.

    Wraps one model (single GWF, or GWF coupled to GWT/GWE) and drives the
    splitting workflow end to end: inspect the simulation :meth:`topology` and pick
    the right FloPy split operation, discover the run :attr:`environment`
    (mf6/mpiexec), then partition the grid into ``nparts`` subdomains via
    :meth:`split_model`, which returns a :class:`ParallelSplitRun` you write, run
    (serial or MPI), and whose results reassemble onto the original grid. Raises
    :class:`ParallelCompatibilityError` for topologies that cannot be split.

    Parameters
    ----------
    model
        The model (or model-like object exposing ``.sim``) to split.
    """

    def __init__(self, model):
        """Bind the parallel workflow to the ``model`` (or ``.sim``-exposing object) to split."""

        self.model = model

    @property
    def environment(self) -> ParallelEnvironment:
        """The discovered MF6/MPI executables as a :class:`ParallelEnvironment`."""

        return ParallelEnvironment(mf6=_find_executable("mf6"), mpiexec=_find_executable("mpiexec"))

    def topology(self) -> dict[str, Any]:
        """Describe the simulation topology and selected FloPy split operation."""

        models = [self.model.sim.get_model(name) for name in self.model.sim.model_names]
        types = [model.model_type[:3].lower() for model in models]
        gwf_count = types.count("gwf")
        unsupported = sorted(set(types) - {"gwf", "gwt", "gwe"})
        if unsupported:
            raise ParallelCompatibilityError(
                f"Model splitting does not currently support model types: {unsupported}."
            )
        if gwf_count != 1:
            raise ParallelCompatibilityError(
                "Model splitting currently requires exactly one foundational GWF model. "
                "Simulations containing multiple connected GWF models are not supported."
            )
        operation = "split_model" if len(models) == 1 else "split_multi_model"
        return {
            "operation": operation,
            "model_names": tuple(model.name for model in models),
            "model_types": tuple(types),
            "coupled": operation == "split_multi_model",
        }

    def split_model(
        self,
        *,
        workspace: str | Path,
        nparts: int | None = None,
        mask=None,
        active_only: bool = True,
        seed: int = 42,
        contiguous: bool = True,
        verbose: bool = False,
        write: bool = True,
        save_mapping: bool = True,
        preserve_source_observations: bool = True,
        preserve_lakes: bool = True,
    ) -> ParallelSplitRun:
        """Split the simulation using one API regardless of model count.

        This method automatically selects FloPy's lower-level ``split_model``
        or ``split_multi_model`` operation based on the simulation topology.
        """

        from flopy.mf6.utils import Mf6Splitter

        topology = self.topology()
        if mask is None:
            if nparts is None or int(nparts) < 2:
                raise ValueError("Provide a partition mask or nparts greater than one.")
            try:
                import pymetis
            except ImportError as error:
                raise ImportError(
                    "Automatic partitioning requires pymetis. Install it with "
                    "'python -m pip install pymetis', or provide mask= explicitly."
                ) from error
            options = pymetis.Options(seed=int(seed), contig=int(bool(contiguous)))
        else:
            options = None

        # FloPy 3.9.5 does not initialize its internal model name when the
        # optional modelname argument is supplied, so start with its default
        # model and switch explicitly only when the foundational GWF is not
        # first in the simulation.
        splitter = Mf6Splitter(self.model.sim)
        if getattr(splitter, "_modelname", None) != self.model.gwf.name:
            splitter.switch_models(modelname=self.model.gwf.name, remap_nodes=True)
        if mask is None:
            mask = splitter.optimize_splitting_mask(
                nparts=int(nparts),
                active_only=active_only,
                options=options,
                verbose=verbose,
            )
        mask = np.asarray(mask, dtype=int)
        if preserve_lakes:
            mask = _preserve_lake_partitions(self.model.gwf, mask)
        operation = getattr(splitter, topology["operation"])
        if preserve_source_observations:
            observation_context = _source_without_observation_packages(self.model.gwf)
        else:
            observation_context = nullcontext()
        with observation_context:
            split_simulation = operation(mask)
        _repair_split_lak_packages(split_simulation)
        run = ParallelSplitRun(self.model, splitter, split_simulation, mask, workspace)
        run.validate()
        if write:
            run.write(save_mapping=save_mapping)
        return run

    def split(self, **kwargs) -> ParallelSplitRun:
        """Short alias for :meth:`split_model`."""

        return self.split_model(**kwargs)

    def prepare(self, **kwargs) -> ParallelSplitRun:
        """Workflow-oriented alias for :meth:`split_model`."""

        return self.split_model(**kwargs)

    def run(
        self,
        *,
        processors: int | None = None,
        validate_serial: bool = True,
        **split_kwargs,
    ) -> ParallelSplitRun:
        """Prepare and run a split simulation through one convenience call.

        When multiple processors are requested, ``validate_serial=True`` first
        runs the split simulation without MPI before launching the parallel run.
        """

        prepared = self.split_model(**split_kwargs)
        if processors is None or int(processors) <= 1:
            prepared.run_serial(write=False)
            return prepared
        if validate_serial:
            prepared.run_serial(write=False)
        prepared.run(processors=int(processors), write=False)
        return prepared
