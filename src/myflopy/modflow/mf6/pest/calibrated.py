"""A calibrated realization's inputs, and the library packages they make (ledger 560).

``IesResults.calibrated(realization="base")`` answers one question: what did that
realization's MODFLOW 6 READ? It answers by REPLAYING the realization's own forward
run, never by re-applying the multipliers in Python:

1. the realization's parameter values -- ``<case>.<iteration>.base.par`` for ``base``
   (16 digits), else that iteration's parameter ensemble (pyEMU writes it at 6
   significant digits, the stated limit for any other realization);
2. a temporary copy of what the forward run reads before MF6 starts -- the
   simulation's inputs (:func:`~myflopy.modflow.mf6.pest.project._simulation_input_files`,
   ledger 383), every file a replayed statement names, and the files a CSV or JSON
   manifest it names lists in turn (pyEMU's ``mult2model_info.csv``, myflopy's
   ``pilot_point_sets.json``);
3. every template in the control file filled at FULL precision (``repr`` of the
   value) -- pyEMU's own filler fits a pilot-point template's 15-character span at 4
   significant digits, 4.5e-4 relative
   (``docs/reproducers/write_back_template_precision.py``);
4. ``forward_run.py`` cut, by its syntax tree, before the first statement calling
   ``pyemu.os_utils.run`` (the model command), and run as a subprocess with the
   absolute ``PYTHONPATH`` (ledger 454) -- the way PESTPP-IES runs it, so physical
   clamps, pilot-point interpolation and the river-bottom carry are the forward
   run's own. In-process, pyEMU's multiprocessing pool broke without a ``__main__``
   guard (measured);
5. each dataset the build's metadata records (``parameters``, ledger 524) read back
   through FloPy with ``load_only=``.

:meth:`CalibratedInputs.packages` then substitutes those values into the packages the
RUN built -- its frozen snapshot (ledger 251), its own library copies (ledger 230) --
never the library as it stands now; :meth:`CalibratedInputs.to_library` registers them
under new keys and forks the run's simulation onto them.
"""

from __future__ import annotations

import ast
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from myflopy._logging import get_logger
from myflopy._pythonpath import pythonpath_for_children
from myflopy.grid_identity import SOURCE_SHA, source_fingerprint
from myflopy.modflow.mf6.pest.ies import IesResults
from myflopy.modflow.mf6.pest.project import _simulation_input_files
from myflopy.project.freeze import thaw_simulation
from myflopy.specs import PackageRef, PackageSpec, package_digest, substitute_values
from myflopy.workspace import Run, _read_package_snapshot

logger = get_logger(__name__)

#: The build-metadata version whose parameter records a write-back reads: each names
#: its model, dataset and written fields (ledger 560).
METADATA_VERSION = 2

#: Where a replay's forward run stops: the first statement calling this (the model
#: command pyEMU writes into ``main()``).
_MODEL_COMMAND = ("pyemu", "os_utils", "run")


@dataclass
class _Dataset:
    """One calibrated FloPy dataset: what the realization's MF6 read, and the run's."""

    model: str
    package: str
    variable: str
    family: str
    fields: tuple[str, ...]
    parameters: list[str]
    additive: bool
    value: Any = None
    base: Any = None
    layers: set = field(default_factory=set)
    periods: set = field(default_factory=set)


def _refuse_old_metadata(results: IesResults) -> list[dict]:
    """The build's parameter records, refused when the template predates them."""

    metadata = results._metadata
    records = metadata.get("parameters")
    if not metadata:
        raise ValueError(
            f"{results.workspace} holds no myflopy_pest_metadata.json, so nothing says "
            "which package and dataset each parameter scaled. A write-back needs a "
            "calibration built by PestProject.build(); rebuild with cal.build() and run it."
        )
    if (
        records is None
        or int(metadata.get("version", 0)) < METADATA_VERSION
        or any("variable" not in record or "model" not in record for record in records)
    ):
        raise ValueError(
            f"The calibration in {results.workspace} was built before its metadata "
            "recorded each parameter's model, dataset and written fields (metadata "
            f"version {metadata.get('version')!r}; a write-back reads version "
            f"{METADATA_VERSION}). Rebuild it with cal.build() and run it again."
        )
    if not records:
        raise ValueError(
            f"The calibration in {results.workspace} parameterizes nothing, so there "
            "is nothing to write back."
        )
    return records


def _iteration(results: IesResults, iteration) -> int:
    """``iteration``, default the posterior, refused when it is not this run's."""

    available = results.iterations
    if not available:
        raise FileNotFoundError(f"No iteration of {results.case} was found in {results.workspace}.")
    if iteration is None:
        return results.posterior_iteration
    if isinstance(iteration, bool) or not isinstance(iteration, (int, np.integer)):
        raise TypeError(f"iteration must be a whole number, not {iteration!r}.")
    if int(iteration) not in available:
        raise ValueError(
            f"Iteration {iteration} is not one of this run's {available} (ledger 507: a "
            "file an earlier run left in the directory is not this run's)."
        )
    return int(iteration)


def _parameter_values(results: IesResults, realization, iteration: int) -> tuple[pd.Series, Path]:
    """The realization's value of every parameter it holds, and the file read.

    The label is checked against the iteration's parameter ensemble through the one
    realization lookup (``IesResults._posterior_row``); ``base`` is then read from
    ``<case>.<iteration>.base.par``, which PESTPP-IES writes at 16 digits, where the
    ensemble holds 6.
    """

    case = results.case
    ensemble = results.par_ensemble(iteration)._df
    row = results._posterior_row(
        realization, ensemble=ensemble,
        what=f"realization in iteration {iteration}'s parameter ensemble",
    )
    if str(realization) != "base":
        path = results._ensemble_file(f"{case}.{iteration}.par")
        return pd.to_numeric(row, errors="raise").astype(float), Path(path)
    path = results.workspace / f"{case}.{iteration}.base.par"
    if not path.exists():
        raise FileNotFoundError(
            f"No base parameter file {path}: PESTPP-IES writes one per iteration for the "
            "'base' realization, and 'base' is read from it alone (16 digits; the "
            "ensemble's 6 would not reproduce its run). Any other realization's label "
            "reads the parameter ensemble."
        )
    frame = results.pyemu.pst_utils.read_parfile(str(path))
    values = pd.to_numeric(frame["parval1"], errors="raise").astype(float)
    values.index = values.index.str.lower()
    return values, path


def _complete_values(pst, values: pd.Series) -> dict[str, float]:
    """A value for every parameter of the control file, or a refusal.

    A FIXED parameter a file leaves out is at its ``parval1``; a TIED one is refused
    (its value is its parent's times a ratio this does not reconstruct), and so is an
    adjustable one the file does not hold.
    """

    data = pst.parameter_data
    tied = sorted(data.index[data["partrans"] == "tied"])
    if tied:
        raise ValueError(
            f"Parameter(s) {tied[:5]} are TIED; a write-back reads each parameter's own "
            "value and does not reconstruct a tie. Untie them and run again."
        )
    complete: dict[str, float] = {}
    missing = []
    for name in data.index:
        key = str(name).lower()
        if key in values.index:
            complete[key] = float(values[key])
        elif data.loc[name, "partrans"] == "fixed":
            complete[key] = float(data.loc[name, "parval1"])
        else:
            missing.append(key)
    if missing:
        raise ValueError(
            f"Adjustable parameter(s) {missing[:5]}{' ...' if len(missing) > 5 else ''} "
            "are not in the realization's parameter file, so what it ran cannot be rebuilt."
        )
    return complete


def _fill_templates(pst, values: Mapping[str, float], source: Path, target: Path) -> None:
    """Fill every template of ``pst`` from ``source`` into ``target`` at full precision.

    The marker span is replaced by ``repr`` of the value -- the shortest text that reads
    back as the same float -- whatever the span's width: the forward run's own readers
    (pyEMU's multiplier files, the pilot-point CSVs) split on separators, not columns.
    pyEMU's ``Pst.write_input_files`` instead fits a pilot-point template's
    15-character span at 4 significant digits, 4.5e-4 relative
    (``docs/reproducers/write_back_template_precision.py``).
    """

    for tpl, model_file in zip(pst.template_files, pst.input_files, strict=True):
        text = (source / tpl).read_text(encoding="utf-8")
        header, _, body = text.partition("\n")
        parts = header.split()
        if len(parts) != 2 or parts[0].lower() != "ptf":
            raise ValueError(f"{source / tpl} does not start with 'ptf <marker>'.")
        marker = re.escape(parts[1])
        pattern = re.compile(f"{marker}([^{marker}]*){marker}")

        def _value(match, tpl=tpl):
            name = match.group(1).strip().lower()
            if name not in values:
                raise ValueError(f"{tpl} names parameter {name!r}, which the control file lacks.")
            return repr(float(values[name]))

        destination = target / model_file
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(pattern.sub(_value, body), encoding="utf-8")


def _is_model_command(node: ast.AST) -> bool:
    """Whether ``node`` contains a call of ``pyemu.os_utils.run``."""

    for inner in ast.walk(node):
        if not isinstance(inner, ast.Call):
            continue
        chain = []
        func = inner.func
        while isinstance(func, ast.Attribute):
            chain.append(func.attr)
            func = func.value
        if isinstance(func, ast.Name):
            chain.append(func.id)
        if tuple(reversed(chain)) == _MODEL_COMMAND:
            return True
    return False


def replay_script(forward_run: Path) -> tuple[str, list[str]]:
    """``forward_run.py`` with ``main()`` cut before the model command, and the lines kept.

    Parameters
    ----------
    forward_run : pathlib.Path
        The calibration's generated forward run.

    Returns
    -------
    tuple of (str, list of str)
        The script's source with ``main()``'s statements from the first one calling
        ``pyemu.os_utils.run`` on removed, and the kept statements as source text.

    Raises
    ------
    RuntimeError
        When the file has no ``main()``, or no statement of it calls
        ``pyemu.os_utils.run``: there is then no point at which the model inputs are
        known to be complete.
    """

    tree = ast.parse(forward_run.read_text(encoding="utf-8"))
    main = next(
        (node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "main"),
        None,
    )
    if main is None:
        raise RuntimeError(f"{forward_run} defines no main(), so it cannot be replayed.")
    cut = next((i for i, statement in enumerate(main.body) if _is_model_command(statement)), None)
    if cut is None:
        raise RuntimeError(
            f"{forward_run}'s main() never calls pyemu.os_utils.run, so nothing marks "
            "where its inputs are complete; it cannot be replayed."
        )
    kept = main.body[:cut]
    main.body = kept or [ast.Pass()]
    return ast.unparse(tree) + "\n", [ast.unparse(statement) for statement in kept]


def _named_paths(statements: list[str], workspace: Path) -> set[Path]:
    """Every file or directory of ``workspace`` a replayed statement names as a string."""

    named: set[Path] = set()
    for source in statements:
        for node in ast.walk(ast.parse(source)):
            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                candidate = _inside(workspace, node.value)
                if candidate is not None:
                    named.add(candidate)
    return named


def _inside(workspace: Path, text: str) -> Path | None:
    """``text`` as a path relative to ``workspace`` when it names something there."""

    text = text.strip()
    if not text or "\n" in text or len(text) > 512:
        return None
    path = Path(text)
    if path.is_absolute():
        return None
    try:
        exists = (workspace / path).exists()
    except OSError:
        return None
    return path if exists else None


def _manifest_paths(path: Path, workspace: Path) -> set[Path]:
    """The files a CSV or JSON manifest in ``workspace`` names (one level, no recursion)."""

    named: set[Path] = set()
    if path.suffix.lower() == ".json":
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return named
        stack = [payload]
        while stack:
            item = stack.pop()
            if isinstance(item, dict):
                stack.extend(item.values())
            elif isinstance(item, list):
                stack.extend(item)
            elif isinstance(item, str):
                found = _inside(workspace, item)
                if found is not None:
                    named.add(found)
    elif path.suffix.lower() == ".csv":
        try:
            with open(path, newline="", encoding="utf-8") as handle:
                for row in csv.reader(handle):
                    for cell in row:
                        found = _inside(workspace, cell)
                        if found is not None:
                            named.add(found)
        except (OSError, csv.Error, UnicodeDecodeError):
            return named
    return named


def _copy(source: Path, target: Path, relative: Path) -> None:
    """Copy ``source/relative`` (a file or a directory) to ``target/relative``."""

    origin, destination = source / relative, target / relative
    if destination.exists():
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    if origin.is_dir():
        shutil.copytree(origin, destination)
    else:
        shutil.copy2(origin, destination)


def _replay(results: IesResults, values: Mapping[str, float], target: Path) -> list[str]:
    """Fill the templates and run the cut forward run in ``target``; the lines replayed."""

    source = results.workspace
    forward_run = source / "forward_run.py"
    if not forward_run.exists():
        raise FileNotFoundError(f"No forward_run.py in {source}, so nothing can be replayed.")
    script, kept = replay_script(forward_run)
    wanted = set(_simulation_input_files(source))
    named = _named_paths(kept, source)
    wanted |= named
    for path in list(named):
        if (source / path).is_file():
            wanted |= _manifest_paths(source / path, source)
    for relative in sorted(wanted):
        _copy(source, target, relative)
    _fill_templates(results.pst, values, source, target)
    (target / "replay_inputs.py").write_text(script, encoding="utf-8")
    environment = dict(os.environ)
    pythonpath = pythonpath_for_children()
    if pythonpath is not None:
        environment["PYTHONPATH"] = pythonpath
    done = subprocess.run(
        [sys.executable, "replay_inputs.py"], cwd=target, env=environment,
        capture_output=True, text=True, check=False,
    )
    if done.returncode != 0:
        raise RuntimeError(
            f"Replaying the forward run of {source} failed (exit {done.returncode}):\n"
            f"{done.stderr[-4000:]}"
        )
    return kept


def _datasets(records: list[dict]) -> dict[tuple[str, str, str], _Dataset]:
    """The (model, package, dataset) triples the parameter records scale, merged."""

    datasets: dict[tuple[str, str, str], _Dataset] = {}
    for record in records:
        for item in record.get("files", []):
            key = (str(record["model"]), str(item["package"]), str(record["variable"]))
            dataset = datasets.get(key)
            if dataset is None:
                dataset = datasets[key] = _Dataset(
                    model=key[0], package=key[1], variable=key[2],
                    family=str(record["family"]), fields=(), parameters=[],
                    additive=bool(record.get("additive")),
                )
            dataset.fields = tuple(dict.fromkeys((*dataset.fields, *record.get("fields", []))))
            if record["name"] not in dataset.parameters:
                dataset.parameters.append(str(record["name"]))
            if item.get("layer") is not None:
                dataset.layers.add(int(item["layer"]))
            if item.get("period") is not None:
                dataset.periods.add(int(item["period"]))
    return datasets


def _read(model, dataset: _Dataset):
    """What ``model``'s package holds of the dataset, as FloPy read it."""

    package = model.get_package(dataset.package)
    if package is None or isinstance(package, list):
        raise ValueError(
            f"Model {dataset.model!r} has no single package named {dataset.package!r} to "
            f"read {dataset.variable!r} from."
        )
    data = getattr(package, dataset.variable, None)
    if data is None:
        raise ValueError(f"Package {dataset.package!r} has no dataset {dataset.variable!r}.")
    if dataset.family == "array":
        value = data.get_data()
        if isinstance(value, dict):
            return {int(per): np.array(array, dtype=float) for per, array in value.items()}
        return np.array(data.array, dtype=float)
    value = data.get_data()
    if isinstance(value, dict):
        return {int(per): records.copy() for per, records in value.items() if records is not None}
    return value.copy()


def _load(workspace: Path, datasets) -> dict[str, Any]:
    """The models of the simulation in ``workspace``, loading only the touched packages."""

    import flopy

    packages = sorted({dataset.package for dataset in datasets})
    simulation = flopy.mf6.MFSimulation.load(
        sim_ws=str(workspace), load_only=packages, verbosity_level=0,
    )
    return {name: simulation.get_model(name) for name in simulation.model_names}


def _name_file_packages(model) -> list[str]:
    """The model's packages by NAME, in its name file's order."""

    rows = model.name_file.packages.get_data()
    return [str(row[2] if len(row) > 2 and row[2] else row[0]).lower() for row in rows]


def _records_frame(value, fields: tuple[str, ...]) -> pd.DataFrame:
    """A list dataset as a table, one row per record, with ``per`` when it has periods."""

    frames = []
    items = value.items() if isinstance(value, dict) else [(None, value)]
    for per, records in items:
        frame = pd.DataFrame.from_records(records, columns=records.dtype.names)
        if per is not None:
            frame.insert(0, "per", int(per))
        frames.append(frame)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _identity(value) -> Any:
    """A record field's identity form: cellids as int tuples, numbers as numbers, text
    LOWER-CASED.

    FloPy lower-cases every string it LOADS (a boundname ``ghb_L3`` reads back
    ``ghb_l3``), and MODFLOW 6 reads a boundname or keyword without regard to case, so
    text compared verbatim refused every package whose run declared a capital letter --
    Ten Trails' GHB carries ``ghb_L3`` (ledger 562).
    """

    if isinstance(value, (tuple, list)):
        return tuple(int(item) for item in value)
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (str, np.str_)):
        return str(value).lower()
    return value


def _substituted_records(base, read, fields: tuple[str, ...], label: str):
    """``base``'s records with ``fields`` taken from ``read``, after checking each cell.

    ``base`` is the package option -- ``{period: [record, ...]}`` or ``[record, ...]``;
    ``read`` what FloPy read of it -- ``{period: recarray}`` or a recarray. Every field
    of ``read`` that is not a float (the cellid, a feature number, a boundname) must
    equal the base record's by position -- text without regard to case
    (:func:`_identity`) -- or the template is not this run. The run's own values of
    those fields are kept, its spelling of a boundname included.
    """

    def _one(records, recarray, where):
        records = list(records)
        names = recarray.dtype.names
        if len(records) != len(recarray):
            raise ValueError(
                f"{label} {where}: the run built {len(records)} record(s) and the "
                f"realization's files hold {len(recarray)}; the calibration's template "
                "is not this run's."
            )
        positions = {name: index for index, name in enumerate(names)}
        missing = [name for name in fields if name not in positions]
        if missing:
            raise ValueError(f"{label} {where}: FloPy read no field(s) {missing}.")
        identity = [
            name for name in names
            if name not in fields and recarray.dtype[name].kind not in "fc"
        ]
        out = []
        for index, (record, row) in enumerate(zip(records, recarray, strict=True)):
            record = list(record)
            if len(record) != len(names):
                raise ValueError(
                    f"{label} {where}: record {index} holds {len(record)} value(s) and "
                    f"FloPy reads {len(names)} ({list(names)}); it cannot be matched by "
                    "position."
                )
            for name in identity:
                ours, theirs = _identity(record[positions[name]]), _identity(row[name])
                if ours != theirs:
                    raise ValueError(
                        f"{label} {where}: record {index}'s {name} is {ours!r} in the run "
                        f"and {theirs!r} in the realization's files; the calibration's "
                        "template is not this run's."
                    )
            for name in fields:
                record[positions[name]] = float(row[name])
            out.append(record)
        return out

    if isinstance(read, dict):
        if not isinstance(base, Mapping):
            raise ValueError(f"{label}: the run declared no period blocks to substitute into.")
        periods = sorted(int(per) for per in base)
        if periods != sorted(read):
            raise ValueError(
                f"{label}: the run declared periods {periods} and the realization's files "
                f"hold {sorted(read)}."
            )
        return {
            per: _one(base[per], read[int(per)], f"period {int(per)}")
            for per in base
        }
    return _one(base, read, "records")


class CalibratedInputs:
    """What one realization's MODFLOW 6 read, rebuilt by replaying its forward run.

    Made by :meth:`IesResults.calibrated <myflopy.modflow.mf6.pest.ies.IesResults.calibrated>`
    (ledger 560), which documents how. Every value is held in memory; the
    replay's directory is gone.

    Attributes
    ----------
    realization : str
        The realization's label.
    iteration : int
        The iteration it was read from.
    case : str
        The control file's stem (the calibration's name).
    workspace : pathlib.Path
        The PESTPP-IES run directory read (the master of a parallel run).
    run_workspace : pathlib.Path
        The run the calibration was built from (its metadata's ``original_workspace``).
    parameter_file : pathlib.Path
        The file the parameter values came from.
    replayed : list of str
        The forward-run statements run before the model command.
    """

    def __init__(
        self, *, results: IesResults, realization: str, iteration: int,
        parameter_file: Path, values: dict[str, float],
        datasets: dict[tuple[str, str, str], _Dataset], name_files: dict[str, list[str]],
        replayed: list[str],
    ):
        """Hold a replay's result; built by :meth:`IesResults.calibrated`, not by hand."""

        self.realization = str(realization)
        self.iteration = int(iteration)
        self.case = results.case
        self.workspace = Path(results.workspace)
        metadata = results._metadata
        self.run_workspace = Path(metadata["original_workspace"])
        self.pest_name = str(metadata.get("project_name") or results.case)
        self.parameter_file = Path(parameter_file)
        self.replayed = list(replayed)
        self._values = dict(values)
        self._datasets = datasets
        self._name_files = name_files
        self._packages: dict[str, PackageSpec] | None = None
        self._labels: dict[tuple[str, str], str] | None = None

    def __repr__(self) -> str:
        """One line: the realization, iteration and the datasets rebuilt."""

        touched = ", ".join(f"{d.package}.{d.variable}" for d in self._datasets.values())
        return (
            f"CalibratedInputs({self.case!r}, realization={self.realization!r}, "
            f"iteration={self.iteration}: {touched})"
        )

    # -- reading -------------------------------------------------------------

    def _dataset(self, package: str, variable: str, model: str | None) -> _Dataset:
        """The one dataset ``(model, package, variable)`` names, or a refusal listing them."""

        matches = [
            dataset for dataset in self._datasets.values()
            if dataset.package == str(package).lower() and dataset.variable == variable
            and (model is None or dataset.model == model)
        ]
        if len(matches) == 1:
            return matches[0]
        known = sorted(f"{d.model}:{d.package}.{d.variable}" for d in self._datasets.values())
        if not matches:
            raise KeyError(
                f"No calibrated dataset {package!r}.{variable!r}"
                + (f" on model {model!r}" if model else "")
                + f"; this realization rebuilt {known}."
            )
        raise KeyError(
            f"{package!r}.{variable!r} is calibrated on several models "
            f"({sorted(d.model for d in matches)}); pass model=."
        )

    def get(self, package: str, variable: str, *, model: str | None = None):
        """The values of one calibrated dataset, as the realization's MODFLOW 6 read them.

        Parameters
        ----------
        package : str
            The FloPy package NAME (``"npf"``, ``"riv"``, ``"crisp"``), any case -- the
            names :meth:`summary` lists in its ``package`` column.
        variable : str
            The FloPy dataset: ``"k"``, ``"k33"``, ``"ss"``, ``"recharge"`` (an
            array RCH), ``"stress_period_data"`` (a list package), ``"packagedata"``
            ... -- the ``variable`` column of :meth:`summary`.
        model : str or None, optional
            The model owning the package, needed only when the same package name and
            dataset are calibrated on two models (a GWF and its GWT). Default
            ``None``.

        Returns
        -------
        numpy.ndarray or dict or pandas.DataFrame
            Griddata (``k``): the whole array, ``(nlay, ncpl)`` on a DISV, every layer
            -- the uncalibrated ones as the model file holds them. A per-period array
            (an array RCH's ``recharge``): ``{period: array}``, one per period block.
            A list dataset: one row per record (with ``per`` for period data) and one
            column per field, the calibrated fields as the realization's files hold
            them -- a ``riv.stage`` offset moves ``rbot`` too (ledger 519). Copies:
            editing one changes nothing here.

        Raises
        ------
        KeyError
            When no calibrated dataset matches, listing those that do, or when two
            models match and ``model`` is ``None``.
        """

        dataset = self._dataset(package, variable, model)
        if dataset.family == "array":
            if isinstance(dataset.value, dict):
                return {per: array.copy() for per, array in dataset.value.items()}
            return dataset.value.copy()
        return _records_frame(dataset.value, dataset.fields)

    def summary(self) -> pd.DataFrame:
        """One row per calibrated dataset: what scaled it and how far it moved from the run.

        Returns
        -------
        pandas.DataFrame
            Columns ``model``, ``package``, ``variable``, ``fields`` (the record
            fields written, empty for an array), ``parameters`` (the
            :meth:`~myflopy.modflow.mf6.pest.project.PestProject.parameterize` names
            scaling it), ``layers`` and ``periods`` (those its files hold; empty when
            the files are not per layer or per period), ``change`` (``"ratio"`` --
            calibrated over the run's value, for a multiplier -- or ``"offset"`` --
            calibrated minus the run's, for an additive elevation, ledger 518) and its
            ``min``, ``median`` and ``max`` over every cell or record, NaN where the
            run's value is 0 (a ratio) or the run's files are gone.
        """

        rows = []
        for dataset in self._datasets.values():
            calibrated, base = self._flat(dataset.value, dataset), self._flat(dataset.base, dataset)
            change = "offset" if dataset.additive else "ratio"
            if base is None or calibrated is None or base.size != calibrated.size:
                moved = np.full(1, np.nan)
            elif dataset.additive:
                moved = calibrated - base
            else:
                with np.errstate(divide="ignore", invalid="ignore"):
                    moved = np.where(base != 0, calibrated / np.where(base != 0, base, 1.0), np.nan)
            finite = moved[np.isfinite(moved)]
            rows.append({
                "model": dataset.model,
                "package": dataset.package,
                "variable": dataset.variable,
                "fields": ", ".join(dataset.fields),
                "parameters": ", ".join(dataset.parameters),
                "layers": ", ".join(str(layer) for layer in sorted(dataset.layers)),
                "periods": ", ".join(str(per) for per in sorted(dataset.periods)),
                "change": change,
                "min": float(finite.min()) if finite.size else np.nan,
                "median": float(np.median(finite)) if finite.size else np.nan,
                "max": float(finite.max()) if finite.size else np.nan,
            })
        return pd.DataFrame(rows, columns=[
            "model", "package", "variable", "fields", "parameters", "layers", "periods",
            "change", "min", "median", "max",
        ])

    @staticmethod
    def _flat(value, dataset: _Dataset) -> np.ndarray | None:
        """Every number of a dataset value the calibration's files hold, flattened in one order.

        Only the layers and periods the parameters' files hold: an uncalibrated layer
        of ``k`` would otherwise pull the median to 1.
        """

        if value is None:
            return None
        if isinstance(value, dict):
            periods = sorted(dataset.periods) or sorted(value)
            items = [value[per] for per in periods if per in value]
        else:
            items = [value]
        if dataset.family == "array":
            arrays = []
            for array in items:
                array = np.asarray(array, dtype=float)
                if dataset.layers and array.ndim >= 2 and not isinstance(value, dict):
                    array = array[sorted(dataset.layers)]
                arrays.append(np.ravel(array))
            return np.concatenate(arrays) if arrays else None
        return np.concatenate([
            np.concatenate([np.asarray(records[name], dtype=float) for name in dataset.fields])
            for records in items
        ]) if items else None

    # -- the run's packages --------------------------------------------------

    def _run(self) -> Run:
        """The run the calibration was built from, reopened from its manifest."""

        try:
            return Run.reopen(self.run_workspace, load=False)
        except FileNotFoundError as error:
            raise FileNotFoundError(
                f"The calibration's model directory {self.run_workspace} holds no run "
                "manifest (run.json): it was not run through a Project, so there are no "
                "packages it built to substitute into. get() and summary() still answer; "
                "packages() and to_library() need a run made by project.prepare_run(...) "
                "and calibrated."
            ) from error

    def packages(self) -> dict[str, PackageSpec]:
        """The packages the RUN built, with this realization's values substituted.

        The base is the run's own record -- its frozen snapshot (ledger 251), whose
        lazy packages are already the arrays and records they resolved to, and its own
        copies of library packages (ledger 230) -- never the library as it stands now,
        which may have moved since the run. Through
        :func:`~myflopy.specs.substitute_values`: an array dataset is replaced whole
        by what the realization's MODFLOW 6 read; a list dataset keeps every record of
        the run's package and takes only its calibrated fields from the realization's
        files, after each record's cellid and other non-float fields are checked
        against the run's, so its grid pin and lake rules stay (ledger 300's trap: a
        correct value edit must not read as "records were edited"). Each package's
        ``metadata["calibrated_from"]`` records the calibration, realization,
        iteration, run, source label and the parameters and forward-run statements
        replayed; ``metadata["source_path"]``/``"source_sha"`` name the parameter file
        read, so ``project.runs.drift()``'s ``source`` column follows it (ledger 276).
        Digests read no metadata, so two write-backs of one realization digest alike.

        Returns
        -------
        dict of str to PackageSpec
            Keyed by the package's LABEL in the run -- its library key for an
            ``mf.ref`` package (``"npf/base"``), else its name (``"ghb"``) -- the keys
            :meth:`to_library` forks from. Computed once.

        Raises
        ------
        FileNotFoundError
            When the calibration's model directory holds no run manifest -- a model
            not run through a Project, or a run deleted since.
        ValueError
            When the run has no frozen snapshot (prepare it again), its packages do not
            line up with its manifest or the realization's name file, a dataset is not
            an option of its package, a list record's cell differs from the run's (the
            template is not this run's), or one label is calibrated on two models.
        """

        if self._packages is not None:
            return dict(self._packages)
        run = self._run()
        spec = thaw_simulation(self.run_workspace, run_name=run.name)
        if spec is None:
            raise ValueError(
                f"The calibration's run {run.name!r} ({self.run_workspace}) has no frozen "
                "snapshot, so what it built cannot be read back. Prepare the run again "
                "(a snapshot is written when a run is built) and calibrate it."
            )
        copies = (
            _read_package_snapshot(self.run_workspace, run.package_snapshot, run_name=run.name)
            if run.package_snapshot else {}
        )
        labels = run.package_keys
        packages: dict[str, PackageSpec] = {}
        owners: dict[str, str] = {}
        self._labels = {}
        for dataset in self._datasets.values():
            label, base = self._base_package(spec, copies, labels, dataset)
            if owners.setdefault(label, dataset.model) != dataset.model:
                raise ValueError(
                    f"Label {label!r} is calibrated on models {owners[label]!r} and "
                    f"{dataset.model!r}; packages() keys by label and cannot hold both."
                )
            self._labels[(dataset.model, dataset.package)] = label
            current = packages.get(label, base)
            packages[label] = self._substitute(current, dataset, label, run)
        self._packages = packages
        return dict(packages)

    def _base_package(self, spec, copies, labels, dataset: _Dataset) -> tuple[str, PackageSpec]:
        """The run's label and built package for the FloPy package ``dataset`` names."""

        model = next((m for m in spec.models if m.name == dataset.model), None)
        if model is None:
            raise ValueError(f"The run's snapshot holds no model {dataset.model!r}.")
        declared = list(labels.get(dataset.model, ()))
        entries = list(model.packages)
        if len(declared) != len(entries):
            raise ValueError(
                f"Model {dataset.model!r}: the run's manifest lists {len(declared)} "
                f"package(s) and its snapshot {len(entries)}; they cannot be paired."
            )
        resolved = []
        for label, entry in zip(declared, entries, strict=True):
            if isinstance(entry, PackageRef):
                if label != entry.key:
                    raise ValueError(
                        f"Model {dataset.model!r}: the manifest's {label!r} sits where the "
                        f"snapshot references {entry.key!r}."
                    )
                if entry.key not in copies:
                    raise ValueError(
                        f"The run holds no copy of library package {entry.key!r} "
                        "(prepared before runs kept one); prepare the run again."
                    )
                entry = copies[entry.key]
            if getattr(entry, "enabled", True):
                resolved.append((label, entry))
        names = self._name_files.get(dataset.model, [])
        if len(names) != len(resolved):
            raise ValueError(
                f"Model {dataset.model!r}: the run built {len(resolved)} package(s) and "
                f"the realization's name file lists {len(names)}; they cannot be paired."
            )
        try:
            index = names.index(dataset.package)
        except ValueError:
            raise ValueError(
                f"Model {dataset.model!r}'s name file has no package {dataset.package!r}."
            ) from None
        return resolved[index]

    def _substitute(self, package: PackageSpec, dataset: _Dataset, label: str, run: Run):
        """``package`` with the dataset's calibrated values, stamped with where they came from."""

        option = dataset.variable
        if option not in package.options and option not in package.inputs:
            raise ValueError(
                f"Package {label!r} declares no {option!r} to substitute; the run's "
                f"package and the calibration's {dataset.package!r} do not line up."
            )
        if dataset.family == "array":
            value = (
                {per: array.copy() for per, array in dataset.value.items()}
                if isinstance(dataset.value, dict) else dataset.value.copy()
            )
        else:
            if option not in package.options:
                raise ValueError(
                    f"Package {label!r}'s {option!r} is a recipe, not records; the run's "
                    "snapshot should hold what it resolved to. Prepare the run again."
                )
            value = _substituted_records(
                package.options[option], dataset.value, dataset.fields, label
            )
        previous = dict((package.metadata or {}).get("calibrated_from") or {})
        provenance = {
            "pest": self.pest_name,
            "case": self.case,
            "workspace": str(self.workspace),
            "iteration": self.iteration,
            "realization": self.realization,
            "run": run.name,
            "run_workspace": str(self.run_workspace),
            "source_key": label,
            "source_digest": (run.package_digests or {}).get(label),
            "parameters": sorted({*previous.get("parameters", []), *dataset.parameters}),
            "datasets": sorted({*previous.get("datasets", []), dataset.variable}),
            "replayed": list(self.replayed),
        }
        declared = package.metadata or {}
        if "source_path" in declared and "declared_source_path" not in previous:
            # An eager `.gpkg` package names the file its records came from; the
            # parameter file now decides its values, so drift's `source` follows that,
            # and the declared file is kept here.
            provenance["declared_source_path"] = declared.get("source_path")
            provenance["declared_source_sha"] = declared.get(SOURCE_SHA)
        for key in ("declared_source_path", "declared_source_sha"):
            if key in previous:
                provenance[key] = previous[key]
        return substitute_values(
            package, {option: value}, provenance=provenance,
            metadata={
                "source_path": str(self.parameter_file),
                SOURCE_SHA: source_fingerprint(self.parameter_file),
            },
        )

    # -- the library ---------------------------------------------------------

    def to_library(
        self, project, *, keys: Mapping[str, str] | None = None, simulation=None,
        overwrite: bool = False, note: str | None = None,
    ) -> dict[str, str]:
        """Register the calibrated packages under NEW keys, and fork the run's simulation onto them.

        Never edits a key: earlier runs, and :meth:`project.runs.drift()
        <myflopy.workspace.RunLibrary.drift>`, keep what they recorded. The arrays save
        through the library's pickle path; on a later grid move a calibrated array
        reads ``unverifiable`` and its route is ``with_grid(remap=True)`` (ledger 300),
        while a list package keeps its pin and reads ``stale`` (re-declare it). Every
        refusal is checked before anything is registered, the simulation fork included.

        The forked simulation is the LIBRARY's entry under the run's simulation name,
        with the calibrated packages swapped in -- so every package it keeps comes from
        the library as it stands now. It is refused when that is not what the run ran
        (ledger 562): a model whose package labels differ from the run's
        manifest, an inline package (or a TDIS/IMS) whose digest differs from the run's
        declared one, or a library key it keeps whose digest differs from the one the
        run recorded (what ``project.runs.drift()`` calls ``edited``/``removed``).
        Otherwise the calibrated values would ride beside packages the calibration
        never ran with, and the fork would not reproduce the realization -- in silence.

        Parameters
        ----------
        project : myflopy.Project
            The project whose run the calibration was built from (its ``runs`` hold a
            run in :attr:`run_workspace`); another raises ``ValueError``.
        keys : mapping of str to str, optional
            New library key per source LABEL (the keys of :meth:`packages`). A label
            left out gets ``f"{label}/{pest}"`` -- ``"npf/base/cal"`` for an
            ``mf.ref("npf/base")`` package of calibration ``cal`` -- with
            ``f"/{realization}"`` appended when the realization is not ``"base"``. A
            label the calibration did not touch raises ``KeyError`` listing those it
            did. Default ``None``: every label takes its default.
        simulation : str or bool or None, optional
            The name of a simulation forked from the run's
            (``project.simulations.fork(<run's simulation>, simulation, swap=...)``),
            each calibrated package swapped in as ``mf.ref(<new key>)`` in its model's
            slot. ``None`` (default) names it ``f"{run's simulation}/{pest}"`` (plus
            the realization suffix); ``False`` forks none. The run's simulation must be
            in ``project.simulations``, else ``KeyError`` naming ``simulation=False``.
        overwrite : bool, optional
            ``False`` (default) raises ``KeyError`` when a new key or the simulation
            name is already registered; ``True`` replaces them.
        note : str or None, optional
            Free text recorded in each package's ``metadata["calibrated_from"]``.
            Default ``None``.

        Returns
        -------
        dict of str to str
            ``{source label: new key}`` for every calibrated package.

        Raises
        ------
        ValueError
            When ``project`` does not own the calibration's run; when the simulation
            to fork is no longer what the run ran (above; pass ``simulation=False``
            to register the packages alone); and what :meth:`packages` raises.
        KeyError
            For a ``keys`` label the calibration did not touch, a key or simulation
            already registered without ``overwrite=True``, or a run whose simulation is
            not in the library while ``simulation`` is not ``False``.
        """

        run = self._owning_run(project)
        packages = self.packages()
        keys = dict(keys or {})
        unknown = sorted(set(keys) - set(packages))
        if unknown:
            raise KeyError(
                f"keys= names {unknown}, which this calibration did not touch; it touched "
                f"{sorted(packages)}."
            )
        suffix = f"/{self.pest_name}" + ("" if self.realization == "base" else f"/{self.realization}")
        chosen = {label: str(keys.get(label, f"{label}{suffix}")) for label in packages}
        if len(set(chosen.values())) != len(chosen):
            raise KeyError(f"keys= gives two packages one key: {chosen}.")
        taken = sorted(key for key in chosen.values() if key in project.packages)
        if taken and not overwrite:
            raise KeyError(
                f"The package library already has {taken}. Pass overwrite=True to replace "
                "them, or keys= to pick others."
            )
        source_simulation = run.spec_name
        name = None
        if simulation is not False:
            name = f"{source_simulation}{suffix}" if simulation in (None, True) else str(simulation)
            if source_simulation not in project.simulations:
                raise KeyError(
                    f"The run's simulation {source_simulation!r} is not in the project's "
                    "simulation library, so it cannot be forked; register it, or pass "
                    "simulation=False to register the packages alone."
                )
            if name in project.simulations and not overwrite:
                raise KeyError(
                    f"The simulation library already has {name!r}. Pass overwrite=True to "
                    "replace it, or simulation= to pick another name."
                )
        swap: dict[str, dict[str, str]] = {}
        for (model, _package), label in (self._labels or {}).items():
            swap.setdefault(model, {})[label] = chosen[label]
        if name is not None:
            self._refuse_a_moved_simulation(project, run, source_simulation, swap)
            # Validated before anything is registered: a swap the fork refuses must
            # not leave the packages behind.
            project.simulations.fork(
                source_simulation, name, swap=swap, overwrite=overwrite, register=False
            )
        for label, package in packages.items():
            if note is not None:
                provenance = dict(package.metadata.get("calibrated_from") or {})
                provenance["note"] = note
                package = package.with_metadata(calibrated_from=provenance)
            project.add_package(chosen[label], package, overwrite=overwrite)
        if name is not None:
            project.simulations.fork(source_simulation, name, swap=swap, overwrite=overwrite)
        return chosen

    def _refuse_a_moved_simulation(self, project, run: Run, source: str, swap) -> None:
        """Refuse to fork ``source`` when the library entry is not what ``run`` ran.

        Compared: each model's package labels against the run's manifest; each inline
        package and simulation-level package against the run's DECLARED spec (its
        manifest recipe -- skipped, with a log line, for a run whose recipe could not
        be stored); and each library key the fork keeps against the digest the run
        recorded (ledger 562).
        """

        library = project.simulations[source]
        declared = self._run()  # the manifest, reopened: labels, digests, recipe
        recorded_keys = declared.package_keys
        recorded_digests = declared.package_digests or {}
        declared_spec = declared.spec if declared.spec_recipe is not None else None
        problems: list[str] = []

        def _label(entry) -> str:
            return entry.key if isinstance(entry, PackageRef) else str(entry.name)

        def _inline_digests(packages) -> dict[str, str | None]:
            return {
                _label(entry): package_digest(entry)
                for entry in packages if not isinstance(entry, PackageRef)
            }

        if declared_spec is not None:
            ran = _inline_digests(declared_spec.packages)
            now = _inline_digests(library.packages)
            problems += [
                f"simulation package {label!r} changed"
                for label in sorted(set(ran) | set(now)) if ran.get(label) != now.get(label)
            ]
        else:
            logger.info(
                "run %r stored no recipe, so the inline packages of simulation %r are "
                "not compared with what it ran", run.name, source,
            )
        declared_models = {m.name: m for m in declared_spec.models} if declared_spec else {}
        for model in library.models:
            labels = tuple(_label(entry) for entry in model.packages)
            if labels != tuple(recorded_keys.get(model.name, ())):
                problems.append(
                    f"model {model.name!r} carries {list(labels)} where the run carried "
                    f"{list(recorded_keys.get(model.name, ()))}"
                )
                continue
            replaced = set((swap or {}).get(model.name, {}))
            if model.name in declared_models:
                ran = _inline_digests(declared_models[model.name].packages)
                now = _inline_digests(model.packages)
                problems += [
                    f"{model.name}: inline package {label!r} changed"
                    for label in sorted(set(ran) | set(now))
                    if label not in replaced and ran.get(label) != now.get(label)
                ]
            for entry in model.packages:
                if not isinstance(entry, PackageRef) or entry.key in replaced:
                    continue
                if entry.key not in project.packages:
                    problems.append(f"{model.name}: library key {entry.key!r} was removed")
                    continue
                then, current = recorded_digests.get(entry.key), package_digest(
                    project.packages[entry.key])
                if then is not None and current is not None and then != current:
                    problems.append(f"{model.name}: library key {entry.key!r} was edited")
        if problems:
            raise ValueError(
                f"The simulation library's {source!r} is not what run {run.name!r} ran, "
                f"so a fork of it would pair the calibrated packages with others the "
                f"calibration never ran with: {'; '.join(problems)}. Pass "
                "simulation=False to register the packages alone (project.runs.drift() "
                "lists the library's edits), then fork the simulation you mean with "
                "project.simulations.fork(..., swap=...)."
            )

    def _owning_run(self, project) -> Run:
        """The run of ``project`` the calibration was built from, or a refusal."""

        target = self.run_workspace.resolve()
        # `project.runs` is a scan cached at first use; a run prepared since is found
        # after one rescan.
        for rescan in (False, True):
            runs = project.refresh_runs() if rescan else project.runs
            for run in runs.values():
                if Path(run.workspace).resolve() == target:
                    return run
        raise ValueError(
            f"Project {getattr(project, 'name', project)!r} owns no run in "
            f"{self.run_workspace}, the run this calibration was built from; pass the "
            "project that prepared it."
        )


def calibrated_inputs(results: IesResults, realization="base", *, iteration=None) -> CalibratedInputs:
    """Replay one realization's forward run and read back what its MODFLOW 6 read.

    The engine of :meth:`IesResults.calibrated
    <myflopy.modflow.mf6.pest.ies.IesResults.calibrated>`, which documents the
    refusals.

    Parameters
    ----------
    results : IesResults
        The PESTPP-IES run to read (its master directory for a parallel run).
    realization : str or int, optional
        The realization's label, compared as text. Default ``"base"``, read from
        ``<case>.<iteration>.base.par``; any other label from the iteration's
        parameter ensemble.
    iteration : int or None, optional
        One of ``results.iterations``. Default ``None``: the posterior iteration.

    Returns
    -------
    CalibratedInputs
        The rebuilt datasets, held in memory; the replay's directory is deleted.
    """

    records = _refuse_old_metadata(results)
    iteration = _iteration(results, iteration)
    values, parameter_file = _parameter_values(results, realization, iteration)
    complete = _complete_values(results.pst, values)
    datasets = _datasets(records)
    with tempfile.TemporaryDirectory(prefix="myflopy-calibrated-") as scratch:
        directory = Path(scratch)
        replayed = _replay(results, complete, directory)
        models = _load(directory, datasets.values())
        name_files = {}
        for dataset in datasets.values():
            model = models.get(dataset.model)
            if model is None:
                raise ValueError(f"The replayed simulation has no model {dataset.model!r}.")
            dataset.value = _read(model, dataset)
            name_files[dataset.model] = _name_file_packages(model)
    run_workspace = Path(results._metadata["original_workspace"])
    if (run_workspace / "mfsim.nam").exists():
        base_models = _load(run_workspace, datasets.values())
        for dataset in datasets.values():
            model = base_models.get(dataset.model)
            if model is not None:
                dataset.base = _read(model, dataset)
    logger.info(
        "rebuilt %d dataset(s) of %s realization %r, iteration %d", len(datasets),
        results.case, realization, iteration,
    )
    return CalibratedInputs(
        results=results, realization=str(realization), iteration=iteration,
        parameter_file=parameter_file, values=complete, datasets=datasets,
        name_files=name_files, replayed=replayed,
    )


__all__ = ["CalibratedInputs", "calibrated_inputs", "replay_script"]
