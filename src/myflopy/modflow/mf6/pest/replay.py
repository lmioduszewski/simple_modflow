"""Replaying a built template's forward run outside PEST++: the steps both replays share.

Two callers run a calibration's ``forward_run.py`` themselves, in a copy, at chosen
parameter values: the write-back (:mod:`~myflopy.modflow.mf6.pest.calibrated`, ledger
560), which stops before MODFLOW 6 to read back what the realization's model read,
and :meth:`PestProject.run_template <myflopy.modflow.mf6.pest.project.PestProject.run_template>`
(ledger 605), which runs it whole -- or without the model command, reading a
completed run's outputs -- and scores what the instruction files read. The steps they
share live here, once:

* :func:`complete_values` -- a value for every parameter of the control file, or a
  refusal naming what is missing, unknown, tied or carries a scale or offset;
* :func:`fill_templates` -- every template filled at FULL precision (``repr``), where
  pyEMU's own filler writes a pilot point at 4 significant digits
  (``docs/reproducers/write_back_template_precision.py``);
* :func:`is_model_command`, :func:`replay_script` and :func:`without_model_command` --
  ``forward_run.py``'s ``main()`` cut at, or without, the statement calling
  ``pyemu.os_utils.run``, found by its syntax tree;
* :func:`run_script` -- the ONE subprocess runner: this interpreter, the copy as working
  directory, and ``PYTHONPATH`` with every relative entry made absolute (ledger 454).
"""

from __future__ import annotations

import ast
import os
import re
import subprocess
import sys
from collections.abc import Mapping
from pathlib import Path

import pandas as pd

from myflopy._logging import get_logger
from myflopy._pythonpath import pythonpath_for_children

logger = get_logger(__name__)

#: The model command pyEMU writes into ``forward_run.py``'s ``main()``: a call of this.
MODEL_COMMAND = ("pyemu", "os_utils", "run")


def complete_values(pst, values: pd.Series, *, source: str, remedy: str) -> dict[str, float]:
    """A value for every parameter of the control file, or a refusal.

    A FIXED parameter ``values`` leaves out is at its ``parval1``; an adjustable one it
    leaves out is refused, and so is a name the control file does not hold -- a
    misspelt fixed parameter would otherwise run at its ``parval1`` in silence (ledger
    605). A TIED parameter is refused (its value is its parent's times a ratio
    this does not reconstruct), and so is one whose ``scale`` is not 1 or whose
    ``offset`` is not 0: PEST writes ``value * scale + offset`` into a template, and
    the templates are filled here with the value alone. myflopy never sets either, so
    one is a hand edit of the control file, and filling past it would run what PEST++
    would not.

    Parameters
    ----------
    pst : pyemu.Pst
        The control file -- anything carrying its ``parameter_data`` (``partrans``,
        ``parval1``, and ``scale``/``offset`` when present), indexed by name.
    values : pandas.Series
        Parameter values indexed by LOWER-CASED pyEMU name.
    source : str
        Where ``values`` came from, for the messages: ``"the realization's parameter
        file"``, ``"values="``.
    remedy : str
        Appended to the refusal of a missing adjustable parameter: what to do instead,
        starting with its own punctuation (``", so what it ran cannot be rebuilt."``).

    Returns
    -------
    dict of str to float
        Every parameter's value by lower-cased name, in the control file's order.

    Raises
    ------
    ValueError
        When a parameter is TIED or carries a scale or offset, when ``values`` names a
        parameter the control file lacks, or when an adjustable parameter has no value
        in ``values``.
    """

    data = pst.parameter_data
    tied = sorted(data.index[data["partrans"] == "tied"])
    if tied:
        raise ValueError(
            f"Parameter(s) {tied[:5]} are TIED; a write-back or a template run reads each "
            "parameter's own value and does not reconstruct a tie. Untie them and run again."
        )
    scale = (pd.to_numeric(data["scale"], errors="coerce") if "scale" in data.columns
             else pd.Series(1.0, index=data.index))
    offset = (pd.to_numeric(data["offset"], errors="coerce") if "offset" in data.columns
              else pd.Series(0.0, index=data.index))
    moved = sorted(data.index[((scale != 1.0) | (offset != 0.0)).to_numpy()])
    if moved:
        first = moved[0]
        raise ValueError(
            f"Parameter(s) {moved[:5]}{' ...' if len(moved) > 5 else ''} carry a scale or "
            f"offset -- {first!r}: scale {float(scale[first]):g}, offset "
            f"{float(offset[first]):g} -- so PEST wrote value * scale + offset into "
            "their templates, which a write-back or a template run, filling each template "
            "with the value, does not reproduce. myflopy sets neither; set scale 1 and "
            "offset 0 in the control file and run again."
        )
    names = [str(name).lower() for name in data.index]
    unknown = [str(key) for key in values.index if str(key) not in set(names)]
    if unknown:
        raise ValueError(
            f"{source} names {len(unknown)} parameter(s) the control file does not hold: "
            f"{unknown[:5]}{' ...' if len(unknown) > 5 else ''}. A misspelling, or the "
            "values of another build of the template; the parameters are "
            "cal.pst.parameter_data.parnme (pyEMU's names, compared without regard to case)."
        )
    complete: dict[str, float] = {}
    missing = []
    for name, key in zip(data.index, names, strict=True):
        if key in values.index:
            complete[key] = float(values[key])
        elif data.loc[name, "partrans"] == "fixed":
            complete[key] = float(data.loc[name, "parval1"])
        else:
            missing.append(key)
    if missing:
        raise ValueError(
            f"Adjustable parameter(s) {missing[:5]}{' ...' if len(missing) > 5 else ''} "
            f"are not in {source}{remedy}"
        )
    return complete


def fill_templates(pst, values: Mapping[str, float], source: Path, target: Path) -> None:
    """Fill every template of ``pst`` from ``source`` into ``target`` at full precision.

    The marker span is replaced by ``repr`` of the value -- the shortest text that reads
    back as the same float -- whatever the span's width: the forward run's own readers
    (pyEMU's multiplier files, the pilot-point CSVs) split on separators, not columns.
    pyEMU's ``Pst.write_input_files`` instead fits a pilot-point template's
    15-character span at 4 significant digits, 4.5e-4 relative
    (``docs/reproducers/write_back_template_precision.py``).

    Parameters
    ----------
    pst : pyemu.Pst
        The control file: its ``template_files`` and the ``input_files`` each fills.
    values : mapping of str to float
        Every parameter's value by lower-cased name (:func:`complete_values`).
    source : pathlib.Path
        The directory the template files are read from.
    target : pathlib.Path
        The directory the filled input files are written into (made where missing).

    Raises
    ------
    ValueError
        When a template does not start ``ptf <marker>``, or names a parameter
        ``values`` lacks.
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


def is_model_command(node: ast.AST) -> bool:
    """Whether ``node`` contains a call of ``pyemu.os_utils.run`` (:data:`MODEL_COMMAND`).

    Parameters
    ----------
    node : ast.AST
        A statement of ``forward_run.py``'s ``main()``.

    Returns
    -------
    bool
        ``True`` when any call inside it is ``pyemu.os_utils.run(...)``.
    """

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
        if tuple(reversed(chain)) == MODEL_COMMAND:
            return True
    return False


def _main_of(forward_run: Path, *, what: str) -> tuple[ast.Module, ast.FunctionDef]:
    """``forward_run``'s syntax tree and its ``main()``, refused when it has none."""

    tree = ast.parse(forward_run.read_text(encoding="utf-8"))
    main = next(
        (node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "main"),
        None,
    )
    if main is None:
        raise RuntimeError(f"{forward_run} defines no main(), so it cannot be {what}.")
    return tree, main


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

    tree, main = _main_of(forward_run, what="replayed")
    cut = next((i for i, statement in enumerate(main.body) if is_model_command(statement)), None)
    if cut is None:
        raise RuntimeError(
            f"{forward_run}'s main() never calls pyemu.os_utils.run, so nothing marks "
            "where its inputs are complete; it cannot be replayed."
        )
    kept = main.body[:cut]
    main.body = kept or [ast.Pass()]
    return ast.unparse(tree) + "\n", [ast.unparse(statement) for statement in kept]


def without_model_command(forward_run: Path) -> str:
    """``forward_run.py`` with the model command removed from ``main()``, every other line kept.

    What :meth:`PestProject.run_template
    <myflopy.modflow.mf6.pest.project.PestProject.run_template>` runs with
    ``outputs_from=``: pyEMU's multipliers and the pilot-point steps before the model,
    and the observation readers after it, read a completed run's outputs copied in.

    Parameters
    ----------
    forward_run : pathlib.Path
        The calibration's generated forward run.

    Returns
    -------
    str
        The script's source.

    Raises
    ------
    RuntimeError
        When the file has no ``main()``, or no statement of it runs the model.
    """

    tree, main = _main_of(forward_run, what="run without its model")
    kept = [statement for statement in main.body if not is_model_command(statement)]
    if len(kept) == len(main.body):
        raise RuntimeError(f"{forward_run}'s main() never runs the model; nothing to cut.")
    main.body = kept or [ast.Pass()]
    return ast.unparse(tree) + "\n"


def run_script(directory: Path, script: str, *, what: str) -> None:
    """Run ``script`` in ``directory`` as a subprocess, the way PEST++ runs the forward run.

    This interpreter (``sys.executable``), ``directory`` as the working directory, and
    ``PYTHONPATH`` with every relative entry made absolute
    (:func:`~myflopy._pythonpath.pythonpath_for_children`, ledger 454), so the script
    imports this process's myflopy. In-process, pyEMU's multiprocessing pool broke
    without a ``__main__`` guard (measured, ledger 560).

    Parameters
    ----------
    directory : pathlib.Path
        The working directory, holding ``script``.
    script : str
        The script's file name in ``directory``.
    what : str
        What is run, for the message: ``"The forward run in <dir>"``.

    Raises
    ------
    RuntimeError
        When the script exits with a nonzero code, quoting the end of its error output;
        or when a relative ``PYTHONPATH`` entry can no longer be resolved (ledger 454).
    """

    environment = dict(os.environ)
    pythonpath = pythonpath_for_children()
    if pythonpath is not None:
        environment["PYTHONPATH"] = pythonpath
    done = subprocess.run(
        [sys.executable, script], cwd=directory, env=environment,
        capture_output=True, text=True, check=False,
    )
    if done.returncode != 0:
        raise RuntimeError(f"{what} failed (exit {done.returncode}):\n{done.stderr[-4000:]}")
    logger.debug("%s finished", what)


__all__ = [
    "MODEL_COMMAND",
    "complete_values",
    "fill_templates",
    "is_model_command",
    "replay_script",
    "run_script",
    "without_model_command",
]
