"""Where a model's output files are: its package directory and the files its OC declares.

Moved here from ``myflopy.project.run_model`` so the PEST build -- an import layer
below it -- reads the same rule (ledger 505): its head and concentration
writers had hard-coded ``<model>.hds``/``.ucn``. Path logic only; it imports nothing
of myflopy's.
"""

from __future__ import annotations

import json
from pathlib import Path


def _model_dir(workspace: Path, model_name: str) -> Path:
    """The directory holding ``model_name``'s package files.

    A simulation with SEVERAL models gives each one its own subdirectory
    (``model_rel_path``; see ``workspace._materialize_models``), so a coupled
    GWF+GWT/GWE/PRT run keeps nothing but ``mfsim.nam`` at the top. Single-model
    runs stay flat. Callers that glob for ``<model>.<suffix>`` must look here, not
    in the workspace root, or they silently discover no packages at all.
    """

    candidate = workspace / model_name
    return candidate if candidate.is_dir() else workspace


#: What the model writes when its OC names no file for a kind: ``<model>.<suffix>``.
#: ``concentration`` is a GWT model's OC (``CONCENTRATION FILEOUT``), read by the PEST
#: build's concentration writer (ledger 505); ``temperature`` a GWE model's
#: (``TEMPERATURE FILEOUT``), read by ``model.temp`` (ledger 598).
_OC_FILEOUT_SUFFIXES = {"budget": "cbc", "head": "hds", "concentration": "ucn", "temperature": "ucn"}


def _model_name_file(workspace: Path, model_name: str) -> Path:
    """The model's name file: where ``mfsim.nam``'s models block says it is.

    MF6 reads ``<type>  <name file>  <name>`` and resolves the name file against the
    SIMULATION directory, so a FloPy ``model_nam_file=`` or ``model_rel_path=`` -- and a
    coupled myflopy run, whose models each sit in a subdirectory -- puts it where the
    old guess (``<model>.nam`` in :func:`_model_dir`) never looked. The guess stays the
    fallback for a workspace with no ``mfsim.nam`` (or one not naming the model).
    """

    for _kind, name_file, name in _declared_model_files(workspace):
        if name == model_name:
            return Path(workspace) / name_file
    return _model_dir(workspace, model_name) / f"{model_name}.nam"


def _name_file_entry(workspace: Path, name_file: Path, entry: str) -> Path:
    """Where a file a model's name file names IS: against the simulation directory first.

    MF6 opens every package file relative to the SIMULATION directory -- FloPy writes a
    coupled run's entries as ``flow/flow.disv`` for that reason -- so joining the entry
    onto the model's own directory doubled it (``<ws>/flow/flow/flow.disv``) and every
    reader of a coupled run's OC, discretization or grid file looked in the wrong place
    (the review of ledger 596). The name file's own directory is the fallback; with
    neither existing, the simulation directory's path is returned, so a refusal names
    where MF6 would have looked.
    """

    for candidate in (Path(workspace) / entry, Path(name_file).parent / entry):
        if candidate.is_file():
            return candidate
    return Path(workspace) / entry


def _name_file_lines(name_file: Path) -> list[list[str]]:
    """The ``#``-stripped, whitespace-split lines of a name file (``[]`` if unreadable)."""

    try:
        text = Path(name_file).read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return []
    return [line.split("#", 1)[0].split() for line in text.splitlines()]


def _name_file_oc(workspace: Path, model_name: str) -> tuple[Path | None, Path | None]:
    """The model's name file and the OC file it names, as ``(name_file, oc_file)``.

    ``name_file`` is ``None`` when there is no name file to read; ``oc_file`` is
    ``None`` when there is none, or when the name file names no OC package. Both are
    found as MF6 finds them (:func:`_model_name_file`, :func:`_name_file_entry`).
    """

    name_file = _model_name_file(workspace, model_name)
    if not name_file.is_file():
        return None, None
    for fields in _name_file_lines(name_file):
        if len(fields) >= 2 and fields[0].upper().startswith("OC"):
            return name_file, _name_file_entry(workspace, name_file, fields[1])
    return name_file, None


def _oc_declaration(workspace: Path, model_name: str, kind: str) -> tuple[Path | None, str | None]:
    """The model's OC file and the file it names for ``kind``, as ``(oc_file, declared)``.

    ``oc_file`` is ``None`` when the name file or the OC file it names cannot be
    read, or the name file names no OC; ``declared`` is ``None`` when the OC file
    names no ``<KIND> FILEOUT``.
    """

    _name_file, oc_file = _name_file_oc(workspace, model_name)
    if oc_file is None or not oc_file.is_file():
        return None, None
    for line in oc_file.read_text(encoding="utf-8", errors="ignore").splitlines():
        fields = line.split("#", 1)[0].split()
        if (
            len(fields) >= 3
            and fields[0].upper() == kind.upper()
            and fields[1].upper() == "FILEOUT"
        ):
            return oc_file, fields[2].strip("'\"")
    return oc_file, None


def _declared_output_path(workspace: Path, model_name: str, kind: str) -> Path:
    """The ``BUDGET``/``HEAD``/``CONCENTRATION``/``TEMPERATURE`` output path the model's OC declares.

    Falls back to ``<model_name>.cbc`` / ``.hds`` / ``.ucn`` when the OC file cannot
    be read or names no such record. A TARGET, and ``model.hds``/``.conc``/``.temp``
    (ledger 598), read their field through :func:`_declared_field_output` instead,
    which refuses the second case.

    A caller that ASSUMES the fallback is wrong whenever the modeller named the
    file themselves: an OC saying ``BUDGET FILEOUT tentrails.cbb`` left every
    cell-budget surface -- ``model.bud(pkg)``, ``model.budget.<term>``,
    ``packages.<pkg>.results.q``, HFB flow -- raising ``FileNotFoundError`` for a
    ``.cbc`` that MF6 never wrote (ledger 301). Parsed from the files rather than
    from the loaded OC package because ``_ensure_packages_loaded`` triggers a
    FULL FloPy load, which is the cost this accessor exists to avoid.
    """

    suffix = _OC_FILEOUT_SUFFIXES[kind]
    fallback = workspace / f"{model_name}.{suffix}"
    _oc_file, declared = _oc_declaration(workspace, model_name, kind)
    if declared is None:
        return fallback
    return _resolved_fileout(workspace, model_name, declared)


def _resolved_fileout(workspace: Path, model_name: str, declared: str) -> Path:
    """Where a declared ``FILEOUT`` name lands: the model's directory, else the simulation's."""

    model_dir = _model_dir(workspace, model_name)
    # MF6 resolves a package's FILEOUT relative to the simulation directory, but
    # a multi-model run keeps each model's files in its own subdirectory. Take
    # whichever of the two actually exists, and the model dir when neither does
    # (so the error names where the run would have written it).
    candidates = [model_dir / declared, workspace / declared]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return candidates[0]


#: The ``mf.oc`` keyword that declares each field file, named by the refusal below.
_FILERECORD_KEYWORDS = {
    "head": "head_filerecord",
    "concentration": "concentration_filerecord",
    "temperature": "temperature_filerecord",
}


def _declared_field_output(workspace: Path, model_name: str, kind: str, *, caller: str) -> Path:
    """The head, concentration or temperature file the OC declares, and no other.

    What a TARGET reads, and the file ``model.hds``/``.conc``/``.temp`` open (ledger
    598; they opened ``<model>.<suffix>`` whatever the OC said). As
    :func:`_declared_output_path`, except that a READABLE OC declaring no
    ``HEAD FILEOUT`` (``CONCENTRATION``/``TEMPERATURE FILEOUT``) is REFUSED (ledger
    576): MF6 refuses ``SAVE HEAD`` without one, so such a run saved no
    head at all, and the ``<model>.hds`` the fallback opened was another run's,
    left in a reused directory -- read in silence by ``compare()`` and by the
    PEST build. A readable name file naming NO OC package is refused the same
    way (ledger 576): without an OC MF6 saves nothing. An OC that
    cannot be read (no name file, a duck-typed model, an OC file the name file
    names but the directory lacks) still falls back to ``<model>.<suffix>``.

    Parameters
    ----------
    workspace : pathlib.Path
        The run directory.
    model_name : str
        The flow (``"head"``), transport (``"concentration"``) or energy
        (``"temperature"``) model's name.
    kind : {"head", "concentration", "temperature"}
        Which file.
    caller : str
        What is reading it, named in the refusal.

    Returns
    -------
    pathlib.Path
        The declared file (which may not exist; the reader says so).

    Raises
    ------
    ValueError
        When the OC is readable and declares no ``<KIND> FILEOUT``, or the name
        file is readable and names no OC.
    """

    oc_file, declared = _oc_declaration(workspace, model_name, kind)
    word = kind.upper()
    if oc_file is None:
        name_file, named_oc = _name_file_oc(workspace, model_name)
        if name_file is not None and named_oc is None:
            # No OC at all saves nothing either: the same stale file, one step
            # earlier (ledger 576).
            raise ValueError(
                f"{caller}: model {model_name!r}'s name file {name_file.name!r} names no "
                f"output control (OC6), so the run saved no {kind} and a "
                f"{model_name}.{_OC_FILEOUT_SUFFIXES[kind]} in the run directory is not this "
                f"run's (ledger 576). Give the model an OC that saves it "
                f"(mf.oc({_FILERECORD_KEYWORDS[kind]}=\"{model_name}.{_OC_FILEOUT_SUFFIXES[kind]}\", "
                f"saverecord=[(\"{word}\", \"LAST\"), ...])) and run the model again."
            )
        return workspace / f"{model_name}.{_OC_FILEOUT_SUFFIXES[kind]}"
    if declared is None:
        raise ValueError(
            f"{caller}: model {model_name!r}'s output control {oc_file.name!r} declares no "
            f"{word} FILEOUT, so the run saved no {kind} -- MF6 refuses SAVE {word} without "
            f"one -- and a {model_name}.{_OC_FILEOUT_SUFFIXES[kind]} in the run directory is "
            f"not this run's (ledger 576). Declare the file "
            f"(mf.oc({_FILERECORD_KEYWORDS[kind]}=\"{model_name}.{_OC_FILEOUT_SUFFIXES[kind]}\", "
            f"saverecord=[(\"{word}\", \"LAST\"), ...])) and run the model again."
        )
    return _resolved_fileout(workspace, model_name, declared)


def recorded_crs(workspace: Path, model_name: str) -> str | None:
    """The CRS a myflopy run recorded for ``model_name`` in its manifest, if it has one.

    ``run.json`` is ``workspace.RUN_MANIFEST_NAME``; the name is spelled here because this
    module imports nothing of myflopy's. Shared by ``mf.load_mf6_run`` (through
    ``project.run_model._recorded_crs``) and the ``mf.Heads`` recipe, which reads a run's
    CRS at declaration.

    Parameters
    ----------
    workspace : pathlib.Path
        The run directory.
    model_name : str
        The model whose CRS to read.

    Returns
    -------
    str or None
        The recorded CRS, or ``None`` with no manifest, an unreadable one, or none
        recorded for the model.
    """

    try:
        recorded = json.loads((Path(workspace) / "run.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return (recorded.get("crs") or {}).get(model_name)


def _declared_model_files(workspace: Path) -> list[tuple[str, str, str]]:
    """``(type, name file, name)`` for each model in ``mfsim.nam``'s models block, in order.

    ``<type>  <namefile>  [<name>]``: the name is optional in MF6 and defaults to the name
    file's stem; the name file is as written, relative to the simulation directory. An
    absent or unreadable ``mfsim.nam`` declares nothing. Read with undecodable bytes
    replaced, as the other name-file readers here are: one in a comment, which MF6
    ignores, raised ``UnicodeDecodeError`` through every reader of the models block.
    """

    namefile = Path(workspace) / "mfsim.nam"
    try:
        lines = namefile.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return []
    models: list[tuple[str, str, str]] = []
    in_models = False
    for line in lines:
        stripped = line.split("#", 1)[0].strip()
        lowered = stripped.lower()
        if lowered.startswith("begin") and "models" in lowered:
            in_models = True
            continue
        if lowered.startswith("end") and "models" in lowered:
            break
        if not in_models or not stripped:
            continue
        fields = stripped.split()
        if len(fields) >= 2:
            name = fields[2] if len(fields) >= 3 else Path(fields[1]).stem
            models.append((fields[0].lower(), fields[1], name))
    return models


def _declared_models(workspace: Path) -> list[tuple[str, str]]:
    """``(type, name)`` for each model in ``mfsim.nam``'s models block, in its order.

    :func:`_declared_model_files` without the name file. An absent or unreadable
    ``mfsim.nam`` declares nothing.
    """

    return [(kind, name) for kind, _name_file, name in _declared_model_files(workspace)]


def discretization_file(workspace: Path, model_name: str) -> tuple[str | None, Path | None]:
    """The model's discretization ``(type, file)`` from its name file: ``("disv6", path)``.

    ``type`` is the name file's ``DIS6``/``DISV6``/``DISU6`` entry lower-cased; ``file`` is
    where MF6 opens it (:func:`_name_file_entry`: the simulation directory first).
    ``(None, None)`` when the name file cannot be read or names none.

    Parameters
    ----------
    workspace : pathlib.Path
        The run (simulation) directory.
    model_name : str
        The model.

    Returns
    -------
    tuple of (str or None, pathlib.Path or None)
        The discretization's type and file.
    """

    name_file = _model_name_file(workspace, model_name)
    for fields in _name_file_lines(name_file):
        if len(fields) >= 2 and fields[0].upper() in ("DIS6", "DISV6", "DISU6"):
            return fields[0].lower(), _name_file_entry(workspace, name_file, fields[1])
    return None, None


def binary_grid_file(workspace: Path, model_name: str) -> tuple[Path | None, str]:
    """The binary grid file MF6 wrote for ``model_name``, and how it was named.

    MF6 names it in the discretization's OPTIONS: ``GRB6 FILEOUT <name>`` (FloPy's
    ``grb_filerecord=``) resolved against the simulation directory, none at all under
    ``NOGRB``, and otherwise the discretization file's own name with ``.grb`` appended.
    The review of ledger 596 measured a ``GRB6 FILEOUT mesh.grb`` run refused for a
    missing ``flow.dis.grb``, with advice about ``NOGRB`` it never used.

    Parameters
    ----------
    workspace : pathlib.Path
        The run (simulation) directory.
    model_name : str
        The model.

    Returns
    -------
    tuple of (pathlib.Path or None, str)
        The grid file (which may not exist; the caller says so) and ``"default"``,
        ``"fileout"`` or ``"nogrb"``; ``(None, "nodis")`` when the model's name file names
        no discretization.
    """

    _kind, dis_file = discretization_file(workspace, model_name)
    if dis_file is None:
        return None, "nodis"
    try:
        handle = Path(dis_file).open("r", encoding="utf-8", errors="ignore")
    except OSError:
        return Path(str(dis_file) + ".grb"), "default"
    in_options = False
    with handle:
        for line in handle:
            fields = line.split("#", 1)[0].split()
            if not fields:
                continue
            head = fields[0].upper()
            if head == "BEGIN" and len(fields) > 1 and fields[1].upper() == "OPTIONS":
                in_options = True
                continue
            if head == "END" and in_options:
                break
            if head == "BEGIN":
                break                       # options come first; past them, stop reading
            if not in_options:
                continue
            if head == "NOGRB":
                return None, "nogrb"
            if head == "GRB6" and len(fields) >= 3 and fields[1].upper() == "FILEOUT":
                name = fields[2].strip("'\"")
                for candidate in (Path(workspace) / name, Path(dis_file).parent / name):
                    if candidate.is_file():
                        return candidate, "fileout"
                return Path(workspace) / name, "fileout"
    return Path(str(dis_file) + ".grb"), "default"
