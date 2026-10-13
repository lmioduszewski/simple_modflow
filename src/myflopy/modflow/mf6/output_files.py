"""Where a model's output files are: its package directory and the files its OC declares.

Moved here from ``myflopy.project.run_model`` so the PEST build -- an import layer
below it -- reads the same rule (ledger 505): its head and concentration
writers had hard-coded ``<model>.hds``/``.ucn``. Path logic only; it imports nothing
of myflopy's.
"""

from __future__ import annotations

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
#: build's concentration writer (ledger 505).
_OC_FILEOUT_SUFFIXES = {"budget": "cbc", "head": "hds", "concentration": "ucn"}


def _name_file_oc(workspace: Path, model_name: str) -> tuple[Path | None, Path | None]:
    """The model's name file and the OC file it names, as ``(name_file, oc_file)``.

    ``name_file`` is ``None`` when there is no name file to read; ``oc_file`` is
    ``None`` when there is none, or when the name file names no OC package.
    """

    model_dir = _model_dir(workspace, model_name)
    name_file = model_dir / f"{model_name}.nam"
    if not name_file.is_file():
        return None, None
    for line in name_file.read_text(encoding="utf-8", errors="ignore").splitlines():
        fields = line.split("#", 1)[0].split()
        if len(fields) >= 2 and fields[0].upper().startswith("OC"):
            return name_file, model_dir / fields[1]
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
    """The ``BUDGET``/``HEAD``/``CONCENTRATION`` output path the model's OC file actually declares.

    Falls back to ``<model_name>.cbc`` / ``.hds`` / ``.ucn`` when the OC file cannot
    be read or names no such record. A TARGET reads its head or concentration
    through :func:`_declared_field_output` instead, which refuses the second case.

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
_FILERECORD_KEYWORDS = {"head": "head_filerecord", "concentration": "concentration_filerecord"}


def _declared_field_output(workspace: Path, model_name: str, kind: str, *, caller: str) -> Path:
    """The head or concentration file a TARGET reads: the one the OC declares, and no other.

    As :func:`_declared_output_path`, except that a READABLE OC declaring no
    ``HEAD FILEOUT`` (``CONCENTRATION FILEOUT``) is REFUSED (ledger
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
        The flow (``"head"``) or transport (``"concentration"``) model's name.
    kind : {"head", "concentration"}
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
