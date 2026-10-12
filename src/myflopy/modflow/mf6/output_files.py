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


def _declared_output_path(workspace: Path, model_name: str, kind: str) -> Path:
    """The ``BUDGET``/``HEAD``/``CONCENTRATION`` output path the model's OC file actually declares.

    Falls back to ``<model_name>.cbc`` / ``.hds`` / ``.ucn`` when the OC file cannot
    be read or names no such record.

    A caller that ASSUMES the fallback is wrong whenever the modeller named the
    file themselves: an OC saying ``BUDGET FILEOUT tentrails.cbb`` left every
    cell-budget surface -- ``model.bud(pkg)``, ``model.budget.<term>``,
    ``packages.<pkg>.results.q``, HFB flow -- raising ``FileNotFoundError`` for a
    ``.cbc`` that MF6 never wrote (ledger 301). Parsed from the files rather than
    from the loaded OC package because ``_ensure_packages_loaded`` triggers a
    FULL FloPy load, which is the cost this accessor exists to avoid.
    """

    suffix = _OC_FILEOUT_SUFFIXES[kind]
    model_dir = _model_dir(workspace, model_name)
    fallback = workspace / f"{model_name}.{suffix}"

    name_file = model_dir / f"{model_name}.nam"
    oc_file = None
    if name_file.is_file():
        for line in name_file.read_text(encoding="utf-8", errors="ignore").splitlines():
            fields = line.split("#", 1)[0].split()
            if len(fields) >= 2 and fields[0].upper().startswith("OC"):
                oc_file = model_dir / fields[1]
                break
    if oc_file is None or not oc_file.is_file():
        return fallback

    declared = None
    for line in oc_file.read_text(encoding="utf-8", errors="ignore").splitlines():
        fields = line.split("#", 1)[0].split()
        if (
            len(fields) >= 3
            and fields[0].upper() == kind.upper()
            and fields[1].upper() == "FILEOUT"
        ):
            declared = fields[2].strip("'\"")
            break
    if declared is None:
        return fallback

    # MF6 resolves a package's FILEOUT relative to the simulation directory, but
    # a multi-model run keeps each model's files in its own subdirectory. Take
    # whichever of the two actually exists, and the model dir when neither does
    # (so the error names where the run would have written it).
    candidates = [model_dir / declared, workspace / declared]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return candidates[0]
