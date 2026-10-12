"""Make ``PYTHONPATH`` absolute, so a child process imports the myflopy this one does.

Python resolves a relative ``PYTHONPATH`` entry against the directory a process STARTS
in. PEST++ starts ``python forward_run.py`` in the template (and in each worker), so
under ``PYTHONPATH=src`` the forward run found no ``src`` there and imported the venv's
editable install -- another checkout -- instead of the code the caller ran (ledger 354
for the test suite, ledger 425 for the reproducers, ledger 454 for
``PestProject.run_ies``/``prior``).

This is the ONE rewrite. It has no dependencies -- not even on the rest of myflopy -- so
``tests/conftest.py`` and ``docs/reproducers/_absolute_pythonpath.py`` load it BY PATH,
before anything is imported, and the library imports it as ``myflopy._pythonpath``.
"""

from __future__ import annotations

import os
import sys


def absolute_pythonpath(value: str, start: str) -> str:
    """``value``, a ``PYTHONPATH``, with every relative entry made absolute against ``start``.

    An empty entry means ``start`` itself, as Python reads one; an absolute entry is
    kept as written. The order is kept, so the import search order does not move.

    Parameters
    ----------
    value : str
        A ``PYTHONPATH`` value, entries separated by ``os.pathsep``.
    start : str
        The directory the relative entries are read against -- the one the process
        that will hand ``value`` on STARTED in.

    Returns
    -------
    str
    """

    return os.pathsep.join(
        entry if os.path.isabs(entry) else os.path.abspath(os.path.join(start, entry))
        for entry in value.split(os.pathsep)
    )


def make_pythonpath_absolute(start: str | None = None) -> str | None:
    """Rewrite this process's ``PYTHONPATH`` with :func:`absolute_pythonpath`, in place.

    Call it before anything changes the working directory: a relative entry means
    what it meant where this process started.

    Parameters
    ----------
    start : str or None, optional
        The directory to read relative entries against. ``None`` (the default) is
        the current working directory.

    Returns
    -------
    str or None
        The rewritten value, or ``None`` when ``PYTHONPATH`` is unset or empty (then
        nothing is set).
    """

    value = os.environ.get("PYTHONPATH")
    if not value:
        return None
    os.environ["PYTHONPATH"] = absolute_pythonpath(value, os.getcwd() if start is None else start)
    return os.environ["PYTHONPATH"]


def pythonpath_for_children() -> str | None:
    """This process's ``PYTHONPATH`` as a child started in ANOTHER directory must get it.

    Each relative entry is made absolute against the current working directory
    (:func:`absolute_pythonpath`) -- which is what Python itself did at this process's
    start, so long as the directory has not changed since. Python put each entry, so
    resolved, on ``sys.path``; a rewritten entry that is NOT there means the working
    directory moved (an ``os.chdir`` in a notebook), and the directory the entry named
    at start cannot be told from it, so that is refused rather than guessed. A process
    started with ``-E``/``-I`` read no ``PYTHONPATH`` and is not checked.

    Returns
    -------
    str or None
        The value to hand a child, or ``None`` when ``PYTHONPATH`` is unset or empty.

    Raises
    ------
    RuntimeError
        When a relative entry, read against the current directory, is not on
        ``sys.path``, naming the entry and both directories.
    """

    value = os.environ.get("PYTHONPATH")
    if not value:
        return None
    start = os.getcwd()
    rewritten = absolute_pythonpath(value, start)
    if sys.flags.ignore_environment:
        return rewritten
    on_path = {os.path.normcase(os.path.abspath(entry or os.curdir)) for entry in sys.path}
    for entry, absolute in zip(value.split(os.pathsep), rewritten.split(os.pathsep), strict=True):
        if os.path.isabs(entry) or os.path.normcase(absolute) in on_path:
            continue
        raise RuntimeError(
            f"PYTHONPATH entry {entry!r} is relative, and read against this process's "
            f"working directory ({start}) it names {absolute}, which is not on sys.path: "
            "the working directory has changed since Python started (or PYTHONPATH has), "
            "so the directory the entry named then cannot be told. PEST++ starts the "
            "forward run in the template, where a relative entry names another directory "
            "and the run imports another myflopy. Start Python with an ABSOLUTE "
            "PYTHONPATH, or set os.environ['PYTHONPATH'] to absolute paths."
        )
    return rewritten
