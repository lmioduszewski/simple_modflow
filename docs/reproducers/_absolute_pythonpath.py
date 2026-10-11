"""Make ``PYTHONPATH`` absolute, so a child process imports the myflopy this one does.

Python resolves a relative ``PYTHONPATH`` entry against the directory a process STARTS
in. PEST++ starts ``python forward_run.py`` in the template (and each worker), so under
``PYTHONPATH=src`` the forward run found no ``src`` there and imported the venv's
editable install -- the main checkout -- instead of the tree the reproducer measures
(ledger 354 for the test suite; ledger 425 for the reproducers).

The ONE rewrite, used by ``tests/conftest.py`` at session start and by every
reproducer whose PEST run starts a forward run that imports myflopy::

    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
    from _absolute_pythonpath import make_pythonpath_absolute  # noqa: E402

    make_pythonpath_absolute()

It has no dependencies, so ``tests/conftest.py`` loads it by path.
"""

from __future__ import annotations

import os


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
