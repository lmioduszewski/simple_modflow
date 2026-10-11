"""A child process started in another directory imports ANOTHER myflopy under a relative PYTHONPATH.

Ledger 354 (ledger 344 found it). Run it from a worktree whose venv has the MAIN
checkout installed editable, with a relative path:

    cd <worktree> && PYTHONPATH=src /home/lukem/python/envs/gw/bin/python \\
        docs/reproducers/relative_pythonpath_subprocess.py

THE SETUP. Python resolves a relative ``PYTHONPATH`` entry against the directory a
process STARTS in. This process starts in the worktree, so ``src`` is the worktree's
``src``. It then starts two children in an empty temporary directory -- as PEST++ starts
``python forward_run.py`` in its template -- one inheriting ``PYTHONPATH`` as given and
one with every relative entry made absolute against this process's start directory (what
``tests/conftest.py`` does at session start since ledger 354), and prints which
``myflopy`` each imports.

WHAT IT SHOWS. Under ``PYTHONPATH=src``, the inheriting child imports the venv's
editable install (the main checkout), not the code this process imports; the rewritten
child imports the same file. Under an absolute ``PYTHONPATH``, or none, all three agree.
"""

from __future__ import annotations

import os
import pathlib
import subprocess
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
# The rewrite tests/conftest.py applies, from the one module both use (ledger 425).
from _absolute_pythonpath import absolute_pythonpath  # noqa: E402


def _child_imports(cwd: str, pythonpath: str | None) -> str:
    env = dict(os.environ)
    if pythonpath is None:
        env.pop("PYTHONPATH", None)
    else:
        env["PYTHONPATH"] = pythonpath
    result = subprocess.run(
        [sys.executable, "-c", "import myflopy; print(myflopy.__file__)"],
        cwd=cwd, env=env, capture_output=True, text=True, check=True,
    )
    return result.stdout.strip().splitlines()[-1]


def main() -> None:
    import myflopy

    given = os.environ.get("PYTHONPATH")
    start = os.getcwd()
    print(f"PYTHONPATH as given          : {given!r}")
    print(f"this process imports         : {myflopy.__file__}")
    with tempfile.TemporaryDirectory() as elsewhere:
        print(f"child in {elsewhere}, inherited: {_child_imports(elsewhere, given)}")
        rewritten = absolute_pythonpath(given, start) if given else None
        print(f"child there, made absolute   : {_child_imports(elsewhere, rewritten)}"
              f"  (PYTHONPATH {rewritten!r})")


if __name__ == "__main__":
    main()
