"""Make ``PYTHONPATH`` absolute, so a child process imports the myflopy this one does.

Python resolves a relative ``PYTHONPATH`` entry against the directory a process STARTS
in. PEST++ starts ``python forward_run.py`` in the template (and each worker), so under
``PYTHONPATH=src`` the forward run found no ``src`` there and imported the venv's
editable install -- the main checkout -- instead of the tree the reproducer measures
(ledger 354 for the test suite; ledger 425 for the reproducers).

The rule is the library's, ``src/myflopy/_pythonpath.py`` -- ``PestProject.run_ies`` and
``prior`` hand PEST++ the same rewrite (ledger 454). It is loaded here BY PATH from
this checkout's ``src`` (it has no dependencies), so nothing is imported before a
reproducer has made the path absolute, and the checkout whose reproducer this is
supplies the rule. Used by ``tests/conftest.py`` at session start and by every
reproducer whose PEST run starts a forward run that imports myflopy::

    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
    from _absolute_pythonpath import make_pythonpath_absolute  # noqa: E402

    make_pythonpath_absolute()
"""

from __future__ import annotations

import importlib.util
import pathlib

_RULE = pathlib.Path(__file__).resolve().parents[2] / "src" / "myflopy" / "_pythonpath.py"
_spec = importlib.util.spec_from_file_location("_myflopy_pythonpath_rule", _RULE)
_module = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_module)

absolute_pythonpath = _module.absolute_pythonpath
make_pythonpath_absolute = _module.make_pythonpath_absolute

__all__ = ["absolute_pythonpath", "make_pythonpath_absolute"]
