"""Every spliced docstring must also be in its SOURCE, where an editor can see it.

Reported as "model.plot.section docstring doesn't show in PyCharm... it shows the
args but the docstring itself is only one line. This is a reoccurring issue."

Measured: `ModelPlots.section` carried ONE line of source docstring and 169 lines
of runtime `__doc__`, assembled at import by `_inherit_verb_docs`. PyCharm and
Pylance are STATIC -- they read the `def` and the literal under it and never run
the module -- so the reader got the one line. It is exactly the defect plan 8.8
fixed for signatures, one field along: `__doc__` set at import reaches `help()`
and reaches no editor.

Ledger 169b had recorded "signature-static, docstring-runtime" as the contract.
That was never a decision anyone made; it was a description of what the code
happened to do, written up as though it were settled. Reversed here.

`scripts/derive_docstrings.py` writes the spliced text back into each source
file. This is the ratchet on it, in both directions: edit `plot.map`'s reference
without regenerating and it fails; hand-edit a generated docstring and it fails.
"""

from __future__ import annotations

import ast
import inspect
import pathlib
import subprocess
import sys

import pytest

REPO = pathlib.Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "derive_docstrings.py"


def _targets():
    """The spliced methods, discovered exactly as the generator discovers them."""

    sys.path.insert(0, str(REPO / "scripts"))
    import derive_docstrings

    return derive_docstrings.targets()


TARGETS = _targets()
IDS = [f"{cls}.{verb}" for _module, cls, verb in TARGETS]


def test_the_sweep_finds_the_spliced_methods_at_all():
    """A guard on the guard: a discovery test that discovers nothing passes
    everything."""

    assert len(TARGETS) >= 25, f"only {len(TARGETS)} spliced methods found"
    assert ("myflopy.plot", "ModelPlots", "section") in TARGETS


def test_the_generated_docstrings_are_committed():
    """`--check` is the whole contract, run the way CI would run it."""

    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--check"],
        capture_output=True, text=True, cwd=REPO, check=False,
    )
    assert result.returncode == 0, (
        "source docstrings have drifted from the splice:\n"
        f"{result.stdout}{result.stderr}"
    )


@pytest.mark.parametrize(("module_name", "cls_name", "verb"), TARGETS, ids=IDS)
def test_a_static_reader_sees_the_whole_docstring(module_name, cls_name, verb):
    """The literal in the file must equal the assembled `__doc__`.

    Asserted per method rather than only through the script so a failure names
    the one that drifted, and compares the SOURCE literal parsed out with `ast`
    -- which is as close as a test gets to what an editor does.
    """

    import importlib

    cls = getattr(importlib.import_module(module_name), cls_name)
    tree = ast.parse(pathlib.Path(inspect.getsourcefile(cls)).read_text())
    literal = None
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == cls_name:
            for fn in node.body:
                if isinstance(fn, ast.FunctionDef) and fn.name == verb:
                    literal = ast.get_docstring(fn)
    assert literal is not None, f"{cls_name}.{verb} has no source docstring"

    runtime = inspect.cleandoc(getattr(cls, verb).__doc__ or "")
    assert literal.strip() == runtime.strip(), (
        f"{cls_name}.{verb}: a static reader sees {len(literal.splitlines())} "
        f"lines, `help()` sees {len(runtime.splitlines())}. "
        f"Run `python scripts/derive_docstrings.py`."
    )


@pytest.mark.parametrize(("module_name", "cls_name", "verb"), TARGETS, ids=IDS)
def test_a_static_docstring_is_worth_reading(module_name, cls_name, verb):
    """Not just present -- substantial, and documenting what it accepts.

    The failure mode this guards is subtler than an empty docstring: a splice
    that silently produced only the method's own one-line summary would satisfy
    the equality test above perfectly, because runtime would be that one line too.
    """

    import importlib

    method = getattr(getattr(importlib.import_module(module_name), cls_name), verb)
    doc = inspect.cleandoc(method.__doc__ or "")
    assert len(doc.splitlines()) > 10, (
        f"{cls_name}.{verb} has a {len(doc.splitlines())}-line docstring"
    )
    assert "Parameters" in doc, f"{cls_name}.{verb} documents no parameters"


def test_the_splice_is_a_fixed_point():
    """Re-splicing generated text must change nothing.

    Load-bearing now that the text lives in the source: the splicers still run at
    every import, so a non-idempotent one would grow each docstring on every
    reload and make the generated file wrong the moment it was read back. It WAS
    non-idempotent -- `ModelPlots.section` went 169 lines to 179 on a second
    pass, duplicating the extended prose, the `**kwargs` entry (a NumPy head with
    no type, invisible to the name filter) and the `Bound form of` footer.
    """

    import myflopy.plot as plot
    from myflopy.layers import StackPlots
    from myflopy.plot import GridPlots, ModelPlots

    for namespace in (ModelPlots, GridPlots, StackPlots):
        for verb in ("map", "section", "surface", "grid", "mosaic", "animate"):
            method = getattr(namespace, verb, None)
            if method is None:
                continue
            before = method.__doc__
            if hasattr(method, "_myflopy_own_doc"):
                del method._myflopy_own_doc
            method.__doc__ = before
            plot._inherit_verb_docs(namespace)
            assert (method.__doc__ or "").strip() == (before or "").strip(), (
                f"{namespace.__name__}.{verb} grew when re-spliced: "
                f"{len((before or '').splitlines())} -> "
                f"{len((method.__doc__ or '').splitlines())} lines"
            )
