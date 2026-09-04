"""Write each spliced docstring back into its own source file.

`model.plot.section` used to carry ONE line of source docstring and 169 lines of
runtime `__doc__`, assembled at import by `_inherit_verb_docs`. PyCharm and
Pylance are static: they read the `def` and the literal beneath it and never run
the module, so a reader hovering the method got the one line. That is the same
defect plan 8.8 fixed for signatures -- `__doc__` set at import reaches `help()`
and reaches no editor -- and ledger 169b recorded the runtime form as the
contract without it ever being chosen.

This closes it the way the project already closes derived state: generate, then
pin. `tests/test_docstrings_are_static.py` re-runs the computation and fails if
the source has drifted, so the splice stays the single author and the source
stays what an editor can see.

Usage::

    python scripts/derive_docstrings.py           # rewrite
    python scripts/derive_docstrings.py --check   # report drift, change nothing

The splice must be IDEMPOTENT for this to be safe -- re-splicing its own output
has to be a fixed point, or every import would grow the text. It was not, and
three duplications had to be fixed first (the extended prose, a `**kwargs` entry
whose head has no type, and the `Bound form of` footer).
"""

from __future__ import annotations

import argparse
import ast
import importlib
import inspect
import pathlib
import sys

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))


def targets() -> list[tuple[str, str, str]]:
    """``(module, class, verb)`` for every method whose docstring is spliced.

    Discovered from the splice wiring itself rather than listed, so a noun added
    to `_inherit_noun_docs` cannot be left behind with a one-line docstring.
    """

    import myflopy.plot as plot  # noqa: F401 - importing runs every splice

    # Nouns ABOVE `myflopy.plot` in the layer graph splice themselves at the
    # bottom of their own module, so they are only stamped once imported.
    import myflopy.project.group.lak  # noqa: F401

    found: list[tuple[str, str, str]] = []
    seen: set[tuple[str, str]] = set()
    for module in list(sys.modules.values()):
        name = getattr(module, "__name__", "")
        if not name.startswith("myflopy.") or "_vendor" in name:
            continue
        for cls_name, cls in vars(module).items():
            if not inspect.isclass(cls) or cls.__module__ != name:
                continue
            for verb in ("map", "section", "surface", "grid", "mosaic", "animate"):
                method = cls.__dict__.get(verb)
                if not inspect.isfunction(method):
                    continue
                # `_myflopy_own_doc` is stamped by the splicers, and only by
                # them: its presence IS the marker that this docstring is
                # assembled rather than hand-written.
                if not hasattr(method, "_myflopy_own_doc"):
                    continue
                if (cls_name, verb) in seen:
                    continue
                seen.add((cls_name, verb))
                found.append((name, cls_name, verb))
    return sorted(found)


def _docstring_span(source: str, cls_name: str, verb: str) -> tuple[int, int, int]:
    """``(start, end, indent)`` 0-based line span of a method's docstring literal."""

    tree = ast.parse(source)
    for node in ast.walk(tree):
        if not (isinstance(node, ast.ClassDef) and node.name == cls_name):
            continue
        for fn in node.body:
            if not (isinstance(fn, ast.FunctionDef) and fn.name == verb):
                continue
            first = fn.body[0]
            if not (isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant)
                    and isinstance(first.value.value, str)):
                raise SystemExit(f"{cls_name}.{verb} has no docstring to replace")
            return first.lineno - 1, first.end_lineno, first.col_offset
    raise SystemExit(f"{cls_name}.{verb} not found")


def _literal(text: str, indent: int) -> list[str]:
    """The docstring as source lines, indented and safely quoted."""

    pad = " " * indent
    if '"""' in text or text.rstrip().endswith("\\"):
        raise SystemExit("generated docstring cannot be written as a \"\"\" literal")
    lines = text.strip().splitlines()
    out = [f'{pad}"""{lines[0]}']
    out += [f"{pad}{line}".rstrip() if line.strip() else "" for line in lines[1:]]
    out.append(f'{pad}"""')
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true",
                        help="report drift without writing")
    args = parser.parse_args()

    stale: list[str] = []
    written = 0
    for module_name, cls_name, verb in targets():
        module = importlib.import_module(module_name)
        method = getattr(module, cls_name).__dict__[verb]
        text = inspect.cleandoc(method.__doc__ or "")
        path = pathlib.Path(inspect.getsourcefile(getattr(module, cls_name)))
        source = path.read_text()
        start, end, indent = _docstring_span(source, cls_name, verb)
        lines = source.splitlines()
        replacement = _literal(text, indent)
        if lines[start:end] == replacement:
            continue
        stale.append(f"{path.relative_to(REPO)}::{cls_name}.{verb}")
        if not args.check:
            path.write_text("\n".join(lines[:start] + replacement + lines[end:]) + "\n")
            written += 1

    if args.check:
        if stale:
            print(f"{len(stale)} docstring(s) out of date:")
            print("\n".join(f"  {s}" for s in stale))
            print("\nRun: python scripts/derive_docstrings.py")
            return 1
        print("all spliced docstrings are written into their source")
        return 0
    print(f"rewrote {written} docstring(s)" if written else "nothing to rewrite")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
