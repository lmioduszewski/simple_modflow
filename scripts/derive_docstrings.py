"""Write each spliced docstring back into its own source file, and pin the shadows.

PyCharm and Pylance are STATIC -- they read the `def` and the literal beneath it
and never run the module -- so a docstring assembled at import reaches `help()`
and reaches no editor. `model.plot.section` carried ONE line of source docstring
against 169 of runtime `__doc__`. This closes it the way the project closes every
other derived thing: generate, then pin.

The first version of this script closed it and BROKE the splice doing so. Once
the assembled text was written into the source, the next import cached THAT as
the method's own docstring, its entries shadowed every inherited one, and the
free verbs stopped reaching anything -- 29 of 29 methods, with `--check` green
throughout, because source and runtime agreed on the same stale text. The fix is
`_FENCE_OPEN_PREFIX`: the generated regions are delimited IN the docstring, so
the splice recovers its own input by deleting them (`strip_generated`) instead of
by remembering it. Idempotency and the ratchet stop being guarded properties and
become structural ones.

Usage::

    python scripts/derive_docstrings.py           # rewrite
    python scripts/derive_docstrings.py --check   # report drift, change nothing

Two artifacts, both pinned by `tests/test_docstrings_are_static.py`: the source
docstrings themselves, and `tests/docstring_shadows.json` -- every name a method
documents locally that its reference ALSO documents. Those are exactly the places
a future edit to a reference will not reach, by design, and the only ones. The
stale `reverse` entries that exposed all of this were three such shadows that
nothing was watching.
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
                # `_myflopy_spliced` is stamped by the splicers, and only by
                # them: its presence IS the marker that this docstring is
                # assembled rather than hand-written. A BOOLEAN, deliberately --
                # its predecessor cached the pre-splice docstring, and once the
                # assembled text was written back into the source the cache held
                # that instead, which is how the splice came to propagate
                # nothing while `--check` stayed green.
                if not getattr(method, "_myflopy_spliced", False):
                    continue
                if (cls_name, verb) in seen:
                    continue
                seen.add((cls_name, verb))
                found.append((name, cls_name, verb))
    return sorted(found)


def source_targets() -> set[tuple[str, str]]:
    """``(class, verb)`` for every method whose SOURCE carries a generated fence.

    The second, independent opinion on which methods the pipeline manages. A
    marker sweep alone would reproduce the bug this design replaces, one level
    up: delete a splice call and the method quietly stops being regenerated
    while its source still looks managed. Runtime discovery alone misses the
    mirror case. `--check` compares the two and names the disagreement.
    """

    from myflopy.plot import _FENCE_OPEN_PREFIX

    found: set[tuple[str, str]] = set()
    for path in (REPO / "src" / "myflopy").rglob("*.py"):
        if "_vendor" in path.parts:
            continue
        try:
            tree = ast.parse(path.read_text())
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            if not isinstance(node, ast.ClassDef):
                continue
            for fn in node.body:
                if isinstance(fn, ast.FunctionDef) and _FENCE_OPEN_PREFIX in (
                    ast.get_docstring(fn) or ""
                ):
                    found.add((node.name, fn.name))
    return found


SHADOWS = REPO / "tests" / "docstring_shadows.json"


def shadow_ledger() -> dict[str, list[str]]:
    """``{"Class.verb": [names]}`` a method documents that its reference also does.

    A shadow is deliberate -- `layer` means something narrower on a record noun
    than on the free verb, and `reverse` is mpl-only on a bare grid. It is also
    permanent: the method's entry wins, so an edit to the reference's entry for
    that name reaches everything EXCEPT the methods listed here.

    Nothing can detect that a shadow's prose has quietly become untrue; that is
    the ceiling, and this does not pretend otherwise. What it does is make the
    set auditable -- a new shadow shows up as a diff in a pinned artifact at
    review time, which is the house pattern (`api_snapshot.json`). The three
    stale `reverse` entries were shadows nobody was watching.
    """

    import myflopy.plot as plot

    out: dict[str, list[str]] = {}
    for module_name, cls_name, verb in targets():
        module = importlib.import_module(module_name)
        method = getattr(module, cls_name).__dict__[verb]
        own = plot.strip_generated(inspect.cleandoc(method.__doc__ or ""))
        free = getattr(plot, verb if hasattr(plot, verb) else "map", None)
        if free is None or not free.__doc__:
            continue
        mine = plot._documented_names(plot._split_sections(own)[1])
        theirs = plot._documented_names(
            plot._split_sections(inspect.cleandoc(free.__doc__))[1]
        )
        shadowed = sorted(mine & theirs)
        if shadowed:
            out[f"{cls_name}.{verb}"] = shadowed
    return out


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
    # A RAW literal when the text carries a backslash. Written plain, a reference
    # that mentions a Windows path, a regex or `"\n"` has its backslash read as
    # an escape by the compiler, so `__doc__` comes back different from the file
    # -- silent, permanent corruption that `--check` then certifies, because the
    # mangled text is the new fixed point. (A raw literal cannot END in a
    # backslash; the guard above already refuses that.)
    prefix = "r" if "\\" in text else ""
    lines = text.strip().splitlines()
    out = [f'{pad}{prefix}"""{lines[0]}']
    out += [f"{pad}{line}".rstrip() if line.strip() else "" for line in lines[1:]]
    out.append(f'{pad}"""')
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true",
                        help="report drift without writing")
    args = parser.parse_args()

    discovered = targets()
    orphans = sorted(source_targets() - {(c, v) for _m, c, v in discovered})
    if orphans:
        # Source says generated, the splice says otherwise. The reverse case --
        # a newly added verb the splice manages whose source has no fences yet
        # -- is NOT an error: the ordinary drift check below catches it and says
        # to regenerate. Treating it as one deadlocks a new method behind a
        # message whose prescribed remedy is the thing it just refused to do.
        print("these carry generated fences but no splice writes them:")
        print("\n".join(f"  {cls}.{verb}" for cls, verb in orphans))
        print("\nEither wire them into `_inherit_verb_docs`/`_inherit_noun_docs`, "
              "or delete the fences and keep the text by hand.")
        return 1

    stale: list[str] = []
    written = 0
    for module_name, cls_name, verb in discovered:
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

    import json

    ledger = json.dumps(shadow_ledger(), indent=1, sort_keys=True) + "\n"
    ledger_stale = not SHADOWS.exists() or SHADOWS.read_text() != ledger
    if ledger_stale and not args.check:
        SHADOWS.write_text(ledger)

    if args.check:
        if ledger_stale:
            stale.append(str(SHADOWS.relative_to(REPO)))
        if stale:
            print(f"{len(stale)} docstring(s) out of date:")
            print("\n".join(f"  {s}" for s in stale))
            print("\nRun: python scripts/derive_docstrings.py")
            print("It REPLACES everything between `.. generated from ...` and "
                  "`.. end generated`, so text hand-written inside a fence is "
                  "discarded. Text outside the fences is yours and is untouched.")
            return 1
        print("all spliced docstrings are written into their source")
        return 0
    note = " and the shadow ledger" if ledger_stale else ""
    print(f"rewrote {written} docstring(s){note}" if written or ledger_stale
          else "nothing to rewrite")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
