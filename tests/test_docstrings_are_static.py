"""Every spliced docstring must be in its SOURCE, and the splice must still work.

Reported as "model.plot.section docstring doesn't show in PyCharm... it shows the
args but the docstring itself is only one line. This is a reoccurring issue."
PyCharm and Pylance are STATIC -- they read the `def` and the literal under it
and never run the module -- so a `__doc__` assembled at import reaches `help()`
and reaches no editor. `scripts/derive_docstrings.py` writes the assembled text
back into each source file.

The FIRST version of that script closed the reported defect and broke the splice
doing it, which is why this file is shaped the way it is. Writing the assembled
text into the source made it the method's own on the next import; its entries
then shadowed every inherited one and the free verbs reached nothing at all --
29 of 29 methods -- while `--check` reported clean, because source and runtime
agreed on the same stale text. Three docstrings told readers that `reverse`
"applies to BOTH branches" for a week after that was measured false.

So a test that only checks source == runtime is not enough: both can be stale
together. `test_an_edit_to_a_reference_reaches_the_bound_copies` is the one that
fails when the mechanism dies, and it is the reason the others can be trusted.
"""

from __future__ import annotations

import ast
import importlib
import inspect
import json
import pathlib
import subprocess
import sys

import pytest

REPO = pathlib.Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "derive_docstrings.py"
SHADOWS = REPO / "tests" / "docstring_shadows.json"


def _targets():
    """The spliced methods, discovered exactly as the generator discovers them."""

    sys.path.insert(0, str(REPO / "scripts"))
    import derive_docstrings

    return derive_docstrings.targets()


TARGETS = _targets()
IDS = [f"{cls}.{verb}" for _module, cls, verb in TARGETS]

#: Pinned, in both directions. Discovery that finds nothing passes everything,
#: and discovery that quietly loses one method is exactly how a docstring goes
#: stale unwatched -- so a change to this number has to be deliberate.
EXPECTED = 29


def test_the_sweep_finds_every_spliced_method_and_no_others():
    assert len(TARGETS) == EXPECTED, (
        f"{len(TARGETS)} spliced methods found, expected {EXPECTED}. If a verb "
        f"or noun was added or removed, update EXPECTED in the same commit."
    )
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

    method = getattr(getattr(importlib.import_module(module_name), cls_name), verb)
    doc = inspect.cleandoc(method.__doc__ or "")
    assert len(doc.splitlines()) > 10, (
        f"{cls_name}.{verb} has a {len(doc.splitlines())}-line docstring"
    )
    assert "Parameters" in doc, f"{cls_name}.{verb} documents no parameters"


@pytest.mark.parametrize(("module_name", "cls_name", "verb"), TARGETS, ids=IDS)
def test_every_generated_region_is_fenced_on_both_sides(module_name, cls_name, verb):
    """An unbalanced fence makes `strip_generated` eat the rest of the text.

    Reads the SOURCE literal, not `__doc__`. `_fence` balances its own output by
    construction, so a runtime check can never fail -- this test was written for
    exactly the case it could not see. Measured on a tree with one closing fence
    deleted by hand: source balance -1, runtime balance 0, test green, and the
    next regeneration silently replaced 31 lines of `GridPlots.section`'s own
    text with the free verb's model-scope entries.
    """

    from myflopy.plot import _FENCE_CLOSE, _FENCE_OPEN_PREFIX

    cls = getattr(importlib.import_module(module_name), cls_name)
    tree = ast.parse(pathlib.Path(inspect.getsourcefile(cls)).read_text())
    doc = None
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == cls_name:
            for fn in node.body:
                if isinstance(fn, ast.FunctionDef) and fn.name == verb:
                    doc = inspect.cleandoc(ast.get_docstring(fn) or "")
    assert doc is not None, f"{cls_name}.{verb} has no source docstring"
    depth = 0
    for line in doc.splitlines():
        stripped = line.strip()
        if stripped.startswith(_FENCE_OPEN_PREFIX):
            assert depth == 0, f"{cls_name}.{verb}: nested opening fence"
            depth = 1
        elif stripped == _FENCE_CLOSE:
            assert depth == 1, f"{cls_name}.{verb}: closing fence with none open"
            depth = 0
    assert depth == 0, f"{cls_name}.{verb}: a fence is never closed"
    assert _FENCE_OPEN_PREFIX in doc, f"{cls_name}.{verb} inherits nothing"


def test_a_marker_in_the_content_itself_is_dropped_not_honoured():
    """The markers are RESERVED, and nothing used to check that.

    A `.. end generated` inside a free verb's Notes is written into a fenced
    region and then closes it early on the next recovery, so the tail is
    reclassified as the method's own, re-inherited beside itself, and the
    docstring grows by two lines on every regeneration. Dropped with a warning
    rather than raised: the splicers run at import, and raising there takes
    `import myflopy` down -- including for `derive_docstrings.py`, which imports
    the module before it can repair anything.
    """

    import myflopy.plot as plot

    scrubbed = plot._scrub_markers(
        "Notes\n-----\nSee below.\n.. end generated\nMore.", "section"
    )
    assert ".. end generated" not in scrubbed
    assert "See below." in scrubbed and "More." in scrubbed

    # And no reference carries one, which is what would put it into the pipeline.
    for verb in ("map", "section", "surface", "grid", "mosaic", "animate"):
        doc = getattr(plot, verb).__doc__ or ""
        assert plot._FENCE_OPEN_PREFIX not in doc, f"plot.{verb} contains a fence marker"
        assert plot._FENCE_CLOSE not in doc, f"plot.{verb} contains a fence marker"


def test_a_backslash_in_a_reference_survives_the_round_trip():
    """Written plain, a backslash in a `\"\"\"` literal is read as an escape.

    So a reference mentioning a Windows path, a regex or ``"\\n"`` would come
    back from the source DIFFERENT from what was written -- and `--check` would
    then certify the corruption, because the mangled text is the new fixed point.
    """

    sys.path.insert(0, str(REPO / "scripts"))
    import derive_docstrings

    text = 'Summary.\n\nUse a path like C:\\Users\\me, or the escape \\n.'
    literal = "\n".join(derive_docstrings._literal(text, 4))
    parsed = ast.get_docstring(ast.parse(f"def f():\n{literal}\n    pass").body[0])
    assert parsed == text, f"round trip corrupted the backslash:\n{parsed!r}"


def test_an_edit_to_a_reference_reaches_the_bound_copies():
    """THE test. It is the one that was missing, and its absence cost a week.

    Everything else here compares the source to the runtime docstring, and the
    defect made both stale together, so all of them passed -- 61 green tests over
    a mechanism that propagated nothing. This edits a free verb's reference in
    memory, re-runs the splice, and demands the change actually arrive.

    Both directions matter. The new entry must appear; the old text must be
    GONE. Under the previous design an edit to a free verb's prose appended a
    second copy beside the stale one rather than replacing it, and three edits
    left `GridPlots.map` at 211 lines carrying three contradictory copies.
    """

    import myflopy.plot as plot

    original = plot.section.__doc__
    sentinel = "PROPAGATION-SENTINEL-9c1f."
    try:
        plot.section.__doc__ = original.replace(
            "The section line as a geometry", sentinel + " The section line as a geometry", 1
        )
        assert sentinel in plot.section.__doc__, "the reference edit did not apply"
        plot._inherit_verb_docs(plot.ModelPlots)
        reached = inspect.cleandoc(plot.ModelPlots.section.__doc__ or "")
        assert sentinel in reached, (
            "an edit to `myflopy.plot.section` did NOT reach ModelPlots.section. "
            "The splice has stopped propagating -- exactly the defect this file "
            "documents. Check that `strip_generated` still recovers the own text."
        )
        assert reached.count(sentinel) == 1, (
            f"the edit arrived {reached.count(sentinel)} times: the splice is "
            f"appending rather than replacing."
        )
    finally:
        plot.section.__doc__ = original
        plot._inherit_verb_docs(plot.ModelPlots)
    assert sentinel not in (plot.ModelPlots.section.__doc__ or ""), (
        "reverting the reference did not revert the bound copy"
    )


def test_hand_written_text_inside_a_fence_does_not_survive():
    """The ratchet direction that had no teeth at all.

    Ledger 178 claimed `--check` "fails in BOTH directions, so hand-editing a
    generated docstring fails too". Measured false: with the source acting as its
    own author, a hand edit to generated text WAS a fixed point, and injecting a
    sentence into `ModelPlots.map`'s inherited block left `--check` at exit 0.
    """

    import myflopy.plot as plot

    doc = inspect.cleandoc(plot.ModelPlots.map.__doc__ or "")
    tampered = doc.replace(
        "Parameters\n----------", "Parameters\n----------\nHAND-EDIT.", 1
    )
    assert plot.strip_generated(tampered) != tampered
    # The tampered line sits inside the fenced run, so recovering the own text
    # discards it -- which is what makes the next regeneration disagree with the
    # source and fail `--check`.
    assert "HAND-EDIT." not in plot.strip_generated(tampered)


def test_the_splice_is_a_fixed_point():
    """Re-splicing generated text must change nothing.

    Load-bearing, because the splicers run at EVERY import: a non-idempotent one
    would grow each docstring on every reload and make the written source wrong
    the moment it was read back. It used to be guarded by containment checks
    that held only while the reference never changed; now it is structural --
    `strip_generated` deletes the previous splice's output before re-splicing,
    so there is nothing left to duplicate.
    """

    import myflopy.plot as plot
    import myflopy.project.group.lak as glak
    from myflopy.layers import StackPlots

    before = {
        (m, c, v): getattr(importlib.import_module(m), c).__dict__[v].__doc__
        for m, c, v in TARGETS
    }
    for _ in range(3):
        for namespace in (plot.ModelPlots, plot.GridPlots, StackPlots):
            plot._inherit_verb_docs(namespace)
        plot._inherit_noun_docs()
        plot.inherit_map_docs(
            glak.GroupLakConnections.map, keep=set(plot.NOUN_MAP_PARAMS)
        )

    grew = [
        f"{c}.{v}"
        for (m, c, v), was in before.items()
        if (getattr(importlib.import_module(m), c).__dict__[v].__doc__ or "") != (was or "")
    ]
    assert not grew, f"re-splicing changed {len(grew)} docstrings: {grew}"


def test_the_shadow_ledger_is_committed():
    """The names a reference edit will NOT reach, pinned so they stay auditable.

    A shadow is deliberate -- `layer` is narrower on a record noun, `reverse` is
    mpl-only on a bare grid -- and permanent: the method's own entry wins. What
    cannot be detected is a shadow whose prose has quietly become untrue, which
    is what the three `reverse` entries were. Pinning the SET at least makes a
    new one visible in review.
    """

    sys.path.insert(0, str(REPO / "scripts"))
    import derive_docstrings

    assert SHADOWS.exists(), "run `python scripts/derive_docstrings.py`"
    committed = json.loads(SHADOWS.read_text())
    assert committed == derive_docstrings.shadow_ledger(), (
        "the shadow ledger has drifted. Run "
        "`python scripts/derive_docstrings.py` and READ the diff: every new "
        "entry is a name whose reference edits will stop reaching that method."
    )


def test_prose_inside_a_parameter_list_is_refused():
    """A paragraph at column 0 inside `Parameters` is read as an entry head.

    `GridPlots.map` carried two, and a NumPy reader showed the first as a
    parameter named ``Requires a CRS`` whose type was ``because the basemap
    does. Use :meth:`grid` for a grid``. The SPLICE created them, by folding the
    method's trailing prose into the merged body -- so the same mistake in a
    REFERENCE would propagate to every copy that inherits it.
    """

    import myflopy.plot as plot

    bad = (
        "Parameters\n----------\n"
        "line : LineString\n    The section line.\n"
        "Requires a CRS, because the basemap does.\n"
    )
    with pytest.raises(ValueError, match="not an entry head"):
        plot._refuse_malformed_entries(bad, "a test")

    ok = "Parameters\n----------\nzmin, zmax : float, optional\n    Limits.\n"
    plot._refuse_malformed_entries(ok, "a test")
    plot._refuse_malformed_entries(
        "Parameters\n----------\n**kwargs\n    Passed through.\n", "a test"
    )


#: Nouns that accept `**trace_kwargs` without documenting it. PRE-EXISTING and
#: left alone -- the noun tier deliberately says less about the Plotly tail than
#: the verb tier does. Pinned so the number cannot grow unnoticed.
UNDOCUMENTED_TAILS = 15


def test_every_spliced_method_documents_what_it_accepts():
    """The direction the bound tier had no test for at all.

    `test_plot_vocabulary` covers the six free verbs and `test_noun_signatures`
    covers the nouns; nothing checked the bound namespaces, so `fill=` could have
    vanished from `model.plot.section`'s docstring with the whole suite green.
    """

    import myflopy.plot as plot

    gaps = {}
    for module_name, cls_name, verb in TARGETS:
        method = getattr(importlib.import_module(module_name), cls_name).__dict__[verb]
        doc = inspect.cleandoc(method.__doc__ or "")
        # `inspect` reports a VAR-KEYWORD as `trace_kwargs`; NumPy documents it
        # as `**trace_kwargs`. Comparing across the stars invents 25 phantom gaps.
        documented = {
            n.lstrip("*")
            for n in plot._documented_names(plot._split_sections(doc)[1])
        }
        accepted = set(inspect.signature(method).parameters) - {"self"}
        if accepted - documented:
            gaps[f"{cls_name}.{verb}"] = sorted(accepted - documented)

    assert sum(len(v) for v in gaps.values()) == UNDOCUMENTED_TAILS, (
        f"accepted-but-undocumented parameters changed: {gaps}"
    )


def test_a_combinators_first_positional_is_not_mistaken_for_a_subject():
    """`mosaic` and `animate` have no subject, and their first argument is real.

    The splice drops the free verb's first positional because binding supplies it
    -- `plot.map(model, ...)` against `model.plot.map(...)`. Deriving that by
    POSITION alone silently undocumented `panels` on `ModelPlots.mosaic` and
    `frames` on `ModelPlots.animate`, which are combinators over pictures and
    take no subject at all. The rule is what the bound method ACCEPTS, not where
    the parameter sits.
    """

    import myflopy.plot as plot

    def documented(method):
        doc = inspect.cleandoc(method.__doc__ or "")
        return plot._documented_names(plot._split_sections(doc)[1])

    assert "panels" in documented(plot.ModelPlots.mosaic)
    assert "frames" in documented(plot.ModelPlots.animate)
    # ...and the genuine subject still goes.
    assert "source" not in documented(plot.ModelPlots.map)
    assert "source" not in documented(plot.ModelPlots.section)
