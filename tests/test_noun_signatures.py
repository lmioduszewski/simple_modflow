"""The 8.8 no-bare-kwargs rule, extended to the NOUN tier (2026-08-30).

`model.hds.map()` and `model.packages.drn.inputs.elev.map()` are nouns: narrower
than the verbs, because the noun already fixes WHAT is drawn. Plan 8.8 required
every picture verb to name its parameters -- editors read the `def` line and
never run the module -- and `tests/test_plot_vocabulary.py` enforces that over
`SCOPES`, which is `{module, model, grid, stack}`. Every enumeration in that file
is a hand-written four-entry whitelist, so the noun tier is unreachable by
construction, in both directions.

It rotted accordingly. Reported as "model.hds.map() docstring doesn't show all
arguments"; measured, 29 of its 39 reachable parameters were invisible, including
`show_mounding`, which does not decorate the map but REWRITES it (z from
153.76/147.67/144.27 to 22.00/16.28/15.79).

Ledger 145 had excused this tier on one question -- "is the signature
`(*args, **kwargs)`?" -- where the rule states four requirements. This file asks
all four, and discovers its subjects instead of listing them.
"""

from __future__ import annotations

import ast
import inspect
import pathlib
import re

import pytest

from myflopy.plot import (
    LAYER_FIELD_MAP_PARAMS,
    NOUN_INERT_PARAMS,
    NOUN_MAP_PARAMS,
    NOUN_REFUSED_PARAMS,
)
from myflopy.plot import map as free_map

SRC = pathlib.Path(__file__).resolve().parent.parent / "src" / "myflopy"

#: Picture verbs, per the vocabulary. Imported shape rather than a literal so a
#: seventh verb reaches this tier too.
VERBS = ("map", "section", "surface", "grid", "mosaic", "animate")

#: Noun classes still carrying the pre-8.8 shape. **EMPTY as of 2026-09-02** --
#: all eighteen are converted (ledger 169b, then 170).
#:
#: Kept rather than deleted, and kept EXACT IN BOTH DIRECTIONS: a class fixed
#: without being removed fails as loudly as a new one added. Emptying it turns
#: the three parametrized tests below from "xfail on the debt list" into "assert
#: on every noun in `src/`", which is the state this file was written to reach.
#: Anything added back here is debt with a ledger entry, not an exemption.
UNCONVERTED: set[str] = set()

#: The bound VERB namespaces. The sweep finds them too -- they define a picture
#: verb with a ``**kwargs`` tail -- but they are the verb tier wearing a noun's
#: shape, and the two tiers are governed by DIFFERENT rules.
#:
#: A noun fixes what is drawn, so `values=`/`type=`/`custom_hover=` contradict
#: it and `NOUN_REFUSED_PARAMS` makes them raise. A scope fixes only the
#: SUBJECT: `vor.plot.map(values=node_ids)` is the whole point of a grid map,
#: and `custom_hover` is the only hover a bare grid has, since the sectioned one
#: needs a model. Applying the noun rules here would delete three working
#: parameters from `vor.plot.map` to satisfy a rule written about something else.
#:
#: They are NOT thereby unchecked. `test_plot_vocabulary.py` holds them to the
#: verb rule -- exact mirror of the free verb for `model`, a subset for `grid`
#: and `stack`, mirrored defaults, and no bare-kwargs forwarder -- and the
#: docstring test below still covers them, because "document what you name"
#: is a requirement of both tiers.
VERB_SCOPES = {"ModelPlots", "GridPlots", "StackPlots"}


def _forwards_to_plot(source: str) -> bool:
    """Whether a method body hands off to a plot function."""

    return any(
        marker in source
        for marker in ("plot.map(", "plot.section(", "_choropleth_factory")
    )


def _discovered() -> list[tuple[str, str, int, ast.FunctionDef, str]]:
    """Every class-level picture verb in src/ that forwards to a plot function.

    A STATIC sweep, deliberately: a runtime walk can only reach the nouns a
    canonical model happens to instantiate, and `PRTPathlineView` needs a
    completed PRT run while `StaticArrayFieldExplorer` is reachable only through
    `declared_fields`. Discovery is the whole point -- a new explorer class must
    not be able to appear uncovered.
    """

    found = []
    for path in sorted(SRC.rglob("*.py")):
        if "_vendor" in path.parts:
            continue
        text = path.read_text()
        try:
            tree = ast.parse(text)
        except SyntaxError:                                   # pragma: no cover
            continue
        for node in ast.walk(tree):
            if not isinstance(node, ast.ClassDef):
                continue
            for fn in node.body:
                if not isinstance(fn, ast.FunctionDef) or fn.name not in VERBS:
                    continue
                source = ast.get_source_segment(text, fn) or ""
                if fn.args.kwarg is not None and _forwards_to_plot(source):
                    rel = str(path.relative_to(SRC.parent))
                    found.append((f"{node.name}.{fn.name}", rel, fn.lineno, fn, source))
    return found


DISCOVERED = _discovered()
IDS = [d[0] for d in DISCOVERED]

#: The noun half of the sweep: everything that is not a bound verb namespace.
NOUNS = [d for d in DISCOVERED if d[0].split(".")[0] not in VERB_SCOPES]
NOUN_IDS = [d[0] for d in NOUNS]


def test_the_sweep_finds_the_noun_tier_at_all():
    """A guard on the guard: a discovery test that discovers nothing passes
    everything."""

    assert len(DISCOVERED) >= 15, f"the sweep found only {len(DISCOVERED)} methods"
    assert "DependentVariableFile.map" in IDS
    # Both halves must be non-empty, or the split silently exempts a whole tier.
    assert NOUNS, "the noun half of the sweep is empty"
    assert len(DISCOVERED) > len(NOUNS), "no verb scope was found to separate"


def test_the_unconverted_list_is_exact_in_both_directions():
    """Debt shrinks or the list is wrong.

    Ledger 145 recorded this tier as measured-clean and eleven months of nouns
    were written to a convention nothing checked. An allowlist that can silently
    absorb a new entry would repeat exactly that.
    """

    discovered = set(IDS)
    stale = UNCONVERTED - discovered
    assert not stale, f"listed as unconverted but no longer found: {sorted(stale)}"
    # The debt is closed. Re-adding a name is allowed -- that is what the list is
    # for -- but it must be a deliberate act with a ledger entry behind it, not
    # the quiet landing place for a noun someone did not want to convert.
    assert not UNCONVERTED, (
        f"the noun tier was fully converted on 2026-09-02; {sorted(UNCONVERTED)} "
        f"is back on the debt list. If that is intended, record it in "
        f"docs/compromises_and_deferrals.md and update this assertion."
    )


def _resolve(name: str, rel: str):
    """The live method for a discovered ``Class.verb``, or None if unimportable."""

    import importlib

    module = importlib.import_module(rel[:-3].replace("/", ".").replace("\\", "."))
    cls_name, verb = name.split(".")
    cls = getattr(module, cls_name, None)
    return None if cls is None else getattr(cls, verb, None)


@pytest.mark.parametrize(("name", "rel", "line", "fn", "source"), DISCOVERED, ids=IDS)
def test_a_noun_verb_documents_every_parameter_it_names(name, rel, line, fn, source):
    """8.8 requirement 1: a parameter you accept is a parameter you document.

    Checked on the LIVE ``__doc__``, not the source docstring, because that is
    the house contract the verb tier already runs on: the reference is written
    once on the free verb and spliced onto every bound copy at import
    (`_inherit_verb_docs`), so no bound method repeats it in source either. What
    a purely static reader gets is the SIGNATURE, which is why the sibling test
    below asserts that separately -- and why it is the load-bearing one.
    """

    if name in UNCONVERTED:
        pytest.xfail(f"{name} is recorded debt (ledger 169b)")

    method = _resolve(name, rel)
    if method is None:                                        # pragma: no cover
        pytest.skip(f"{name} is not importable from {rel}")
    doc = method.__doc__ or ""
    named = [a.arg for a in fn.args.kwonlyargs]
    assert "Parameters" in doc, f"{rel}:{line} {name} documents nothing it accepts"

    # An ENTRY, not a mention. `name in doc` passes when the parameter merely
    # appears in the prose or a Raises clause -- which is how a docstring
    # carrying one real entry (`backend`) and 30 names in passing read as fully
    # documented. Require `name : type` at column 0 with a description under it.
    # NumPy allows a GROUPED head -- `zmin, zmax : float, optional` documents two
    # parameters in one entry -- so split the head rather than assuming one name.
    entries = {}
    for match in re.finditer(r"^([\w, ]+?) :.*\n((?:[ \t]+\S.*\n?)*)", doc, re.M):
        body = match.group(2).strip()
        for param in (n.strip() for n in match.group(1).split(",")):
            if param:
                entries[param] = body
    missing = [a for a in named if a not in entries]
    assert not missing, f"{rel}:{line} {name} names but never documents {missing}"
    bodyless = [a for a in named if not entries[a]]
    assert not bodyless, f"{rel}:{line} {name} lists {bodyless} with no description"


@pytest.mark.parametrize(("name", "rel", "line", "fn", "source"), NOUNS, ids=NOUN_IDS)
def test_a_noun_verb_offers_the_measured_parameter_set(name, rel, line, fn, source):
    """8.8 requirement 2: name what the destination accepts.

    "A narrower scope may offer FEWER parameters, never other ones" -- so a noun
    must not silently forward a drawing parameter it never names.
    """

    if name in UNCONVERTED:
        pytest.xfail(f"{name} is recorded debt (ledger 169b)")

    named = {a.arg for a in fn.args.kwonlyargs}
    missing = [p for p in NOUN_MAP_PARAMS if p not in named]
    assert not missing, f"{rel}:{line} {name} forwards but never names {missing}"


@pytest.mark.parametrize(("name", "rel", "line", "fn", "source"), NOUNS, ids=NOUN_IDS)
def test_a_noun_verb_offers_nothing_meaningless(name, rel, line, fn, source):
    """The rule pointed the other way: naming a parameter that cannot act is the
    same defect as hiding one that can.

    `values`/`type` contradict what a noun is; the legacy hover trio is dead
    against a noun's own `hover_spec`; and the inert four were measured to change
    nothing on either a results noun or a record noun.
    """

    if name in UNCONVERTED:
        pytest.xfail(f"{name} is recorded debt (ledger 169b)")

    named = {a.arg for a in fn.args.kwonlyargs}
    offered = named & (set(NOUN_REFUSED_PARAMS) | set(NOUN_INERT_PARAMS))
    assert not offered, f"{rel}:{line} {name} offers parameters it cannot honour: {sorted(offered)}"


def test_the_tiers_partition_the_free_verb():
    """The tiers are a statement ABOUT `plot.map`, so they must stay true to it."""

    verb = {
        n for n, p in inspect.signature(free_map).parameters.items()
        if p.kind is p.KEYWORD_ONLY
    }
    declared = (
        set(NOUN_MAP_PARAMS) | set(LAYER_FIELD_MAP_PARAMS)
        | set(NOUN_REFUSED_PARAMS) | set(NOUN_INERT_PARAMS)
    )
    # `per`/`layer` are the noun's own selectors, named by every noun already;
    # `values` is positional-or-keyword on the verb, so it is not in `verb`.
    unaccounted = verb - declared - {"per", "layer"}
    assert not unaccounted, f"plot.map parameters in no tier: {sorted(unaccounted)}"
    invented = declared - verb - {"values"}
    assert not invented, f"tiers name parameters plot.map does not have: {sorted(invented)}"


@pytest.mark.parametrize(("name", "rel", "line", "fn", "source"), DISCOVERED, ids=IDS)
def test_a_noun_verbs_examples_name_paths_that_exist(name, rel, line, fn, source,
                                                     canonical_run):
    """An Examples block is a promise that those calls work. Check the paths.

    Written after two of the fifteen conversions shipped an Examples block
    addressed to a noun that does not exist at that path:
    ``model.packages.lak.inputs.connections`` (the explorer hangs off the
    PACKAGE, not off ``.inputs``) and ``model.surface_water.results.exchange``
    (it is ``model.packages.surface_water.results.q``). Both were plausible,
    both were wrong, and neither the signature tests nor a docstring-completeness
    check could see them -- an example is prose to every one of those.

    Only the ATTRIBUTE CHAIN is walked, not the call: resolving
    ``model.packages.drn.inputs.elev`` proves the path, where actually drawing
    eighteen maps with every parameter combination in every docstring would
    turn this file into the slowest in the suite for a much smaller gain.
    """

    method = _resolve(name, rel)
    if method is None:                                        # pragma: no cover
        pytest.skip(f"{name} is not importable from {rel}")
    doc = method.__doc__ or ""
    assert "Examples" in doc, f"{rel}:{line} {name} shows no example calls"

    broken = []
    for example in re.findall(r"^\s*>>> (.+)$", doc, re.M):
        chain = re.match(r"(model|m)((?:\.[a-z_0-9]+)+)\(", example)
        if not chain:
            continue                      # `vor.`/`run.`/`stack.` need other subjects
        subject, walked = canonical_run, []
        for part in chain.group(2).lstrip(".").split("."):
            if part in VERBS or part in ("plot", "get", "summary"):
                break                     # the verb itself; the path is what we check
            walked.append(part)
            try:
                subject = getattr(subject, part)
            except AttributeError as error:
                broken.append(f"model.{'.'.join(walked)} ({error})")
                break
    assert not broken, (
        f"{rel}:{line} {name} shows examples addressed to paths that do not "
        f"exist: {broken}"
    )


#: NumPy's section order. A docstring that lists them out of order still parses,
#: but no renderer lays it out the way its author meant.
_CANONICAL_ORDER = (
    "Parameters", "Other Parameters", "Returns", "Yields", "Raises", "Warns",
    "See Also", "Notes", "References", "Examples",
)


def _sections(doc: str) -> list[str]:
    """Section headings in the order they appear (a title with a dashed rule)."""

    lines = doc.splitlines()
    return [
        line.strip()
        for i, line in enumerate(lines[:-1])
        if line.strip() in _CANONICAL_ORDER
        and lines[i + 1].strip()
        and set(lines[i + 1].strip()) == {"-"}
    ]


@pytest.mark.parametrize(("name", "rel", "line", "fn", "source"), DISCOVERED, ids=IDS)
def test_a_spliced_docstring_is_one_well_formed_document(name, rel, line, fn, source):
    """The splice must MERGE the two docstrings, not concatenate them.

    `_merge_parameter_sections` folded the `Parameters` blocks and then appended
    everything else from both sides, so 7 of 8 nouns carried two `Returns`, two
    `See Also` and two `Examples`. Worse than untidy: the LAST examples a reader
    reaches were the free verb's, so `help(model.packages.ghb.inputs.cond.map)`
    ended on `>>> vor.plot.map(values=node_ids)` -- a different subject, and a
    call this noun REFUSES. The fix also has to order the result, because
    own-then-inherited put an inherited `Returns` after the method's own
    `Examples`.
    """

    method = _resolve(name, rel)
    if method is None:                                        # pragma: no cover
        pytest.skip(f"{name} is not importable from {rel}")
    heads = _sections(method.__doc__ or "")

    repeated = sorted({h for h in heads if heads.count(h) > 1})
    assert not repeated, f"{rel}:{line} {name} has duplicate sections: {repeated}"

    ranks = [_CANONICAL_ORDER.index(h) for h in heads]
    assert ranks == sorted(ranks), (
        f"{rel}:{line} {name} sections are out of NumPy order: {heads}"
    )


@pytest.mark.parametrize(("name", "rel", "line", "fn", "source"), NOUNS, ids=NOUN_IDS)
def test_a_noun_shows_examples_of_ITSELF(name, rel, line, fn, source):
    """The last thing a reader sees must be a call to THIS noun.

    A noun that writes no `Examples` of its own inherits `plot.map`'s, which are
    verb-scope calls (`model.plot.map(values=...)`, `vor.plot.map(...)`) -- and
    `values=` is a parameter every noun raises on. Documenting a call that
    raises is worse than documenting nothing.
    """

    method = _resolve(name, rel)
    if method is None:                                        # pragma: no cover
        pytest.skip(f"{name} is not importable from {rel}")
    doc = method.__doc__ or ""
    assert "Examples" in _sections(doc), f"{rel}:{line} {name} shows no examples"

    examples = [ln.strip() for ln in doc.splitlines() if ln.strip().startswith(">>> ")]
    assert examples, f"{rel}:{line} {name} has an empty Examples block"
    verb_scope = [e for e in examples if ".plot.map(" in e or ".plot.section(" in e]
    assert not verb_scope, (
        f"{rel}:{line} {name} shows VERB-scope examples, which belong to "
        f"`myflopy.plot`, not to this noun: {verb_scope}"
    )
    # A WORD boundary, not a substring: `connection_type=` and
    # `lak_connection_type=` are real parameters of two of these nouns and both
    # end in `type=`. And a refused name the method actually ACCEPTS is not
    # refused -- `HfbPackageExplorer.map` takes `values` positionally on purpose,
    # because its subject is the barriers and the cells are a backdrop.
    accepted = set(inspect.signature(method).parameters)
    refusable = [p for p in NOUN_REFUSED_PARAMS if p not in accepted]
    pattern = re.compile(r"(?<![\w.])(" + "|".join(refusable) + r")\s*=") if refusable else None
    refused = [e for e in examples if pattern and pattern.search(e)]
    assert not refused, (
        f"{rel}:{line} {name} documents a call it would raise on: {refused}"
    )
