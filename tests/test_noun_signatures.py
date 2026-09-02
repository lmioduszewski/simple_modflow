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

#: Noun classes still carrying the pre-8.8 shape. EXACT IN BOTH DIRECTIONS: a
#: class fixed without being removed fails just as loudly as a new one added, so
#: the list can only shrink and cannot quietly absorb a regression.
#:
#: Every entry is debt, not an exemption. See ledger 169.
UNCONVERTED = {
    "CellPackageInputFieldExplorer.map",
    "StaticArrayFieldExplorer.map",
    "UzfFieldInputsExplorer.map",
    "HfbResultsExplorer.map",
    "HfbPackageExplorer.map",
    "CellBudgetResultsExplorer.map",
    "StageResultsExplorer.map",
    "LakBudgetResultsExplorer.map",
    "LakConnectionsExplorer.map",
    "SfrBudgetResultsExplorer.map",
    "SurfaceWaterExchangeResultsExplorer.map",
    "SurfaceWaterInputFieldExplorer.map",
    "PRTPathlineView.map",
    "UzfInput.map",
    "GroupLakConnections.map",
    "GridPlots.map",
    "StackPlots.section",
}


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


def test_the_sweep_finds_the_noun_tier_at_all():
    """A guard on the guard: a discovery test that discovers nothing passes
    everything."""

    assert len(DISCOVERED) >= 15, f"the sweep found only {len(DISCOVERED)} methods"
    assert "DependentVariableFile.map" in IDS


def test_the_unconverted_list_is_exact_in_both_directions():
    """Debt shrinks or the list is wrong.

    Ledger 145 recorded this tier as measured-clean and eleven months of nouns
    were written to a convention nothing checked. An allowlist that can silently
    absorb a new entry would repeat exactly that.
    """

    discovered = set(IDS)
    stale = UNCONVERTED - discovered
    assert not stale, f"listed as unconverted but no longer found: {sorted(stale)}"


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
        pytest.xfail(f"{name} is recorded debt (ledger 169)")

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


@pytest.mark.parametrize(("name", "rel", "line", "fn", "source"), DISCOVERED, ids=IDS)
def test_a_noun_verb_offers_the_measured_parameter_set(name, rel, line, fn, source):
    """8.8 requirement 2: name what the destination accepts.

    "A narrower scope may offer FEWER parameters, never other ones" -- so a noun
    must not silently forward a drawing parameter it never names.
    """

    if name in UNCONVERTED:
        pytest.xfail(f"{name} is recorded debt (ledger 169)")

    named = {a.arg for a in fn.args.kwonlyargs}
    missing = [p for p in NOUN_MAP_PARAMS if p not in named]
    assert not missing, f"{rel}:{line} {name} forwards but never names {missing}"


@pytest.mark.parametrize(("name", "rel", "line", "fn", "source"), DISCOVERED, ids=IDS)
def test_a_noun_verb_offers_nothing_meaningless(name, rel, line, fn, source):
    """The rule pointed the other way: naming a parameter that cannot act is the
    same defect as hiding one that can.

    `values`/`type` contradict what a noun is; the legacy hover trio is dead
    against a noun's own `hover_spec`; and the inert four were measured to change
    nothing on either a results noun or a record noun.
    """

    if name in UNCONVERTED:
        pytest.xfail(f"{name} is recorded debt (ledger 169)")

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
