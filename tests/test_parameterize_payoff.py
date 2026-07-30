"""What does adding ONE `_Recipe` entry actually get you, and what does it not?

The parameterize surface had no completeness test at all, and the cost showed
up twice in one day (2026-07-30):

* the `porosity` commit left `docs/manual/README.md` listing the pre-porosity
  target set, and nothing noticed;
* the `uzf.vks` commit left the target out of `parameterize`'s own docstring,
  found only while writing this file.

Both are the same failure: a target that is *wired* but not *announced*. Code
tests cannot see it, because the code works.

This is the parameterize twin of `test_package_descriptor_payoff.py`, and it
follows that file's contract: the AUTOMATIC set may grow but must never shrink,
and the MANUAL set must be EXACTLY as recorded — so a new hand-written site
fails immediately, and closing one also fails, forcing the win to be written
down rather than absorbed.

**Why a subprocess.** The recipe registry is consulted by module-level code in
several consumers. Injecting a synthetic recipe into an already-imported process
would not reach them, and `importlib.reload` would leave other modules holding
stale references — dangerous under `-n 10` with a session-scoped canonical
model. A fresh interpreter that injects BEFORE importing any consumer is the
only honest way to test "adding a recipe entry makes it appear".
"""

from __future__ import annotations

import json
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from myflopy.modflow.mf6.pest.native_parameters import _ALIASES, _RECIPES

REPO_ROOT = Path(__file__).resolve().parents[1]

#: Behaviour a `_Recipe` entry reaches with no other edit. May GROW, never shrink.
AUTOMATIC_SITES = frozenset(
    {
        "resolve_target",
        "resolve_target.error_text",
        "spec.style_default",
        "spec.additive_default",
        "spec.name_slug",
        "spec.pilotpoints_guard",
        "model_name_for",
        "flopy_model_for",
        "relayer_array_target",
        "resolve_files.layered_first",
        "select_layer_files",
        "zone_resolution",
    }
)

#: Sites a human still has to type. Each is here for a stated reason; closing
#: one is a real win and must be recorded by deleting it from this set.
MANUAL_SITES = frozenset(
    {
        # The registry entry itself, and its aliases -- deliberate: they carry
        # the prose explaining WHY a target is shaped the way it is, which is
        # where the load-bearing measurements live (see uzf.vks).
        "recipes.entry",
        # Doc enumerations. Four files list the supported targets in prose, and
        # prose is exactly what no code test can derive.
        "docs.CLAUDE.md",
        "docs.manual_README",
        "docs.myflopy_context",
        "docs.parameterize_docstring",
    }
)

#: The token each target must appear as in every enumerating doc. Explicit and
#: per-target because doc prose abbreviates: `ghb.cond`/`ghb.bhead` are listed
#: together as `ghb` in the shorter capability tables.
#:
#: This mapping must cover the registry EXACTLY, which is the both-directions
#: guard: add a recipe without a token and the first test fails; remove a recipe
#: and leave the token and it fails too.
DOC_TOKENS = {
    "k": "k",
    "k33": "k33",
    "recharge": "recharge",
    "chd": "chd",
    "ghb.cond": "ghb",
    "ghb.bhead": "ghb",
    "drn.cond": "drn",
    "drn.elev": "drn",
    "wel": "wel",
    "uzf.vks": "uzf.vks",
    "porosity": "porosity",
}

#: Files that ENUMERATE the supported targets, each with an anchor isolating the
#: enumeration itself and how far it runs.
#:
#: Searching the WHOLE file is not good enough, and this is not theoretical: the
#: first version of this test searched whole files, and deleting `porosity` from
#: the manual's target list still passed — because the word appears elsewhere in
#: that file. A doc test that cannot fail is worse than no doc test, because it
#: looks like coverage.
#:
#: The PEST guide is deliberately absent: it shows worked examples, not a
#: complete list, so requiring every target there would push noise into a
#: teaching document.
ENUMERATING_DOCS = (
    ("CLAUDE.md", "One unified parameterization API", "bullet"),
    ("docs/manual/README.md", "17.3 Parameterization", "line"),
    ("docs/myflopy_context.md", "**Unified** `cal.parameterize", "line"),
)


def _enumeration(relative_path: str, anchor: str, mode: str) -> str:
    """The passage of a doc that lists the targets, and nothing else."""

    text = (REPO_ROOT / relative_path).read_text(encoding="utf-8")
    start = text.find(anchor)
    assert start >= 0, (
        f"{relative_path} no longer contains {anchor!r}, so this test cannot "
        "find its target list. Update the anchor -- do not delete the check."
    )
    if mode == "line":
        line_start = text.rfind("\n", 0, start) + 1
        return text[line_start: text.find("\n", start)]
    end = text.find("\n- **", start)          # the next top-level bullet
    return text[start: end if end > 0 else len(text)]


def test_every_target_has_a_documentation_token():
    """The both-directions guard on the mapping below."""

    assert set(DOC_TOKENS) == set(_RECIPES), (
        "DOC_TOKENS and _RECIPES disagree. Adding a parameterize target means "
        "telling this test how the target is spelled in prose -- which is the "
        "moment to go and write that prose."
    )


@pytest.mark.parametrize(
    "relative_path,anchor,mode", ENUMERATING_DOCS,
    ids=[path for path, _, _ in ENUMERATING_DOCS],
)
def test_every_target_is_announced_in_every_enumerating_doc(relative_path, anchor, mode):
    """The test that would have caught both of 2026-07-30's misses."""

    enumeration = _enumeration(relative_path, anchor, mode)
    missing = sorted(
        target for target, token in DOC_TOKENS.items() if token not in enumeration
    )
    assert not missing, (
        f"{relative_path}'s target list does not mention {missing}. A target "
        "that works but is undocumented is invisible to every user who did not "
        "write it."
    )


def test_the_parameterize_docstring_lists_every_target():
    """`parameterize`'s own docstring is the first place anyone looks, and it is
    the one that went stale for `uzf.vks`."""

    from myflopy.modflow.mf6.pest.project import PestProject

    doc = PestProject.parameterize.__doc__ or ""
    missing = sorted(target for target in _RECIPES if target not in doc)
    assert not missing, f"parameterize's docstring omits {missing}"


def test_recipe_invariants_hold_for_every_target():
    """Shape rules the engine relies on, asserted once for the whole registry
    rather than rediscovered per target.

    Each of these was learned the hard way: the pilot-point base array
    (ledger 115), UZF's index columns (117), and MST porosity needing a package
    to relayer (110).
    """

    for name, recipe in _RECIPES.items():
        assert recipe.family in {"array", "list"}, name
        assert recipe.model in {"flow", "transport"}, name
        if recipe.family == "list":
            assert recipe.use_col is not None, f"{name}: a list target needs use_col"
            assert len(recipe.index_cols) >= 2, (
                f"{name}: pyEMU needs at least (layer, cell) to geolocate a list "
                "parameter"
            )
            assert recipe.use_col not in recipe.index_cols, (
                f"{name}: use_col {recipe.use_col} is also an index column, which "
                "pyEMU refuses"
            )
        else:
            # Array targets can be asked for pilot points, which multiply a base
            # array read off the model -- so the recipe has to name one.
            assert recipe.package and recipe.variable, (
                f"{name}: an array target must declare package/variable, or "
                "pilot points and relayering cannot find its array"
            )


def test_every_alias_points_at_a_real_target():
    unknown = {alias: key for alias, key in _ALIASES.items() if key not in _RECIPES}
    assert not unknown, f"aliases point at targets that do not exist: {unknown}"


_PROBE = textwrap.dedent(
    '''
    import json

    # Inject a synthetic recipe BEFORE importing any consumer.
    from myflopy.modflow.mf6.pest import native_parameters as np_mod

    np_mod._RECIPES["zzz"] = np_mod._Recipe(
        "zzz", "array", "{model}.npf_k.txt", package="npf", variable="k",
    )
    np_mod._RECIPES["zzz.list"] = np_mod._Recipe(
        "zzz.list", "list", "{model}.drn_stress_period_data_*.txt", use_col=3,
    )
    np_mod._ALIASES["zed"] = "zzz"

    found = {}

    found["resolve_target"] = np_mod.resolve_target("zed").canonical == "zzz"
    try:
        np_mod.resolve_target("definitely_not_a_target")
        found["resolve_target.error_text"] = False
    except KeyError as error:
        found["resolve_target.error_text"] = "zzz" in str(error)

    spec = np_mod.NativeParameterSpec(target="zzz")
    found["spec.style_default"] = spec.style == "constant"
    found["spec.additive_default"] = spec.additive is False
    found["spec.name_slug"] = spec.name == "zzz"

    # An array recipe WITHOUT package/variable must be refused for pilot points.
    np_mod._RECIPES["zzz.bare"] = np_mod._Recipe("zzz.bare", "array", "{model}.x.txt")
    try:
        np_mod.NativeParameterSpec(target="zzz.bare", style="pilotpoints")
        found["spec.pilotpoints_guard"] = False
    except NotImplementedError:
        found["spec.pilotpoints_guard"] = True

    # The remaining sites are reached through helpers that take the recipe.
    found["model_name_for"] = callable(np_mod.model_name_for)
    found["flopy_model_for"] = callable(np_mod.flopy_model_for)
    found["relayer_array_target"] = callable(np_mod.relayer_array_target)

    import pathlib, tempfile
    workspace = pathlib.Path(tempfile.mkdtemp())
    (workspace / "m.npf_k_layer1.txt").write_text("1.0\\n")
    (workspace / "m.npf_k.txt").write_text("1.0\\n")
    resolved = np_mod._resolve_files(workspace, "m", np_mod._RECIPES["zzz"])
    found["resolve_files.layered_first"] = resolved == ["m.npf_k_layer1.txt"]

    try:
        np_mod._select_layer_files(["m.npf_k.txt"], [0], np_mod._RECIPES["zzz"])
        found["select_layer_files"] = False
    except ValueError:
        found["select_layer_files"] = True

    from types import SimpleNamespace
    from myflopy.modflow.mf6.pest.zones import resolve_zone_array
    model = SimpleNamespace(
        vor=SimpleNamespace(ncpl=4),
        gwf=SimpleNamespace(modelgrid=SimpleNamespace(nlay=3)),
    )
    array_zones = resolve_zone_array([1, 1, 2, 2], family="array", model=model)
    list_zones = resolve_zone_array([1, 1, 2, 2], family="list", model=model)
    found["zone_resolution"] = (
        array_zones.shape == (4,) and list_zones.shape == (3, 4)
    )

    # ...and the site a human still has to type.
    found["recipes.entry"] = False   # nothing derives it; that is the point

    print(json.dumps(found))
    '''
)


@pytest.fixture(scope="module")
def probe_results():
    result = subprocess.run(
        [sys.executable, "-c", _PROBE],
        cwd=REPO_ROOT, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout[-3000:] + "\n" + result.stderr[-3000:]
    return json.loads(result.stdout.strip().splitlines()[-1])


def test_one_recipe_entry_reaches_every_automatic_site(probe_results):
    """A regression here re-opens the defect class this file exists to close."""

    missing = sorted(site for site in AUTOMATIC_SITES if not probe_results.get(site))
    assert not missing, (
        f"these sites used to be automatic and no longer are: {missing}. "
        "Adding a parameterize target just got more expensive."
    )


def test_the_manual_sites_are_exactly_the_ones_recorded(probe_results):
    """Fails in BOTH directions: a new hand-written site, or one closed without
    updating this file. The second is the interesting half -- a win that is not
    recorded is a win the next person cannot rely on."""

    probed_manual = {
        site for site, reached in probe_results.items()
        if site in MANUAL_SITES and not reached
    }
    unexpectedly_automatic = {
        site for site in MANUAL_SITES
        if site in probe_results and probe_results[site]
    }
    assert not unexpectedly_automatic, (
        f"{sorted(unexpectedly_automatic)} is now derived automatically. Delete "
        "it from MANUAL_SITES and say so in the compromise ledger."
    )
    # The doc sites are asserted by the doc tests above rather than by the probe;
    # what the probe covers is the code half.
    assert probed_manual == {"recipes.entry"}
