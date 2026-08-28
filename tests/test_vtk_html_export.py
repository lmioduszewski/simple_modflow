"""A written 3-D scene must open from disk, not only over HTTP.

PyVista emits its vtk.js bundle as ``<script type="module">``. A module script is
not executed when the page is opened from ``file://``, and the bundle's last act
is ``window.OfflineLocalView = {...}``, which a following CLASSIC script calls.
So the viewer shell runs, the loader never does, and vtk.js shows its "Drop File
/ Explore Scene" placeholder with

    Uncaught ReferenceError: OfflineLocalView is not defined

in the console. Served over HTTP the identical file is fine, which is what makes
it look like a browser bug rather than a packaging one -- it was reported as
"doesn't display correctly in Firefox" (2026-08-27).
"""

from __future__ import annotations

import re

import numpy as np
import pytest

import myflopy as mf


@pytest.fixture(scope="module")
def scene(tmp_path_factory):
    pytest.importorskip("pyvista")
    ws = tmp_path_factory.mktemp("vtk_html")
    tri = mf.TriangleGrid(model_ws=str(ws), angle=30)
    tri.set_domain_rectangle(x_dist=400, y_dist=300, origin=(0, 0), max_area=2000)
    tri.build()
    vor = mf.VoronoiGridPlus(tri, crs="2927")
    stack = mf.LayerStack(
        vor, top=mf.Surface.from_array(np.full(vor.ncpl, 100.0))
    ).add("a", thickness=20)
    return stack.plot.surface("all", backend="vtk")


@pytest.mark.slow
def test_the_written_page_has_no_module_script(scene, tmp_path):
    """The one marker that makes a page unopenable from disk."""

    path = scene.html(tmp_path / "scene.html")
    assert 'type="module"' not in path.read_text(encoding="utf-8")


@pytest.mark.slow
def test_the_loader_is_still_defined_before_it_is_used(scene, tmp_path):
    """Stripping the marker must not reorder anything.

    A classic script runs at parse time -- BEFORE the consumer below it -- where
    the module was merely deferred and happened to win a ``setTimeout(..., 0)``
    race. So this ordering is more robust after the change, not less.
    """

    html = scene.html(tmp_path / "scene.html").read_text(encoding="utf-8")
    assign = html.find("window.OfflineLocalView")
    use = html.find("OfflineLocalView.load")
    assert assign != -1 and use != -1, "the export no longer looks like vtk.js"
    assert assign < use


@pytest.mark.slow
def test_the_bundle_really_is_classic_script_safe(scene, tmp_path):
    """Guards the assumption the fix rests on.

    Dropping ``type="module"`` is only valid because this bundle uses no module
    syntax. If PyVista ever ships one that does, this fails rather than writing
    a page that silently throws.
    """

    html = scene.html(tmp_path / "scene.html").read_text(encoding="utf-8")
    bodies = [
        html[m.end():html.find("</script>", m.end())]
        for m in re.finditer(r"<script([^>]*)>", html)
    ]
    joined = "\n".join(bodies)
    for pattern, what in (
        (r"^\s*import\s+", "top-level import"),
        (r"^\s*export\s+", "top-level export"),
        (r"\bimport\.meta\b", "import.meta"),
        (r"\bimport\s*\(", "dynamic import()"),
    ):
        assert not re.search(pattern, joined, re.M), f"bundle now uses {what}"

    # Top-level await also requires a module; count awaits at brace depth 0.
    depth = toplevel = 0
    for token in re.finditer(r"[{}]|\bawait\b", joined):
        text = token.group(0)
        if text == "{":
            depth += 1
        elif text == "}":
            depth -= 1
        elif depth == 0:
            toplevel += 1
    assert toplevel == 0, "bundle now uses top-level await; it needs type=module"


@pytest.mark.slow
def test_the_page_is_still_self_contained(scene, tmp_path):
    """The fix must not have traded offline use for disk use."""

    html = scene.html(tmp_path / "scene.html").read_text(encoding="utf-8")
    external = re.findall(r'(?:src|href)\s*=\s*["\'](https?://[^"\']+)', html)
    assert not external, f"page now fetches from the network: {external[:3]}"
