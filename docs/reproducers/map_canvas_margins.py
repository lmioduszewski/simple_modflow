"""A plotly map is drawn into a SMALLER canvas than its declared margins say.

Ledger 311 (amending ledger 307). ``viz.map_canvas`` fitted a map's zoom to the
figure's size less its DECLARED margins -- ``r=0`` on a standalone map, the
template's ``l=50``/``r=0`` on an animation. At draw time plotly.js GROWS those
margins to fit what it adds: the colorbar (the right margin measured 81 px) and an
animation's slider and play buttons (80 px each side). The fit was therefore for a
canvas wider than the one drawn, and a wide window lost both edges -- measured on the
canonical model with ``extent=(0, 900, 2100, 1200)``: a standalone map at 700x450
showed 97.3% of the extent's width, an animation 91.4%.

``map_canvas(..., colorbar=True, controls=True)`` now reserves
``COLORBAR_MARGIN`` (100) and ``CONTROLS_MARGIN`` (80). After the fix, measured in a
browser, every picture CONTAINS its window: standalone 700x450 -> 1.135x the extent's
width visible (canvas 619 px), 1100x600 -> 1.854x, animation -> 1.142x (canvas 540
px), each mosaic panel -> 1.113x, the default grid map -> 1.100x on its tighter axis.

The check needs a real browser, so this writes the HTML and the JavaScript to run in
it rather than asserting anything in CI::

    PYTHONPATH=src python docs/reproducers/map_canvas_margins.py <out-dir>
    python -m http.server --directory <out-dir> 8791
    # open http://127.0.0.1:8791/single_700.html, run MEASURE_JS in the console

``MEASURE_JS`` reads MapLibre's own visible bounds -- the only ground truth; a test
that takes ``map_canvas`` as the truth cannot see this defect by construction, which
is how ledger 307's tests passed over it.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

MEASURE_JS = """
const boxes = await fetch('boxes.json').then(r => r.json());
const gd = document.querySelector('.js-plotly-plot');
const key = location.pathname.includes('grid') ? 'grid' : 'extent';
Object.keys(gd._fullLayout).filter(k => /^map\\d*$/.test(k)).map(sid => {
  const sub = gd._fullLayout[sid]._subplot.map, b = sub.getBounds();
  const [w, s, e, n] = boxes[key];
  return {sid, canvas: [sub.getContainer().clientWidth, sub.getContainer().clientHeight],
          widthShare: (b.getEast() - b.getWest()) / (e - w),
          heightShare: (b.getNorth() - b.getSouth()) / (n - s),
          contains: b.getWest() <= w && b.getEast() >= e && b.getSouth() <= s && b.getNorth() >= n};
});
"""


def main(out: Path) -> None:
    from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, build_canonical_model

    out.mkdir(parents=True, exist_ok=True)
    workspace = out / "model"
    shutil.rmtree(workspace, ignore_errors=True)
    model = build_canonical_model(workspace, config=CanonicalModelConfig.testing())
    assert model.run_simulation()[0]

    wide = (0.0, 900.0, 2100.0, 1200.0)
    single = model.plot.map(layer=0, extent=wide)
    html = {"include_plotlyjs": True, "full_html": True}
    for width, height in ((700, 450), (1100, 600)):
        single.fig.write_html(
            out / f"single_{width}.html", default_width=f"{width}px",
            default_height=f"{height}px", **html,
        )
    model.hds.animate(extent=wide).write_html(
        out / "anim_700.html", default_width="700px", default_height="450px", **html
    )
    model.hds.mosaic(by="layer", extent=wide).write_html(
        out / "mosaic.html", default_width="700px", **html
    )
    grid = model.plot.map(layer=0)
    grid.fig.write_html(
        out / "grid_700.html", default_width="700px", default_height="450px", **html
    )
    json.dump(
        {"extent": list(single.extent.wgs84), "grid": list(grid.latlon_bounds)},
        (out / "boxes.json").open("w"),
    )
    (out / "measure.js").write_text(MEASURE_JS)
    print(f"wrote {out}; serve it and run measure.js in each page's console")


if __name__ == "__main__":
    main(Path(sys.argv[1] if len(sys.argv) > 1 else "map_canvas_margins"))
