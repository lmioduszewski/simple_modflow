"""`plot.map(values=)` of the wrong shape drew a garbled map, without a word.

Found 2026-10-15 by the recalibration-helpers audit (defect NEW1): the obvious way to map
a calibrated K, `model.plot.map(values=ies.calibrated().get("npf", "k"))`, hands the map
the WHOLE `(nlay, ncpl)` array.

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/map_values_shape.py

THE MECHANISM. Both verbs that take `values=` -- the free `plot.map` (and `model.plot.map`,
which calls it) and `vor.plot.map`, which builds its picture itself -- did
`custom_zs = list(values)` and nothing else, and `Choro.__init__` stored that list past
the one place that checked a length (the `custom_zs` setter's `assert`, which nothing
calls). `list()` of a `(nlay, ncpl)` array is nlay rows, of a dict its KEYS, of a
DataFrame its column NAMES. Plotly pairs `z` with `locations` by position, so it drew
whatever it was given; Matplotlib's branch assigned the list to a GeoDataFrame column
and raised pandas' "Length of values (4) does not match length of index (441)", a
message naming neither `values=` nor a layer.

THE SETUP. `CanonicalModelConfig.testing()` (21x21 Voronoi, 441 cells, 4 layers), built
and NOT run -- `values=` needs no results. `k` is the model's own `npf.k`, the shape
`calibrated().get("npf", "k")` returns. Each shape is passed to each scope at
`layer=1`; a picture that builds is reported by the length of its cell trace's `z`
against its 441 `locations` (plotly) or as a `Figure` (mpl).

WHAT IT FOUND at 55947911: every wrong shape was ACCEPTED on plotly at every scope --
`(4, 441)` as z of length 4, a flattened `(1764,)` as z of length 1764 whose first 441
are LAYER 0 under `layer=1`, a `(1, 441)` row as z of length 1, a dict as its keys, a
wrong length as itself -- and raised pandas' error on mpl; a single number raised
"'float' object is not iterable". A masked array's masked cells were painted 0.0. After
the fix (ledger 601) every wrong shape raises at the call on both
backends, naming the expected `(441,)` and `values=array[1], layer=1`; a `(1, 441)` row
draws as that row and a masked cell draws as a gap.
"""

from __future__ import annotations

import tempfile
import warnings
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
warnings.simplefilter("ignore")

from myflopy import plot  # noqa: E402
from myflopy.modflow.mf6.canonical_example import (  # noqa: E402
    CanonicalModelConfig,
    build_canonical_model,
)


def describe(build):
    try:
        picture = build()
    except Exception as exc:  # noqa: BLE001 - a reproducer reports whatever was raised
        return f"raises {type(exc).__name__}: {str(exc).splitlines()[0][:150]}"
    if hasattr(picture, "fig"):
        trace = picture.fig.data[0]
        return f"accepted: z {len(trace.z)} vs {len(trace.locations)} locations"
    return f"accepted: {type(picture).__name__}"


with tempfile.TemporaryDirectory() as tmp:
    model = build_canonical_model(Path(tmp) / "canonical", config=CanonicalModelConfig.testing())
    k = np.asarray(model.gwf.npf.k.array, dtype=float)
    ncpl = model.vor.ncpl
    print(f"grid: {ncpl} cells, k shape {k.shape}\n")

    shapes = {
        "(nlay, ncpl)": k,
        "(nlay*ncpl,)": k.reshape(-1),
        "(1, ncpl)": k[1:2],
        "(ncpl,)": k[1],
        "100 values": k[1][:100],
        "ncpl + 1": np.r_[k[1], 1.0],
        "dict {per: array}": {0: k[0], 1: k[1]},
        "DataFrame": pd.DataFrame({"k": k[1]}),
        "3.0": 3.0,
    }
    scopes = {
        "plot.map(model)": lambda v, b: plot.map(model, v, layer=1, backend=b),
        "model.plot.map": lambda v, b: model.plot.map(values=v, layer=1, backend=b),
        "vor.plot.map": lambda v, b: model.vor.plot.map(values=v, backend=b),
    }
    for label, values in shapes.items():
        print(label)
        for scope, call in scopes.items():
            for backend in ("plotly", "mpl"):
                print(f"  {scope:16s} {backend:6s} {describe(lambda: call(values, backend))}")

    print("\nflattened at layer=1 -- which layer did plotly paint?")
    try:
        z = np.asarray(model.plot.map(values=k.reshape(-1), layer=1).fig.data[0].z, float)
        print(f"  first {ncpl} z equal layer 0: {np.array_equal(z[:ncpl], k[0])}; "
              f"layer 1: {np.array_equal(z[:ncpl], k[1])}")
    except ValueError as exc:
        print(f"  raises: {str(exc).splitlines()[0][:150]}")

    print("\na masked array (every 7th cell masked)")
    mask = np.arange(ncpl) % 7 == 0
    z = np.asarray(
        model.plot.map(values=np.ma.masked_array(k[1], mask=mask), layer=1).fig.data[0].z,
        dtype=float,
    )
    print(f"  {mask.sum()} masked cells: NaN (a gap) {np.isnan(z[mask]).sum()}, "
          f"0.0 {(z[mask] == 0.0).sum()}")
