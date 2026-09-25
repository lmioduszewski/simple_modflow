"""A typed ``ievt`` beside an omitted ET ``surface``, and a numpy-integer index (ledger 324).

    /home/lukem/python/envs/gw/bin/python docs/reproducers/lazy_typed_ievt.py

No MF6 needed; writes the input files only. About 3 s.

PART 1 -- ``mf.evt.array(rate=, depth=, ievt=1)`` with no ``context=`` and ``surface``
omitted goes LAZY for the surface (ledger 285). Before ledger 324 the lazy branch built
its options from the value fields and the other keywords, and a typed ``ievt`` is
neither: the written EVT file had NO ``IEVT`` block, which MF6 reads as layer 1 in every
column. Printed: whether the spec carries ``ievt``, what FloPy holds, and the lines of
the written file naming IEVT. On a tree before the fix the last two are ``None`` and
nothing.

PART 2 -- ``mf.rch.array(recharge=, irch=np.int64(1))``, the EAGER form: FloPy refused a
numpy integer scalar ("object of type 'numpy.int64' has no len()", an MFDataException).
Printed: whether it builds.

PART 3 -- how many stored manifests hold a lazy array spec, the declarations whose
digest the carry could move: every ``*.json`` under the repo and the model tree naming
``_build_lazy_array``. Read only.
"""

from __future__ import annotations

import pathlib
import subprocess
import tempfile
import warnings

import numpy as np

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid

REPO = pathlib.Path(__file__).resolve().parents[2]


def _simulation(grid, package):
    stack = mf.LayerStack(top=mf.Flat(100.0)).add("fill", thickness=5.0).add("sand", thickness=30.0)
    model = mf.gwf("m", context=mf.ModelContext(grid=grid), packages=[
        mf.disv(stack), mf.npf(k=10.0), mf.ic(strt=99.0), package,
    ])
    return mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])


def main() -> None:
    warnings.simplefilter("ignore")
    grid = irregular_voronoi_grid(nrow=4, ncol=4, cell_size=50.0)
    with tempfile.TemporaryDirectory() as tmp:
        workspace = pathlib.Path(tmp) / "typed"
        spec = mf.evt.array(rate=1e-3, depth=5.0, ievt=1)
        built = _simulation(grid, spec).build_flopy(workspace)
        held = built.built_model("m").model.get_package("evt").ievt.get_data()
        built.simulation.write_simulation(silent=True)
        written = []
        for path in workspace.glob("*.evt*"):
            lines = path.read_text().splitlines()
            for i, line in enumerate(lines):
                if line.strip().lower() == "ievt":     # the block's header, then its data
                    written.append(f"{line.strip()} {' '.join(lines[i + 1].split()[:4])}")
        print(f"PART 1  spec carries ievt: {'ievt' in spec.options}; inputs {dict(spec.inputs)}")
        print(f"        FloPy holds ievt: {None if held is None else np.unique(held[0]).tolist()}")
        print(f"        written IEVT (1-based in the file): {written}")

        try:
            _simulation(grid, mf.rch.array(recharge=1e-3, irch=np.int64(1))).build_flopy(
                pathlib.Path(tmp) / "numpy"
            )
            print("PART 2  irch=np.int64(1): builds")
        except Exception as error:  # noqa: BLE001 - a reproducer reports whatever FloPy raises
            print(f"PART 2  irch=np.int64(1): {type(error).__name__}")

    roots = [REPO, pathlib.Path("/home/lukem/models")]
    found = []
    for root in roots:
        if root.exists():
            result = subprocess.run(
                ["grep", "-rl", "--include=*.json", "_build_lazy_array", str(root)],
                capture_output=True, text=True, check=False,
            )
            found += [line for line in result.stdout.splitlines() if "/.claude/worktrees/" not in line]
    print(f"PART 3  stored JSON naming a lazy array spec: {len(found)}")


if __name__ == "__main__":
    main()
