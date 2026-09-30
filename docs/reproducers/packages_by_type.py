"""Several packages of one type: found by TYPE, and never one read as the model's in silence.

Ledger 352. FloPy's ``get_package(name)`` answers a package NAMED like its argument
before it looks at types; myflopy's ``model.fl.<pkg>``, ``model.package(<type>)`` and
every LAK/SFR/UZF view asked it by a bare type string. And the budget-file reader
under every cell-budget view kept the FIRST record of each saved time.

    PYTHONPATH=<abs path to a tree's src> /home/lukem/python/envs/gw/bin/python \\
        docs/reproducers/packages_by_type.py

It measures whichever tree is on ``PYTHONPATH`` (give it ABSOLUTE: ledger 344). The
"before" column below is the parent commit's ``src`` (89364c2, the main checkout on
2026-09-30); "after" is this tree. Needs ``mf6``. The models are
``tests/test_packages_by_type.py``'s (imported from it, so the two cannot drift): the
6 x 6 irregular Voronoi grid of ``tests/test_obs_route_declared_identity.py`` (ncpl 36,
one layer, CHD on the west and east edges, three 10-day periods of one step) carrying
two packages each of LAK (cells 8-9, 25-26), SFR (13-16, 31-32), RIV (20; 21-22) and
UZF (1-2, 3-4), the first of each named like its type (``lak``, ``sfr``, ``riv``,
``uzf``) or none so (``lak_a`` ...).

MEASURED (before -> after):

A  ``lak`` + ``lak_2``: ``model.fl.lak`` / ``model.package("lak")`` -> ``lak`` ALONE ->
   ``['lak', 'lak_2']``; ``model.fl.riv`` -> ``riv`` alone -> ``['riv', 'crisp']``.
   ``model.lak.results.stage.get()`` -> ``lak``'s lake (cells 8, 9) as the model's,
   ``lak_2``'s absent, no word -> ``SeveralPackagesError`` naming both and
   ``model.lak(package='lak_2')``; ``model.lak(package="lak_2")`` reads cells 25, 26.
B  ``lak_a`` + ``lak_b``: the same views -> ``AttributeError: 'list' object has no
   attribute 'nlakes'`` (SFR: ``'nreaches'``) -> ``SeveralPackagesError``.
C  ``riv`` + ``crisp``: ``model.riv.results.q.get()`` period 0 -> ``riv``'s cell 20 alone
   (the budget file holds cells 20, 21, 22) -> all three, ``pname`` RIV / CRISP;
   ``model.budget.riv`` likewise.
D  UZF declared BEFORE RCH: MF6 still wrote ``RCH`` before ``UZF-GWRCH`` in the budget
   file, so FloPy's first-substring match of "RCH" found RCH on both trees -- the
   hazard the old ``paknam=`` comment named was NOT reached here.
E  ``get_data(text="RIV", paknam="riv")`` raises ``ValueError``: MF6 writes the MODEL's
   name in ``paknam`` (``M``) and the package's in ``paknam2``, so the old reader's
   ``paknam=`` never matched and always fell back to the text alone.
F  Ten Trails (``/home/lukem/models/mf6/tentrails/runs/hm_mounding_new_lak``, read-only,
   skipped when absent): 38 periods of RIV ``riv`` (500 cells) + ``crisp`` (521 cells).
   ``model.riv.results.q.get()`` period 0 -> 500 rows, -208,680.73 ft3/d -> 1,021 rows,
   -448,427.27 (``crisp`` -239,746.54: 53.46%); all 38,798 entries equal the budget file.
"""

from __future__ import annotations

import pathlib
import sys
import tempfile
import warnings

import flopy
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
warnings.simplefilter("ignore")

import myflopy as mf  # noqa: E402 -- after the path
from tests.test_obs_route_declared_identity import _run  # noqa: E402
from tests.test_packages_by_type import APART, LIKE_THE_TYPE, _two_of_each  # noqa: E402

print("myflopy from", mf.__file__)


def attempt(label, read):
    try:
        out = read()
    except Exception as error:  # noqa: BLE001 -- a reproducer reports whatever it hit
        print(f"  {label}: {type(error).__name__}: {str(error)[:160]}")
        return None
    print(f"  {label}: {out}")
    return out


def names_of(found):
    return [p.package_name for p in found] if isinstance(found, list) else [found.package_name]


def cells_of(frame):
    return sorted(set(int(c) for c in frame["cell"]))


with tempfile.TemporaryDirectory() as tmp:
    root = pathlib.Path(tmp)
    like = _two_of_each(root / "like", LIKE_THE_TYPE)
    apart = _two_of_each(root / "apart", APART)

    print("\nA  lak + lak_2 (named like the type)")
    attempt("model.fl.lak", lambda: names_of(like.fl.lak))
    attempt('model.package("lak")', lambda: names_of(like.package("lak")))
    attempt("model.fl.riv", lambda: names_of(like.fl.riv))
    attempt("lak.results.stage.get() cells", lambda: cells_of(like.lak.results.stage.get()))
    attempt('lak(package="lak_2") cells', lambda: cells_of(like.lak(package="lak_2").results.stage.get()))

    print("\nB  lak_a + lak_b (none named like the type)")
    attempt("lak.results.stage.get()", lambda: cells_of(apart.lak.results.stage.get()))
    attempt("sfr.results.stage.get()", lambda: cells_of(apart.sfr.results.stage.get()))

    print("\nC  riv + crisp: the budget file against model.riv.results.q")
    cbc = flopy.utils.CellBudgetFile(f"{like.workspace}/m.cbc", precision="double")
    for record_index in np.flatnonzero(
        (np.char.strip(cbc.recordarray["text"].astype(str)) == "RIV") & (cbc.recordarray["kper"] == 1)
    ):
        record = cbc.get_record(int(record_index))
        print(f"  file: {cbc.recordarray['paknam2'][record_index].decode().strip():<6}"
              f" cells {list(np.asarray(record['node']) - 1)}")
    attempt("riv.results.q.get(per=0) cells", lambda: cells_of(like.riv.results.q.get(per=0)))
    attempt("budget.riv.get() cells", lambda: cells_of(like.budget.riv.get()))

    print("\nE  paknam= against MF6's budget file")
    attempt('get_data(text="RIV", paknam="riv")',
            lambda: len(cbc.get_data(text="RIV", paknam="riv", kstpkper=(0, 0))))
    print(f"  paknam (from) names: {[n.decode().strip() for n in cbc.get_unique_package_names()]}")

    print("\nD  UZF declared BEFORE RCH")
    uzf = mf.uzf.flopy(
        packagedata=[[i, (0, c), 1, -1, 0.1, 0.1, 0.05, 0.30, 0.15, 4.0] for i, c in enumerate([1, 2])],
        perioddata={0: [[i, 0.05, 0.0, 1.0, 0.1, 0.0, 0.0, 0.0] for i in range(2)]}, save_flows=True)
    rch = mf.rch.flopy(stress_period_data={0: [[(0, 3), 0.002], [(0, 4), 0.003]]}, save_flows=True)
    first, _ = _run(root / "uzf_first", [uzf, rch])
    order = flopy.utils.CellBudgetFile(f"{first.workspace}/m.cbc", precision="double").textlist
    print(f"  budget record order: {[t.decode().strip() for t in order]}")
    attempt("rch.results.q.get(per=0) cells", lambda: cells_of(first.rch.results.q.get(per=0)))

TEN_TRAILS = pathlib.Path("/home/lukem/models/mf6/tentrails/runs/hm_mounding_new_lak")
print("\nF  Ten Trails (read-only)")
if not (TEN_TRAILS / "tentrails.cbb").exists():
    print(f"  skipped: {TEN_TRAILS} is not here")
else:
    model = mf.load_mf6_run(str(TEN_TRAILS), verbosity_level=0)
    attempt("model.fl.riv", lambda: names_of(model.fl.riv))
    frame = model.riv.results.q.get(per=0)
    by = "pname" if "pname" in frame.columns else "package"
    print(f"  riv.results.q.get(per=0): {len(frame)} rows, "
          f"{frame.groupby(by)['q_gwf'].sum().round(2).to_dict()}")
