"""PstFrom's ``chunk_len`` decides whether every realization spawns a process pool (ledger 511).

A forward run's first step is pyEMU's ``apply_list_and_array_pars``, which multiplies
every parameterized model file. ``PstFrom`` writes its ``chunk_len`` into that call, and
pyEMU's default is 50: at 50 files or more it splits them into chunks and applies the
chunks in a ``multiprocessing`` ``spawn`` pool of ``min(cpu_count, 60)`` processes, each
of which re-imports the forward run's modules. That happens in EVERY realization, beside
every other IES worker's. myflopy builds with ``chunk_len=10**6``
(``PYEMU_APPLY_CHUNK_LEN``), so one process applies every file.

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/pstfrom_chunk_len_pool.py

Needs ``mf6`` and pyEMU; writes only under a temporary directory, which it removes.

THE SETUP. A 6 x 6 irregular Voronoi grid of 100 ft cells, one layer from 110 to 0 ft,
K 10 ft/d, a CHD on the west and east edge cells declared in each of 61 one-day steady
periods, so FloPy writes 61 CHD list files; the calibration is ``chd`` constant (61
multiplier files) and one head target. The template's apply step is then run as its own
Python process -- ``apply_list_and_array_pars(arr_par_file="mult2model_info.csv",
chunk_len=N)`` in the template, the call ``forward_run.py`` makes -- once with pyEMU's
50 and once with the library's 10**6, sampling the process tree every 20 ms: the most
child processes alive at once, and the most proportional set size (PSS) the tree held.
Wall time is the apply process's. Each setting is run three times; the table shows the
middle of the three.

MEASURED (pyEMU 1.4.0, FloPy 3.10.0, 20 cores, 2026-10-12, load average ~1.5-3)
  chunk_len   child processes (peak)   tree PSS (peak)   wall (3 runs)
  50          21                       732 MB            2.92-2.97 s
  10**6       0                        201 MB            1.71-1.79 s
The 21 are 20 pool workers (``min(cpu_count, 60)``) and the resource tracker. The PSS
is sampled every 20 ms, so it is a floor, and grows with the files: on Ten Trails'
cells it was ~3.4 GB (below).

The scale finding this records (scale-4, on a 61-period copy of Ten Trails' steady
model with GHB/RIV/DRN repeated each period): 21 child processes at peak (20 pool
workers and the resource tracker), ~3.4 GB summed PSS and 4.4 s at ``chunk_len=50``,
against 0 children and 1.3 s at ``10**6``. The threshold is 50 files, not 51: pyEMU
builds ``main_chunks + [remainder]``, so exactly 50 files give two chunks, the second
empty, and still open the pool.
"""

from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _absolute_pythonpath import make_pythonpath_absolute  # noqa: E402

make_pythonpath_absolute()

import statistics  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402
import warnings  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import psutil  # noqa: E402

import myflopy as mf  # noqa: E402
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid  # noqa: E402
from myflopy.modflow.mf6.pest.project import PYEMU_APPLY_CHUNK_LEN  # noqa: E402

NPER = 61
APPLY = (
    "import pyemu\n"
    "pyemu.helpers.apply_list_and_array_pars(arr_par_file='mult2model_info.csv', "
    "chunk_len={chunk_len})\n"
)


def build(root: pathlib.Path) -> pathlib.Path:
    """The model in the docstring, run and built; its template directory."""

    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    chd = {per: [[(0, int(c)), 102.0 - 0.01 * per] for c in edge] for per in range(NPER)}
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=110.0, botm=[0.0]),
        mf.npf(k=10.0, icelltype=0), mf.ic(strt=100.0),
        mf.chd(stress_period_data=chd),
        mf.oc(head_filerecord="m.hds", saverecord=[("HEAD", "LAST")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=NPER, perioddata=[(1.0, 1, 1.0)] * NPER)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = mf.Project(root / "p", name="p").run("r", simulation)
        cal = run.model().pest("chunks", workspace=root / "chunks", start_datetime="2024-01-01")
        cal.parameterize("chd", style="constant", bounds=(0.9, 1.1))
        cal.observe(mf.HeadTargets(
            locations=pd.DataFrame({"name": ["W1"], "cell": [14], "layer": [0]}),
            values=pd.DataFrame({"time": [NPER - 1], "W1": [100.0]}),
        ))
        cal.build("chunks.pst", noptmax=0)
    return cal.template_workspace


def measure(template: pathlib.Path, chunk_len: int) -> tuple[int, float, float]:
    """(peak child processes, peak tree PSS in MB, wall seconds) of one apply step."""

    start = time.perf_counter()
    process = subprocess.Popen(
        [sys.executable, "-c", APPLY.format(chunk_len=chunk_len)], cwd=template,
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    parent = psutil.Process(process.pid)
    children = pss = 0.0
    while process.poll() is None:
        try:
            tree = [parent, *parent.children(recursive=True)]
            children = max(children, len(tree) - 1)
            total = 0
            for member in tree:
                try:
                    total += member.memory_full_info().pss
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    continue
            pss = max(pss, total / 2**20)
        except psutil.NoSuchProcess:
            break
        time.sleep(0.02)
    wall = time.perf_counter() - start
    if process.returncode:
        raise RuntimeError(f"the apply step failed at chunk_len={chunk_len}")
    return int(children), pss, wall


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="chunk_len_") as directory:
        template = build(pathlib.Path(directory))
        files = len(pd.read_csv(template / "mult2model_info.csv"))
        written = "chunk_len=" + (template / "forward_run.py").read_text().split(
            "chunk_len=")[1].split(")")[0]
        print(f"{files} multiplier files; this tree's forward_run.py applies them with "
              f"{written}; cpu_count {psutil.cpu_count()}")
        print(f"{'chunk_len':>10} {'children':>9} {'PSS MB':>8} {'wall s':>7}")
        for chunk_len in (50, PYEMU_APPLY_CHUNK_LEN):
            runs = [measure(template, chunk_len) for _ in range(3)]
            middle = sorted(runs, key=lambda run: run[2])[1]
            print(f"{chunk_len:>10} {middle[0]:>9} {middle[1]:>8.0f} {middle[2]:>7.2f}"
                  f"   (walls {', '.join(f'{run[2]:.2f}' for run in runs)};"
                  f" median {statistics.median(run[2] for run in runs):.2f})")


if __name__ == "__main__":
    main()
