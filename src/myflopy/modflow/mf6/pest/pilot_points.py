"""Native pilot-point K parameterization (inverse-distance weighting).

pyEMU's built-in pilot points fall back to pure-Python kriging on unstructured
(DISV/Voronoi) grids -- minutes per build -- and its fast ``pypestutils`` path
either does not engage for unstructured grids or explodes into a non-stationary
hyperparameter setup. So pilot points on Voronoi are done the way the (now
retired) legacy path did: place points on the grid and **interpolate to cells
with inverse-distance weighting at forward-run time**. This module is the
build-time half (placement + support files + parameter registration); the
forward-run half is ``forward_run.apply_pilotpoints_to_array``.

Exposed through the one unified API::

    cal.parameterize("k", style="pilotpoints", pp_space=6, correlation=...,
                     layers=[0, 1], bounds=(0.05, 20), physical=(1e-3, 300))
    cal.parameterize("k", style="pilotpoints", pp_points=my_gdf)   # explicit
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from myflopy.modflow.mf6.pest.native_parameters import _flatten_array_file, _resolve_files


def _cell_centers(model) -> tuple[np.ndarray, np.ndarray]:
    """The model grid's cell-center ``(x, y)`` coordinate arrays, flattened."""

    grid = model.gwf.modelgrid
    return (
        np.asarray(grid.xcellcenters, dtype=float).reshape(-1),
        np.asarray(grid.ycellcenters, dtype=float).reshape(-1),
    )


def place_pilot_points(model, *, pp_space=None, pp_points=None) -> pd.DataFrame:
    """Return a pilot-point table (``ppname, x, y``) for the model grid.

    ``pp_points`` (a GeoDataFrame, an array, or a list of ``(x, y)``) places
    points explicitly. Otherwise ``pp_space`` lays a regular net at a spacing of
    ``pp_space`` median-cell-widths and keeps the cell nearest each net node.
    """

    xc, yc = _cell_centers(model)
    if pp_points is not None:
        if hasattr(pp_points, "geometry"):
            pts = np.array([(geom.x, geom.y) for geom in pp_points.geometry], dtype=float)
        else:
            pts = np.asarray(pp_points, dtype=float).reshape(-1, 2)
        frame = pd.DataFrame({"x": pts[:, 0], "y": pts[:, 1]})
    else:
        space = int(pp_space) if pp_space else 5
        areas = np.asarray([g.area for g in model.vor.gdf_vorPolys.geometry], dtype=float)
        step = float(np.sqrt(np.median(areas))) * space
        bx = np.floor((xc - xc.min()) / step).astype(int)
        by = np.floor((yc - yc.min()) / step).astype(int)
        chosen = {}
        for cell, key in enumerate(zip(bx.tolist(), by.tolist(), strict=False)):
            # nearest cell to the net-node centre of its bin
            node = (xc.min() + (key[0] + 0.5) * step, yc.min() + (key[1] + 0.5) * step)
            d2 = (xc[cell] - node[0]) ** 2 + (yc[cell] - node[1]) ** 2
            if key not in chosen or d2 < chosen[key][1]:
                chosen[key] = (cell, d2)
        cells = sorted(item[0] for item in chosen.values())
        frame = pd.DataFrame({"x": xc[cells], "y": yc[cells]})
    frame.insert(0, "ppname", [f"pp{i:04d}" for i in range(len(frame))])
    return frame


def _write_template(frame: pd.DataFrame, csv_path: Path) -> Path:
    """Write a pp value CSV plus its PEST template (one parameter per point)."""

    tpl_path = csv_path.with_suffix(csv_path.suffix + ".tpl")
    frame.loc[:, ["parnme", "value"]].to_csv(csv_path, index=False)
    with tpl_path.open("w", encoding="utf-8") as handle:
        handle.write("ptf ~\n")
        templated = frame.loc[:, ["parnme", "value"]].copy()
        templated["value"] = templated["parnme"].map(lambda name: f"~  {name}  ~")
        templated.to_csv(handle, index=False)
    return tpl_path


def add_pilot_point_parameter(project, spec) -> pd.DataFrame:
    """Compile a ``style='pilotpoints'`` K spec onto the native PEST workspace.

    Places pilot points, writes the pp value template + support tables, injects
    the IDW apply as a pre-model command in the native forward run, and records
    the parameter frame so :meth:`PestProject.build` can register the template
    parameters after ``build_pst``.
    """

    files = _resolve_files(project.template_workspace, project.model.name, spec.recipe)
    if spec.layers is not None:
        wanted = {int(layer) for layer in spec.layers}
        files = [
            name for name in files
            if (m := re.search(r"_layer(\d+)\.txt$", name)) and (int(m.group(1)) - 1) in wanted
        ] or files
    spec.resolved_files = list(files)
    for name in files:
        _flatten_array_file(Path(project.template_workspace) / name)

    pp = place_pilot_points(project.model, pp_space=spec.pp_space, pp_points=spec.pp_points)
    xc, yc = _cell_centers(project.model)
    kdata = np.asarray(project.model.gwf.npf.k.get_data(), dtype=float)
    helper = str(Path(__file__).with_name("forward_run.py"))
    lo, hi = (spec.physical if spec.physical is not None else (-1.0e30, 1.0e30))
    frames = []
    for name in files:
        match = re.search(r"_layer(\d+)\.txt$", name)
        layer = int(match.group(1)) - 1 if match else 0
        base = f"{spec.name}l{layer}"
        frame = pd.DataFrame({
            "parnme": [f"{base}pp{i:04d}" for i in range(len(pp))],
            "pargp": base,
            "partrans": spec.resolved_transform,
            "parval1": 1.0,
            "value": 1.0,
            "parlbnd": float(spec.bounds[0]),
            "parubnd": float(spec.bounds[1]),
            "x": pp["x"].to_numpy(),
            "y": pp["y"].to_numpy(),
            "layer": layer,
        })
        pp_csv = project.template_workspace / f"{base}_pp.csv"
        _write_template(frame, pp_csv)
        points_csv = project.template_workspace / f"{base}_pp.points.csv"
        cells_csv = project.template_workspace / f"{base}_pp.cells.csv"
        frame.loc[:, ["parnme", "x", "y", "layer"]].to_csv(points_csv, index=False)
        pd.DataFrame({
            "cell": np.arange(kdata.shape[1]),
            "x": xc, "y": yc, "base_k": kdata[layer], "layer": layer,
        }).to_csv(cells_csv, index=False)
        project.pf.add_py_function(
            helper,
            "apply_pilotpoints_to_array("
            f"parameter_csv='{pp_csv.name}', points_meta_csv='{points_csv.name}', "
            f"cells_meta_csv='{cells_csv.name}', k_file='{name}', "
            f"lower_limit={lo!r}, upper_limit={hi!r})",
            is_pre_cmd=True,
        )
        frames.append((frame, pp_csv.with_suffix(pp_csv.suffix + ".tpl")))
    project._pilot_point_frames[spec.name] = frames
    return pd.concat([f for f, _ in frames], ignore_index=True)


def register_pilot_point_parameters(project) -> None:
    """Register pilot-point template parameters after ``pf.build_pst``."""

    pst = project.pst
    for frames in project._pilot_point_frames.values():
        for frame, tpl_path in frames:
            pst.add_parameters(str(tpl_path), pst_path=".")
            names = frame["parnme"].tolist()
            data = pst.parameter_data
            data.loc[names, "parnme"] = names
            data.loc[names, "pargp"] = frame["pargp"].to_numpy()
            data.loc[names, "partrans"] = frame["partrans"].to_numpy()
            data.loc[names, "parval1"] = frame["parval1"].to_numpy()
            data.loc[names, "parlbnd"] = frame["parlbnd"].to_numpy()
            data.loc[names, "parubnd"] = frame["parubnd"].to_numpy()
