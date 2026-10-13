"""Native pilot points on any MF6 grid: placement, interpolation factors, apply step.

pyEMU's built-in pilot points fall back to pure-Python kriging on unstructured
(DISV/Voronoi) grids -- minutes per build -- and its fast ``pypestutils`` path
either does not engage for unstructured grids or explodes into a non-stationary
hyperparameter setup. So myflopy places the points itself and **interpolates them
to cells by inverse-distance weighting** in the forward run. This module is the
build-time half (placement, interpolation factors, template files, parameter
registration); the forward-run half is ``forward_run.apply_pilotpoints_to_array``.

A pilot-point SET is one ``parameterize`` call's points on one array file (one
model layer) and one zone. Each cell of the set's zone is interpolated from the
``pp_nearest`` nearest points of THAT set, by inverse distance squared, in log10
space for a log-transformed parameter, with the weights computed ONCE at build and
stored in the template (``<set>_pp.factors.npz``). Cells outside every set keep
their value. Every set on one array file is applied by ONE forward-run command,
run after pyEMU's ``apply_list_and_array_pars``, which MULTIPLIES the file as
pyEMU just wrote it (or, for a file pyEMU does not parameterize, the base captured
at build) by the product of every set covering each cell, so pilot points stack
with constant/zone/grid multipliers and with each other (ledger 525).
The sets, their files and the base each file starts from are listed in
:data:`PILOT_POINT_MANIFEST` in the template.

Exposed through the one unified API -- three alternatives (several calls on one
project each need their own ``name=``, ledger 523)::

    cal.parameterize("k", style="pilotpoints", pp_spacing=1500.0, correlation=3000.0,
                     layers=[0, 1], bounds=(0.05, 20), physical=(1e-3, 300))
    cal.parameterize("k", style="pilotpoints", pp_points=my_gdf)   # explicit
    cal.parameterize("k", style="pilotpoints", pp_spacing=500.0, zones=units)  # a net per unit
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from myflopy._logging import get_logger
from myflopy.modflow.mf6.grid.geometry import cell_shares, model_cell_polygons, where_they_are
from myflopy.modflow.mf6.pest.native_parameters import (
    DEFAULT_PP_NEAREST,
    _flatten_array_file,
    array_row_length,
    staged_model,
    target_files,
)
from myflopy.modflow.mf6.pest.zones import resolve_zone_array

logger = get_logger(__name__)

#: The template file naming every pilot-point set and, per array file, the base the
#: forward run multiplies (ledger 525). Read by
#: ``forward_run.apply_pilotpoints_to_array``; a write-back finds a set's points,
#: zone, layer and interpolation factors here.
PILOT_POINT_MANIFEST = "pilot_point_sets.json"


def _cell_centers(gwf) -> tuple[np.ndarray, np.ndarray]:
    """A FloPy model's grid's cell-center ``(x, y)`` coordinate arrays, flattened."""

    grid = gwf.modelgrid
    return (
        np.asarray(grid.xcellcenters, dtype=float).reshape(-1),
        np.asarray(grid.ycellcenters, dtype=float).reshape(-1),
    )


def _point_array(pp_points) -> np.ndarray:
    """Explicit pilot-point locations as an ``(n, 2)`` float array."""

    geometry = getattr(pp_points, "geometry", None)
    if geometry is not None:
        return np.array([(geom.x, geom.y) for geom in geometry], dtype=float).reshape(-1, 2)
    return np.asarray(pp_points, dtype=float).reshape(-1, 2)


def _locate_points(gwf, points: np.ndarray, ncells: int) -> np.ndarray:
    """The cell each explicit pilot point lies in, by the shared point-in-cell rule.

    :func:`~myflopy.modflow.mf6.grid.geometry.cell_shares` on the cell polygons of
    ``gwf``'s ``modelgrid`` -- the build's STAGED model, never the caller's object
    (ledger 356) -- so a point on a face or vertex goes to the lowest-numbered cell it
    touches. A point in no cell is refused, naming it, with both extents.
    """

    import shapely

    polys = model_cell_polygons(gwf)
    if len(polys) != ncells:
        raise ValueError(
            f"Explicit pilot points are located on the model's {len(polys)} cell "
            f"polygons, but the array they scale has {ncells} values per file, so a "
            "point's cell cannot be numbered in it. Use pp_spacing= to lay a net instead."
        )
    cells = np.full(len(points), -1, dtype=int)
    for index, (x, y) in enumerate(points):
        shares = cell_shares(polys, shapely.Point(float(x), float(y)))
        if shares:
            cells[index] = int(shares[0][0])
    off = np.flatnonzero(cells < 0)
    if off.size:
        import geopandas as gpd

        first = int(off[0])
        missed = gpd.GeoSeries(shapely.points(points[off]), crs=polys.crs)
        raise ValueError(
            f"{off.size} of {len(points)} explicit pilot point(s) lie in no cell of the "
            f"grid -- point {first} at ({points[first, 0]:g}, {points[first, 1]:g}) is "
            "the first. A pilot point outside the grid moves no cell it is meant to. "
            + where_they_are(missed, polys, noun="pilot points", source_crs=polys.crs)
        )
    return cells


def place_pilot_points(gwf, *, pp_spacing=None, pp_points=None, zones=None) -> pd.DataFrame:
    """Return a pilot-point table (``ppname, x, y, zone, cell``) for the model grid.

    Give exactly one of ``pp_spacing`` (a regular net) and ``pp_points`` (explicit
    locations). A net is laid over the cell centres of ``gwf`` (the build's STAGED
    flow model, ledger 356) in squares ``pp_spacing`` wide, from the lowest cell
    centre's ``x`` and ``y``; in each square holding a cell centre, the cell whose
    centre is nearest the square's middle becomes a point (the lowest-numbered cell
    on a tie). With ``zones``, each zone label ``> 0`` gets its own net over ITS
    cells on the same squares, so a zone holding any cell holds a point and zoned
    nets are subsets of one lattice; cells labelled ``<= 0`` get no point
    (ledger 527).

    Parameters
    ----------
    gwf : flopy.mf6.ModflowGwf
        The FloPy flow model whose ``modelgrid`` is read -- the build's STAGED copy
        (:func:`~myflopy.modflow.mf6.pest.native_parameters.staged_model`), whose
        cells are the template's as written: its cell centres, which a net is
        snapped to, and with ``pp_points`` its cell polygons
        (:func:`~myflopy.modflow.mf6.grid.geometry.model_cell_polygons`), which
        each explicit point is located in. The polygons were read off the CALLER's
        model object until every grid the build reads came from the staged copy
        (ledger 356). Their CRS is FloPy's ``modelgrid.crs`` -- ``None`` when the
        run's files declare none -- and is read only to word the refusal of a point
        off the grid; the points are never reprojected.
    pp_spacing : float or None, optional
        Net spacing in model length units (the grid's own, e.g. feet), finite and
        greater than zero. ``None`` (default) requires ``pp_points``.
    pp_points : geopandas.GeoDataFrame or array_like or None, optional
        Explicit point locations -- Point geometries, or ``(x, y)`` pairs
        reshapable to ``(n, 2)`` -- in the grid's own coordinates; a GeoDataFrame's
        CRS is not read or reprojected. Each point takes the zone of the cell it
        lies in (the lowest-numbered cell on a face or vertex). Default ``None``
        requires ``pp_spacing``.
    zones : array_like of int or None, optional
        One integer zone label per cell of a layer (flat, in the order of
        ``gwf``'s cell centres). ``None`` (default) makes one net over every cell,
        each point's ``zone`` reading ``0``.

    Returns
    -------
    pandas.DataFrame
        Columns ``ppname`` (``pp0000``, ... numbered within each zone), ``x``,
        ``y``, ``zone`` and ``cell`` (the cell the point is at or in), one row per
        point; a net's points in zone, then cell, order, explicit points in the
        order given within each zone.

    Raises
    ------
    ValueError
        When neither or both of ``pp_spacing``/``pp_points`` are given, the spacing
        is not a finite number above zero, ``zones`` is not one label per cell; when
        an explicit point lies in no cell, in a cell labelled ``<= 0``, or at the
        location of another; or when a zone ``> 0`` is left with no explicit point.
    """

    xc, yc = _cell_centers(gwf)
    if zones is None:
        labels = np.zeros(xc.size, dtype=int)
    else:
        labels = np.asarray(zones, dtype=int).reshape(-1)
        if labels.size != xc.size:
            raise ValueError(
                f"zones holds {labels.size} labels and the grid {xc.size} cells; "
                "give one integer label per cell."
            )
    if (pp_spacing is None) == (pp_points is None):
        raise ValueError(
            "Pilot points are placed by pp_spacing= (a regular net, in model length "
            "units) OR pp_points= (explicit locations): give exactly one."
        )
    if pp_points is not None:
        points = _point_array(pp_points)
        if not len(points):
            raise ValueError("pp_points= holds no point.")
        _, first, counts = np.unique(points, axis=0, return_index=True, return_counts=True)
        if (counts > 1).any():
            twice = int(first[np.flatnonzero(counts > 1)[0]])
            raise ValueError(
                f"Explicit pilot point {twice} at ({points[twice, 0]:g}, "
                f"{points[twice, 1]:g}) is given more than once; two parameters at one "
                "location are one parameter twice (their prior correlation is 1)."
            )
        cells = _locate_points(gwf, points, xc.size)
        zone = labels[cells]
        if zones is not None:
            outside = np.flatnonzero(zone <= 0)
            if outside.size:
                bad = int(outside[0])
                raise ValueError(
                    f"Explicit pilot point {bad} at ({points[bad, 0]:g}, "
                    f"{points[bad, 1]:g}) lies in cell {int(cells[bad])}, whose zone is "
                    f"{int(zone[bad])}: a label <= 0 is not parameterized, so the point "
                    "would move nothing. Move it into a zone, or drop it."
                )
            empty = sorted(set(labels[labels > 0].tolist()) - set(zone.tolist()))
            if empty:
                raise ValueError(
                    f"Zone {empty[0]} ({int((labels == empty[0]).sum())} cells) holds none "
                    "of the explicit pilot points, so its cells have nothing to be "
                    "interpolated from. Add a point in it, give it label 0 to leave it "
                    "alone, or lay a net with pp_spacing=."
                )
        frame = pd.DataFrame({"x": points[:, 0], "y": points[:, 1], "zone": zone, "cell": cells})
        frame = frame.sort_values("zone", kind="stable").reset_index(drop=True)
    else:
        step = float(pp_spacing)
        if not np.isfinite(step) or step <= 0:
            raise ValueError(f"pp_spacing={pp_spacing!r} must be a finite length above zero.")
        x0, y0 = xc.min(), yc.min()
        bx = np.floor((xc - x0) / step).astype(np.int64)
        by = np.floor((yc - y0) / step).astype(np.int64)
        d2 = (xc - (x0 + (bx + 0.5) * step)) ** 2 + (yc - (y0 + (by + 0.5) * step)) ** 2
        table = pd.DataFrame({"cell": np.arange(xc.size), "zone": labels, "bx": bx, "by": by,
                              "d2": d2})
        if zones is not None:
            table = table.loc[table["zone"] > 0]
        chosen = (
            table.sort_values(["d2", "cell"], kind="stable")
            .drop_duplicates(["zone", "bx", "by"])
            .sort_values(["zone", "cell"], kind="stable")
        )
        cells = chosen["cell"].to_numpy()
        frame = pd.DataFrame({"x": xc[cells], "y": yc[cells],
                              "zone": chosen["zone"].to_numpy(), "cell": cells})
    frame.insert(0, "ppname", [f"pp{i:04d}" for i in frame.groupby("zone").cumcount()])
    return frame


def idw_factors(cell_x, cell_y, point_x, point_y, nearest: int = DEFAULT_PP_NEAREST):
    """Inverse-distance-squared weights of each cell's ``nearest`` nearest points.

    The interpolation factors of one pilot-point set, computed once at build: each
    cell takes the ``min(nearest, n_points)`` points nearest its centre, weighted
    by ``1 / d**2`` and normalized to sum to one; a cell whose centre IS a point's
    takes that point alone. The forward run applies them as
    ``sum(weight * log10(value[index]))`` per cell (linear for a parameter that is
    not log-transformed) -- one gather and one row sum, no distances.

    Parameters
    ----------
    cell_x, cell_y : array_like of float
        The cell centres to interpolate to, model length units.
    point_x, point_y : array_like of float
        The set's pilot-point locations, model length units; at least one.
    nearest : int, optional
        How many nearest points each cell is interpolated from; ``>= 1``. Default
        :data:`DEFAULT_PP_NEAREST` (8).

    Returns
    -------
    tuple of (numpy.ndarray, numpy.ndarray)
        ``index`` (``int32``, shape ``(n_cells, k)``): each cell's points, nearest
        first, as positions in ``point_x``; ``weight`` (``float64``, same shape):
        their weights, each row summing to one.
    """

    from scipy.spatial import cKDTree

    cells = np.column_stack([np.asarray(cell_x, dtype=float), np.asarray(cell_y, dtype=float)])
    points = np.column_stack([np.asarray(point_x, dtype=float), np.asarray(point_y, dtype=float)])
    k = max(1, min(int(nearest), len(points)))
    distance, index = cKDTree(points).query(cells, k=k)
    distance = np.asarray(distance, dtype=float).reshape(len(cells), k)
    index = np.asarray(index).reshape(len(cells), k)
    d2 = distance**2
    exact = d2[:, 0] == 0.0
    weight = np.zeros_like(d2)
    weight[~exact] = 1.0 / d2[~exact]
    weight[exact, 0] = 1.0
    weight /= weight.sum(axis=1, keepdims=True)
    return index.astype(np.int32), weight


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
    """Compile a ``style='pilotpoints'`` array spec onto the native PEST workspace.

    Places the spec's points, then for each of its layer files and each zone writes
    a SET: the pp value template, ``<set>_pp.points.csv`` (``parnme, x, y, layer,
    zone``; ``layer`` blank for a set on period files, as the manifest's ``None``)
    and ``<set>_pp.factors.npz`` (the cells it covers, their centres and
    :func:`idw_factors`), and records it on ``project._pilot_point_sets`` so that
    :func:`emit_pilot_point_commands` adds one forward-run command per array file
    and :func:`register_pilot_point_parameters` registers the parameters after
    ``build_pst``. A set is named ``<name>l<layer>`` -- ``<name>l<layer>z<zone>``
    with ``zones=`` -- lower-cased, which is its PEST group and its parameters'
    stem (``...pp0000``).

    Parameters
    ----------
    project : PestProject
        The calibration being built; its template workspace, staged model and
        ``_pilot_point_sets`` are read and written.
    spec : NativeParameterSpec
        An array-target ``style="pilotpoints"`` spec. ``spec.resolved_files`` and
        ``spec.pilot_point_sets`` (one summary per set: ``name``, ``layer``,
        ``zone``, ``points``, ``cells``) are filled in.

    Returns
    -------
    pandas.DataFrame
        Every set's parameter rows (``parnme, pargp, partrans, parval1, value,
        parlbnd, parubnd, x, y, layer, zone``); ``layer`` is ``None`` for a set on
        period files.

    Raises
    ------
    ValueError
        When ``zones=`` labels no cell above zero; when an array file does not hold
        one value per cell of the grid; when two sets take one name (which
        :meth:`PestProject.parameterize` already prevents by refusing a repeated
        ``name=``, ledger 523); or what :func:`place_pilot_points` refuses.
    """

    # FloPy's own records of the staged model, by package TYPE, narrowed to
    # `layers=` (ledger 517); sets `spec.resolved_files`.
    files = [target.filename for target in target_files(project, spec)]

    # The cell centres the template's model holds, as written to its files: the
    # build's staged copy, not the caller's object (ledger 356).
    gwf = staged_model(project, project.model.name)
    row_length = array_row_length(gwf.modelgrid)
    template = Path(project.template_workspace)
    for name in files:
        _flatten_array_file(template / name, row_length)

    labels = None
    if spec.zones is not None:
        labels = np.asarray(
            resolve_zone_array(spec.zones, family="array", model=gwf)
        ).reshape(-1)
        # The labels the zone array holds, as `style="zone"` records them for the
        # build's metadata; the sets are the labels above 0.
        spec.resolved_zone_labels = sorted(int(label) for label in np.unique(labels))
        if not (labels > 0).any():
            raise ValueError(
                f"zones= labels no cell above 0 for {spec.recipe.canonical!r} pilot points "
                f"(name {spec.name!r}), so it parameterizes nothing: a label <= 0 is "
                "left alone. Give the cells to calibrate a positive label."
            )
    pp = place_pilot_points(gwf, pp_spacing=spec.pp_spacing, pp_points=spec.pp_points,
                            zones=labels)
    xc, yc = _cell_centers(gwf)
    for name in files:
        n_values = len((template / name).read_text().split())
        if n_values != xc.size:
            raise ValueError(
                f"Pilot points on {spec.recipe.canonical!r} interpolate to the grid's "
                f"{xc.size} cell centres, and {name} holds {n_values} values, not one "
                "per cell; such a file cannot take pilot points."
            )
    lower, upper = spec.physical if spec.physical is not None else (None, None)
    interpolation = "log" if spec.resolved_transform == "log" else "linear"
    periodic = spec.recipe.periodic
    taken = {item["record"]["name"] for item in project._pilot_point_sets}
    frames, summaries = [], []
    # A set per (period group, layer, zone). A LAYER file is a set of its own; the
    # PERIOD files of one group share their sets -- one field interpolated once per
    # run and multiplied into each file (ledger 547). The layer is
    # FloPy's record of the file, never parsed from its name (ledger 517).
    for group, pargps in zip(spec.period_groups, spec.period_group_pargps, strict=True):
        by_layer: dict[int | None, list[str]] = {}
        for target in group.targets:
            by_layer.setdefault(target.layer, []).append(target.filename)
        for file_layer, array_files in by_layer.items():
            # A whole-grid griddata file is layer 0, as it always was; period files
            # have no layer and carry none in their name.
            layer = file_layer if file_layer is not None else (None if periodic else 0)
            token = "" if layer is None else f"l{layer}"
            for zone, points in pp.groupby("zone", sort=True):
                zone = int(zone)
                base = (f"{spec.name}{group.name_suffix}{token}"
                        + ("" if labels is None else f"z{zone}")).lower()
                if base in taken:
                    raise ValueError(
                        f"Two pilot-point sets are named {base!r}: two style='pilotpoints' "
                        f"calls on {array_files[0]} under name {spec.name!r}. PEST would see "
                        "one set's parameters and the forward run apply both; give each call "
                        "its own name= (each set is named <name>[_<period group>]"
                        "[l<layer>][z<zone>])."
                    )
                taken.add(base)
                pargps.append(base)
                cells = np.arange(xc.size) if labels is None else np.flatnonzero(labels == zone)
                frame = pd.DataFrame({
                    "parnme": [f"{base}pp{i:04d}" for i in range(len(points))],
                    "pargp": base,
                    "partrans": spec.resolved_transform,
                    "parval1": 1.0,
                    "value": 1.0,
                    "parlbnd": float(spec.bounds[0]),
                    "parubnd": float(spec.bounds[1]),
                    "x": points["x"].to_numpy(dtype=float),
                    "y": points["y"].to_numpy(dtype=float),
                    # The manifest's value: None (a blank in the points CSV) for
                    # period files, which have no layer -- never a -1 sentinel.
                    "layer": layer,
                    "zone": zone if labels is not None else 0,
                })
                pp_csv = template / f"{base}_pp.csv"
                tpl_path = _write_template(frame, pp_csv)
                points_csv = template / f"{base}_pp.points.csv"
                frame.loc[:, ["parnme", "x", "y", "layer", "zone"]].to_csv(points_csv, index=False)
                index, weight = idw_factors(xc[cells], yc[cells], frame["x"], frame["y"],
                                            spec.pp_nearest)
                factors = template / f"{base}_pp.factors.npz"
                np.savez(factors, cells=cells.astype(np.int64), x=xc[cells], y=yc[cells],
                         index=index, weight=weight)
                record = {
                    "name": base,
                    "spec": spec.name,
                    "target": spec.recipe.canonical,
                    "array_files": list(array_files),
                    "layer": layer,
                    "zone": zone if labels is not None else None,
                    "period_group": group.label,
                    "periods": [int(period) for period in group.periods] if periodic else None,
                    "parameter_csv": pp_csv.name,
                    "template": tpl_path.name,
                    "points_csv": points_csv.name,
                    "factors": factors.name,
                    "interpolation": interpolation,
                    "nearest": int(index.shape[1]),
                    "lower_limit": None if lower is None else float(lower),
                    "upper_limit": None if upper is None else float(upper),
                    "points": len(frame),
                    "cells": int(cells.size),
                }
                project._pilot_point_sets.append(
                    {"record": record, "frame": frame, "tpl": tpl_path, "spec": spec}
                )
                summaries.append({key: record[key] for key in
                                  ("name", "layer", "zone", "period_group", "points", "cells")})
                frames.append(frame)
                logger.info(
                    "Pilot-point set %s: %d points over %d cells of %d file(s) (%s ...)", base,
                    len(frame), cells.size, len(array_files), array_files[0],
                )
    spec.pilot_point_sets = summaries
    return pd.concat(frames, ignore_index=True)


def _pyemu_model_files(pf) -> set[str]:
    """The model files pyEMU's ``apply_list_and_array_pars`` rewrites, by name."""

    if not getattr(pf, "par_dfs", None):
        return set()
    relations = pf.parfile_relations
    return {Path(str(name)).name for name in relations["model_file"]}


def _pyemu_limits(project, array_file: str) -> tuple[float | None, float | None]:
    """The tightest ``physical=`` of every pyEMU multiplier on ``array_file``.

    pyEMU clamps its own product to its ``ult_lbound``/``ult_ubound`` before the
    pilot-point step multiplies the file again, so a pilot-point factor above 1 would
    carry a cell past a constant's ``physical=`` unless the step clamps to it too
    (ledger 528). ``(None, None)`` when no such multiplier declares one.
    """

    lower = upper = None
    for spec in project._native_parameter_specs:
        if (spec.style == "pilotpoints" or spec.physical is None
                or spec.recipe.family != "array"
                or array_file not in {str(name) for name in spec.resolved_files}):
            continue
        low, high = float(spec.physical[0]), float(spec.physical[1])
        lower = low if lower is None else max(lower, low)
        upper = high if upper is None else min(upper, high)
    return lower, upper


def _refuse_limits_that_exclude_each_other(array_file: str, records: list[dict], template: Path,
                                           file_limits=(None, None)):
    """Refuse a covered cell whose ``physical=`` ranges -- every covering set's and
    the file's pyEMU multipliers' (``file_limits``) -- share no value."""

    covered = {record["name"]: np.load(template / record["factors"])["cells"] for record in records}
    size = 1 + max(int(cells.max()) for cells in covered.values())
    lower, upper = np.full(size, -np.inf), np.full(size, np.inf)
    if file_limits[0] is not None:
        lower[:] = file_limits[0]
    if file_limits[1] is not None:
        upper[:] = file_limits[1]
    for record in records:
        cells = covered[record["name"]]
        if record["lower_limit"] is not None:
            lower[cells] = np.maximum(lower[cells], record["lower_limit"])
        if record["upper_limit"] is not None:
            upper[cells] = np.minimum(upper[cells], record["upper_limit"])
    clash = np.flatnonzero(lower > upper)
    if clash.size:
        cell = int(clash[0])
        names = [record["name"] for record in records if cell in set(covered[record["name"]].tolist())]
        also = ("" if file_limits == (None, None) else
                f", with the physical= {file_limits} of the pyEMU multipliers on the file,")
        raise ValueError(
            f"Pilot-point sets {names} all cover cell {cell} of {array_file}, and their "
            f"physical= limits{also} share no value (at least {lower[cell]:g}, at most "
            f"{upper[cell]:g}); a cell is clamped to every covering limit at once."
        )


def emit_pilot_point_commands(project) -> None:
    """Add ONE forward-run command per pilot-point array file, and write the manifest.

    Called once by :meth:`PestProject.build` after every parameter is declared and
    before ``build_pst``. Each command runs after pyEMU's
    ``apply_list_and_array_pars`` (which ``PstFrom`` puts first) and multiplies every
    set on its file in. A file pyEMU ALSO parameterizes -- a constant, zone or grid
    multiplier on it -- is multiplied as pyEMU has just rewritten it from its
    ``org/`` copy (``base_from_file``); any other is rebuilt from
    ``<file>.pp_base``, the file as it stood at build, copied here. Either way a
    worker directory reused across realizations gets the same file every run
    (ledger 525). The file's entry also carries the tightest
    ``physical=`` of the pyEMU multipliers on it (``lower_limit``/``upper_limit``),
    which the forward run clamps every covered cell to as well, so a pilot-point
    factor cannot carry a cell past a constant's limits (ledger 528).

    Parameters
    ----------
    project : PestProject
        The calibration being built; its ``_pilot_point_sets``, ``pf`` and template
        workspace are read, and :data:`PILOT_POINT_MANIFEST` is written there.

    Raises
    ------
    ValueError
        When the ``physical=`` limits covering one cell -- its pilot-point sets' and
        the pyEMU multipliers' on its file -- share no value.
    """

    sets = [item["record"] for item in project._pilot_point_sets]
    if not sets:
        return
    template = Path(project.template_workspace)
    managed = _pyemu_model_files(project.pf)
    by_file: dict[str, list[dict]] = {}
    for record in sets:
        for array_file in record["array_files"]:
            by_file.setdefault(array_file, []).append(record)
    helper = str(Path(__file__).with_name("forward_run.py"))
    files = []
    # Files covered by the SAME sets share one command, so each set is interpolated
    # once per run however many files it multiplies -- every period file of an array
    # RCH (ledger 547); a K layer file has sets of its own and a command
    # of its own, as before.
    commands: dict[tuple[str, ...], list[str]] = {}
    for array_file, records in by_file.items():
        from_file = array_file in managed
        limits = _pyemu_limits(project, array_file) if from_file else (None, None)
        _refuse_limits_that_exclude_each_other(array_file, records, template, limits)
        base_file = None
        if not from_file:
            base_file = f"{array_file}.pp_base"
            shutil.copyfile(template / array_file, template / base_file)
        files.append({"array_file": array_file, "base_from_file": from_file,
                      "base_file": base_file, "lower_limit": limits[0],
                      "upper_limit": limits[1],
                      "sets": [record["name"] for record in records]})
        commands.setdefault(tuple(record["name"] for record in records), []).append(array_file)
    for array_files in commands.values():
        project.pf.add_py_function(
            helper,
            f"apply_pilotpoints_to_array(array_files={array_files!r}, "
            f"manifest={PILOT_POINT_MANIFEST!r})",
            is_pre_cmd=True,
        )
    (template / PILOT_POINT_MANIFEST).write_text(
        json.dumps({"version": 1, "files": files, "sets": sets}, indent=1), encoding="utf-8"
    )


def register_pilot_point_parameters(project) -> None:
    """Register every pilot-point set's template parameters after ``pf.build_pst``."""

    pst = project.pst
    for item in project._pilot_point_sets:
        frame = item["frame"]
        pst.add_parameters(str(item["tpl"]), pst_path=".")
        names = frame["parnme"].tolist()
        data = pst.parameter_data
        data.loc[names, "parnme"] = names
        data.loc[names, "pargp"] = frame["pargp"].to_numpy()
        data.loc[names, "partrans"] = frame["partrans"].to_numpy()
        data.loc[names, "parval1"] = frame["parval1"].to_numpy()
        data.loc[names, "parlbnd"] = frame["parlbnd"].to_numpy()
        data.loc[names, "parubnd"] = frame["parubnd"].to_numpy()
