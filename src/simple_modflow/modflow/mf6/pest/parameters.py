"""Parameter support-file builders for the first ``simple_modflow`` PEST slice."""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

from simple_modflow.modflow.mf6.pest.gis import derive_bounds, load_vector_parameter_source
from simple_modflow.modflow.mf6.pest.geostats import build_geostruct
from simple_modflow.modflow.mf6.pest.specs import (
    DrainConductanceParameter,
    DrainElevationParameter,
    KPilotPointParameter,
)
from simple_modflow.modflow.mf6.drn import DRN as DRNFromVector


def _derive_parameter_value_bounds(
    *,
    base: pd.Series,
    bounds: tuple[float, float] | None,
    bounds_mode: str,
    lower_col: pd.Series | None = None,
    upper_col: pd.Series | None = None,
    parameter_space: str = "absolute",
) -> tuple[pd.Series, pd.Series]:
    """Derive bounds in the same space as the adjustable parameter values.

    For multiplier-style parameters, the adjustable value is a factor (often
    starting from ``1.0``), not the absolute model property value. In that case
    bounds must also live in multiplier space.
    """

    mode = str(bounds_mode).strip().lower()
    space = str(parameter_space).strip().lower()
    if space != "multiplier":
        return derive_bounds(
            base,
            bounds,
            bounds_mode,
            lower_bound_column=lower_col,
            upper_bound_column=upper_col,
        )

    if mode == "multiplier":
        if bounds is None:
            raise ValueError("Multiplier bounds mode requires a bounds tuple.")
        lower = pd.Series(float(bounds[0]), index=base.index)
        upper = pd.Series(float(bounds[1]), index=base.index)
    elif mode == "multiplier_from_columns":
        if lower_col is None or upper_col is None:
            raise ValueError(
                "multiplier_from_columns mode requires lower_bound_column and upper_bound_column."
            )
        lower = pd.to_numeric(lower_col, errors="coerce")
        upper = pd.to_numeric(upper_col, errors="coerce")
    elif mode == "absolute":
        if bounds is None:
            raise ValueError("Absolute bounds mode requires a bounds tuple.")
        lower = pd.Series(float(bounds[0]), index=base.index) / base
        upper = pd.Series(float(bounds[1]), index=base.index) / base
    elif mode == "from_columns":
        if lower_col is None or upper_col is None:
            raise ValueError("from_columns mode requires lower_bound_column and upper_bound_column.")
        lower = pd.to_numeric(lower_col, errors="coerce") / base
        upper = pd.to_numeric(upper_col, errors="coerce") / base
    else:
        raise ValueError(f"Unsupported bounds_mode {bounds_mode!r}.")

    invalid = lower.isna() | upper.isna() | (base == 0)
    if invalid.any():
        raise ValueError("Derived parameter bounds contain missing or invalid values.")
    swapped = lower > upper
    if swapped.any():
        raise ValueError("Derived parameter bounds contain lower values greater than upper values.")
    return lower.astype(float), upper.astype(float)


def _write_parameter_template(frame: pd.DataFrame, csv_path: Path) -> Path:
    """Write a CSV plus matching PEST template file for row-based parameters."""

    csv_path.parent.mkdir(parents=True, exist_ok=True)
    template_path = csv_path.with_suffix(csv_path.suffix + ".tpl")
    out = frame.loc[:, ["parnme", "value"]].copy()
    out.to_csv(csv_path, index=False)
    templated = out.copy()
    templated["value"] = templated["parnme"].map(lambda name: f"~   {name}   ~")
    with template_path.open("w", encoding="utf-8") as handle:
        handle.write("ptf ~\n")
        templated.to_csv(handle, index=False)
    return template_path


def _register_template_parameters(project, frame: pd.DataFrame, template_path: Path):
    """Add template-file parameters to the built ``pyemu.Pst`` object."""

    pst = project.pst
    new_pars = pst.add_parameters(str(template_path), pst_path=".")
    parnames = frame["parnme"].tolist()
    par_data = pst.parameter_data
    par_data.loc[parnames, "parnme"] = parnames
    par_data.loc[parnames, "pargp"] = frame["pargp"].to_numpy()
    par_data.loc[parnames, "parval1"] = frame["parval1"].astype(float).to_numpy()
    par_data.loc[parnames, "parlbnd"] = frame["parlbnd"].astype(float).to_numpy()
    par_data.loc[parnames, "parubnd"] = frame["parubnd"].astype(float).to_numpy()
    par_data.loc[parnames, "partrans"] = frame["partrans"].to_numpy()
    return new_pars


def _build_k_support(project, spec: KPilotPointParameter):
    """Build pilot-point parameter rows plus cell/point support metadata."""

    source = load_vector_parameter_source(spec.source)
    centers = project.model.vor.gdf_vorPolys.copy()
    centers["cell"] = centers.index.astype(int)
    centers["x"] = centers.geometry.centroid.x
    centers["y"] = centers.geometry.centroid.y
    points = gpd.GeoDataFrame(
        centers.loc[:, ["cell", "x", "y"]],
        geometry=gpd.points_from_xy(centers["x"], centers["y"]),
        crs=centers.crs,
    )
    joined = gpd.sjoin(points, source, how="inner", predicate="intersects")
    if joined.empty:
        raise ValueError("No Voronoi cell centers intersect the K parameter source polygons.")

    spacing = float(spec.pp_spacing or 1.0)
    xmin = float(joined["x"].min())
    ymin = float(joined["y"].min())
    joined["pp_i"] = np.floor((joined["x"] - xmin) / spacing).astype(int)
    joined["pp_j"] = np.floor((joined["y"] - ymin) / spacing).astype(int)
    if spec.source.zone_column and spec.source.zone_column in joined.columns:
        group_columns = [spec.source.zone_column, "pp_i", "pp_j"]
    else:
        group_columns = ["pp_i", "pp_j"]
    grouped = joined.sort_values(["y", "x"]).groupby(group_columns, as_index=False).first()

    base = pd.to_numeric(grouped[spec.source.value_column], errors="coerce")
    lower_col = (
        grouped[spec.source.lower_bound_column]
        if spec.source.lower_bound_column and spec.source.lower_bound_column in grouped.columns
        else None
    )
    upper_col = (
        grouped[spec.source.upper_bound_column]
        if spec.source.upper_bound_column and spec.source.upper_bound_column in grouped.columns
        else None
    )
    lower, upper = _derive_parameter_value_bounds(
        base=base,
        bounds=spec.bounds,
        bounds_mode=spec.bounds_mode,
        lower_col=lower_col,
        upper_col=upper_col,
        parameter_space=spec.parameter_space,
    )
    zone_values = (
        grouped[spec.source.zone_column].astype(str).tolist()
        if spec.source.zone_column and spec.source.zone_column in grouped.columns
        else ["all"] * len(grouped)
    )
    frame = pd.DataFrame(
        {
            "parnme": [f"{spec.name}_pp_{i:04d}" for i in range(len(grouped))],
            "pargp": [f"{spec.name}_pp"] * len(grouped),
            "partrans": [spec.transform] * len(grouped),
            "parval1": [1.0 if spec.parameter_space == "multiplier" else value for value in base],
            "parlbnd": lower.to_numpy(),
            "parubnd": upper.to_numpy(),
            "value": [1.0 if spec.parameter_space == "multiplier" else value for value in base],
            "x": grouped["x"].to_numpy(),
            "y": grouped["y"].to_numpy(),
            "zone": zone_values,
            "base_value": base.to_numpy(),
            "layer": [int(spec.layers[0])] * len(grouped),
        }
    )

    cells_meta = joined.loc[:, ["cell", "x", "y"]].copy()
    cells_meta["zone"] = (
        joined[spec.source.zone_column].astype(str).to_numpy()
        if spec.source.zone_column and spec.source.zone_column in joined.columns
        else "all"
    )
    cells_meta["base_k"] = pd.to_numeric(joined[spec.source.value_column], errors="coerce").to_numpy()
    cells_meta["layer"] = int(spec.layers[0])
    cells_meta = cells_meta.drop_duplicates(subset=["layer", "cell"]).reset_index(drop=True)
    points_meta = frame.loc[:, ["parnme", "x", "y", "zone", "base_value", "layer"]].rename(
        columns={"base_value": "base_k"}
    )
    return frame, points_meta, cells_meta


def build_k_pilotpoint_frame(project, spec: KPilotPointParameter) -> pd.DataFrame:
    """Build a pilot-point parameter table from Voronoi cells and GIS zones."""

    frame, _, _ = _build_k_support(project, spec)
    return frame


def build_drain_parameter_frame(
    spec: DrainElevationParameter | DrainConductanceParameter,
    *,
    parameter_kind: str,
) -> pd.DataFrame:
    """Build a drain-feature or drain-group parameter table."""

    source = load_vector_parameter_source(spec.source)
    value_col = spec.source.value_column
    base = pd.to_numeric(source[value_col], errors="coerce")
    lower_col = (
        source[spec.source.lower_bound_column]
        if spec.source.lower_bound_column and spec.source.lower_bound_column in source.columns
        else None
    )
    upper_col = (
        source[spec.source.upper_bound_column]
        if spec.source.upper_bound_column and spec.source.upper_bound_column in source.columns
        else None
    )

    if "group" in spec.parameter_style and spec.source.group_column and spec.source.group_column in source.columns:
        group_key = source[spec.source.group_column].astype(str)
        grouped = pd.DataFrame(
            {
                "group": group_key,
                "base_value": base,
                "lower_src": lower_col if lower_col is not None else np.nan,
                "upper_src": upper_col if upper_col is not None else np.nan,
            }
        ).groupby("group", as_index=False).first()
        ids = grouped["group"].astype(str)
        base = grouped["base_value"]
        lower, upper = _derive_parameter_value_bounds(
            base=base,
            bounds=spec.bounds,
            bounds_mode=spec.bounds_mode,
            lower_col=grouped["lower_src"],
            upper_col=grouped["upper_src"],
            parameter_space="multiplier"
            if isinstance(spec, DrainConductanceParameter) and "multiplier" in spec.parameter_style
            else "absolute",
        )
    else:
        if spec.source.feature_id_column is None:
            raise ValueError("Drain feature parameters require feature_id_column when not grouping.")
        ids = source[spec.source.feature_id_column].astype(str)
        lower, upper = _derive_parameter_value_bounds(
            base=base,
            bounds=spec.bounds,
            bounds_mode=spec.bounds_mode,
            lower_col=lower_col,
            upper_col=upper_col,
            parameter_space="multiplier"
            if isinstance(spec, DrainConductanceParameter) and "multiplier" in spec.parameter_style
            else "absolute",
        )

    if isinstance(spec, DrainConductanceParameter) and "multiplier" in spec.parameter_style:
        parval1 = np.ones(len(ids), dtype=float)
        value = np.ones(len(ids), dtype=float)
    else:
        parval1 = base.to_numpy(dtype=float)
        value = parval1.copy()

    return pd.DataFrame(
        {
            "parnme": [f"{spec.name}_{i}" for i in ids],
            "pargp": [spec.name] * len(ids),
            "partrans": [spec.transform] * len(ids),
            "parval1": np.zeros(len(ids), dtype=float)
            if isinstance(spec, DrainElevationParameter)
            else parval1,
            "parlbnd": lower.to_numpy(dtype=float),
            "parubnd": upper.to_numpy(dtype=float),
            "value": np.zeros(len(ids), dtype=float)
            if isinstance(spec, DrainElevationParameter)
            else value,
            "feature_id": ids.to_list(),
            "parameter_kind": [parameter_kind] * len(ids),
        }
    )


def _build_drain_support(project, spec: DrainElevationParameter | DrainConductanceParameter):
    """Build drain parameter rows plus stress-period support metadata."""

    frame = build_drain_parameter_frame(
        spec,
        parameter_kind="drn_elev" if isinstance(spec, DrainElevationParameter) else "drn_cond",
    )
    support_key = "feature_id"
    if "group" in spec.parameter_style and spec.source.group_column:
        support_key = "group"

    source = load_vector_parameter_source(spec.source)
    builder = DRNFromVector(
        model=project.model,
        shp_gpkg=Path(spec.source.path),
        uid=spec.source.feature_id_column or "name",
        crs=spec.source.crs or 2927,
        idomain=getattr(project.model, "_idomain", None),
    )
    fields = {
        "name": spec.source.feature_id_column or "name",
        "height_over_btm": spec.source.value_column,
        "conductance": spec.source.value_column if isinstance(spec, DrainConductanceParameter) else "cond",
        "layer": spec.source.layer_column or "layer",
        "min_elev": "min_elev",
    }
    region_rows = []
    nper = project.model.nper
    data = project.model.gwf.drn.stress_period_data.data
    for per in range(nper):
        records = pd.DataFrame(data[per]).copy()
        records["per"] = int(per)
        records["row_index"] = np.arange(len(records), dtype=int)
        layers = []
        cells = []
        for cellid in records["cellid"]:
            layers.append(int(cellid[0]))
            cells.append(int(cellid[1]))
        records["layer"] = layers
        records["cell"] = cells
        for name, row, active_cells in builder.iter_polygon_boundary_features(
            name_field=fields["name"],
            edges_only=False,
        ):
            row_layer = int(row[fields["layer"]]) - 1
            group_value = (
                str(row[spec.source.group_column])
                if spec.source.group_column and spec.source.group_column in row.index
                else None
            )
            matched = records[(records["layer"] == row_layer) & (records["cell"].isin(active_cells))]
            for match in matched.itertuples(index=False):
                region_rows.append(
                    {
                        "per": int(match.per),
                        "row_index": int(match.row_index),
                        "layer": int(match.layer),
                        "cell": int(match.cell),
                        "feature_id": str(name),
                        "group": group_value,
                        "base_elev": float(match.elev),
                        "base_cond": float(match.cond),
                    }
                )
    support = pd.DataFrame(region_rows).drop_duplicates(subset=["per", "row_index"])
    if support.empty:
        raise ValueError("No drain package rows could be matched to the provided drain parameter source.")
    key_column = "group" if support_key == "group" else "feature_id"
    par_map = frame.loc[:, ["parnme", "feature_id"]].rename(columns={"feature_id": key_column})
    support = support.merge(par_map, on=key_column, how="left")
    if support["parnme"].isna().any():
        raise ValueError("Some drain rows could not be mapped to parameter names.")
    return frame, support


def register_parameter_spec(project, spec):
    """Register one parameter specification with a built ``PestProject``."""
    prepared = project._prepared_parameters[spec.name]
    frame = prepared["frame"]
    template_path = prepared["template_path"]
    _register_template_parameters(project, frame, template_path)
    if spec.geostruct is not None:
        project._registered_geostructs[spec.name] = build_geostruct(spec.geostruct)
    project._parameter_frames[spec.name] = frame
    return frame


def prepare_parameter_spec(project, spec):
    """Write template/support files for a parameter spec before ``build_pst()``."""

    if isinstance(spec, KPilotPointParameter):
        frame, points_meta, cells_meta = _build_k_support(project, spec)
        csv_path = project.template_workspace / f"{spec.name}_pilot_points.csv"
        template_path = _write_parameter_template(frame, csv_path)
        points_meta_path = project.template_workspace / f"{spec.name}_pilot_points.meta.csv"
        cells_meta_path = project.template_workspace / f"{spec.name}_pilot_points.cells.csv"
        points_meta.to_csv(points_meta_path, index=False)
        cells_meta.to_csv(cells_meta_path, index=False)
        config = {
            "kind": "k_pilotpoints",
            "parameter_csv": csv_path.name,
            "points_meta_csv": points_meta_path.name,
            "cells_meta_csv": cells_meta_path.name,
            "parameter_space": spec.parameter_space,
        }
    elif isinstance(spec, DrainElevationParameter):
        frame, support = _build_drain_support(project, spec)
        csv_path = project.template_workspace / f"{spec.name}_drain_elevation.csv"
        template_path = _write_parameter_template(frame, csv_path)
        support_path = project.template_workspace / f"{spec.name}_drain_elevation.support.csv"
        support.to_csv(support_path, index=False)
        config = {
            "kind": "drn_elev",
            "parameter_csv": csv_path.name,
            "support_csv": support_path.name,
        }
    elif isinstance(spec, DrainConductanceParameter):
        frame, support = _build_drain_support(project, spec)
        csv_path = project.template_workspace / f"{spec.name}_drain_conductance.csv"
        template_path = _write_parameter_template(frame, csv_path)
        support_path = project.template_workspace / f"{spec.name}_drain_conductance.support.csv"
        support.to_csv(support_path, index=False)
        config = {
            "kind": "drn_cond",
            "parameter_csv": csv_path.name,
            "support_csv": support_path.name,
        }
    else:
        raise TypeError(f"Unsupported parameter spec type: {type(spec).__name__}")

    prepared = {"frame": frame, "template_path": template_path, "config": config}
    project._prepared_parameters[spec.name] = prepared
    return prepared
