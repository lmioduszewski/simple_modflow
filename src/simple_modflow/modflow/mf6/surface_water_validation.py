"""Validation helpers for coupled MF6 surface-water package workflows.

This module provides a small reporting layer for validating generated
``SFR``, ``LAK``, and ``MVR`` inputs before or after they are attached to a
model. The goal is not to replace MODFLOW 6's own parsing, but to catch the
most common geometry, indexing, and cross-package mistakes closer to the
builder layer with clearer error messages.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.sfr import SFR
    from simple_modflow.modflow.mf6.simulation.base import SimulationBase


_VALID_LAK_CONNECTION_TYPES = {"VERTICAL", "HORIZONTAL", "EMBEDDEDH", "EMBEDDEDV"}
_VALID_MVR_TYPES = {"FACTOR", "EXCESS", "THRESHOLD", "UPTO"}


@dataclass(slots=True)
class SurfaceWaterValidationIssue:
    """One validation finding from a surface-water configuration check.

    Parameters
    ----------
    level
        Either ``"error"`` or ``"warning"``.
    code
        Stable short code for the finding type.
    message
        Human-readable description of the issue.
    context
        Optional structured details such as package names, lake ids, or reach
        numbers that help pinpoint the issue.
    """

    level: str
    code: str
    message: str
    context: dict[str, Any] = field(default_factory=dict)


@dataclass
class SurfaceWaterValidationReport:
    """Collection of validation findings for one coupled surface-water setup.

    The report is intentionally lightweight:

    - use :meth:`add_error` and :meth:`add_warning` while validating
    - inspect :attr:`ok`, :attr:`errors`, and :attr:`warnings`
    - call :meth:`raise_for_errors` to fail fast with a concise message
    - call :meth:`to_frame` or :meth:`summary_frame` for notebook-friendly
      inspection
    """

    issues: list[SurfaceWaterValidationIssue] = field(default_factory=list)

    @property
    def errors(self) -> list[SurfaceWaterValidationIssue]:
        """Return only error-level findings."""

        return [issue for issue in self.issues if issue.level == "error"]

    @property
    def warnings(self) -> list[SurfaceWaterValidationIssue]:
        """Return only warning-level findings."""

        return [issue for issue in self.issues if issue.level == "warning"]

    @property
    def ok(self) -> bool:
        """Whether the report contains no error-level findings."""

        return len(self.errors) == 0

    def add_error(self, code: str, message: str, **context: Any) -> None:
        """Record one error-level finding."""

        self.issues.append(
            SurfaceWaterValidationIssue(
                level="error",
                code=code,
                message=message,
                context=context,
            )
        )

    def add_warning(self, code: str, message: str, **context: Any) -> None:
        """Record one warning-level finding."""

        self.issues.append(
            SurfaceWaterValidationIssue(
                level="warning",
                code=code,
                message=message,
                context=context,
            )
        )

    def extend(self, other: "SurfaceWaterValidationReport") -> None:
        """Append findings from another report."""

        self.issues.extend(other.issues)

    def summary(self) -> dict[str, int | bool]:
        """Return a compact count summary."""

        return {
            "ok": self.ok,
            "num_errors": len(self.errors),
            "num_warnings": len(self.warnings),
            "num_issues": len(self.issues),
        }

    def summary_frame(self) -> pd.DataFrame:
        """Return a one-row DataFrame summary for notebook display."""

        return pd.DataFrame([self.summary()])

    def to_frame(self) -> pd.DataFrame:
        """Return findings as a DataFrame.

        Returns
        -------
        pandas.DataFrame
            A tidy table with ``level``, ``code``, ``message``, and
            ``context`` columns. ``context`` is preserved as a dictionary so it
            can still be inspected programmatically.
        """

        return pd.DataFrame(
            [
                {
                    "level": issue.level,
                    "code": issue.code,
                    "message": issue.message,
                    "context": issue.context,
                }
                for issue in self.issues
            ]
        )

    def raise_for_errors(self, prefix: str = "Surface-water validation failed.") -> None:
        """Raise ``ValueError`` when any error-level findings are present."""

        if self.ok:
            return
        error_lines = [f"- [{issue.code}] {issue.message}" for issue in self.errors]
        raise ValueError(prefix + "\n" + "\n".join(error_lines))


def _coerce_int(value: Any) -> int | None:
    """Return ``value`` as ``int`` when possible, otherwise ``None``."""

    try:
        if isinstance(value, bool):
            return None
        return int(value)
    except (TypeError, ValueError):
        return None


def _coerce_float(value: Any) -> float | None:
    """Return ``value`` as ``float`` when possible, otherwise ``None``."""

    try:
        if value is None:
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def _normalize_lak_connection_rows(connectiondata: list | None) -> list[list]:
    """Accept both flat and legacy extra-nested LAK connection rows."""

    if connectiondata is None:
        return []
    if len(connectiondata) == 1 and isinstance(connectiondata[0], list):
        first = connectiondata[0]
        if len(first) > 0 and isinstance(first[0], list):
            return first
    return connectiondata


def _model_ncpl(model: "SimulationBase") -> int | None:
    """Return the number of plan-view cells for ``model`` when available."""

    if getattr(model, "vor", None) is not None:
        return int(model.vor.ncpl)
    grid = getattr(model, "modelgrid", None)
    if grid is not None and getattr(grid, "xcellcenters", None) is not None:
        return int(np.asarray(grid.xcellcenters).reshape(-1).size)
    return None


def _model_nlay(model: "SimulationBase") -> int | None:
    """Return the number of layers for ``model`` when available."""

    try:
        return int(model.nlay)
    except Exception:
        return None


def _cell_top_bottom(model: "SimulationBase", layer: int, cell: int) -> tuple[float | None, float | None]:
    """Return the top and bottom elevation for one zero-based cellid."""

    vor = getattr(model, "vor", None)
    if vor is not None and getattr(vor, "gdf_topbtm", None) is not None:
        df = vor.gdf_topbtm
        top_col = 0 if layer == 0 else layer
        botm_col = layer + 1
        if top_col in df.columns and botm_col in df.columns:
            try:
                return float(df.loc[cell, top_col]), float(df.loc[cell, botm_col])
            except Exception:
                pass

    grid = getattr(model, "modelgrid", None)
    if grid is None:
        return None, None

    top = getattr(grid, "top", None)
    botm = getattr(grid, "botm", None)
    try:
        if top is not None:
            top_arr = np.asarray(top, dtype=float).reshape(-1)
            cell_top = float(top_arr[cell] if layer == 0 else np.asarray(botm, dtype=float)[layer - 1].reshape(-1)[cell])
        else:
            cell_top = None
        if botm is not None:
            cell_botm = float(np.asarray(botm, dtype=float)[layer].reshape(-1)[cell])
        else:
            cell_botm = None
        return cell_top, cell_botm
    except Exception:
        return None, None


def _validate_zero_based_cellid(
    report: SurfaceWaterValidationReport,
    *,
    model: "SimulationBase",
    cellid: Any,
    code_prefix: str,
    context: dict[str, Any],
) -> tuple[int | None, int | None]:
    """Validate a zero-based MF6 cellid tuple and return layer/cell when valid."""

    if not isinstance(cellid, tuple) or len(cellid) != 2:
        report.add_error(
            f"{code_prefix}_cellid_shape",
            "Cell ids must be zero-based `(layer, cell)` tuples.",
            cellid=cellid,
            **context,
        )
        return None, None

    layer = _coerce_int(cellid[0])
    cell = _coerce_int(cellid[1])
    if layer is None or cell is None:
        report.add_error(
            f"{code_prefix}_cellid_type",
            "Cell ids must contain integer-like zero-based layer and cell values.",
            cellid=cellid,
            **context,
        )
        return None, None

    if layer < 0 or cell < 0:
        report.add_error(
            f"{code_prefix}_cellid_negative",
            "Cell ids must stay zero-based and non-negative.",
            cellid=cellid,
            **context,
        )
        return None, None

    nlay = _model_nlay(model)
    ncpl = _model_ncpl(model)
    if nlay is not None and layer >= nlay:
        report.add_error(
            f"{code_prefix}_layer_range",
            "Layer index is outside the model layer range.",
            cellid=cellid,
            nlay=nlay,
            **context,
        )
    if ncpl is not None and cell >= ncpl:
        report.add_error(
            f"{code_prefix}_cell_range",
            "Cell index is outside the model grid range.",
            cellid=cellid,
            ncpl=ncpl,
            **context,
        )
    return layer, cell


def _infer_mvr_package_count(model: "SimulationBase", package_name: str) -> int | None:
    """Infer how many mover-addressable entries a package exposes."""

    gwf = getattr(model, "gwf", None)
    if gwf is None:
        return None
    package = getattr(gwf, package_name, None)
    if package is None:
        return None

    if hasattr(package, "packagedata") and getattr(package.packagedata, "array", None) is not None:
        try:
            return int(len(package.packagedata.array))
        except Exception:
            return None
    return None


def validate_lak_configuration(
    model: "SimulationBase",
    *,
    nlakes: int,
    packagedata: list | None,
    connectiondata: list | None,
    perioddata: dict | None = None,
) -> SurfaceWaterValidationReport:
    """Validate LAK package inputs against the current model grid.

    Parameters
    ----------
    model
        Model providing grid and surface information.
    nlakes
        Number of lakes declared for the LAK package.
    packagedata
        LAK packagedata rows.
    connectiondata
        LAK connectiondata rows.
    perioddata
        Optional LAK period data dictionary.
    """

    report = SurfaceWaterValidationReport()
    rows = _normalize_lak_connection_rows(connectiondata)
    if nlakes < 1:
        report.add_error("lak_nlakes", "LAK must declare at least one lake.", nlakes=nlakes)
    if not rows:
        report.add_error("lak_connectiondata_empty", "LAK connectiondata is empty.")

    connections_per_lake: dict[int, int] = {}
    connected_cells_by_lake: dict[int, list[tuple[int, int]]] = {}
    seen_conn_ids: set[tuple[int, int]] = set()

    for row_index, row in enumerate(rows):
        if not isinstance(row, list | tuple) or len(row) < 9:
            report.add_error(
                "lak_connection_row_shape",
                "Each LAK connection row must include at least 9 values.",
                row_index=row_index,
                row=row,
            )
            continue

        lake_no = _coerce_int(row[0])
        conn_id = _coerce_int(row[1])
        cellid = row[2]
        conn_type = str(row[3]).upper()
        leakance = _coerce_float(row[4])
        belev = _coerce_float(row[5])
        telev = _coerce_float(row[6])
        conn_len = _coerce_float(row[7])
        conn_width = _coerce_float(row[8])

        context = {"row_index": row_index, "lake_no": lake_no, "conn_id": conn_id}
        if lake_no is None or not (0 <= lake_no < nlakes):
            report.add_error(
                "lak_lake_index",
                "LAK connection rows must reference valid zero-based lake ids.",
                row_index=row_index,
                lake_no=row[0],
                nlakes=nlakes,
            )
            continue

        if conn_id is None or conn_id < 0:
            report.add_error(
                "lak_connection_index",
                "LAK connection rows must use non-negative connection indices.",
                row_index=row_index,
                conn_id=row[1],
                lake_no=lake_no,
            )
        elif (lake_no, conn_id) in seen_conn_ids:
            report.add_error(
                "lak_connection_duplicate",
                "LAK connection indices must be unique within each lake.",
                row_index=row_index,
                lake_no=lake_no,
                conn_id=conn_id,
            )
        else:
            seen_conn_ids.add((lake_no, conn_id))

        layer, cell = _validate_zero_based_cellid(
            report,
            model=model,
            cellid=cellid,
            code_prefix="lak",
            context=context,
        )
        if layer is not None and cell is not None:
            connected_cells_by_lake.setdefault(lake_no, []).append((layer, cell))

        if conn_type not in _VALID_LAK_CONNECTION_TYPES:
            report.add_warning(
                "lak_connection_type",
                "LAK connection type is not one of the expected MF6 connection types.",
                conn_type=conn_type,
                **context,
            )
        if leakance is None or leakance <= 0:
            report.add_error(
                "lak_leakance",
                "Lake-bed leakance must be positive.",
                leakance=row[4],
                **context,
            )
        if belev is not None and telev is not None and belev > telev:
            report.add_error(
                "lak_elevation_order",
                "LAK connection bottom elevation cannot exceed the top elevation.",
                belev=belev,
                telev=telev,
                **context,
            )
        if conn_type == "HORIZONTAL":
            if conn_len is None or conn_len <= 0:
                report.add_error(
                    "lak_horizontal_length",
                    "Horizontal LAK connections must have positive connection length.",
                    conn_len=row[7],
                    **context,
                )
            if conn_width is None or conn_width <= 0:
                report.add_error(
                    "lak_horizontal_width",
                    "Horizontal LAK connections must have positive connection width.",
                    conn_width=row[8],
                    **context,
                )

        connections_per_lake[lake_no] = connections_per_lake.get(lake_no, 0) + 1

    if not packagedata or len(packagedata) != nlakes:
        report.add_error(
            "lak_packagedata_count",
            "LAK packagedata row count must match `nlakes`.",
            nlakes=nlakes,
            num_rows=0 if packagedata is None else len(packagedata),
        )
    else:
        seen_lake_ids: set[int] = set()
        for row_index, row in enumerate(packagedata):
            if not isinstance(row, list | tuple) or len(row) < 3:
                report.add_error(
                    "lak_packagedata_shape",
                    "Each LAK packagedata row must contain at least lake id, stage, and connection count.",
                    row_index=row_index,
                    row=row,
                )
                continue

            lake_no = _coerce_int(row[0])
            if lake_no is None or not (0 <= lake_no < nlakes):
                report.add_error(
                    "lak_packagedata_lake_index",
                    "LAK packagedata rows must reference valid zero-based lake ids.",
                    row_index=row_index,
                    lake_no=row[0],
                    nlakes=nlakes,
                )
                continue
            if lake_no in seen_lake_ids:
                report.add_error(
                    "lak_packagedata_duplicate_lake",
                    "LAK packagedata rows must use unique lake ids.",
                    row_index=row_index,
                    lake_no=lake_no,
                )
            seen_lake_ids.add(lake_no)

            declared_connections = _coerce_int(row[2])
            actual_connections = connections_per_lake.get(lake_no, 0)
            if declared_connections != actual_connections:
                report.add_error(
                    "lak_packagedata_connection_count",
                    "LAK packagedata connection count does not match the provided connectiondata.",
                    lake_no=lake_no,
                    declared_connections=declared_connections,
                    actual_connections=actual_connections,
                )

    if perioddata is not None:
        nper = getattr(model, "nper", None)
        for per, settings in perioddata.items():
            per_idx = _coerce_int(per)
            if per_idx is None or per_idx < 0 or (nper is not None and per_idx >= nper):
                report.add_error(
                    "lak_period_index",
                    "LAK perioddata keys must be valid zero-based stress periods.",
                    period=per,
                    nper=nper,
                )
                continue
            for setting_index, setting in enumerate(settings):
                if not isinstance(setting, list | tuple) or len(setting) < 3:
                    report.add_error(
                        "lak_perioddata_shape",
                        "Each LAK perioddata row must include lake id plus at least one setting pair.",
                        period=per,
                        setting_index=setting_index,
                        setting=setting,
                    )
                    continue
                normalized_setting = list(setting)
                if isinstance(normalized_setting[0], str) and _coerce_int(normalized_setting[1]) is not None:
                    normalized_setting = [normalized_setting[1], normalized_setting[0], *normalized_setting[2:]]

                lake_no = _coerce_int(normalized_setting[0])
                if lake_no is None or not (0 <= lake_no < nlakes):
                    report.add_error(
                        "lak_perioddata_lake_index",
                        "LAK perioddata rows must reference valid zero-based lake ids.",
                        period=per,
                        setting_index=setting_index,
                        lake_no=setting[0],
                        nlakes=nlakes,
                    )
                    continue
                normalized_keys = [str(item).lower() for item in normalized_setting[1::2]]
                normalized_values = normalized_setting[2::2]
                if "stage" in normalized_keys and connected_cells_by_lake.get(lake_no):
                    stage_idx = normalized_keys.index("stage")
                    stage_value = _coerce_float(normalized_values[stage_idx])
                    tops: list[float] = []
                    bottoms: list[float] = []
                    for layer, cell in connected_cells_by_lake[lake_no]:
                        cell_top, cell_botm = _cell_top_bottom(model, layer, cell)
                        if cell_top is not None:
                            tops.append(cell_top)
                        if cell_botm is not None:
                            bottoms.append(cell_botm)
                    if stage_value is not None and tops and stage_value > max(tops):
                        report.add_warning(
                            "lak_stage_above_top",
                            "Lake stage is above the maximum connected-cell top elevation.",
                            period=per,
                            lake_no=lake_no,
                            stage=stage_value,
                            max_top=max(tops),
                        )
                    if stage_value is not None and bottoms and stage_value < min(bottoms):
                        report.add_warning(
                            "lak_stage_below_bottom",
                            "Lake stage is below the minimum connected-cell bottom elevation.",
                            period=per,
                            lake_no=lake_no,
                            stage=stage_value,
                            min_bottom=min(bottoms),
                        )

    return report


def validate_sfr_configuration(sfr: "SFR") -> SurfaceWaterValidationReport:
    """Validate an SFR builder's derived reach and connection data."""

    report = SurfaceWaterValidationReport()
    model = sfr.model
    packagedata = sfr.packagedata
    connectiondata = sfr.connectiondata
    perioddata = sfr.perioddata

    if sfr.total_nreaches < 1:
        report.add_error("sfr_no_reaches", "SFR requires at least one reach.")
        return report

    if len(packagedata) != sfr.total_nreaches:
        report.add_error(
            "sfr_packagedata_count",
            "SFR packagedata row count must match `total_nreaches`.",
            total_nreaches=sfr.total_nreaches,
            num_rows=len(packagedata),
        )
    if len(connectiondata) != sfr.total_nreaches:
        report.add_error(
            "sfr_connectiondata_count",
            "SFR connectiondata row count must match `total_nreaches`.",
            total_nreaches=sfr.total_nreaches,
            num_rows=len(connectiondata),
        )

    seen_reaches: set[int] = set()
    seen_cells: dict[tuple[int, int], int] = {}
    for row_index, row in enumerate(packagedata):
        if not isinstance(row, list | tuple) or len(row) < 10:
            report.add_error(
                "sfr_packagedata_shape",
                "Each SFR packagedata row must contain the standard MF6 reach fields.",
                row_index=row_index,
                row=row,
            )
            continue

        rno = _coerce_int(row[0])
        cellid = row[1]
        rlen = _coerce_float(row[2])
        rwid = _coerce_float(row[3])
        rgrd = _coerce_float(row[4])
        rtp = _coerce_float(row[5])
        rbth = _coerce_float(row[6])
        rhk = _coerce_float(row[7])
        man = _coerce_float(row[8])
        ncon = _coerce_int(row[9])

        if rno is None or rno < 0:
            report.add_error(
                "sfr_reach_index",
                "SFR packagedata rows must use non-negative zero-based reach numbers.",
                row_index=row_index,
                rno=row[0],
            )
            continue
        if rno in seen_reaches:
            report.add_error(
                "sfr_reach_duplicate",
                "SFR packagedata rows must use unique reach numbers.",
                row_index=row_index,
                rno=rno,
            )
        seen_reaches.add(rno)

        layer, cell = _validate_zero_based_cellid(
            report,
            model=model,
            cellid=cellid,
            code_prefix="sfr",
            context={"row_index": row_index, "rno": rno},
        )
        if layer is not None and cell is not None:
            if (layer, cell) in seen_cells:
                report.add_warning(
                    "sfr_duplicate_cell",
                    "Multiple SFR reaches map to the same grid cell.",
                    row_index=row_index,
                    rno=rno,
                    previous_rno=seen_cells[(layer, cell)],
                    cellid=cellid,
                )
            else:
                seen_cells[(layer, cell)] = rno

            cell_top, cell_botm = _cell_top_bottom(model, layer, cell)
            if rtp is not None and cell_top is not None and rtp > cell_top:
                report.add_warning(
                    "sfr_top_above_cell_top",
                    "SFR reach top is above the host-cell top elevation.",
                    rno=rno,
                    reach_top=rtp,
                    cell_top=cell_top,
                    cellid=cellid,
                )
            if rtp is not None and rbth is not None and cell_botm is not None and (rtp - rbth) < cell_botm:
                report.add_warning(
                    "sfr_bottom_below_cell_bottom",
                    "SFR streambed bottom is below the host-cell bottom elevation.",
                    rno=rno,
                    streambed_bottom=rtp - rbth,
                    cell_bottom=cell_botm,
                    cellid=cellid,
                )

        if rlen is None or rlen <= 0:
            report.add_error("sfr_reach_length", "SFR reach lengths must be positive.", rno=rno, rlen=row[2])
        if rwid is None or rwid <= 0:
            report.add_error("sfr_reach_width", "SFR reach widths must be positive.", rno=rno, rwid=row[3])
        if rgrd is None or rgrd < 0:
            report.add_error("sfr_reach_gradient", "SFR reach gradients must be non-negative.", rno=rno, rgrd=row[4])
        if rbth is None or rbth <= 0:
            report.add_error(
                "sfr_streambed_thickness",
                "SFR streambed thickness must be positive.",
                rno=rno,
                rbth=row[6],
            )
        if rhk is None or rhk < 0:
            report.add_error(
                "sfr_streambed_k",
                "SFR streambed hydraulic conductivity must be non-negative.",
                rno=rno,
                rhk=row[7],
            )
        if man is None or man <= 0:
            report.add_error(
                "sfr_mannings",
                "SFR Manning's n values must be positive.",
                rno=rno,
                man=row[8],
            )

        if row_index < len(connectiondata):
            conn_row = connectiondata[row_index]
            if not isinstance(conn_row, list | tuple) or len(conn_row) < 1:
                report.add_error(
                    "sfr_connection_row_shape",
                    "Each SFR connection row must contain a reach id followed by connection ids.",
                    row_index=row_index,
                    connection_row=conn_row,
                )
            else:
                conn_rno = _coerce_int(conn_row[0])
                if conn_rno != rno:
                    report.add_error(
                        "sfr_connection_row_reach",
                        "SFR connection rows must start with the matching reach number.",
                        row_index=row_index,
                        rno=rno,
                        connection_rno=conn_row[0],
                    )
                actual_ncon = len(conn_row) - 1
                if ncon != actual_ncon:
                    report.add_error(
                        "sfr_connection_count",
                        "SFR packagedata `ncon` must match the number of listed connections.",
                        rno=rno,
                        declared_ncon=ncon,
                        actual_ncon=actual_ncon,
                    )

    expected_reaches = set(range(sfr.total_nreaches))
    if seen_reaches != expected_reaches:
        missing = sorted(expected_reaches - seen_reaches)
        if missing:
            report.add_error(
                "sfr_missing_reaches",
                "SFR packagedata reach numbering is not contiguous from 0 to `total_nreaches - 1`.",
                missing_reaches=missing[:10],
                num_missing=len(missing),
            )

    nper = getattr(model, "nper", None)
    for per, settings in perioddata.items():
        per_idx = _coerce_int(per)
        if per_idx is None or per_idx < 0 or (nper is not None and per_idx >= nper):
            report.add_error(
                "sfr_period_index",
                "SFR perioddata keys must be valid zero-based stress periods.",
                period=per,
                nper=nper,
            )
            continue
        for setting_index, setting in enumerate(settings):
            if not isinstance(setting, list | tuple) or len(setting) < 3:
                report.add_error(
                    "sfr_perioddata_shape",
                    "Each SFR perioddata row must include a reach id and setting/value fields.",
                    period=per,
                    setting_index=setting_index,
                    setting=setting,
                )
                continue
            reach = _coerce_int(setting[0])
            keyword = str(setting[1]).upper()
            if reach is None or not (0 <= reach < sfr.total_nreaches):
                report.add_error(
                    "sfr_perioddata_reach",
                    "SFR perioddata rows must reference valid zero-based reach ids.",
                    period=per,
                    setting_index=setting_index,
                    reach=setting[0],
                    total_nreaches=sfr.total_nreaches,
                )
            if keyword == "INFLOW":
                inflow = _coerce_float(setting[2])
                if inflow is not None and inflow < 0:
                    report.add_warning(
                        "sfr_negative_inflow",
                        "SFR inflow is negative.",
                        period=per,
                        reach=reach,
                        inflow=inflow,
                    )

    return report


def validate_mvr_configuration(
    model: "SimulationBase",
    *,
    maxmvr: int,
    maxpackages: int,
    packages: list,
    perioddata: dict,
) -> SurfaceWaterValidationReport:
    """Validate MVR package definitions and zero-based mover references."""

    report = SurfaceWaterValidationReport()
    if not isinstance(packages, list) or len(packages) == 0:
        report.add_error("mvr_packages_empty", "MVR `packages` must list at least one source/target package.")
        return report

    normalized_packages: list[str] = []
    for package_index, package_row in enumerate(packages):
        if not isinstance(package_row, list | tuple) or len(package_row) < 1:
            report.add_error(
                "mvr_package_row_shape",
                "Each MVR package row must include at least the package name.",
                package_index=package_index,
                package_row=package_row,
            )
            continue
        package_name = str(package_row[0]).lower()
        normalized_packages.append(package_name)

    if len(normalized_packages) > maxpackages:
        report.add_error(
            "mvr_maxpackages",
            "MVR `maxpackages` is smaller than the number of declared package rows.",
            maxpackages=maxpackages,
            declared_packages=len(normalized_packages),
        )

    gwf_package_names = {
        str(name).lower()
        for name in getattr(getattr(model, "gwf", None), "package_name_dict", {}).keys()
    }

    max_records_in_period = 0
    for per, records in perioddata.items():
        max_records_in_period = max(max_records_in_period, len(records))
        for record_index, record in enumerate(records):
            if not isinstance(record, list | tuple) or len(record) != 6:
                report.add_error(
                    "mvr_perioddata_shape",
                    "Each MVR perioddata row must contain 6 values: source package/id, target package/id, type, and value.",
                    period=per,
                    record_index=record_index,
                    record=record,
                )
                continue

            src_pkg = str(record[0]).lower()
            src_id = _coerce_int(record[1])
            dst_pkg = str(record[2]).lower()
            dst_id = _coerce_int(record[3])
            mover_type = str(record[4]).upper()
            mover_value = _coerce_float(record[5])

            if src_pkg not in normalized_packages or dst_pkg not in normalized_packages:
                report.add_error(
                    "mvr_package_reference",
                    "MVR perioddata must reference packages declared in the `packages` block.",
                    period=per,
                    record_index=record_index,
                    source_package=src_pkg,
                    target_package=dst_pkg,
                    declared_packages=normalized_packages,
                )

            if src_pkg not in gwf_package_names:
                report.add_error(
                    "mvr_missing_source_package",
                    "Source MVR package is not attached to the model.",
                    period=per,
                    record_index=record_index,
                    source_package=src_pkg,
                )
            if dst_pkg not in gwf_package_names:
                report.add_error(
                    "mvr_missing_target_package",
                    "Target MVR package is not attached to the model.",
                    period=per,
                    record_index=record_index,
                    target_package=dst_pkg,
                )

            if mover_type not in _VALID_MVR_TYPES:
                report.add_error(
                    "mvr_type",
                    "MVR transfer type must be one of the standard MF6 mover types.",
                    period=per,
                    record_index=record_index,
                    mover_type=mover_type,
                )
            if src_id is None or src_id < 0:
                report.add_error(
                    "mvr_source_id",
                    "MVR source ids must be zero-based and non-negative.",
                    period=per,
                    record_index=record_index,
                    source_id=record[1],
                )
            if dst_id is None or dst_id < 0:
                report.add_error(
                    "mvr_target_id",
                    "MVR target ids must be zero-based and non-negative.",
                    period=per,
                    record_index=record_index,
                    target_id=record[3],
                )

            source_count = _infer_mvr_package_count(model, src_pkg)
            target_count = _infer_mvr_package_count(model, dst_pkg)
            if source_count is not None and src_id is not None and src_id >= source_count:
                report.add_error(
                    "mvr_source_range",
                    "MVR source id is outside the available package index range.",
                    period=per,
                    record_index=record_index,
                    source_package=src_pkg,
                    source_id=src_id,
                    source_count=source_count,
                )
            if target_count is not None and dst_id is not None and dst_id >= target_count:
                report.add_error(
                    "mvr_target_range",
                    "MVR target id is outside the available package index range.",
                    period=per,
                    record_index=record_index,
                    target_package=dst_pkg,
                    target_id=dst_id,
                    target_count=target_count,
                )

            if mover_value is None:
                report.add_error(
                    "mvr_value_type",
                    "MVR mover values must be numeric.",
                    period=per,
                    record_index=record_index,
                    mover_value=record[5],
                )
            elif mover_value < 0:
                report.add_error(
                    "mvr_value_negative",
                    "MVR mover values must be non-negative.",
                    period=per,
                    record_index=record_index,
                    mover_value=mover_value,
                )
            elif mover_type == "FACTOR" and mover_value > 1.0:
                report.add_warning(
                    "mvr_factor_gt_one",
                    "MVR FACTOR transfers are usually between 0 and 1.",
                    period=per,
                    record_index=record_index,
                    mover_value=mover_value,
                )

    if max_records_in_period > maxmvr:
        report.add_error(
            "mvr_maxmvr",
            "MVR `maxmvr` is smaller than the maximum number of mover records in one stress period.",
            maxmvr=maxmvr,
            max_records_in_period=max_records_in_period,
        )

    return report


def validate_surface_water_configuration(
    model: "SimulationBase",
    *,
    nlakes: int | None = None,
    lak_packagedata: list | None = None,
    lak_connectiondata: list | None = None,
    lak_perioddata: dict | None = None,
    sfr: "SFR" | None = None,
    maxmvr: int | None = None,
    maxpackages: int | None = None,
    mvr_packages: list | None = None,
    mvr_perioddata: dict | None = None,
) -> SurfaceWaterValidationReport:
    """Validate a coupled ``LAK``/``SFR``/``MVR`` configuration together.

    Parameters
    ----------
    model
        Model providing grid, layer, and attached-package context.
    nlakes, lak_packagedata, lak_connectiondata, lak_perioddata
        Optional LAK inputs to validate.
    sfr
        Optional :class:`~simple_modflow.modflow.mf6.sfr.SFR` builder instance.
    maxmvr, maxpackages, mvr_packages, mvr_perioddata
        Optional MVR inputs to validate.
    """

    report = SurfaceWaterValidationReport()
    if nlakes is not None or lak_packagedata is not None or lak_connectiondata is not None or lak_perioddata is not None:
        report.extend(
            validate_lak_configuration(
                model,
                nlakes=1 if nlakes is None else nlakes,
                packagedata=lak_packagedata,
                connectiondata=lak_connectiondata,
                perioddata=lak_perioddata,
            )
        )
    if sfr is not None:
        report.extend(validate_sfr_configuration(sfr))
    if mvr_packages is not None or mvr_perioddata is not None:
        report.extend(
            validate_mvr_configuration(
                model,
                maxmvr=0 if maxmvr is None else maxmvr,
                maxpackages=0 if maxpackages is None else maxpackages,
                packages=[] if mvr_packages is None else mvr_packages,
                perioddata={} if mvr_perioddata is None else mvr_perioddata,
            )
        )
    return report
