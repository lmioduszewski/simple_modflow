"""Head observation targets (`HeadTargets` + bound views/plots)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, ClassVar

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point

from myflopy._logging import get_logger
from myflopy.grid_identity import excavation_cells
from myflopy.modflow.mf6.grid.domain import model_domain
from myflopy.modflow.mf6.grid.geometry import cell_shares, domain_active, miss_reason
from myflopy.modflow.mf6.observations._shared import (
    _build_values_from_cells,
    _coerce_label_list,
    _coerce_scalar_or_list,
    _copy_frame,
    _default_obs_csv_name,
    _ensure_unique_observation_names,
    _in_grid_crs,
    _load_locations,
    _normalize_identifier_series,
    _normalize_values_frame,
    _placement_note,
    _report_missed_features,
    _simulated_field_by_period,
)
from myflopy.viz import mpl_axes

logger = get_logger(__name__)


def _dug_cells(model) -> dict[int, str]:
    """``{cell: excavation name}`` for whatever records this model's context or grid carry.

    Only used to NAME why a target's cell is not active (ledger 321): a cell inside an
    excavation is lake water, and the head there belongs to the cell under the bed.
    """

    import types

    found: dict[int, str] = {}
    for context in (
        getattr(getattr(model, "gwf", None), "myflopy_context", None),
        getattr(model, "myflopy_context", None),
        types.SimpleNamespace(grid=getattr(model, "vor", None)),
    ):
        for cell, name in excavation_cells(context).items():
            found.setdefault(cell, name)
    return found


_DUG_ADVICE = (
    " A cell inside an excavation is lake water, not aquifer -- observe the head in the "
    "cell under its bed, the column's top active layer."
)


@dataclass
class HeadTargets:
    """Observation locations plus target heads through time or stress period.

    Parameters
    ----------
    locations
        Path to a vector dataset or an in-memory ``GeoDataFrame`` describing
        observation locations.
    values
        CSV path or DataFrame containing target heads. Both long and wide forms
        are supported.
    name_column
        Name of the observation identifier column in ``locations``.
    layer_column
        Optional layer column in ``locations``. Defaults to layer 0 when absent.
    group_column
        Optional observation-group column in ``locations``.
    weight_column
        Optional observation-weight column in ``locations``.
    time_column
        Column in ``values`` representing stress period, time stamp, or another
        row identifier shared with simulated results.
    value_column
        Value column name used when ``values`` is supplied in long form.

    Notes
    -----
    A target keyed by stress period reads the value at that period's END (ledger
    336) -- the record that ends the period on the head (or concentration) file's
    own clock, the one rule the stage and flow targets read too. A period whose
    final step output control did not save (``SAVE HEAD FIRST``, ``STEPS`` short
    of it) is refused: its end-of-period value does not exist in the file, and it
    used to be read at the last step saved, wherever that fell.
    """

    #: model attribute holding the full value table for this field
    _TABLE_ATTRIBUTE: ClassVar[str] = "all_heads"
    #: the value column inside that table
    _STORE_COLUMN: ClassVar[str] = "elev"
    #: human-readable field name, used in errors and plot labels
    _FIELD_LABEL: ClassVar[str] = "head"

    locations: str | Path | gpd.GeoDataFrame | pd.DataFrame | pd.Series | dict | list | tuple
    values: str | Path | pd.DataFrame | pd.Series | dict | list | tuple
    name_column: str = "name"
    layer_column: str | None = "layer"
    group_column: str | None = "group"
    weight_column: str | None = "weight"
    time_column: str = "time"
    value_column: str = "head"
    times: list | tuple | pd.Index | pd.Series | None = None

    def __post_init__(self):
        """Load and normalize locations and target tables."""

        self._locations = _load_locations(
            self.locations,
            name_column=self.name_column,
            layer_column=self.layer_column,
            group_column=self.group_column,
            weight_column=self.weight_column,
        )
        self._values = _normalize_values_frame(
            self.values,
            time_column=self.time_column,
            value_column=self.value_column,
            times=self.times,
        )
        if "name" not in self._values.columns:
            if len(self._locations["name_key"].unique()) != 1:
                raise ValueError(
                    "HeadTargets values without explicit names require exactly one target location."
                )
            only_name = str(self._locations["name"].iloc[0])
            self._values["name"] = only_name
            self._values["name_key"] = _normalize_identifier_series(self._values["name"])
        missing = sorted(
            set(self._values["name_key"].unique()) - set(self._locations["name_key"].unique())
        )
        if missing:
            raise ValueError(
                "Head target values reference locations not found in the locations layer: "
                + ", ".join(missing[:10])
            )

    @classmethod
    def from_cells(
        cls,
        *,
        cells,
        values,
        times=None,
        names=None,
        layers=0,
        groups=None,
        weights=None,
        time_column: str = "time",
        value_column: str = "head",
    ):
        """Create head targets directly from explicit model cell ids."""

        cell_list = [int(value) for value in list(cells)]
        count = len(cell_list)
        if count == 0:
            raise ValueError("from_cells(...) requires at least one cell.")
        names_list = _coerce_label_list(names, count=count, default_prefix="OBS")
        layer_list = [int(value) for value in _coerce_scalar_or_list(layers, count=count, default=0)]
        group_list = _coerce_scalar_or_list(groups, count=count, default=pd.NA)
        weight_list = _coerce_scalar_or_list(weights, count=count, default=np.nan)
        locations = pd.DataFrame(
            {
                "name": names_list,
                "layer": layer_list,
                "cell": cell_list,
                "group": group_list,
                "weight": weight_list,
            }
        )
        values_frame = _build_values_from_cells(
            values,
            names=names_list,
            times=times,
            time_column=time_column,
            value_column=value_column,
        )
        return cls(
            locations=locations,
            values=values_frame,
            name_column="name",
            layer_column="layer",
            group_column="group",
            weight_column="weight",
            time_column=time_column,
            value_column=value_column,
        )

    @classmethod
    def from_records(
        cls,
        records,
        *,
        time_column: str = "time",
        value_column: str = "head",
    ):
        """Create head targets from long-format record dictionaries/dataframes."""

        frame = pd.DataFrame(records).copy()
        if frame.empty:
            raise ValueError("from_records(...) requires at least one record.")
        if "name" not in frame.columns:
            raise ValueError("from_records(...) requires a 'name' column.")
        if time_column not in frame.columns:
            raise ValueError(f"from_records(...) requires the time column {time_column!r}.")
        if value_column not in frame.columns:
            raise ValueError(f"from_records(...) requires the value column {value_column!r}.")
        if "cell" not in frame.columns and "geometry" not in frame.columns and "node" not in frame.columns:
            raise ValueError("from_records(...) requires 'cell', 'node', or 'geometry' columns.")
        locations_columns = [
            column
            for column in ["name", "layer", "cell", "node", "group", "weight", "geometry", "x", "y"]
            if column in frame.columns
        ]
        locations = frame.loc[:, locations_columns].drop_duplicates(subset=["name"]).copy()
        if "geometry" in locations.columns:
            locations = gpd.GeoDataFrame(locations, geometry="geometry")
        values_frame = frame.loc[:, [time_column, "name", value_column]].copy()
        return cls(
            locations=locations,
            values=values_frame,
            name_column="name",
            layer_column="layer",
            group_column="group",
            weight_column="weight",
            time_column=time_column,
            value_column=value_column,
        )

    @property
    def locations_gdf(self) -> gpd.GeoDataFrame:
        """Normalized observation locations table."""

        return self._locations.copy()

    def get(self) -> pd.DataFrame:
        """Return the merged long-format target table."""

        if isinstance(self._locations, gpd.GeoDataFrame):
            location_frame = self._locations.drop(columns=self._locations.geometry.name)
        else:
            location_frame = self._locations.copy()
        merged = self._values.merge(
            location_frame,
            on="name_key",
            how="left",
            suffixes=("", "_loc"),
        )
        if "name_loc" in merged.columns:
            merged = merged.drop(columns=["name"]).rename(columns={"name_loc": "name"})
        return _copy_frame(merged)

    def to_long(self) -> pd.DataFrame:
        """Return the long-format target table."""

        return self.get()

    def to_wide(self) -> pd.DataFrame:
        """Return the target table in wide format with one column per location."""

        frame = self.get().pivot_table(
            index="time",
            columns="name",
            values="head_target",
            aggfunc="first",
        ).reset_index()
        frame.columns.name = None
        return frame

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of the target dataset."""

        frame = self.get()
        summary = {
            "n_locations": int(frame["name_key"].nunique()),
            "n_rows": int(len(frame)),
            "n_groups": int(frame["group"].dropna().nunique()) if "group" in frame else 0,
            "n_weighted": int(frame["weight"].notna().sum()) if "weight" in frame else 0,
        }
        return pd.DataFrame([summary])

    def plot_locations(self, *, ax=None, **kwargs):
        """Plot target locations and return the matplotlib axis."""

        if not isinstance(self._locations, gpd.GeoDataFrame):
            raise ValueError("plot_locations() requires geospatial target locations.")
        if ax is None:
            _, ax = mpl_axes()
        self._locations.plot(ax=ax, **kwargs)
        ax.set_title("Head Target Locations")
        return ax

    def match_to_model(self, model) -> pd.DataFrame:
        """Map target locations to zero-based model cell ids.

        Parameters
        ----------
        model
            ``SimulationBase`` or another model-like object exposing
            ``vor.gdf_vorPolys``.

        Returns
        -------
        pandas.DataFrame
            One row per MATCHED target location, with its ``cell`` and ``layer``.

        Notes
        -----
        A point resolves to the ONE cell that holds it: on a shared face or
        vertex, the lowest-numbered cell active in the target's own layer
        (ledger 233). It used to reach every cell it touched -- two on a face,
        four on a vertex -- so one well became four observations, and PEST's
        build and its forward run read different cells for it.

        Targets are reprojected to the grid's CRS first. One that lands on no
        active cell is left out and named in a warning; if NONE of them land,
        that raises, rather than handing back an observation set that is quietly
        about somewhere else. A target given as a cell id is checked instead: an
        id this grid does not have, or one inactive in its layer, raises.
        """

        if "cell" in self._locations.columns:
            matches = self._validated_cells(model)
            if matches["x"].isna().any() or matches["y"].isna().any():
                vor = model.vor.gdf_vorPolys.reset_index().rename(columns={"index": "cell"})
                centroids = pd.DataFrame(
                    {
                        "cell": vor["cell"].astype(int),
                        "x": vor.geometry.centroid.x,
                        "y": vor.geometry.centroid.y,
                    }
                )
                matches = matches.drop(columns=["x", "y"], errors="ignore").merge(centroids, on="cell", how="left")
            return _copy_frame(
                matches.loc[:, ["name", "name_key", "layer", "group", "weight", "cell", "x", "y"]]
            )

        return self._matched_from_geometry(model)

    def _validated_cells(self, model) -> pd.DataFrame:
        """The locations, with their explicit cell ids checked against the grid.

        An id is an assertion about the model, so a wrong one raises rather than
        being warned about the way a missed point is. Before this, an id the grid
        did not have reached the value table as NaN and dropped out of ``stats``
        silently, and an inactive id was read where MF6 writes no value.
        """

        matches = self._locations.copy()
        ncpl = len(model.vor.gdf_vorPolys)
        domain = model_domain(model, ncpl)
        cells = matches["cell"].astype(int)
        layers = matches["layer"].astype(int)
        outside = [
            f"{name!r} (cell {cell})"
            for name, cell in zip(matches["name"], cells, strict=True)
            if not 0 <= cell < ncpl
        ]
        if outside:
            raise ValueError(
                f"These {self._FIELD_LABEL} targets name cells this grid does not "
                f"have (it has {ncpl}, numbered from 0): {', '.join(outside)}."
            )
        dug = _dug_cells(model)
        inactive = [
            f"{name!r} (cell {cell}, layer {layer}"
            + (f", inside excavation {dug[cell]!r}" if cell in dug else "")
            + ")"
            for name, cell, layer in zip(matches["name"], cells, layers, strict=True)
            if not domain_active(domain, layer, cell)
        ]
        if inactive:
            raise ValueError(
                f"These {self._FIELD_LABEL} targets name cells that are not active "
                f"in their own layer, where MF6 writes no {self._FIELD_LABEL}: "
                f"{', '.join(inactive)}."
                + (_DUG_ADVICE if any("inside excavation" in entry for entry in inactive) else "")
            )
        return matches

    def _matched_from_geometry(self, model) -> pd.DataFrame:
        """Resolve each target's point to the one cell it sits in."""

        polys = model.vor.gdf_vorPolys
        locations = self._locations_in_grid_crs(model)
        domain = model_domain(model, len(polys))
        columns = ["name", "name_key", "layer", "group", "weight", "cell", "x", "y"]
        rows: list[dict[str, Any]] = []
        missed: list[str] = []
        dug: dict[int, str] | None = None
        for _, row in locations.iterrows():
            geometry = row.geometry
            name = str(row["name"])
            if geometry is not None and not geometry.is_empty and not isinstance(geometry, Point):
                raise ValueError(
                    f"{self._FIELD_LABEL.capitalize()} target {name!r} is a "
                    f"{geometry.geom_type}. A {self._FIELD_LABEL} is measured at a "
                    f"POINT: give its coordinates, or the cell it sits in."
                )
            layer = int(row["layer"])
            hits = cell_shares(
                polys,
                geometry,
                # Eligibility only. ``cell_shares`` sorts whatever this returns,
                # so the lowest-numbered usable cell wins a tie without a second
                # copy of that rule here.
                usable=lambda cells, layer=layer: [
                    cell for cell in cells if domain_active(domain, layer, cell)
                ],
            )
            if not hits:
                reason = miss_reason(polys, geometry, target=f"active cell in layer {layer}")
                if dug is None:
                    dug = _dug_cells(model)
                inside = [cell for cell, _share in cell_shares(polys, geometry)] if dug else []
                if inside and inside[0] in dug:
                    reason += (
                        f", inside excavation {dug[inside[0]]!r} -- the head there is the "
                        "cell under its bed, in the column's top active layer"
                    )
                missed.append(f"{name!r} ({reason})")
                continue
            rows.append({column: row.get(column) for column in columns} | {"cell": hits[0][0]})
        if missed:
            self._report_missed(model, missed, total=len(locations))
        frame = pd.DataFrame(rows, columns=columns)
        frame["cell"] = frame["cell"].astype(int)
        frame["layer"] = frame["layer"].astype(int)
        return _copy_frame(frame)

    def _locations_in_grid_crs(self, model) -> gpd.GeoDataFrame:
        """The locations in the GRID's CRS, with x/y recomputed from them.

        Targets are read in whatever CRS their file carries, and the grid decides
        where a cell is. Before this, targets in degrees over a State Plane grid
        were located by their raw degrees: two wells whose cells were 5 and 10
        both resolved to cell 4.
        """

        locations = self._locations
        reprojected = _in_grid_crs(locations, model)
        if reprojected is locations:
            return locations
        reprojected["x"] = reprojected.geometry.x
        reprojected["y"] = reprojected.geometry.y
        return reprojected

    def _report_missed(self, model, missed: list[str], *, total: int) -> None:
        """Raise when no target lands on an active cell; warn, by name, when only some do.

        A partial miss can be deliberate -- one target file serving several grids
        -- so it warns. None of them landing means the targets and the model are
        in different places, which is what the extents and CRSs in the message
        are read for. Every consumer (``simulated_heads``, ``compare``, ``stats``,
        ``to_flopy_obs``) calls ``match_to_model`` again, so the warning is said
        once per set of missed targets rather than once per call.
        """

        _report_missed_features(
            self,
            missed,
            total=total,
            subject=f"{self._FIELD_LABEL} target",
            nothing_left="compare against",
            where=lambda: self._where(model),
            log=logger,
        )

    def _where(self, model) -> str:
        """Where the targets are and where the grid is -- what an empty match is acted on with."""

        return _placement_note(self._locations, model)

    def observation_labels(self, model) -> pd.DataFrame:
        """The matched targets, each with the label its observation is known by.

        ONE definition of that label, shared by the MF6 ``continuous`` records,
        the simulated CSV's columns, the mapping the forward run reads, and the
        pst rows. A name used once IS its label; a name used in two layers becomes
        ``<name>_lay<layer>`` (then ``_c<cell>``), which is what the FloPy records
        have always done.

        The PEST side keyed on the bare name instead, so a well screened in two
        layers was one column: the build pivoted ``first`` and read layer 0, the
        forward run wrote each mapping row in turn and left layer 1, and the two
        disagreed by whatever the layers differ by (5.0 against 1005.0 measured).
        """

        matches = self.match_to_model(model).copy()
        matches["obsname"] = _ensure_unique_observation_names(matches)
        return matches

    def simulated_heads(self, model) -> pd.DataFrame:
        """Return simulated heads for the target locations by stress period.

        The result is a wide table suitable for pyEMU list-style observation
        setup: one row per stress period and one column per observation LABEL
        (see :meth:`observation_labels`, which a repeated name disambiguates).
        Each period's value is the one at its END (ledger 336), read on the
        head file's own clock.

        Parameters
        ----------
        model : SimulationBase
            A run model holding the head (or, for ``ConcTargets``, the
            concentration) file, as ``mf.load_mf6_run(ws)`` or ``run.model()``
            returns.

        Returns
        -------
        pandas.DataFrame
            ``per`` plus one column per observation label.

        Raises
        ------
        ValueError
            When a period the file saved lacks its final step -- naming the
            period, the steps saved and the ``SAVE HEAD LAST`` setting that
            fixes it; or when no target lands on an active cell (ledger 211).
        """

        matches = self.observation_labels(model)
        sim = _simulated_field_by_period(
            model,
            table_attribute=self._TABLE_ATTRIBUTE,
            store_column=self._STORE_COLUMN,
        )
        merged = matches.merge(sim, on=["layer", "cell"], how="left")
        wide = merged.pivot_table(
            index="per",
            # The label, not the bare name: one column per (name, layer), so a
            # two-screen well is two observations rather than one column whose
            # value depended on which of its rows got there first.
            columns="obsname",
            values="sim_head",
            aggfunc="first",
        ).reset_index()
        wide.columns.name = None
        return wide

    def compare(self, model) -> pd.DataFrame:
        """Compare target heads against simulated heads by stress period.

        The first implementation assumes the target table's ``time_column`` can
        be interpreted directly as a zero-based stress period integer. This
        matches the initial Cumberland calibration use case where targets are
        defined per stress period. ``sim_head`` is the head at each period's END
        (ledger 336), as :meth:`simulated_heads` reads it.

        Parameters
        ----------
        model : SimulationBase
            A run model holding the head (or concentration) file.

        Returns
        -------
        pandas.DataFrame
            One row per target row: the target's columns plus ``per``,
            ``cell``, ``sim_head``, ``residual`` (simulated minus target) and
            ``abs_residual``.

        Raises
        ------
        ValueError
            When a target's ``time`` is not a zero-based stress-period integer,
            or a period the file saved lacks its final step (as
            :meth:`simulated_heads`).
        """

        targets = self.get().copy()
        per = pd.to_numeric(targets["time"], errors="coerce")
        if per.isna().any():
            raise ValueError(
                "HeadTargets.compare(...) currently expects the time column to "
                "contain zero-based stress period integers."
            )
        targets["per"] = per.astype(int)
        sim = _simulated_field_by_period(
            model,
            table_attribute=self._TABLE_ATTRIBUTE,
            store_column=self._STORE_COLUMN,
        )
        if "cell" in targets.columns and targets["cell"].notna().all():
            merged = targets.merge(sim, on=["per", "layer", "cell"], how="left")
            if (targets["x"].isna().any() or targets["y"].isna().any()) and hasattr(model, "vor"):
                matches = self.match_to_model(model).loc[:, ["name", "name_key", "layer", "cell", "x", "y"]]
                merged = merged.drop(columns=["x", "y"], errors="ignore").merge(
                    matches,
                    on=["name", "name_key", "layer", "cell"],
                    how="left",
                )
        else:
            matches = self.match_to_model(model).loc[:, ["name", "name_key", "layer", "cell", "x", "y"]]
            merged = (
                targets.merge(matches, on=["name", "name_key", "layer"], how="left")
                .merge(sim, on=["per", "layer", "cell"], how="left")
            )
        merged["residual"] = merged["sim_head"] - merged["head_target"]
        merged["abs_residual"] = merged["residual"].abs()
        return _copy_frame(merged)

    def residuals(self, model) -> pd.DataFrame:
        """Return a compact residual table."""

        frame = self.compare(model)
        return frame.loc[
            :,
            [
                "name",
                "group",
                "time",
                "per",
                "layer",
                "cell",
                "head_target",
                "sim_head",
                "residual",
                "abs_residual",
                "weight",
            ],
        ]

    def stats(self, model) -> pd.DataFrame:
        """Return aggregate residual statistics for the target set."""

        frame = self.compare(model).dropna(subset=["residual"])
        if frame.empty:
            return pd.DataFrame(
                [{"n": 0, "mean_error": np.nan, "mae": np.nan, "rmse": np.nan}]
            )
        residual = frame["residual"].to_numpy(dtype=float)
        stats = {
            "n": int(len(frame)),
            "mean_error": float(np.mean(residual)),
            "mae": float(np.mean(np.abs(residual))),
            "rmse": float(np.sqrt(np.mean(residual**2))),
        }
        return pd.DataFrame([stats])

    def to_flopy_obs(
        self,
        model,
        *,
        csv_name: str = "head_targets.csv",
        kind: str = "HEAD",
    ) -> dict[str, list[tuple[str, str, tuple[int, int]]]]:
        """Return a FloPy ``continuous`` dict for MF6 head observations."""

        matches = self.observation_labels(model)
        records = [
            (str(row.obsname), str(kind).upper(), (int(row.layer), int(row.cell)))
            for row in matches.itertuples(index=False)
        ]
        return {str(csv_name): records}

    def attach_flopy_obs(
        self,
        model,
        *,
        pname: str = "gwf_obs",
        filename: str = "head_targets.obs",
        csv_name: str | None = None,
        kind: str = "HEAD",
        package=None,
    ):
        """Attach an MF6 utility observation package for these head targets."""

        import flopy

        owner = model.gwf if package is None else package
        csv_name = _default_obs_csv_name(filename, "head_targets.csv") if csv_name is None else str(csv_name)
        continuous = self.to_flopy_obs(model, csv_name=csv_name, kind=kind)
        return flopy.mf6.modflow.mfutlobs.ModflowUtlobs(
            owner,
            pname=pname,
            continuous=continuous,
            filename=str(filename),
        )

    def calibration_plot(self, model, *, type: str = "calibration"):
        """Build a calibration plot directly from these targets and one model."""

        from myflopy.modflow.calcs.calibration import CalibrationPlot

        return CalibrationPlot.from_targets(self, model=model, type=type)


class BoundHeadTargets:
    """Model-bound head-target helper returned from ``model.targets.heads``."""

    def __init__(self, model, targets: HeadTargets):
        """Bind ``targets`` to ``model`` so its methods default to that model."""

        self.model = model
        self.targets = targets
        self._plot = None

    def get(self) -> pd.DataFrame:
        """The target definitions (delegates to the wrapped :class:`HeadTargets`)."""

        return self.targets.get()

    def to_long(self) -> pd.DataFrame:
        """The long-format target table (delegates to the wrapped targets)."""

        return self.targets.to_long()

    def to_wide(self) -> pd.DataFrame:
        """The wide-format target table (delegates to the wrapped targets)."""

        return self.targets.to_wide()

    def summary(self) -> pd.DataFrame:
        """A compact summary of the target set (delegates to the wrapped targets)."""

        return self.targets.summary()

    def plot_locations(self, *, ax=None, **kwargs):
        """Plot the target locations on a map axis (delegates to the wrapped targets)."""

        return self.targets.plot_locations(ax=ax, **kwargs)

    def match_to_model(self, model=None) -> pd.DataFrame:
        """Match each target location to its grid cell (defaults to the bound model)."""

        return self.targets.match_to_model(self.model if model is None else model)

    def simulated_heads(self, model=None) -> pd.DataFrame:
        """Simulated heads at the target locations (defaults to the bound model)."""

        return self.targets.simulated_heads(self.model if model is None else model)

    def compare(self, model=None) -> pd.DataFrame:
        """Compare observed vs simulated heads (defaults to the bound model)."""

        return self.targets.compare(self.model if model is None else model)

    def residuals(self, model=None) -> pd.DataFrame:
        """The residual table for the head targets (defaults to the bound model)."""

        return self.targets.residuals(self.model if model is None else model)

    def stats(self, model=None) -> pd.DataFrame:
        """Aggregate residual statistics (defaults to the bound model)."""

        return self.targets.stats(self.model if model is None else model)

    def to_flopy_obs(
        self,
        *,
        csv_name: str = "head_targets.csv",
        kind: str = "HEAD",
    ):
        """A FloPy ``continuous`` observation dict for the head targets on the bound model."""

        return self.targets.to_flopy_obs(self.model, csv_name=csv_name, kind=kind)

    def attach_flopy_obs(
        self,
        *,
        pname: str = "gwf_obs",
        filename: str = "head_targets.obs",
        csv_name: str | None = None,
        kind: str = "HEAD",
        package=None,
    ):
        """Attach an MF6 head observation package to the bound model."""

        return self.targets.attach_flopy_obs(
            self.model,
            pname=pname,
            filename=filename,
            csv_name=csv_name,
            kind=kind,
            package=package,
        )

    def calibration_plot(self, *, type: str = "calibration"):
        """A calibration plot for the head targets on the bound model."""

        return self.targets.calibration_plot(self.model, type=type)

    @property
    def plot(self):
        """Plotting helpers for this model-bound target set."""

        if self._plot is None:
            self._plot = BoundHeadTargetPlots(self)
        return self._plot


class BoundHeadTargetPlots:
    """Plotting facade for model-bound head targets."""

    def __init__(self, bound_targets: BoundHeadTargets):
        """Wrap a :class:`BoundHeadTargets` to expose its plotting helpers."""

        self.bound_targets = bound_targets

    @property
    def model(self):
        """The model the underlying targets are bound to."""

        return self.bound_targets.model

    @property
    def targets(self):
        """The underlying :class:`HeadTargets` being plotted."""

        return self.bound_targets.targets

    def locations(self, *, ax=None, **kwargs):
        """Plot target locations on a map axis."""

        return self.targets.plot_locations(ax=ax, **kwargs)

    def calibration(self, *, type: str = "calibration"):
        """Return a calibration/heads plot for the bound targets."""

        return self.bound_targets.calibration_plot(type=type)

    def obs_vs_sim(self, *, baseline=None, ax=None, backend: str = "plotly"):
        """Return a target-vs-simulated cross plot for one or two models.

        ``backend`` selects ``"plotly"`` (default) or ``"matplotlib"``.
        """

        from myflopy.modflow.calcs.calibration import CalibrationPlot

        current = self.targets.compare(self.model)
        baseline_frame = None if baseline is None else self.targets.compare(baseline)
        return CalibrationPlot.from_obs_vs_sim(current, baseline_compare=baseline_frame, backend=backend)

    def timeseries(self, name: str | None = None, *, baseline=None, ax=None, backend: str = "plotly"):
        """Return a time-series plot for one target location.

        ``backend`` selects ``"plotly"`` (default) or ``"matplotlib"``.
        """

        from myflopy.modflow.calcs.calibration import CalibrationPlot

        current = self.targets.compare(self.model)
        baseline_frame = None if baseline is None else self.targets.compare(baseline)
        return CalibrationPlot.from_timeseries(
            current,
            name=name,
            baseline_compare=baseline_frame,
            backend=backend,
        )

    def residuals_by_period(self, *, baseline=None, backend: str = "plotly"):
        """Return a by-period residual summary plot.

        ``backend`` selects ``"plotly"`` (default) or ``"matplotlib"``.
        """

        from myflopy.modflow.calcs.calibration import CalibrationPlot

        current = self.targets.compare(self.model)
        baseline_frame = None if baseline is None else self.targets.compare(baseline)
        return CalibrationPlot.from_residuals_by_period(
            current,
            baseline_compare=baseline_frame,
            backend=backend,
        )


