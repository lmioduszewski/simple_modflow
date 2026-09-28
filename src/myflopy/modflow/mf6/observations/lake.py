"""Lake stage observation targets (`LakeStageTargets` + bound view)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar

import numpy as np
import pandas as pd

from myflopy.modflow.mf6.observations._shared import (
    _aggregate_named_series,
    _copy_frame,
    _default_compare_stats,
    _default_obs_csv_name,
    _find_named_observation_output,
    _normalize_compare_time,
    _normalize_lake_stage_values,
    _numeric_periods,
    _observation_label,
)
from myflopy.modflow.mf6.observations.plots import BoundNamedSeriesTargetPlots
from myflopy.modflow.mf6.period_end import period_ends


def _lake_stage_table(model) -> pd.DataFrame:
    """``packages.lak.results.stage``, refused when a period's end was not saved (ledger 336).

    The table's ``per`` is the saved record's INDEX (ledger 301's stage clock),
    so it cannot be reduced to each period's end the way the SFR table is -- that
    stays the open limit ledger 331 pins. What it can do is refuse a schedule that
    saved no END for some period (``SAVE HEAD FIRST``, ``STEPS`` short of the end):
    there, one record per period made the index look like the period while every
    value was read before the period ended. Checked on the LAK stage file's own
    clock, which MF6 writes when it saves HEAD.

    Parameters
    ----------
    model : SimulationBase
        A run model with a LAK ``stage_filerecord``.

    Returns
    -------
    pandas.DataFrame
        ``packages.lak.results.stage.get()``, a copy.

    Raises
    ------
    ValueError
        When a period the LAK stage file saved lacks its final step -- naming the
        period, the steps saved and ``SAVE HEAD LAST``.
    """

    frame = model.packages.lak.results.stage.get().copy()
    # Called for the refusal alone: the record clock keeps every record.
    period_ends(
        model,
        reader=lambda: model.lak.output.stage(),
        file="the LAK stage file",
        saves="HEAD",
    )
    return frame


@dataclass
class LakeStageTargets:
    """Observed lake-stage time series tied to named lakes, for review or PEST.

    The lake-stage analogue of :class:`HeadTargets`: it pairs named lakes
    (``locations`` -- a ``name -> lake`` mapping, zero-based MF6 lake ids) with
    observed stage ``values`` through time, then aligns them with the model's
    simulated LAK stage (read from the LAK stage output, or a workspace
    observation CSV when present). Pass an instance directly to
    ``cal.observe(...)`` for calibration, or call ``.compare(model)`` /
    ``.stats(model)`` to review residuals without a PEST workflow. Omitting
    ``values`` yields an observation *template* useful for capturing simulated
    series only.

    A lake target on the stage file is keyed by SAVED RECORD, not by period
    (ledgers 301, 331) -- the one family the period-end rule cannot reduce yet --
    but a schedule that saves no period's END (``SAVE HEAD FIRST``, ``STEPS``
    short of it) is refused as it is for heads, stages and flows (ledger 336).

    Parameters
    ----------
    locations : dict, pandas.Series, pandas.DataFrame, list, tuple or None, optional
        The named lakes: a ``{name: lake}`` dict, a ``Series`` indexed by name,
        or a ``DataFrame`` / list of record dicts with ``name`` and ``lake``
        columns (``ValueError`` when either is missing). ``lake`` is MF6's lake
        number, ZERO-based (:meth:`to_flopy_obs` writes it one-based); each is
        passed through ``int`` and each name through ``str``. Several names may
        share a lake. Default ``None``: no targets.
    values : str, pathlib.Path, DataFrame, Series, dict, list, tuple or None, optional
        Observed lake stage through time, in the model's length units. A path is
        read with ``pandas.read_csv``. LONG form has ``name`` and
        ``value_column`` columns (headers matched ignoring case); WIDE form has
        one column per target name beside ``time_column``; a table of only
        ``time_column`` and ``value_column``, or a ``Series`` indexed by time,
        is one series and needs exactly one named lake (``ValueError``
        otherwise). Names are NOT checked against ``locations`` and are joined
        to them exactly (case counts): a row naming no target gets no ``lake``
        and no simulated stage in :meth:`compare`.
        Default ``None`` makes a template that only captures simulated series.
    time_column : str, optional
        Column of ``values`` holding each row's time label, matched exactly.
        Default ``"time"``; ``ValueError`` when ``values`` lacks it after
        ``times`` is applied. :meth:`compare` reads it as a zero-based integer
        -- see there for what it is matched against on each route.
    value_column : str, optional
        Column of a LONG ``values`` table holding the observed stage, matched
        ignoring case. Default ``"stage"``.
    times : list, tuple, pandas.Index, pandas.Series or None, optional
        Time labels inserted as ``time_column`` when ``values`` is a ``dict`` or
        list/tuple without that column; a list/tuple of another length raises
        ``ValueError``. Ignored for a ``DataFrame``, ``Series`` or CSV path.
        Default ``None``.

    Examples
    --------
    >>> LakeStageTargets(locations={"valley_lake": 0},
    ...                  values=observed_stage_df)        # -> cal.observe(...)
    >>> LakeStageTargets(locations={"valley_lake": 0}).compare(model)

    See Also
    --------
    HeadTargets, SfrStageTargets, SfrFlowTargets, DrnFlowTargets
    """

    locations: dict[str, int] | pd.Series | pd.DataFrame | list | tuple | None = None
    values: str | Path | pd.DataFrame | pd.Series | dict | list | tuple | None = None
    time_column: str = "time"
    value_column: str = "stage"
    times: list | tuple | pd.Index | pd.Series | None = None

    #: What a lake with no stage is WRITTEN as in the CSV PEST reads (ledger 329),
    #: by both writers -- the base run and ``write_named_series_targets``, which
    #: reads it off the instance because pyEMU copies that function's source and
    #: not its module. The series itself says ``NaN``; PEST++ fails the run on an
    #: empty field and on ``nan`` (``docs/reproducers/pest_no_value_fields.py``),
    #: so a lake dry in a period nobody even weights would fail every realization
    #: in which it is dry. MF6's own no-value magnitude (the stage file's inactive
    #: lake, the head file's inactive cell); the sign is not kept because the mask
    #: has already merged dry and inactive, and a squared residual has none.
    _PEST_NO_STAGE: ClassVar[float] = 1.0e30

    def __post_init__(self):
        """Normalize the ``name -> lake`` locations and observed values into internal tables."""

        normalized: dict[str, int] = {}
        source = self.locations if self.locations is not None else {}
        if isinstance(source, pd.Series):
            source = source.to_dict()
        elif isinstance(source, pd.DataFrame):
            if {"name", "lake"} - set(source.columns):
                raise ValueError("LakeStageTargets locations DataFrame requires 'name' and 'lake' columns.")
            source = source.set_index("name")["lake"].to_dict()
        elif isinstance(source, (list, tuple)):
            source_frame = pd.DataFrame(source)
            if {"name", "lake"} - set(source_frame.columns):
                raise ValueError("LakeStageTargets locations records require 'name' and 'lake' keys.")
            source = source_frame.set_index("name")["lake"].to_dict()
        for name, lake_id in dict(source).items():
            normalized[str(name)] = int(lake_id)
        self._series = normalized
        self._values = _normalize_lake_stage_values(
            self.values,
            time_column=self.time_column,
            value_column=self.value_column,
            times=self.times,
        )
        if not self._values.empty and "name" not in self._values.columns:
            if len(self._series) != 1:
                raise ValueError(
                    "LakeStageTargets values without explicit names require exactly one named lake target."
                )
            only_name = next(iter(self._series))
            self._values["name"] = only_name

    @classmethod
    def from_series(
        cls,
        *,
        lake: str,
        lake_id: int,
        values,
        times=None,
        time_column: str = "time",
        value_column: str = "stage",
    ):
        """Create one named lake-stage target from a Series/list-like input.

        Parameters
        ----------
        lake : str
            The target's name.
        lake_id : int
            MF6's lake number for it, ZERO-based (passed through ``int``).
        values : pandas.Series or sequence of float
            Observed stages, in the model's length units. A ``Series`` carries
            its time labels as its index; any other sized sequence is one stage
            per entry of ``times``.
        times : sequence or None, optional
            Time labels for a non-``Series`` ``values``, one per value (another
            length raises ``ValueError`` from pandas). Default ``None``:
            ``0, 1, 2, ...``. Ignored for a ``Series``.
        time_column : str, optional
            Name of the time column in the table built here, stored on the
            result. Default ``"time"``.
        value_column : str, optional
            Name of the stage column in the table built from a ``Series``,
            stored on the result. Default ``"stage"``.

        Returns
        -------
        LakeStageTargets
            A target set holding the one lake.
        """

        if isinstance(values, pd.Series):
            series = values.copy()
            series.index.name = time_column
            frame = pd.DataFrame(
                {
                    time_column: series.index.to_list(),
                    "name": [lake] * len(series),
                    value_column: series.to_list(),
                }
            )
            return cls(locations={lake: lake_id}, values=frame, time_column=time_column, value_column=value_column)
        if times is None:
            times = list(range(len(values)))
        frame = pd.DataFrame({time_column: list(times), lake: list(values)})
        return cls(locations={lake: lake_id}, values=frame, time_column=time_column, value_column=value_column)

    @classmethod
    def from_records(
        cls,
        records,
        *,
        time_column: str = "time",
        value_column: str = "stage",
    ):
        """Create lake-stage targets from long-format records.

        Parameters
        ----------
        records : list of dict or pandas.DataFrame
            One row per target per time, as ``pandas.DataFrame(records)`` reads
            it, with columns ``name``, ``lake`` (MF6's lake number, ZERO-based),
            ``time_column`` and ``value_column``, all matched exactly.
            ``ValueError`` when there are no records or a column is missing. A
            name's lake is taken from its FIRST row.
        time_column : str, optional
            Column holding each record's time label. Default ``"time"``.
        value_column : str, optional
            Column holding each record's observed stage. Default ``"stage"``.

        Returns
        -------
        LakeStageTargets
            The targets, one per distinct ``name``.
        """

        frame = pd.DataFrame(records).copy()
        if frame.empty:
            raise ValueError("from_records(...) requires at least one record.")
        if {"name", "lake", time_column, value_column} - set(frame.columns):
            missing = {"name", "lake", time_column, value_column} - set(frame.columns)
            raise ValueError("from_records(...) is missing required columns: " + ", ".join(sorted(missing)))
        locations = (
            frame.loc[:, ["name", "lake"]]
            .drop_duplicates(subset=["name"])
            .set_index("name")["lake"]
            .to_dict()
        )
        values_frame = frame.loc[:, [time_column, "name", value_column]].copy()
        return cls(locations=locations, values=values_frame, time_column=time_column, value_column=value_column)

    def get(self) -> pd.DataFrame:
        """Return the named lake-stage definitions or merged values."""

        definitions = pd.DataFrame(
            [{"name": name, "lake": lake_id} for name, lake_id in self._series.items()]
        )
        if self._values.empty:
            return definitions
        return _copy_frame(self._values.merge(definitions, on="name", how="left"))

    def to_long(self) -> pd.DataFrame:
        """Return the long-format lake-stage target table."""

        return self.get()

    def to_wide(self) -> pd.DataFrame:
        """Return the target table in wide format with one column per lake series."""

        if self._values.empty:
            return pd.DataFrame(columns=["time", *self._series.keys()])
        frame = self.get().pivot_table(
            index="time",
            columns="name",
            values="stage",
            aggfunc="first",
        ).reset_index()
        frame.columns.name = None
        return frame

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of the lake-stage target set."""

        return pd.DataFrame(
            [{"n_lakes": int(len(self._series)), "n_rows": int(len(self._values))}]
        )

    def _obs_column_map(self) -> dict[str, str]:
        """Map each target name to its MF6 observation label (for matching output CSVs)."""

        return {
            str(name): _observation_label(name)
            for name in self._series.keys()
        }

    def simulated_series(self, model) -> pd.DataFrame:
        """Return one wide simulated-stage table: a ``time`` column and one column per target.

        Read from the model's LAK stage observations when one workspace CSV holds
        every target's column (MF6's observation labels, as :meth:`to_flopy_obs`
        writes them), else from ``packages.lak.results.stage``. The two routes
        label ``time`` differently: the observation route carries MF6's simulated
        time, the stage-file route the zero-based index of the saved stage record.

        **That index is the stress period only when output control saves one step
        per period, and that step ENDS it** (ledgers 331, 336). Saving several --
        ``SAVE HEAD ALL`` at ``nstp=3`` -- writes three records per period, so a
        target keyed by stress period silently reads a step inside an EARLIER
        period, where the ``SfrStageTargets`` and ``HeadTargets`` beside it read
        each period's end. Key a lake target by saved-record index there. The open
        fix is ``build_lak_stage_result_table`` carrying the stage file's own
        ``kstpkper`` (ledger 301). A schedule that saves no period's end at all
        (``SAVE HEAD FIRST``, ``STEPS`` short of it) made the index look like the
        period with every value read early; that is refused now, as it is for
        heads, stages and flows.

        **A period with no stage is a row holding ``NaN``, never a missing row**
        (ledger 329). An INACTIVE lake is ``NaN`` on both routes: the observation
        file's 3e30 had been read as a stage of 0.0 ft, and the stage file's +1e30
        as itself. The routes differ at a DRY lake because MF6's two outputs do: its
        observation reports the stage it solved, which sits at the lake bottom
        (100 ft on the reproducer's lake, measured), while the stage file writes
        -1e30. So a dry period is the bottom on the observation route and ``NaN`` on
        the stage-file route.

        The stage-file route keeps one row per saved record and one column per
        target, even a target dry in every record, so the shape does not depend on
        which lakes are dry. That is what PEST needs: pyEMU builds the instruction
        file from the BASE run's copy of this table, one line advance per row and
        each value found by counting commas past a skipped ``time``, so a
        realization one row short fails at end of file and one with an extra
        column reads the wrong column without error (both measured,
        ``docs/reproducers/pest_no_value_fields.py``). The file PEST reads carries
        1e30 in place of ``NaN``, because PEST++ refuses an empty field and ``nan``
        alike (see ``write_named_series_targets``).

        **A lake reads its ONE stage on both routes** (ledger 331). The stage
        table repeats a lake's stage on every cell the lake connects to, and the
        stage-file route used to ADD those copies: the user's pond, 763 connected
        cells, read 401,339.97 ft for a stage of 526.00. It takes the one number
        per saved record now, as a stage is not added (ledger 297), and refuses
        copies that disagree rather than averaging them.

        Parameters
        ----------
        model : SimulationBase
            A run model with a LAK package, as ``mf.load_mf6_run(ws)`` or
            ``run.model()`` returns. The observation route needs its workspace;
            the stage-file route needs a LAK ``stage_filerecord``.

        Returns
        -------
        pandas.DataFrame
            ``time`` plus one float column per target name, ``NaN`` where the
            lake has no stage. On the stage-file route a target naming a lake the
            model does not have gets no column. Empty (those columns, no rows)
            when the stage file holds no record.

        Raises
        ------
        AttributeError
            When no workspace CSV holds the targets and there is no stage file
            to fall back on: the model has no LAK package ("Package 'lak' is not
            attached to model ...") or its LAK declares no ``stage_filerecord``.
        ValueError
            On the stage-file route: when a period the stage file saved lacks its
            final step -- naming the period, the steps saved and the ``SAVE HEAD
            LAST`` setting that fixes it (ledger 336); or when the rows of one
            lake in one saved record are not copies of one stage -- naming the
            lake, the record and the spread. MF6's stage table cannot produce the
            second; it guards a table that is not the one this reads.

        Examples
        --------
        >>> targets = mf.LakeStageTargets(locations={"pond": 0})
        >>> series = targets.simulated_series(model)   # one row per saved record
        >>> series["pond"].isna()                       # True where the pond is dry
        """

        observed = _find_named_observation_output(
            model, expected_columns=self._obs_column_map(), no_value=np.nan
        )
        if observed is not None:
            return observed
        frame = _lake_stage_table(model)
        if frame.empty:
            return pd.DataFrame(columns=["time", *self._series.keys()])
        frame = _aggregate_named_series(
            frame, id_column="lake", value_column="stage", kind="elevation"
        )
        definitions = pd.DataFrame([{"name": name, "lake": lake_id} for name, lake_id in self._series.items()])
        # A target naming a lake the model does not have matches no row; drop it
        # here, as the pivot always did, rather than as a NaN ``time`` row.
        merged = definitions.merge(frame, on="lake", how="left").dropna(subset=["per"])
        merged["per"] = merged["per"].astype(int)
        # ``dropna=False`` keeps a record in which every target is dry, and a
        # target dry in EVERY record -- the default drops both, and the second
        # would give a later realization one column more than PEST's base run.
        wide = merged.pivot_table(
            index="per", columns="name", values="stage", aggfunc="first", dropna=False
        ).reset_index()
        wide.columns.name = None
        return wide.rename(columns={"per": "time"})

    def compare(self, model) -> pd.DataFrame:
        """Compare target lake stages against simulated lake stages.

        ``sim_stage`` comes from the same two routes as :meth:`simulated_series`,
        and is the lake's ONE stage on both: the stage-file route used to add the
        copies the stage table holds for each connected cell (ledger 331). It is
        ``NaN`` where MF6 reports no stage, and so are ``residual`` and
        ``abs_residual``; which lakes that covers depends on the route (ledger
        329). Targets are matched on ``per``, read from their ``time`` as a
        zero-based integer -- which on the observation route is compared with
        MF6's simulated time, and on the stage-file route with the index of the
        saved stage record. That index is the stress period only at one saved step
        per period: with several, a target keyed by period silently reads a step
        inside an earlier period (measured at ``nstp=3``: periods 0, 1, 2 read
        steps 1, 2, 3 of period 0), unlike ``SfrStageTargets`` and
        ``HeadTargets``, which read each period's last step. Key it by
        saved-record index there until ``build_lak_stage_result_table`` carries
        the stage file's ``kstpkper`` (ledgers 331, 301). A schedule that saves no
        period's end is refused (ledger 336).

        Parameters
        ----------
        model : SimulationBase
            A run model with a LAK package, as ``mf.load_mf6_run(ws)`` or
            ``run.model()`` returns.

        Returns
        -------
        pandas.DataFrame
            One row per target row (per simulated row when the targets hold no
            values): ``name``, ``lake``, ``time``, ``per``, ``stage_target``,
            ``sim_stage``, ``residual`` (simulated minus target) and
            ``abs_residual``.

        Raises
        ------
        AttributeError
            As :meth:`simulated_series`: no workspace CSV holds the targets and
            the model has no LAK stage file.
        ValueError
            On the stage-file route, when a period the stage file saved lacks its
            final step or a lake's copies in one record disagree (as
            :meth:`simulated_series`), or when a target's ``time`` is not a
            number.
        """

        definitions = pd.DataFrame([{"name": name, "lake": lake_id} for name, lake_id in self._series.items()])
        observed = _find_named_observation_output(
            model, expected_columns=self._obs_column_map(), no_value=np.nan
        )
        if observed is not None:
            simulated = observed.melt(id_vars=["time"], var_name="name", value_name="sim_stage")
            simulated["per"] = simulated["time"]
            simulated = definitions.merge(simulated, on="name", how="left")
        else:
            simulated = _lake_stage_table(model)
            simulated = _aggregate_named_series(
                simulated, id_column="lake", value_column="stage", kind="elevation"
            )
            simulated = definitions.merge(simulated, on="lake", how="left").rename(columns={"stage": "sim_stage"})
        if self._values.empty:
            compare = simulated.rename(columns={"per": "time"}).copy()
            compare["stage_target"] = np.nan
        else:
            targets = self.get().copy()
            targets["per"] = _numeric_periods(targets, caller="LakeStageTargets.compare(...)")
            if "stage" in targets.columns and "stage_target" not in targets.columns:
                targets = targets.rename(columns={"stage": "stage_target"})
            compare = targets.merge(simulated, on=["name", "lake", "per"], how="left")
        compare = _normalize_compare_time(compare)
        compare["residual"] = compare["sim_stage"] - compare["stage_target"]
        compare["abs_residual"] = compare["residual"].abs()
        return _copy_frame(compare)

    def residuals(self, model) -> pd.DataFrame:
        """Return a compact residual table for lake-stage targets.

        Parameters
        ----------
        model : SimulationBase
            A run model with a LAK package; passed to :meth:`compare`, whose
            routes and ``Raises`` apply.

        Returns
        -------
        pandas.DataFrame
            :meth:`compare`'s rows, reduced to ``name``, ``lake``, ``time``,
            ``per``, ``stage_target``, ``sim_stage``, ``residual`` and
            ``abs_residual``.
        """

        frame = self.compare(model)
        return frame.loc[:, ["name", "lake", "time", "per", "stage_target", "sim_stage", "residual", "abs_residual"]]

    def stats(self, model) -> pd.DataFrame:
        """Return aggregate residual statistics for the lake-stage target set.

        Parameters
        ----------
        model : SimulationBase
            A run model with a LAK package; passed to :meth:`compare`, whose
            routes and ``Raises`` apply.

        Returns
        -------
        pandas.DataFrame
            One row: ``n`` (rows with a residual -- a lake with no stage has
            none), ``mean_error``, ``mae`` and ``rmse`` of simulated minus
            target stage, in model length units. ``n`` is 0 and the rest
            ``NaN`` when no row has a residual.
        """

        return _default_compare_stats(self.compare(model))

    def to_flopy_obs(
        self,
        *,
        csv_name: str = "lake_stage_targets.csv",
        kind: str = "STAGE",
    ) -> dict[str, list[tuple[str, str, int]]]:
        """Return a FloPy ``continuous`` dict for MF6 lake-stage observations.

        Parameters
        ----------
        csv_name : str, optional
            The CSV file MF6 writes these observations to -- the dict's one key.
            Default ``"lake_stage_targets.csv"``. :meth:`simulated_series` reads
            observations back from a workspace CSV holding every target's
            column, and skips any file whose name contains ``simulated_`` or
            ``_target_``.
        kind : str, optional
            The MF6 LAK observation type written on every record, upper-cased.
            Default ``"STAGE"``. The observation route of :meth:`simulated_series`
            and :meth:`compare` reads whatever that CSV holds as the stage.

        Returns
        -------
        dict of str to list of tuple
            ``{csv_name: [(label, kind, lake + 1), ...]}``, one record per
            target: its name with every character that is not a letter or digit
            made ``_`` and leading and trailing ``_`` stripped (``"obs"`` if
            nothing is left), and its lake number ONE-based. The shape
            ``ModflowUtlobs(continuous=)`` and a LAK package's ``observations=``
            take.
        """

        records = [
            (_observation_label(name), str(kind).upper(), int(lake_id) + 1)
            for name, lake_id in self._series.items()
        ]
        return {str(csv_name): records}

    def attach_flopy_obs(
        self,
        model,
        *,
        pname: str = "lak_obs",
        filename: str = "lake_stage_targets.obs",
        csv_name: str | None = None,
        kind: str = "STAGE",
        package=None,
    ):
        """Attach an MF6 utility observation package for lake-stage targets.

        The package is added to the FloPy objects in memory; writing and running
        the simulation afterwards is the caller's.

        Parameters
        ----------
        model : SimulationBase
            The model whose FloPy LAK package, ``model.gwf.lak``, receives the
            OBS package unless ``package`` is given (a model with no LAK raises
            ``AttributeError`` there). The records themselves need no model.
        pname : str, optional
            FloPy package name. Default ``"lak_obs"``.
        filename : str, optional
            The OBS input file name FloPy writes, in the owning model's
            directory. Default ``"lake_stage_targets.obs"``.
        csv_name : str or None, optional
            The CSV file MF6 writes the observations to. Default ``None`` takes
            ``filename``'s stem plus ``.csv`` (``"lake_stage_targets.csv"`` for
            the default ``filename``).
        kind : str, optional
            The MF6 LAK observation type, as in :meth:`to_flopy_obs`. Default
            ``"STAGE"``.
        package : flopy object or None, optional
            The FloPy package (or model) to attach the OBS package to instead of
            ``model.gwf.lak`` -- for a model with several LAK packages, say.
            Default ``None``.

        Returns
        -------
        flopy.mf6.modflow.mfutlobs.ModflowUtlobs
            The observation package.
        """

        import flopy

        owner = model.gwf.lak if package is None else package
        csv_name = _default_obs_csv_name(filename, "lake_stage_targets.csv") if csv_name is None else str(csv_name)
        continuous = self.to_flopy_obs(csv_name=csv_name, kind=kind)
        return flopy.mf6.modflow.mfutlobs.ModflowUtlobs(
            owner,
            pname=pname,
            continuous=continuous,
            filename=str(filename),
        )

    def calibration_plot(self, model, *, type: str = "calibration"):
        """Build a calibration plot directly from these targets and one model.

        Parameters
        ----------
        model : SimulationBase
            A run model with a LAK package; the plot is drawn from
            :meth:`compare`, whose ``Raises`` apply.
        type : str, optional
            Which plot:

            - ``"calibration"`` (default) or ``"obs_vs_sim"``: target stage
              against simulated stage, one point per row, titled "Observed vs
              simulated lake stage", without fit statistics;
            - ``"residuals_by_period"``: mean absolute residual per ``time``.

            ``"heads"`` RAISES ``ValueError`` (the lake comparison has no
            ``layer`` column, which that plot requires), and any other string
            draws the cross plot WITH fit statistics.

        Returns
        -------
        myflopy.modflow.calcs.calibration.CalibrationPlot
            A Plotly figure (a ``viz.Fig``).
        """

        from myflopy.modflow.calcs.calibration import CalibrationPlot

        compare = self.compare(model)
        return CalibrationPlot.from_compare(
            compare,
            type="obs_vs_sim" if type == "calibration" else type,
            target_column="stage_target",
            simulated_column="sim_stage",
            title="Observed vs simulated lake stage" if type != "heads" else "Lake stage",
            yaxis_title="Simulated stage",
        )


class BoundLakeStageTargets:
    """Model-bound lake-stage helper returned from ``model.targets.lake_stage``."""

    def __init__(self, model, targets: LakeStageTargets):
        """Bind ``targets`` to ``model`` so its methods default to that model."""

        self.model = model
        self.targets = targets
        self._plot = None

    def get(self) -> pd.DataFrame:
        """The target definitions (delegates to the wrapped targets)."""

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

    def simulated_series(self, model=None) -> pd.DataFrame:
        """The wide simulated series (defaults to the bound model)."""

        return self.targets.simulated_series(self.model if model is None else model)

    def compare(self, model=None) -> pd.DataFrame:
        """Compare observed vs simulated values (defaults to the bound model)."""

        return self.targets.compare(self.model if model is None else model)

    def residuals(self, model=None) -> pd.DataFrame:
        """The residual table (defaults to the bound model)."""

        return self.targets.residuals(self.model if model is None else model)

    def stats(self, model=None) -> pd.DataFrame:
        """Aggregate residual statistics (defaults to the bound model)."""

        return self.targets.stats(self.model if model is None else model)

    def to_flopy_obs(
        self,
        *,
        csv_name: str = "lake_stage_targets.csv",
        kind: str = "STAGE",
    ):
        """A FloPy ``continuous`` observation dict for the lake-stage targets."""

        return self.targets.to_flopy_obs(csv_name=csv_name, kind=kind)

    def attach_flopy_obs(
        self,
        *,
        pname: str = "lak_obs",
        filename: str = "lake_stage_targets.obs",
        csv_name: str | None = None,
        kind: str = "STAGE",
        package=None,
    ):
        """Attach an MF6 lake-stage observation package to the bound model."""

        return self.targets.attach_flopy_obs(
            self.model,
            pname=pname,
            filename=filename,
            csv_name=csv_name,
            kind=kind,
            package=package,
        )

    def calibration_plot(self, *, type: str = "calibration"):
        """A calibration plot for the lake-stage targets on the bound model."""

        return self.targets.calibration_plot(self.model, type=type)

    @property
    def plot(self):
        """Plotting helpers (obs-vs-sim, time series) for these lake-stage targets."""

        if self._plot is None:
            self._plot = BoundNamedSeriesTargetPlots(
                self,
                target_column="stage_target",
                simulated_column="sim_stage",
                title="Observed vs simulated lake stage",
                yaxis_title="Stage",
            )
        return self._plot


