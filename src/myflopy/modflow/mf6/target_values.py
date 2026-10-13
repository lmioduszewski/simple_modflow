"""A head or concentration target's simulated value at each stress period's END: the ONE reader.

``HeadTargets.compare()``/``simulated_heads()`` and PEST's forward-run writer
(``pest/forward_run._write_head_target_csv``, and its concentration twin) read the
values a target set is compared on through :func:`read_target_values`
(ledger 557). They were two derivations: ``compare()`` built a
``HeadsPlus`` twice per call and the whole ``kstpkper x layer x cell`` table
(7.0-9.5 s and 1.9 GB per call on Ten Trails at 38 periods), and opened
``<model>.hds`` whatever the OC declared, while the forward run read the
period-end records of the OC-declared file straight from the binary (ledger 505).

What the reader does, in order:

1. picks each period's END record on the file's own clock -- the one rule every
   target reads, :func:`~myflopy.modflow.mf6.period_end.period_end_steps`
   (ledger 336), refusing a period whose end the file does not hold;
2. reads each mapping row's ``(layer, cell)`` value from it -- the file's own
   number, MODFLOW's no-value ``1e30`` included (a target is not a view, ledgers
   351 and 355);
3. combines the rows of a SCREENED label -- several rows under one ``obsname``,
   one per layer its screen opens into -- into one value,
   :func:`~myflopy.modflow.mf6.hydraulics.screened_head` (ledger
   556), with the K it is handed: the model's NPF K for
   ``compare()``, the realization's K files for the forward run
   (:func:`k_from_sources`).

Light on purpose -- numpy, pandas and two light myflopy modules at import; FloPy
when a file is read -- because PEST's forward run imports it on every model run.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from myflopy.modflow.mf6.hydraulics import screened_head
from myflopy.modflow.mf6.period_end import period_end_steps, saved_steps

__all__ = [
    "DIFFERENCE_KINDS",
    "SCREEN_COLUMNS",
    "k_from_sources",
    "period_end_records",
    "read_target_values",
    "target_differences",
    "target_value_table",
]

#: The columns a target map carries for a SCREENED row (ledger 556):
#: the screen's interval inside the row's layer, the layer's NPF ``ICELLTYPE`` at
#: the cell, and the open length a THICKSTRT layer conducts through (``NaN`` when
#: the length follows the head). A map without them -- every map written before
#: screens, and every map of an unscreened set -- reads one cell per label.
SCREEN_COLUMNS = ("open_bottom", "open_top", "icelltype", "fixed_b")


def period_end_records(
    path,
    *,
    perlen: Sequence[float] | None,
    file: str,
    saves: str,
    text: str = "head",
) -> list[tuple[int, np.ndarray]]:
    """Each stress period's END record of a binary head or concentration file.

    Parameters
    ----------
    path : str or pathlib.Path
        The binary file (the one the OC declares, ``HEAD FILEOUT`` or
        ``CONCENTRATION FILEOUT``).
    perlen : sequence of float or None
        Each stress period's length, from the TDIS. With it a period the file
        saved without its final step is REFUSED (ledger 336). ``None`` -- a
        forward run generated before ledger 336 -- takes each period's last saved
        step unchecked.
    file : str
        What the file is called in the refusal (its name relative to the run).
    saves : str
        The OC word that schedules it (``"HEAD"``, ``"CONCENTRATION"``), for the
        refusal.
    text : str, default "head"
        FloPy's record-text selector (``"head"``, ``"concentration"``).

    Returns
    -------
    list of tuple of (int, numpy.ndarray)
        ``(per, values)`` per period the file saved, zero-based, in period order;
        ``values`` is ``(nlay, ncells)`` -- a DISV layer's cells, a DIS layer's
        ``row * ncol + col``, a DISU's nodes as one layer -- the file's own
        numbers.

    Raises
    ------
    ValueError
        With ``perlen``, when a period the file saved lacks its final step.
    FileNotFoundError
        When the file does not exist (from FloPy).
    """

    import flopy

    reader = flopy.utils.HeadFile(str(path), text=text)
    try:
        ends = period_end_steps(saved_steps(reader), perlen=perlen, file=file, saves=saves)
        records = []
        for per, kstpkper in sorted(ends.items()):
            data = np.asarray(reader.get_data(kstpkper=kstpkper), dtype=float)
            records.append((int(per), data.reshape(data.shape[0], -1)))
        return records
    finally:
        reader.close()


def _labels(mapping: pd.DataFrame) -> pd.Series:
    """Each mapping row's label: its ``obsname``, else (a map written before labels) its ``name``."""

    if "obsname" in mapping.columns:
        named = mapping["obsname"].astype(str)
        if "name" in mapping.columns:
            named = named.where(named.str.len() > 0, mapping["name"].astype(str))
        return named
    return mapping["name"].astype(str)


def target_value_table(
    records: Iterable[tuple[int, np.ndarray]],
    mapping: pd.DataFrame,
    *,
    k=None,
) -> pd.DataFrame:
    """One value per (period, label) from period-end records: the cell's, or a screen's.

    Parameters
    ----------
    records : iterable of tuple of (int, numpy.ndarray)
        ``(per, values)`` with ``values`` ``(nlay, ncells)``, as
        :func:`period_end_records` returns them.
    mapping : pandas.DataFrame
        The target map: ``obsname`` (else ``name``), ``layer``, ``cell``, one row
        per (label, layer) read. A label on ONE row reads that cell. A label on
        several rows that carry finite ``open_bottom``/``open_top`` is a screened
        well (:data:`SCREEN_COLUMNS`): its rows' values are combined by
        :func:`~myflopy.modflow.mf6.hydraulics.screened_head`, ordered by layer,
        with ``icelltype`` and ``fixed_b`` when present. A row with no screen
        interval sharing a label with another row is refused -- a map that would
        read one column twice.
    k : array_like or None, optional
        Horizontal K, ``(nlay, ncells)`` (or anything indexable as ``k[layer]``
        then reshaped to ``ncells``). Required when any label is screened; unused
        otherwise. Default ``None``.

    Returns
    -------
    pandas.DataFrame
        ``per``, then one column per label in the order its first row appears in
        ``mapping``; one row per record.

    Raises
    ------
    ValueError
        When a screened label is mapped but ``k`` is ``None``; or when one label
        holds several rows that are not all screened.
    """

    mapping = mapping.reset_index(drop=True)
    labels = _labels(mapping)
    order = list(dict.fromkeys(labels))
    layers = mapping["layer"].astype(int).to_numpy()
    cells = mapping["cell"].astype(int).to_numpy()
    screened_rows = (
        np.isfinite(pd.to_numeric(mapping["open_top"], errors="coerce").to_numpy(dtype=float))
        if "open_top" in mapping.columns else np.zeros(len(mapping), dtype=bool)
    )
    singles: dict[str, int] = {}
    screens: dict[str, np.ndarray] = {}
    for label, rows in pd.Series(np.arange(len(mapping))).groupby(labels.to_numpy(), sort=False):
        rows = rows.to_numpy()
        if len(rows) == 1 and not screened_rows[rows[0]]:
            singles[str(label)] = int(rows[0])
            continue
        if not screened_rows[rows].all():
            raise ValueError(
                f"The target map reads label {label!r} from {len(rows)} rows (layers "
                f"{sorted(int(layers[row]) for row in rows)}) without a screen interval "
                "on each; one label is one value, so only a SCREENED well's rows are "
                "combined (open_bottom/open_top, ledger 556). Rebuild the map."
            )
        screens[str(label)] = rows[np.argsort(layers[rows], kind="stable")]
    if screens and k is None:
        raise ValueError(
            f"The target map holds screened wells ({', '.join(sorted(screens)[:5])}), "
            "whose head is weighted by each layer's K, and no K was given."
        )
    if screens:
        k_array = np.asarray(k, dtype=float)
        k_array = k_array.reshape(k_array.shape[0], -1)
        bottoms = mapping["open_bottom"].astype(float).to_numpy()
        tops = mapping["open_top"].astype(float).to_numpy()
        celltypes = (
            pd.to_numeric(mapping["icelltype"], errors="coerce").fillna(1).astype(int).to_numpy()
            if "icelltype" in mapping.columns else np.ones(len(mapping), dtype=int)
        )
        fixed = (
            pd.to_numeric(mapping["fixed_b"], errors="coerce").to_numpy(dtype=float)
            if "fixed_b" in mapping.columns else np.full(len(mapping), np.nan)
        )
    rows_out = []
    for per, values in records:
        values = np.asarray(values, dtype=float)
        values = values.reshape(values.shape[0], -1)
        row: dict[str, float] = {"per": int(per)}
        for label in order:
            if label in singles:
                index = singles[label]
                row[label] = float(values[layers[index], cells[index]])
                continue
            index = screens[label]
            row[label] = float(screened_head(
                values[layers[index], cells[index]],
                k_array[layers[index], cells[index]],
                bottoms[index],
                tops[index],
                icelltype=celltypes[index],
                fixed_thickness=fixed[index],
            ))
        rows_out.append(row)
    return pd.DataFrame(rows_out, columns=["per", *order])


def k_from_sources(k_sources, *, base=".") -> np.ndarray:
    """The horizontal K a forward run's model read: one array per layer, from the files MF6 reads.

    The build records, per layer, the external file the staged model's NPF reads
    ``K`` from -- the one pyEMU's multipliers and the pilot points rewrite before
    the model runs, so the file holds THIS realization's K -- with that layer's
    ``FACTOR``; a layer written ``CONSTANT`` is recorded as its constant
    (ledger 556).

    Parameters
    ----------
    k_sources : list of dict
        One entry per layer, in layer order: ``{"layer": int, "file": str,
        "factor": float}`` or ``{"layer": int, "constant": float}``. ``file`` is
        relative to ``base`` (an absolute path is read as it is).
    base : str or pathlib.Path, default "."
        The directory the files are relative to: the forward run's working
        directory.

    Returns
    -------
    numpy.ndarray
        ``(nlay, ncells)``; a constant layer is broadcast to the file layers'
        cell count.

    Raises
    ------
    ValueError
        When the layers are not ``0 .. nlay-1`` in order, or two file layers hold
        different cell counts.
    FileNotFoundError
        When a recorded file is missing.
    """

    entries = sorted(k_sources, key=lambda entry: int(entry["layer"]))
    if [int(entry["layer"]) for entry in entries] != list(range(len(entries))):
        raise ValueError(
            f"k_sources name layers {[int(entry['layer']) for entry in entries]}; "
            "they must be every layer from 0, one entry each."
        )
    arrays: list[np.ndarray | float] = []
    for entry in entries:
        if entry.get("file") is not None:
            path = Path(base) / str(entry["file"])
            # Every number in the file, whatever its line layout: FloPy wraps an
            # external array at 20 values a line, pyEMU writes it its own way.
            values = np.array(path.read_text().split(), dtype=float)
            arrays.append(values * float(entry.get("factor", 1.0)))
        else:
            arrays.append(float(entry["constant"]))
    sizes = {array.size for array in arrays if isinstance(array, np.ndarray)}
    if len(sizes) > 1:
        raise ValueError(f"The K files hold different cell counts {sorted(sizes)}.")
    if not sizes:
        return np.asarray(arrays, dtype=float).reshape(len(arrays), 1)
    ncells = sizes.pop()
    return np.vstack([
        array if isinstance(array, np.ndarray) else np.full(ncells, array) for array in arrays
    ])


def read_target_values(
    path,
    mapping: pd.DataFrame,
    *,
    perlen: Sequence[float] | None,
    file: str,
    saves: str,
    text: str = "head",
    k=None,
) -> pd.DataFrame:
    """The simulated value of each target label at each period's END, from one binary file.

    :func:`period_end_records` then :func:`target_value_table`: what
    ``HeadTargets.compare()``/``simulated_heads()`` and the PEST forward-run
    writers read (ledger 557).

    Parameters
    ----------
    path : str or pathlib.Path
        The OC-declared binary head or concentration file.
    mapping : pandas.DataFrame
        The target map (see :func:`target_value_table`).
    perlen : sequence of float or None
        Each period's length (see :func:`period_end_records`).
    file : str
        The file's name for the refusal.
    saves : str
        The OC word scheduling it (``"HEAD"``/``"CONCENTRATION"``).
    text : str, default "head"
        FloPy's record-text selector.
    k : array_like or None, optional
        Horizontal K per layer and cell; required when the map holds a screened
        well. Default ``None``.

    Returns
    -------
    pandas.DataFrame
        ``per`` plus one column per label, in map order.
    """

    records = period_end_records(path, perlen=perlen, file=file, saves=saves, text=text)
    return target_value_table(records, mapping, k=k)


#: The two derived observations :func:`target_differences` makes (ledger
#: 565): ``"mean"`` -- each value minus its label's mean over its
#: measured periods -- and ``"previous"`` -- each measured period minus the previous
#: one of that label.
DIFFERENCE_KINDS = ("mean", "previous")


def target_differences(
    table: pd.DataFrame, periods: pd.DataFrame, *, kind: str, undefined: float = np.nan
) -> pd.DataFrame:
    """Each label's DIFFERENCES over its own measured periods: the ONE rule (ledger 565).

    What ``HeadDifferences`` compares -- applied to the targets and to the simulated
    values alike -- and what PEST's forward run writes
    (``forward_run._write_target_difference_csv``), so a calibration and the
    off-PEST review read one definition. A period the label was not measured in, and
    a label's first measured period under ``kind="previous"``, has no difference:
    ``undefined`` there (the forward run passes ``0.0``, which the build weights 0,
    because PEST++ fails a run on an empty field, ledger 351). A NaN VALUE at a
    measured period is not that: it stays ``NaN`` whatever ``undefined`` is, so a
    realization whose model wrote no number fails its run rather than scoring a flat
    hydrograph (ledger 565).

    Parameters
    ----------
    table : pandas.DataFrame
        ``per`` (zero-based stress period) plus one column of values per label, as
        :func:`read_target_values` returns.
    periods : pandas.DataFrame
        ``obsname`` (a label of ``table``) and ``per``: each label's MEASURED periods.
        A label it does not name has no difference anywhere.
    kind : str
        ``"mean"``: value minus the mean of the label's values at its measured
        periods (a ``NaN`` among them makes every difference of the label ``NaN``);
        ``"previous"``: value minus the value at the label's previous measured
        period.
    undefined : float, default numpy.nan
        What a row with NO difference holds -- a period the label was not measured
        in, its first measured period under ``"previous"``, every period of a label
        ``periods`` does not name. Never applied where a difference is defined and
        its values are ``NaN``.

    Returns
    -------
    pandas.DataFrame
        ``table``'s shape and columns: ``per`` unchanged, each label's differences.

    Raises
    ------
    ValueError
        When ``kind`` is not one of :data:`DIFFERENCE_KINDS`, or ``periods`` names a
        period ``table`` holds no row for.
    """

    if kind not in DIFFERENCE_KINDS:
        raise ValueError(f"kind={kind!r} is not a difference; pass 'mean' or 'previous'.")
    out = table.copy()
    rows = {int(per): position for position, per in enumerate(table["per"].astype(int))}
    for label in table.columns[1:]:
        measured = sorted({int(per) for per in periods.loc[
            periods["obsname"].astype(str) == str(label), "per"]})
        missing = [per for per in measured if per not in rows]
        if missing:
            raise ValueError(
                f"Label {label!r} is measured in period(s) {missing}, which the simulated "
                f"table holds no row for (it has {sorted(rows)})."
            )
        values = table[label].to_numpy(dtype=float)
        result = np.full(len(table), float(undefined))
        positions = [rows[per] for per in measured]
        if positions:
            picked = values[positions]
            if kind == "mean":
                result[positions] = picked - np.mean(picked)
            else:
                result[positions[1:]] = np.diff(picked)
        out[label] = result
    return out
