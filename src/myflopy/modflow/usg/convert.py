"""Map an imported MODFLOW-USG model onto a MODFLOW 6 :class:`SimulationSpec`.

Every judgment this module makes is recorded in :func:`conversion_notes`, which
is what :meth:`UsgModel.report` prints. The rule is that a conversion may be
lossy but never silently lossy: anything approximated or omitted is counted and
named, with the cells it touched.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

import myflopy.package_api as api
from myflopy._logging import get_logger
from myflopy.modflow.mf6.boundary_support import bottom_clearance, river_pair_lift
from myflopy.specs import ModelContext, PackageSpec, SimulationSpec

if TYPE_CHECKING:
    from myflopy.modflow.usg.model import UsgModel

logger = get_logger(__name__)

__all__ = ["Finding", "conversion_notes", "mf6_findings", "to_mf6"]

#: USG list packages and the MODFLOW 6 factory that receives them.
_LIST_TARGETS = {"CHD": "chd", "DRN": "drn", "GHB": "ghb", "RIV": "riv", "WEL": "wel"}

#: The packages ``to_mf6(include=...)`` chooses among. Everything else it writes --
#: DISV, NPF, IC, STO, OC -- is written whatever ``include`` says.
_OPTIONAL_PACKAGES = (*_LIST_TARGETS.values(), "rch", "evt", "hfb")

#: The elevation fields MODFLOW 6 measures against the cell bottom, and the column each
#: one occupies in a MODFLOW-USG record (column 0 is the node). A USG ``RIV`` record is
#: ``node, stage, cond, rbot`` -- the same order MODFLOW 6 writes -- and ``rbot`` is the
#: field MODFLOW 6 actually checks. A river whose BED is below its cell bottom while its
#: water surface is not is the ordinary case, and a stage-only check cannot see it at
#: all: measured, stage 55 / rbot 40 in a cell bottomed at 48 is refused by name at
#: ``icelltype = 1`` and was reported by nothing (ledger 291).
_BOTTOM_CHECKED: dict[str, tuple[tuple[str, int], ...]] = {
    "DRN": (("elev", 1),),
    "GHB": (("bhead", 1),),
    "RIV": (("rbot", 3), ("stage", 1)),
}

#: Where each family keeps its conductance, for saying how much of the package the
#: offending records carry.
_CONDUCTANCE_COLUMN = {"DRN": 2, "GHB": 2, "RIV": 2}

#: What each family's clamp says when it runs. The GHB wording is deliberately
#: unchanged: it is accurate, and two notebooks carry it in their saved output.
_CLAMP_REPORT = {
    "GHB": "%s: raised %d boundary head(s) to their cell bottom so MODFLOW 6 accepts them "
           "(fix_for_mf6=True)",
    "DRN": "%s: raised %d drain elevation(s) to their cell bottom so MODFLOW 6 accepts "
           "them (fix_for_mf6=True)",
    "RIV": "%s: raised %d river bed(s) to their cell bottom, carrying a stage up with its "
           "bed only where the stage was below the bottom too (fix_for_mf6=True)",
}

#: What ``fix_for_mf6=True`` does about each one, in the finding's own words.
_BOTTOM_FIXES = {
    ("DRN", "elev"): (
        "to_mf6(fix_for_mf6=True) raises each elev to just above its cell bottom"
    ),
    ("GHB", "bhead"): (
        "to_mf6(fix_for_mf6=True) raises each bhead to just above its cell bottom"
    ),
    ("RIV", "rbot"): (
        "to_mf6(fix_for_mf6=True) raises each rbot to just above its cell bottom, and "
        "carries the stage up by the same delta where the stage was below the bottom "
        "too, so the declared river depth survives; a stage already inside its cell "
        "is left where it was drawn"
    ),
    ("RIV", "stage"): (
        "to_mf6(fix_for_mf6=True) carries each stage up with its own river bottom where "
        "that bottom is lifted; where the bottom is already inside the cell, a stage "
        "below it means the stage is under its own bed -- see riv_stage_below_river_bottom"
    ),
}


@dataclass(frozen=True, slots=True)
class Finding:
    """One way in which MODFLOW 6 is stricter than MODFLOW-USG about this model.

    MODFLOW-USG accepts inputs MODFLOW 6 refuses -- most often a head-type
    boundary placed below the bottom of the cell it sits in, which USG treats as
    a strong sink. Finding these before the run is worth more than discovering
    them from a failed run, because the listing file names six of them and stops.

    **Whether MODFLOW 6 actually refuses depends on the layer**, so ``severity``
    is read from the model, never assumed (ledger 291). ``"reject"`` means
    MODFLOW 6 will not start; ``"allow"`` means it starts and the converted model
    behaves differently from the one the modeller drew -- silently, which is the
    reason a non-blocking finding is worth reporting at all.
    """

    code: str
    severity: str
    count: int
    detail: str
    fix: str

    @property
    def blocks_run(self) -> bool:
        """True when MODFLOW 6 will refuse to run until this is resolved.

        False does NOT mean "harmless": an ``"allow"`` finding is one MODFLOW 6
        accepts and runs, having changed what the boundary does.
        """

        return self.severity == "reject"


def mf6_findings(model: UsgModel) -> list[Finding]:
    """Return everything MODFLOW 6 will object to in this model, before running it."""

    return [
        *_below_bottom_findings(model),
        *_river_pair_findings(model),
        *_chd_findings(model),
    ]


def _convertible(model: UsgModel, layers: np.ndarray) -> bool:
    """True when any of these layers is convertible, which is the whole severity split."""

    return bool(np.asarray(model.icelltype)[layers].any())


def _below_bottom_findings(model: UsgModel) -> list[Finding]:
    """Findings for every elevation MODFLOW 6 measures against its own cell's bottom.

    **What MODFLOW 6 does with one depends on ``icelltype``, and a converted USG model
    lands on either side of that split** -- ``icelltype`` is the LPF ``LAYTYP``, carried
    straight across. Where the layer is convertible MODFLOW 6 REFUSES to start and names
    the record, and the refusal is a per-period data read, so a model whose first
    offending period is late dies after MODFLOW 6 has solved the earlier ones. Where the
    layer is confined MODFLOW 6 accepts the record, terminates normally, and drains
    against a head below the cell it sits in -- measured 2.0x the discharge of the same
    drain placed inside its cell, with nothing said in the listing file.

    So the severity is read from the model rather than asserted (ledger 291). Both halves
    used to be asserted, and asserted wrongly: RIV was inspected on ``stage`` where
    MODFLOW 6 checks ``rbot``, and DRN was exempted outright on the premise that a drain
    below its cell bottom is simply inactive, which neither branch supports.

    The test here is ``< bottom``, while the clamp lifts anything below ``bottom +
    clearance``. The asymmetry is deliberate and in the safe direction: a value sitting
    exactly ON its cell bottom is not *below* it, so saying so would be false, but it
    can still round below in the written file -- so the correction covers it while the
    report does not claim it.
    """

    findings: list[Finding] = []
    for ftype, checks in _BOTTOM_CHECKED.items():
        records = model.boundaries.get(ftype)
        if records is None:
            continue
        for label, column in checks:
            # Every number in the sentence comes from ONE period -- the one with the most
            # offenders. Taking each as a max over periods independently reported a
            # share/total pair that occurred in no period at all (a constructed case
            # understated the share 50x), and a count that read as a record total when it
            # was a per-period maximum.
            worst_period = None
            convertible = False
            conductance = _CONDUCTANCE_COLUMN[ftype]
            for block in records.periods.values():
                if not len(block) or block.shape[1] <= column:
                    continue
                cellid = model.to_cellid(block[:, 0].astype(np.int64))
                valid = cellid[:, 0] >= 0
                if not valid.any():
                    continue
                layers = cellid[valid, 0]
                bottom = model.botm[layers, cellid[valid, 1]]
                below = block[valid, column] < bottom
                if not below.any():
                    continue
                convertible = convertible or _convertible(model, layers[below])
                share = total = 0.0
                if block.shape[1] > conductance:
                    share = float(block[valid, conductance][below].sum())
                    total = float(block[valid, conductance].sum())
                here = (
                    int(below.sum()),
                    float((bottom - block[valid, column])[below].max()),
                    share,
                    total,
                )
                if worst_period is None or here[0] > worst_period[0]:
                    worst_period = here
            if worst_period is None:
                continue
            offenders, worst, share, total = worst_period
            portion = f", carrying {share:.3g} of the period's {total:.3g} conductance" if total else ""
            verdict = (
                "MODFLOW-USG allows this; MODFLOW 6 refuses to start and names the "
                "record. The refusal is read per stress period, so a model whose first "
                "offending period is late fails only after the earlier ones have solved."
                if convertible
                else "Every layer holding one is CONFINED (icelltype 0), where MODFLOW 6 "
                "accepts the record and terminates normally -- but the boundary then "
                "works against a head below its own cell and exchanges more water with "
                "the model than the one that was drawn, with nothing said in the "
                "listing file."
            )
            findings.append(
                Finding(
                    code=f"{ftype.lower()}_{label}_below_cell_bottom",
                    severity="reject" if convertible else "allow",
                    count=offenders,
                    detail=(
                        f"{offenders} {ftype} record(s) have {label} below the bottom of "
                        f"their cell in one stress period, the deepest by {worst:.2f}"
                        f"{portion}. {verdict}"
                    ),
                    fix=_BOTTOM_FIXES[ftype, label],
                )
            )
    return findings


def _river_pair_findings(model: UsgModel) -> list[Finding]:
    """A river whose stage sits below its own bed, which MODFLOW 6 refuses either way.

    Unlike everything else here this does not depend on ``icelltype``: MODFLOW 6 checks
    ``stage >= rbot`` as a property of the record, and refuses a confined model for it
    too. ``fix_for_mf6`` deliberately does not touch it -- which of the two numbers is
    wrong is not something a converter can know, and a pair lift preserves the inversion
    rather than repairing it.
    """

    records = model.boundaries.get("RIV")
    if records is None:
        return []
    offenders = 0
    worst = 0.0
    for block in records.periods.values():
        if not len(block) or block.shape[1] < 4:
            continue
        inverted = block[:, 1] < block[:, 3]
        if not inverted.any():
            continue
        offenders = max(offenders, int(inverted.sum()))
        worst = max(worst, float((block[:, 3] - block[:, 1])[inverted].max()))
    if not offenders:
        return []
    return [
        Finding(
            code="riv_stage_below_river_bottom",
            severity="reject",
            count=offenders,
            detail=(
                f"{offenders} RIV record(s) give a stage below their own river bottom "
                f"(by up to {worst:.2f}). MODFLOW 6 refuses to start at any icelltype: a "
                "river's stage is its water surface and its rbot the bed beneath it."
            ),
            fix=(
                "correct the source RIV records -- fix_for_mf6 will not guess which of "
                "the two numbers is the wrong one"
            ),
        )
    ]


def _chd_findings(model: UsgModel) -> list[Finding]:
    """Constant heads below their cell bottom, whose meaning also splits on icelltype."""

    chd = model.boundaries.get("CHD")
    if chd is None:
        return []
    offenders = 0
    convertible = False
    for block in chd.periods.values():
        if not len(block):
            continue
        cellid = model.to_cellid(block[:, 0].astype(np.int64))
        valid = cellid[:, 0] >= 0
        if not valid.any():
            continue
        layers = cellid[valid, 0]
        bottom = model.botm[layers, cellid[valid, 1]]
        head = block[valid, 2] if block.shape[1] > 2 else block[valid, 1]
        below = head < bottom
        if not below.any():
            continue
        offenders = max(offenders, int(below.sum()))
        convertible = convertible or _convertible(model, layers[below])
    if not offenders:
        return []
    fix = (
        "to_mf6(fix_for_mf6=True) omits each such record in the periods where it is "
        "below the bottom, and keeps it in the periods where it is not"
    )
    if convertible:
        return [
            Finding(
                code="chd_head_below_cell_bottom",
                severity="reject",
                count=offenders,
                detail=(
                    f"up to {offenders} CHD record(s) per period hold a head below their "
                    "cell bottom. MODFLOW 6 refuses to start and names the record. "
                    "MODFLOW-USG accepts it, but a convertible cell then has zero "
                    "saturated thickness and passes no water -- so the boundary was "
                    "already inert wherever this happens."
                ),
                fix=fix,
            )
        ]
    return [
        Finding(
            code="chd_head_below_cell_bottom",
            severity="allow",
            count=offenders,
            detail=(
                f"up to {offenders} CHD record(s) per period hold a head below their cell "
                "bottom. Every layer holding one is CONFINED (icelltype 0), where MODFLOW "
                "6 accepts the record and runs -- and where the cell keeps its full "
                "thickness, so the boundary is NOT inert: it pins the head below the "
                "cell and pulls water out (measured -1200 ft3/d on a one-cell bench). "
                "Dropping it is a real change to the water balance, not a tidy-up."
            ),
            fix=fix,
        )
    ]


def to_mf6(
    model: UsgModel,
    *,
    name: str = "usg",
    crs: str | None = None,
    newton: bool | None = None,
    under_relaxation: bool = True,
    start_date_time: str | None = None,
    save_budget: bool = True,
    complexity: str = "COMPLEX",
    include: str | tuple[str, ...] | None = None,
    fix_for_mf6: bool = False,
    local_origin: bool = True,
) -> SimulationSpec:
    """Build a MODFLOW 6 simulation spec from a USG model. See :meth:`UsgModel.to_mf6`."""

    if model.grid is None:
        raise ValueError(
            "conversion needs the grid geometry: a USG DISU carries no coordinates, so "
            "read_usg() must be given the model's .gsf file"
        )
    wanted = _included(include)

    if not fix_for_mf6:
        for finding in mf6_findings(model):
            if finding.blocks_run:
                logger.warning(
                    "MODFLOW 6 will reject this model: %s %s", finding.detail, finding.fix
                )
            else:
                # An `allow` finding is the one nothing else says out loud: MODFLOW 6
                # runs, and the boundary does something other than what was drawn.
                # Leaving it to validate()/report() alone is ledger 288's defect one
                # module over -- a degradation reported through a channel nobody hears.
                logger.warning(
                    "MODFLOW 6 will run this model but change what a boundary does: %s %s",
                    finding.detail,
                    finding.fix,
                )

    keep = lambda pkg: wanted is None or pkg in wanted  # noqa: E731

    grid_props = model.grid.get_disv_gridprops()
    vertices, cell2d, xorigin, yorigin = _localize(grid_props, shift=local_origin)
    packages: list[PackageSpec] = [
        api.disv(
            nlay=model.nlay,
            ncpl=grid_props["ncpl"],
            nvert=len(vertices),
            vertices=vertices,
            cell2d=cell2d,
            top=model.top,
            botm=model.botm,
            idomain=model.idomain,
            xorigin=xorigin,
            yorigin=yorigin,
        ),
        # PRT reads specific discharge and saturation from the flow model's budget;
        # without them every particle-tracking run on a converted model failed.
        api.npf(
            k=model.k,
            k33=model.k33,
            icelltype=model.icelltype,
            save_flows=True,
            save_specific_discharge=True,
            save_saturation=True,
        ),
        api.ic(strt=model.strt),
    ]

    # Pass a period map only when it has entries. An EMPTY dict is not the same
    # as omitting it: FloPy treats `steady_state={}` as a supplied period map and
    # writes `BEGIN period / END period` with no keyword inside, which silently
    # loses the TRANSIENT flags as well -- MODFLOW 6 then defaults every period
    # to steady state, and the model solves to NaN rather than failing loudly.
    steady_flags = np.asarray(model.disu.steady, dtype=bool)
    transient = {i: True for i, flag in enumerate(steady_flags) if not flag}
    steady = {i: True for i, flag in enumerate(steady_flags) if flag}
    packages.append(
        api.sto(
            ss=model.ss,
            sy=model.sy,
            iconvert=np.broadcast_to(model.icelltype[:, None], (model.nlay, model.ncpl)).copy(),
            transient=transient or None,
            steady_state=steady or None,
        )
    )

    for ftype, target in _LIST_TARGETS.items():
        records = model.boundaries.get(ftype)
        if records is None or not keep(target):
            continue
        data = _list_period_data(model, records, fix_for_mf6=fix_for_mf6)
        if not data:
            logger.debug("%s has no groundwater records after CLN filtering; skipped", ftype)
            continue
        factory = getattr(api, target)
        packages.append(factory.flopy(stress_period_data=data, name=target))

    if model.rch is not None and keep("rch"):
        packages.append(_recharge_spec(model))
    if model.ets is not None and keep("evt"):
        packages.append(_evt_spec(model))
    if model.hfb.size and keep("hfb"):
        spec = _hfb_spec(model)
        if spec is not None:
            packages.append(spec)

    packages.append(
        api.oc(
            head_filerecord=f"{name}.hds",
            budget_filerecord=f"{name}.cbc" if save_budget else None,
            saverecord=[("HEAD", "ALL")] + ([("BUDGET", "ALL")] if save_budget else []),
        )
    )

    use_newton = newton if newton is not None else bool(model.sms and model.sms.is_newton)
    newton_options = None
    if use_newton:
        newton_options = "NEWTON UNDER_RELAXATION" if under_relaxation else "NEWTON"

    flow = api.gwf(
        name,
        packages=packages,
        context=ModelContext(grid=model.grid),
        newtonoptions=newton_options,
        save_flows=True,
    )

    perioddata = [
        (float(p), int(n), float(t))
        for p, n, t in zip(model.disu.perlen, model.disu.nstp, model.disu.tsmult, strict=True)
    ]
    simulation_packages = [
        api.tdis(
            nper=len(perioddata),
            perioddata=perioddata,
            time_units=model.time_units,
            start_date_time=start_date_time or model.start_date_time,
        ),
        _ims_spec(model, name=name, complexity=complexity),
    ]

    logger.info(
        "converted %s to MODFLOW 6: %d packages on a %d x %d DISV grid, %d periods",
        model.name_file.path.name,
        len(packages),
        model.nlay,
        model.ncpl,
        len(perioddata),
    )
    if crs:
        model.grid.crs = crs
    return SimulationSpec(name, models=(flow,), packages=tuple(simulation_packages))


# ---------------------------------------------------------------- pieces


def _included(include) -> set[str] | None:
    """``to_mf6``'s ``include=`` as a set of package names, or ``None`` for all of them.

    A bare string is ONE name: ``include="chd"`` keeps CHD. Iterated as given it was
    read letter by letter, and since no package is called ``"c"`` it kept none, in
    silence -- as did a misspelt name in a tuple. A name outside
    :data:`_OPTIONAL_PACKAGES`, or an item that is not a string, is refused naming the
    names it can be (ledger 429).
    """

    if include is None:
        return None
    names = (include,) if isinstance(include, str) else tuple(include)
    valid = ", ".join(repr(name) for name in _OPTIONAL_PACKAGES)
    for name in names:
        if not isinstance(name, str):
            raise TypeError(
                f"to_mf6(include={include!r}): {name!r} is not a package name; include= "
                f"takes names among {valid}"
            )
    unknown = [name for name in names if name.lower() not in _OPTIONAL_PACKAGES]
    if unknown:
        raise ValueError(
            f"to_mf6(include={include!r}): {', '.join(repr(n) for n in unknown)} "
            f"{'is not an optional package' if len(unknown) == 1 else 'are not optional packages'}. "
            f"include= chooses among {valid}; DISV, NPF, IC, STO and OC are always written"
        )
    return {name.lower() for name in names}


def _localize(grid_props: dict, *, shift: bool):
    """Move the mesh next to the origin and declare where it came from.

    MODFLOW 6 builds DISV connection geometry -- shared edge lengths and the
    centre-to-edge distances that become conductances -- from the raw vertex
    coordinates. On a State Plane grid those are around 1.34 MILLION feet, and
    the differences that matter are a few feet, so the subtraction loses most of
    its significant digits. For the smallest cells here (under 3 square feet)
    the result degrades to NaN, and MODFLOW 6 then reports a NaN budget and
    "Normal termination" rather than an error: the run looks like it worked.

    Measured on the Ten Trails grid: identical input, run twice, differing only
    in this shift -- 0 of 9405 cells finite as-is, 9405 of 9405 finite shifted.

    Shifting the mesh to its own lower-left corner and passing that corner as
    ``xorigin``/``yorigin`` keeps the model georeferenced (FloPy's model grid and
    every export put the cells back where they belong) while the arithmetic
    happens near zero.
    """

    vertices = grid_props["vertices"]
    cell2d = grid_props["cell2d"]
    if not shift:
        return vertices, cell2d, 0.0, 0.0

    xorigin = float(min(v[1] for v in vertices))
    yorigin = float(min(v[2] for v in vertices))
    moved_vertices = [[v[0], v[1] - xorigin, v[2] - yorigin] for v in vertices]
    moved_cells = [
        [r[0], r[1] - xorigin, r[2] - yorigin, r[3], *r[4 : 4 + r[3]]] for r in cell2d
    ]
    logger.debug("shifted the mesh by (%.1f, %.1f) and declared it as the DISV origin", xorigin, yorigin)
    return moved_vertices, moved_cells, xorigin, yorigin


def _list_period_data(model: UsgModel, records, *, fix_for_mf6: bool = False) -> dict[int, list[Any]]:
    """Convert one list boundary condition to a MODFLOW 6 period dict.

    A USG period whose ``ITMP`` was negative reused the previous period. MODFLOW
    6 expresses exactly that by omitting the period from the dict, so a reused
    period is simply absent here -- no data is duplicated and no period is
    accidentally emptied.

    Every period that DEFINED its own records is a key, even when none of them
    reaches MODFLOW 6: ``ITMP = 0`` (the package switched off), a period whose
    records are all on CLN nodes, and a period whose every CHD record
    ``fix_for_mf6`` omitted each become an EMPTY list, which MODFLOW 6 writes as
    an empty period block and reads as "delete every record". Skipping such a
    period left the previous period's boundaries running through it (ledger
    367). A package with no record in ANY period returns an empty dict, so
    the caller writes no package at all.
    """

    ftype = records.ftype
    data: dict[int, list[Any]] = {}
    clamped = [0]
    dropped = [0]
    dropped_convertible = False
    for period, block in sorted(records.periods.items()):
        if len(block) == 0:
            # ITMP = 0, or every record of the period was on a CLN node: the period
            # is DEFINED and holds nothing, which is not the same as reusing the last.
            data[period] = []
            continue
        nodes = block[:, 0].astype(np.int64)
        cellid = model.to_cellid(nodes)
        valid = cellid[:, 0] >= 0
        values_block = block[valid, 1:]
        if fix_for_mf6 and ftype == "CHD":
            # A constant head below its cell bottom is fatal in MODFLOW 6 where the
            # layer is CONVERTIBLE, and there it is inert in MODFLOW-USG too (zero
            # saturated thickness), so omitting it reproduces the USG behaviour
            # exactly. In a CONFINED layer neither half holds: the cell keeps its full
            # thickness, MODFLOW 6 accepts the record, and dropping it removes a
            # boundary that was moving water (measured -1200 ft3/d on a one-cell
            # bench). The drop is unconditional all the same -- see ledger 291, where
            # gating it is recorded as open rather than guessed at.
            # It has to be decided PER PERIOD: these heads follow a seasonal cycle, so
            # the same record is below the bottom in a low-stage month and above it in
            # a high one.
            bottom = model.botm[cellid[valid, 0], cellid[valid, 1]]
            head = values_block[:, 1] if values_block.shape[1] > 1 else values_block[:, 0]
            keep_row = head > bottom + bottom_clearance(bottom)
            dropped_convertible = dropped_convertible or (
                not keep_row.all() and _convertible(model, cellid[valid, 0][~keep_row])
            )
            if not keep_row.all():
                dropped[0] += int((~keep_row).sum())
                cellid = cellid[valid][keep_row]
                values_block = values_block[keep_row]
                valid = np.ones(len(cellid), dtype=bool)

        # `values_block` drops the node, so its columns are the record's trailing
        # values: (bhead, cond) / (elev, cond) / (stage, cond, rbot).
        column = 2 if ftype == "RIV" else 0
        floorable = fix_for_mf6 and ftype in ("GHB", "DRN", "RIV")
        if floorable and values_block.shape[1] <= column:
            # A record short of its own fields is malformed upstream. Leaving it as
            # read keeps the failure where it belongs, rather than turning it into an
            # IndexError inside the clamp -- or, worse, dropping the period.
            logger.debug(
                "%s: %d value column(s) in period %d, too few to floor %s; left as read",
                ftype, values_block.shape[1], period, "rbot" if column else "the elevation",
            )
            floorable = False
        if floorable:
            bottom = model.botm[cellid[valid, 0], cellid[valid, 1]]
            # Clamp to just ABOVE the bottom, not to it. MODFLOW 6 compares the two
            # numbers as it read them back from text, and both are written with
            # finite precision -- so a value set exactly equal to the bottom still
            # trips "IS LESS THAN CELL BOTTOM" once the last digit rounds the wrong
            # way. The nudge is relative so it stays negligible at any elevation
            # (5e-5 ft at 50 ft); `bottom_clearance` is the same rule the GeoPackage
            # boundary sources floor by.
            floor = bottom + bottom_clearance(bottom)
            below = values_block[:, column] < floor
            if below.any():
                values_block = values_block.copy()
                if ftype == "RIV":
                    # A RIV record is a PAIR. MODFLOW 6 checks `rbot` against the cell
                    # bottom and `stage >= rbot` separately, so lifting `rbot` alone
                    # turns a model that ran into a hard failure WHERE THE STAGE WAS
                    # BELOW THE BOTTOM TOO, and clamping both to the floor collapses
                    # the declared depth to zero. A stage already inside its cell is a
                    # legal, physical number and is left exactly where it was drawn --
                    # raising it unconditionally put one 3 ft above land surface and
                    # reversed the river's sign (ledger 291).
                    values_block[below, 0], values_block[below, 2] = river_pair_lift(
                        values_block[below, 0], values_block[below, 2], floor[below]
                    )
                else:
                    values_block[below, column] = floor[below]
                clamped[0] += int(below.sum())
        rows = []
        for (layer, cell), values in zip(cellid[valid], values_block, strict=True):
            rows.append(_bc_record(ftype, (int(layer), int(cell)), values))
        # An empty `rows` is kept: the period was defined, and omitting it would
        # tell MODFLOW 6 to reuse the previous period's records (ledger 367).
        data[period] = rows
    if dropped[0]:
        # The convertible wording is byte-identical to what this always said, and is
        # true there. The confined one is new, because there the record was neither
        # inert in MODFLOW-USG nor rejected by MODFLOW 6.
        logger.warning(
            (
                "%s: omitted %d record-period(s) whose head sits below the cell bottom -- "
                "inert in MODFLOW-USG, rejected by MODFLOW 6 (fix_for_mf6=True)"
                if dropped_convertible
                else "%s: omitted %d record-period(s) whose head sits below the cell "
                "bottom in a CONFINED layer, where MODFLOW 6 would have accepted them and "
                "the cell keeps its full thickness -- this removes water the model was "
                "moving (fix_for_mf6=True)"
            ),
            ftype,
            dropped[0],
        )
    if clamped[0]:
        # One message per family, each naming the field it actually moved. A single
        # generic line would have to call a drain invert and a river bed a "head".
        logger.warning(_CLAMP_REPORT[ftype], ftype, clamped[0])
    if not any(data.values()):
        # Not one record in any period: there is no package to write.
        return {}
    return data


def _bc_record(ftype: str, cellid: tuple[int, int], values: np.ndarray) -> list[Any]:
    """Build one MODFLOW 6 boundary record from a USG record's trailing values.

    ``CHD`` is the one that loses information: MODFLOW-USG ramps a constant head
    from ``shead`` to ``ehead`` across the stress period, while MODFLOW 6 holds
    one value for the whole period. The end-of-period head is taken, because
    that is the value MODFLOW-USG itself reaches by the end of the period and so
    the one a following period continues from.
    """

    if ftype == "CHD":
        return [cellid, float(values[1] if values.size > 1 else values[0])]
    return [cellid, *(float(v) for v in values)]


def _recharge_spec(model: UsgModel) -> PackageSpec:
    """Build MODFLOW 6 ``RCHA`` from the USG ``RCH`` arrays.

    USG recharge is already one array per period, so the array form is the
    faithful shape as well as the compact one -- on the Ten Trails model, 0.88 s
    end to end against 5.35 s for the equivalent list.

    ``NRCHOP = 1`` (top grid layer) and ``NRCHOP = 3`` (highest ACTIVE cell) both
    resolve through ``irch``: MODFLOW 6 without an IRCH array applies recharge to
    layer 1 unconditionally and silently skips any column whose layer 1 is
    inactive, which is not what either option means.
    """

    recharge = {period: np.asarray(values) for period, values in sorted(model.rch.rech.items())}
    return api.rch.array(
        recharge=recharge,
        context=ModelContext(grid=model.grid, domain=model.idomain),
        name="rch",
    )


def _evt_spec(model: UsgModel) -> PackageSpec:
    """Build MODFLOW 6 ``EVT`` from the USG ``ETS`` package.

    Two mismatches are resolved here, both forced:

    ``NETSOP = 3`` applies ET to the highest *active* cell in each column.
    MODFLOW 6 has no such option, so the cell is resolved once from IBOUND --
    which is static for the whole run, making this exact rather than an
    approximation.

    ``NETSEG > 1`` (segmented ET) cannot use the array-based ``EVTA`` at all,
    because ``READASARRAYS`` and segments are mutually exclusive in MODFLOW 6.
    The package therefore becomes list-based, one record per active column per
    period, which is a large but honest translation of the same physics.
    """

    ets = model.ets
    columns = np.flatnonzero(model.has_active_column)
    layers = model.uppermost_active[columns]
    nseg = max(int(ets.netseg), 1)

    def latest(store: dict[int, Any], period: int, default=None):
        """Value for ``period``, falling back to the most recent one defined."""

        if period in store:
            return store[period]
        earlier = [p for p in store if p < period]
        return store[max(earlier)] if earlier else default

    data: dict[int, list[Any]] = {}
    for period in range(model.nper):
        rate = latest(ets.rate, period)
        if rate is None:
            continue
        surf = latest(ets.surf, period, np.zeros(model.ncpl))
        depth = latest(ets.depth, period, np.ones(model.ncpl))
        pxdp = latest(ets.pxdp, period, [])
        petm = latest(ets.petm, period, [])
        rows = []
        for layer, cell in zip(layers, columns, strict=True):
            record: list[Any] = [
                (int(layer), int(cell)),
                float(surf[cell]),
                float(rate[cell]),
                float(depth[cell]),
            ]
            record += [float(a[cell]) for a in pxdp[: nseg - 1]]
            record += [float(a[cell]) for a in petm[: nseg - 1]]
            rows.append(record)
        data[period] = rows

    logger.info(
        "ETS -> list-based EVT: %d periods x %d columns = %d records (nseg=%d)",
        len(data),
        columns.size,
        sum(len(v) for v in data.values()),
        nseg,
    )
    return api.evt.flopy(stress_period_data=data, nseg=nseg, name="evt")


def _hfb_spec(model: UsgModel) -> PackageSpec | None:
    """Build MODFLOW 6 ``HFB`` from USG node pairs.

    Vertical barriers ARE legal on DISV as of MODFLOW 6 6.7.0, so a cross-layer
    pair is no longer dropped on sight -- it is handed to the same validator every
    other route uses, which accepts a pair in adjacent layers and rejects one that
    skips a layer. (FloPy 3.10's embedded definition still asserts the old
    same-layer rule and is out of date.)

    ``mf.hfb`` also dedupes: MODFLOW-USG tolerates a repeated barrier face, while
    MODFLOW 6 silently applies its series formula twice and then leaves that face
    permanently wrong.
    """

    if not len(model.hfb):
        return None

    first = model.to_cellid(model.hfb[:, 0].astype(np.int64))
    second = model.to_cellid(model.hfb[:, 1].astype(np.int64))
    on_grid = (first[:, 0] >= 0) & (second[:, 0] >= 0)
    if not on_grid.any():
        logger.warning("every HFB barrier references a CLN node; HFB omitted")
        return None
    if not on_grid.all():
        logger.warning(
            "%d of %d HFB barriers reference CLN nodes and were dropped with the CLN network",
            int((~on_grid).sum()),
            on_grid.size,
        )

    pairs = [
        ((int(a[0]), int(a[1])), (int(b[0]), int(b[1])))
        for a, b in zip(first[on_grid], second[on_grid], strict=True)
    ]
    return api.hfb(
        pairs=pairs,
        hydchr=[float(v) for v in model.hfb[on_grid, 2]],
        context=ModelContext(grid=model.grid),
        name="hfb",
    )


def _ims_spec(model: UsgModel, *, name: str, complexity: str) -> PackageSpec:
    """Translate the USG ``SMS`` solver onto MODFLOW 6 ``IMS``.

    Only the convergence criteria and iteration caps carry over literally. The
    rest of SMS -- its under-relaxation flavour, backtracking and linear
    acceleration -- has no term-for-term MODFLOW 6 equivalent, so a complexity
    preset supplies those and the result is a *comparable*, not identical,
    solver.

    The preset defaults to ``COMPLEX`` rather than ``MODERATE`` deliberately. A
    USG model that reached for SMS with a Newton formulation is a hard nonlinear
    problem -- convertible layers, wetting and drying, thin cells -- and on the
    Ten Trails model ``MODERATE`` does not merely converge poorly: MODFLOW 6
    dies with SIGFPE during the first solve. ``COMPLEX`` runs the same model to
    normal termination, with or without the SMS-derived tolerances.
    """

    sms = model.sms
    if sms is None:
        return api.ims(models=(name,), complexity=complexity)
    # Carry only what means the same thing in both solvers: the OUTER (nonlinear)
    # head-change closure and the iteration caps. SMS's HICLOSE is not MODFLOW 6's
    # INNER_DVCLOSE -- the inner solvers differ, and SMS's remaining inner controls
    # (IACL, NORDER, LEVEL, NORTH, RCLOSEPCGU) have no MODFLOW 6 counterpart at
    # all. Copying half of SMS's numbers on top of half a preset produces a solver
    # that is neither, so the preset owns the inner solve outright.
    values: dict[str, Any] = {
        "complexity": complexity,
        "outer_dvclose": sms.hclose,
        "outer_maximum": sms.mxiter,
        "inner_maximum": sms.iter1,
    }
    if sms.linear_acceleration:
        values["linear_acceleration"] = sms.linear_acceleration
    if sms.is_newton and sms.theta is not None:
        # SMS's delta-bar-delta under-relaxation and backtracking map field for
        # field onto MODFLOW 6's, and they are not decoration: they are how the
        # original model was made to converge. Dropping them and keeping only the
        # tolerances leaves MODFLOW 6 taking 100,000-foot Newton steps on a model
        # whose heads span 400 feet.
        values.update(
            under_relaxation="DBD",
            under_relaxation_theta=sms.theta,
            under_relaxation_kappa=sms.akappa,
            under_relaxation_gamma=sms.gamma,
            under_relaxation_momentum=sms.amomentum,
        )
        if sms.numtrack:
            values.update(
                backtracking_number=sms.numtrack,
                backtracking_tolerance=sms.btol,
                backtracking_reduction_factor=sms.breduc,
                backtracking_residual_limit=sms.reslim,
            )
    return api.ims(models=(name,), **values)


# ---------------------------------------------------------------- report


def conversion_notes(model: UsgModel) -> str:
    """Return the human-readable account of the conversion. See :meth:`UsgModel.report`."""

    lines: list[str] = []
    add = lines.append

    add(f"MODFLOW-USG model: {model.name_file.path.name}")
    add(f"  grid       {model.nlay} layers x {model.ncpl} cells = {model.nodes} nodes")
    add(f"  time       {model.nper} stress periods, {model.disu.perlen.sum():.0f} {model.time_units}")
    add(f"  units      length={model.length_units}, time={model.time_units}")
    if model.length_units == "unknown":
        add("             ^ DISU LENUNI is 0; MODFLOW 6 will be told 'unknown' unless you say otherwise")
    if model.start_date_time:
        add(f"  starts     {model.start_date_time} (inferred from the DISU period labels)")
    active = int((model.idomain != 0).sum())
    add(f"  active     {active} of {model.idomain.size} cells")
    add("")

    add("CONVERTED")
    add("  DISU  -> DISV      layered grid confirmed; geometry from the .gsf")
    add(f"  BAS6  -> IDOMAIN+IC  starting heads {model.strt.min():.1f}-{model.strt.max():.1f}")
    convertible = "convertible" if model.icelltype.any() else "confined"
    add(f"  LPF   -> NPF+STO   icelltype {convertible}; K {model.k.min():.3g}-{model.k.max():.3g}")
    if model.sms:
        add("  SMS   -> IMS       approximate: criteria carry over, solver internals do not")
    for ftype in ("CHD", "DRN", "GHB", "WEL", "RIV"):
        records = model.boundaries.get(ftype)
        if records is None:
            continue
        if records.n_records == 0 and records.n_cln_records:
            add(f"  {ftype:<5} -> (none)     every record is on a CLN node; see below")
            continue
        note = f"{records.n_defined} period(s) defined, {len(records.reused)} reused"
        switched_off = sum(1 for block in records.periods.values() if len(block) == 0)
        if switched_off:
            note += f", {switched_off} with no groundwater record (switched off)"
        add(f"  {ftype:<5} -> {ftype:<9} {records.n_records} records, {note}")
    if model.rch is not None:
        add(f"  RCH   -> RCHA      {len(model.rch.rech)} arrays, NRCHOP={model.rch.nrchop}")
    if model.ets is not None:
        add(f"  ETS   -> EVT       segmented (NETSEG={model.ets.netseg}), list-based")
    if model.hfb.size:
        add(f"  HFB6  -> HFB       {len(model.hfb)} barriers")
    add("")

    findings = mf6_findings(model)
    if findings:
        # Not every finding is a refusal. Whether MODFLOW 6 rejects a below-bottom
        # elevation or accepts it and quietly changes what the boundary does depends on
        # the layer's icelltype (ledger 291), so a single "IS STRICTER" heading would be
        # wrong half the time; the marker on each line says which this one is.
        add("MODFLOW 6 DIFFERS")
        add("  REJECTS = will not start; allows = starts, and the boundary behaves differently")
        for finding in findings:
            marker = "REJECTS" if finding.blocks_run else "allows "
            add(f"  {marker}  {finding.detail}")
            add(f"           fix: {finding.fix}")
        add("")

    approximations, omissions = _caveats(model)
    if approximations:
        add("APPROXIMATED")
        for line in approximations:
            add(f"  - {line}")
        add("")
    if omissions:
        add("NOT CONVERTED")
        for line in omissions:
            add(f"  - {line}")
        add("")
    return "\n".join(lines)


def _caveats(model: UsgModel) -> tuple[list[str], list[str]]:
    """Return the approximations and the omissions, each as report lines."""

    approximations: list[str] = []
    omissions: list[str] = []

    chd = model.boundaries.get("CHD")
    if chd is not None and chd.periods:
        drift = [
            float(np.abs(block[:, 2] - block[:, 1]).max())
            for block in chd.periods.values()
            if block.shape[1] > 2 and len(block)
        ]
        if drift and max(drift) > 0:
            approximations.append(
                f"CHD ramps from shead to ehead within each period (max {max(drift):.3f} "
                "of head); MODFLOW 6 holds one value, so the end-of-period head is used"
            )

    if model.ets is not None:
        if model.ets.netsop == 3:
            approximations.append(
                "ETS applied to the highest ACTIVE cell (NETSOP=3), which MODFLOW 6 cannot "
                "express; resolved once from IBOUND, which is static -- so this is exact"
            )
        if model.ets.netseg > 1:
            columns = int(model.has_active_column.sum())
            approximations.append(
                f"ETS is segmented (NETSEG={model.ets.netseg}), and MODFLOW 6 cannot combine "
                f"segments with array input, so EVT becomes list-based: about "
                f"{columns * model.nper:,} records"
            )

    if model.sms is not None:
        approximations.append(
            "SMS -> IMS keeps the convergence criteria; under-relaxation, backtracking and "
            "linear acceleration have no term-for-term equivalent and come from the preset"
        )

    if "DRAWDOWN SAVE" in model.oc_words:
        omissions.append("OC requested DRAWDOWN output, which MODFLOW 6 does not produce")

    if model.cln is not None:
        cells = model.cln.cells_touched(model.ncpl)
        bodies = model.cln.waterbodies
        streams = model.cln.streams
        omissions.append(
            f"CLN ({model.cln.nclnnds} nodes) has no MODFLOW 6 counterpart and is NOT "
            f"converted. It touches {cells.size} groundwater cells in "
            f"{len(model.cln.features)} features "
            f"({len(bodies)} waterbody, {len(streams)} stream). "
            "Rebuild them as LAK/SFR from model.cln_polygons()."
        )
        for feature in model.cln.features:
            omissions.append(f"    [{feature.index}] {feature.describe()}")

    for ftype, records in model.boundaries.items():
        if records.n_cln_records:
            omissions.append(
                f"{ftype}: {records.n_cln_records} records address CLN nodes, not groundwater "
                f"cells, and are dropped with the CLN network"
            )
    return approximations, omissions
