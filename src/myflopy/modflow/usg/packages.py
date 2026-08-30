"""Per-package readers for MODFLOW-USG input.

Each function takes a resolved :class:`~myflopy.modflow.usg._io.NameFile` and
returns plain data -- arrays, dicts and dataclasses -- with no MODFLOW 6 opinion
baked in. The MF6 mapping lives in :mod:`myflopy.modflow.usg.convert`, so that
the reader stays usable for anything else you would want to do with a USG model
(re-gridding it, plotting it, comparing two of them).

``DISU`` is read through FloPy's ``flopy.mfusg`` loader, which handles it
correctly. ``BAS6`` and ``LPF`` are read here instead, because FloPy 3.10 cannot
load them when starting heads come from a binary unit -- see
:mod:`myflopy.modflow.usg._io` for the specifics.

**Period-block discrimination.** A USG list boundary condition alternates a
period header with the records for that period, and neither is self-describing:
a header may carry trailing free text (``71  Stress Period 11  10  October``) or
trailing numbers, so counting fields cannot tell them apart. What does tell them
apart is that every record here is an integer node followed by *floating point*
values, while a header is integers (and words) only. That is the discriminator
:func:`_looks_like_record` applies, and it is why the readers do not trust
``ITMP`` for the record count -- only for the reuse decision.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from myflopy._logging import get_logger
from myflopy.modflow.usg._io import ArrayCursor, NameFile, free_floats, free_ints

logger = get_logger(__name__)

__all__ = [
    "BasData",
    "DisuData",
    "EtsData",
    "LpfData",
    "PeriodRecords",
    "RchData",
    "SmsData",
    "read_bas",
    "read_disu",
    "read_ets",
    "read_hfb",
    "read_list_bc",
    "read_lpf",
    "read_rch",
    "read_sms",
]

#: MODFLOW time-unit codes (DISU ITMUNI).
TIME_UNITS = {0: "unknown", 1: "seconds", 2: "minutes", 3: "hours", 4: "days", 5: "years"}

#: MODFLOW length-unit codes (DISU LENUNI). 0 means the model never said.
LENGTH_UNITS = {0: "unknown", 1: "feet", 2: "meters", 3: "centimeters"}


@dataclass(slots=True)
class DisuData:
    """Discretization and timing read from a USG ``DISU``."""

    nodes: int
    nlay: int
    njag: int
    ivsd: int
    nper: int
    itmuni: int
    lenuni: int
    nodelay: np.ndarray
    top: np.ndarray
    bot: np.ndarray
    area: np.ndarray
    iac: np.ndarray
    ja: np.ndarray
    perlen: np.ndarray
    nstp: np.ndarray
    tsmult: np.ndarray
    steady: np.ndarray

    @property
    def ncpl(self) -> int:
        """Nodes per layer; meaningful only when the grid is layered."""

        return int(self.nodelay[0])

    @property
    def is_layered(self) -> bool:
        """True when every layer holds the same number of nodes."""

        return bool(np.all(self.nodelay == self.nodelay[0]))

    @property
    def time_units(self) -> str:
        """The DISU time unit as a MODFLOW 6 word."""

        return TIME_UNITS.get(self.itmuni, "unknown")

    @property
    def length_units(self) -> str:
        """The DISU length unit as a MODFLOW 6 word."""

        return LENGTH_UNITS.get(self.lenuni, "unknown")


@dataclass(slots=True)
class BasData:
    """``BAS6`` contents: the active-cell mask and the starting heads."""

    ibound: np.ndarray
    hnoflo: float
    strt: np.ndarray
    options: tuple[str, ...] = ()


@dataclass(slots=True)
class LpfData:
    """``LPF`` contents: flow properties, one row per layer."""

    ilpfcb: int
    hdry: float
    options: tuple[str, ...]
    laytyp: np.ndarray
    layavg: np.ndarray
    chani: np.ndarray
    layvka: np.ndarray
    laywet: np.ndarray
    hk: np.ndarray
    vka: np.ndarray
    ss: np.ndarray
    sy: np.ndarray

    @property
    def is_convertible(self) -> np.ndarray:
        """Per-layer convertible flag; anything non-zero is water-table."""

        return self.laytyp != 0


@dataclass(slots=True)
class SmsData:
    """``SMS`` solver settings, as read.

    The nonlinear block (item 2, present when ``NONLINMETH`` is non-zero) is the
    part worth carrying: SMS's delta-bar-delta under-relaxation and its
    backtracking controls have exact MODFLOW 6 counterparts, and they are the
    tuning the original modeller did to make a hard unconfined model converge.
    """

    hclose: float
    hiclose: float
    mxiter: int
    iter1: int
    iprsms: int
    nonlinmeth: int
    linmeth: int
    options: tuple[str, ...] = ()
    extra: tuple[float, ...] = ()
    # item 2 -- delta-bar-delta under-relaxation and backtracking
    theta: float | None = None
    akappa: float | None = None
    gamma: float | None = None
    amomentum: float | None = None
    numtrack: int | None = None
    btol: float | None = None
    breduc: float | None = None
    reslim: float | None = None
    # item 3 -- linear solver
    iacl: int | None = None

    @property
    def is_newton(self) -> bool:
        """True when SMS was configured for a Newton (non-Picard) formulation."""

        return self.nonlinmeth != 0

    @property
    def linear_acceleration(self) -> str | None:
        """MODFLOW 6's word for SMS's ``IACL`` linear acceleration choice."""

        return {0: "CG", 1: "BICGSTAB", 2: "BICGSTAB"}.get(self.iacl)


@dataclass(slots=True)
class PeriodRecords:
    """Per-period records of one list boundary condition.

    ``periods`` maps a zero-based stress period to its records. A period that
    reused the previous one (``ITMP < 0``) is *not* a key -- the reuse is
    recorded in :attr:`reused` instead, because MODFLOW 6 expresses reuse the
    same way, by omitting the period from its period dict.
    """

    ftype: str
    periods: dict[int, np.ndarray] = field(default_factory=dict)
    cln_periods: dict[int, np.ndarray] = field(default_factory=dict)
    reused: tuple[int, ...] = ()
    options: tuple[str, ...] = ()
    maxbound: int = 0

    @property
    def n_defined(self) -> int:
        """How many periods carry their own records."""

        return len(self.periods)

    @property
    def n_records(self) -> int:
        """Total GWF records across every period that defines its own."""

        return int(sum(len(v) for v in self.periods.values()))

    @property
    def n_cln_records(self) -> int:
        """Total records addressed to CLN nodes across all periods."""

        return int(sum(len(v) for v in self.cln_periods.values()))


@dataclass(slots=True)
class RchData:
    """``RCH`` contents: one recharge array per period that defines a new one."""

    nrchop: int
    irchcb: int
    rech: dict[int, np.ndarray] = field(default_factory=dict)
    reused: tuple[int, ...] = ()


@dataclass(slots=True)
class EtsData:
    """``ETS`` (segmented evapotranspiration) contents.

    ``pxdp``/``petm`` hold ``netseg - 1`` arrays each for the periods that
    define them -- the segment breakpoints as a fraction of extinction depth and
    the corresponding fraction of the maximum rate.
    """

    netsop: int
    ietscb: int
    npets: int
    netseg: int
    surf: dict[int, np.ndarray] = field(default_factory=dict)
    rate: dict[int, np.ndarray] = field(default_factory=dict)
    depth: dict[int, np.ndarray] = field(default_factory=dict)
    pxdp: dict[int, list[np.ndarray]] = field(default_factory=dict)
    petm: dict[int, list[np.ndarray]] = field(default_factory=dict)
    ietsl: dict[int, np.ndarray] = field(default_factory=dict)


# ---------------------------------------------------------------- DISU


def read_disu(name_file: NameFile) -> DisuData:
    """Read the ``DISU`` package via FloPy's USG loader.

    FloPy reads DISU correctly, so this delegates rather than reimplementing
    ~400 lines of array plumbing. Only BAS/LPF need a hand-rolled reader.
    """

    from flopy.mfusg import MfUsg

    entry = name_file.package("DISU")
    if entry is None:
        raise ValueError(f"{name_file.path.name}: no DISU package -- not a MODFLOW-USG model")

    model = MfUsg.load(
        name_file.path.name,
        model_ws=str(name_file.workspace),
        load_only=["DISU"],
        check=False,
        verbose=False,
    )
    disu = model.disu
    arr = lambda x: np.asarray(x.array if hasattr(x, "array") else x)  # noqa: E731

    return DisuData(
        nodes=int(disu.nodes),
        nlay=int(disu.nlay),
        njag=int(disu.njag),
        ivsd=int(disu.ivsd),
        nper=int(disu.nper),
        itmuni=int(disu.itmuni),
        lenuni=int(disu.lenuni),
        nodelay=arr(disu.nodelay).astype(int),
        top=arr(disu.top).astype(float),
        bot=arr(disu.bot).astype(float),
        area=arr(disu.area).astype(float),
        iac=arr(disu.iac).astype(int),
        ja=arr(disu.ja).astype(int),
        perlen=arr(disu.perlen).astype(float),
        nstp=arr(disu.nstp).astype(int),
        tsmult=arr(disu.tsmult).astype(float),
        steady=arr(disu.steady).astype(bool),
    )


def layered_report(disu: DisuData) -> tuple[bool, str]:
    """Classify a DISU as layered (DISV-convertible) or not, and say why.

    A layered grid is one where every layer holds the same nodes in the same
    plan positions, so node ``n`` in layer ``k`` sits directly above node ``n``
    in layer ``k+1``. MODFLOW 6's DISV requires exactly that. The test is on the
    connectivity rather than on ``IVSD``, because ``IVSD`` states the modeller's
    intent while the connections state what was actually built: every connection
    must be either same-layer or a strict vertical step of ``ncpl``.
    """

    if not disu.is_layered:
        return False, (
            "layers hold different node counts "
            f"({', '.join(str(int(n)) for n in disu.nodelay[:6])}...)"
        )

    ncpl = disu.ncpl
    iac, ja = disu.iac, disu.ja
    ptr = np.concatenate([[0], np.cumsum(iac)])
    lateral = vertical = other = 0
    for node in range(disu.nodes):
        neighbours = ja[ptr[node] + 1 : ptr[node + 1]] - 1
        same_layer = (neighbours // ncpl) == (node // ncpl)
        step = np.abs(neighbours - node)
        lateral += int(same_layer.sum())
        vertical += int((step == ncpl).sum())
        other += int((~same_layer & (step != ncpl)).sum())

    if other:
        return False, (
            f"{other} of {lateral + vertical + other} connections are neither "
            "same-layer nor a strict vertical step (nested or ghost-node refinement)"
        )
    return True, f"{lateral} lateral + {vertical} vertical connections, none irregular"


# ---------------------------------------------------------------- BAS6


def read_bas(name_file: NameFile, disu: DisuData) -> BasData:
    """Read ``BAS6``: options line, IBOUND per layer, HNOFLO, then STRT per layer."""

    entry = name_file.package("BAS6") or name_file.package("BAS")
    if entry is None:
        raise ValueError(f"{name_file.path.name}: no BAS6 package")

    cursor = ArrayCursor(entry.path, name_file)
    options = tuple(w.upper() for w in cursor.next_line().split())
    counts = [int(n) for n in disu.nodelay]

    ibound = np.concatenate([cursor.read_array(n, int) for n in counts])
    hnoflo = cursor.read_floats()[0]
    strt = np.concatenate([cursor.read_array(n, float) for n in counts])

    logger.debug(
        "BAS6: %d active of %d nodes, HNOFLO=%g", int((ibound != 0).sum()), ibound.size, hnoflo
    )
    return BasData(ibound=ibound, hnoflo=hnoflo, strt=strt, options=options)


# ---------------------------------------------------------------- LPF


def read_lpf(name_file: NameFile, disu: DisuData) -> LpfData:
    """Read ``LPF``: the header, five per-layer flag rows, then the property arrays.

    Which property arrays exist per layer depends on the flags: ``LAYVKA`` picks
    whether the vertical entry is a conductivity or a ratio, and storage arrays
    appear only when at least one period is transient.
    """

    entry = name_file.package("LPF")
    if entry is None:
        raise ValueError(f"{name_file.path.name}: no LPF package")

    cursor = ArrayCursor(entry.path, name_file)
    header = cursor.next_line()
    numbers = free_floats(header)
    ilpfcb, hdry = int(numbers[0]), numbers[1]
    options = tuple(w.upper() for w in header.split() if not _is_number(w))

    nlay = disu.nlay
    laytyp = np.asarray(cursor.read_ints(nlay))
    layavg = np.asarray(cursor.read_ints(nlay))
    chani = np.asarray(cursor.read_floats(nlay))
    layvka = np.asarray(cursor.read_ints(nlay))
    laywet = np.asarray(cursor.read_ints(nlay))

    transient = bool(np.any(~disu.steady))
    counts = [int(n) for n in disu.nodelay]
    hk, vka, ss, sy = [], [], [], []
    for layer in range(nlay):
        n = counts[layer]
        hk.append(cursor.read_array(n, float))
        if chani[layer] <= 0:
            cursor.read_array(n, float)  # HANI, unused: MF6 expresses this as k22
        vka.append(cursor.read_array(n, float))
        if transient:
            ss.append(cursor.read_array(n, float))
            if laytyp[layer] != 0:
                sy.append(cursor.read_array(n, float))
            else:
                sy.append(np.zeros(n))
        if laywet[layer] != 0 and laytyp[layer] != 0:
            cursor.read_array(n, float)  # WETDRY, no MF6 equivalent (Newton replaces it)

    return LpfData(
        ilpfcb=ilpfcb,
        hdry=hdry,
        options=options,
        laytyp=laytyp,
        layavg=layavg,
        chani=chani,
        layvka=layvka,
        laywet=laywet,
        hk=np.concatenate(hk),
        vka=np.concatenate(vka),
        ss=np.concatenate(ss) if ss else np.zeros(disu.nodes),
        sy=np.concatenate(sy) if sy else np.zeros(disu.nodes),
    )


# ---------------------------------------------------------------- SMS


def read_sms(name_file: NameFile) -> SmsData:
    """Read the ``SMS`` solver package header."""

    entry = name_file.package("SMS")
    if entry is None:
        raise ValueError(f"{name_file.path.name}: no SMS package")

    cursor = ArrayCursor(entry.path, name_file)
    line = cursor.next_line()
    numbers = free_floats(line)
    options = tuple(w.upper() for w in line.split() if not _is_number(w))
    extra: list[float] = []
    while not cursor.at_end():
        extra.extend(free_floats(cursor.next_line()))

    nonlinmeth = int(numbers[5])
    data = SmsData(
        hclose=numbers[0],
        hiclose=numbers[1],
        mxiter=int(numbers[2]),
        iter1=int(numbers[3]),
        iprsms=int(numbers[4]),
        nonlinmeth=nonlinmeth,
        linmeth=int(numbers[6]),
        options=options,
        extra=tuple(extra),
    )

    # Item 2 follows only for a non-Picard formulation:
    # THETA AKAPPA GAMMA AMOMENTUM NUMTRACK BTOL BREDUC RESLIM
    if nonlinmeth != 0 and len(extra) >= 8:
        (
            data.theta,
            data.akappa,
            data.gamma,
            data.amomentum,
            numtrack,
            data.btol,
            data.breduc,
            data.reslim,
        ) = extra[:8]
        data.numtrack = int(numtrack)
        # Item 3 (linear) starts after item 2; its first field is IACL.
        if len(extra) >= 9:
            data.iacl = int(extra[8])
    elif nonlinmeth == 0 and extra:
        data.iacl = int(extra[0])
    return data


# ---------------------------------------------------------------- list BCs


def _is_number(token: str) -> bool:
    """True when ``token`` parses as a MODFLOW free-format number."""

    try:
        float(token.replace("D", "E").replace("d", "e"))
    except ValueError:
        return False
    return True


def _looks_like_record(line: str) -> bool:
    """True when ``line`` is a boundary record rather than a period header.

    Every list record is an integer node followed by at least one *floating
    point* value; a period header is integers and words only. That distinction
    survives the trailing free text and trailing period numbers that real files
    put on their header lines.
    """

    tokens = line.split()
    if len(tokens) < 2 or not _is_number(tokens[0]):
        return False
    return any(_is_number(t) and ("." in t or "e" in t.lower()) for t in tokens[1:])


def read_list_bc(
    name_file: NameFile,
    ftype: str,
    *,
    nper: int,
    n_gwf_nodes: int,
    n_values: int,
) -> PeriodRecords:
    """Read a USG list boundary condition into per-period record arrays.

    Parameters
    ----------
    ftype
        Package type in the name file, e.g. ``"DRN"``.
    nper
        Stress periods the model actually runs. Files routinely carry MORE
        period blocks than this -- the Ten Trails WEL holds 612 and its CHD 792
        for a 72-period run -- and MODFLOW reads the first ``nper`` and ignores
        the rest, so this truncates rather than trusting the file's length.
    n_gwf_nodes
        Node count of the groundwater grid. Records above it address CLN nodes
        and are separated into :attr:`PeriodRecords.cln_periods`.
    n_values
        Numeric fields per record *including* the node, e.g. 3 for ``DRN``
        (node, elevation, conductance).

    Returns
    -------
    PeriodRecords
    """

    entry = name_file.package(ftype)
    if entry is None:
        raise ValueError(f"{name_file.path.name}: no {ftype} package")

    cursor = ArrayCursor(entry.path, name_file)
    header = cursor.next_line()
    options = tuple(w.upper() for w in header.split() if not _is_number(w))
    maxbound = free_ints(header)[0] if free_ints(header) else 0

    periods: dict[int, np.ndarray] = {}
    cln_periods: dict[int, np.ndarray] = {}
    reused: list[int] = []

    for period in range(nper):
        if cursor.at_end():
            break
        head_line = cursor.next_line()
        itmp = free_ints(head_line)
        if itmp and itmp[0] < 0:
            reused.append(period)
            continue

        rows: list[list[float]] = []
        while not cursor.at_end() and _looks_like_record(cursor.peek()):
            rows.append(free_floats(cursor.next_line())[:n_values])
            cursor.pos += 0  # records consumed by peek/next_line above

        if not rows:
            periods[period] = np.empty((0, n_values))
            continue

        block = np.asarray(rows, dtype=float)
        nodes = block[:, 0].astype(np.int64)
        is_cln = nodes > n_gwf_nodes
        # A package whose header declares a separate CLN count (WEL does) puts
        # CLN wells in the same block but numbers them from 1 within the CLN
        # grid, so "above n_gwf_nodes" does not catch them. The header's third
        # count is what says how many trailing records are CLN.
        counts = [n for n in free_ints(head_line)]
        if len(counts) >= 3 and counts[0] == 0 and counts[2] == len(rows):
            cln_periods[period] = block
            periods[period] = np.empty((0, n_values))
            continue

        periods[period] = block[~is_cln]
        if is_cln.any():
            cln_periods[period] = block[is_cln]

    logger.debug(
        "%s: %d periods defined, %d reused, %d GWF records, %d CLN records",
        ftype,
        len(periods),
        len(reused),
        sum(len(v) for v in periods.values()),
        sum(len(v) for v in cln_periods.values()),
    )
    return PeriodRecords(
        ftype=ftype.upper(),
        periods=periods,
        cln_periods=cln_periods,
        reused=tuple(reused),
        options=options,
        maxbound=maxbound,
    )


def read_hfb(name_file: NameFile) -> np.ndarray:
    """Read ``HFB6`` barriers as an ``(n, 3)`` array of ``node1, node2, hydchr``."""

    entry = name_file.package("HFB6") or name_file.package("HFB")
    if entry is None:
        return np.empty((0, 3))

    cursor = ArrayCursor(entry.path, name_file)
    counts = free_ints(cursor.next_line())
    nhfbnp = counts[2] if len(counts) >= 3 else counts[-1]
    rows = [free_floats(cursor.next_line())[:3] for _ in range(nhfbnp)]
    return np.asarray(rows, dtype=float)


# ---------------------------------------------------------------- RCH / ETS


def read_rch(name_file: NameFile, *, nper: int, ncpl: int) -> RchData:
    """Read ``RCH``: a header, then per period a flag and (if non-negative) an array."""

    entry = name_file.package("RCH")
    if entry is None:
        raise ValueError(f"{name_file.path.name}: no RCH package")

    cursor = ArrayCursor(entry.path, name_file)
    header = free_ints(cursor.next_line())
    nrchop, irchcb = header[0], (header[1] if len(header) > 1 else 0)

    rech: dict[int, np.ndarray] = {}
    reused: list[int] = []
    for period in range(nper):
        if cursor.at_end():
            break
        inrech = free_ints(cursor.next_line())[0]
        if inrech < 0:
            reused.append(period)
            continue
        rech[period] = cursor.read_array(ncpl, float)

    logger.debug("RCH: %d periods with a new array, %d reused", len(rech), len(reused))
    return RchData(nrchop=nrchop, irchcb=irchcb, rech=rech, reused=tuple(reused))


def read_ets(name_file: NameFile, *, nper: int, ncpl: int) -> EtsData:
    """Read ``ETS`` (segmented evapotranspiration).

    Per period the flags are ``INETSS INETSR INETSX [INIETS] [INSGDF]``: the ET
    surface, the maximum rate, the extinction depth, the layer indicator (only
    when ``NETSOP == 2``) and the segment definition (only when
    ``NETSEG > 1``). A negative flag reuses the previous period's array.
    """

    entry = name_file.package("ETS")
    if entry is None:
        raise ValueError(f"{name_file.path.name}: no ETS package")

    cursor = ArrayCursor(entry.path, name_file)
    header = free_ints(cursor.next_line())
    netsop, ietscb = header[0], header[1]
    npets = header[2] if len(header) > 2 else 0
    netseg = header[3] if len(header) > 3 else 1

    data = EtsData(netsop=netsop, ietscb=ietscb, npets=npets, netseg=netseg)
    for period in range(nper):
        if cursor.at_end():
            break
        flags = free_ints(cursor.next_line())
        inetss, inetsr, inetsx = flags[0], flags[1], flags[2]
        index = 3
        iniets = -1
        if netsop == 2:
            iniets = flags[index] if len(flags) > index else -1
            index += 1
        insgdf = flags[index] if len(flags) > index else -1

        if inetss >= 0:
            data.surf[period] = cursor.read_array(ncpl, float)
        if inetsr >= 0:
            data.rate[period] = cursor.read_array(ncpl, float)
        if inetsx >= 0:
            data.depth[period] = cursor.read_array(ncpl, float)
        if netsop == 2 and iniets >= 0:
            data.ietsl[period] = cursor.read_array(ncpl, int)
        if netseg > 1 and insgdf >= 0:
            data.pxdp[period] = [cursor.read_array(ncpl, float) for _ in range(netseg - 1)]
            data.petm[period] = [cursor.read_array(ncpl, float) for _ in range(netseg - 1)]

    logger.debug(
        "ETS: netsop=%d netseg=%d; %d surf / %d rate / %d depth arrays over %d periods",
        netsop,
        netseg,
        len(data.surf),
        len(data.rate),
        len(data.depth),
        nper,
    )
    return data


def read_oc_words(name_file: NameFile) -> tuple[str, ...]:
    """Return the distinct output-control words a USG ``OC`` file asks for."""

    entry = name_file.package("OC")
    if entry is None:
        return ()
    words: list[str] = []
    for raw in Path(entry.path).read_text(errors="replace").splitlines():
        upper = raw.upper()
        for word in ("HEAD SAVE", "DRAWDOWN SAVE", "SAVE BUDGET", "PRINT BUDGET",
                     "SAVE HEAD", "PRINT HEAD", "COMPACT BUDGET"):
            if word in upper and word not in words:
                words.append(word)
    return tuple(words)
