"""Semantic helpers for building MODFLOW 6 MVR package specifications."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

from myflopy.advanced import mvr_spec
from myflopy.specs import PackageSpec


@dataclass(frozen=True, slots=True)
class MoverConnection:
    """One MVR endpoint: a feature in an advanced package that can give/take water.

    Identifies a single provider or receiver for a :class:`Move` by package name
    and zero-based feature id (an SFR reach number, a LAK lake number, a UZF cell
    number, etc.). Construct it directly when you know the id, or -- for SFR/LAK --
    let myflopy resolve it from geometry via ``mf.sfr_connection(...)`` /
    ``mf.lak_connection(...)``.

    Parameters
    ----------
    package : str
        The name the package was declared under on the model (``name=`` of
        ``mf.sfr``/``mf.lak``/``mf.uzf``/..., ``"sfr"``/``"lak"``/``"uzf"`` by
        default). It is lowercased; an empty name raises ``ValueError``. The
        :class:`MVRBuilder` spec ``requires=`` it, so a name that matches no package
        declared before the mover is refused (``ValueError``) when the model spec
        is created, or at build for library references. The package must also be
        declared with ``mover=True`` for MF6 to move water through it; myflopy
        does not check that.
    index : int
        Zero-based feature number within that package -- an SFR reach, a LAK lake,
        a UZF object (FloPy adds one on write). Converted with ``int()``; a negative
        value raises ``ValueError``.

    Examples
    --------
    >>> MoverConnection("sfr", 0)          # the first SFR reach
    >>> mf.sfr_connection(sfr, "main_stem")  # the stream outlet, resolved for you
    """

    package: str
    index: int

    def __post_init__(self) -> None:
        """Lowercase the package name and validate a non-negative zero-based ``index``."""

        if not self.package:
            raise ValueError("package is required.")
        if int(self.index) < 0:
            raise ValueError("index must be zero-based and non-negative.")
        object.__setattr__(self, "package", str(self.package).lower())
        object.__setattr__(self, "index", int(self.index))


@dataclass(frozen=True, slots=True)
class Move:
    """Route water from one advanced-package feature to another via MVR.

    A move takes available water leaving the ``source`` feature and delivers it to
    the ``receiver`` -- e.g. a stream's outflow into a lake, or rejected UZF
    infiltration into a stream. Pass moves to ``mf.mvr(moves=(...))``.

    Parameters
    ----------
    source : MoverConnection
        The provider: the feature whose available water is moved. Anything else
        raises ``ValueError``.
    receiver : MoverConnection
        The feature that receives the moved water. Anything else raises
        ``ValueError``.
    method : {"FACTOR", "EXCESS", "THRESHOLD", "UPTO"}, default "FACTOR"
        MF6's ``MVRTYPE``, how ``value`` is read; case-insensitive (stored
        uppercased), and any other word raises ``ValueError``. ``"FACTOR"`` moves
        that fraction of the available water; ``"EXCESS"`` moves the water above
        ``value``; ``"THRESHOLD"`` moves ``value`` once at least that much is
        available and nothing below it; ``"UPTO"`` moves ``value``, or all of the
        available water when there is less (MF6's MVR rules, which mirror SFR's
        ``cprior``).
    value : float, default 1.0
        The factor (dimensionless, for ``"FACTOR"``) or rate (model length^3/time,
        for the other three) that ``method`` applies. Converted with ``float()``.
        The default with ``"FACTOR"`` moves all available water.

    Examples
    --------
    >>> Move(MoverConnection("sfr", 5), MoverConnection("lak", 0))          # 100% of reach 5 -> lake 0
    >>> Move(mf.sfr_connection(sfr, "main_stem"),
    ...      mf.lak_connection(lak, "valley_lake"), value=0.5)              # half the outflow
    """

    source: MoverConnection
    receiver: MoverConnection
    method: str = "FACTOR"
    value: float = 1.0

    def __post_init__(self) -> None:
        """Validate the endpoints are connections and uppercase/validate the transfer ``method``."""

        if not isinstance(self.source, MoverConnection):
            raise ValueError("source must be a MoverConnection.")
        if not isinstance(self.receiver, MoverConnection):
            raise ValueError("receiver must be a MoverConnection.")
        method = self.method.upper()
        if method not in {"FACTOR", "EXCESS", "THRESHOLD", "UPTO"}:
            raise ValueError("method must be FACTOR, EXCESS, THRESHOLD, or UPTO.")
        object.__setattr__(self, "method", method)
        object.__setattr__(self, "value", float(self.value))

    def to_record(self) -> list[Any]:
        """Return one FloPy/MODFLOW MVR perioddata record."""

        return [
            self.source.package,
            self.source.index,
            self.receiver.package,
            self.receiver.index,
            self.method,
            self.value,
        ]


@dataclass(frozen=True, slots=True)
class MVRBuilder:
    """Engine that turns semantic :class:`Move` declarations into an MF6 MVR package.

    Collects the provider/receiver packages from the moves, validates them, and
    emits a :class:`~myflopy.specs.PackageSpec` with the MVR packagedata +
    perioddata. ``moves`` may be a flat sequence (applied every period) or a
    ``{period: [Move, ...]}`` mapping.

    This is the engine under the package-first ``mf.mvr(...)`` facade -- prefer that
    front door for new work. The moved packages must be declared on the model with
    ``mover=True`` and ordered before the mover.

    Parameters
    ----------
    nper : int
        Number of stress periods; a value below 1 raises ``ValueError``. Every
        period ``0 .. nper - 1`` gets a ``PERIOD`` block, empty when no move is
        declared for it.
    moves : sequence of Move, or dict of int to sequence of Move
        The transfers. A flat sequence is applied in every period. A mapping is
        keyed by zero-based period, and a period it leaves out gets NO moves (an
        empty block, not a repeat of the previous one). A key outside
        ``range(nper)``, an item that is not a :class:`Move`, and no move at all
        each raise ``ValueError`` at :meth:`build`, not here. The MVR ``packages``
        block is derived from the moves, in order of first appearance.
    name : str, default "mvr"
        Package name: FloPy's ``pname``, the input file ``<model>.<name>``, the
        output files ``<model>.<name>.bud`` and ``<model>.<name>.csv``, and the
        built spec's name.
    print_input : bool, default False
        Write MF6's ``PRINT_INPUT`` option (echo the mover records to the listing
        file).
    print_flows : bool, default False
        Write MF6's ``PRINT_FLOWS`` option (print the moved flow rates to the
        listing file on the time steps where Output Control prints the budget).
    options : Mapping[str, Any], optional
        Further keywords for :func:`~myflopy.advanced.mvr_spec`, which passes any
        it does not take itself (``maxmvr``/``maxpackages`` it does) on to
        :class:`flopy.mf6.ModflowGwfmvr`. Checked against that class's signature at
        :meth:`build`; an unknown name raises ``TypeError``. ``modelnames=True``
        raises ``ValueError`` there: a model's mover takes no MF6 ``MODELNAMES``
        (ledger 413). Defaults to none.
    """

    nper: int
    moves: Mapping[int, Sequence[Move]] | Sequence[Move]
    name: str = "mvr"
    print_input: bool = False
    print_flows: bool = False
    options: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate that ``nper`` is at least 1."""

        if self.nper < 1:
            raise ValueError("nper must be at least 1.")

    @property
    def perioddata(self) -> dict[int, list[list[Any]]]:
        """Return normalized MVR perioddata for every stress period."""

        if isinstance(self.moves, Mapping):
            by_period = {int(period): tuple(moves) for period, moves in self.moves.items()}
        else:
            by_period = {period: tuple(self.moves) for period in range(self.nper)}

        result: dict[int, list[list[Any]]] = {period: [] for period in range(self.nper)}
        for period, moves in by_period.items():
            if period not in result:
                raise ValueError(f"MVR period {period} is outside nper={self.nper}.")
            for move in moves:
                if not isinstance(move, Move):
                    raise ValueError("moves must contain Move objects.")
                result[period].append(move.to_record())
        return result

    @property
    def packages(self) -> list[list[str]]:
        """Return package declarations needed by the MVR package."""

        names: list[str] = []
        for records in self.perioddata.values():
            for record in records:
                for package in (record[0], record[2]):
                    if package not in names:
                        names.append(package)
        return [[name] for name in names]

    def validate(self) -> None:
        """Raise clear errors for invalid semantic MVR inputs."""

        if not self.packages:
            raise ValueError("MVRBuilder requires at least one move.")
        for _period, records in self.perioddata.items():
            for record in records:
                for index in (record[1], record[3]):
                    if int(index) < 0:
                        raise ValueError("MVR mover indexes must be zero-based and non-negative.")

    def build(self) -> PackageSpec:
        """Validate stored configuration and return its package spec."""

        self.validate()
        spec = mvr_spec(
            self.packages,
            self.perioddata,
            name=self.name,
            print_input=self.print_input,
            print_flows=self.print_flows,
            **dict(self.options),
        )
        return spec.with_metadata(
            builder="MVRBuilder",
            package_names=[record[0] for record in self.packages],
            move_count=sum(len(records) for records in self.perioddata.values()),
        )


__all__ = ["MVRBuilder", "Move", "MoverConnection"]
