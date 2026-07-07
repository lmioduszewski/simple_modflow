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

    Attributes
    ----------
    package
        The advanced package name (e.g. ``"sfr"``, ``"lak"``, ``"uzf"``). The
        package must be declared with ``mover=True``.
    index
        Zero-based feature id within that package (reach / lake / cell number).

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

    Attributes
    ----------
    source, receiver
        The provider and receiver :class:`MoverConnection` endpoints.
    method
        How ``value`` is interpreted: ``"FACTOR"`` (move this fraction of available
        water, default), ``"EXCESS"`` (water above ``value``), ``"THRESHOLD"``, or
        ``"UPTO"`` (at most ``value``).
    value
        The factor/rate for ``method`` (default ``1.0`` = all available water).

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
    """

    nper: int
    moves: Mapping[int, Sequence[Move]] | Sequence[Move]
    name: str = "mvr"
    print_input: bool = False
    print_flows: bool = False
    modelnames: bool = False
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
            modelnames=self.modelnames,
            **dict(self.options),
        )
        return spec.with_metadata(
            builder="MVRBuilder",
            package_names=[record[0] for record in self.packages],
            move_count=sum(len(records) for records in self.perioddata.values()),
        )


__all__ = ["MVRBuilder", "Move", "MoverConnection"]
