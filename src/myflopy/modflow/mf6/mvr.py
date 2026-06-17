"""Semantic helpers for building MODFLOW 6 MVR package specifications."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

from myflopy.advanced import mvr_spec
from myflopy.specs import PackageSpec


@dataclass(frozen=True, slots=True)
class MoverConnection:
    """One package endpoint that can send or receive moved water."""

    package: str
    index: int

    def __post_init__(self) -> None:
        if not self.package:
            raise ValueError("package is required.")
        if int(self.index) < 0:
            raise ValueError("index must be zero-based and non-negative.")
        object.__setattr__(self, "package", str(self.package).lower())
        object.__setattr__(self, "index", int(self.index))


@dataclass(frozen=True, slots=True)
class Move:
    """Move water from one MODFLOW package endpoint to another."""

    source: MoverConnection
    receiver: MoverConnection
    method: str = "FACTOR"
    value: float = 1.0

    def __post_init__(self) -> None:
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
    """Prepare an MVR package from semantic package-to-package moves."""

    nper: int
    moves: Mapping[int, Sequence[Move]] | Sequence[Move]
    name: str = "mvr"
    print_input: bool = False
    print_flows: bool = False
    modelnames: bool = False
    options: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
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
