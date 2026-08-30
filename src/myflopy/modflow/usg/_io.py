"""Low-level readers for MODFLOW-USG text input.

MODFLOW's array input is a *control record* followed (sometimes) by the values
themselves, and the same four words appear in every package::

    CONSTANT   <value>
    INTERNAL   <multiplier> <fmtin> <iprn>   [free-text label]
    EXTERNAL   <unit> <multiplier> <fmtin> <iprn>   [free-text label]
    OPEN/CLOSE <fname> <multiplier> <fmtin> <iprn>  [free-text label]

``EXTERNAL`` names a *unit number*, which only the name file can resolve to a
path -- which is why :class:`NameFile` has to be read before anything else.

A **negative** unit number means the record is unformatted (binary). MODFLOW's
``U2DREL`` reads such a record as a header (``KSTP KPER PERTIM TOTIM TEXT NCOL
NROW ILAY``) followed by the values, which is exactly the layout of a head-save
file -- so a model can, and the Ten Trails model does, seed its starting heads
from a previous run's ``.hds``. FloPy 3.10 cannot read that: it routes the
record to :func:`flopy.utils.util_array.Util2d.load_txt`, which tests ``"," in
line`` against ``bytes`` and raises ``TypeError``. Its BAS loader therefore
fails outright on such a model, and its LPF loader then fails too because it
reaches for the BAS that never loaded. :class:`ArrayCursor` handles the binary
case by delegating to :class:`flopy.utils.HeadUFile`, which reads the same file
correctly.

Nothing here is USG-specific except the assumption that arrays are addressed by
node count rather than (nrow, ncol); a structured MODFLOW file would read the
same way given the right count.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from myflopy._logging import get_logger

logger = get_logger(__name__)

__all__ = ["ArrayCursor", "NameFile", "free_floats", "free_ints"]

#: Matches one free-format number, including Fortran ``D`` exponents and the
#: bare-exponent form (``-9.99000e+002``) that Groundwater Vistas writes.
_NUMBER = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?")

#: Array control words, in the order the parser tests them.
_CONTROL_WORDS = ("CONSTANT", "INTERNAL", "EXTERNAL", "OPEN/CLOSE")


def _to_float(token: str) -> float:
    """Parse one MODFLOW numeric token, accepting Fortran ``D`` exponents."""

    return float(token.replace("D", "E").replace("d", "e"))


def free_floats(text: str) -> list[float]:
    """Return every free-format number on ``text`` as a float."""

    return [_to_float(t) for t in _NUMBER.findall(text)]


def free_ints(text: str) -> list[int]:
    """Return every free-format number on ``text`` as an int (truncating)."""

    return [int(_to_float(t)) for t in _NUMBER.findall(text)]


def _is_comment(line: str) -> bool:
    """True for a MODFLOW comment line (``#`` in column 1) or a blank line."""

    stripped = line.strip()
    return not stripped or stripped.startswith("#")


@dataclass(frozen=True, slots=True)
class NameFileEntry:
    """One row of a MODFLOW name file."""

    ftype: str
    unit: int
    path: Path
    binary: bool = False
    replace: bool = False


@dataclass(slots=True)
class NameFile:
    """A parsed MODFLOW-USG name file: package rows plus the unit-number map.

    ``EXTERNAL`` array records address files by unit number, so every package
    reader needs this map. Both ``DATA`` and ``DATA(BINARY)`` rows are recorded;
    ``binary`` distinguishes them.

    Attributes
    ----------
    path
        The name file itself.
    workspace
        Directory the name file lives in; every relative path resolves against it.
    entries
        Every row, in file order.
    """

    path: Path
    workspace: Path
    entries: tuple[NameFileEntry, ...] = ()
    #: unit number -> entry, for resolving ``EXTERNAL``.
    units: dict[int, NameFileEntry] = field(default_factory=dict)

    @classmethod
    def read(cls, path: str | Path) -> NameFile:
        """Parse a name file.

        Rows are ``FTYPE UNIT FNAME [REPLACE]``. A path is resolved
        case-insensitively against the workspace, because these files are
        routinely written on Windows and read on Linux -- the Ten Trails name
        file says ``SYQvr.dat`` where the directory holds ``SyQvr.dat``.
        """

        path = Path(path)
        workspace = path.parent
        # One case-folded index of the directory, so each lookup is a dict hit
        # rather than a directory scan.
        on_disk = {p.name.lower(): p for p in workspace.iterdir() if p.is_file()}

        entries: list[NameFileEntry] = []
        for raw in path.read_text(errors="replace").splitlines():
            if _is_comment(raw):
                continue
            parts = raw.split()
            if len(parts) < 3:
                continue
            ftype = parts[0].upper()
            try:
                unit = int(parts[1])
            except ValueError:
                logger.debug("name-file row has no unit number, skipped: %r", raw.strip())
                continue
            fname = parts[2]
            resolved = on_disk.get(Path(fname).name.lower(), workspace / fname)
            entries.append(
                NameFileEntry(
                    ftype=ftype,
                    unit=unit,
                    path=resolved,
                    binary=ftype == "DATA(BINARY)",
                    replace=any(p.upper() == "REPLACE" for p in parts[3:]),
                )
            )

        units = {e.unit: e for e in entries}
        return cls(path=path, workspace=workspace, entries=tuple(entries), units=units)

    def package(self, ftype: str) -> NameFileEntry | None:
        """Return the row for ``ftype`` (e.g. ``"DISU"``), or ``None``."""

        wanted = ftype.upper()
        for entry in self.entries:
            if entry.ftype == wanted:
                return entry
        return None

    @property
    def package_types(self) -> tuple[str, ...]:
        """Every package FTYPE present, excluding the DATA/LIST plumbing rows."""

        skip = {"DATA", "DATA(BINARY)", "LIST", "GLOBAL"}
        return tuple(e.ftype for e in self.entries if e.ftype not in skip)


class ArrayCursor:
    """A position in one MODFLOW text file, able to read array control records.

    The cursor owns the file's lines and an index into them. Every ``read_*``
    advances past what it consumed, so a package reader is a straight-line
    sequence of reads in input order -- which is how the MODFLOW input format is
    defined and the only way to parse it, since records carry no self-describing
    length.

    Parameters
    ----------
    path
        The file to read.
    name_file
        Supplies the unit-number map for ``EXTERNAL`` records. Optional only for
        files that contain no ``EXTERNAL`` record.
    """

    def __init__(self, path: str | Path, name_file: NameFile | None = None) -> None:
        self.path = Path(path)
        self.name_file = name_file
        self.lines: list[str] = self.path.read_text(errors="replace").splitlines()
        self.pos: int = 0
        #: Cache of binary EXTERNAL payloads, keyed by unit; each read of a
        #: negative unit takes the next layer from the same file.
        self._binary_cache: dict[int, list[np.ndarray]] = {}
        self._binary_taken: dict[int, int] = {}

    # -- cursor movement ---------------------------------------------------

    def peek(self) -> str:
        """Return the next significant line without consuming it."""

        while self.pos < len(self.lines) and _is_comment(self.lines[self.pos]):
            self.pos += 1
        return self.lines[self.pos] if self.pos < len(self.lines) else ""

    def next_line(self) -> str:
        """Consume and return the next significant (non-comment, non-blank) line."""

        line = self.peek()
        self.pos += 1
        return line

    def at_end(self) -> bool:
        """True once no significant line remains."""

        return not self.peek()

    # -- scalar records ----------------------------------------------------

    def read_ints(self, count: int | None = None) -> list[int]:
        """Read one line as integers; if ``count`` is given, keep reading lines."""

        if count is None:
            return free_ints(self.next_line())
        out: list[int] = []
        while len(out) < count:
            out.extend(free_ints(self.next_line()))
        return out[:count]

    def read_floats(self, count: int | None = None) -> list[float]:
        """Read one line as floats; if ``count`` is given, keep reading lines."""

        if count is None:
            return free_floats(self.next_line())
        out: list[float] = []
        while len(out) < count:
            out.extend(free_floats(self.next_line()))
        return out[:count]

    # -- array control records ---------------------------------------------

    def read_array(self, count: int, dtype: type = float) -> np.ndarray:
        """Read one array control record and return ``count`` values.

        Handles ``CONSTANT`` / ``INTERNAL`` / ``EXTERNAL`` / ``OPEN/CLOSE``,
        applies the record's multiplier, and follows a negative ``EXTERNAL``
        unit into its binary file.
        """

        line = self.next_line()
        word = _control_word(line)
        if word is None:
            raise ValueError(
                f"{self.path.name}: expected an array control record "
                f"(one of {', '.join(_CONTROL_WORDS)}) at line {self.pos}, got: {line.strip()!r}"
            )

        if word == "CONSTANT":
            value = _token_float(line, 1, default=0.0)
            return np.full(count, value, dtype=dtype)

        if word == "INTERNAL":
            multiplier = _multiplier(line, token=1)
            values = self.read_floats(count)
            return (np.asarray(values, dtype=float) * multiplier).astype(dtype)

        if word == "OPEN/CLOSE":
            fname = line.split()[1].strip("'\"")
            multiplier = _multiplier(line, token=2)
            target = self._resolve(fname)
            values = _read_text_values(target, count)
            return (values * multiplier).astype(dtype)

        # EXTERNAL
        unit = int(_token_float(line, 1, default=0.0))
        multiplier = _multiplier(line, token=2)
        entry = self._unit(unit)
        if unit < 0 or (entry is not None and entry.binary):
            values = self._read_binary_layer(unit, count)
        else:
            if entry is None:
                raise ValueError(
                    f"{self.path.name}: EXTERNAL unit {unit} is not declared in the name file"
                )
            values = _read_text_values(entry.path, count, unit=unit)
        return (values * multiplier).astype(dtype)

    # -- helpers -----------------------------------------------------------

    def _unit(self, unit: int) -> NameFileEntry | None:
        """Look up a unit number, tolerating the sign that marks binary."""

        if self.name_file is None:
            return None
        return self.name_file.units.get(abs(unit)) or self.name_file.units.get(unit)

    def _resolve(self, fname: str) -> Path:
        """Resolve an ``OPEN/CLOSE`` name against the workspace, case-insensitively."""

        base = self.name_file.workspace if self.name_file else self.path.parent
        candidate = base / fname
        if candidate.is_file():
            return candidate
        wanted = Path(fname).name.lower()
        for found in base.iterdir():
            if found.is_file() and found.name.lower() == wanted:
                return found
        return candidate

    def _read_binary_layer(self, unit: int, count: int) -> np.ndarray:
        """Take the next layer's worth of values from a binary EXTERNAL unit.

        Successive ``EXTERNAL -<unit>`` records walk successive layers of the
        same file, so the payload is read once and handed out in order.
        """

        key = abs(unit)
        if key not in self._binary_cache:
            entry = self._unit(unit)
            if entry is None:
                raise ValueError(
                    f"{self.path.name}: binary EXTERNAL unit {unit} is not declared "
                    "in the name file"
                )
            self._binary_cache[key] = _read_binary_layers(entry.path)
            self._binary_taken[key] = 0
            logger.debug(
                "read %d binary layer(s) from unit %d (%s)",
                len(self._binary_cache[key]),
                key,
                entry.path.name,
            )

        taken = self._binary_taken[key]
        layers = self._binary_cache[key]
        if taken >= len(layers):
            raise ValueError(
                f"{self.path.name}: binary unit {unit} holds {len(layers)} layer(s) "
                f"but a {taken + 1}th was requested"
            )
        self._binary_taken[key] = taken + 1
        values = np.asarray(layers[taken], dtype=float).ravel()
        if values.size != count:
            raise ValueError(
                f"binary unit {unit} layer {taken + 1} has {values.size} values, expected {count}"
            )
        return values


def _control_word(line: str) -> str | None:
    """Return the array control word on ``line``, or ``None`` if there is none."""

    first = line.strip().split()[0].upper() if line.strip() else ""
    for word in _CONTROL_WORDS:
        if first == word:
            return word
    return None


def _token_float(line: str, index: int, *, default: float) -> float:
    """Return whitespace-token ``index`` of ``line`` as a float.

    Control records must be read by TOKEN position, never by position among the
    numbers the line happens to contain: ``OPEN/CLOSE`` carries a *filename* in
    token 1, and a name like ``tops2.dat`` would otherwise contribute a phantom
    number and shift every field after it. That is not hypothetical -- indexing
    by number position silently read ``IPRN = -1`` as the multiplier and negated
    the entire ET surface.
    """

    tokens = line.split()
    if len(tokens) <= index:
        return default
    try:
        return float(tokens[index].replace("D", "E").replace("d", "e"))
    except ValueError:
        return default


def _multiplier(line: str, *, token: int) -> float:
    """Return the ``CNSTNT`` multiplier of a control record, defaulting to 1.0.

    ``token`` is the whitespace-token index of ``CNSTNT``: 1 for ``INTERNAL``,
    2 for ``EXTERNAL`` (which leads with its unit) and ``OPEN/CLOSE`` (which
    leads with its filename). A multiplier of 0 means 1.0, per MODFLOW.
    """

    value = _token_float(line, token, default=1.0)
    return 1.0 if value == 0.0 else value


def _read_text_values(path: Path, count: int, *, unit: int | None = None) -> np.ndarray:
    """Read the first ``count`` free-format numbers from a text file."""

    values: list[float] = []
    with open(path, errors="replace") as handle:
        for line in handle:
            values.extend(_to_float(t) for t in _NUMBER.findall(line))
            if len(values) >= count:
                break
    if len(values) < count:
        where = f"unit {unit} ({path.name})" if unit is not None else path.name
        raise ValueError(f"{where} holds {len(values)} values, expected {count}")
    return np.asarray(values[:count], dtype=float)


def _read_binary_layers(path: Path) -> list[np.ndarray]:
    """Return the last time step's per-layer arrays from a MODFLOW binary file.

    A negative ``EXTERNAL`` unit is an unformatted record with a head-file
    header, so the file is read as one. MODFLOW reads such a unit *sequentially*
    -- each array record consumes the next record in the file -- so the FIRST
    time step is the one a set of per-layer reads lands on, not the last. A
    restart-style ``initial_*.hds`` holds a single time step either way.
    """

    from flopy.utils import HeadUFile

    heads = HeadUFile(str(path))
    data = heads.get_data(kstpkper=heads.get_kstpkper()[0])
    return [np.asarray(layer, dtype=float).ravel() for layer in data]
