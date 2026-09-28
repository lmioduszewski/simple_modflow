"""Every public parameter is documented: a ``name : type`` entry with a description.

PyCharm and Pylance are static readers. What they show for a call is the ``def``
line and the docstring written in the source, so a parameter with no entry is one
the user can see and cannot look up. Until this file, only the plot verbs were
checked (`test_plot_vocabulary`, `test_noun_signatures`, `test_docstrings_are_static`).
Everything else was left to discipline. On 2026-09-27 the rest of the public
surface was missing 2,142 entries or types:

- ``mf.disv`` described only its hand-written form, and 10 of its 19 parameters
  were undocumented (the ledger-270 stack form).
- All seven ``.gpkg`` helpers gained ``gaps=`` and ``snap_within=`` without an
  entry for either.
- ``LakeTableBuilder`` named all 11 of its parameters but typed none of them,
  and its own example raised.

That last case is why the rule is an ENTRY and not a mention. A name alone
passes a docstring that says nothing about the parameter.

The rule, per callable:

1. Every named parameter has an entry in ``Parameters`` or ``Other Parameters``
   whose head reads ``name : type`` and whose body is not empty. A grouped head
   ``a, b : float`` documents both names. For a DATACLASS, an ``Attributes`` entry
   counts too, because its fields are its constructor.
2. A ``*args`` or ``**kwargs`` tail has an entry (``**options``) with a body
   saying where the extra arguments go. No type is required for it.
3. A ``Parameters`` entry names something the callable accepts. The exception is
   a callable with a ``**kwargs`` tail, which may document the keywords it forwards.
   A line at column 0 inside the block that is not a valid entry head is flagged,
   because ``help()`` renders it as a parameter.

Surface: every name in ``myflopy.__all__`` and ``myflopy.__engine__``, the public
methods on the helper objects (``mf.chd.gpkg``, ``mf.rch.array``, ...), and the
public methods myflopy defines on those exported classes. Names starting with ``_``
are private. A class whose ``__init__`` belongs to another library (the FloPy model
subclasses, the ``dict`` libraries) is not measured as a constructor, but its own
methods are.

NOT covered yet (ledger 338): the read-grammar view classes reached as
``model.packages.<pkg>.<inputs|results>.<noun>``. They are not exported, so they
cannot be discovered from ``__all__``. Their plot verbs are already checked by
`test_noun_signatures`.

Run as a script to list what is still missing, optionally only for some files
(paths relative to ``src/myflopy``, matched exactly):

    python tests/test_public_parameters_documented.py specs.py modflow/mf6/grid/voronoi.py
"""

from __future__ import annotations

import dataclasses
import enum
import inspect
import sys
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / "src"

#: Section headings a NumPy docstring may carry (numpydoc's own list).
SECTION_HEADS = frozenset({
    "Parameters", "Other Parameters", "Returns", "Yields", "Receives", "Raises",
    "Warns", "Warnings", "See Also", "Notes", "References", "Examples",
    "Attributes", "Methods",
})
PARAMETER_SECTIONS = ("Parameters", "Other Parameters")
NAMED = (
    inspect.Parameter.POSITIONAL_ONLY,
    inspect.Parameter.POSITIONAL_OR_KEYWORD,
    inspect.Parameter.KEYWORD_ONLY,
)

#: Floors, not exact counts: a discovery that silently finds nothing must fail,
#: but a new public method must not. Measured 2026-09-27: 171 / 24 / 28 / 492, each
#: function counted once however many exported classes inherit it.
MIN_CALLABLES = {"all": 165, "engine": 22, "helper": 26, "method": 470}


# ------------------------------------------------------------------ parsing
@dataclasses.dataclass(frozen=True)
class Entry:
    section: str
    names: tuple[str, ...]      # stars stripped
    type: str
    body: str
    head: str


def _sections(doc: str) -> dict[str, list[str]]:
    """Section title -> its lines, for a cleaned NumPy docstring."""

    lines = doc.splitlines()
    out: dict[str, list[str]] = {}
    current = None
    i = 0
    while i < len(lines):
        line = lines[i]
        nxt = lines[i + 1].strip() if i + 1 < len(lines) else ""
        if (
            line.strip() in SECTION_HEADS
            and not line[:1].isspace()
            and nxt
            and set(nxt) == {"-"}
            and len(nxt) >= len(line.strip())
        ):
            current = line.strip()
            out.setdefault(current, [])
            i += 2
            continue
        if current is not None:
            out[current].append(line)
        i += 1
    return out


def entries(doc: str | None) -> list[Entry]:
    """Every entry of the Parameters, Other Parameters and Attributes sections."""

    if not doc:
        return []
    found: list[Entry] = []
    for section, lines in _sections(inspect.cleandoc(doc)).items():
        if section not in (*PARAMETER_SECTIONS, "Attributes"):
            continue
        head = None
        body: list[str] = []

        def close():
            if head is None:
                return
            name_part, _, type_part = head.partition(" :")
            if ":" in name_part:          # `name: type`, no space before the colon
                name_part, _, type_part = head.partition(":")
            names = tuple(n.strip().lstrip("*") for n in name_part.split(",") if n.strip())
            found.append(Entry(section, names, type_part.strip(),
                               "\n".join(body).strip(), head))

        for line in lines:
            if line.strip() and not line[:1].isspace():
                close()
                head, body = line.rstrip(), []
            elif head is not None:
                body.append(line)
        close()
    return found


# ----------------------------------------------------------------- discovery
@dataclasses.dataclass
class Subject:
    qualname: str
    tier: str
    where: str
    signature: inspect.Signature
    docs: list[str]
    is_dataclass: bool = False


def _where(obj) -> str:
    try:
        path = Path(inspect.getsourcefile(obj)).resolve()
        line = inspect.getsourcelines(obj)[1]
    except (TypeError, OSError):
        return "?"
    try:
        rel = path.relative_to(SRC / "myflopy")
    except ValueError:
        rel = path
    return f"{rel}:{line}"


def _own_class_doc(cls) -> str:
    """The class's OWN docstring; a dataclass's generated ``Name(field: ...)`` is not one."""

    doc = cls.__dict__.get("__doc__") or ""
    if dataclasses.is_dataclass(cls) and doc.startswith(f"{cls.__name__}("):
        return ""
    return doc


def _myflopy(obj) -> bool:
    return (getattr(obj, "__module__", "") or "").startswith("myflopy")


def _init_owner(cls):
    for klass in cls.__mro__:
        if klass is object:
            return None
        if "__init__" in klass.__dict__:
            return klass
    return None


def _class_subject(name, tier, cls):
    """A class measured as a constructor, or the reason it is not."""

    if issubclass(cls, enum.Enum):
        return "an enum"
    if issubclass(cls, BaseException):
        return "an exception"
    owner = _init_owner(cls)
    if owner is not None and not _myflopy(owner) and not dataclasses.is_dataclass(cls):
        return f"__init__ belongs to {owner.__module__}.{owner.__name__}"
    try:
        sig = inspect.signature(cls)
    except (TypeError, ValueError) as exc:
        return f"no signature ({exc})"
    docs = [_own_class_doc(cls)]
    if owner is not None and _myflopy(owner):
        docs.append(owner.__dict__["__init__"].__doc__ or "")
    return Subject(name, tier, _where(cls), sig, docs, dataclasses.is_dataclass(cls))


def _is_helper(obj) -> bool:
    cls = type(obj)
    return (
        not inspect.isclass(obj)
        and not inspect.isroutine(obj)
        and not inspect.ismodule(obj)
        and "__call__" in cls.__dict__
        and _myflopy(cls)
    )


def discover():
    """(subjects, skipped) over the whole public surface, each function once."""

    import myflopy as mf

    subjects: list[Subject] = []
    skipped: dict[str, str] = {}
    seen: set[int] = set()
    classes = []

    def add_routine(qualname, tier, func, doc=None):
        key = id(getattr(func, "__func__", func))
        if key in seen:
            return
        seen.add(key)
        try:
            sig = inspect.signature(func)
        except (TypeError, ValueError) as exc:
            skipped[qualname] = f"no signature ({exc})"
            return
        subjects.append(Subject(qualname, tier, _where(func), sig,
                                [inspect.getdoc(func) if doc is None else doc]))

    for tier, names in (("all", mf.__all__), ("engine", mf.__engine__)):
        for name in names:
            obj = getattr(mf, name)
            if inspect.isclass(obj):
                result = _class_subject(name, tier, obj)
                if isinstance(result, Subject):
                    subjects.append(result)
                else:
                    skipped[name] = result
                if not issubclass(obj, (enum.Enum, BaseException)):
                    classes.append(obj)
            elif inspect.isroutine(obj):
                if _myflopy(obj):
                    add_routine(name, tier, obj)
                else:
                    skipped[name] = f"defined in {obj.__module__}"
            elif _is_helper(obj):
                cls = type(obj)
                call = cls.__dict__["__call__"]
                subjects.append(Subject(f"mf.{name}", tier, _where(call),
                                        inspect.signature(obj),
                                        [call.__doc__ or _own_class_doc(cls)]))
                seen.add(id(call))
                for meth in sorted(m for m in dir(cls) if not m.startswith("_")):
                    bound = getattr(obj, meth)
                    if inspect.isroutine(bound):
                        add_routine(f"mf.{name}.{meth}", "helper", bound)
            else:
                skipped[name] = f"not callable ({type(obj).__name__})"

    for cls in classes:
        for klass in cls.__mro__:
            if not _myflopy(klass):
                continue
            for mname, raw in klass.__dict__.items():
                if mname.startswith("_"):
                    continue
                if isinstance(raw, (classmethod, staticmethod)):
                    func = getattr(klass, mname)
                elif inspect.isfunction(raw):
                    func = raw
                else:
                    continue          # properties and plain attributes
                add_routine(f"{klass.__name__}.{mname}", "method", func)
    return subjects, skipped


# -------------------------------------------------------------------- rule
def complaints(subject: Subject) -> list[str]:
    """What is wrong with one callable's parameter documentation."""

    params = [p for p in subject.signature.parameters.values()
              if p.name not in {"self", "cls"} and not p.name.startswith("_")]
    named = [p.name for p in params if p.kind in NAMED]
    var = [p for p in params if p.kind not in NAMED]
    has_varkw = any(p.kind is inspect.Parameter.VAR_KEYWORD for p in var)

    documented: dict[str, Entry] = {}
    bad_heads: list[str] = []
    param_names: set[str] = set()
    for doc in subject.docs:
        for entry in entries(doc):
            if entry.section == "Attributes" and not subject.is_dataclass:
                continue
            if not entry.names or not all(n.isidentifier() for n in entry.names):
                if entry.section in PARAMETER_SECTIONS:
                    bad_heads.append(entry.head)
                continue
            for n in entry.names:
                if entry.section in PARAMETER_SECTIONS:
                    param_names.add(n)
                    documented[n] = entry
                else:
                    documented.setdefault(n, entry)

    out = []
    missing = [n for n in named if n not in documented]
    if missing:
        out.append(f"no entry for {missing}")
    untyped = [n for n in named if n in documented and not documented[n].type]
    if untyped:
        out.append(f"no `name : type` head for {untyped}")
    bodyless = [n for n in named if n in documented and not documented[n].body]
    if bodyless:
        out.append(f"no description for {bodyless}")
    tails = [("**" if p.kind is inspect.Parameter.VAR_KEYWORD else "*") + p.name
             for p in var if p.name not in documented or not documented[p.name].body]
    if tails:
        out.append(f"no entry saying where {tails} go")
    if not has_varkw:
        accepted = {p.name for p in params} | {"self", "cls"}
        extra = sorted(n for n in param_names - accepted if not n.startswith("_"))
        if extra:
            out.append(f"documents {extra}, which it does not accept")
    if bad_heads:
        out.append(f"column-0 lines in Parameters that are not entries: {bad_heads}")
    return out


SUBJECTS, SKIPPED = discover()
IDS = [s.qualname for s in SUBJECTS]


def test_discovery_found_the_public_surface():
    """A discovery that found nothing would pass every case below."""

    counts = {tier: sum(s.tier == tier for s in SUBJECTS) for tier in MIN_CALLABLES}
    short = {t: (counts[t], floor) for t, floor in MIN_CALLABLES.items() if counts[t] < floor}
    assert not short, f"discovered fewer callables than measured: {short}"


def test_only_known_kinds_are_skipped():
    """A name is skipped for a stated reason, never because discovery could not read it."""

    unexpected = {n: r for n, r in SKIPPED.items() if r.startswith("no signature")}
    assert not unexpected, unexpected


@pytest.mark.parametrize("subject", SUBJECTS, ids=IDS)
def test_every_public_parameter_is_documented(subject):
    problems = complaints(subject)
    assert not problems, f"{subject.where} {subject.qualname}: " + "; ".join(problems)


def _report(files: list[str]) -> int:
    """Print every callable with a gap, optionally only those defined in ``files``.

    ``files`` are paths relative to ``src/myflopy`` and match EXACTLY, so
    ``specs.py`` does not also select ``modflow/mf6/pest/specs.py``.
    """

    total = 0
    for s in sorted(SUBJECTS, key=lambda s: (s.where.split(":")[0], s.qualname)):
        if files and s.where.split(":")[0] not in files:
            continue
        problems = complaints(s)
        if problems:
            total += 1
            print(f"{s.where}  {s.qualname}")
            for p in problems:
                print(f"    {p}")
    print(f"\n{total} callables with gaps")
    return 1 if total else 0


if __name__ == "__main__":
    sys.exit(_report(sys.argv[1:]))
