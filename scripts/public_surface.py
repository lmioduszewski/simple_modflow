"""The public surface of myflopy: one discovery for every guard that walks it.

Two guards measure the SAME set of callables, so they read it from here rather
than each keeping a walk of its own:

- ``tests/test_public_parameters_documented.py`` -- every parameter of every
  callable found here has a numpydoc entry (ledger 338);
- ``scripts/derive_api_snapshot.py`` -- every callable found here has its
  signature pinned in ``tests/api_snapshot.json`` (ledger 443), so a default
  or a keyword changing on a public method shows in the snapshot diff.

The surface: every name in ``myflopy.__all__`` and ``myflopy.__engine__``, the
public methods on the helper objects (``mf.chd.gpkg``, ``mf.rch.array``, ...), and
the public methods myflopy defines on those exported classes. Names starting with
``_`` are private. A class whose ``__init__`` belongs to another library (the FloPy
model subclasses, the ``dict`` libraries) is not measured as a constructor, but its
own methods are. Each function is found ONCE, under the class that defines it,
however many exported classes inherit it.

The REACHED tier (ledger 436) adds the classes a user gets to from that surface
without importing them -- the read grammar (``model.<pkg>.<inputs|results>.<noun>``,
``model.outputs.<pkg>``, ``group.packages.<pkg>``, ``diff.packages.<pkg>``), the
bound ``plot`` namespaces, ``model.targets`` and what their methods return --
found by `reached_classes` with no list. Their public methods, plus ``__call__`` and
``__getitem__``, are measured by both guards; their constructors are not, since
myflopy builds every one of them.
"""

from __future__ import annotations

import ast
import dataclasses
import enum
import functools
import importlib
import inspect
import pkgutil
import re
import sys
import textwrap
import types
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src"


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
    helper_types = []

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
                helper_types.append(cls)
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

    reached = reached_classes([*classes, *helper_types])
    for cls in sorted(reached, key=lambda c: (c.__module__, c.__qualname__)):
        for klass in cls.__mro__:
            if not _ours(klass):
                continue
            for mname, raw in klass.__dict__.items():
                if mname.startswith("_") and mname not in _CALLED:
                    continue
                if isinstance(raw, (classmethod, staticmethod)):
                    func = getattr(klass, mname)
                elif inspect.isfunction(raw):
                    func = raw
                else:
                    continue
                add_routine(f"{klass.__qualname__}.{mname}", "reached", func)
    return subjects, skipped


# ------------------------------------------------------ the reached classes
#: Dunders that ARE the public surface of a view: ``model.lak(package=...)`` is
#: ``LakPackageExplorer.__call__``, ``model.regions["x"]`` a ``__getitem__``.
_CALLED = ("__call__", "__getitem__")
#: Members whose bodies NAVIGATE: reading ``obj.<name>`` (a property, or
#: ``__getattr__`` for the registry's fields), indexing or calling the object.
_NAVIGATING = ("__getattr__", *_CALLED)
_TOKEN = re.compile(r"[A-Za-z_][A-Za-z_0-9]*")


def _ours(obj) -> bool:
    """Written in myflopy, not vendored into it.

    ``myflopy._vendor.figs`` is a snapshot of the author's ``figs`` project that its
    README forbids editing by hand (``scripts/sync_vendored_figs.py`` rewrites it), so
    a gap there is closed upstream, not here. Reached through ``GridMesh.fig``.
    """

    return _myflopy(obj) and not obj.__module__.startswith("myflopy._vendor")


def _canonical(cls) -> bool:
    """``cls`` is the object its module names, not a ghost.

    ``@dataclass(slots=True)`` builds a NEW class, and the one it replaced stays in
    its base's ``__subclasses__()`` until collected (measured: every ``DataSourceSpec``
    subclass is listed twice), so a subclass walk would measure each one again.
    """

    obj = sys.modules.get(cls.__module__)
    for part in cls.__qualname__.split("."):
        obj = getattr(obj, part, None)
    return obj is cls


def _class_index():
    """Every myflopy class by its ``__name__``, after importing every myflopy module.

    A view module imports its classes under ``TYPE_CHECKING`` or inside a function,
    so the name a body or annotation uses is often not in the module's globals.
    Importing every module costs ~1.5 s, once.
    """

    import myflopy

    for info in pkgutil.walk_packages(myflopy.__path__, "myflopy."):
        importlib.import_module(info.name)
    index: dict[str, set] = {}
    for module in list(sys.modules.values()):
        if not (getattr(module, "__name__", "") or "").startswith("myflopy"):
            continue
        for value in vars(module).values():
            if inspect.isclass(value) and _myflopy(value) and _canonical(value):
                index.setdefault(value.__name__, set()).add(value)
    return index


def _code_names(code) -> set[str]:
    names = set(code.co_names)
    for const in code.co_consts:
        if isinstance(const, types.CodeType):
            names |= _code_names(const)
    return names


def _annotation_names(func) -> set[str]:
    ann = getattr(func, "__annotations__", {}).get("return")
    if ann is None:
        return set()
    text = ann if isinstance(ann, str) else getattr(ann, "__qualname__", repr(ann))
    return set(_TOKEN.findall(text))


def _documented_return_names(func) -> set[str]:
    """The type named on the first line of a numpydoc ``Returns`` section."""

    found = re.search(r"^Returns\n-+\n(.+)$", inspect.getdoc(func) or "", re.M)
    if not found:
        return set()
    line = found.group(1)
    return set(_TOKEN.findall(line.split(" : ", 1)[1] if " : " in line else line))


def _constructed_return_names(func) -> set[str]:
    """The class a ``return X(...)`` statement builds, read from the source."""

    try:
        tree = ast.parse(textwrap.dedent(inspect.getsource(func)))
    except (OSError, TypeError, SyntaxError):
        return set()
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Return) and isinstance(node.value, ast.Call):
            callee = node.value.func
            if isinstance(callee, ast.Name):
                names.add(callee.id)
            elif isinstance(callee, ast.Attribute):
                names.add(callee.attr)
    return names


def _functions_of(raw):
    if isinstance(raw, property):
        return [f for f in (raw.fget,) if f is not None]
    if isinstance(raw, functools.cached_property):
        return [raw.func]
    if isinstance(raw, (classmethod, staticmethod)):
        return [raw.__func__]
    return [raw] if inspect.isfunction(raw) else []


def reached_classes(exported):
    """Every myflopy class a user reaches from ``exported``, and is not one of them.

    This is how the read grammar is found without listing it: ``model.lak`` is a
    property, ``model.lak.results.stage`` two more, ``model.npf.inputs.k`` a
    ``__getattr__``, ``model.lak(package="lak_2")`` a ``__call__``, ``group.hds`` an
    attribute ``ModelGroup.__init__`` sets, ``model.diff(other)`` a method whose
    ``Returns`` names ``ModelDiff``. From each class reached, an edge goes to every
    myflopy class named by

    - the BODY or return annotation of a public property, ``__getattr__``,
      ``__getitem__`` or ``__call__``, following the module functions and the
      class's own private methods it calls (``model.hds`` is ``get_hds(self)``, which
      builds a ``HeadsPlus``);
    - the value ``__init__`` assigns to ANY ``self.<name>``, private too:
      ``group.packages.rch.inputs`` is the ``self._rch`` ``ModelGroup.__init__``
      builds and ``GroupPackages.rch`` hands on;
    - the return annotation or ``Returns`` type of a public method, or the class a
      ``return X(...)`` in it builds (``stack.plot.map()`` returns a
      ``LayerThicknessMap`` while its generated ``Returns`` says ``Choro``) -- NOT the
      rest of its body, which names every class it merely USES: following whole
      bodies instead reached 27 classes these rules do not (196 against 169), all
      internal (``Boundaries`` and its legacy subclasses, hover blocks,
      grid-building records, ``_Bundler``).

    And from every class reached to its subclasses, because an accessor typed as a
    base (``CellBudgetResultsExplorer``) answers a subclass
    (``LakBudgetResultsExplorer``).

    Blind spot: an attribute assigned OUTSIDE ``__init__`` is not an edge
    (``vor.tessellation_report``, set by the repair, is a ``TessellationReport``
    found by nothing here); its one public method takes no parameter, so nothing
    measurable is missed today.

    Measured 2026-10-11 against walking a run canonical model's attributes, and a
    ``ModelGroup`` and a ``ModelDiff`` of it, by ``getattr``: each of the 107 classes
    that walk met is found here or exported, and the PRT views (the canonical model
    has no PRT) are found besides. Return annotations ALONE missed 56 of the 107
    (``model.hds``, ``model.budget`` and every ``group.*`` and ``diff.*`` namespace
    are unannotated), and ``SpatialView``'s subclasses alone miss every namespace
    that draws nothing. The walk is static: a model built at collection would cost
    an MF6 run in every xdist worker and still miss what its packages do not
    exercise. Vendored classes are left out (`_ours`).
    """

    index = _class_index()
    exported = set(exported)

    def classes_named(names, globals_):
        out = set()
        for name in names:
            value = globals_.get(name)
            candidates = (
                {value} if inspect.isclass(value) else index.get(name, set())
                if value is None else set()
            )
            out |= {c for c in candidates
                    if _ours(c) and _canonical(c) and not issubclass(c, BaseException)}
        return out

    def body(func, owner, done):
        func = inspect.unwrap(func)
        code = getattr(func, "__code__", None)
        if code is None or id(func) in done:
            return set()
        done.add(id(func))
        names = _code_names(code) | _annotation_names(func)
        out = classes_named(names, func.__globals__)
        for name in names:
            value = func.__globals__.get(name)
            if inspect.isfunction(value) and _myflopy(value):
                out |= body(value, None, done)
            elif owner is not None and name.startswith("_") and not name.startswith("__"):
                raw = inspect.getattr_static(owner, name, None)
                for inner in _functions_of(raw):
                    out |= body(inner, owner, done)
        return out

    def init_attributes(cls, init):
        try:
            tree = ast.parse(textwrap.dedent(inspect.getsource(init)))
        except (OSError, TypeError, SyntaxError):
            return set()
        names: set[str] = set()
        for node in ast.walk(tree):
            targets = (node.targets if isinstance(node, ast.Assign)
                       else [node.target] if isinstance(node, ast.AnnAssign) else [])
            if node.__class__ in (ast.Assign, ast.AnnAssign) and node.value is not None and any(
                isinstance(t, ast.Attribute) and isinstance(t.value, ast.Name)
                and t.value.id == "self"
                for t in targets
            ):
                for sub in ast.walk(node.value):
                    if isinstance(sub, ast.Name):
                        names.add(sub.id)
                    elif isinstance(sub, ast.Attribute):
                        names.add(sub.attr)
        out = classes_named(names, init.__globals__)
        done: set[int] = set()
        for name in names:
            value = init.__globals__.get(name)
            if inspect.isfunction(value) and _myflopy(value):
                out |= body(value, None, done)
            elif name.startswith("_") and not name.startswith("__"):
                for inner in _functions_of(inspect.getattr_static(cls, name, None)):
                    out |= body(inner, cls, done)
        return out

    def edges(cls):
        out = set()
        for klass in cls.__mro__:
            if not _myflopy(klass):
                continue
            for name, raw in klass.__dict__.items():
                if name == "__init__" and inspect.isfunction(raw):
                    out |= init_attributes(cls, raw)
                    continue
                navigating = (
                    isinstance(raw, (property, functools.cached_property))
                    and not name.startswith("_")
                ) or (name in _NAVIGATING and inspect.isfunction(raw))
                for func in _functions_of(raw):
                    if navigating:
                        out |= body(func, cls, set())
                    elif not name.startswith("_"):
                        out |= classes_named(
                            _annotation_names(func) | _documented_return_names(func)
                            | _constructed_return_names(func),
                            func.__globals__,
                        )
        return out

    def subclasses(cls):
        out, stack = set(), [cls]
        while stack:
            for sub in stack.pop().__subclasses__():
                if _ours(sub) and _canonical(sub) and sub not in out:
                    out.add(sub)
                    stack.append(sub)
        return out

    seen: set = set()
    stack = list(exported)
    while stack:
        cls = stack.pop()
        if cls in seen:
            continue
        seen.add(cls)
        found = edges(cls)
        for target in list(found):
            found |= subclasses(target)
        stack.extend(found - seen)
    return seen - exported
