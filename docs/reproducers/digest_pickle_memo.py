"""`package_digest` is sensitive to string IDENTITY, not only to value.

Ledger 284 found this while giving a lazy areal builder a ``package="rch"`` keyword.

`package_digest` builds a payload of six elements and pickles it ONCE::

    payload = (
        _canonical(package.name),                                    # a bare str
        _callable_ref(package.builder),                              # a dict
        pickle.dumps(_canonical(package.options), protocol=5),       # pre-pickled
        pickle.dumps(_canonical(package.inputs), protocol=5),        # pre-pickled
        _canonical(package.enabled),
        _canonical(tuple(sorted(package.requires))),
    )
    return hashlib.sha256(pickle.dumps(payload, protocol=5)).hexdigest()[:16]

`options` and `inputs` are pickled SEPARATELY first, so their internal memo tables are
self-contained -- which is why this has never bitten. The remaining elements share one
memo table. If the builder reference contains a string EQUAL to the package name, a
freshly declared spec has the same interned object in both places and pickle emits a
4-byte back-reference for the second occurrence; a spec revived from JSON holds two
distinct equal strings and pickle writes the string out twice. Equal values, different
bytes, different digest -- and `runs.drift()` then reports a package nobody touched as
EDITED, permanently, because the reloaded value never goes back.

Measured here: 168 bytes against 172.

**Not fixed, deliberately.** Pre-pickling the remaining four elements would make the
digest memo-independent, and would move EVERY stored digest -- a migration, not a
side effect of a feature commit. Ledger 284 avoids the collision instead, by keeping the
lazy areal partial's ``keywords`` empty.

**The exposure is narrow**, and that is what makes avoidance reasonable: the collision
needs a second string equal to `name` among the four un-pre-pickled elements, i.e. inside
the builder reference or `requires`. Nothing in the library puts one there today.

Run: `python docs/reproducers/digest_pickle_memo.py`
"""

from __future__ import annotations

import functools
import hashlib
import json
import pickle


def _build(package_class, model, **kwargs):  # pragma: no cover - never called
    """Stand-in for a lazy builder bound to its FloPy class."""

    raise NotImplementedError


def _callable_ref(builder) -> dict:
    """The shape `specs._callable_ref` produces for a `functools.partial`."""

    return {
        "partial": f"{builder.func.__module__}:{builder.func.__qualname__}",
        "args": [{"$callable": "flopy.mf6.modflow.mfgwfrch:ModflowGwfrch"}],
        "keywords": dict(builder.keywords),
    }


def _digest(name: str, builder) -> tuple[str, int]:
    payload = (name, _callable_ref(builder))
    blob = pickle.dumps(payload, protocol=5)
    return hashlib.sha256(blob).hexdigest()[:16], len(blob)


def main() -> None:
    declared_name = "rch"
    # A DECLARED spec: `name` and the partial's keyword are the same interned literal.
    declared = functools.partial(_build, package=declared_name)

    # A REVIVED spec: both strings come back from JSON as distinct objects.
    revived_name = json.loads(json.dumps(declared_name))
    revived = functools.partial(_build, package=json.loads(json.dumps(declared_name)))

    with_keyword = [
        _digest(declared_name, declared),
        _digest(revived_name, revived),
    ]
    print("WITH a `package=` keyword in the partial")
    print(f"  declared : {with_keyword[0][0]}  ({with_keyword[0][1]} bytes)")
    print(f"  revived  : {with_keyword[1][0]}  ({with_keyword[1][1]} bytes)")
    print(f"  stable   : {with_keyword[0][0] == with_keyword[1][0]}")
    print(f"  same string object when declared: {declared_name is declared.keywords['package']}")
    print(f"  same string object when revived : {revived_name is revived.keywords['package']}")

    # The fix ledger 284 took: no keywords at all, so nothing can collide with `name`.
    bare_declared = functools.partial(_build)
    bare_revived = functools.partial(_build)
    without = [_digest(declared_name, bare_declared), _digest(revived_name, bare_revived)]
    print()
    print("WITHOUT any keyword (what the library ships)")
    print(f"  declared : {without[0][0]}  ({without[0][1]} bytes)")
    print(f"  revived  : {without[1][0]}  ({without[1][1]} bytes)")
    print(f"  stable   : {without[0][0] == without[1][0]}")

    assert not with_keyword[0][0] == with_keyword[1][0], "the defect did not reproduce"
    assert without[0][0] == without[1][0], "the shipped shape should be stable"
    print("\nBoth claims reproduced.")


if __name__ == "__main__":
    main()
