"""LAK ``IMPLICIT``: MF6's option to solve each lake stage in the flow matrix (ledger 326).

MODFLOW 6's LAK package gained an ``IMPLICIT`` option on its ``develop`` branch (commit
``e2dfc84c``, PR #2860, July 2026; in the 6.8 line, in no release before it): each lake
stage becomes one more unknown of the groundwater-flow matrix, solved WITH the heads,
instead of by MF6's default substitution iteration, which solves the stage and the
heads by turns. FloPy 3.10's ``ModflowGwflak`` has no keyword for it -- its class is
generated from an older DFN -- so ``mf.lak(implicit=True)`` needs a SECOND WRITER for
the LAK file, and this module is it, in one place:

* :class:`ModflowGwflakImplicit` -- ``ModflowGwflak`` whose ``write`` inserts ``IMPLICIT``
  into the file's ``OPTIONS`` block after FloPy writes it (:func:`insert_implicit`).
  Every route that WRITES the package -- ``Run.write``/``execute``, a
  ``BuiltSimulation``, ``sim.write_simulation()`` on a built model, a frozen run rebuilt
  on reopen through the builder -- reaches its own ``write``. The flag IS the class, so
  a route that REBUILDS a LAK package loses it unless it carries the class, and three
  did, silently: FloPy's model splitter (``ParallelModelWorkflow.split_model`` re-classes
  the split packages), a package artifact (it records ``implicit``) and FloPy's loader
  (below). ``model.config`` reads the class too, as the lake's ``implicit`` setting.
* :func:`keep_implicit` -- FloPy LOADS a file holding ``IMPLICIT`` without a word and
  DROPS it (measured: no warning; a rewrite of the loaded simulation wrote no
  ``IMPLICIT``). myflopy's two load paths (``Run.open``, ``mf.load_mf6_run``) call this to
  turn such a loaded package into a :class:`ModflowGwflakImplicit`, so a rewrite keeps it.
  A plain ``flopy.mf6.MFSimulation.load`` + ``write_simulation`` still loses it.
* :func:`implicit_support` / :func:`warn_if_implicit_unsupported` -- an MF6 executable
  that predates the option stops reading the LAK file. The executable is probed ONCE
  (a two-cell model with ``IMPLICIT``, cached per path, modification time and size) when
  a simulation holding an implicit lake is run, and a build without it is WARNED about
  -- the run still goes ahead, and MF6 stops with its own message.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import warnings
from dataclasses import dataclass
from pathlib import Path

import flopy

from myflopy._logging import get_logger
from myflopy._warnings import stacklevel_outside_package

logger = get_logger(__name__)

#: The keyword, as MF6 reads it (any case) and as it is written here.
IMPLICIT = "IMPLICIT"


def insert_implicit(path: str | os.PathLike) -> bool:
    """Put ``IMPLICIT`` in a LAK file's ``OPTIONS`` block, once; ``True`` when it wrote.

    Idempotent: a block that already holds the keyword (any case, as its own line) is
    left alone, so a second write, or a file FloPy wrote from a package that kept it,
    never carries it twice. The keyword goes on its own line directly above ``END
    OPTIONS``; a file with no ``OPTIONS`` block gets one, above its first block.

    Parameters
    ----------
    path : str or os.PathLike
        The written LAK input file.

    Returns
    -------
    bool
        Whether the file was changed.

    Raises
    ------
    OSError
        The file cannot be read or written.

    Examples
    --------
    >>> import pathlib, tempfile
    >>> lak = pathlib.Path(tempfile.mkdtemp()) / "m.lak"
    >>> _ = lak.write_text("BEGIN options\\n  BOUNDNAMES\\nEND options\\n")
    >>> insert_implicit(lak), insert_implicit(lak)
    (True, False)
    >>> print(lak.read_text())
    BEGIN options
      BOUNDNAMES
      IMPLICIT
    END options
    <BLANKLINE>
    """

    path = Path(path)
    lines = path.read_text().splitlines(keepends=True)
    start = end = None
    for index, line in enumerate(lines):
        words = line.split()
        if len(words) >= 2 and words[0].upper() == "BEGIN" and words[1].upper() == "OPTIONS":
            start = index
        elif start is not None and len(words) >= 2 and words[0].upper() == "END" and words[1].upper() == "OPTIONS":
            end = index
            break
    if start is not None and end is not None:
        if any(line.split()[:1] == [IMPLICIT] or (line.split() and line.split()[0].upper() == IMPLICIT)
               for line in lines[start + 1:end]):
            return False
        lines.insert(end, f"  {IMPLICIT}\n")
    else:
        first = next(
            (i for i, line in enumerate(lines) if line.split()[:1] and line.split()[0].upper() == "BEGIN"),
            len(lines),
        )
        lines[first:first] = ["BEGIN options\n", f"  {IMPLICIT}\n", "END options\n", "\n"]
    path.write_text("".join(lines))
    return True


def has_implicit(path: str | os.PathLike) -> bool:
    """Whether a written LAK file's ``OPTIONS`` block holds ``IMPLICIT``.

    Parameters
    ----------
    path : str or os.PathLike
        A LAK input file.

    Returns
    -------
    bool
        ``False`` also for a file that does not exist.

    Examples
    --------
    >>> import pathlib, tempfile
    >>> lak = pathlib.Path(tempfile.mkdtemp()) / "m.lak"
    >>> _ = lak.write_text("BEGIN options\\n  implicit\\nEND options\\n")
    >>> has_implicit(lak), has_implicit(lak.with_suffix(".none"))
    (True, False)
    """

    try:
        lines = Path(path).read_text().splitlines()
    except OSError:
        return False
    inside = False
    for line in lines:
        words = [word.upper() for word in line.split()]
        if words[:2] == ["BEGIN", "OPTIONS"]:
            inside = True
        elif words[:2] == ["END", "OPTIONS"]:
            return False
        elif inside and words[:1] == [IMPLICIT]:
            return True
    return False


class ModflowGwflakImplicit(flopy.mf6.ModflowGwflak):
    """FloPy's ``ModflowGwflak``, writing MF6's LAK ``IMPLICIT`` option it has no keyword for.

    Built by ``mf.lak(implicit=True)`` (and a :class:`~myflopy.LAKBuilder` with
    ``implicit=True``); rarely constructed by hand. It IS a ``ModflowGwflak`` --
    ``package_type`` ``"lak"``, every FloPy keyword, ``isinstance`` holds -- and differs
    only in :meth:`write`, which lets FloPy write the file and then inserts ``IMPLICIT``
    into its ``OPTIONS`` block (:func:`insert_implicit`, once).

    Parameters
    ----------
    model : flopy.mf6.MFModel
        The GWF model the package belongs to.
    loading_package : bool, default False
        FloPy's own flag (``True`` while a file is being loaded).
    **kwargs : Any
        Every ``flopy.mf6.ModflowGwflak`` keyword, unchanged: the constructor IS
        FloPy's.

    Raises
    ------
    flopy.mf6.mfbase.MFDataException
        Anything FloPy's own ``ModflowGwflak`` constructor refuses (a malformed
        ``packagedata`` row, say): the constructor is FloPy's, unchanged.

    See Also
    --------
    insert_implicit : The one edit made to the written file.
    keep_implicit : Turns a LOADED package whose file holds ``IMPLICIT`` into this class.

    Examples
    --------
    >>> import flopy
    >>> sim = flopy.mf6.MFSimulation()
    >>> gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
    >>> lak = ModflowGwflakImplicit(gwf, nlakes=1, packagedata=[[0, 10.0, 0]])
    >>> isinstance(lak, flopy.mf6.ModflowGwflak), lak.package_type
    (True, 'lak')
    """

    #: Read by ``grid_identity.implicit_solver_complaints``, which cannot import this
    #: module without lifting the spec layer above it (tests/import_layers.json).
    writes_lak_implicit = True

    def __init__(self, model, loading_package=False, **kwargs) -> None:
        # FloPy's constructor, unchanged; named here so the three parameters the class
        # documents are the signature an editor shows.
        super().__init__(model, loading_package=loading_package, **kwargs)

    def write(self, ext_file_action=None):
        """Write the LAK file as FloPy does, then put ``IMPLICIT`` in its ``OPTIONS`` block.

        Writing the simulation (``sim.write_simulation()``, ``Run.write``) calls every
        package's ``write``, this one's included, so it is rarely called by hand.

        Parameters
        ----------
        ext_file_action : flopy.mf6.mfbase.ExtFileAction, optional
            Passed to FloPy's ``write``; ``None`` keeps FloPy's default.

        Raises
        ------
        OSError
            The written file cannot be read back or rewritten (:func:`insert_implicit`).
        flopy.mf6.mfbase.MFDataException
            Anything FloPy's own ``write`` refuses.

        Examples
        --------
        >>> import tempfile, flopy
        >>> sim = flopy.mf6.MFSimulation(sim_ws=tempfile.mkdtemp())
        >>> gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
        >>> lak = ModflowGwflakImplicit(gwf, nlakes=1, packagedata=[[0, 10.0, 0]],
        ...                             filename="m.lak")
        >>> lak.write()
        >>> has_implicit(lak.get_file_path())
        True
        """

        if ext_file_action is None:
            super().write()
        else:
            super().write(ext_file_action=ext_file_action)
        insert_implicit(self.get_file_path())


# FloPy's `MFPackage.__init_subclass__` registers EVERY subclass as the class its LOADER
# builds for that package type (`PackageContainer.packages_by_abbr`), so defining the
# class above made every LAK file loaded afterwards -- with or without IMPLICIT in it --
# a `ModflowGwflakImplicit`, and writing it back inserted the keyword (measured: a file
# without it loaded as the subclass). The loader keeps FloPy's own class; a loaded file
# that DOES hold IMPLICIT is re-classed explicitly, by `keep_implicit`.
flopy.mf6.mfbase.PackageContainer.packages_by_abbr[
    ModflowGwflakImplicit.package_abbr
] = flopy.mf6.ModflowGwflak


def keep_implicit(simulation) -> list[str]:
    """Make every loaded LAK package whose FILE holds ``IMPLICIT`` write it again.

    FloPy's loader reads ``IMPLICIT`` without a word and keeps no trace of it, so a
    loaded simulation written back out silently loses the option. Each such package is
    re-classed as :class:`ModflowGwflakImplicit` (FloPy's own object, one method
    different), so writing the simulation writes the keyword again. Run by ``Run.open``
    and ``mf.load_mf6_run``; a bare ``flopy.mf6.MFSimulation.load`` does not.

    Parameters
    ----------
    simulation : flopy.mf6.MFSimulation
        A loaded simulation. Packages not yet loaded (FloPy's ``load_only``) are skipped.

    Returns
    -------
    list of str
        The file names of the packages re-classed.

    Examples
    --------
    Write a lake with ``IMPLICIT``, load it with FloPy (which drops the option), keep it:

    >>> import tempfile, flopy
    >>> workspace = tempfile.mkdtemp()
    >>> sim = flopy.mf6.MFSimulation(sim_ws=workspace)
    >>> _ = flopy.mf6.ModflowTdis(sim)
    >>> ims = flopy.mf6.ModflowIms(sim)
    >>> gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
    >>> sim.register_ims_package(ims, ["m"])
    >>> _ = flopy.mf6.ModflowGwfdis(gwf, nrow=1, ncol=2)
    >>> _ = ModflowGwflakImplicit(
    ...     gwf, nlakes=1, packagedata=[[0, 0.5, 1]], filename="m.lak",
    ...     connectiondata=[[0, 0, (0, 0, 0), "vertical", 1.0, 0.0, 0.0, 0.0, 0.0]])
    >>> sim.write_simulation(silent=True)
    >>> loaded = flopy.mf6.MFSimulation.load(sim_ws=workspace, verbosity_level=0)
    >>> type(loaded.get_model("m").get_package("lak")).__name__
    'ModflowGwflak'
    >>> keep_implicit(loaded)
    ['m.lak']
    >>> type(loaded.get_model("m").get_package("lak")).__name__
    'ModflowGwflakImplicit'
    """

    kept = []
    for name in getattr(simulation, "model_names", None) or ():
        model = simulation.get_model(name)
        for package in getattr(model, "packagelist", ()) or ():
            if str(getattr(package, "package_type", "")) != "lak":
                continue
            if isinstance(package, ModflowGwflakImplicit):
                continue
            if type(package) is not flopy.mf6.ModflowGwflak:
                continue
            if has_implicit(package.get_file_path()):
                package.__class__ = ModflowGwflakImplicit
                kept.append(str(package.filename))
    return kept


@dataclass(frozen=True)
class ImplicitSupport:
    """What probing one MF6 executable found: whether it accepts LAK ``IMPLICIT``.

    Parameters
    ----------
    executable : str
        The executable as asked for (a name looked up on ``PATH``, or a path).
    path : str or None
        Where it was found; ``None`` when it was not.
    version : str or None
        What ``<exe> --version`` printed, e.g. ``"mf6: 6.8.0.dev0+0418efe4"``.
    supported : bool or None
        ``True`` when a model with ``IMPLICIT`` ran to Normal termination, ``False``
        when it did not, ``None`` when nothing could be run.
    detail : str
        MF6's own error lines when it refused, or why nothing was probed.

    Examples
    --------
    >>> ImplicitSupport("mf6", None, None, None, "not found").supported is None
    True
    """

    executable: str
    path: str | None
    version: str | None
    supported: bool | None
    detail: str = ""


#: One verdict per executable file: (real path, mtime in ns, size) -> ImplicitSupport.
_PROBED: dict[tuple[str, int, int], ImplicitSupport] = {}

#: A two-cell model whose LAK file carries IMPLICIT, written as plain MF6 input so the
#: probe needs nothing but the executable.
_PROBE_FILES = {
    "mfsim.nam": "BEGIN options\nEND options\nBEGIN timing\n  TDIS6 p.tdis\nEND timing\n"
                 "BEGIN models\n  gwf6 p.nam p\nEND models\nBEGIN exchanges\nEND exchanges\n"
                 "BEGIN solutiongroup 1\n  ims6 p.ims p\nEND solutiongroup\n",
    "p.tdis": "BEGIN options\nEND options\nBEGIN dimensions\n  NPER 1\nEND dimensions\n"
              "BEGIN perioddata\n  1.0 1 1.0\nEND perioddata\n",
    "p.ims": "BEGIN options\nEND options\nBEGIN nonlinear\n  OUTER_MAXIMUM 50\nEND nonlinear\n"
             "BEGIN linear\n  INNER_MAXIMUM 50\n  LINEAR_ACCELERATION BICGSTAB\nEND linear\n",
    "p.nam": "BEGIN options\nEND options\nBEGIN packages\n  DIS6 p.dis\n  NPF6 p.npf\n"
             "  IC6 p.ic\n  CHD6 p.chd\n  LAK6 p.lak\nEND packages\n",
    "p.dis": "BEGIN options\nEND options\nBEGIN dimensions\n  NLAY 1\n  NROW 1\n  NCOL 2\n"
             "END dimensions\nBEGIN griddata\n  delr\n    CONSTANT 10.0\n  delc\n    CONSTANT 10.0\n"
             "  top\n    CONSTANT 10.0\n  botm\n    CONSTANT 0.0\nEND griddata\n",
    "p.npf": "BEGIN options\nEND options\nBEGIN griddata\n  icelltype\n    CONSTANT 0\n"
             "  k\n    CONSTANT 1.0\nEND griddata\n",
    "p.ic": "BEGIN options\nEND options\nBEGIN griddata\n  strt\n    CONSTANT 9.0\nEND griddata\n",
    "p.chd": "BEGIN options\nEND options\nBEGIN dimensions\n  MAXBOUND 1\nEND dimensions\n"
             "BEGIN period 1\n  1 1 2 9.0\nEND period\n",
    "p.lak": "BEGIN options\n  IMPLICIT\nEND options\nBEGIN dimensions\n  NLAKES 1\n"
             "  NOUTLETS 0\n  NTABLES 0\nEND dimensions\nBEGIN packagedata\n  1 9.5 1\n"
             "END packagedata\nBEGIN connectiondata\n  1 1 1 1 1 VERTICAL 1.0 0.0 0.0 0.0 0.0\n"
             "END connectiondata\n",
}


def _resolve(executable: str) -> str | None:
    """The executable's path: a file given by path, else ``PATH``'s answer."""

    candidate = Path(os.path.expanduser(str(executable)))
    if candidate.is_file():
        return str(candidate)
    return shutil.which(str(executable))


def implicit_support(executable: str = "mf6") -> ImplicitSupport:
    """Probe an MF6 executable, once, for LAK ``IMPLICIT``.

    Runs ``<exe> --version``, then a two-cell model whose LAK file carries ``IMPLICIT``
    in a temporary directory, and caches the verdict per executable FILE -- its real
    path, modification time and size -- so a rebuilt MF6 is probed again and an
    unchanged one never is. Takes well under a second.

    Parameters
    ----------
    executable : str, default "mf6"
        A name looked up on ``PATH``, or a path to the executable.

    Returns
    -------
    ImplicitSupport
        ``supported`` is ``None`` (and nothing is cached) when the executable is not
        found or cannot be started. Nothing is raised: an unreadable, unstartable or
        timed-out executable is a verdict (``detail`` says which), not an error.

    Examples
    --------
    An executable that cannot be found gets no verdict; on a 6.8 build, ``supported``
    is ``True`` and ``version`` reads like ``"mf6: 6.8.0.dev0+0418efe4"``:

    >>> verdict = implicit_support("/no/such/mf6")
    >>> verdict.supported is None, verdict.detail
    (True, 'not found on PATH or as a file')
    """

    path = _resolve(executable)
    if path is None:
        return ImplicitSupport(str(executable), None, None, None, "not found on PATH or as a file")
    try:
        stat = os.stat(path)
    except OSError as error:
        return ImplicitSupport(str(executable), path, None, None, f"cannot be read: {error}")
    key = (os.path.realpath(path), stat.st_mtime_ns, stat.st_size)
    if key in _PROBED:
        return _PROBED[key]
    try:
        version = subprocess.run(
            [path, "--version"], capture_output=True, text=True, timeout=60
        ).stdout.strip() or None
        with tempfile.TemporaryDirectory(prefix="lak_implicit_probe_") as root:
            for name, text in _PROBE_FILES.items():
                Path(root, name).write_text(text)
            finished = subprocess.run([path], cwd=root, capture_output=True, text=True, timeout=120)
            listing = Path(root, "mfsim.lst")
            text = finished.stdout + (listing.read_text() if listing.exists() else "")
    except (OSError, subprocess.SubprocessError) as error:
        logger.debug("the MF6 IMPLICIT probe could not run %s: %s", path, error)
        return ImplicitSupport(str(executable), path, None, None, f"could not be run: {error}")
    supported = "Normal termination" in text
    detail = "" if supported else " / ".join(
        line.strip() for line in text.splitlines()
        if "ERROR" in line.upper() or "IMPLICIT" in line.upper() or "UNKNOWN" in line.upper()
    )[:500]
    verdict = ImplicitSupport(str(executable), path, version, supported, detail)
    _PROBED[key] = verdict
    return verdict


def implicit_packages(simulation) -> list[str]:
    """The file names of every LAK package in ``simulation`` that writes ``IMPLICIT``.

    Parameters
    ----------
    simulation : flopy.mf6.MFSimulation
        A built or loaded simulation.

    Returns
    -------
    list of str
        One file name per package, in model and package order; empty when none asks.

    Examples
    --------
    >>> implicit_packages(flopy.mf6.MFSimulation())
    []
    """

    found = []
    for name in getattr(simulation, "model_names", None) or ():
        model = simulation.get_model(name)
        for package in getattr(model, "packagelist", ()) or ():
            if getattr(package, "writes_lak_implicit", False):
                found.append(str(package.filename))
    return found


def warn_if_implicit_unsupported(simulation, executable: str | None = None) -> ImplicitSupport | None:
    """Warn when ``simulation`` asks for LAK ``IMPLICIT`` and its MF6 cannot read it.

    Nothing happens for a simulation without an implicit lake. Otherwise the executable
    -- ``executable``, else the simulation's own ``exe_name`` -- is probed
    (:func:`implicit_support`, cached) and a ``UserWarning`` names it and its version
    when it refused the option or could not be found. It never refuses: the run goes
    ahead, and an MF6 without the option stops reading the LAK file with its own
    message, which ends up in the run's report.

    Parameters
    ----------
    simulation : flopy.mf6.MFSimulation
        The simulation about to run.
    executable : str, optional
        The executable it will run with; ``None`` reads ``simulation.exe_name``.

    Returns
    -------
    ImplicitSupport or None
        The probe's verdict, or ``None`` when no package asks for ``IMPLICIT``.

    Examples
    --------
    >>> warn_if_implicit_unsupported(flopy.mf6.MFSimulation()) is None
    True
    """

    packages = implicit_packages(simulation)
    if not packages:
        return None
    exe = executable if executable is not None else getattr(simulation, "exe_name", None) or "mf6"
    verdict = implicit_support(exe)
    listing = ", ".join(repr(name) for name in packages)
    if verdict.supported is None:
        message = (
            f"LAK package(s) {listing} ask for IMPLICIT, and the MF6 executable {exe!r} "
            f"{verdict.detail}, so whether it accepts the option could not be checked. "
            "IMPLICIT is in MF6's develop branch and the 6.8 line only; an MF6 without it "
            "stops reading the LAK file."
        )
    elif not verdict.supported:
        message = (
            f"LAK package(s) {listing} ask for IMPLICIT, and the MF6 executable at "
            f"{verdict.path} ({verdict.version or 'version unknown'}) does not accept it"
            + (f" -- it said: {verdict.detail}" if verdict.detail else "")
            + ". The run goes ahead and MF6 will stop reading the LAK file with its own "
            "message. IMPLICIT is in MF6's develop branch (commit e2dfc84c) and the 6.8 "
            "line, in no earlier release: point the run at such a build, or declare the "
            "lake with implicit=False."
        )
    else:
        return verdict
    logger.warning(message)
    warnings.warn(message, UserWarning, stacklevel=stacklevel_outside_package())
    return verdict


__all__ = [
    "IMPLICIT",
    "ImplicitSupport",
    "ModflowGwflakImplicit",
    "has_implicit",
    "implicit_packages",
    "implicit_support",
    "insert_implicit",
    "keep_implicit",
    "warn_if_implicit_unsupported",
]
