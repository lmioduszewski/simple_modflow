"""Interpolate an elevation raster surface from vector contours using GRASS GIS.

This wraps GRASS ``r.surf.contour`` to build a raster elevation surface from
drawn elevation contours -- the contour-following interpolation myflopy's
``Surface.from_contours`` relies on.

GRASS works in ONE frame, its location's CRS, and imports with ``-o`` (no
projection check of its own). So every input is brought into that frame first
(ledger 602): :func:`contour_frame` reads each file's CRS and the
location's, :meth:`ContourFrame.refuse` refuses what cannot be reconciled, and
:func:`stage_contour_inputs` writes the contours and a region polygon,
reprojected, as the GeoPackage copies GRASS is actually given.

GRASS is an optional **system** dependency (it is not installed via pip), so its
Python modules are imported lazily. Importing this module never requires GRASS;
only actually running an interpolation does.
"""

from __future__ import annotations

import contextlib
import os
import shutil
import subprocess
import sys
import tempfile
import threading
from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
import shapely
from pyproj import CRS, Transformer
from pyproj.exceptions import ProjError

from myflopy._defaults import DEFAULT_EPSG
from myflopy._logging import get_logger
from myflopy.areas import _same_crs
from myflopy.modflow.utils.datatypes.readers import read_vector_layer

logger = get_logger(__name__)

#: Executable names a POSIX GRASS install puts on ``PATH``. Most distributions
#: ship the unversioned ``grass``; the versioned names are the fallbacks for
#: installs that keep several side by side.
_GRASS_EXECUTABLES = ("grass", "grass84", "grass83", "grass8", "grass78")

#: The GRASS location every contour surface shares unless ``location=`` names
#: another. One location has ONE CRS, set when it is created.
DEFAULT_LOCATION = "myflopy_interp"


def _default_grassdata() -> Path:
    """The GRASS database directory used when ``grassdata=`` is not given."""

    return Path.home() / "grassdata"


def _import_grass_modules():
    """The bare imports, so the retry after a ``sys.path`` fix is one call."""

    import grass.script.setup as gsetup
    from grass.pygrass.modules.shortcuts import general as g
    from grass.pygrass.modules.shortcuts import raster as r
    from grass.pygrass.modules.shortcuts import vector as v

    return r, g, v, gsetup


def _grass_python_path(grass_bin: Path | str) -> Path | None:
    """Ask a GRASS launcher where its Python bindings live, or ``None``.

    ``grass --config python_path`` prints the directory directly; older
    launchers only answer ``--config path`` (the install prefix), whose
    ``etc/python`` is the same place.
    """

    for option, suffix in (("python_path", None), ("path", "etc/python")):
        try:
            result = subprocess.run(
                [str(grass_bin), "--config", option], capture_output=True, text=True
            )
        except (OSError, subprocess.SubprocessError):
            logger.debug("could not run GRASS launcher %s", grass_bin, exc_info=True)
            return None
        if result.returncode != 0 or not result.stdout.strip():
            continue
        path = Path(result.stdout.strip().splitlines()[0])
        if suffix:
            path = path / suffix
        if path.is_dir():
            return path
    return None


def _add_grass_python_path(grass_bin: Path | str | None) -> Path | None:
    """Put the GRASS Python bindings on ``sys.path``; return the path added."""

    try:
        launcher = Path(grass_bin) if grass_bin else _default_grass_bin()
    except ValueError:
        logger.debug("no GRASS launcher to ask for the bindings path", exc_info=True)
        return None
    path = _grass_python_path(launcher)
    if path is None:
        return None
    if str(path) not in sys.path:
        sys.path.append(str(path))
    return path


def _grass_modules(grass_bin: Path | str | None = None):
    """Import GRASS Python modules lazily, with a clear error if unavailable.

    The bindings ship *inside* the GRASS install (``<prefix>/etc/python``) and
    are not on a normal interpreter's ``sys.path``, so a first failure is not
    the answer: the launcher is asked where they live and the import retried.
    That is what spares callers a hand-set ``PYTHONPATH``.
    """

    try:
        return _import_grass_modules()
    except Exception as error:  # noqa: BLE001 - see below  # pragma: no cover
        # Broad on purpose, and it does not swallow: every failure to load the
        # optional GRASS SYSTEM dependency is re-raised as one actionable
        # ImportError with the cause attached. GRASS's own import chain runs
        # ctypes bindings and shells out to the launcher, so it raises far more
        # than ImportError -- OSError, RuntimeError, CalledProcessError -- and
        # the caller needs the same message for all of them.
        last_error = error

    added = _add_grass_python_path(grass_bin)  # pragma: no cover - needs GRASS
    if added is not None:  # pragma: no cover - needs GRASS
        logger.debug("added GRASS Python bindings at %s to sys.path", added)
        try:
            return _import_grass_modules()
        except Exception as error:  # noqa: BLE001 - same reason as above
            last_error = error

    raise ImportError(  # pragma: no cover - needs GRASS
        "GRASS GIS Python modules are required for contour interpolation. "
        "Install GRASS GIS and run from a GRASS-enabled environment. The "
        "bindings live in <grass prefix>/etc/python and are found via the "
        "launcher, so set GRASS_BIN (or pass grass_bin=...) if GRASS is "
        "installed somewhere this could not discover."
    ) from last_error


def _grass_search_dirs() -> list[Path]:
    """Return likely locations of a GRASS launcher (OSGeo4W and QGIS bundles)."""

    dirs = [
        Path.home() / "AppData" / "Local" / "Programs" / "OSGeo4W" / "bin",
        Path("C:/OSGeo4W/bin"),
        Path("C:/OSGeo4W64/bin"),
    ]
    dirs.extend(sorted(Path("C:/Program Files").glob("QGIS */bin"), reverse=True))
    return dirs


def _find_grass_launcher(search_dirs) -> Path | None:
    """Return the highest-versioned ``grass*.bat`` launcher found, if any.

    Windows-shaped: OSGeo4W and the QGIS bundles ship ``.bat`` launchers in a
    known directory. The ``python-grass*.bat`` wrapper is skipped in favor of
    the main launcher. POSIX installs are found by :func:`_find_grass_on_path`.
    """

    for directory in search_dirs:
        directory = Path(directory)
        if not directory.is_dir():
            continue
        for candidate in sorted(directory.glob("grass*.bat"), reverse=True):
            if candidate.is_file() and not candidate.name.startswith("python"):
                return candidate
    return None


def _find_grass_on_path() -> Path | None:
    """Return the first GRASS executable on ``PATH``, if any (POSIX installs).

    A Linux or macOS install is a plain executable -- ``/usr/bin/grass``,
    ``/opt/homebrew/bin/grass`` -- never a ``.bat``, and it is already on
    ``PATH``, so ``which`` is the whole search.
    """

    for name in _GRASS_EXECUTABLES:
        found = shutil.which(name)
        if found:
            return Path(found)
    return None


def _default_grass_bin() -> Path:
    """Resolve a GRASS launcher: ``GRASS_BIN`` env var, then auto-discovery.

    ``GRASS_BIN`` always wins. Otherwise the Windows bundle directories are
    globbed for a ``grass*.bat``, and on every other platform ``PATH`` is
    searched for a GRASS executable.
    """

    env = os.environ.get("GRASS_BIN")
    if env:
        return Path(env)
    found = _find_grass_launcher(_grass_search_dirs())
    if found is None and os.name != "nt":
        found = _find_grass_on_path()
    if found is not None:
        return found
    if os.name == "nt":
        example, install = "the path to grass84.bat", "install GRASS via OSGeo4W"
    else:
        example, install = (
            "/usr/bin/grass",
            "install GRASS GIS with your package manager (it is a SYSTEM "
            "dependency, not a pip one) -- PATH was searched for "
            + ", ".join(_GRASS_EXECUTABLES),
        )
    raise ValueError(
        f"Could not find a GRASS launcher. Pass grass_bin=... (e.g. {example}), "
        f"set the GRASS_BIN environment variable, or {install}."
    )


@contextlib.contextmanager
def _teed_console(enabled: bool):
    """Route the C-level output of GRASS subprocesses into Python's stdout.

    GRASS modules -- ``v.to.rast``, ``r.surf.contour``, ``r.out.gdal`` -- are
    subprocesses that write their progress percentages to the file descriptors
    they inherit. In a terminal that already lands on screen, which is why this
    is off by default and a no-op there in every way that matters.

    In Jupyter it does NOT: fd 1 belongs to the kernel's console, not the cell,
    and `contextlib.redirect_stdout` cannot help because it rebinds
    ``sys.stdout`` while the writing happens below Python entirely. So the fds
    themselves are swapped for a pipe and a reader thread re-emits what arrives
    through ``sys.stdout``, which in a kernel is the object that reaches the cell.

    Read in small chunks rather than by line: GRASS separates percentages with
    carriage returns and no newline, so line buffering would hold the whole run
    back and deliver it at the end, which is the opposite of progress.
    """

    if not enabled:
        yield
        return

    # Where to re-emit. In a terminal `sys.stdout` IS fd 1, which is about to
    # become the pipe -- writing there would feed the reader its own output --
    # so the duplicate of the original fd is used instead. A Jupyter OutStream
    # has no real fileno and raises, which is exactly how it is recognised.
    try:
        writes_to_fd1 = sys.stdout.fileno() == 1
    except (AttributeError, OSError, ValueError):
        writes_to_fd1 = False

    read_fd, write_fd = os.pipe()
    saved_out, saved_err = os.dup(1), os.dup(2)
    sink = os.fdopen(os.dup(saved_out), "w", buffering=1) if writes_to_fd1 else sys.stdout

    def pump() -> None:
        with os.fdopen(read_fd, "rb", 0) as reader:
            while chunk := reader.read(128):
                sink.write(chunk.decode("utf-8", "replace"))
                sink.flush()

    reader_thread = threading.Thread(target=pump, name="grass-progress", daemon=True)
    reader_thread.start()
    try:
        os.dup2(write_fd, 1)
        os.dup2(write_fd, 2)
        os.close(write_fd)
        yield
    finally:
        # Restoring drops the last references to the pipe's write end, which is
        # what lets the reader see EOF and finish.
        os.dup2(saved_out, 1)
        os.dup2(saved_err, 2)
        os.close(saved_out)
        os.close(saved_err)
        reader_thread.join(timeout=10)
        if sink is not sys.stdout:
            sink.close()


def _location_crs(location_path: Path | str) -> str | None:
    """The CRS a GRASS location was created in, or ``None`` when it does not exist.

    Read from the files GRASS writes into ``PERMANENT``: ``PROJ_SRID`` (GRASS 8,
    ``EPSG:2927``), then ``PROJ_WKT``, then GRASS 7's ``PROJ_EPSG`` (``epsg: 2927``).
    A location holding none of them (an XY location, or one not created yet) has no
    CRS to compare against and reads ``None``.
    """

    permanent = Path(location_path) / "PERMANENT"
    srid = permanent / "PROJ_SRID"
    if srid.is_file() and srid.read_text().strip():
        return srid.read_text().strip()
    wkt = permanent / "PROJ_WKT"
    if wkt.is_file() and wkt.read_text().strip():
        return wkt.read_text().strip()
    epsg = permanent / "PROJ_EPSG"
    if epsg.is_file():
        for line in epsg.read_text().splitlines():
            key, _, value = line.partition(":")
            if key.strip().lower() == "epsg" and value.strip():
                return f"EPSG:{value.strip()}"
    return None


#: Two CRSs are ONE frame here when carrying a point from one to the other moves it by
#: no more than this fraction of the target CRS's area of use (ledger 603):
#: about 0.019 ft over EPSG:2927's 1.9 million ft. pyproj's ``equals`` -- the house
#: ``_same_crs`` -- is exact, and the user's own lidar DEM and ground-surface contours
#: carry an EPSG:2927 written with no authority code, a false easting of
#: 1640416.66666667 where EPSG's definition reads 1640416.667, and its standard
#: parallels in the other order: ``equals`` says no, and the points move 0.00033 ft.
#: A real datum change moves feet (EPSG:2286 -> 2927: 4.24 ft) and stays a different
#: frame.
_SAME_FRAME_FRACTION = 1e-8


def _same_frame(crs, target) -> bool:
    """Whether ``crs`` is ``target``'s frame: the same CRS, or one ``to_crs`` would not move.

    Measured, not judged by name: a 5 x 5 lattice over ``target``'s area of use, in
    ``target``'s coordinates, is read as ``crs`` coordinates and carried into ``target``
    by the transformation ``to_crs`` would use; the two are one frame when no point
    moves more than ``_SAME_FRAME_FRACTION`` of the area's span. A pair the transform
    cannot carry (degrees against feet), or a target with no area of use, is not.
    """

    if _same_crs(crs, target):
        return True
    target_crs = CRS.from_user_input(target)
    area = target_crs.area_of_use
    if area is None:
        return False
    lon, lat = np.meshgrid(np.linspace(area.west, area.east, 5),
                           np.linspace(area.south, area.north, 5))
    try:
        xs, ys = Transformer.from_crs("EPSG:4326", target_crs, always_xy=True).transform(
            lon.ravel(), lat.ravel())
        moved_x, moved_y = Transformer.from_crs(
            CRS.from_user_input(crs), target_crs, always_xy=True).transform(xs, ys)
    except ProjError:
        logger.debug("no transformation from %s to %s: not one frame", crs, target,
                     exc_info=True)
        return False
    moved = np.hypot(np.asarray(moved_x) - xs, np.asarray(moved_y) - ys)
    span = max(np.ptp(xs), np.ptp(ys))
    return bool(np.isfinite(moved).all() and moved.max() <= _SAME_FRAME_FRACTION * span)


def _crs_label(crs) -> str:
    """A CRS as a message names it: ``EPSG:n`` when it has a code, else its name."""

    parsed = CRS.from_user_input(crs)
    code = parsed.to_epsg()
    return f"EPSG:{code}" if code is not None else parsed.name


@dataclass(frozen=True)
class ContourInput:
    """One file a contour surface reads, and the CRS it carries.

    Attributes
    ----------
    role : {'contours', 'region_vector', 'region_raster'}
        What the file is to the interpolation; messages name it by this.
    path : pathlib.Path
        The file, as the surface declared it.
    crs : str or None
        Its CRS as ``pyproj``'s ``to_string()`` spells it (``"EPSG:4326"``), or
        ``None`` when the file carries none -- or does not exist yet.
    exists : bool
        Whether the file exists. A missing file is never refused here: a cache is
        reused under the stale rule, and GRASS names the path if it must run.
    """

    role: str
    path: Path
    crs: str | None
    exists: bool

    @property
    def is_vector(self) -> bool:
        """True for the contours and a region polygon, the two inputs that are reprojected."""

        return self.role != "region_raster"


@dataclass(frozen=True)
class ContourFrame:
    """The frame a GRASS contour interpolation works in, and each input's own CRS.

    Built by :func:`contour_frame`, which reads CRSs and nothing else. GRASS works in
    its LOCATION's CRS and imports with ``-o``, so it never compares the two itself;
    this is where they are compared (ledger 602).

    Attributes
    ----------
    target : str
        ``"EPSG:<epsg>"``, the CRS the inputs are reprojected into and the raster is
        written in.
    location : pathlib.Path
        The GRASS location directory, ``<grassdata>/<location>``.
    location_crs : str or None
        The CRS that location was created in, as GRASS recorded it; ``None`` while
        the location does not exist (it is created in ``target`` on first use).
    inputs : tuple of ContourInput
        The contours, then the region (a polygon or a raster), each with its CRS.
    """

    target: str
    location: Path
    location_crs: str | None
    inputs: tuple[ContourInput, ...]

    def _input(self, role: str) -> ContourInput | None:
        """The input playing ``role``, or ``None``."""

        return next((item for item in self.inputs if item.role == role), None)

    @property
    def reprojected(self) -> dict[str, str]:
        """``{role: crs}`` for each existing vector input in a CRS other than ``target``.

        Empty when every input is already in the location's frame -- the same CRS, or
        one reprojecting would not move (ledger 603) -- or carries no CRS,
        which is READ as that frame, the rule every myflopy vector reader follows; so a
        surface whose files never move keeps its cache key unchanged.
        """

        return {
            item.role: item.crs
            for item in self.inputs
            if item.is_vector and item.exists and item.crs is not None
            and not _same_frame(item.crs, self.target)
        }

    def refuse(self) -> None:
        """Raise for a frame no reprojection can reconcile; otherwise return ``None``.

        Run before any cache is read and before GRASS starts, so a refusal needs no
        GRASS and a cache made under such a frame is never reused in silence.

        Raises
        ------
        ValueError
            The GRASS location already exists in a CRS other than ``target`` -- GRASS
            reads ``epsg=`` only when it CREATES a location, so the raster would be
            written in the location's CRS whatever ``epsg=`` says (measured: ``"32610"``
            on an EPSG:2927 location wrote a raster labelled EPSG:2927 holding metre
            coordinates). Or a ``region_raster`` is in another CRS: a raster region is
            read as it is, never reprojected.
        """

        if self.location_crs is not None and not _same_frame(self.location_crs, self.target):
            have = _crs_label(self.location_crs)
            code = self.target.split(":", 1)[1]
            # `epsg=` takes a bare code, so it is offered only when the location has one.
            or_epsg = (f", or epsg={have.split(':', 1)[1]!r}" if have.startswith("EPSG:")
                       else "")
            raise ValueError(
                f"The GRASS location '{self.location}' already exists in {have}, and "
                f"epsg={code!r} asks for {self.target}. GRASS works in the location's own "
                "CRS and never reads epsg= again, so the raster would be written in "
                f"{have}. Pass location='{self.location.name}_{code}' (a location of its "
                f"own, created on first use){or_epsg}."
            )
        raster = self._input("region_raster")
        if (
            raster is not None and raster.exists and raster.crs is not None
            and not _same_frame(raster.crs, self.target)
        ):
            raise ValueError(
                f"region_raster '{raster.path.name}' is in {_crs_label(raster.crs)}, and "
                f"GRASS works in its location's {self.target}. A raster region is read as "
                "it is, never reprojected, so its extent would land in the wrong place. "
                "Pass region_vector=<the domain polygon>, which is reprojected for you, or "
                f"reproject the raster to {self.target} first."
            )


def _vector_crs(path: Path) -> str | None:
    """A vector file's CRS as ``to_string()`` spells it, reading one feature."""

    crs = gpd.read_file(path, rows=1).crs
    return None if crs is None else crs.to_string()


def _raster_crs(path: Path) -> str | None:
    """A raster's CRS as ``to_string()`` spells it, or ``None`` when it carries none."""

    with rasterio.open(path) as handle:
        return None if handle.crs is None else handle.crs.to_string()


def contour_frame(
    contours: Path | str,
    *,
    region_vector: Path | str | None = None,
    region_raster: Path | str | None = None,
    epsg: str | int = str(DEFAULT_EPSG),
    grassdata: Path | str | None = None,
    location: str | None = None,
) -> ContourFrame:
    """Read the CRS of a contour interpolation's inputs and of its GRASS location.

    Reads CRSs only -- one feature of each vector file, a raster's header, the
    location's ``PERMANENT`` projection files -- so it is cheap enough to run on every
    resolution of a contour surface, and it raises nothing for a file that does not
    exist yet (``cache_status()`` is asked about paths not yet written).

    Parameters
    ----------
    contours : str or pathlib.Path
        The contour-line vector file.
    region_vector : str or pathlib.Path, optional
        The region polygon, when the region is a vector. Default None.
    region_raster : str or pathlib.Path, optional
        The region raster, when the region is a raster. Default None.
    epsg : str or int, default "2927"
        The EPSG code (bare, as ``"2927"`` or ``2927``) of the GRASS location: the
        CRS the vector inputs are reprojected into and the raster is written in.
    grassdata : str or pathlib.Path, optional
        The GRASS database directory. ``None`` (the default) is ``~/grassdata``, as
        :class:`ContourSurfaceInterpolator` uses.
    location : str, optional
        The GRASS location's name inside ``grassdata``. ``None`` (the default) is
        ``"myflopy_interp"``.

    Returns
    -------
    ContourFrame
        The target CRS, the location and its CRS (``None`` while it does not exist),
        and each input with its own CRS.
    """

    grassdata = Path(grassdata) if grassdata is not None else _default_grassdata()
    location_path = grassdata / (location if location is not None else DEFAULT_LOCATION)
    inputs = []
    for role, path, reader in (
        ("contours", contours, _vector_crs),
        ("region_vector", region_vector, _vector_crs),
        ("region_raster", region_raster, _raster_crs),
    ):
        if path is None:
            continue
        path = Path(path)
        exists = path.exists()
        inputs.append(ContourInput(role, path, reader(path) if exists else None, exists))
    return ContourFrame(
        target=f"EPSG:{epsg}",
        location=location_path,
        location_crs=_location_crs(location_path),
        inputs=tuple(inputs),
    )


def _in_target(frame: ContourFrame, item: ContourInput) -> gpd.GeoDataFrame:
    """A vector input read and brought into ``frame.target``.

    A file in another CRS is reprojected (``to_crs``, the rule every ``.gpkg`` source
    follows); one carrying NO CRS is read as already in the location's frame -- the
    rule every myflopy vector reader follows (decided 2026-10-16, ledger
    602) -- and labelled so; and one in a CRS that reprojecting would not
    move (:func:`_same_frame`, ledger 603) keeps its coordinates and is
    labelled ``frame.target`` too.
    """

    frame_gdf = read_vector_layer(item.path)
    if frame_gdf.crs is None:
        return frame_gdf.set_crs(frame.target)
    if not _same_frame(frame_gdf.crs, frame.target):
        try:
            return frame_gdf.to_crs(frame.target)
        except ProjError as error:
            # pyproj's own message ("Error creating Transformer from CRS.") names
            # neither CRS nor the file -- a CAD export's engineering CRS does this.
            raise ValueError(
                f"{item.role} '{item.path.name}' is in {_crs_label(frame_gdf.crs)}, which "
                f"cannot be carried into the GRASS location's {frame.target} (ledger "
                "604). Give the file the CRS it was drawn in: "
                "gpd.read_file(path).set_crs(<its EPSG>, allow_override=True)"
                ".to_file(path)."
            ) from error
    return frame_gdf.set_crs(frame.target, allow_override=True)


def _keep_columns(gdf: gpd.GeoDataFrame, names) -> gpd.GeoDataFrame:
    """``gdf`` holding only its geometry and the columns ``names`` match, ignoring case."""

    wanted = {str(name).lower() for name in names}
    keep = [column for column in gdf.columns
            if column != gdf.geometry.name and str(column).lower() in wanted]
    return gdf[[*keep, gdf.geometry.name]]


def stage_contour_inputs(
    frame: ContourFrame, scratch: Path | str, *, z_field: str | None = None
) -> tuple[Path, Path | None]:
    """Write the contours and a region polygon, in the location's CRS, for GRASS.

    The staged GeoPackage copies are the only files GRASS imports. A copy of a file
    already in the frame imports identically to the file itself (measured: 17,688
    pixels, max |diff| 0.0), so staging always happens, one code path. A contour file
    holding several layers is read as its FIRST, as every myflopy vector reader reads
    one. A copy carries only what GRASS reads -- the region its geometry, the contours
    theirs and ``z_field`` -- because a GeoPackage refuses attributes another format
    holds: a shapefile's ``fid`` column that is not a unique integer (what QGIS's
    "Merge vector layers" writes) is GeoPackage's own key (ledger 604).

    Parameters
    ----------
    frame : ContourFrame
        From :func:`contour_frame`, already passed through :meth:`ContourFrame.refuse`.
    scratch : str or pathlib.Path
        An existing directory to write into; the caller removes it.
    z_field : str, optional
        The contours' elevation attribute, matched IGNORING CASE as GRASS's SQLite
        attribute table matches it (``z="ELEV"`` reads a column ``elev``). Given, the
        staged contours carry it and no other attribute, and a file without it is
        refused here, before GRASS starts. ``None`` (the default) keeps every
        attribute of the contours. The region keeps no attribute either way.

    Returns
    -------
    tuple of (pathlib.Path, pathlib.Path or None)
        The staged contours, and the staged region polygon (``None`` when the region
        is a raster, which GRASS reads as it is).

    Raises
    ------
    ValueError
        No contour reaches the region's extent -- GRASS's computational region is
        that box, so ``r.surf.contour`` would have nothing to interpolate from and
        runs at 0% indefinitely (ledger 149). The message gives both extents in the
        location's CRS and, when the contour file was moved or carries no CRS, what
        it was read as. Or the contours carry no ``z_field`` attribute (the message
        lists those they carry). Or a vector input's CRS cannot be carried into the
        location's at all (an engineering CRS, as a CAD export writes; ledger
        604).
    """

    scratch = Path(scratch)
    contours_item = frame._input("contours")
    contours = _in_target(frame, contours_item)
    if z_field is not None:
        staged = _keep_columns(contours, [z_field])
        if staged.shape[1] == 1:
            have = [str(column) for column in contours.columns
                    if column != contours.geometry.name]
            raise ValueError(
                f"Contours '{contours_item.path.name}' carry no attribute {z_field!r}, so "
                "GRASS has no elevation to rasterize. They carry: "
                f"{', '.join(map(repr, have)) or 'no attributes'}. Pass z=<that column> "
                "(z_field= on ContourSurfaceInterpolator)."
            )
        contours = staged
    region_vector = frame._input("region_vector")
    region_raster = frame._input("region_raster")
    staged_region = None
    # The raster first, as `_set_region` reads them: given both, GRASS's region is
    # the raster's, so the overlap is checked against it and no polygon is staged.
    if region_raster is not None:
        bounds = region_path = None
        if region_raster.exists:
            with rasterio.open(region_raster.path) as handle:
                bounds, region_path = tuple(handle.bounds), region_raster.path
    elif region_vector is not None:
        region = _keep_columns(_in_target(frame, region_vector), [])
        bounds, region_path = tuple(region.total_bounds), region_vector.path
        staged_region = scratch / "region.gpkg"
        region.to_file(staged_region, driver="GPKG")
    else:
        bounds = region_path = None
    if bounds is not None:
        box = shapely.box(*bounds)
        if not contours.geometry.intersects(box).any():
            def extent(values) -> str:
                return "(" + ", ".join(f"{value:,.6g}" for value in values) + ")"

            if contours_item.crs is None:
                read_as = (f" The contours' file carries no CRS and was read as "
                           f"{frame.target}.")
            elif not _same_frame(contours_item.crs, frame.target):
                read_as = f" The contours' file is in {_crs_label(contours_item.crs)}."
            else:
                read_as = ""
            raise ValueError(
                f"None of the {len(contours)} contours in '{contours_item.path.name}' "
                f"reach the extent of region '{region_path.name}', so GRASS has no "
                "contour to interpolate from -- r.surf.contour then runs at 0% "
                f"indefinitely (ledger 149). The contours cover "
                f"{extent(contours.total_bounds)} and the region "
                f"{extent(bounds)}, both in {frame.target}.{read_as} Check they are the "
                "same place in the same CRS."
            )
    staged_contours = scratch / "contours.gpkg"
    contours.to_file(staged_contours, driver="GPKG")
    return staged_contours, staged_region


class ContourSurfaceInterpolator:
    """Create an interpolated raster surface from vector elevation contours.

    Parameters mirror the GRASS workflow: import the contours, set a region from
    a raster or vector, rasterize the contours by their elevation attribute,
    interpolate with ``r.surf.contour``, optionally clip, and write a GeoTIFF.

    GRASS works in its LOCATION's CRS. :meth:`run` first brings the contours and a
    region polygon into it (:func:`contour_frame`, :func:`stage_contour_inputs`),
    and refuses before GRASS starts what it cannot reconcile (ledger
    602).

    Parameters
    ----------
    contours : str or pathlib.Path
        Vector file of elevation contour lines, read as its first layer. In a CRS
        other than ``epsg``'s it is reprojected into it; one carrying no CRS is read
        as already in it.
    out : str or pathlib.Path
        Where the GeoTIFF is written, in the location's CRS.
    region_raster : str or pathlib.Path, optional
        Raster whose extent is the GRASS region. It must already be in the
        location's CRS (or carry none, and be read as in it): a raster is read as it
        is, so one in another CRS is refused. One of this or ``region_vector`` is
        needed.
    region_vector : str or pathlib.Path, optional
        Polygon whose extent (and shape, with ``clip``) is the GRASS region,
        reprojected like ``contours``. Used when ``region_raster`` is not given.
    z_field : str, default "Elev"
        The contours' elevation attribute, matched ignoring case as GRASS matches
        it. Contours without it are refused before GRASS starts.
    resolution : float, default 4
        Cell size, in the units of the location's CRS.
    epsg : str, default "2927"
        Bare EPSG code of the GRASS location, used to create it on first use. An
        existing location keeps its own CRS, so one that disagrees is refused,
        naming ``location=``.
    grass_bin : str or pathlib.Path, optional
        The GRASS launcher; ``None`` finds it (``GRASS_BIN``, then discovery).
    grassdata : str or pathlib.Path, optional
        The GRASS database directory; ``None`` is ``~/grassdata``.
    location : str, default "myflopy_interp"
        The GRASS location's name inside ``grassdata``.
    clip : bool, default True
        Null the result outside the region (the polygon itself, for a vector).
    progress : bool, default False
        Re-emit GRASS's own progress through ``sys.stdout`` (see
        :func:`_teed_console`).
    """

    def __init__(
        self,
        contours: Path | str,
        *,
        out: Path | str,
        region_raster: Path | str | None = None,
        region_vector: Path | str | None = None,
        z_field: str = "Elev",
        resolution: float = 4,
        epsg: str = str(DEFAULT_EPSG),
        grass_bin: Path | str | None = None,
        grassdata: Path | str | None = None,
        location: str = DEFAULT_LOCATION,
        clip: bool = True,
        progress: bool = False,
    ):
        """Configure a GRASS contour-to-raster interpolation (see the class docstring for parameters)."""

        self.contours = Path(contours)
        self.out = Path(out)
        self.region_raster = region_raster
        self.region_vector = region_vector
        self.z_field = z_field
        self.resolution = resolution
        self.epsg = str(epsg)
        self.grass_bin = Path(grass_bin) if grass_bin else None
        self.grassdata = Path(grassdata) if grassdata else _default_grassdata()
        self.location = location
        self.clip = clip
        self.progress = progress
        self.session = None
        self._vect = "vectContours"
        self._rast = "rastContours"

    @property
    def location_path(self) -> Path:
        """The GRASS location directory for this interpolation (``grassdata/location``)."""

        return self.grassdata / self.location

    def frame(self) -> ContourFrame:
        """The frame this interpolation works in: :func:`contour_frame` of its inputs."""

        return contour_frame(
            self.contours, region_vector=self.region_vector, region_raster=self.region_raster,
            epsg=self.epsg, grassdata=self.grassdata, location=self.location,
        )

    def run(self) -> Path:
        """Run the full interpolation and return the written GeoTIFF path.

        The inputs' frame is checked and the vector inputs staged in the location's
        CRS before GRASS is imported or started, so every refusal needs no GRASS.

        Raises
        ------
        ValueError
            From :meth:`ContourFrame.refuse` (an existing location in another CRS than
            ``epsg``, a region raster in another CRS) or :func:`stage_contour_inputs`
            (no contour reaches the region's extent, the contours carry no
            ``z_field``, or an input's CRS cannot be carried into the location's).
        """

        frame = self.frame()
        frame.refuse()
        with tempfile.TemporaryDirectory(prefix="myflopy-contours-") as scratch:
            contours, region_vector = stage_contour_inputs(
                frame, scratch, z_field=self.z_field
            )
            r, g, v, gsetup = _grass_modules(self.grass_bin)
            with _teed_console(self.progress):
                return self._run_grass(
                    r, g, v, gsetup, contours=contours, region_vector=region_vector
                )

    def _run_grass(self, r, g, v, gsetup, *, contours: Path, region_vector: Path | None) -> Path:
        """The GRASS pipeline itself, with the console already routed if asked.

        ``contours``/``region_vector`` are the STAGED copies, already in the
        location's CRS: ``-o`` skips GRASS's own projection check, which names
        neither CRS, because the frame was checked before GRASS started.
        """

        self._start_session(gsetup)
        v.in_ogr(input=Path(contours).as_posix(), output=self._vect, overwrite=True, flags="o")
        self._set_region(r, g, v, region_vector=region_vector)
        v.to_rast(
            overwrite=True,
            input=self._vect,
            type="line",
            output=self._rast,
            use="attr",
            attribute_column=self.z_field,
        )
        r.surf_contour(input=self._rast, output="interpdSurface", overwrite=True)
        if self.clip:
            r.mapcalc(
                expression="interpdSurface = if(Region, interpdSurface, null())",
                overwrite=True,
            )
        r.out_gdal(
            input="interpdSurface",
            output=self.out.as_posix(),
            format="GTiff",
            overwrite=True,
        )
        self.session.finish()
        return self.out

    def _start_session(self, gsetup) -> None:
        """Open a GRASS session in the location, creating the location first if needed."""

        try:
            self.session = gsetup.init(self.grassdata, self.location)
        except (ValueError, subprocess.CalledProcessError):
            # ValueError is the NORMAL first-run path, not an error: GRASS 8.4's
            # `init` raises it for "Path <grassdata> does not exist" and for a
            # location that has not been created yet -- and `grassdata` defaults
            # to ~/grassdata, which most users will not have. CalledProcessError
            # comes from `get_install_path` querying the GRASS binary as a
            # subprocess, which fires on a half-configured install. Creating the
            # location and retrying is the answer to both.
            logger.debug(
                "no GRASS session at %s; creating the location first",
                self.location_path, exc_info=True,
            )
            self._create_location()
            self.session = gsetup.init(self.grassdata, self.location)

    def _create_location(self) -> None:
        """Create the GRASS location for this EPSG via the GRASS binary (raises on failure)."""

        grass_bin = self.grass_bin or _default_grass_bin()
        cmd = [str(grass_bin), "-c", "epsg:" + self.epsg, "-e", str(self.location_path)]
        result = subprocess.run(cmd, capture_output=True)
        if result.returncode != 0:
            raise RuntimeError(
                "Failed to create GRASS location: "
                + result.stderr.decode(errors="ignore")
            )

    def _set_region(self, r, g, v, *, region_vector: Path | None) -> None:
        """Set the GRASS computational region + mask from the region raster or vector.

        ``region_vector`` is the staged polygon; the raster is read as it is, which
        :meth:`ContourFrame.refuse` made safe by refusing one in another CRS.
        """

        if self.region_raster:
            r.in_gdal(overwrite=True, input=str(self.region_raster), output="Region", flags="o")
        elif region_vector is not None:
            v.in_ogr(overwrite=True, input=Path(region_vector).as_posix(), output="VRegion", flags="o")
            g.region(vector="VRegion", res=self.resolution, flags="p")
            v.to_rast(overwrite=True, input="VRegion", output="Region", use="value", value=1)
        else:
            raise ValueError(
                "Provide region_raster or region_vector to define the GRASS region."
            )
        g.region(raster="Region", res=self.resolution, flags="p")


def interpolate_contours_to_raster(
    contours: Path | str, *, out: Path | str, **kwargs
) -> Path:
    """Interpolate an elevation raster from vector contours via GRASS.

    Convenience wrapper around :class:`ContourSurfaceInterpolator`: ``contours``,
    ``out`` and every keyword (``**kwargs``) go to its constructor, documented there,
    and its :meth:`~ContourSurfaceInterpolator.run` is returned.
    """

    return ContourSurfaceInterpolator(contours, out=out, **kwargs).run()
