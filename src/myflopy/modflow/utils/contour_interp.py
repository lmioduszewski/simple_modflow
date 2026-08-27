"""Interpolate an elevation raster surface from vector contours using GRASS GIS.

This wraps GRASS ``r.surf.contour`` to build a raster elevation surface from
drawn elevation contours -- the contour-following interpolation myflopy's
``Surface.from_contours`` relies on.

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
import threading
from pathlib import Path

from myflopy._logging import get_logger

logger = get_logger(__name__)

#: Executable names a POSIX GRASS install puts on ``PATH``. Most distributions
#: ship the unversioned ``grass``; the versioned names are the fallbacks for
#: installs that keep several side by side.
_GRASS_EXECUTABLES = ("grass", "grass84", "grass83", "grass8", "grass78")


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


class ContourSurfaceInterpolator:
    """Create an interpolated raster surface from vector elevation contours.

    Parameters mirror the GRASS workflow: import the contours, set a region from
    a raster or vector, rasterize the contours by their elevation attribute,
    interpolate with ``r.surf.contour``, optionally clip, and write a GeoTIFF.
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
        epsg: str = "2927",
        grass_bin: Path | str | None = None,
        grassdata: Path | str | None = None,
        location: str = "myflopy_interp",
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
        self.grassdata = Path(grassdata) if grassdata else Path.home() / "grassdata"
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

    def run(self) -> Path:
        """Run the full interpolation and return the written GeoTIFF path."""

        r, g, v, gsetup = _grass_modules(self.grass_bin)
        with _teed_console(self.progress):
            return self._run_grass(r, g, v, gsetup)

    def _run_grass(self, r, g, v, gsetup) -> Path:
        """The GRASS pipeline itself, with the console already routed if asked."""

        self._start_session(gsetup)
        v.in_ogr(input=self.contours.as_posix(), output=self._vect, overwrite=True, flags="o")
        self._set_region(r, g, v)
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

    def _set_region(self, r, g, v) -> None:
        """Set the GRASS computational region + mask from the region raster or vector."""

        if self.region_raster:
            r.in_gdal(overwrite=True, input=str(self.region_raster), output="Region", flags="o")
        elif self.region_vector:
            v.in_ogr(overwrite=True, input=str(self.region_vector), output="VRegion", flags="o")
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

    Convenience wrapper around :class:`ContourSurfaceInterpolator`.
    """

    return ContourSurfaceInterpolator(contours, out=out, **kwargs).run()
