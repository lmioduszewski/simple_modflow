"""Interpolate an elevation raster surface from vector contours using GRASS GIS.

This wraps GRASS ``r.surf.contour`` to build a raster elevation surface from
drawn elevation contours -- the contour-following interpolation myflopy's
``Surface.from_contours`` relies on.

GRASS is an optional **system** dependency (it is not installed via pip), so its
Python modules are imported lazily. Importing this module never requires GRASS;
only actually running an interpolation does.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path


def _grass_modules():
    """Import GRASS Python modules lazily, with a clear error if unavailable."""

    try:
        import grass.script.setup as gsetup
        from grass.pygrass.modules.shortcuts import general as g
        from grass.pygrass.modules.shortcuts import raster as r
        from grass.pygrass.modules.shortcuts import vector as v
    except Exception as error:  # pragma: no cover - exercised only with GRASS
        raise ImportError(
            "GRASS GIS Python modules are required for contour interpolation. "
            "Install GRASS GIS and run from a GRASS-enabled environment."
        ) from error
    return r, g, v, gsetup


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

    The ``python-grass*.bat`` wrapper is skipped in favor of the main launcher.
    """

    for directory in search_dirs:
        directory = Path(directory)
        if not directory.is_dir():
            continue
        for candidate in sorted(directory.glob("grass*.bat"), reverse=True):
            if candidate.is_file() and not candidate.name.startswith("python"):
                return candidate
    return None


def _default_grass_bin() -> Path:
    """Resolve a GRASS launcher: GRASS_BIN env var, then auto-discovery."""

    env = os.environ.get("GRASS_BIN")
    if env:
        return Path(env)
    found = _find_grass_launcher(_grass_search_dirs())
    if found is not None:
        return found
    raise ValueError(
        "Could not find a GRASS launcher. Pass grass_bin=... (e.g. the path to "
        "grass84.bat), set the GRASS_BIN environment variable, or install GRASS "
        "via OSGeo4W."
    )


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
        self.session = None
        self._vect = "vectContours"
        self._rast = "rastContours"

    @property
    def location_path(self) -> Path:
        """The GRASS location directory for this interpolation (``grassdata/location``)."""

        return self.grassdata / self.location

    def run(self) -> Path:
        """Run the full interpolation and return the written GeoTIFF path."""

        r, g, v, gsetup = _grass_modules()
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
        except Exception:
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
