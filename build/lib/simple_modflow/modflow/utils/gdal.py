"""from osgeo import gdal, ogr
from shapely.wkt import loads
from shapely.geometry import Polygon
import geopandas as gpd


def get_contours_as_polygons(raster_path, contour_interval=10, union=True):

    # Open the raster file
    raster_ds = gdal.Open(raster_path)
    if not raster_ds:
        raise FileNotFoundError(f"Could not open raster file: {raster_path}")

    # Get the raster band (first band)
    band = raster_ds.GetRasterBand(1)

    # Create an in-memory data source for contours
    driver = ogr.GetDriverByName("MEM")
    contour_ds = driver.CreateDataSource("in_memory")
    contour_layer = contour_ds.CreateLayer("contours", geom_type=ogr.wkbLineString)

    # Add a field to store contour elevation
    field_defn = ogr.FieldDefn("elevation", ogr.OFTReal)
    contour_layer.CreateField(field_defn)

    # Generate contours
    gdal.ContourGenerate(
        band,
        contourInterval=contour_interval,  # Contour interval
        contourBase=0,  # Base level for contours
        fixedLevelCount=[],  # No fixed contour levels
        useNoData=False,  # Ignore NoData value
        noDataValue=0,  # NoData value
        dstLayer=contour_layer,  # Output OGR layer
        idField=-1,  # No ID field
        elevField=0  # Field index for elevation
    )

    # Extract contours as Shapely polygons
    contours = {}
    for feature in contour_layer:
        geom = feature.GetGeometryRef()  # Get OGR geometry
        elevation = feature.GetField("elevation")  # Get elevation value
        shapely_geom = Polygon(loads(geom.ExportToWkt()))  # Convert to Shapely
        assert isinstance(shapely_geom, Polygon)  # Check for Polygon
        if elevation in contours.keys():
            contours[elevation].append(shapely_geom)
        else:
            contours[elevation] = [shapely_geom]
    if union:
        union_contours = {}
        for elev, polys in contours.items():
            if len(polys) > 1:
                polys = [gpd.GeoDataFrame(geometry=polys).union_all()]
            union_contours[elev] = polys[0]
            union_contours = {key: union_contours[key] for key in sorted(union_contours.keys())}
        return union_contours

    return contours
"""

from osgeo import gdal, ogr
from shapely.wkt import loads as wkt_loads
from shapely.geometry import Polygon, MultiPolygon, LineString, MultiLineString
from shapely.ops import unary_union, polygonize, snap
from shapely.validation import explain_validity
import geopandas as gpd  # optional; we won't depend on it for unions now

# Shapely >= 2.0 has make_valid; older versions can fall back to buffer(0)
try:
    from shapely import make_valid  # Shapely 2.x
except Exception:  # pragma: no cover
    make_valid = None


def _fix_to_polygons(geom, min_area=0.0):
    """
    Return a list of valid Polygon parts >= min_area from any polygonal Shapely geometry.
    """
    if geom.is_empty:
        return []
    g = make_valid(geom) if make_valid is not None else geom.buffer(0)
    if g.is_empty:
        return []
    parts = []
    if isinstance(g, Polygon):
        if g.area >= min_area:
            parts.append(Polygon(g.exterior))
    elif isinstance(g, MultiPolygon):
        for p in g.geoms:
            if p.area >= min_area:
                parts.append(Polygon(p.exterior))
    return parts


def get_contours_as_polygons(raster_path, contour_interval=10, union=True, *,
                             min_area=0.0, snap_tol=0.0):
    """
    Extract contours from a raster and return them as polygons built from polygonized contour lines.

    Parameters
    ----------
    raster_path : str
        Path to the input raster file.
    contour_interval : float
        Interval between contour levels.
    union : bool
        If True, return a single (unioned) polygon per elevation. If False, return a list of polygons per elevation.
    min_area : float, optional
        Minimum polygon area to keep (in raster units^2). Default 0 keeps everything.
    snap_tol : float, optional
        Snap tolerance applied before polygonization to help close nearly-closed lines. 0 disables snapping.

    Returns
    -------
    dict[float, Polygon | list[Polygon]]
        Mapping elevation -> polygon(s). Keys are sorted ascending.
    """
    # --- Open raster
    raster_ds = gdal.Open(raster_path)
    if not raster_ds:
        raise FileNotFoundError(f"Could not open raster file: {raster_path}")
    band = raster_ds.GetRasterBand(1)

    # --- Create in-memory OGR layer for contour lines (vector driver "Memory")
    # Note: The GDAL *raster* memory driver is "MEM"; the OGR *vector* memory driver is "Memory".
    vdriver = ogr.GetDriverByName("Memory")
    contour_ds = vdriver.CreateDataSource("")  # empty name is fine for in-mem
    contour_layer = contour_ds.CreateLayer("contours", geom_type=ogr.wkbLineString)

    # Add elevation field
    field_defn = ogr.FieldDefn("elevation", ogr.OFTReal)
    contour_layer.CreateField(field_defn)

    # --- Generate contours as lines
    # If you have a NoData value on the band, replace noDataValue accordingly and set useNoData=True
    gdal.ContourGenerate(
        band,
        contourInterval=float(contour_interval),
        contourBase=0.0,
        fixedLevelCount=[],    # no fixed levels
        useNoData=False,       # set True + your noDataValue if needed
        noDataValue=0.0,
        dstLayer=contour_layer,
        idField=-1,
        elevField=0
    )

    # --- Group lines by elevation
    by_elev = {}
    for feat in contour_layer:
        elev = float(feat.GetField("elevation"))
        wkt = feat.GetGeometryRef().ExportToWkt()
        g = wkt_loads(wkt)  # Shapely geometry: usually LineString / MultiLineString
        by_elev.setdefault(elev, []).append(g)

    # --- Build polygons per elevation by polygonizing the merged linework
    out = {}
    for elev, geoms in by_elev.items():
        # Flatten/normalize to a list of LineStrings
        lines = []
        for g in geoms:
            if g.is_empty:
                continue
            if isinstance(g, LineString):
                if len(g.coords) >= 4:  # skip degenerate 2-3 point segments
                    lines.append(g)
            elif isinstance(g, MultiLineString):
                for ls in g.geoms:
                    if len(ls.coords) >= 4:
                        lines.append(ls)
            else:
                # If a polygon slipped in (rare), convert its boundary to lines
                if hasattr(g, "boundary"):
                    b = g.boundary
                    if isinstance(b, LineString):
                        lines.append(b)
                    elif isinstance(b, MultiLineString):
                        lines.extend([ls for ls in b.geoms if len(ls.coords) >= 4])

        if not lines:
            # Nothing usable at this elevation
            continue

        merged = unary_union(lines)
        if snap_tol and not merged.is_empty:
            # snapping can help close "almost closed" rings due to floating precision
            merged = snap(merged, merged, snap_tol)

        # polygonize returns an iterable of Polygons constructed from the line network
        poly_iter = polygonize(merged)
        polys = []
        for p in poly_iter:
            polys.extend(_fix_to_polygons(p, min_area=min_area))

        if not polys:
            # Helpful debugging: report why the first geometry was invalid, if any
            reason = ""
            try:
                reason = explain_validity(geoms[0])
            except Exception:
                pass
            # Don’t hard-crash; just skip this elevation if nothing polygonized
            # You can switch to raising an error if this must succeed:
            # raise RuntimeError(f"No valid polygons at elev {elev}: {reason}")
            continue

        if union:
            out[elev] = unary_union(polys)
        else:
            out[elev] = polys

    # Sort by elevation
    out = dict(sorted(out.items(), key=lambda kv: kv[0]))

    # Clean up GDAL/OGR handles
    contour_layer = None
    contour_ds = None
    raster_ds = None

    return out
