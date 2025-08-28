from osgeo import gdal, ogr
from shapely.wkt import loads
from shapely.geometry import Polygon
import geopandas as gpd


def get_contours_as_polygons(raster_path, contour_interval=10, union=True):
    """
    Extract contours from a raster and return them as Shapely polygons.

    Parameters:
        raster_path (str): Path to the input raster file.
        contour_interval (float): Interval between contour levels.
        union (bool): If True, return union of all contours at each level.

    Returns:
        list[dict]: A list of dictionaries with "elevation" and "geometry" (Shapely polygons).
    """
    # Open the raster file
    raster_ds = gdal.Open(raster_path)
    if not raster_ds:
        raise FileNotFoundError(f"Could not open raster file: {raster_path}")

    # Get the raster band (first band)
    band = raster_ds.GetRasterBand(1)

    # Create an in-memory data source for contours
    driver = ogr.GetDriverByName("Memory")
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
