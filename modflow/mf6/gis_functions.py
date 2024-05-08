from osgeo import gdal, ogr, osr
import numpy as np
from pathlib import Path
from scipy.interpolate import griddata
import geopandas as gpd

def get_intersecting_polygons(point_shapefile, polygon_shapefile):
    """find intersection of a point shapefile with respect to a 
    polygon shapefile. Return the intersection

    Args:
        point_shapefile (Path): Path of point file
        polygon_shapefile (Path): Path of polygon file

    Returns:
        Geopandas: geopandas dataframe of intersection
    """
    
    # Read the shapefiles into GeoDataFrames
    points = gpd.read_file(point_shapefile)
    polygons = gpd.read_file(polygon_shapefile)
    
    # Ensure that the GeoDataFrames have the same Coordinate Reference System (CRS)
    if points.crs != polygons.crs:
        points = points.to_crs(polygons.crs)
    
    # Use sjoin to find which points lie inside which polygons
    intersections = gpd.sjoin(points, polygons, how="inner", op="within")

    # Return a list of tuples containing point IDs and the corresponding polygon IDs they intersect
    # Assuming 'id' is the column containing the ID for both points and polygons
    return intersections.sort_index()

def contour_to_raster(
    shapefile_path: Path = None, 
    raster_resolution: int = None, 
    output_raster_path: Path = None, 
    epsg = 2927 # Default to NAD83 Harn Washington South
    ) -> gdal.RasterizeLayer:
    """Rasterizes a vector contour shapefile to a given resolution, epsg code,
    and output path

    Args:
        shapefile_path (Path): Path to shapefile
        raster_resolution (int): output raster resolution
        output_raster_path (Path): Path for output file
        epsg (int, optional): Integer of EPSG code. Defaults to 2927#DefaulttoNAD83HarnWashingtonSouth.

    Returns:
        raster: returns the rasterized contours
    """
    
    # Open the shapefile
    shp = ogr.Open(shapefile_path.as_posix())
    layer = shp.GetLayer()

    # Get the extent of the shapefile
    x_min, x_max, y_min, y_max = layer.GetExtent()

    # Create the raster dataset
    x_res = int((x_max - x_min) / raster_resolution)
    y_res = int((y_max - y_min) / raster_resolution)

    raster_ds = gdal.GetDriverByName("GTiff").Create(
        output_raster_path.as_posix(), 
        x_res, 
        y_res, 
        1, 
        gdal.GDT_Float32
    )
    raster_ds.SetGeoTransform(
        (x_min, 
         raster_resolution, 
         0, 
         y_max, 
         0, 
         -raster_resolution)
    )
    raster_band = raster_ds.GetRasterBand(1)
    raster_band.SetNoDataValue(-9999)

    # Setup the spatial reference
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)  
    raster_ds.SetProjection(srs.ExportToWkt())

    # Rasterize the shapefile layer to our new dataset
    return gdal.RasterizeLayer(
        raster_ds, 
        [1], 
        layer, 
        options = ["ATTRIBUTE=Elev"]
        )
    
def interpolate_raster_griddata(
        input_raster_path: Path, 
        output_raster_path: Path, 
        epsg_code: int
        ):
    """Quick and dirty interpolation of an elevation raster from
    a rasterized contour file. Uses griddata from the python module
    scipy

    Args:
        input_raster_path (Path): input raster path
        output_raster_path (Path): output raster path
        epsg_code (int): epsg projection code
    """    
    
    # Open the raster file
    ds = gdal.Open(input_raster_path.as_posix())
    band = ds.GetRasterBand(1)

    # Read the data into numpy arrays
    data = band.ReadAsArray()
    gt = ds.GetGeoTransform()

    # Get the coordinates for each cell
    x = np.arange(gt[0], gt[0] + gt[1] * ds.RasterXSize, gt[1])
    y = np.arange(gt[3], gt[3] + gt[5] * ds.RasterYSize, gt[5])

    # Create a grid for the coordinates
    X, Y = np.meshgrid(x, y)

    # Mask out the no data values
    valid_data = data != band.GetNoDataValue()
    coords = np.column_stack((X[valid_data], Y[valid_data]))
    values = data[valid_data]

    # Perform the interpolation
    interpolated_data = griddata(coords, values, (X, Y), method="cubic")

    # Save the interpolated data to a new raster
    output_ds = gdal.GetDriverByName("GTiff").Create(
        output_raster_path.as_posix(), 
        ds.RasterXSize, 
        ds.RasterYSize, 
        1, 
        gdal.GDT_Float32
    )
    output_ds.SetGeoTransform(gt)

    # Set the coordinate reference system to the provided EPSG code
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg_code)
    output_ds.SetProjection(srs.ExportToWkt())

    return output_ds.GetRasterBand(1).WriteArray(interpolated_data)

    """# Close datasets
    output_ds = None
    ds = None"""
    """# Close datasets
    raster_ds = None
    shp = None"""
    
class Contours:
    """Generic class for a vector contour object based on a input
    shapefile Path.
    """
    def __init__(
        self,
        path: Path
        ):
        
        self.shp_name = path.name[:-4]
        self.path = path #path to contour shapefile
        
    def get_rasterized_contours(
        self,
        out_path: Path = None,
        resolution: int = 4,
        epsg: int = 2927 # Default to NAD83 Harn Washington South
        ):
        """Contours class-specific method to save a rasterized version of the
        Contour vector object

        Args:
            out_path (Path): Path for output raster. Defaults to None.
            resolution (int): Resolution of output raster. Defaults to 4.
            epsg (int): Projection of data as an EPSG code. 
            Defaults to 2927 # NAD83 Harn Washington South.
        """
        
        contour_to_raster(
            self.path,
            raster_resolution = resolution,
            output_raster_path = out_path,
            epsg = epsg)
        
    def get_interpolated_raster_griddata(
        self,
        rasterized_contours_path: Path = None,
        resolution: int = 4,
        epsg: int = 2927,
        interpolated_contours_path: Path = None
        ):
        """Contour class-specific method to generate a rasterized version of the
        vector contours and then save another raster of interpolated data based on 
        the rasterized contours, uses the griddata method from scipy.

        Args:
            rasterized_contours_path (Path, optional): Path to save rasterized contours. Default name based on the shapefile name.
            resolution (int, optional): resolution of rasters. Defaults to 4.
            epsg (int, optional): EPSG Code. Defaults to 2927.
            interpolated_contours_path (Path, optional): Path to save interpolated raster. Default name based on the shapefile nam.
        """
        
        if rasterized_contours_path == None:
            rasterized_contours_path = Path().cwd().joinpath(f'{self.shp_name}.tif')
        if interpolated_contours_path == None:
            interpolated_contours_path = Path().cwd().joinpath(f'interpolated {self.shp_name}.tif')
        
        self.get_rasterized_contours(
            out_path = rasterized_contours_path,
            resolution = resolution,
            epsg = epsg)
        interpolate_raster_griddata(
            input_raster_path = rasterized_contours_path,
            output_raster_path = interpolated_contours_path,
            epsg_code = epsg)