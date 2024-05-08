import subprocess
import sys
import os
from pathlib import Path
os.environ['GISBASE'] = Path(r'C:\OSGeo4W\apps\grass\grass83').as_posix()

from grass.pygrass.modules.shortcuts import raster as r
from grass.pygrass.modules.shortcuts import general as g
from grass.pygrass.modules.shortcuts import vector as v
import grass.script.setup as gsetup

class SurfaceInterpFromShp:

    def __init__(
            self,
            location: str = 'temp_loc',
            mapset: str = 'MAPSET',
            grassdata:str = 'grassdata',
            epsg: str = '2927',
            shp_path: Path = None,
            region_dimensions_raster:Path = None,
            write_interpolated_surface_and_finish: bool = False,
            output_resolution = 4,
            shp_attribute_for_z = 'Elev',
            data_dir: Path = Path.home() / 'Python',
            surf_out = 'interp_surface.tif'
    ):
        """
        Class used to create an interpolated raster surface from a set of vector contours.
        As in, an elevation surface from elevation contours.
        :param location:
        :param mapset:
        :param grassdata:
        :param epsg:
        :param shp_path:
        :param region_dimensions_raster:
        :param write_interpolated_surface_and_finish:
        :param output_resolution:
        :param shp_attribute_for_z:
        :param data_dir:
        :param surf_out:

        Example Usage:

            top_of_qpf_shp = Path('shapefile path').as_posix()
            region_raster = Path('raster to define region').as_posix()
            interp_qpf = SurfaceInterpFromShp(
                shp_path=top_of_qpf_shp,
                region_dimensions_raster= region_raster,
                shp_attribute_for_z='Elevation',
                output_resolution=15
            )
            interp_qpf.write_surf()
        """
        self.grass8bin = Path(r'C:\OSGeo4W\bin\grass83.bat')
        self.grassdata = Path.home().joinpath(grassdata)
        self.epsg_code = epsg
        self.session = None

        self.location = location
        self.location_path = self.grassdata / self.location

        self.mapset = mapset
        self.mapset_path = self.grassdata.joinpath(self.location, self.mapset)

        #  default names to use for the various working rasters and vectors in the grass session
        self.grassname_vect_cont = 'vectContours'
        self.grassname_rast_cont = 'rastContours'
        # define the full surface output path
        self.surf_out = (data_dir / surf_out).as_posix()
        self.region_dimensions_raster = region_dimensions_raster

        #  if writing an interpolated surface on object init is set to True, run self.write_surf()
        self.shp_path = shp_path
        self.shp_attribute_for_z = shp_attribute_for_z
        self.write_interpolated_surface_and_finish = write_interpolated_surface_and_finish
        self.output_resolution = output_resolution
        if self.shp_path is not None and self.write_interpolated_surface_and_finish:
            self.write_surf()

    def write_surf(self, finish=True):
        """create a grass location, start a grass session, import the shapefile of contours,
        rasterize the contours, generate an interpolated raster surface, then write as a
        geotiff, and finally close the session if 'finish' arg is True"""

        print('starting grass session')
        self.start_grass_session()
        print('importing vector contours')
        self.import_vector_contours(self.shp_path)

        #  then define the region using a raster (such as a clip of lidar or other raster
        #  that on the same projection as the vector contours
        print('setting region from raster')
        self.set_region_from_raster(self.region_dimensions_raster)
        print('rasterizing contours')
        self.rasterize_vector_contours()
        print('interpolating surface')
        self.interpolate_surface_from_rasterized_contours()
        print('done, closing grass session')

        if finish is True:
            self.session.finish()

    def create_grass_location(self):
        """Create location on the computer's hard drive (a file folder) that will
        contain the working files for the grass gis session"""
        grass_cmd = [self.grass8bin, '-c', 'epsg:' + self.epsg_code, '-e', self.location_path]
        p = subprocess.Popen(grass_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        if p.returncode == 0:
            print("Location successfully created.")
        else:
            print("Error creating location.")
            print(err)

    def start_grass_session(self):
        """start the grass gis session based on provided location. The location will
        be created if it does not exist on the hard drive in the grassdata directory."""
        grassdata = self.grassdata
        location = self.location
        try: #start the grass session
            self.session = gsetup.init(grassdata, location,
                #mapset
                    )
        except: #if location doesn't exist, create that first then start session
            self.create_grass_location()
            self.session = gsetup.init(grassdata, location,
                #mapset
                    )

    @property
    def gisenv(self):
        return g.gisenv()

    def import_vector_contours(self, shp_path=None, grass_vector_name=None, info=True):
        """
        imports a shapfile of vector contours to the grass gis session
        :param shp_path: path to the shapefile
        :return: None
        """
        if grass_vector_name is not None:
            self.grassname_vect_cont = grass_vector_name
        if shp_path is not None:
            self.shp_path = shp_path
        v.in_ogr(
            input=self.shp_path,
            output=self.grassname_vect_cont,
            overwrite=True,
            flags='o'
        )
        if info:
            v.info(
                map=self.grassname_vect_cont,
                flags='c'
            )

    def set_region_from_raster(self, raster=None, resolution=None):
        """
        sets the grass gis region from a provided raster. The raster should be on the
        same projection as the contour shapefile.
        :param raster: raster used to define the region
        :param resolution: resolution to use for the raster defined region
        :return: None
        """
        raster = raster if raster is not None else self.region_dimensions_raster
        resolution = resolution if resolution is not None else self.output_resolution
        r.in_gdal(
            overwrite=True,
            input=raster,
            output='regionRaster',
            flags='o'
        )
        g.region(
            raster='regionRaster',
            res=resolution,
            flags='p'
        )

    def rasterize_vector_contours(self, attribute_for_z=None):
        """Create a rasterized version of the shapefile contours provided as an arguement
        to the function. The rasterized contours are then used to interpolate a surface in
        a separate function: interpolate_surface_from_rasterized_contours"""
        attribute_for_z = attribute_for_z if attribute_for_z is not None else self.shp_attribute_for_z
        v.to_rast(
            overwrite=True,
            input=self.grassname_vect_cont,
            type='line',
            output=self.grassname_rast_cont,
            use='attr',
            attribute_column=attribute_for_z
        )
        r.out_gdal(
            input=self.grassname_rast_cont,
            output='rast_contours.tif',
            format='GTiff',
            overwrite=True,
            verbose=True
        )

    def interpolate_surface_from_rasterized_contours(self, surf_out=None):
        """
        interpolate a surface from rasterized contours
        :param surf_out: file (geotiff) name to output
        :return: writes a file to surf_out
        """
        surf_out = surf_out if surf_out is not None else self.surf_out
        r.surf_contour(
            input=self.grassname_rast_cont,
            output='interpdSurface',
            overwrite=True,
            verbose=True
        )
        print('writing to geotiff')
        r.out_gdal(
            input='interpdSurface',
            output=surf_out,
            format='GTiff',
            overwrite=True,
            verbose=True
        )


if __name__ == '__main__':

    top_of_qpf_shp = Path('C:/Users/lukem/Python/MODFLOW/LakePointe/inputs/surfaces/top_of_qpff_bot.shp').as_posix()
    region_raster = Path('C:/Users/lukem/Python/MODFLOW/LakePointe/inputs/lakepointe_lidar.tif').as_posix()
    interp_qpf = SurfaceInterpFromShp(
        shp_path=top_of_qpf_shp,
        region_dimensions_raster= region_raster,
        shp_attribute_for_z='Elevation',
        output_resolution=15,
        surf_out='qpfftop_model_bottom.tif'
    )
    interp_qpf.write_surf()