import subprocess
import sys
import os
from grass_session import Session
import grass_session
import grass.pygrass as pyg
from pathlib import Path
from grass.script import core as gcore
import grass.pygrass.modules.shortcuts as shortCuts
from grass.pygrass.modules.shortcuts import raster as r
from grass.pygrass.modules.shortcuts import general as g
from grass.pygrass.modules.shortcuts import vector as v
import grass.script as gs
import grass.script.setup as gsetup

location = "vector"
mapset = "PERMANENT"
QvtRasterOutput = Path().cwd().joinpath("raster", "QvtRasterfromVectorOutput.tif")
QvtRasterOutputPosix = Path().cwd().joinpath("raster", "QvtRasterfromVectorOutput.tif").as_posix()
gisdb = Path.home().joinpath("grassdata")
mpath = Path.home().joinpath("grassdata", location, mapset)
lpath = Path.home().joinpath("grassdata", location)
grassdata = gisdb
location = location
#mapset = mpath

epsg_code = '2927'
grass8bin = Path("C:/").joinpath('OSGeo4W', 'bin', 'grass83.bat')

qvtshp = Path().cwd().joinpath('shp', 'Qvt contours', 'Qvt Contours Geophys v3.shp').as_posix()
qvtoutput = 'qvtRast_V3.tif'

cdf_lidar = Path().cwd().joinpath('raster','cdf_lidar.tif').as_posix()

def create_grass_location(
    grass8bin=grass8bin,
    epsg_code=epsg_code,
    location_path=lpath
    ):

    grass_cmd = [grass8bin, '-c', 'epsg:' + epsg_code, '-e', lpath]
    #grass_cmd = [grass8bin, '-c', cdf_lidar, '-e', lpath]

    p = subprocess.Popen(grass_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    if p.returncode == 0:
        print("Location successfully created.")
    else:
        print("Error creating location.")
        print(err)

try: #start the grass session
    session = gsetup.init(
        grassdata, 
        location, 
        #mapset
            )
except: #if location doesn't exist, create that first then start session
    create_grass_location()
    session = gsetup.init(
        grassdata, 
        location, 
        #mapset
            )
    
def run_rSurfContour():

    r.surf_contour(
        input='qvtRastOut',
        output='interpedRastQvt',
        overwrite=True,
        verbose=True
    )
    r.out_gdal(
        input='interpedRastQvt',
        output=qvtoutput,
        format='GTiff',
        overwrite=True,
        verbose=True
    )


g.gisenv()
v.in_ogr(
    input=qvtshp,
    output='qvtcontV',
    overwrite=True,
    flags='o'
)
r.in_gdal(
    overwrite=True,
    input=cdf_lidar,
    output='cdflidar',
    flags='o'
)

g.region(
    raster='cdflidar',
    res=4,
    flags='p'
)
v.info(
    map='qvtcontV',
    flags='c'
)
v.to_rast(
    overwrite=True,
    input='qvtcontV',
    type='line',
    output='qvtRastOut',
    use='attr',
    attribute_column='Elev'
)

run_rSurfContour()

session.finish()