import sys
from pathlib import Path

qgis_path = Path('C:/OSGeo4W/apps/qgis/python/').as_posix()
sys.path.append(qgis_path)

import qgis.core