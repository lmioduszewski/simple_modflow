from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import shapely as shp
from flopy.utils.geospatial_utils import GeoSpatialUtil
from flopy.utils.triangle import Triangle
from shapely.geometry import Polygon

from simple_modflow.modflow.mf6.grid.helpers import densify_poly
from simple_modflow.modflow.utils.datatypes.readers import read_shp_gpkg


class TriangleGrid(Triangle):

    def __init__(self, angle=32, *args, **kwargs):
        super().__init__(angle=angle, *args, **kwargs)

    def add_circle(
        self,
        radius: float = 100.0,
        center_coords: tuple = (0, 0),
        polygon_to_add=None,
        point_to_add=None,
        point_region_size_max=1,
        return_only: bool = False,
        radians_step: int = 0.1,
    ):
        theta = np.arange(0.0, 2 * np.pi, radians_step)
        x = radius * np.cos(theta) + center_coords[0]
        y = radius * np.sin(theta) + center_coords[1]
        circle_poly = [(x_coord, y_coord) for x_coord, y_coord in zip(x, y)]
        if not return_only:
            self.add_polygon(circle_poly)

        if polygon_to_add:
            polygon_to_add = Polygon(shell=polygon_to_add)
            self.add_polygon(polygon=polygon_to_add)
        if point_to_add:
            self.add_region(point=point_to_add, maximum_area=point_region_size_max)

        return shp.Polygon(circle_poly)

    def add_rectangle(
        self,
        x_dist=100,
        y_dist=100,
        origin=(0, 0),
        return_only=False,
        max_area=None,
    ):
        x_min, y_min = origin[0], origin[1]
        x_max, y_max = x_min + x_dist, y_min + y_dist
        polygon_coords = ((x_min, y_min), (x_min, y_max), (x_max, y_max), (x_max, y_min))
        polygon = shp.Polygon(polygon_coords)
        if not return_only:
            self.add_polygon(polygon)
            if max_area:
                representative_point = polygon.representative_point().coords[0]
                self.add_region(representative_point, maximum_area=max_area)

        return polygon

    def add_regions(self, points, attributes=None, maximum_areas=None):
        if attributes is None:
            attributes = [0 for _ in points]
        if maximum_areas is None:
            maximum_areas = [None for _ in points]

        regions = zip(points, attributes, maximum_areas)
        for region in regions:
            self._regions.append(region)
        return

    def add_poly_regions(
        self,
        shp_gpkg: list | Path,
        points: list | tuple = None,
        use_representative_point: bool = True,
        maximum_areas: list | int | float = None,
        *args,
        **kwargs,
    ):
        if isinstance(shp_gpkg, Path):
            shp_gpkg = [shp_gpkg]
        polys = read_shp_gpkg(shp_gpkg).geometry
        self.add_region(*args, **kwargs)
        return NotImplementedError

    def add_polygon(
        self,
        polygon: shp.Polygon | Path,
        domain: shp.Polygon | Path = None,
        buffer: int | float = 0,
        simplify_tolerance=None,
        ignore_holes=True,
        max_area=None,
        densify_dist: int = None,
    ):
        if isinstance(polygon, Path):
            polygon = read_shp_gpkg(polygon).union_all()
        if domain:
            if isinstance(domain, Path):
                domain = read_shp_gpkg(domain).union_all()
            assert isinstance(domain, shp.Polygon), 'domain must be a shapely polygon'
            if not domain.contains(polygon):
                print('clipping to domain!')
                if isinstance(polygon, shp.Polygon):
                    polygon = gpd.GeoDataFrame(geometry=[polygon])
                polygon = polygon.clip(domain).union_all()
        if buffer != 0:
            polygon = polygon.buffer(buffer)
        if domain:
            assert domain.contains(polygon), 'polygon not fully within domain, try adding a negative buffer'
        if simplify_tolerance:
            polygon = polygon.simplify(simplify_tolerance)
        if densify_dist:
            polygon = densify_poly(polygon, densify_dist)

        super().add_polygon(polygon, ignore_holes=ignore_holes)

        if max_area:
            point = (polygon.representative_point().x, polygon.representative_point().y)
            self.add_region(point, maximum_area=max_area)

    def add_points(self, points: shp.MultiPoint):
        if not isinstance(points, shp.MultiPoint):
            try:
                points = shp.MultiPoint(points)
            except TypeError:
                print('points must be of type shp.MultiPoint')

        geom = GeoSpatialUtil(points).points
        self._polygons.append(geom)

    def add_line_buffer(
        self,
        line: Path,
        buffer: int = 10,
        simplify_tolerance: int = 10,
        densify_dist: int = None,
        domain: shp.Polygon | Path = None,
        max_area: int = None,
        negative_buffer_after_clipping: float | int = 0,
    ):
        if negative_buffer_after_clipping:
            assert negative_buffer_after_clipping <= 0, (
                f'clipping buffer must be less than or equal to zero, not {negative_buffer_after_clipping}'
            )
        line_geom = read_shp_gpkg(line).union_all()
        line_buffer = line_geom.buffer(buffer)
        if simplify_tolerance:
            line_buffer = line_buffer.simplify(simplify_tolerance)
        if densify_dist:
            line_buffer = densify_poly(line_buffer, densify_dist)
        self.add_polygon(
            line_buffer,
            domain=domain,
            max_area=max_area,
            buffer=negative_buffer_after_clipping,
        )

    @staticmethod
    def generate_dissipating_point_cloud(
        polygon: shp.Polygon | Path = None,
        buffer_dist: int | float = 1000,
        num_buffers: int = 5,
        min_area: int | float = None,
        max_area: int | float = None,
        min_spacing: int | float = 200,
        max_spacing: int | float = 1000,
        method: str = 'power',
        exponent: int | float = 2,
    ):
        if isinstance(polygon, Path):
            polygon = read_shp_gpkg(polygon).union_all()
        assert isinstance(polygon, shp.Polygon), 'provided argument cannot be converted to a polygon'

        min_spacing = np.sqrt(min_area) if min_spacing is None else min_spacing
        max_spacing = np.sqrt(max_area) if max_spacing is None else max_spacing

        if method == 'exponential':
            linear_space = np.linspace(0, 1, num_buffers)
            curve = np.exp(linear_space) - 1
            curve /= curve[-1]
        elif method == 'power':
            linear_space = np.linspace(0, 1, num_buffers)
            curve = linear_space ** exponent
        else:
            raise ValueError("Method must be 'exponential' or 'power'.")

        buffer_distances = min_spacing + curve * (buffer_dist - min_spacing)
        spacings = np.linspace(min_spacing, max_spacing, num_buffers)

        pols = [polygon]
        for dist, spac in zip(buffer_distances, spacings):
            pnts = []
            for interp_len in np.arange(0, polygon.buffer(dist).exterior.length, spac):
                pnts.append(polygon.buffer(dist).exterior.line_interpolate_point(interp_len))
            pols.append(shp.Polygon(pnts))

        return pols

    def add_dissipating_buffer_zones(
        self,
        max_area: int | float = 100,
        **kwargs,
    ):
        polys = self.generate_dissipating_point_cloud(**kwargs)
        for i, geom in enumerate(polys):
            if i == 0:
                self.add_polygon(geom, max_area=max_area)
            else:
                self.add_polygon(geom)

    def add_polygons_with_multiregions(
        self,
        shp_gpkg,
        buff_num=5,
        buff_sep=200,
        min_area=50,
        max_area=5_000,
    ):
        poly = read_shp_gpkg(shp_gpkg).union_all()
        rep = poly.representative_point().xy
        rep = [rep[0][0], rep[1][0]]
        y_center = poly.centroid.xy[1][0]
        right_x_edge = poly.bounds[2]

        buff_start = right_x_edge + (buff_sep / 2)
        buff_end = buff_start + (buff_sep * (buff_num - 1))

        x_buffs = np.linspace(buff_start, buff_end, buff_num).tolist()
        y_buffs = [y_center for _ in range(buff_num)]
        buff_points = [rep] + list(zip(x_buffs, y_buffs))
        buff_points = [shp.Point(point) for point in buff_points]

        pols = [poly] + [poly.buffer(buff_sep * mult) for mult in range(1, buff_num + 1)]
        buff_areas = np.linspace(min_area, max_area, buff_num + 1)

        assert len(buff_areas) == len(buff_points) == len(pols)
        for i, pol in enumerate(pols):
            self.add_polygon(pol, simplify_tolerance=10, densify_dist=50)
            self.add_region(buff_points[i], maximum_area=buff_areas[i])
