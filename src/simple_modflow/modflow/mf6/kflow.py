"""Vector-driven hydraulic-conductivity builders for MF6 workflows."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from simple_modflow.modflow.mf6.boundaries import Boundaries

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor
    from simple_modflow.modflow.mf6.simulation.base import SimulationBase


class KFromVector(Boundaries):
    """Map polygon attributes from a vector file onto model-cell K arrays."""

    def __init__(
        self,
        model: "SimulationBase" = None,
        vor: "Vor" = None,
        shp_gpkg: Path = None,
        uid: str = "name",
        crs: int = 2927,
        idomain: list[int] | pd.Series = None,
        idomain_path: Path = None,
    ):
        """Parameters
        ----------
        model
            Model associated with the conductivity data.
        vor
            Voronoi grid helper used to intersect polygons with model cells.
        shp_gpkg
            Polygon or geopackage path describing K zones.
        uid
            Unique-id field in the feature attributes.
        crs
            EPSG code for the feature geometry.
        idomain, idomain_path
            Optional active-domain definition. When provided, inactive cells are
            skipped by default.
        """

        super().__init__(model, vor, shp_gpkg, uid, crs, idomain=idomain, idomain_path=idomain_path)
        self.bound_type = "k"

    def to_frame(
        self,
        fields: dict | None = None,
        nlay: int = 1,
        defaults: list | None = None,
        include_inactive: bool = False,
    ) -> pd.DataFrame:
        """Return conductivity values as a ``(layer, cell)`` indexed dataframe.

        Parameters
        ----------
        fields
            Mapping of logical field names to geometry attribute names. Defaults
            to ``{"name": "name", "k": "k", "layer": "layer"}``.
        nlay
            Number of model layers to include in the returned frame.
        defaults
            Optional per-layer fallback values used to fill NaNs.
        include_inactive
            When ``False``, inactive cells are skipped if an idomain definition
            was provided.
        """

        vor = self._require_vor()
        if fields is None:
            fields = {"name": "name", "k": "k", "layer": "layer"}

        k_index = pd.MultiIndex.from_product(
            [list(range(nlay)), list(range(vor.ncpl))],
            names=["layer", "cell"],
        )
        k_df = pd.DataFrame(index=k_index, columns=["k"], dtype=float)

        gdf = self._indexed_gdf(fields["name"])
        for name, cell_nums in self._candidate_cells_by_name(edges_only=False).items():
            row = gdf.loc[name]
            layer_idx = int(row[fields["layer"]]) - 1
            cells = [int(cell) for cell in cell_nums]
            if not include_inactive:
                cells = [cell for cell in cells if self.inactive_cells is None or cell not in self.inactive_cells]
            for cell in cells:
                k_df.loc[(layer_idx, cell), "k"] = float(row[fields["k"]])

        if defaults is not None:
            if len(defaults) != nlay:
                raise ValueError("defaults must be the same length as nlay")
            for layer_idx, default in enumerate(defaults):
                layer_mask = k_df.index.get_level_values("layer") == layer_idx
                k_df.loc[layer_mask, "k"] = k_df.loc[layer_mask, "k"].fillna(default)

        return k_df

    def from_polygons(
        self,
        *,
        return_array: bool = True,
        fields: dict | None = None,
        nlay: int = 1,
        defaults: list | None = None,
        include_inactive: bool = False,
    ):
        """Build conductivity data from polygon features.

        Parameters
        ----------
        return_array
            When ``True``, return an ``(nlay, ncpl)`` numpy array. When
            ``False``, return a ``(layer, cell)`` indexed dataframe.
        fields, nlay, defaults, include_inactive
            Passed through to :meth:`to_array` / :meth:`to_frame`.
        """

        if return_array:
            return self.to_array(
                fields=fields,
                nlay=nlay,
                defaults=defaults,
                include_inactive=include_inactive,
            )
        return self.to_frame(
            fields=fields,
            nlay=nlay,
            defaults=defaults,
            include_inactive=include_inactive,
        )

    def from_vector(self, **kwargs):
        """Alias for :meth:`from_polygons` for shapefile/geopackage workflows."""

        return self.from_polygons(**kwargs)

    def to_array(
        self,
        fields: dict | None = None,
        nlay: int = 1,
        defaults: list | None = None,
        include_inactive: bool = False,
    ) -> np.ndarray:
        """Return conductivity values as an ``(nlay, ncpl)`` numpy array."""

        k_df = self.to_frame(
            fields=fields,
            nlay=nlay,
            defaults=defaults,
            include_inactive=include_inactive,
        )
        return np.array([k_df.loc[layer_idx, "k"].to_numpy(dtype=float) for layer_idx in range(nlay)])
