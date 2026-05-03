"""Compatibility exports for the legacy Voronoi module."""

from __future__ import annotations

from simple_modflow.modflow.mf6.grid.triangle import TriangleGrid
from simple_modflow.modflow.mf6.grid.voronoi import VoronoiGridPlus

__all__ = [
    "TriangleGrid",
    "VoronoiGridPlus",
]
