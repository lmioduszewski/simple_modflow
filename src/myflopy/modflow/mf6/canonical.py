"""Authoritative contract and grid helpers for the myflopy canonical model."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import Voronoi
from shapely.geometry import Polygon, box

from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus

CANONICAL_TARGET_NAMES = frozenset(
    {"heads", "lake_stage", "sfr_stage", "sfr_flow", "drn_flow"}
)
CANONICAL_PACKAGE_NAMES = frozenset(
    {
        "disv",
        "npf",
        "sto",
        "chd",
        "ghb",
        "drn",
        "rch",
        "wel",
        "lak",
        "sfr",
        "mvr",
        "uzf",
        "oc",
    }
)


@dataclass(frozen=True, slots=True)
class CanonicalModelContract:
    """The invariants every build of the canonical valley model must satisfy.

    The canonical alluvial-valley model is reused as the single fixture across
    examples, integration tests, and PEST notebooks, so its shape must stay stable.
    This frozen dataclass pins those expectations -- layer count and roles,
    icelltype profile, minimum cell/period/SFR-reach counts (with stricter
    ``full_profile`` thresholds), and the required package and target sets -- and
    :meth:`validate` raises ``AssertionError`` listing every violation on a built
    model. The shared instance is :data:`CANONICAL_MODEL_CONTRACT`.

    Attributes are the thresholds/role lists above; see the field defaults for the
    canonical values.
    """

    nlay: int = 4
    minimum_full_ncpl: int = 10_000
    minimum_nper: int = 6
    minimum_sfr_reaches: int = 40
    minimum_full_sfr_reaches: int = 150
    icelltype: tuple[int, ...] = (1, 1, 0, 0)
    layer_roles: tuple[str, ...] = (
        "upper unconfined aquifer",
        "lower unconfined aquifer",
        "aquitard",
        "confined aquifer",
    )
    required_packages: frozenset[str] = CANONICAL_PACKAGE_NAMES
    required_targets: frozenset[str] = CANONICAL_TARGET_NAMES

    def validate(self, model, *, full_profile: bool = False) -> None:
        """Raise ``AssertionError`` when a built model violates the contract."""

        grid = model.gwf.modelgrid
        errors: list[str] = []
        packages = {str(name).lower() for name in model.gwf.package_names}
        targets = set(model.targets.keys())

        if int(grid.nlay) != self.nlay:
            errors.append(f"expected {self.nlay} layers, found {grid.nlay}")
        if full_profile and int(grid.ncpl) < self.minimum_full_ncpl:
            errors.append(
                f"full profile requires >= {self.minimum_full_ncpl:,} cells, found {grid.ncpl:,}"
            )
        if int(model.nper) < self.minimum_nper:
            errors.append(f"expected >= {self.minimum_nper} stress periods, found {model.nper}")
        if missing := self.required_packages - packages:
            errors.append(f"missing packages: {sorted(missing)}")
        if missing := self.required_targets - targets:
            errors.append(f"missing targets: {sorted(missing)}")

        sfr_packagedata = model.gwf.sfr.packagedata.get_data()
        if len(sfr_packagedata) < self.minimum_sfr_reaches:
            errors.append(
                f"SFR requires >= {self.minimum_sfr_reaches} reaches, "
                f"found {len(sfr_packagedata)}"
            )
        if full_profile and len(sfr_packagedata) < self.minimum_full_sfr_reaches:
            errors.append(
                f"full-profile SFR requires >= {self.minimum_full_sfr_reaches} reaches, "
                f"found {len(sfr_packagedata)}"
            )

        icelltype_data = np.asarray(model.gwf.npf.icelltype.get_data())
        if icelltype_data.size == self.nlay:
            icelltype = tuple(int(value) for value in icelltype_data.reshape(-1))
        else:
            icelltype = tuple(
                int(np.asarray(icelltype_data[layer]).reshape(-1)[0])
                for layer in range(self.nlay)
            )
        if icelltype != self.icelltype:
            errors.append(f"expected NPF icelltype={self.icelltype}, found {icelltype}")

        conductivity = np.asarray(model.gwf.npf.k.get_data(), dtype=float)
        if conductivity.shape[0] != self.nlay:
            errors.append(f"expected conductivity for {self.nlay} layers, found shape {conductivity.shape}")
        elif np.nanmedian(conductivity[2]) >= np.nanmedian(conductivity[1]) * 0.05:
            errors.append("layer 3 must be a distinct low-K aquitard")

        well_records = model.gwf.wel.stress_period_data.get_data()
        well_layers = {
            int(record["cellid"][0])
            for records in well_records.values()
            for record in records
        }
        if not ({1, 3} <= well_layers):
            errors.append(
                "WEL must pump the lower unconfined aquifer (layer 2) and confined aquifer (layer 4)"
            )

        areas = np.asarray([geometry.area for geometry in model.vor.gdf_vorPolys.geometry])
        if areas.size > 1 and np.allclose(areas, areas[0]):
            errors.append("DISV cells are uniform rectangles, not an irregular Voronoi mesh")
        for region_name in ("all_lakes", "all_streams"):
            cells = model.get_region_cells(region_name)
            if cells and np.nanmedian(areas[cells]) >= np.nanmedian(areas) * 0.90:
                errors.append(f"{region_name} cells must be measurably refined")
        for region_name in (
            "north_seepage_springs",
            "south_seepage_springs",
            "infiltration_pond",
        ):
            if not model.get_region_cells(region_name):
                errors.append(f"required canonical region {region_name!r} is empty")

        if errors:
            raise AssertionError("Canonical model contract failed:\n- " + "\n- ".join(errors))


CANONICAL_MODEL_CONTRACT = CanonicalModelContract()


def canonical_head_signals(model) -> dict[str, object]:
    """Summarize the head response of a run canonical model into a few scalars.

    Reads the head output of a built+run canonical model and reduces it to the
    per-layer signals the examples/tests assert on: how many time frames were
    written, the spatial head range in each layer at the final step, and the
    temporal change/maximum drawdown per layer (with MF6 dry/inactive sentinels
    masked out). Use it to confirm the model produced the intended hydraulic
    behavior; not a general post-processor.

    Parameters
    ----------
    model
        A built and *run* canonical model (exposing ``gwf.output.head()``).

    Returns
    -------
    dict
        ``frame_count`` plus per-layer ``spatial_range_by_layer``,
        ``temporal_range_by_layer``, and ``maximum_drawdown_by_layer`` arrays.
    """

    reader = model.gwf.output.head()
    frames = np.stack(
        [np.asarray(reader.get_data(kstpkper=key), dtype=float) for key in reader.get_kstpkper()]
    ).reshape(-1, int(model.gwf.modelgrid.nlay), int(model.gwf.modelgrid.ncpl))
    valid = np.isfinite(frames) & (np.abs(frames) < 1.0e29)
    frames = np.where(valid, frames, np.nan)
    temporal_change = frames - frames[0]
    return {
        "frame_count": int(frames.shape[0]),
        "spatial_range_by_layer": np.nanmax(frames[-1], axis=1) - np.nanmin(frames[-1], axis=1),
        "temporal_range_by_layer": np.nanmax(temporal_change, axis=(0, 2))
        - np.nanmin(temporal_change, axis=(0, 2)),
        "maximum_drawdown_by_layer": -np.nanmin(temporal_change, axis=(0, 2)),
    }


def canonical_feature_signals(model) -> dict[str, object]:
    """Summarize the canonical model's head-change and seepage responses after a run.

    Pulls the simulated series at the canonical model's observation features and
    reduces them to two compact signals the examples/tests check: the head-change
    range at each head target (capturing pond/pumping influence) and the peak
    magnitude of each DRN seepage series. Expects the canonical model's target
    set (``model.targets.heads`` / ``model.targets.drn_flow``).

    Parameters
    ----------
    model
        A built and *run* canonical model with the canonical target set attached.

    Returns
    -------
    dict
        ``head_change`` (per head target) and ``seepage_peak`` (per DRN series).
    """

    heads = model.targets.heads.targets.simulated_heads(model)
    head_change = {
        column: float(heads[column].max() - heads[column].min())
        for column in heads.columns
        if column != "per"
    }
    seepage = model.targets.drn_flow.targets.simulated_series(model)
    seepage_peak = {
        column: float(np.nanmax(np.abs(seepage[column].to_numpy(dtype=float))))
        for column in seepage.columns
        if column != "time"
    }
    return {"head_change": head_change, "seepage_peak": seepage_peak}


def canonical_sfr_signals(model, *, per: int = 0) -> dict[str, object]:
    """Summarize the canonical stream's routing, stage, and GW-exchange after a run.

    Walks the SFR long-profile and budget for one stress period and reduces them
    to the stream signals the examples/tests assert on: reach counts (total/wet/
    dry), wetted-depth extremes, the split of losing vs gaining reaches (using the
    MF6 sign convention where stream->aquifer flow is positive), and the
    routed-flow range. Confirms the canonical valley produces the intended
    gaining/losing stream behavior.

    Parameters
    ----------
    model
        A built and *run* canonical model with an SFR package.
    per
        Stress period to summarize (default 0).

    Returns
    -------
    dict
        Reach counts, min/max wetted depth, losing/gaining reach counts, and
        min/max routed flow for the period.
    """

    profile = model.packages.sfr.results.long_profile(per=per)
    stage = profile["stage"].to_numpy(dtype=float)
    streambed_top = profile["streambed_top"].to_numpy(dtype=float)
    exchange = profile["q"].to_numpy(dtype=float)
    wet = np.isfinite(stage) & (stage > -1.0e20)
    depth = np.where(wet, stage - streambed_top, np.nan)

    routing = model.packages.sfr.budget.flow_ja_face.get(per=per)
    routed_flow = np.asarray([], dtype=float)
    if not routing.empty:
        last_kstpkper = max(routing["kstpkper"], key=lambda key: int(key[0]))
        routing = routing.loc[
            routing["kstpkper"].apply(lambda key: tuple(key) == tuple(last_kstpkper))
        ]
        routed_flow = routing.loc[routing["q"] < 0.0, "q"].abs().to_numpy(dtype=float)

    tolerance = 1.0e-6
    return {
        "reach_count": int(len(profile)),
        "wet_reach_count": int(wet.sum()),
        "dry_reach_count": int((~wet).sum()),
        "minimum_depth": float(np.nanmin(depth)) if wet.any() else np.nan,
        "maximum_depth": float(np.nanmax(depth)) if wet.any() else np.nan,
        # MF6 SFR GWF flow is positive from the stream to groundwater.
        "losing_reach_count": int((exchange > tolerance).sum()),
        "gaining_reach_count": int((exchange < -tolerance).sum()),
        "minimum_routed_flow": float(np.nanmin(routed_flow)) if routed_flow.size else np.nan,
        "maximum_routed_flow": float(np.nanmax(routed_flow)) if routed_flow.size else np.nan,
    }


def canonical_partition_mask(model, nparts: int) -> np.ndarray:
    """Build a domain-decomposition mask for the canonical model that keeps lakes whole.

    A purpose-built partitioner for splitting the canonical valley model for
    parallel runs: it lays out ``nparts`` balanced vertical bands across the grid,
    then adjusts them so each physical lake stays within a single subdomain and
    every partition remains spatially contiguous. Returns the per-cell subdomain
    assignment expected by the parallel split workflow.

    Parameters
    ----------
    model
        The canonical model to partition.
    nparts
        Number of subdomains (must be >= 2).

    Returns
    -------
    numpy.ndarray
        Per-cell integer subdomain ids.

    Raises
    ------
    ValueError
        If ``nparts`` is less than 2.
    """

    if nparts < 2:
        raise ValueError("nparts must be at least 2")
    from myflopy.modflow.mf6.parallel import (
        _preserve_lake_partitions,
        _repair_partition_contiguity,
    )

    x = np.asarray(model.gwf.modelgrid.xcellcenters, dtype=float).reshape(-1)
    order = np.argsort(x, kind="stable")
    mask = np.empty(x.size, dtype=int)
    for partition, columns in enumerate(np.array_split(order, nparts)):
        mask[columns] = partition
    mask = _preserve_lake_partitions(model.gwf, mask)
    return _repair_partition_contiguity(model.gwf, mask)


def irregular_voronoi_grid(
    *,
    nrow: int,
    ncol: int,
    cell_size: float,
    seed: int = 2026,
    jitter_fraction: float = 0.32,
    refined: bool = True,
) -> VoronoiGridPlus:
    """Return a deterministic bounded irregular Voronoi DISV helper.

    Seeds retain row-major ordering so stable feature and boundary selections can
    be shared by the full and validation profiles.
    """

    rng = np.random.default_rng(seed)
    x_fraction = (np.arange(ncol, dtype=float) + 0.5) / ncol
    y_fraction = (np.arange(nrow, dtype=float) + 0.5) / nrow
    if refined:
        # Allocate more seeds to the mountain front, surface-water corridor,
        # pumping centers, and downgradient capture boundary.
        x_fraction = np.interp(
            x_fraction,
            [0.0, 0.20, 0.45, 0.90, 1.0],
            [0.0, 0.20, 0.60, 0.90, 1.0],
        )
        y_fraction = np.interp(
            y_fraction,
            [0.0, 0.15, 0.45, 0.85, 1.0],
            [0.0, 0.15, 0.45, 0.85, 1.0],
        )
    x = x_fraction * ncol * cell_size
    y = y_fraction * nrow * cell_size
    points = np.asarray([(xv, yv) for yv in y for xv in x], dtype=float)
    if refined:
        width = ncol * cell_size
        height = nrow * cell_size
        normalized_x = points[:, 0] / width
        normalized_y = points[:, 1] / height
        # Concentrate seeds along the converging stream network so the SFR
        # corridor and terminal lake resolve on finer cells: two tributaries
        # meeting the valley axis at the (0.40, 0.50) confluence, then a main
        # stem running east to the lake.
        trib_fraction = np.clip((normalized_x - 0.03) / 0.37, 0.0, 1.0)
        north_y = 0.85 - 0.35 * trib_fraction
        south_y = 0.15 + 0.35 * trib_fraction
        main_y = np.full_like(normalized_x, 0.5)
        far = 10.0
        distances = np.vstack(
            [
                np.where(normalized_x <= 0.41, np.abs(normalized_y - north_y), far),
                np.where(normalized_x <= 0.41, np.abs(normalized_y - south_y), far),
                np.where(normalized_x >= 0.39, np.abs(normalized_y - main_y), far),
            ]
        )
        candidates = np.vstack([north_y, south_y, main_y])
        choice = np.argmin(distances, axis=0)
        target_y = candidates[choice, np.arange(len(normalized_x))]
        near_stream = np.min(distances, axis=0) < 0.15
        points[near_stream, 1] += (
            target_y[near_stream] - normalized_y[near_stream]
        ) * height * 0.55
    points += rng.uniform(-jitter_fraction, jitter_fraction, size=points.shape) * np.asarray(
        [np.median(np.diff(x)), np.median(np.diff(y))]
    )

    width = ncol * cell_size
    height = nrow * cell_size
    domain = box(0.0, 0.0, width, height)
    regions, vertices = _finite_voronoi_regions(Voronoi(points), radius=max(width, height) * 4)

    vertex_ids: dict[tuple[float, float], int] = {}
    verts: list[tuple[float, float]] = []
    iverts: list[list[int]] = []
    centers: list[tuple[float, float]] = []
    for point, region in zip(points, regions):
        polygon = Polygon(vertices[region]).intersection(domain)
        if polygon.geom_type != "Polygon":
            polygon = max(polygon.geoms, key=lambda geometry: geometry.area)
        coords = list(polygon.exterior.coords)[:-1]
        cell_vertices: list[int] = []
        for coord in coords:
            key = (round(float(coord[0]), 8), round(float(coord[1]), 8))
            if key not in vertex_ids:
                vertex_ids[key] = len(verts)
                verts.append(key)
            cell_vertices.append(vertex_ids[key])
        iverts.append(cell_vertices)
        centers.append((float(point[0]), float(point[1])))

    return VoronoiGridPlus(
        verts=np.asarray(verts, dtype=float),
        iverts=iverts,
        xcyc=np.asarray(centers, dtype=float),
    )


def _finite_voronoi_regions(vor: Voronoi, radius: float) -> tuple[list[list[int]], np.ndarray]:
    """Reconstruct infinite 2-D SciPy Voronoi regions as finite polygons."""

    center = vor.points.mean(axis=0)
    new_regions: list[list[int]] = []
    new_vertices = vor.vertices.tolist()
    ridges: dict[int, list[tuple[int, int, int]]] = {}
    for (point_a, point_b), (vertex_a, vertex_b) in zip(vor.ridge_points, vor.ridge_vertices):
        ridges.setdefault(point_a, []).append((point_b, vertex_a, vertex_b))
        ridges.setdefault(point_b, []).append((point_a, vertex_a, vertex_b))

    for point_id, region_id in enumerate(vor.point_region):
        region = vor.regions[region_id]
        if region and all(vertex >= 0 for vertex in region):
            new_regions.append(region)
            continue

        rebuilt = [vertex for vertex in region if vertex >= 0]
        for neighbor, vertex_a, vertex_b in ridges[point_id]:
            if vertex_b < 0:
                vertex_a, vertex_b = vertex_b, vertex_a
            if vertex_a >= 0:
                continue
            tangent = vor.points[neighbor] - vor.points[point_id]
            tangent /= np.linalg.norm(tangent)
            normal = np.asarray([-tangent[1], tangent[0]])
            midpoint = vor.points[[point_id, neighbor]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, normal)) * normal
            new_vertices.append((vor.vertices[vertex_b] + direction * radius).tolist())
            rebuilt.append(len(new_vertices) - 1)

        polygon = np.asarray([new_vertices[vertex] for vertex in rebuilt])
        centroid = polygon.mean(axis=0)
        angles = np.arctan2(polygon[:, 1] - centroid[1], polygon[:, 0] - centroid[0])
        new_regions.append(np.asarray(rebuilt)[np.argsort(angles)].tolist())
    return new_regions, np.asarray(new_vertices)
