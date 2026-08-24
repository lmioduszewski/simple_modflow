"""Layer-centric authoring facade for model discretization.

`LayerStack` is a thin, friendly front door over :class:`~myflopy.surfaces.LayerSurfaces`.
You declare a model top, then ``add`` named layers -- each carrying its own bottom
definition and (optionally) its own ``min_thickness`` / ``pinch`` policy::

    from myflopy.layers import LayerStack, Raster, Flat, Contours

    stack = LayerStack(vor, top=Raster("ground.tif"), length_units="feet")
    stack.add("alluvium", bottom=Raster("allu_bot.tif"))
    stack.add("clay",     thickness=30, min_thickness=1, pinch="passthrough")
    stack.add("bedrock",  bottom=Contours("bedrock.gpkg", z="elev"), fill="propagate")

    result = stack.build()          # result.top / .botm / .idomain / .report()
    disv   = stack.to_disv(vor)     # ready-to-use mf.disv spec

This module adds **no new geometry logic**: it translates "top + N named layers"
into the "N+1 surfaces" list and delegates sampling, reconcile, and pinch-out to
the existing, tested `LayerSurfaces` engine.
"""

from __future__ import annotations

from dataclasses import dataclass
from dataclasses import replace as _dc_replace
from pathlib import Path
from typing import Any

import numpy as np

from myflopy._optional import require
from myflopy.modflow.mf6.grid.plotting import GridPlots
from myflopy.surfaces import LayerSurfaces, Surface
from myflopy.viz import Fig, MplPicture, Picture, VtkScene, mpl_axes

# Convenience source aliases -- the existing Surface constructors under friendlier
# names for layer authoring (NOT new implementations).
Raster = Surface.raster
Flat = Surface.flat
Contours = Surface.from_contours
Points = Surface.from_points
Array = Surface.from_array
Isopach = Surface.isopach
# Surface algebra (lower/upper envelopes, clamping, zone selection).
Min = Surface.minimum
Max = Surface.maximum
Clamp = Surface.clamp
Where = Surface.where

_UNSET = object()


def _coerce_surface(x) -> Surface:
    """Accept a Surface, or coerce a path -> Raster and a number -> Flat."""

    if isinstance(x, Surface):
        return x
    if isinstance(x, (str, Path)):
        return Surface.raster(x)
    if isinstance(x, (int, float)):
        return Surface.flat(float(x))
    raise TypeError(f"Expected a Surface, path, or number; got {type(x).__name__}.")


def _make_surface(bottom, thickness, fill) -> Surface:
    """Build a layer's bottom surface from ``bottom=`` or ``thickness=``."""

    if (bottom is None) == (thickness is None):
        raise ValueError("Provide exactly one of bottom= or thickness=.")
    if thickness is not None:
        surface = (
            thickness
            if isinstance(thickness, Surface)
            else Surface.constant_thickness(float(thickness))
        )
    else:
        surface = _coerce_surface(bottom)
    if fill is not None:
        surface = _dc_replace(surface, fill=fill)
    return surface


def _reconcile_args(reconcile) -> tuple[bool, str]:
    """Map the facade ``reconcile`` argument to ``LayerSurfaces.sample`` kwargs."""

    if reconcile in (False, None):
        return False, "bottom"
    if reconcile is True:
        return True, "bottom"
    if reconcile in ("bottom", "top"):
        return True, reconcile
    raise ValueError("reconcile must be 'bottom', 'top', or False.")


def modflow_surfaces(source, *, resample: bool = True):
    """Read an existing MODFLOW model's surfaces as :class:`Surface` objects.

    ``source`` is a flopy model (anything exposing ``.modelgrid``) or a flopy
    modelgrid directly. Returns ``(top_surface, [bottom_surface, ...])``. With
    ``resample=True`` (default) each surface is interpolated from the source cell
    centres, so it transfers onto a *different* grid; with ``resample=False`` the
    arrays are used verbatim (the target grid must match cell-for-cell).
    """
    mg = getattr(source, "modelgrid", source)
    top = np.asarray(mg.top, dtype=float).ravel()
    botm = np.asarray(mg.botm, dtype=float)
    botm = botm.reshape(botm.shape[0], -1)
    if resample:
        xc = np.asarray(mg.xcellcenters, dtype=float).ravel()
        yc = np.asarray(mg.ycellcenters, dtype=float).ravel()
        top_s = Surface.from_points(xc, yc, top)
        botm_s = [Surface.from_points(xc, yc, botm[k]) for k in range(botm.shape[0])]
    else:
        top_s = Surface.from_array(top)
        botm_s = [Surface.from_array(botm[k]) for k in range(botm.shape[0])]
    return top_s, botm_s


@dataclass
class _Layer:
    name: str
    surface: Surface
    min_thickness: float | None = None
    pinch: str | None = None


@dataclass
class LayerQCReport:
    """Diagnostics for a built layer stack -- the problems to fix before MF6 runs.

    The structured result of QC-ing a :class:`LayerBuildResult` (via
    ``result.qc()`` or :meth:`LayerStack.qc`). It counts the geometry pathologies
    that make a DISV grid fail or behave oddly -- cells with no source coverage
    (NaN top/bottom), active cells bounded by a NaN surface or with non-positive
    thickness (fatal), overly thin cells, pinched-out cells, isolated active cells
    with no connection, and how many disconnected active components exist. When
    produced by ``LayerStack.qc`` it also reports how much top-down reconciliation
    moved each surface. Check :attr:`ok` for a fatal/clean verdict, or ``str(report)``
    for a human-readable per-layer breakdown.

    The integer/list fields hold counts (per layer where noted); see the inline
    field comments. ``isolated_active`` lists ``(layer, cell)`` pairs.
    """

    nlay: int
    ncpl: int
    names: list[str]
    nan_top: int                       # cells with no top elevation (no coverage)
    nan_botm: list[int]                # per layer: cells with no bottom elevation
    nan_active_cells: int              # active cells bounded by a NaN surface (fatal)
    nonpositive_active: list[int]      # per layer: active cells with thickness <= 0
    thin: list[int]                    # per layer: cells thinner than min_thickness
    pinched: list[int]                 # per layer: cells removed (idomain != 1)
    isolated_active: list[tuple]       # (layer, cell) active cells with no connection
    n_active_components: int           # connected components of active cells
    component_sizes: list[int]         # sizes, largest first
    # reconcile diagnostics -- only filled by LayerStack.qc (needs the surfaces)
    reconcile_adjusted: list[int] | None = None    # per layer: cells reconcile moved
    reconcile_max_shift: list[float] | None = None  # per layer: largest move

    @property
    def ok(self) -> bool:
        """True when nothing fatal was found (no NaN-bounded or isolated active cells,
        no non-positive active thickness). Thin/pinched cells are expected, not fatal."""
        return (
            self.nan_active_cells == 0
            and sum(self.nonpositive_active) == 0
            and len(self.isolated_active) == 0
        )

    def __str__(self) -> str:
        """Render a human-readable multi-line QC summary (status, warnings, per-layer stats)."""

        head = "OK" if self.ok else "PROBLEMS FOUND"
        lines = [
            f"LayerStack QC [{head}]: {self.nlay} layers, {self.ncpl} cells, "
            f"{self.n_active_components} active component(s)"
        ]
        if self.nan_active_cells:
            lines.append(f"  ** {self.nan_active_cells} active cell(s) bounded by a NaN "
                         "surface (no source coverage) -- fix the source or fill='propagate'")
        if self.isolated_active:
            lines.append(f"  ** {len(self.isolated_active)} isolated active cell(s) with no "
                         "connection -- prune_isolated() deactivates them")
        if sum(self.nonpositive_active):
            lines.append(f"  ** {sum(self.nonpositive_active)} active cell(s) with thickness <= 0")
        if self.n_active_components > 1:
            shown = ", ".join(str(s) for s in self.component_sizes[:5])
            lines.append(f"  note: active cells split into {self.n_active_components} components "
                         f"(sizes: {shown}{'...' if self.n_active_components > 5 else ''})")
        for i, name in enumerate(self.names):
            extra = ""
            if self.reconcile_adjusted is not None:
                extra = (f"  reconcile_moved={self.reconcile_adjusted[i]} "
                         f"(max {self.reconcile_max_shift[i]:.2f})")
            lines.append(
                f"  [{i}] {name:<14} nan_bottom={self.nan_botm[i]} "
                f"thin={self.thin[i]} pinched={self.pinched[i]} "
                f"thickness<=0(active)={self.nonpositive_active[i]}{extra}"
            )
        return "\n".join(lines)


@dataclass
class LayerBuildResult:
    """The disv-ready arrays produced by :meth:`LayerStack.build`.

    Bundles the discretization arrays a layer stack resolves to -- ``top``
    ``(ncpl,)``, ``botm`` ``(nlay, ncpl)``, ``idomain`` ``(nlay, ncpl)`` (1 active,
    -1 pass-through, 0 inactive), and the derived ``thickness`` -- alongside the
    per-layer metadata (``names``, ``min_thickness``, ``pinch`` policy, units) and a
    back-reference to the ``vor`` grid. Feed ``.top``/``.botm``/``.idomain`` straight
    into ``mf.disv(...)``. Inspect quality with :meth:`report` (a per-layer
    thickness/pinch summary) or ``.qc()`` (a :class:`LayerQCReport`), and publish the
    surfaces onto the grid for the GIS-aware builders with :meth:`attach_to_grid`
    (or ``LayerStack.build(attach=True)``).

    The array shapes and idomain encoding are noted in the inline field comments.
    """

    top: np.ndarray          # (ncpl,)
    botm: np.ndarray         # (nlay, ncpl)
    idomain: np.ndarray      # (nlay, ncpl): 1 active, -1 pass-through, 0 inactive
    thickness: np.ndarray    # (nlay, ncpl)
    names: list[str]
    min_thickness: list[float]
    pinch: list[str]
    length_units: str
    time_units: str
    vor: Any = None          # grid the result was built on (for views)

    @property
    def nlay(self) -> int:
        """Number of layers (rows of ``botm``)."""

        return self.botm.shape[0]

    @property
    def n_pinched(self) -> int:
        """Cells removed from the solution (idomain != 1)."""
        return int((self.idomain != 1).sum())

    def report(self) -> str:
        """Per-layer thickness + thin/pinched-cell summary."""
        lines = [
            f"LayerStack: {self.nlay} layers, {self.top.size} cells "
            f"[{self.length_units}/{self.time_units}], {self.n_pinched} pinched cells"
        ]
        for i, name in enumerate(self.names):
            t = self.thickness[i]
            thin = int((t < float(self.min_thickness[i])).sum())
            pinched = int((self.idomain[i] != 1).sum())
            lines.append(
                f"  [{i}] {name:<14} thk min={t.min():.2f} mean={t.mean():.2f} "
                f"max={t.max():.2f}  thin(<{self.min_thickness[i]})={thin} "
                f"pinched={pinched} ({self.pinch[i]})"
            )
        return "\n".join(lines)

    def attach_to_grid(self, vor=None):
        """Publish ``top``/``botm`` onto ``vor.gdf_topbtm`` for grid-aware builders.

        Writes a centroid GeoDataFrame with integer columns (``0`` = model top,
        ``1..nlay`` = layer bottoms) -- the format the surface-aware builders read
        (SFR reach tops, LAK lake-cell layering). Returns the grid.

        Usually you do not call this directly: pass ``attach=True`` to
        :meth:`LayerStack.build`.
        """
        import geopandas as gpd

        grid = self.vor if vor is None else vor
        if grid is None:
            raise ValueError("No grid to attach to; pass vor= or build the stack with a grid.")
        columns = {0: np.asarray(self.top, dtype=float)}
        for i in range(self.nlay):
            columns[i + 1] = np.asarray(self.botm[i], dtype=float)
        grid.gdf_topbtm = gpd.GeoDataFrame(
            {"geometry": grid.gdf_vorPolys.geometry, **columns},
            geometry="geometry", crs=grid.crs,
        )
        return grid

    # -- QC / validation -------------------------------------------------- #
    def _active_components(self):
        """Union-find over active cells (idomain == 1). Returns ``(labels, sizes)``:
        ``labels`` is an ``(nlay, ncpl)`` array of component roots (-1 where not
        active) and ``sizes`` maps a root to its cell count.

        Connections follow MF6: horizontal between active plan-neighbors, and
        vertical down a column where ``idomain == -1`` (pass-through) bridges
        active cells while ``idomain == 0`` (inactive) blocks them."""
        from collections import Counter

        nlay, ncpl = self.idomain.shape
        parent = list(range(nlay * ncpl))

        def find(a):
            """Union-find root of flat cell index ``a``, with path compression."""

            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        def union(a, b):
            """Merge the union-find sets containing flat cell indices ``a`` and ``b``."""

            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb

        active = self.idomain == 1
        adj = self.vor.adjacent_cells_idx
        for k in range(nlay):
            base = k * ncpl
            for i in range(ncpl):
                if not active[k, i]:
                    continue
                for j in adj[i]:
                    if 0 <= j < ncpl and active[k, j]:
                        union(base + i, k * ncpl + j)
        for i in range(ncpl):
            last = None
            for k in range(nlay):
                d = self.idomain[k, i]
                if d == 1:
                    if last is not None:
                        union(last * ncpl + i, k * ncpl + i)
                    last = k
                elif d == 0:
                    last = None      # inactive blocks vertical flow
                # d == -1: pass-through keeps ``last`` reachable
        labels = np.full((nlay, ncpl), -1, dtype=int)
        sizes = Counter()
        for k in range(nlay):
            for i in range(ncpl):
                if active[k, i]:
                    r = find(k * ncpl + i)
                    labels[k, i] = r
                    sizes[r] += 1
        return labels, sizes

    def qc(self) -> LayerQCReport:
        """Check the built stack for problems MODFLOW 6 would choke on.

        Catches **NaN-bounded active cells** (a surface had no source coverage),
        **non-positive active thickness**, and **isolated active cells** (no
        connection to any neighbour) -- plus thin/pinched counts and the number
        of connected active components. Returns a :class:`LayerQCReport` whose
        ``.ok`` is False if anything fatal was found. ``LayerStack.qc`` adds
        reconcile diagnostics on top of this."""
        nlay, ncpl = self.idomain.shape
        active = self.idomain == 1
        surf_top = np.vstack([self.top[None, :], self.botm[:-1]])  # top of each layer
        nan_bounded = np.isnan(surf_top) | np.isnan(self.botm)
        labels, sizes = self._active_components()
        isolated = [
            (k, i)
            for k in range(nlay)
            for i in range(ncpl)
            if active[k, i] and sizes[labels[k, i]] == 1
        ]
        return LayerQCReport(
            nlay=nlay,
            ncpl=ncpl,
            names=list(self.names),
            nan_top=int(np.isnan(self.top).sum()),
            nan_botm=[int(np.isnan(self.botm[k]).sum()) for k in range(nlay)],
            nan_active_cells=int((nan_bounded & active).sum()),
            nonpositive_active=[int(((self.thickness[k] <= 0) & active[k]).sum())
                                for k in range(nlay)],
            thin=[int((self.thickness[k] < float(self.min_thickness[k])).sum())
                  for k in range(nlay)],
            pinched=[int((self.idomain[k] != 1).sum()) for k in range(nlay)],
            isolated_active=isolated,
            n_active_components=len(sizes),
            component_sizes=sorted(sizes.values(), reverse=True),
        )

    def validate(self) -> LayerBuildResult:
        """Raise ``ValueError`` if :meth:`qc` finds anything fatal; else return self."""
        report = self.qc()
        if not report.ok:
            raise ValueError("Layer stack failed QC:\n" + str(report))
        return self

    def prune_isolated(self) -> LayerBuildResult:
        """Deactivate (``idomain -> 0``) active cells that have no connection.

        Returns a new result. Isolated cells have no edges, so removing them
        never disconnects anything else -- a single pass is enough."""
        isolated = self.qc().isolated_active
        if not isolated:
            return self
        idomain = self.idomain.copy()
        for k, i in isolated:
            idomain[k, i] = 0
        return _dc_replace(self, idomain=idomain)

    # -- views (thin wrappers over flopy / the surface API) --------------- #
    def _resolve_line(self, line, x, y) -> dict:
        """Build a flopy cross-section line from explicit points, x=, y=, or a
        default West-East line through the grid centre."""
        if line is not None:
            return {"line": [tuple(pt) for pt in line]}
        minx, miny, maxx, maxy = (float(v) for v in self.vor.gdf_vorPolys.total_bounds)
        if x is not None:
            return {"line": [(float(x), miny), (float(x), maxy)]}
        if y is not None:
            return {"line": [(minx, float(y)), (maxx, float(y))]}
        ymid = (miny + maxy) / 2.0  # default: West-East section through the centre
        return {"line": [(minx, ymid), (maxx, ymid)]}

    def vertex_grid(self):
        """A flopy ``VertexGrid`` carrying this result's geometry."""
        import flopy

        p = self.vor.get_disv_gridprops()
        return flopy.discretization.VertexGrid(
            vertices=p["vertices"], cell2d=p["cell2d"],
            top=self.top, botm=self.botm, idomain=self.idomain,
            nlay=self.nlay, ncpl=p["ncpl"], crs=str(getattr(self.vor, "crs", None)),
        )

    def _draw_cross_section(
        self, line=None, *, x=None, y=None, color_by="layer", ax=None,
        cmap="tab10", show_grid=True, legend=True, title=None,
    ):
        """Render the layer cross-section into ``ax`` and return it.

        Private: this is the RENDERER. ``stack.plot.section(y=300)`` is the verb,
        and it returns a Picture that answers ``.show()``/``.save()``/``.html()``
        like every other. Matplotlib-native by necessity -- see
        :class:`LayerSection`.

        Layer coloring delegates to the shared
        :func:`~myflopy.modflow.mf6.cross_section_plotting.plot_layered_cross_section`
        renderer (the same core behind the model-aware ``plot_model_cross_section``)."""
        import matplotlib.pyplot as plt

        from myflopy.modflow.mf6.cross_section_plotting import (
            ModelCrossSectionStyle,
            plot_layered_cross_section,
        )

        vg = self.vertex_grid()
        line_spec = self._resolve_line(line, x, y)

        if color_by == "thickness":
            import flopy

            if ax is None:
                _, ax = mpl_axes(figsize=(9, 4))
            xsec = flopy.plot.PlotCrossSection(modelgrid=vg, line=line_spec, ax=ax)
            xsec.plot_array(self.thickness, cmap="viridis")
            if show_grid:
                xsec.plot_grid(lw=0.25, color="0.3")
            ax.set_title(title or "Layer cross-section")
            ax.set_xlabel(f"distance along section [{self.length_units}]")
            ax.set_ylabel(f"elevation [{self.length_units}]")
            return ax

        cm = plt.get_cmap(cmap)
        colors = [cm(i % 10) for i in range(self.nlay)]  # matplotlib RGBA tuples
        style = ModelCrossSectionStyle(
            figsize=(9, 4), grid_linewidth=0.25, grid_color="0.3",
            layer_alpha=0.75, title_fontsize=12, label_fontsize=10,
            legend_loc="upper right", legend_frameon=True, legend_fontsize=8,
        )
        _, ax = plot_layered_cross_section(
            vg, line_spec, ax=ax, style=style,
            layer_colors=colors, layer_labels=self.names,
            show_grid=show_grid, show_layers=True, show_head=False,
            show_legend=legend, title=title or "Layer cross-section",
            xlabel=f"distance along section [{self.length_units}]",
            ylabel=f"elevation [{self.length_units}]",
        )
        return ax

    @property
    def surface_names(self) -> list[str]:
        """Plottable surface names: the model top plus each layer's bottom contact."""
        return ["top"] + list(self.names)

    def _surface_z(self, name) -> np.ndarray:
        """Elevation array for a surface name (``"top"`` or a layer's bottom)."""
        if name == "top":
            return self.top
        return self.botm[self.names.index(name)]

    def _resolve_surface_names(self, layer) -> list[str]:
        """Normalize the ``layer`` selector to a list of valid surface names."""
        valid = self.surface_names
        if isinstance(layer, str) and layer == "all":
            return valid
        names = list(layer) if isinstance(layer, (list, tuple)) else [layer]
        bad = [n for n in names if n not in valid]
        if bad:
            raise KeyError(f"no surface named {bad!r}; choose from {valid} or 'all'.")
        return names

    @staticmethod
    def _rgb(rgba) -> str:
        """Format a 0-1 RGBA tuple as a Plotly ``rgb(r,g,b)`` string (0-255)."""

        r, g, b = (int(round(255 * c)) for c in rgba[:3])
        return f"rgb({r},{g},{b})"

    def _surface_fig(
        self, layer="top", *, resolution=120, colorscale="Earth_r",
        color_by=None, opacity=None, height=None,
    ):
        """The 3-D surface figure. Reached as ``stack.plot.surface(...)``.

        Private because it is Layer 1 only: it BUILDS the figure and nothing
        else. Writing it out and opening it are `.html(path)` / `.show()` on the
        Picture that wraps it -- this used to take `html_path=`/`browser=` and do
        both itself.

        Interactive 3D surface(s) of one or more layers (plotly).

        ``layer`` selects which surface(s) to draw and may be:

        * ``"top"`` or any layer name -- that layer's *bottom* contact (default ``"top"``),
        * a list of names, e.g. ``["top", "clay", "bedrock"]``, overlaid in one scene,
        * ``"all"`` -- the model top plus every layer bottom.

        Discover the choices with ``result.surface_names``. A single surface with
        relief is shaded by elevation (``colorscale``); a flat surface or several
        surfaces get distinct solid colors. Override with
        ``color_by="elevation"`` / ``"surface"`` and tune ``opacity``.

        ``height`` sets the figure height in pixels (default ``None`` = fill the
        container / browser window)."""

        from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface

        names = self._resolve_surface_names(layer)
        crs = str(getattr(self.vor, "crs", None))
        xs = np.asarray(self.vor.centroids[0])
        ys = np.asarray(self.vor.centroids[1])

        # Each surface is built by the shared InterpolatedSurface.surface_trace
        # so there is one go.Surface builder across the toolkit; we keep the
        # interpolated z-array to set a common elevation scale below.
        grids = []
        for name in names:
            isurf = InterpolatedSurface(
                xs=xs, ys=ys, zs=np.asarray(self._surface_z(name)),
                surf_type="lyr", resolution=resolution, crs=crs,
            )
            grids.append((name, isurf, isurf.surface))

        # Overall elevation span -> colour mode, shared colour range and z-axis.
        zmin = min(float(np.nanmin(zz)) for *_, zz in grids)
        zmax = max(float(np.nanmax(zz)) for *_, zz in grids)
        flat = zmax - zmin < 1e-6
        if color_by is None:
            color_by = "elevation" if (len(names) == 1 and not flat) else "surface"
        if opacity is None:
            opacity = 1.0 if len(names) == 1 else 0.85
        zaxis = dict(title=f"elev [{self.length_units}]")
        aspectmode = "manual"
        if flat:
            # A perfectly flat surface has zero vertical extent; a lone flat
            # go.Surface then fails to render in hardware-WebGL viewers (its faces
            # get a 0/0 colour and the trace has no z-depth), while multi-surface
            # scenes dodge this because their combined z-range is non-zero. Give
            # the sheet a faint, non-planar relief and a real z-axis so it always
            # draws; ``cube`` aspect keeps the near-flat sheet prominent.
            mid = 0.5 * (zmin + zmax)
            pad = max(1.0, abs(mid) * 0.01)
            gx0, gy0 = grids[0][1].xy_meshgrid
            xspan = float(np.nanmax(gx0) - np.nanmin(gx0)) or 1.0
            yspan = float(np.nanmax(gy0) - np.nanmin(gy0)) or 1.0
            ripple = (
                ((gx0 - np.nanmin(gx0)) / xspan - 0.5)
                + ((gy0 - np.nanmin(gy0)) / yspan - 0.5)
            ) * (0.5 * pad)                      # ~+/-pad/2 of imperceptible relief
            grids = [(n, s, zz + ripple) for (n, s, zz) in grids]
            zmin, zmax = mid - pad, mid + pad
            zaxis["range"] = [zmin, zmax]
            aspectmode = "cube"

        fig = Fig()
        if color_by == "elevation":
            for i, (name, isurf, zz) in enumerate(grids):
                fig.add_trace(isurf.surface_trace(
                    surface=zz, colorscale=colorscale, name=name,
                    cmin=zmin, cmax=zmax, opacity=opacity, showscale=(i == 0),
                    colorbar=dict(title=f"elev [{self.length_units}]"),
                ))
        else:  # one solid color per surface, distinguished by a legend
            import matplotlib.pyplot as plt

            cmap = plt.get_cmap("tab10")
            for i, (name, isurf, zz) in enumerate(grids):
                c = self._rgb(cmap(i % 10))
                fig.add_trace(isurf.surface_trace(
                    surface=zz, colorscale=[[0, c], [1, c]], name=name,
                    cmin=zmin, cmax=zmax, opacity=opacity,
                    showscale=False, showlegend=True,
                ))

        multi = len(names) > 1
        scene = dict(aspectmode=aspectmode, zaxis=zaxis)
        if aspectmode == "manual":
            scene["aspectratio"] = dict(x=1, y=0.6, z=0.45)
        fig.update_layout(
            title=("layer surfaces" if multi else f"{names[0]} surface"),
            autosize=True, height=height, margin=dict(l=0, r=0, t=40, b=0),
            showlegend=(color_by == "surface" and multi), scene=scene,
        )
        return fig

    def _resolve_layer_indices(self, layers) -> list[int]:
        """Normalize a layer selector to sorted, unique layer indices.

        Accepts ``None``/``"all"`` (every layer), a single layer name or integer
        index, or a list mixing names and indices."""
        if layers is None or (isinstance(layers, str) and layers == "all"):
            return list(range(self.nlay))
        items = list(layers) if isinstance(layers, (list, tuple)) else [layers]
        idx = []
        for it in items:
            if isinstance(it, (int, np.integer)) and not isinstance(it, bool):
                i = int(it)
                if not 0 <= i < self.nlay:
                    raise IndexError(f"layer index {i} out of range [0, {self.nlay}).")
                idx.append(i)
            elif it in self.names:
                idx.append(self.names.index(it))
            else:
                raise KeyError(f"no layer named {it!r}; choose from {self.names} or an index.")
        return sorted(set(idx))

    def _vtk_plotter(
        self, layers=None, *, color_by="layer", scale=8, cmap="tab10",
        width=900, height=580,
    ):
        """Build the PyVista plotter for the layered grid volume.

        Private: this is Layer 1 only -- it BUILDS the scene. Displaying it and
        writing it out belong to the :class:`~myflopy.viz.VtkScene` that wraps it.
        Reached as ``stack.plot.grid(backend="vtk")``.

        ``layers`` selects which layers to show: ``None``/``"all"`` (default) for
        every layer, a single layer name or index, or a list mixing them
        (e.g. ``["sand", "clay"]`` or ``[0, 2]``). Colors stay keyed to each
        layer's position, so a subset keeps the same colors it has in the full
        stack."""
        import tempfile
        from pathlib import Path

        import flopy
        from flopy.export.vtk import Vtk

        pv = require("pyvista", feature="interactive 3-D scenes")

        sel = self._resolve_layer_indices(layers)
        p = self.vor.get_disv_gridprops()
        ws = Path(tempfile.mkdtemp(prefix="layer_view_"))
        sim = flopy.mf6.MFSimulation(sim_name="layerview", sim_ws=str(ws))
        flopy.mf6.ModflowTdis(sim)
        flopy.mf6.ModflowIms(sim)
        gwf = flopy.mf6.ModflowGwf(sim, modelname="layers")
        flopy.mf6.ModflowGwfdisv(
            gwf, nlay=self.nlay, ncpl=p["ncpl"], nvert=len(p["vertices"]),
            vertices=p["vertices"], cell2d=p["cell2d"],
            top=self.top, botm=self.botm, idomain=self.idomain,
        )
        vtk = Vtk(model=gwf, vertical_exageration=1, binary=True, smooth=False)
        vtk.add_model(gwf)
        # Always attach a per-cell layer index so a subset can be selected.
        vtk.add_array(np.repeat(np.arange(self.nlay)[:, None], p["ncpl"], axis=1), "layer")
        mesh = vtk.to_pyvista()
        if isinstance(mesh, pv.MultiBlock):
            mesh = mesh.combine()
        if len(sel) != self.nlay:  # keep only the requested layers' cells
            layer_cell = np.asarray(mesh.cell_data["layer"]).astype(int)
            mesh = mesh.extract_cells(np.isin(layer_cell, sel))
        if color_by == "layer":
            scalars, use_cmap = "layer", cmap
        else:
            mesh["elevation"] = mesh.points[:, 2]
            scalars, use_cmap = "elevation", "terrain"
        mesh_kwargs = dict(show_edges=True, cmap=use_cmap)
        if color_by == "layer":
            # Discrete colors keyed to each layer's global index; label only the
            # layers actually shown.
            mesh_kwargs.update(
                n_colors=self.nlay,
                clim=[-0.5, self.nlay - 0.5],
                annotations={float(i): self.names[i] for i in sel},
                scalar_bar_args=dict(title="layer", n_labels=0),
            )
        plotter = pv.Plotter(off_screen=True, window_size=(width - 40, height - 40))
        plotter.add_mesh(mesh, scalars=scalars, **mesh_kwargs)
        plotter.set_scale(zscale=scale)
        plotter.add_axes()
        plotter.camera_position = "yz"
        return plotter

    def _thickness_values(self, layer=None):
        """Per-cell thickness (total, or a single named layer) and its label."""

        if layer is None:
            return self.thickness.sum(axis=0), "total thickness"
        return self.thickness[self.names.index(layer)], f"{layer!r} thickness"

    def _draw_thickness_map(self, *, layer=None, ax=None):
        """Render the thickness choropleth into ``ax`` (geopandas, no basemap)."""

        values, title = self._thickness_values(layer)
        gdf = self.vor.gdf_vorPolys.copy().assign(_thickness=values)
        if ax is None:
            _, ax = mpl_axes()
        gdf.plot(column="_thickness", ax=ax, legend=True)
        ax.set_title(f"{title} [{self.length_units}]")
        ax.set_aspect("equal")
        return ax

    @property
    def plot(self) -> StackPlots:
        """The plotting verbs for this layer geometry: ``map``, ``section``,
        ``surface``, ``grid``.

        Same verbs as ``model.plot`` and ``vor.plot``, and everything returned is
        a :class:`~myflopy.viz.Picture`. A stack has no results and no time, so
        there is no ``animate``; ``qc()`` stays a report rather than becoming a
        picture, because it is text you read.

        ``surface`` is a height field -- one contact as ``z(x, y)``. The layered
        cell VOLUME is ``grid(backend="vtk")``, a different shape entirely.

        Replaces ``thickness_map()``/``preview()`` (now ``plot.map()``),
        ``cross_section()`` (``plot.section()``), ``surface_3d()``
        (``plot.surface()``), ``vtk_3d()`` (``plot.grid()``) and ``views()``
        (compose what you want with ``myflopy.plot.mosaic``).

        Examples
        --------
        >>> stack.plot.map()                          # total thickness
        >>> stack.plot.map("sand", basemap=True)      # one layer, on a basemap
        >>> stack.plot.section(y=300)
        >>> stack.plot.surface("all")
        >>> stack.plot.grid(["sand", "clay"], scale=12)
        """

        return StackPlots(self)


class LayerSection(MplPicture):
    """A filled, layer-coloured cross-section through a built stack.

    Matplotlib by necessity, not by preference: the renderer is FloPy's
    ``PlotCrossSection``, and the Plotly ``GridSection`` draws cell outlines
    without layer fills. Rather than exempt it from the picture grammar, it is an
    :class:`~myflopy.viz.MplPicture` -- same verbs, Axes underneath.
    """

    def __init__(self, result, line=None, *, x=None, y=None, **kwargs):
        """Bind a section of ``result`` along ``line`` (or ``x=``/``y=``)."""

        self._result = result
        self._line, self._x, self._y = line, x, y
        self._kwargs = kwargs
        self.title = kwargs.get("title") or "Layer cross-section"

    def draw(self, ax=None, **kwargs):
        """Render the section into ``ax`` (or a new one) and return the Axes."""

        return self._result._draw_cross_section(
            self._line, x=self._x, y=self._y, ax=ax, **{**self._kwargs, **kwargs}
        )


class LayerThicknessMap(MplPicture):
    """Per-cell layer thickness, drawn on the grid without a basemap.

    Deliberately NOT a :class:`Choro`: a stack under construction is often on
    synthetic or local coordinates, and a web basemap would put it in the ocean.
    ``stack.plot.map(basemap=True)`` gives the georeferenced choropleth when the
    grid really is where it says it is.
    """

    def __init__(self, result, layer=None):
        """Bind a thickness map of ``result`` (total, or one named layer)."""

        self._result = result
        self._layer = layer
        _, self.title = result._thickness_values(layer)

    def draw(self, ax=None, **kwargs):
        """Render the thickness map into ``ax`` (or a new one)."""

        return self._result._draw_thickness_map(layer=self._layer, ax=ax, **kwargs)


class LayerSurface(Picture):
    """One or more layer contacts as an interactive 3-D Plotly surface."""

    def __init__(self, result, layer="top", **kwargs):
        """Bind a 3-D surface of ``result`` for the chosen contact(s)."""

        self._result = result
        self._layer = layer
        self._kwargs = kwargs
        self._built = None

    @property
    def fig(self) -> Fig:
        """The assembled 3-D figure (built once, then cached)."""

        if self._built is None:
            self._built = self._result._surface_fig(self._layer, **self._kwargs)
        return self._built


#: `_vtk_plotter`'s own defaults, so `grid(backend="plotly")` can tell a value
#: the caller chose from one it merely inherited. Mirrored, not imported, because
#: the signature is the public contract and a test pins the two equal.
_VTK_GRID_DEFAULTS = {
    "color_by": "layer", "scale": 8, "cmap": "tab10", "width": 900, "height": 580,
}


class StackPlots:
    """The plotting verbs for layer geometry -- ``stack.plot.map()`` (plan 8.5a).

    Three verbs, because a layer stack can answer three questions: how thick is
    it (:meth:`map`), what does it look like in section (:meth:`section`), and
    what shape is a given contact (:meth:`surface`). No ``animate`` -- a stack has
    no stress periods. ``qc()`` stays a method on the stack, because a QC report
    is text you read, not a picture you look at.

    Everything returned is a :class:`~myflopy.viz.Picture`, so it renders inline
    and answers ``.show()`` / ``.save(path)`` / ``.html(path)``. Two of the three
    are Matplotlib underneath and answer those over ``.axes``; only
    :meth:`surface` has a Plotly ``.fig``.
    """

    def __init__(self, result):
        """Bind the plotting verbs to a built :class:`LayerBuildResult`."""

        self.result = result

    def __repr__(self):
        """Name the verbs, since tab-completion is how this gets found."""

        return f"StackPlots({self.result.nlay} layers: map, section, surface, grid)"

    def map(self, layer=None, *, basemap: bool = False, **kwargs):
        """Per-cell thickness -- total, or one named layer.

        Draws on the grid with no basemap by default, because a stack is often
        still on synthetic coordinates. ``basemap=True`` routes through the
        shared choropleth instead, for a grid that really is georeferenced.

        Parameters
        ----------
        layer : str or int, optional
            One layer by name or index. With none, total thickness.
        basemap : bool, default False
            Route through the shared choropleth (``vor.plot.map``) so the cells
            sit on a web basemap. Requires a real CRS.
        **kwargs
            Choropleth styling, forwarded to :meth:`~myflopy.modflow.mf6.grid
            .plotting.GridPlots.map`. **Only meaningful with ``basemap=True``**
            -- the default renderer is Matplotlib and takes none of them.

        Raises
        ------
        TypeError
            If styling arguments are given without ``basemap=True``. They used
            to be accepted and silently discarded, which quietly produced an
            unstyled picture.
        """

        if basemap:
            values, _ = self.result._thickness_values(layer)
            # `GridPlots`, not `myflopy.plot`: that module is a layer ABOVE this
            # one in the import graph, so reaching it would need a deferred
            # import and the exact-match ratchet only moves down. Same function
            # either way -- `vor.plot.map` is what the front door calls too.
            return GridPlots(self.result.vor).map(values=list(values), **kwargs)
        if kwargs:
            raise TypeError(
                f"{', '.join(sorted(kwargs))} style the choropleth, which is "
                f"only drawn with basemap=True; the default thickness map is "
                f"Matplotlib and ignores them."
            )
        return LayerThicknessMap(self.result, layer=layer)

    def section(
        self,
        line=None,
        *,
        x=None,
        y=None,
        color_by: str = "layer",
        cmap: str = "tab10",
        show_grid: bool = True,
        legend: bool = True,
        title: str | None = None,
        **kwargs,
    ) -> LayerSection:
        """A filled, layer-coloured cross-section: ``stack.plot.section(y=300)``.

        Parameters
        ----------
        line : LineString or Path, optional
            The section line. Alternatively give ``x=`` or ``y=`` for an
            axis-aligned slice.
        x, y : float, optional
            Draw the section along a constant x or constant y.
        color_by : str, default 'layer'
            Cell scalar the fill is keyed to.
        cmap : str, default 'tab10'
            Colormap for that scalar.
        show_grid : bool, default True
            Draw cell edges over the fill.
        legend : bool, default True
            Include the layer legend.
        title : str, optional
            Plot title. Defaults to "Layer cross-section".
        **kwargs
            Forwarded to :class:`LayerSection`.

        Returns
        -------
        LayerSection
            A Matplotlib :class:`~myflopy.viz.Picture`; ``.axes`` rather than
            ``.fig``.
        """

        return LayerSection(
            self.result, line, x=x, y=y, color_by=color_by, cmap=cmap,
            show_grid=show_grid, legend=legend, title=title, **kwargs,
        )

    def surface(
        self,
        layer="top",
        *,
        resolution: int = 120,
        colorscale: str = "Earth_r",
        color_by: str | None = None,
        opacity: float | None = None,
        height: int | None = None,
        **kwargs,
    ) -> LayerSurface:
        """One or more contacts as an interactive 3-D surface.

        ``layer`` is a name, a list of names, or ``"all"``; discover them with
        ``result.surface_names``. That is a height field ``z(x, y)``; the layered
        grid VOLUME is :meth:`grid`, a different shape entirely.

        Parameters
        ----------
        layer : str or list of str, default 'top'
            Which contact(s) to draw. ``"all"`` draws every one.
        resolution : int, default 120
            Interpolation grid size per axis.
        colorscale : str, default 'Earth_r'
            Plotly colorscale for the height field.
        color_by : str, optional
            Colour by a scalar other than elevation.
        opacity : float, optional
            Surface opacity, useful when stacking several contacts.
        height : int, optional
            Figure height in pixels.
        **kwargs
            Forwarded to :class:`LayerSurface`.

        Returns
        -------
        LayerSurface
            A Plotly :class:`~myflopy.viz.Picture`.
        """

        return LayerSurface(
            self.result, layer, resolution=resolution, colorscale=colorscale,
            color_by=color_by, opacity=opacity, height=height, **kwargs,
        )

    def grid(
        self,
        layers=None,
        *,
        backend: str = "vtk",
        color_by: str = "layer",
        scale: float = 8,
        cmap: str = "tab10",
        width: int = 900,
        height: int = 580,
    ):
        """The layered grid mesh itself.

        ``backend="vtk"`` (the default here) renders the cell VOLUME in 3-D,
        coloured by layer -- the picture the old ``vtk_3d()`` drew, minus its
        habit of writing an HTML file into the working directory on every call.
        ``backend="plotly"`` gives the flat 2-D mesh instead, the same picture as
        ``vor.plot.grid()``.

        A ``backend`` switch is honest here because both branches draw the SAME
        subject -- this grid -- and differ only in renderer. That is why the 3-D
        volume is `grid`, not `surface`: `surface` means a height field.

        Parameters
        ----------
        layers : str or int or list, optional
            *(vtk only)* Which layers to show: a name, an index, or a list mixing
            them. Colours stay keyed to each layer's position, so a subset looks
            the same as it does in the full stack.
        backend : {'vtk', 'plotly'}, default 'vtk'
            ``'vtk'`` renders the 3-D volume and needs the ``viz3d`` extra;
            ``'plotly'`` draws the flat 2-D mesh.
        color_by : str, default 'layer'
            *(vtk only)* Cell scalar to colour by.
        scale : float, default 8
            *(vtk only)* Vertical exaggeration.
        cmap : str, default 'tab10'
            *(vtk only)* Colormap for ``color_by``.
        width, height : int, default 900, 580
            *(vtk only)* Scene size in pixels.

        Returns
        -------
        VtkScene or GridMesh
            A :class:`~myflopy.viz.Picture` either way. ``VtkScene`` exposes
            ``.scene`` instead of ``.fig``.

        Raises
        ------
        ValueError
            If ``backend`` is neither value, or if a vtk-only argument is given
            with ``backend="plotly"``.
        """

        scene_args = {"color_by": color_by, "scale": scale, "cmap": cmap,
                      "width": width, "height": height}
        if backend == "plotly":
            stray = [n for n, v in scene_args.items() if v != _VTK_GRID_DEFAULTS[n]]
            if layers is not None:
                stray.insert(0, "layers")
            if stray:
                raise ValueError(
                    f"{', '.join(stray)} configure the 3-D scene; the flat mesh "
                    f"has no use for them. Drop backend='plotly', or drop these."
                )
            return GridPlots(self.result.vor).grid()
        if backend != "vtk":
            raise ValueError(
                f"backend must be 'vtk' or 'plotly', not {backend!r}."
            )
        return VtkScene(
            self.result._vtk_plotter(layers, **scene_args),
            title="layered grid",
        )


class LayerStack:
    """Author a model's vertical layering from a top surface + named layers.

    ``LayerStack`` is the **friendly facade** for building MODFLOW layer geometry.
    You declare the model top, then ``.add(...)`` one named layer at a time (by its
    *bottom* surface or its *thickness*), and ``.build()`` resolves the stack into
    DISV-ready ``top`` / ``botm`` / ``idomain`` arrays. It is a thin layer over the
    :class:`~myflopy.surfaces.LayerSurfaces` engine -- it compiles to it (see
    :meth:`_layer_surfaces`) and never re-implements the sampling / reconcile /
    pinch-out logic. Each surface is an atomic :class:`~myflopy.surfaces.Surface`
    (use the aliases :class:`Raster`, :class:`Contours`, :class:`Points`,
    :class:`Flat`, :class:`Array`).

    What :meth:`build` does for you: samples every surface onto the grid cells
    (area-weighted by default), resolves them **top-down** so a flat/relative
    surface stays flat where it fits and is lowered only where it would intrude on
    the surface above ("flat where possible, fit between"), and turns sub-minimum
    or inverted thicknesses into **pinch-outs** (idomain) per your per-layer policy.

    Parameters
    ----------
    vor
        The grid (a ``VoronoiGridPlus``/vertex grid) whose cells the surfaces are
        sampled onto.
    top
        The model-top surface -- a :class:`Surface` or any value the aliases accept
        (a raster path via :class:`Raster`, an array via :class:`Array`, a constant
        via :class:`Flat`, ...).
    length_units, time_units
        Model units (default feet / days); ``length_units`` flows to DISV and is
        used to convert any surface declaring different ``units=``.

    Examples
    --------
    >>> from myflopy import LayerStack
    >>> from myflopy.layers import Raster, Contours
    >>> stack = (
    ...     LayerStack(vor, top=Raster("ground.tif"))   # land surface from a DEM
    ...     .add("alluvium",   thickness=25)            # 25-ft upper aquifer
    ...     .add("aquitard",   thickness=8, pinch="inactive")   # pinches out where thin
    ...     .add("bedrock",    bottom=Contours("bedrock_top.shp"))
    ... )
    >>> layers = stack.build()          # -> layers.top, layers.botm, layers.idomain, layers.nlay
    >>> print(stack.qc())               # NaN/thickness/connectivity report
    >>> layers.cross_section(x=500)     # quick W-E section to eyeball it

    Feed the result straight into ``mf.disv`` and the model context::

        gp = vor.get_disv_gridprops()
        ctx = mf.ModelContext(grid=vor, domain=layers.idomain)
        flow = mf.gwf("flow", context=ctx, packages=[
            mf.disv(nlay=layers.nlay, ncpl=gp["ncpl"], nvert=len(gp["vertices"]),
                    vertices=gp["vertices"], cell2d=gp["cell2d"],
                    top=layers.top, botm=layers.botm, idomain=layers.idomain),
            ...,
        ])

    See Also
    --------
    from_modflow : Seed an editable stack from an existing model's top/botm.
    build : Resolve the stack into DISV arrays (the per-layer pinch options).
    qc : Geometry quality-control report (NaN, thin/pinched, connectivity).
    """

    def __init__(self, vor, top, *, length_units: str = "feet", time_units: str = "days"):
        """Start an editable layer stack on grid ``vor`` with the model-top surface ``top``.

        Layers are added below the top with :meth:`add`; ``length_units`` /
        ``time_units`` are carried onto the built :class:`LayerBuildResult`.
        """

        self.vor = vor
        self._top = _coerce_surface(top)
        self._layers: list[_Layer] = []
        self.length_units = length_units
        self.time_units = time_units

    @classmethod
    def from_modflow(
        cls,
        vor,
        source,
        *,
        names: list[str] | None = None,
        resample: bool = True,
        length_units: str = "feet",
        time_units: str = "days",
    ) -> LayerStack:
        """Build a stack on ``vor`` from an existing MODFLOW model's top/botm.

        ``source`` is a flopy model or modelgrid; its top and per-layer bottoms
        become this stack's surfaces (interpolated onto ``vor`` when
        ``resample=True``, used verbatim when ``False``). Name the layers with
        ``names`` (defaults to ``layer1..N``). Edit the returned stack like any
        other -- e.g. ``.replace(...)`` a bottom, then ``.build()`` / ``.to_disv()``.
        """
        top_s, botm_s = modflow_surfaces(source, resample=resample)
        if names is None:
            names = [f"layer{i + 1}" for i in range(len(botm_s))]
        if len(names) != len(botm_s):
            raise ValueError(f"{len(names)} names given for {len(botm_s)} layers.")
        stack = cls(vor, top=top_s, length_units=length_units, time_units=time_units)
        for name, surface in zip(names, botm_s, strict=False):
            stack.add(name, bottom=surface)
        return stack

    # -- authoring -------------------------------------------------------- #
    def add(
        self,
        name: str,
        *,
        bottom=None,
        thickness=None,
        min_thickness: float | None = None,
        pinch: str | None = None,
        fill: str | None = None,
    ) -> LayerStack:
        """Append a named layer beneath the current bottom. Returns ``self`` (chainable).

        Define the layer either by its **bottom** surface or its **thickness**
        (exactly one). Thickness is measured down from the surface above, so layers
        compose naturally as you stack them.

        Parameters
        ----------
        name
            Layer name (used in QC, plots, and as the surface label).
        bottom
            The layer's bottom as a :class:`Surface` / alias (e.g.
            ``Contours("base.shp")``, ``Array(values)``, ``Flat(90.0)``). Mutually
            exclusive with ``thickness``.
        thickness
            Constant or per-cell thickness below the surface above. Mutually
            exclusive with ``bottom``.
        min_thickness
            Minimum layer thickness; thinner cells are handled per ``pinch``.
            Defaults to the stack-wide value passed to :meth:`build`.
        pinch
            What to do where a layer is thinner than ``min_thickness``:
            ``"passthrough"`` (keep the cell active, default), ``"inactive"``
            (idomain 0 -- a true pinch-out), or ``"floor"`` (clamp to the minimum).
        fill
            For raster/derived bottoms, how to fill cells with no source data
            (e.g. ``"propagate"`` to inherit the surface above -> pinch).

        Examples
        --------
        >>> stack.add("sand", thickness=20)                       # 20-ft layer
        >>> stack.add("clay", thickness=5, pinch="inactive")      # pinches out where thin
        >>> stack.add("bedrock", bottom=Contours("bedrock.shp"))  # bottom from contours
        """
        if any(layer.name == name for layer in self._layers):
            raise ValueError(f"layer {name!r} already exists.")
        surface = _make_surface(bottom, thickness, fill)
        self._layers.append(_Layer(name, surface, min_thickness, pinch))
        return self

    def insert_below(
        self, name: str, new_name: str, *, bottom=None, thickness=None,
        min_thickness: float | None = None, pinch: str | None = None, fill: str | None = None,
    ) -> LayerStack:
        """Insert a new layer directly below the existing layer ``name``."""
        if any(layer.name == new_name for layer in self._layers):
            raise ValueError(f"layer {new_name!r} already exists.")
        surface = _make_surface(bottom, thickness, fill)
        self._layers.insert(
            self._index(name) + 1, _Layer(new_name, surface, min_thickness, pinch)
        )
        return self

    def replace(
        self, name: str, *, bottom=None, thickness=None,
        min_thickness=_UNSET, pinch=_UNSET, fill=None,
    ) -> LayerStack:
        """Update an existing layer in place; unspecified fields are kept."""
        idx = self._index(name)
        layer = self._layers[idx]
        if bottom is not None or thickness is not None:
            surface = _make_surface(bottom, thickness, fill)
        else:
            surface = layer.surface if fill is None else _dc_replace(layer.surface, fill=fill)
        self._layers[idx] = _Layer(
            name,
            surface,
            layer.min_thickness if min_thickness is _UNSET else min_thickness,
            layer.pinch if pinch is _UNSET else pinch,
        )
        return self

    def remove(self, name: str) -> LayerStack:
        """Remove a layer by name."""
        del self._layers[self._index(name)]
        return self

    @property
    def names(self) -> list[str]:
        """The layer names, top-to-bottom (excluding the model top)."""

        return [layer.name for layer in self._layers]

    def _index(self, name: str) -> int:
        """The position of the layer named ``name`` (raises ``KeyError`` if absent)."""

        for i, layer in enumerate(self._layers):
            if layer.name == name:
                return i
        raise KeyError(f"no layer named {name!r} (have {self.names}).")

    # -- compilation ------------------------------------------------------ #
    def _layer_surfaces(self) -> LayerSurfaces:
        """Compile to the :class:`LayerSurfaces` engine: the top plus each layer bottom, labeled."""

        surfaces = [self._top] + [layer.surface for layer in self._layers]
        labels = ["top"] + self.names
        return LayerSurfaces(surfaces, labels=labels)

    def _per_layer_config(self, default_min_thickness, default_pinch):
        """Resolve each layer's ``(min_thickness, pinch)``, filling unset ones with the defaults."""

        min_thk = [
            default_min_thickness if layer.min_thickness is None else layer.min_thickness
            for layer in self._layers
        ]
        pinch = [
            default_pinch if layer.pinch is None else layer.pinch
            for layer in self._layers
        ]
        return min_thk, pinch

    def refresh(self) -> LayerStack:
        """Rebuild the cached raster of every derived (contour) surface now."""
        for surface in [self._top] + [layer.surface for layer in self._layers]:
            if getattr(surface, "is_derived", False):
                surface.resolve_source(refresh=True, warn=False)
        return self

    def cache_status(self) -> dict[str, str]:
        """Map ``name -> "missing"/"fresh"/"stale"`` for each derived surface."""
        named = [("top", self._top)] + [(layer.name, layer.surface) for layer in self._layers]
        return {
            name: surface.cache_status()
            for name, surface in named
            if getattr(surface, "is_derived", False)
        }

    def build(
        self,
        *,
        default_min_thickness: float = 1.0,
        default_pinch: str = "passthrough",
        reconcile="bottom",
        min_sep: float = 0.1,
        method: str = "area",
        refresh: bool = False,
        attach: bool = False,
    ) -> LayerBuildResult:
        """Resolve, reconcile, and pinch out the stack into DISV-ready arrays.

        Samples every surface onto the grid, reconciles crossing surfaces, applies
        the per-layer pinch policy, and returns the ``top``/``botm``/``idomain``
        arrays ready for ``mf.disv(...)``. A stale derived (contour) surface is
        reused with a warning; pass ``refresh=True`` to rebuild its cache first.

        Parameters
        ----------
        default_min_thickness : float, default 1.0
            Minimum layer thickness used where a layer does not set its own.
        default_pinch : str, default "passthrough"
            Default thin-layer policy: ``"passthrough"`` (idomain -1),
            ``"inactive"`` (idomain 0, a true pinch-out), or ``"floor"`` (clamp to
            the minimum).
        reconcile : str, default "bottom"
            How crossing surfaces are reconciled (e.g. push conflicts to the
            ``"bottom"``).
        min_sep : float, default 0.1
            Minimum vertical separation enforced between reconciled surfaces.
        method : str, default "area"
            Raster sampling method: ``"area"`` (area-weighted) or ``"centroid"``.
        refresh : bool, default False
            Rebuild any stale derived (contour) surface caches before sampling.
        attach : bool, default False
            Also publish the result onto ``vor.gdf_topbtm`` (see
            :meth:`LayerBuildResult.attach_to_grid`) so surface-aware builders --
            ``mf.sfr`` reach tops, ``mf.lak`` lake-cell layering -- can read the
            elevations. Use it when your model has SFR/LAK on this grid.

        Returns
        -------
        LayerBuildResult
            Bundles ``top`` ``(ncpl,)``, ``botm`` / ``idomain`` ``(nlay, ncpl)``,
            derived ``thickness``, and per-layer metadata.

        Raises
        ------
        ValueError
            If no layers have been added with :meth:`add`.

        Examples
        --------
        >>> layers = (mf.LayerStack(vor, top=mf.Raster("ground.tif"))
        ...           .add("sand", thickness=20)
        ...           .add("clay", bottom=mf.Contours("base.shp"), pinch="inactive")
        ...           .build(attach=True))
        >>> mf.disv(nlay=layers.nlay, ..., top=layers.top, botm=layers.botm,
        ...         idomain=layers.idomain)
        """
        if not self._layers:
            raise ValueError("Add at least one layer with .add(...) before build().")
        ls = self._layer_surfaces()
        rec_on, which = _reconcile_args(reconcile)
        gdf = ls.sample(
            self.vor, reconcile=rec_on, which=which, min_sep=min_sep,
            method=method, length_units=self.length_units, refresh=refresh,
        )
        top, botm = ls._split_top_botm(gdf)
        thickness = ls._thickness(gdf)
        min_thk, pinch = self._per_layer_config(default_min_thickness, default_pinch)
        ls._validate_pinch_invariant(
            min_thk, pinch, thickness.shape[0],
            {"reconcile": rec_on, "min_sep": min_sep},
        )
        idomain = ls._idomain_from_thickness(thickness, min_thk, pinch)
        result = LayerBuildResult(
            top=top, botm=botm, idomain=idomain, thickness=thickness,
            names=self.names, min_thickness=min_thk, pinch=pinch,
            length_units=self.length_units, time_units=self.time_units,
            vor=self.vor,
        )
        if attach:
            result.attach_to_grid(self.vor)
        return result

    def qc(
        self,
        *,
        default_min_thickness: float = 1.0,
        default_pinch: str = "passthrough",
        reconcile="bottom",
        min_sep: float = 0.1,
        method: str = "area",
        refresh: bool = False,
    ) -> LayerQCReport:
        """Build the stack and run QC, including reconcile diagnostics.

        On top of :meth:`LayerBuildResult.qc` (NaN coverage, isolated active
        cells, thickness), this samples the surfaces *without* reconciling and
        reports, per layer, how many cells reconcile had to move and the largest
        move -- showing where surfaces were crossing before reconcile fixed them.
        Returns a :class:`LayerQCReport`."""
        result = self.build(
            default_min_thickness=default_min_thickness, default_pinch=default_pinch,
            reconcile=reconcile, min_sep=min_sep, method=method, refresh=refresh,
        )
        report = result.qc()
        ls = self._layer_surfaces()
        raw = ls.sample(
            self.vor, reconcile=False, method=method,
            length_units=self.length_units, refresh=refresh,
        )
        _, raw_botm = ls._split_top_botm(raw)
        adjusted, max_shift = [], []
        for k in range(result.nlay):
            shift = np.abs(np.asarray(raw_botm[k], float) - np.asarray(result.botm[k], float))
            finite = shift[np.isfinite(shift)]
            adjusted.append(int((finite > 1e-6).sum()))
            max_shift.append(float(finite.max()) if finite.size else 0.0)
        report.reconcile_adjusted = adjusted
        report.reconcile_max_shift = max_shift
        return report

    # -- views: build with defaults, then view the result ----------------- #
    @property
    def plot(self) -> StackPlots:
        """The plotting verbs for this stack: ``map``, ``section``, ``surface``.

        Builds the stack with default options and returns the result's namespace,
        so ``stack.plot.map()`` is ``stack.build().plot.map()``. Build with
        non-default options first if you need them.
        """

        return self.build().plot


    def to_disv(
        self,
        vor=None,
        *,
        name: str = "disv",
        attach: bool = True,
        default_min_thickness: float = 1.0,
        default_pinch: str = "passthrough",
        reconcile="bottom",
        min_sep: float = 0.1,
        method: str = "area",
        refresh: bool = False,
    ):
        """Build a ready-to-use ``mf.disv`` spec (with pinch-out idomain)."""
        if not self._layers:
            raise ValueError("Add at least one layer before to_disv().")
        vor = self.vor if vor is None else vor
        ls = self._layer_surfaces()
        rec_on, which = _reconcile_args(reconcile)
        min_thk, pinch = self._per_layer_config(default_min_thickness, default_pinch)
        return ls.to_disv(
            vor,
            pinch_out=True,
            minimum_thickness=min_thk,
            pinch=pinch,
            length_units=self.length_units,
            name=name,
            attach=attach,
            reconcile=rec_on,
            which=which,
            min_sep=min_sep,
            method=method,
            refresh=refresh,
        )


__all__ = [
    "LayerStack",
    "LayerBuildResult",
    "LayerQCReport",
    "modflow_surfaces",
    "Raster",
    "Flat",
    "Contours",
    "Points",
    "Array",
    "Isopach",
    "Min",
    "Max",
    "Clamp",
    "Where",
]
