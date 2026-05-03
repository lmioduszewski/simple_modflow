from __future__ import annotations
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from .voronoiplus import VoronoiGridPlus as Vor
    from simple_modflow.modflow.utils.surfaces import InterpolatedSurface
    from simple_modflow.modflow.mf6.sfr import SFR

import flopy
from pathlib import Path
from simple_modflow.modflow.mf6.sfr import SFR
from functools import cached_property
from simple_modflow.modflow.mf6.simulation.accessors import (
    build_choro,
    build_xsection,
    get_all_heads,
    get_budget,
    get_budget_cumulative,
    get_budget_incremental,
    get_hds,
    get_inputs,
    get_kstpkper,
    get_lak_output,
    get_sfr_output,
    get_surface,
    get_uzf_output,
)
from simple_modflow.modflow.mf6.simulation.discretization import (
    DisuGrid,
    DisvGrid,
    TemporalDiscretization,
)
from simple_modflow.modflow.mf6.simulation.indexing import (
    build_idomain,
    build_model_times,
    build_ncpl_arr,
    build_node_to_lni,
    build_offsets,
    coerce_per_dates,
)
from simple_modflow.modflow.mf6.simulation.packages import (
    CHD,
    GHB,
    LAK,
    UZF,
    Drains,
    InitialConditions,
    KFlow,
    OutputControl,
    Recharge,
    Storage,
)
from simple_modflow.modflow.mf6.simulation.runtime import run_simulation as run_model_simulation


class SimulationBase:

    def __init__(
            self,
            name: str = 'mf6_model',
            mf_folder_path: Path = Path().home().joinpath('mf6'),
            nper: int = 1,
            vor: Vor = None,
            per_dates: list | pd.DatetimeIndex = None,
            idomain_path: Path = None,
            newton: bool = True,
            complexity: str = 'MODERATE',
    ):
        """
        Initializes a MODFLOW 6 model with specified configurations. This class
        sets up the environment for running a simulation with flopy and provides
        options for specifying simulation parameters, workspace paths, and solver characteristics.

        :param name: Name of the MODFLOW 6 model.
        :type name: str
        :param mf_folder_path: Path to the folder where model files should be stored.
        :type mf_folder_path: Path
        :param nper: Number of stress periods in the model simulation.
        :type nper: int
        :param vor: Optional argument for additional configurations. The type should be defined elsewhere in the code.
        :type vor: Vor
        :param per_dates: List or DatetimeIndex defining stress period start and end dates.
        :type per_dates: list | pd.DatetimeIndex
        :param idomain_path: Path to the idomain configuration file.
        :type idomain_path: Path
        :param newton: Boolean indicating whether Newton's solver options should be enabled.
        :type newton: bool
        :param complexity: Specifies the complexity level of the IMS (Iterative Model Solver) package. Possible values are 'SIMPLE', 'MODERATE', or 'COMPLEX'.
        :type complexity: str
        """
        self.name = name
        self.vor = vor
        self.nlay = None
        self.nper = nper
        self.num_steps = None
        self.per_len = None
        self._obs = None
        self._times = None
        self.model_output_folder_path = mf_folder_path.joinpath(f'{name}')
        self._master_celld = {}
        self._hds = None
        self._idomain_path = None
        self._bud = None
        self._sfr_input = None
        self._per_dates = None
        self._idomain_gdf = None
        self._inactive_cells = None
        self._node_to_lni = None
        self._lni_to_node = None
        self._kstpkper = None
        self._idomain = None

        self.per_dates = per_dates
        self.idomain_path = idomain_path

        if newton:
            newtonoptions = 'under_relaxation'
        elif not newton:
            newtonoptions = None


        self.sim = flopy.mf6.MFSimulation(
            sim_name=self.name,
            exe_name="mf6",
            version="mf6",
            sim_ws=self.model_output_folder_path,
        )
        self.gwf = flopy.mf6.ModflowGwf(
            self.sim,
            modelname=self.name,
            model_nam_file=f"{self.name}.nam",
            print_flows=False,
            newtonoptions=newtonoptions,
            save_flows=True
        )
        """self.ims = flopy.mf6.modflow.mfims.ModflowIms(
            self.sim,
            print_option='ALL',
            csv_inner_output_filerecord='ims_inner_convergence.csv',
            csv_outer_output_filerecord='ims_outer_convergence.csv',
            pname="ims",
            complexity="COMPLEX",
            under_relaxation="DBD",
            under_relaxation_theta=0.72,
            under_relaxation_kappa=0.1,
            # under_relaxation_gamma=0.2,
            under_relaxation_momentum=0.001,
            backtracking_number=20,
            backtracking_tolerance=1.1,
            backtracking_reduction_factor=0.2,
            backtracking_residual_limit=100,
            outer_maximum=300,
            inner_maximum=300,
            outer_dvclose=0.01,
            inner_dvclose=0.001,
            # rcloserecord=[0.5, 'strict'],
            relaxation_factor=0.97,
            linear_acceleration='BICGSTAB',
        )"""

        self.ims = flopy.mf6.modflow.mfims.ModflowIms(
            self.sim,
            print_option="SUMMARY",
            csv_inner_output_filerecord="ims_inner_convergence.csv",
            csv_outer_output_filerecord="ims_outer_convergence.csv",
            pname="ims_gwf",
            complexity=complexity,

            # Nonlinear controls
            under_relaxation="DBD",
            under_relaxation_theta=0.7,  # 0.7–0.8 typical
            under_relaxation_kappa=0.20,
            under_relaxation_momentum=0.001,  # disable momentum for stability
            backtracking_number=50,
            backtracking_tolerance=1.1,
            backtracking_reduction_factor=0.2,
            backtracking_residual_limit=10,

            # Iteration budgets
            outer_maximum=800,
            inner_maximum=300,

            # Convergence criteria
            outer_dvclose=1e-5,
            inner_dvclose=1e-6,
            # rcloserecord=[1e-4, 'relative_rclose'],

            # Linear solver + conditioning
            linear_acceleration="BICGSTAB",
            scaling_method="L2NORM",
            reordering_method="RCM",
            preconditioner_levels=5,
            preconditioner_drop_tolerance=0,

            relaxation_factor=0.97,
        )

        self.sim.register_ims_package(self.ims, [self.name])

    """def obs_wells(self, well_dict: dict, output_filename: str = 'well_obs.csv'):
        
        self._obs = flopy.mf6.ModflowUtlobs(
            self.gwf,
            continuous={output_filename:
                    [
                        ((k, 'HEAD', (0, v[0])) for k, v in well_dict.items())
                    ]                
            }
        )"""

    @property
    def inactive_cells(self):
        """returns a list of inactive cells in the model"""
        if self._inactive_cells is None:

            icells = self.vor.get_vor_cells_as_series(self.idomain_gdf.geometry).to_list()
            self._inactive_cells = icells

        return self._inactive_cells

    @property
    def idomain(self):
        """
        Property method that retrieves or computes the idomain array for the Voronoi grid. The idomain
        represents the active/inactive cell status for each layer in the grid. This property constructs
        the idomain by reading a GeoPackage or shapefile, assigning Voronoi cells to respective layers,
        and determining inactive cells for each layer.

        :rtype: list[list[int]]
        :return: The idomain as a list of lists, where each sublist represents a layer and contains the
                 activity status (1 for active, 0 for inactive) for each cell in that layer.
        """
        if self._idomain is None:
            self._idomain, self._idomain_gdf = build_idomain(self.vor, self.idomain_path)

        return self._idomain

    @property
    def idomain_path(self):
        """returns the path to the idomain shapefile/geopackage"""
        return self._idomain_path

    @idomain_path.setter
    def idomain_path(self, idomain_path):
        if idomain_path is not None:
            assert isinstance(idomain_path, Path), 'idomain_path must be a Path object'
            self._idomain_path = idomain_path

    """@property
    def idom_vor(self):
        return self.vor.gdf_vorPolys.loc[
            self.modelgrid.idomain.T.flatten() == 1]"""

    @property
    def per_dates(self):
        return self._per_dates

    @per_dates.setter
    def per_dates(self, per_dates):
        self._per_dates = coerce_per_dates(per_dates)

    @property
    def times(self):
        """dict where keys are stress periods and time steps (tuple with
        time step then stress period) and values are corresponding model times"""
        if self._times is None:
            self._times = build_model_times(self.gwf)
        return self._times

    @property
    def modelgrid(self) -> flopy.discretization.vertexgrid.VertexGrid:
        modelgrid: flopy.discretization.vertexgrid.VertexGrid = self.gwf.modelgrid
        return modelgrid

    @cached_property
    def node_to_lni(self) -> dict[int, tuple[int, int]]:
        """
        Computes a mapping from node indices to layer-node indices. This is an efficient
        way to associate a node in a flattened multi-layered structure with its corresponding
        layer and local node index. Works with DISV grids.

        :return: Dictionary mapping node indices to tuples, where each tuple contains the
            layer index and the node's position within that layer.
        :rtype: dict[int, tuple[int, int]]
        """
        return build_node_to_lni(self._ncpl_arr, self._offsets)

    @cached_property
    def _ncpl_arr(self) -> np.ndarray:
        return build_ncpl_arr(self.modelgrid)

    @cached_property
    def _offsets(self) -> np.ndarray:
        return build_offsets(self._ncpl_arr)

    def node_from_lni(self, layer: int, idx_in_layer: int) -> int:
        """
        Returns the global node index computed from the given layer index and the
        index within the specified layer. Works with DISV grids.

        :param layer: The layer number of the node in the structure.
        :type layer: int
        :param idx_in_layer: The index of the node within the given layer.
        :type idx_in_layer: int
        :return: The global node index calculated from the given layer and index.
        :rtype: int
        """
        return int(self._offsets[layer] + idx_in_layer)

    @property
    def cellids(self):
        return list(self.node_to_lni.keys())

    @property
    def sfr_input(self) -> SFR:
        if self._sfr_input is not None:
            assert isinstance(self._sfr_input, SFR), 'sfr input not SFR class'
            sfr_input: SFR = self._sfr_input
            return sfr_input
        else:
            return None

    @property
    def hds(self):
        return get_hds(self)

    @property
    def all_heads(self):
        return get_all_heads(self)

    @property
    def srf(self):
        """instaniates a ModelSurface object for the model"""
        return get_surface(self)

    def cor(
            self,
            kstpkper: tuple = None,
            per: int = None,
            layer: int = 0,
            type: str = 'hds',
            custom_hover: dict = None,
            custom_zs: list = None,
            zmin: float | int = None,
            zmax: float | int = None,
            zoom: int = 13,
            show_layer_elevs: bool = True,
            show_mounding: bool = False,
            hover_heads: bool = True,
            hover_ks: bool = False,
            locs=None,
            rch_scale=None,
            bgs=False,
            hillshade_path: Path = None,
            colorscale: str = None,
            logscale: bool = False,
            **kwargs
    ) -> Choro:
        """
        kwargs can be any allowable keyword arguments from the Choro class
        Class defining the basic choropleth plots generated from a modflow model.

        :param kstpkper: tuple of stress period and time step to plot
        :param per: can just provide stress period. appropriate kstpkper tuple will be determined, will throw an
        error if more than one valid kstpkper in the model output exists with the provided per
        :param layer: what layer to plot, a zero-index. 0 equals layer 1.
        :param type: type of choropleth to plot - options are 'hds', 'ks', 'input_rch', 'output_rch'
        :param custom_hover: custom dictionary of hover labels to use for choropleth.
        Must be same length as no. of cells in model.
        :param custom_zs: custom list of z values to use for choropleth. Must be same length as no. of cells in model.
        :param zmin: minimum z value to use for choropleth colorscale.
        :param zmax: maximum z value to use for choropleth colorscale.
        :param zoom: zoom level for choropleth map. Default is 13.
        :param show_layer_elevs: Default is True. To show elevations of all layers on hover
        :param show_mounding: if True, colorscale will be mounding over given Layer
        :param hover_heads: boolean to show heads on hover.
        :param hover_ks: boolean to show Kh on hover.
        :param locs: specify the path of a shapefile or geopackage with location points to show on choropleth
        :param rch_scale: value to scale z-values in choropleth. For example to convert units in recharge to another L/T
        :param bgs: if True, will take precedence, and will plot water levels in given Layer relative to top of model
        :param hillshade_path: path to hillshade raster.
        :param colorscale: colorscale to use for choropleth. Default is 'earth'.
        :param logscale: if True, will use logarithmic scale for colorscale.
        """
        return build_choro(
            self,
            kstpkper=kstpkper,
            per=per,
            layer=layer,
            type=type,
            custom_hover=custom_hover,
            custom_zs=custom_zs,
            zmin=zmin,
            zmax=zmax,
            zoom=zoom,
            show_layer_elevs=show_layer_elevs,
            show_mounding=show_mounding,
            hover_heads=hover_heads,
            hover_ks=hover_ks,
            locs=locs,
            rch_scale=rch_scale,
            bgs=bgs,
            hillshade_path=hillshade_path,
            colorscale=colorscale,
            logscale=logscale,
            **kwargs
        )

    def xs(
            self,
            per: int = None,
            kstpkper: tuple = None,
            layer: int = 0,
            cells: int | list[int] = None,
            x_or_y: str = None,
            spacing: int = 10,
            num_points: int = 100,
            extrapolate_beyond_section_ends: bool = False,
            interpolate: bool = False,
            use_rbf: bool = False,
            show_model_top = True,
            show_model_btm = False,
    ):
        """
        Use to plot a cross-section of heads through a model. Can be used to create an animation
        of head changes for all stress periods. The cross-section line can be defined by providing
        one cell (the 'cells' parameter) or as two ends by providing two cells to the 'cells'
        parameter. If just one cell is given, 'x_or_y' parameter defines whether the cross-section
        is vertical (along 'y' axis) or horizontal (along 'x' axis).

        Examples:

            Show an animated cross-section of all stress periods:

                XSection(model, cells=[1653, 651, 1241]).ani.show() ...OR...
                XSection(model, cells=69, layer=2, x_or_y='y').ani.show()

            Show just a cross-section of one stress period, no animation:

                XSection(model, cells=[1653, 651, 1241], kstpkper=(9, 50)).show()

        :param model: model (SimulationBase object) instance
        :param per: stress period number, O-based index; will take presedence over kstpkper if provided
        :param kstpkper: defaults to the first model stress period if not provided
        :param layer: defaults to 0
        :param cells: defines cross-section location. Can provide any number of cells
        :param x_or_y: only used if one cell is given, defines whether
        the cross-section is vertical (along 'y' axis) or horizontal (along 'x' axis).
        :param spacing: x distance between points on the plot
        :param num_points: number of points in the cross-section plot
        :param extrapolate_beyond_section_ends: not implemented
        :param surf_type: can be hds (default) or lyr (for model layers)

        """
        return build_xsection(
            self,
            per=per,
            kstpkper=kstpkper,
            layer=layer,
            cells=cells,
            x_or_y=x_or_y,
            spacing=spacing,
            num_points=num_points,
            extrapolate_beyond_section_ends=extrapolate_beyond_section_ends,
            interpolate=interpolate,
            use_rbf=use_rbf,
            show_model_top=show_model_top,
            show_model_btm=show_model_btm,
        )

    @property
    def inputs(self):
        return get_inputs(self)

    @property
    def lak(self):
        return get_lak_output(self)

    @property
    def uzf(self):
        return get_uzf_output(self)

    @property
    def sfr(self):
        return get_sfr_output(self)

    def bud(self, package: str = None):
        return get_budget(self, package)

    @property
    def budget_cumulative(self):
        return get_budget_cumulative(self)

    @property
    def budget_incremental(self):
        return get_budget_incremental(self)

    @property
    def master_celld(self):
        return self._master_celld

    @property
    def kstpkper(self):
        return get_kstpkper(self)

    def run_simulation(self):
        return run_model_simulation(self)



    """def plot_hds(self, kstpkper, zoom=13, plot_mounding=False, layer=0, zmin=None, zmax=None):
        layer_nums = self.vor.gdf_topbtm.columns[2:].to_list()
        hover = {"": ["" for cell in range(self.vor.ncpl)]}
        hover.update({
            f'Top of Model': self.vor.gdf_topbtm.loc[:, 0].to_list()
        })
        hover.update({
            f'Layer {layer} Bottom': self.vor.gdf_topbtm.loc[:, layer_num].to_list() for layer_num in layer_nums
        })

        self.hds.plot_choropleth(
            kstpkper=kstpkper,
            zoom=zoom,
            plot_mounding=plot_mounding,
            custom_hover=hover,
            layer=layer,
            zmax=zmax,
            zmin=zmin
        )"""
