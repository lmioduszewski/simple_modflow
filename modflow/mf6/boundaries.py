from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as vgp

class Boundaries:
    
    def __init__(
        self,
        vor: vgp = None,
        bound_type: str = None
    ):
        
        self.vor = vor
        self.bound_type = bound_type

    def get_drn_stress_period_data(
        self, 
        cells: list,
        bottom_addition: float = 0,
        conductance: float = 100,
        disMf: str = 'disu',
        bottoms: dict = None,
        ) -> list:
        """Returns a list of lists. Each nested list corresponds to the DRN package
        boundary data for a particular voronoi cell in the grid, which includes cell
        ID, elevation of drain, and conductance. Can be passed to the flopy DRN package.

        Args:
            cells (list): list of cell IDs in this drain
            bottoms (dict): dict of bottoms, keys are the cell indices, values are the bottom elevations
            bottom_addition (float): height above the bottom of cell for the drain. This is added to the bottom of cell elevation derived from the Voronoi grid object.
            conductance (float): conductance for this drain cell
            disMf (str, optional): Either 'disu' or 'disv' works, and refers to the MODFLOW6 discretization package being used. Defaults to 'disu'.

        Returns:
            list: List of lists that contain the data for this drain and can be passed to the flopy DRN package
        """
        
        if self.vor is None:
            return print("No voronoi grid defined")
            
        drn_values = []
        for cell in cells:
            if bottoms:
                thisdrn = [cell, (bottoms[cell] + bottom_addition), conductance]
            elif self.vor.gdf_topbtm is not None:
                thisdrn = [cell, (self.vor.gdf_topbtm.loc[cell, "bottom"] + bottom_addition), conductance]
            else:
                thisdrn = [cell, bottom_addition, conductance]
            if disMf == "disv":
                thisdrn = [0] + thisdrn  # add layer num for disv grid
            drn_values.append(thisdrn)
        return drn_values