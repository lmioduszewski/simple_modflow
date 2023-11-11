import numpy as np
import pandas as pd
import matplotlib as mpl
from voronoiplus import VoronoiGridPlus as vgp
import flopy.utils.binaryfile as bf
import math
from pathlib import Path
import flopy
import shapely as shp
import plotly.graph_objs as go
import plotly

idxx = pd.IndexSlice  # for easy index slicing in a MultiIndex DataFrame


class BudgetPlus(bf.CellBudgetFile):

    def __init__(
            self,
            cbb_path: bf.CellBudgetFile,
            vor: vgp = None
    ):

        super().__init__(cbb_path)
        self.cbb = bf.CellBudgetFile(cbb_path)
        self.nperstp = len(self.get_kstpkper())
        self.vor = vor
        self.mdx = self.mdxVor(vor=self.vor)
        self.bud_recs = self.get_buds_recs()
        self.df_jaFlows = self.get_jaFlows()
        self.fig_layout = go.Layout(  # default fig layout
            dragmode='pan'
        )
        self.fig_config = {  # default fig config
            'scrollZoom': True
        }

    def mdxVor(self,
               cells: list = None,
               vor: vgp = None,
               ) -> pd.MultiIndex:
        """Generates a Pandas MultiIndex object, based number of cells
                and number of stress periods

                Args:
                cells (list, optional): Number of cells in the grid. Defaults to None.
                vor (Vor, optional): VoronoiGridPlus object to find number of cells.
                        This will take precedence if both cells and vor are defined.
                        Defaults to None.

                Returns:
                pd.MultiIndex: Pandas MultiIndex, with levels corresponding to stress 
                                        periods and all cells in the Voronoi grid
                """

        splist = list(range(self.nperstp))
        mdx = None
        if self.vor != None:
            vor = self.vor
        if cells != None:
            celllist = list(range(cells))
            mdx = pd.MultiIndex.from_product(
                iterables=[
                    splist,
                    celllist
                ],
                names=['sp', 'cells']
            )
        if vor != None:
            vorlist = list(range(len(vor.iverts)))
            mdx = pd.MultiIndex.from_product(
                iterables=[
                    splist,
                    vorlist
                ],
                names=['sp', 'cells']
            )
        return mdx

    def get_buds_recs(self):
        """Returns a dict of all budget records

                Returns:
                    dict: dict of budget records
                """

        # make list of the names of unique records in the budget file
        num_records = len(self.get_unique_record_names())
        record_list = []
        for record in range(num_records):
            thisrecord = str(self.get_unique_record_names()[record]).split()[1][0:-1]
            record_list.append(thisrecord)
        """make dict of DataFrames for each budget record"""
        dictBudRecords = {}
        for record in record_list:
            dictBudRecords[record] = {}
            thisRecordDict = dictBudRecords[record]
            for stpper in self.get_kstpkper():
                thisstp = stpper[0]
                thisper = stpper[1]
                thisrec = self.get_data(
                    kstpkper=(
                        thisstp,
                        thisper
                    ),
                    text=record,
                )
                """make Pandas DataFrame or Series of this record"""
                try:
                    """DataFrame will work for most records, except as below"""
                    thisrec = pd.DataFrame(thisrec[0])
                except:
                    """Series only works for storage and face-to-face cell flow"""
                    thisrec = pd.DataFrame(pd.Series(thisrec[0][0][0]))
                thisRecordDict[f'{record}-per{thisper}-stp{thisstp}'] = thisrec
        return dictBudRecords

    def get_jaFlows(
            self,
            record: str = None
    ) -> pd.DataFrame:
        """Get cell-to-cell flow, adjacent cell flows (jaFlows)
                
                Args:
                    record (str): name of record to get flows from

                Returns:
                    pd.DataFrame: DataFrame of adjacent cell flows
                """

        if record == None:
            record = list(self.get_buds_recs()['FLOW-JA-FACE'].keys())[0]

        thisRec = self.get_buds_recs()['FLOW-JA-FACE'][record]
        jaSlice = 0
        budJa = self.vor.ja
        jaFlowDict = {}
        jaFlowList = []

        for cell, num_iac in enumerate(self.vor.iac):
            theseFlows = list(thisRec.iloc[jaSlice:jaSlice + num_iac][0])
            theseCells = budJa[jaSlice:jaSlice + num_iac]
            jaSlice += num_iac
            cells = [cell for i in range(num_iac)]
            theseJaFlows = list(zip(cells, theseCells, theseFlows))
            thisCellDict = {
                cell: theseJaFlows,
            }
            jaFlowDict.update(thisCellDict)
            jaFlowList += theseJaFlows

        df_jaFlows = pd.DataFrame().from_records(jaFlowList)
        df_jaFlows.columns = ['cell1', 'cell2', 'flows']
        df_jaFlows.index.names = [record]

        return df_jaFlows

    def get_jaFlowVectors(
            self,
            df_jaFlows: pd.DataFrame = None,
            quiver_scaling_factor: int = 1,
    ) -> pd.DataFrame:

        x_list = []
        y_list = []
        u_list = []
        v_list = []
        norm_list = []
        u_new_list = []
        v_new_list = []
        x_new_list = []
        y_new_list = []

        for vector in range(len(df_jaFlows)):
            """get the flow value for this vector"""
            thisFlow = df_jaFlows.iloc[vector].loc['flows']

            """calculate shared face length for this vector"""
            cell1 = int(df_jaFlows.iloc[vector].loc['cell1'])
            cell2 = int(df_jaFlows.iloc[vector].loc['cell2'])
            cell1 = self.vor.gdf_vorPolys['geometry'].iloc[cell1]
            cell2 = self.vor.gdf_vorPolys['geometry'].iloc[cell2]
            shared_face_length = self.vor.shared_face_length(cell1, cell2)

            thisvector_idx1 = int(df_jaFlows.iloc[vector].loc['cell1'])
            thisvector_idx2 = int(df_jaFlows.iloc[vector].loc['cell2'])
            thisvector_x1 = self.vor.centroids_x[thisvector_idx1]
            thisvector_y1 = self.vor.centroids_y[thisvector_idx1]
            thisvector_x2 = self.vor.centroids_x[thisvector_idx2]
            thisvector_y2 = self.vor.centroids_y[thisvector_idx2]
            this_u = thisvector_x2 - thisvector_x1
            this_v = thisvector_y2 - thisvector_y1
            this_norm = math.sqrt((this_u ** 2) + (this_v ** 2))
            this_x_new, this_y_new, this_u_new, this_v_new = self.calculate_new_vector_with_same_center(
                x=thisvector_x1,
                y=thisvector_y1,
                u=this_u,
                v=this_v,
                new_magnitude=abs(thisFlow / (quiver_scaling_factor * shared_face_length))  # normalize by face length
            )

            x_list += [thisvector_x1]
            y_list += [thisvector_y1]
            u_list += [this_u]
            v_list += [this_v]
            norm_list += [this_norm]
            x_new_list += [this_x_new]
            y_new_list += [this_y_new]
            u_new_list += [this_u_new]
            v_new_list += [this_v_new]

        df_jaFlows['x'] = x_list
        df_jaFlows['y'] = y_list
        df_jaFlows['u'] = u_list
        df_jaFlows['v'] = v_list
        df_jaFlows['norm'] = norm_list
        df_jaFlows['x_new'] = x_new_list
        df_jaFlows['y_new'] = y_new_list
        df_jaFlows['u_new'] = u_new_list
        df_jaFlows['v_new'] = v_new_list

        return df_jaFlows

    def plot_quiver(
            self,
            df_jaFlows: pd.DataFrame,
            flow_filter: int = -100,
            quiver_scaling_factor: int = 1
    ) -> go.Figure:

        if df_jaFlows.empty:
            df_jaFlows = self.df_jaFlows

        df_jaFlows = self.get_jaFlowVectors(
            df_jaFlows=df_jaFlows,
            quiver_scaling_factor=quiver_scaling_factor)

        filtered_jaFlows = df_jaFlows[df_jaFlows.loc[:, 'flows'] < flow_filter]
        fig_vorquiv = self.vor.plot2d()
        quiver_fig = plotly.figure_factory.create_quiver(
            x=filtered_jaFlows['x_new'],
            y=filtered_jaFlows['y_new'],
            u=filtered_jaFlows['u_new'],
            v=filtered_jaFlows['v_new'],
            scale=1,
            hoverinfo='text',
            hovertext=filtered_jaFlows['flows']

        )
        quiver_gl = go.Scattergl(
            x=quiver_fig.data[0].x,
            y=quiver_fig.data[0].y,
            mode=quiver_fig.data[0].mode,
            hoverinfo='text',
            hovertext=filtered_jaFlows['flows'],
            line_width=1,
        )
        quiver_fig = go.Figure(
            data=[
                fig_vorquiv.data[0],
                quiver_gl
            ]
        )
        quiver_fig.layout = self.fig_layout

        return quiver_fig.show(renderer='browser', config=self.fig_config)

    def get_rch_bud(
            self,
            kstpkper_just_one: tuple = None
    ):
        """Method to get a Dataframe of the recharge budget for the model.
                Recharge is given for every cell of the model for every time step
                and stress period saved by MODFLOW. Can define what stress period
                and time step for which to get the recharge budget, called
                (kstpkper_just_one) to differentiate from kstpkper in this 
                function, which includes all time steps and stress periods."""

        cell_budget = self.cbb
        num_rch_cells = len(self.vor.gdf_vorPolys)
        kstpkper = self.get_kstpkper()

        """set generic MultiIndex for all stress periods and all cells"""
        multiIdxRch = pd.MultiIndex.from_product(
            iterables=[
                kstpkper,
                list(range(num_rch_cells))
            ],
            names=['kstpkper', 'cell']
        )
        """Set up a MultiIndex DataFrame to hold the recharge budget
                for all cells and stress periods"""
        df_rch_budget = pd.DataFrame(
            index=multiIdxRch,
            columns=['node', 'node2', 'q']
        )
        """get data for each saved time step and stress period"""
        for stp_per in kstpkper:
            spRch = cell_budget.get_data(
                kstpkper=stp_per,
                text='RCH'
            )
            spRch = pd.DataFrame(spRch[0])
            """copy and paste this stress period data to the MultiIndex DataFrame"""
            df_rch_budget.loc[idxx[stp_per, :len(spRch) - 1], idxx[:]] = spRch.values

        if kstpkper_just_one:
            df_rch_budget = self.get_rch_bud().loc[idxx[kstpkper_just_one, :], :]

        return df_rch_budget

    def get_rch_bud_for_cells(
            self,
            cell_list: list = None,
            kstpkper: tuple = None,
    ):
        """Method to get a DataFrame of the recharge budget of a list
                of model cells passed to the function for all time step and 
                stress periods saved by the MODFLOW model. Can also define desired
                stress period and time step (kstpkper).

                Args:
                    cell_list (list, optional): A list of cells that defines the
                    region for which to get the recharge budget. Defaults to [], an
                    empty list.
                    
                    kstpkper (tuple, optional): A tuple of the stress period and time
                    step for which to get a recharge budget. Defaults to (), an empty
                    tuple.

                Returns:
                    pd.DataFrame: Returns the Pandas DataFrame or Series (if defining kstpkper).
                """

        dfRechargeBudget = self.get_rch_bud()
        dfRechargeBudget = dfRechargeBudget.loc[idxx[:, cell_list], 'q']

        if kstpkper:
            dfRechargeBudget = dfRechargeBudget.loc[kstpkper]

        return dfRechargeBudget

    @staticmethod
    def toGal(cft):
        """converter from cubic feet to gallons

                Args:
                    cft (float): input cubic feet

                Returns:
                    float: returns gallons
                """
        return cft * 7.48052

    @staticmethod
    def calculate_new_vector_with_same_center(x, y, u, v, new_magnitude):
        """calculate the new vector x,y coordintes and u,v (of x,y) components
                based on the new desired magnitude. New vector will have the same center
                as the input vector, but with the new magnitude

                Args:
                    x (float): starting x-coordinate for the input vector
                    y (float): starting y-coordinate for the input vector
                    u (float): x-component of the vector direction and magnitude
                    v (float): y-component of the vector direction and magnitude
                    new_magnitude (float): magnitude of vector to calculate

                Returns:
                    floats: returns x, y, u, and v of new vector
                """

        # Compute the end point of the original vector
        end_x = x + u
        end_y = y + v

        # Compute the midpoint (center) of the original vector
        center_x = (x + end_x) / 2.0
        center_y = (y + end_y) / 2.0

        # Compute the direction of the original vector
        direction_x = u
        direction_y = v

        # Compute the norm (magnitude) of the original vector
        original_magnitude = math.sqrt(direction_x ** 2 + direction_y ** 2)

        if original_magnitude == 0:
            return None, None, None, None

        # Calculate the ratio of the new magnitude to the original magnitude
        ratio = new_magnitude / original_magnitude

        # Compute the new direction
        new_direction_x = ratio * direction_x
        new_direction_y = ratio * direction_y

        # Compute the new vector endpoints
        new_u = new_direction_x
        new_v = new_direction_y
        new_x = center_x - new_direction_x / 2.0
        new_y = center_y - new_direction_y / 2.0

        return new_x, new_y, new_u, new_v

    @staticmethod
    def calculate_shared_face_center(Polygon1: shp.Polygon, Polygon2: shp.Polygon):
        # calculate the intersection of the two polygons
        shared_face = Polygon1.intersection(Polygon2)

        # if the intersection is a LineString (which should be the case for Voronoi cells)
        if shared_face.geom_type == 'LineString':
            # calculate the center of the shared line segment
            center_x = (shared_face.coords[0][0] + shared_face.coords[1][0]) / 2
            center_y = (shared_face.coords[0][1] + shared_face.coords[1][1]) / 2

            return (center_x, center_y)

        else:
            print("Shared face is not a line segment.")
            return None
