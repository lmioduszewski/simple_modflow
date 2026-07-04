import sys
import os
import shutil
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
import pandas as pd
import plotly
import plotly.express as px
import dash
import plotly.graph_objects as go
import mlti
from matplotlib.figure import Figure
from pathlib import Path
import flopy

def get_heads(folder_path=None, 
              stp=None, 
              per=None, 
              lyr=None, 
              headfile=None):
    # Read the binary head file and plot the results
    # We can use the existing Flopy HeadFile class because
    # the format of the headfile for MODFLOW 6 is the same
    # as for previous MODFLOW verions
        
    fname = os.path.join(folder_path, headfile)
       
    hds = flopy.utils.binaryfile.HeadFile(fname)
    hds_per = hds.get_data(kstpkper=(stp, per))
    
    if lyr is None:
        heads_to_return=hds_per
    else:
        heads_to_return=hds_per[lyr]
    return heads_to_return

def df_xsec(lyr=0, 
            headfile=None, 
            nrow=None, 
            ncol=None, 
            nper=None, 
            trans_tstp=None, 
            folder_path=None):
    """Creates a cross-section of y-values along center 
    column of a defined layer by reading the head file from MODFLOW"""
    
    df_rows = nrow * (nper - 1 ) * trans_tstp
    df_xsec = pd.DataFrame(np.ones((df_rows,4)), 
                           columns=["x-distance", 
                                    "y-head", 
                                    "time steps", 
                                    "stress period"]
                           )
    """for loop to cycle through every cell along a x-section 
    through the middle (ncol/2) of the model grid 
    for every time step 'i' in every stress period 'j'. """
    for j in np.arange(1, nper, 1): #starting at period 1, the first transient stress period
              
        for i in np.arange(0, trans_tstp, 1): #starting at timestep 0 of stress period j
            print(f"Stress period = {j} | Time step = {i}", end='\r')
            tstp_strt = (nrow * i) + ((j - 1) * 1000)
            tstp_end = tstp_strt + nrow

            y_heads = pd.DataFrame(get_heads(folder_path, 
                                             stp=i, 
                                             per=j, 
                                             lyr=lyr,
                                             headfile=headfile)
                                   ) #get heads for time step i in period j
            y_heads = y_heads[ncol/2] #get heads along center column of model grid

            df_xsec.iloc[tstp_strt:tstp_end, 0] = np.arange(0, nrow, 1) #set x-distance for time step i
            df_xsec.iloc[tstp_strt:tstp_end, 1] = y_heads.iloc[0:(nrow)] #set y-head for time step i
            df_xsec.iloc[tstp_strt:tstp_end, 2] = (i + ((j - 1) * trans_tstp)) #set current 'i' time step value
            df_xsec.iloc[tstp_strt:tstp_end, 3] = j #set current 'j' stress period value
    return df_xsec

def plot_heads(hds_plt: pd.DataFrame) -> go.Figure.show:
    """Plots heads from the MODFLOW6 output .hds file"""
    
    plyt = go.Figure(data=
                     [go.Surface(
                         z=hds_plt,
                         contours = {
                             "z": {
                                 "show": True, 
                                 "start": 0.5, 
                                 "end": 5, 
                                 "size": 0.5
                                 }
                             }
                         )
                      ]
                     )
    
    """Update the layout and general characteristics of the plotly go.Surface 3D plot."""
    plyt.update_layout(
        title='Mounding', 
        autosize=False,
        width=750, 
        height=750,
        margin=dict(l=65, r=50, b=65, t=90),
        )
    plyt.update_traces(contours_z=dict(show=True, usecolormap=True,
                                  highlightcolor="limegreen", project_z=True))
    plyt.update_scenes(zaxis_range=(0,20))
    
    """Make the plot visible"""
    plyt.show()
    