import numpy as np
import pandas as pd
from pathlib import Path
import plotly.express as px
import plotly.graph_objs as go

# Define paths to excel files, import to DataFrames, and set indexes and column names
pump_times_filename = Path.home() / 'Documents' / 'Tehaleh' / 'CDF pump run times.xlsx'
df_pump_times = pd.read_excel(pump_times_filename)
df_columns = pd.Index(
    ['Date', '1A', '1B', '2A', '2B', '3A', '3B',
     '4A', '4B', '5A', '5B']
)
df_pump_times.columns = df_columns
df_pump_times = df_pump_times.set_index('Date')
df_rates = pd.read_excel(pump_times_filename, 1)
df_rates = df_rates.set_index('pump')


def calc_flows_from_pump_run_times(
        df_times: pd.DataFrame,
        df_rates: pd.DataFrame,
):
    df = df_times.copy(deep=True)
    for lobe in range(1, 6):
        lobeA = f'{lobe}A'
        lobeB = f'{lobe}B'
        more_zero_A = (df[lobeA] > 0)
        more_zero_B = (df[lobeB] > 0)
        is_zero_A = (df[lobeA] == 0)
        is_zero_B = (df[lobeB] == 0)
        more_zero_both = more_zero_A & more_zero_B
        diff_b_a = df[lobeB] - df[lobeA]
        just_A = more_zero_A & is_zero_B
        just_B = more_zero_B & is_zero_A

        flow_cases = [
            diff_b_a.loc[df.loc[more_zero_both, lobeA:lobeB].loc[diff_b_a > 0].index] * df_rates.loc[lobeB][0],
            diff_b_a.loc[df.loc[more_zero_both, lobeA:lobeB].loc[diff_b_a < 0].index] * df_rates.loc[lobeA][0] * -1,
            df.loc[more_zero_both, lobeA:lobeB].loc[diff_b_a > 0].loc[:, lobeA] * df_rates.loc['2_pumps'][0],
            df.loc[more_zero_both, lobeA:lobeB].loc[diff_b_a < 0].loc[:, lobeB] * df_rates.loc['2_pumps'][0],
            df.loc[just_B, lobeB] * df_rates.loc[lobeB][0],
            df.loc[just_A, lobeA] * df_rates.loc[lobeA][0]
        ]
        flow_cases[0].name = lobeB
        flow_cases[1].name = lobeA

        for flows in flow_cases:
            df.update(flows)
        return df


def add_pumpA_pumpB_flows(df: pd.DataFrame):
    df_to_concat = []
    for lobe in range(1, 6):
        lobeA = f'{lobe}A'
        lobeB = f'{lobe}B'
        df1 = df[lobeA] + df[lobeB]
        df1.name = lobe
        df_to_concat += [df1]
    df_flows = pd.concat(df_to_concat, axis=1)
    return df_flows


df_flows = calc_flows_from_pump_run_times(df_pump_times, df_rates)
df_flows = add_pumpA_pumpB_flows(df_flows)
df_percentages = df_flows.div(df_flows.sum(axis=1), axis=0)
df_percentages.columns = ('1_per', '2_per', '3_per', '4_per', '5_per')
df_final = pd.concat([df_flows, df_percentages], axis=1)

fig = go.Figure()
fig.update_layout(
    yaxis2=dict(
        title="Percentage of total flow",
        overlaying="y",
        side="right"
    ),
    yaxis_title='Flow (gpd)',
    dragmode='pan'
)
for col in df_final:
    if 'per' in str(df_final[col].name):
        y_axis = 'y2'
        line_mode = 'dash'
    else:
        y_axis = 'y1'
        line_mode = 'solid'
    fig.add_scattergl(
        x=df_final[col].index,
        y=df_final[col].tolist(),
        yaxis=y_axis,
        name=df_final[col].name,
        line_dash=line_mode
    )

    fig.show(renderer='browser', config={'scrollZoom': True})
    #fig.write_html('run_times.html', config={'scrollZoom': True})
