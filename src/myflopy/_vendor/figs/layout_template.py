import plotly.graph_objects as go
from pathlib import Path
import base64
from ._layout_presets import LayoutPresets


# AESI logo, embedded at sync time (see scripts/sync_vendored_figs.py)
from ._logo_data import AESI_LOGO_B64 as aesi_logo

### MODEBAR ###
modebar = go.layout.Modebar(
    add=[
        'togglespikelines',
        'hovercompare',
        'togglehover',
        'drawline',
        'drawopenpath',
        'drawclosedpath',
        'drawcircle',
        'drawrect',
        'eraseshape',
        'hoverclosest',
    ]
)

### MARGIN ###
margin = go.layout.Margin(autoexpand=True, b=120, t=20, l=50, r=0)

### X-AXIS ###
xaxis_template = go.layout.XAxis(
    showticklabels=True,
    gridcolor='lightgray',
    griddash='dot',
    spikethickness=1,
    spikemode='across',
    ticks='outside'
)

### Y-AXIS ###
yaxis_template = go.layout.YAxis(
    gridcolor='lightgray',
    griddash='dot',
    spikethickness=1,
    spikemode='across',
    ticks='outside'
)

### DEFAULT FONT ###
font_template = go.layout.Font(
    family='Calibri, Open Sans',
    size=12,
    color='black'
)
# add a plot border shape
plot_border = go.layout.Shape(
    xref="paper", yref="paper",
    type='rect',
    x0=0, y0=0,
    x1=1, y1=1,
    line=dict(color="black", width=1),
    layer="above",
    visible=True,
    name='border'
)

### MASTER LAYOUT OBJECT ###
template_layout = go.Layout(
    dragmode='pan',
    annotations=[],
    shapes=[plot_border],
    modebar=modebar,
    paper_bgcolor='white',
    font=font_template,
    plot_bgcolor='white',
    xaxis=xaxis_template,
    #xaxis2=xaxis_template,
    yaxis=yaxis_template,
    #yaxis2=yaxis_template,
    margin=margin,
)

template = go.layout.Template(layout=template_layout)
layout = go.Layout(template=template)

# Scattergl template
scattergl_template = go.Scattergl(
    marker=go.scattergl.Marker(
        size=4,
        symbol='circle',
        line=go.scattergl.marker.Line(  # the line edge of the trace marker
            color='black',
            width=0.5,
        ),
    ),
)

# Titles
titles = {
            'water levels': 'Water Level Data',
            'map': 'Monitoring Locations',
            'rain': '',
            'flow': 'Stream Flow Data'
            }
