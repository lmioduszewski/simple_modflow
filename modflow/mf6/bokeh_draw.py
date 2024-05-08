from bokeh.models import ColumnDataSource, FreehandDrawTool, Div, CustomJS, BoxEditTool, Slider
from bokeh.plotting import figure, output_file, show
from bokeh.layouts import column

# Output to an HTML file
output_file("freehand_draw.html")

# Create a new plot with a title
p = figure(width=400, height=400, title="Draw on this!")

# Create a data source for the freehand draw tool with initial empty lists.
source = ColumnDataSource(data={'x': [[1,2,3]], 'y': [[1,2,3]]})
box_source = ColumnDataSource(data={})

# Use a multi_line glyph for drawing.
line_glyph = p.multi_line('x', 'y', source=source, line_width=2, alpha=0.8)
# Rect glyph for box drawing and editing
box_glyph = p.rect(source=box_source,angle_units='deg',angle=0)

slider = Slider(start=0,end=360,value=0,step=1)

# Add the freehand draw tool to the plot and activate it.
draw_tool = FreehandDrawTool(renderers=[line_glyph], num_objects=1)
box_tool = BoxEditTool(renderers=[box_glyph],empty_value=1)
p.add_tools(draw_tool,box_tool)
p.toolbar.active_drag = box_tool

# Create a Div to show the drawn data as JSON.
div = Div(text="JSON representation of your shape will appear here.", width=400, height=100)

# Use a JavaScript callback to update the Div when the data changes.
callback = CustomJS(args=dict(source=box_source, div=div), code="""
    const data = source.data;
    div.text = JSON.stringify(data, undefined, 2);
""")
box_source.js_on_change('data', callback)
slider.js_link('value',box_glyph,'angle')

# Arrange the plot and the Div in a layout.
layout = column(p, slider, div)

# Show the layout.
show(layout)
