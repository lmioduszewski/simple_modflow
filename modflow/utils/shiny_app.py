from shiny.express import input, ui, render
from shinywidgets import render_plotly
from shiny import reactive

ui.input_selectize(
    "var", "Select variable",
    choices=["bill_length_mm", "body_mass_g"]
)
ui.input_checkbox(
    'check', 'Baby',

)

@render_plotly
def hist():
    import plotly.express as px
    from palmerpenguins import load_penguins
    df = load_penguins()
    return px.histogram(df, x=input.var())

reactive.value()
