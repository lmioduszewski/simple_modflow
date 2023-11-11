import plotly.graph_objs as go
import plotly.subplots as subplots

class MfBaseFigure(go.Figure):
    def __init__(self):
        super().__init__()

        # FIG TEMPLATE
        self._config = {
            'scrollZoom': True,
            'displaylogo': False
        }
        self._modebar = go.layout.Modebar(
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
        self._xaxis_template = go.layout.XAxis(gridcolor='lightgray', griddash='dot')
        self._yaxis_template = go.layout.YAxis(gridcolor='lightgray', griddash='dot')
        self._font_template = go.layout.Font()

        self._template_layout = go.Layout(
            dragmode='pan',
            modebar=self._modebar,
            paper_bgcolor='white',
            plot_bgcolor='white',
            xaxis=self._xaxis_template,
            yaxis=self._yaxis_template,
        )
        self._template = go.layout.Template(layout=self._template_layout)
        self.layout = go.Layout(template=self._template)

    def show(self, *args, **kwargs):
        super().show(config=self._config, renderer='browser', *args, **kwargs)

