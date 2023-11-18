import plotly.graph_objs as go
import plotly.subplots as subplots
import plotly


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

    def show(self, renderer='browser', *args, **kwargs):
        super().show(config=self._config, renderer=renderer, *args, **kwargs)

    @staticmethod
    def create_hover(name_list: dict = None):
        """
        Function to return a hover template for a Plotly figure
        :return: hovertemplate
        """
        names = list(name_list.keys())
        lists = list(name_list.values())
        list_len = len(lists[0])
        custom_data = []

        for i in range(list_len):
            this_list = [param[i] for param in lists]
            custom_data.append(this_list)

        hover_template_list = [
            f'<b>{name}: </b>%{{customdata[{i}]}}<br>' for i, name in enumerate(names)
        ]
        hover_template_list.append('<extra></extra>')
        hover_template = ''.join(hover_template_list)

        return custom_data, hover_template

if __name__ == "__main__":

    d = {'model': [1, 2, 3], 'version': [4, 5, 6], 'user': [7, 8, 9]}
    fig = MfBaseFigure()
    hover = fig.create_hover(d)
    print(hover[0])
    print(hover[1])
