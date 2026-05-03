from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
    from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor


class Animation:
    """
    Represents an animation for a simulation model.

    This class provides mechanisms to generate interactive animations using
    update menus and sliders. It is designed to integrate with simulation
    models and enable user interaction to play, pause, and navigate frames
    of the animation.

    :ivar model: The simulation model associated with the animation.
    :type model: SimulationBase
    """
    def __init__(self, model: SimulationBase):
        self.model = model

    @property
    def updatemenus(self):
        updatemenus = [{
            'type': 'buttons',
            'buttons': [
                {'args': [None, {'frame': {'duration': 125, 'redraw': False},
                                 'transition': {'duration': 0, 'easing': 'quad-in'},
                                 'fromcurrent': True,
                                 'mode': 'afterall'}],
                 'label': 'Play',
                 'method': 'animate'},
                {'args': [[None], {'mode': 'immediate',
                                   'frame': {'duration': 0, 'redraw': False}
                                   }],
                 'label': 'Pause',
                 'method': 'animate'}
            ]
        }]

        return updatemenus

    @property
    def sliders(self):
        sliders = [{
            'active': 0,
            'yanchor': 'top',
            'xanchor': 'left',
            'currentvalue': {
                'font': {'size': 20},
                'prefix': 'Frame:',
                'visible': True,
                'xanchor': 'right'
            },
            'transition': {'duration': 0,
                           'easing': 'linear'},
            'pad': {'b': 10, 't': 50},
            'len': 0.9,
            'x': 0.1,
            'y': 0,
            'steps': [{
                'args': [[f'{per}'],
                         {'frame': {'duration': 100, 'redraw': False},
                          'mode': 'immediate',
                          'fromcurrent': True,
                          'transition': {
                              'duration': 0,
                              'easing': 'linear'
                          }}],
                'label': f'{per}',
                'method': 'animate'} for per in self.model.kstpkper]}
        ]

        return sliders
