from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.mfsimbase import SimulationBase


def create_custom_modelgrid(model: SimulationBase, BaseGridClass):
    class ModelGrid(model.gwf.modelgrid):

        def __init__(self):
            # super().__init__()  # Initialize with inherited arguments

            self._node_to_lni = None
            self._lni_to_node = None
            self.model = model

        @property
        def node_to_lni(self) -> dict:
            """
            create a dict where keys are model nodes and values are corresponding layer node indices,
            i.e. layer specific index of node
            :return: dict
            """
            if self._node_to_lni is None:

                node_to_lni = {}
                for node in range(self.nnodes):
                    node_to_lni[node] = self.get_lni([node])[0]
                self._node_to_lni = node_to_lni

            return self._node_to_lni

        @property
        def lni_to_node(self) -> dict:
            """
            create a dict where keys are layer node indices and values are corresponding model nodes
            :return: dict
            """
            if self._lni_to_node is None:

                lni_to_node = {}
                for node in range(self.nnodes):
                    lni_to_node[self.get_lni([node])[0]] = node
                self._lni_to_node = lni_to_node

            return self._lni_to_node

    return ModelGrid
