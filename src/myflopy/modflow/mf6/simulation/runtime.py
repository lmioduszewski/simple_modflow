from __future__ import annotations

import pickle

from myflopy._logging import get_logger

logger = get_logger(__name__)


def run_simulation(model):
    """Write, pickle, and run ``model``'s MF6 simulation, returning ``(success, output)``.

    Writes the input files, pickles the model to a ``.model`` file for reloading,
    then runs MF6 and reports success.
    """

    model.sim.write_simulation()

    model_file_path = model.model_output_folder_path / f'{model.name}.model'
    try:
        with open(model_file_path, 'wb') as file:
            pickle.dump(model, file)
        logger.info("saved the model object to %s", model_file_path)
    except Exception as exc:  # noqa: BLE001 - see below; the graph is unbounded
        # Deliberately broad, on an asymmetry worth stating. `pickle.dump(model)`
        # serializes an UNBOUNDED third-party object graph -- the flopy
        # simulation, the Voronoi grid, geopandas/shapely, plus whatever the
        # user hung on the model -- and any `__reduce__` in it may raise a type
        # of its own choosing (a ctypes handle, for one, raises ValueError). The
        # `.model` pickle is a convenience snapshot nothing in this repo reads,
        # while an escaped exception here would abort before
        # `run_simulation()` on the next line and kill the MF6 run itself.
        logger.warning(
            "could not save the .model snapshot to %s: %s", model_file_path, exc,
        )

    success, buff = model.sim.run_simulation(silent=False, report=True)
    print("\nSuccess is: ", success)
    return success, buff
