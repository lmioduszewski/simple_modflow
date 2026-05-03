from __future__ import annotations

import pickle


def run_simulation(model):
    model.sim.write_simulation()

    model_file_path = model.model_output_folder_path / f'{model.name}.model'
    try:
        with open(model_file_path, 'wb') as file:
            pickle.dump(model, file)
        print(f"\nSaved model object to .model file: {model_file_path}\n")
    except Exception as exc:
        print(f"\nError saving model object to .model file: {exc}\n")

    success, buff = model.sim.run_simulation(silent=False, report=True)
    print("\nSuccess is: ", success)
    return success, buff
