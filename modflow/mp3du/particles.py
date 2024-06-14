import os
import json
import subprocess
import flopy
from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
from pathlib import Path


class ParticleTrackingInput:
    def __init__(
            self,
            model: SimulationBase = None,
            writep3dgsf_path: Path = None,
            mp3du_path: Path = None,
            model_output_files: dict = None,
            output_path: Path = None,
            porosities_by_layer = None,
            particle_shp: Path = None
    ):
        self.model = model
        mp3du_default = Path(r"C:\Users\lukem\Python\Projects\simple_modflow\modflow\mp3du\bin\mp3du.exe")
        writegsf_default = Path(r"C:\Users\lukem\Python\Projects\simple_modflow\modflow\mp3du\bin\writep3dgsf.exe")
        self.writep3dgsf_path = writegsf_default if writep3dgsf_path is None else writep3dgsf_path
        self.mp3du_path = mp3du_default if mp3du_path is None else mp3du_path
        self._model_output_files = model_output_files
        self._output_path = output_path
        self._porosities_by_layer = porosities_by_layer
        self.path_file_path = self.output_path.joinpath('mp3du.p3d')
        self.gsf_output_path = self.output_path.joinpath(f'{self.model.name}.gsf')
        self.particle_shp = particle_shp
        self._variables = None


    @property
    def model_output_files(self):
        if self._model_output_files is None:
            model_name = self.model.name
            model_output_files = {
                'grb': f'{model_name}.grb',
                'tdis': f'{model_name}.tdis',
                'hds': f'{model_name}.hds',
                'cbc': f'{model_name}.cbc'
            }
            self._model_output_files = model_output_files
        return self._model_output_files

    @property
    def output_path(self) -> Path:
        if self._output_path is None:
            self._output_path = self.model.model_output_folder_path
        return self._output_path

    @property
    def porosities_by_layer(self):
        return self._porosities_by_layer

    @porosities_by_layer.setter
    def porosities_by_layer(self, val):
        assert isinstance(val, list), 'porosities_by_layer must be a list of porosities, one per layer'
        assert len(val) == self.model.gwf.modelgrid.nlay, 'number of porosities must match number of layers'
        for por in val:
            assert isinstance(por, (int, float)), f'porosity value: {por} must be an integer or a float (decimal)'
        self._porosities_by_layer = val

    @property
    def variables(self):
        if self._variables is None:
            porosities = self.porosities_by_layer
            vars = {
                'VELOCITY METHOD LAYER': 3,
                'POROSITY': porosities,
                'RETARDATION': 1,
                'DispH': 0,
                'DISPT': 0,
                'DISPV': 0
            }
            self._variables = vars
        return self._variables

    def create_gsf_file(self):
        # Use the provided grb file with writeP3DGSF.exe to create the GSF file
        gsf_output_path = self.gsf_output_path
        cmd = [self.writep3dgsf_path, self.model_output_files['grb'], gsf_output_path]
        run = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        if run.returncode != 0:
            print(f"Error running writeP3DGSF.exe: {run.stderr}")
        else:
            print(f"GSF file created successfully at {gsf_output_path}")

    def create_modflow_input_files(self):
        cwd = Path.cwd()
        # Copy files from the model path to the output path if necessary
        for typ, file in self.model_output_files.items():
            if not os.path.exists(cwd / file):
                original_file = self.model.model_output_folder_path / file
                if os.path.exists(original_file):
                    os.symlink(original_file, file)
                else:
                    raise FileNotFoundError(f"Required file {original_file} not found.")

    def create_path_file(self):
        # Create the PATH file with per-cell properties
        path_file_path = self.path_file_path
        porosities_by_layer = self.porosities_by_layer
        with open(path_file_path, 'w') as f:
            f.write("# PATH3D input file\n")
            for variable in self.variables:
                for layer in range(self.model.nlay):
                    if variable == 'POROSITY':
                        f.write(f"  CONSTANT    {variables[variable][layer]}   POROSITY {layer + 1}\n")
                    else:
                        f.write(f"  CONSTANT    {variables[variable]}   {variable} {layer + 1}\n")
        print(f"PATH file created at {path_file_path}")
        return path_file_path

    def run_mp3du(self, json_file_path):
        # Run the mp3du.exe with the created JSON file
        cmd = [self.mp3du_path, '-i', json_file_path]

        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        if result.returncode != 0:
            print(f"Error running mp3du.exe: {result.stderr}")
        else:
            print(f"MP3DU ran successfully. Output:\n{result.stdout}")

    def create_json_file(self):
        # Create the JSON configuration file
        json_data = {
            "FLOW_MODEL_TYPE": {
                "USGS_HFWK": {
                    "GRB_FILE": self.model_output_files['grb'],
                    "TDIS_FILE": self.model_output_files['tdis'],
                    "PATH_FILE": self.path_file_path,
                    "HDS_FILE": self.model_output_files['hds'],
                    "CBB_FILE": self.model_output_files['cbc'],
                    "GSF_FILE": {
                        "TYPE": "GSF_V.1.1.0",
                        "FILE_NAME": self.gsf_output_path
                    },
                    "OUTPUT_PRECISION": "DOUBLE",
                    "IFACE": [{"GHB": 7}, {"RCH": 6}, {"DRN": 7}, {"SFR": 6}, {"LAK": 7}],
                    "THREAD_COUNT": 1
                }
            },
            "SIMULATIONS": [
                {
                    "PATHLINE": {
                        "NAME": self.model.name,
                        "DIRECTION": "FORWARD",
                        "THREAD_COUNT": 1,
                        "INITIAL_STEPSIZE": 0.1,
                        "EULER_DT": 1.0e-7,
                        "ADAPTIVE_STEP_ERROR": 1.000000e-06,
                        "CAPTURE_RADIUS": 10.000000,
                        "OPTIONS": ["TRACK_TO_TERMINATION"],
                        "PARTICLE_START_LOCATIONS": {
                            "SHAPEFILE": {
                                "FILE_NAME": self.particle_shp,
                                "CELLID_ATTR": "Node",
                                "TIME_ATTR": "TimeRel",
                                "ZLOC_ATTR": "ZLoc",
                                "ADDTL_ATTR": ["LocName"]
                            }
                        }
                    }
                }
            ]
        }

        json_file_path = self.output_path / 'mp3du_input.json'
        with open(json_file_path, 'w') as f:
            json.dump(json_data, f, indent=4)

        print(f"JSON configuration file created at {json_file_path}")
        return json_file_path

    def run(self):
        self.create_gsf_file()
        self.create_modflow_input_files()
        path_file_path = self.create_path_file()
        json_file_path = self.create_json_file()
        self.run_mp3du(json_file_path)


if __name__ == "__main__":

    # Example usage
    model_name = 'your_model_name'
    model_path = 'path_to_your_model'
    output_path = 'path_to_output_directory'
    writep3dgsf_path = 'path_to_writeP3DGSF_executable'
    mp3du_path = 'path_to_mp3du_executable'
    grb_file_path = 'path_to_your_grb_file'

    pti = ParticleTrackingInput(
        model,
        porosities_by_layer=None,
        particle_shp=None
    )
    pti.run()
