import os
import json
import subprocess
import flopy
from simple_modflow import VoronoiGridPlus
from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
import pandas as pd
from simple_modflow.modflow.utils.datatypes.readers import read_shp_gpkg
from pathlib import Path
import pickle
import shutil
from flopy.utils import CellBudgetFile
from collections import OrderedDict
import numpy as np


class ParticleTrackingInput:
    def __init__(
            self,
            model: SimulationBase = None,
            writep3dgsf_path: Path = None,
            mp3du_path: Path = None,
            model_output_files: dict = None,
            output_path: Path = None,
            porosities_by_layer=None,
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
        self.particle_shp = particle_shp
        self._variables = None

    @property
    def model_output_files(self):
        if self._model_output_files is None:
            model_name = self.model.name
            model_output_files = {
                'grb': f'{model_name}.disv.grb',
                'tdis': f'{model_name}.tdis',
                'hds': f'{model_name}.hds',
                'cbc': f'{model_name}.cbc',
                'gsf': f'{model_name}.gsf'
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
        gsf_json = {
            "FLOW_MODEL_TYPE": {
                "USGS_HFWK": {
                    "GRB_FILE": self.model_output_files['grb'],
                    "GSF_FILE": {
                        "TYPE": "HFWK_GRB_V.1.0.0"
                    }
                }
            },
            "OUTPUT_FILENAME": f"{self.model.name}.gsf"
        }
        gsf_json_file_path = self.output_path / 'grb_to_gsf.json'
        with open(gsf_json_file_path, 'w') as f:
            json.dump(gsf_json, f, indent=4)

        cmd = [self.writep3dgsf_path.as_posix(), gsf_json_file_path.as_posix(), 'colorcode']
        run = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        if run.returncode != 0:
            print(f"Error running writeP3DGSF.exe: {run.stderr}")
        else:
            print("GSF file created successfully")

    def create_modflow_input_files(self):
        """
        Creates necessary MODFLOW input files by copying them from the model's output folder to the
        current working directory if they do not already exist.

        If the required files are missing in the current working directory, they are copied from the
        model's output folder path. Raises an error if the files do not exist in the source folder.

        :raises FileNotFoundError: If a required file does not exist in the model's output folder path.
        """
        cwd = Path.cwd()
        # Copy files from the model path to the output path if necessary
        for typ, file in self.model_output_files.items():
            if not os.path.exists(cwd / file):
                original_file = self.model.model_output_folder_path / file
                if os.path.exists(original_file):
                    # os.symlink(original_file, file)
                    shutil.copyfile(original_file, file)  # Copies the file instead of linking it
                elif typ == 'gsf':
                    print(f"GSF file not found. Creating a new one from the grb file.")
                    continue
                else:
                    raise FileNotFoundError(f"Required file {original_file} not found.")

    def create_path_file(self):
        # Create the PATH file with per-cell properties
        path_file_path = self.path_file_path
        porosities_by_layer = self.porosities_by_layer
        with open(path_file_path, 'w') as f:
            f.write("# PATH3D input file\n\n")
            for variable in self.variables:
                for layer in range(self.model.gwf.modelgrid.nlay):
                    if variable == 'POROSITY':
                        f.write(f"  CONSTANT    {self.variables[variable][layer]}   POROSITY {layer + 1}\n")
                    else:
                        f.write(f"  CONSTANT    {self.variables[variable]}   {variable} {layer + 1}\n")
        print(f"PATH file created at {path_file_path}")
        return path_file_path

    def run_mp3du(self, json_file_path):
        # Run the mp3du.exe with the created JSON file
        cmd = [self.mp3du_path.as_posix(), json_file_path, 'colorcode']
        print(cmd)
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
                    "PATH_FILE": self.path_file_path.as_posix(),
                    "HDS_FILE": self.model_output_files['hds'],
                    "CBB_FILE": self.model_output_files['cbc'],
                    "GSF_FILE": {
                        "TYPE": "GSF_V.1.1.0",
                        "FILE_NAME": self.model_output_files['gsf']
                    },
                    "OUTPUT_PRECISION": "DOUBLE",
                    "IFACE": [{"RCH": 6}, {"DRN": 7}, {"SFR": 6}],
                    "THREAD_COUNT": 4
                }
            },
            "SIMULATIONS": [
                {
                    "PATHLINE": {
                        "NAME": self.model.name,
                        "DIRECTION": "BACKWARD",
                        "THREAD_COUNT": 4,
                        "INITIAL_STEPSIZE": 0.1,
                        "EULER_DT": 1.0e-7,
                        "ADAPTIVE_STEP_ERROR": 1.000000e-04,
                        "CAPTURE_RADIUS": 10.000000,
                        "SIMULATION_END_TIME": 0,
                        "OPTIONS": ["TRACK_TO_TERMINATION"],
                        "PARTICLE_START_LOCATIONS": {
                            "SHAPEFILE": {
                                "FILE_NAME": self.particle_shp.as_posix(),
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
        self.create_modflow_input_files()
        self.create_gsf_file()
        path_file_path = self.create_path_file()
        json_file_path = self.create_json_file().as_posix()
        self.run_mp3du(json_file_path)

    def get_output_json(self):

        output_json = {
            "MP3DU_BIN": "cumb_v2b_PATHLINE.bin",
            "OUTPUTS": [
                {"SUMMARY": {
                }},
                {"DBF_TABLE": {
                    "FILE_NAME": "NAME_OF_OUTPUT.dbf"
                }},
                {"PATHLINE_WHOLE": {
                    "FILE_NAME": "NAME_01_OF_OUTPUT.shp"
                }},
                {"PATHLINE_PARTS": {
                    "FILE_NAME": "NAME_02_OF_OUTPUT.shp"
                }},
                {"POINTS_IN_TIME": {
                    "FILE_NAME": "NAME_03_OF_OUTPUT.shp"
                }},
                {"ENDPOINT": {
                    "FILE_NAME": "NAME_04_OF_OUTPUT.shp"
                }}
            ]
        }
        with open('P3DOutput_json.json', 'w') as f:
            json.dump(output_json, f, indent=4)


class PRT:

    def __init__(
            self,
            model: SimulationBase = None,
            name: str = None,
            mf_folder_path: Path = Path().home().joinpath('mf6'),
    ):
        """
        Initializes an instance of a class responsible for managing a MODFLOW 6 model
        simulation setup using FloPy, setting up the simulation environment and
        configuring the PRT (particle tracking) model. The initializer associates a
        simulation model, configures a unique name for the model, and prepares the folder
        paths for storing simulation outputs. It also creates a FloPy MFSimulation
        object and a ModflowPRT model instance based on the provided configuration.

        :param model: The associated simulation model, implementing a basic simulation
            representation to integrate with FloPy.
        :type model: SimulationBase
        :param name: The unique name for the MODFLOW simulation. If not provided, it
            defaults to the name of the given model.
        :type name: str, optional
        :param mf_folder_path: The base folder path where the MODFLOW files and all
            related outputs will be stored. The user can specify an alternative path,
            or it defaults to the user's home directory appended with 'mf6'.
        :type mf_folder_path: Path, optional
        """
        self.model = model
        self.name = self.model.name if name is None else name
        self.model_output_folder_path = mf_folder_path.joinpath(f'{name}')
        self.nper = self.model.gwf.modeltime.nper
        self.bud_file = self.model.model_output_folder_path.joinpath(f'{model.name}.cbc')
        self.hds_file = self.model.model_output_folder_path.joinpath(f'{model.name}.hds')

        self.prt_sim = flopy.mf6.MFSimulation(
            sim_name='sim_prt',
            exe_name="mf6",
            version="mf6",
            sim_ws=self.model_output_folder_path
        )
        self.tdis = flopy.mf6.modflow.ModflowTdis(
            simulation=self.prt_sim,
            time_units=self.model.gwf.modeltime.time_units,
            # pname="tdis-prt",
            nper=1,
            perioddata=[(1, 1, 1)],
            # perioddata=self.generate_tdis_from_cbc(self.bud_file),
        )
        self.prt = flopy.mf6.modflow.ModflowPrt(
            simulation=self.prt_sim,
            modelname=f'{self.name}.prt',
            model_nam_file=f'{self.name}.prt.nam',
            version="mf6",
            exe_name="mf6",
            print_input=True,
            print_flows=False,
            save_flows=True,
        )
        self.ims = flopy.mf6.ModflowIms(
            simulation=self.prt_sim,
            pname="ims",
            complexity="COMPLEX",
        )

    @staticmethod
    def generate_tdis_from_cbc(cbc_path: str) -> list[tuple[float, int, float]]:
        """
        Generate TDIS perioddata using only the last saved time step of each stress period
        in a MODFLOW 6 .cbc budget file. Intended for use with transport models that
        read only the saved (last) time step per period.

        Parameters
        ----------
        cbc_path : str
            Path to the MODFLOW 6 .cbc budget file

        Returns
        -------
        list of (perlen, nstp, tsmult) tuples suitable for TDIS
        """
        cbc = CellBudgetFile(cbc_path, precision="double")
        kstpkper = cbc.get_kstpkper()
        times = cbc.get_times()

        last_time_by_kper = OrderedDict()
        prev_time = 0.0

        for (kstp, kper), totim in zip(kstpkper, times):
            if kper not in last_time_by_kper or kstp > last_time_by_kper[kper][0]:
                dt = totim - prev_time
                last_time_by_kper[kper] = (kstp, dt)
            prev_time = totim

        # Format: (perlen, nstp, tsmult)
        perioddata = [(dt, 1, 1.0) for (kstp, dt) in last_time_by_kper.values()]
        print(perioddata)
        return perioddata





class PrtMip:

    def __init__(
            self,
            prt_model: PRT = None,
            porosity: int | list = 0.2,
            retfactor: int = 1,
            izone: int = 0,
    ):
        self.prt_mip = flopy.mf6.modflow.ModflowPrtmip(
            model=prt_model.prt,
            porosity=porosity,
            retfactor=retfactor,
            izone=izone,
            filename=f'prt-{prt_model.name}.mip',
            pname='mip'
        )


class PrtOc:

    def __init__(
            self,
            prt_model: PRT = None,
            save_record=("BUDGET", "LAST"),
            print_record=None
    ):
        self.prt_oc = flopy.mf6.ModflowPrtoc(
            model=prt_model.prt,
            pname='oc',
            filename=f'prt-{prt_model.name}.oc',
            budget_filerecord=f'prt-{prt_model.name}.cbc',
            track_filerecord=f'prt-{prt_model.name}.trk',
            saverecord=save_record,
            printrecord=print_record,
        )


class PrtPrp:
    def __init__(
            self,
            prt_model: PRT = None,
            exit_solve_tolerance=0.00001,
            # dev_exit_solve_method=1,  # Brent method = 1
            stoptime=None,
            stoptraveltime=None,
            istopzone=0,
            shp_gpkg_path: Path = None,
            vor: VoronoiGridPlus = None,
            local_z=0.5
    ):

        particle_data = read_shp_gpkg(shp_gpkg_path).geometry
        pnts_vor = vor.get_vor_cells_as_series(particle_data).to_list()
        pnt_cells = vor.gdf_vorPolys.loc[pnts_vor].geometry
        nprt = len(pnt_cells)
        packagedata = []
        for i, pnt in enumerate(pnt_cells):
            packagedata.append(
                [
                    i, (0, pnts_vor[i]), pnt.centroid.x, pnt.centroid.y, local_z
                ]
            )

        self.prt_prp = flopy.mf6.modflow.ModflowPrtprp(
            model=prt_model.prt,
            pname='prtprp',
            filename=f'prt-{prt_model.name}.prp',
            print_input=True,
            # dev_exit_solve_method=dev_exit_solve_method,
            # exit_solve_tolerance=exit_solve_tolerance,
            local_z=True,
            # track_filerecord=f'prt-{prt_model.name}.trk',
            # stoptime=stoptime,
            # stoptraveltime=stoptraveltime,
            stop_at_weak_sink=True,
            # istopzone=istopzone,
            drape=True,
            nreleasepts=nprt,
            packagedata=packagedata,
            # perioddata=["first"] * prt_model.nper
        )


class PrtDisv:

    def __init__(
            self,
            gwf_model: SimulationBase = None,
            prt_model: PRT = None,
    ):
        modelgrid = gwf_model.modelgrid
        vor = gwf_model.vor
        self.prt_disv = flopy.mf6.modflow.ModflowPrtdisv(
            model=prt_model.prt,
            pname='disv',
            filename=f'prt-{prt_model.name}.disv',
            length_units='feet',
            export_array_ascii=False,
            nlay=modelgrid.nlay,
            ncpl=modelgrid.ncpl,
            nvert=modelgrid.nvert,
            top=modelgrid.top,
            botm=modelgrid.botm,
            idomain=modelgrid.idomain,
            vertices=vor.get_disv_gridprops()['vertices'],
            cell2d=modelgrid.cell2d
        )


class PrtFmi:

    def __init__(
            self,
            prt_model: PRT = None,
            gwf_model: SimulationBase = None,
    ):

        self.prt_fmi = flopy.mf6.modflow.ModflowPrtfmi(
            save_flows=True,
            model=prt_model.prt,
            filename=f'prt-{prt_model.name}.fmi',
            pname='fmi',
            packagedata=[
                ['GWFBUDGET', prt_model.bud_file.as_posix()],
                ['GWFHEAD', prt_model.hds_file.as_posix()],
            ]
        )


if __name__ == "__main__":
    with open(Path(r"C:\Users\lukem\mf6\cum8cNoET\cum8cNoET.model"), 'rb') as file:
        model: SimulationBase = pickle.load(file)
    particles = Path(r"C:\Users\lukem\mf6\Cumberland general\particles\edge_particles.shp")
    pti = ParticleTrackingInput(
        model=model,
        porosities_by_layer=[0.25],
        particle_shp=particles
    )
    pti.run()
    pti.get_output_json()
