import json
import os
from src.make_runs import *
# from src.make_runs import run_simulation
from src.launch_model import MA_model_scipy
from src.models.ModelProperties import ModelProperties
from src.read_netcdf import *

workdir = os.getenv('INPUT_DESTINATION')

model_properties = ModelProperties(os.getenv('INPUT_SOURCE'), os.getenv('INPUT_DESTINATION'))
try:
    model_properties.parse(os.getenv('INPUT_MODEL_PROPERTIES_JSON'))
except:
    raise RuntimeError('Cannot parse the value of the parameters_json environment variable')

if not model_properties.isDataDownloadTaskCompleted():
    raise RuntimeError('Data not downloaded')

full_json = model_properties.attrs


if os.getenv('RUN_SIMULATION_WITH_FARMS') is not None:
    farm_pos_file = '/_farmdistribution/opt_posfarms.txt'
else:
    farm_pos_file = None


if full_json['metadata']['zone'] != "Europe":

    dict_dataCmd = full_json['dataset_parameters']['datasets']

    dict_to_AllData = dataCmd_to_AllData(dict_dataCmd, model_properties.file_template)
    algaeData = AllData(dict_to_AllData)
    time_spent = run_simulation(
        f"{workdir}/concat.nc",
        full_json,
        algaeData,
        farm_pos_file=farm_pos_file
    )

else:
    for area_name in ['IBI', 'NWS', 'MED', 'Baltic', 'BS', 'Arctic']:
        with open(f'/media/share/reference_data/{area_name}/parameters.json') as f:
            dict_dataCmd_area = json.load(f)['datasets']

        europe_json = full_json.copy()

        if area_name == 'Baltic': # No good data for 2021 Baltic so 2020 is used
            europe_json['dataset_parameters']['year'] = 2020
            dict_to_AllData = dataCmd_to_AllData(dict_dataCmd_area,
                                '/media/share/reference_data/{Place}/_pretreated/{Parameter}/{Parameter}{Place}modelNetCDF2020-01to2021-01.nc')
        else:
            europe_json['dataset_parameters']['year'] = 2021
            dict_to_AllData = dataCmd_to_AllData(dict_dataCmd_area,
                                '/media/share/reference_data/{Place}/_pretreated/{Parameter}/{Parameter}{Place}modelNetCDF2021-01to2022-01.nc')

        algaeData = AllData(dict_to_AllData)
        time_spent = run_simulation(
            f"{workdir}/concat_{area_name}.nc", 
            europe_json, 
            algaeData, 
            farm_pos_file=farm_pos_file
        )

        print(f'AREA {area_name} IS DONE. TIME SPENT: {time_spent/60} minutes.')
