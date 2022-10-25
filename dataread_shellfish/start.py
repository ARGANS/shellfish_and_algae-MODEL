import pandas as pd
import numpy as np
import datetime
import time
import os
import multiprocessing as mp
from src.make_runs import *
from src.launch_model import SF_model_scipy
from src.models.ModelProperties import ModelProperties
from src.read_netcdf import *

workdir = os.getenv('INPUT_DESTINATION')


model_properties = ModelProperties(os.getenv('INPUT_SOURCE'), os.getenv('INPUT_DESTINATION'))
# self, source_path, destination_path) 
try:
    model_properties.parse(os.getenv('INPUT_MODEL_PROPERTIES_JSON'))
except:
    raise RuntimeError('Cannot parse the value of the parameters_json environment variable')

if not model_properties.isDataDownloadTaskCompleted():
    raise RuntimeError('Data not downloaded')

full_json = model_properties.attrs

if full_json['metadata']['zone'] != "Europe":

    dict_dataCmd = full_json['dataset_parameters']['datasets']

    dict_to_AllData = dataCmd_to_AllData(dict_dataCmd, model_properties.file_template)
    shellData = AllData(dict_to_AllData)

    time_spent = run_simulation(f"{workdir}/concat.nc", full_json, shellData)

else:
    for area_name in ['IBI', 'NWS', 'MED', 'Baltic', 'BS', 'Arctic']:
        with open(f'/media/share/reference_data/{area_name}/parameters.json') as f:
            dict_dataCmd_area = json.load(f)['datasets']

        dict_to_AllData = dataCmd_to_AllData(dict_dataCmd_area,
                            '/media/share/reference_data_shellfish/{Place}/_pretreated/{Parameter}/{Parameter}{Place}modelNetCDF2021-01to2022-01.nc')
        shellData = AllData(dict_to_AllData)

        # Ensure that we call the correct year in case the user changed it.
        europe_json = full_json.copy()
        europe_json['dataset_parameters']['year'] = 2021

        time_spent = run_simulation(f"{workdir}/concat_{area_name}.nc", full_json, shellData)

        print(f'AREA {area_name} IS DONE. TIME SPENT: {time_spent/60} minutes.')
