import pandas as pd
import numpy as np
import datetime
import time
import os
import multiprocessing as mp
from src.make_runs import *
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

dict_dataCmd = full_json['dataset_parameters']['datasets']
dict_to_AllData = dataCmd_to_AllData(dict_dataCmd, model_properties.file_template)

algaeData = AllData(dict_to_AllData)


time_spent = run_simulation(f"{workdir}/concat.nc", full_json, algaeData)#,'/media/global/ZEE_saccharina_scA_farms-for-5MT-FW_surf1.00000_depth-30_dist10_posfarms.txt')