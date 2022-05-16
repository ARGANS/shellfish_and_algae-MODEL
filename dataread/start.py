from make_runs import open_data_input, initialize_result, run_scenario_a_monthly
import pandas as pd
import numpy as np
import time
import os
import multiprocessing as mp
from launch_model import MA_model_scipy
from models.ModelProperties import ModelProperties

# TODO get dataCmd.csv from a mounted volume
dataRef: pd.DataFrame = pd.read_csv('./dataCmd.csv', delimiter=';')

model_properties = ModelProperties(os.getenv('DATASET_ID'), os.getenv('TASK_ID'))
try:
    model_properties.parse(os.getenv('PARAMETERS_JSON'))
except:
    raise RuntimeError('Cannot parse the value of the parameters_json environment variable')

if not model_properties.isDataDownloadTaskCompleted():
    raise RuntimeError('Data not downloaded')

input_args = {
    'zone' : model_properties.attrs['metadata']['zone'],
    'file_adress' : model_properties.file_template,
    'dataRef' : dataRef,
    'paramNames' : ['Ammonium', 'Nitrate', 'Temperature', 'northward_Water_current', 'eastward_Water_current']
}

### Initialize the netcdf reading interface
algaeData = open_data_input(**input_args)

### get the copernicus grid and mask

# TODO parametrize sim_area
sim_area = {
    'longitude': (-4, -3),
    'latitude': (48.5, 49),
    # 'longitude': (-180, 180),
    # 'latitude': (-90, 90),
    'time_index': 0,
    'depth': 3
}

longitudes, _ = algaeData.parameterData['Temperature'].getVariable('longitude', **sim_area)
latitudes, _ = algaeData.parameterData['Temperature'].getVariable('latitude', **sim_area)

mask1 = algaeData.parameterData['Temperature'].getVariable(**sim_area)[0].mask
mask2 = algaeData.parameterData['eastward_Water_current'].getVariable(**sim_area)[0].mask
mask3 = algaeData.parameterData['northward_Water_current'].getVariable(**sim_area)[0].mask
mask = np.logical_or(mask1, np.logical_or(mask2, mask3))
n_slices = 10
lon_split = np.array_split(longitudes, n_slices)
mask_split = np.array_split(mask, n_slices, 1)

model = MA_model_scipy(model_properties.parameters)

### Create datasets for monthly sim
for i, (lon_arr, mask_arr) in enumerate(zip(lon_split, mask_split)):
    initialize_result(
        model_properties.getMonthlySimulationsPath(i), 
        np.array(range(1,13)), latitudes, lon_arr, model.names, mask_arr)

pool = mp.Pool(10)
y0 = np.array([0, 0, 0, 1000, 0], dtype=np.float64)
t0 = time.time()

n_cells = pool.starmap_async(run_scenario_a_monthly, [(
        model_properties.getMonthlySimulationsPath(i),
        model_properties.parameters, 
        y0, 
        input_args, 
        model_properties.year, 
        True, 
        True
    ) for i in range(n_slices)]).get()

with open(model_properties.results_dir_path + '/stats.log', 'w') as file:
    file.write(str(n_cells) + '\n')
    file.write(str((time.time() - t0) / sum(n_cells)) + '\n')
