from pprint import pprint
from make_runs import open_data_input, initialize_result, run_scenario_a_monthly
import pandas as pd
import numpy as np
import time
import json
import os
import multiprocessing as mp
from launch_model import MA_model_scipy

# TODO provide dataCmd.csv from mounted volume
dataRef: pd.DataFrame = pd.read_csv('./dataCmd.csv', delimiter=';')

parameters_json_value:str = os.getenv('parameters_json')
try:
    input_parameters:dict = json.loads(parameters_json_value)
except:
    raise RuntimeError('Cannot parse the value of the parameters_json environment variable')

# TODO safe in a manifest file all information about zone and paramNames
input_args = {
    'zone' : 'IBI',
    'file_adress' : '/media/share/data/{zone}/{param}/{param}{zone}modelNetCDF2021-01to2022-01.nc',
    'dataRef' : dataRef,
    'paramNames' : ['Ammonium', 'Nitrate', 'Temperature', 'northward_Water_current', 'eastward_Water_current']
}

### Initialize the netcdf reading interface
algaeData = open_data_input(**input_args)
print('algaeData')
pprint(algaeData)

### get the copernicus grid and mask

# TODO parametrize sim_area
sim_area = {
    'longitude': (-4, -3),
    'latitude': (48.5, 49),
    #'longitude': (-180, 180),
    #'latitude': (-90, 90),
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

model = MA_model_scipy(input_parameters['parameters'])

### Create datasets for monthly sim
for i, (lon_arr, mask_arr) in enumerate(zip(lon_split, mask_split)):
    initialize_result(f"/media/share/results/simulations/monthly/monthly_simulations_{i:03d}.nc",
    np.array(range(1,13)), latitudes, lon_arr, model.names, mask_arr)

pool = mp.Pool(10)
y0 = np.array([0, 0, 0, 1000, 0], dtype=np.float64)
t0 = time.time()

n_cells = pool.starmap_async(run_scenario_a_monthly,
    [(f"/media/share/results/simulations/monthly/monthly_simulations_{i:03d}.nc", 
        input_parameters['parameters'], y0, input_args, 2021, True, True) for i in range(n_slices)]).get()
# n_cells = run_scenario_a_monthly("/media/share/results/complete_simulations_monthly_test_0.nc", model, y0, input_args, 2021)

print('RESULTS')
pprint(n_cells)
pprint((time.time()-t0)/sum(n_cells))