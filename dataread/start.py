import pandas as pd
import numpy as np
import datetime
import time
import os
import multiprocessing as mp
from src.make_runs import open_data_input, initialize_result, run_scenario_a_monthly
from src.launch_model import MA_model_scipy
from src.models.ModelProperties import ModelProperties
from src.read_netcdf import *

dataRef: pd.DataFrame = pd.read_csv('/media/global/dataCmd.csv', delimiter=';')

model_properties = ModelProperties(os.getenv('INPUT_SOURCE'), os.getenv('INPUT_DESTINATION'))
try:
    model_properties.parse(os.getenv('INPUT_MODEL_PROPERTIES_JSON'))
except:
    raise RuntimeError('Cannot parse the value of the parameters_json environment variable')

if not model_properties.isDataDownloadTaskCompleted():
    raise RuntimeError('Data not downloaded')

input_args = {
    'zone' : model_properties.attrs['metadata']['zone'],
    'file_adress' : model_properties.file_template,
    'dataRef' : dataRef,
    'paramNames' : ['Ammonium', 'Nitrate', 'Phosphate', 'Temperature', 'northward_Water_current', 'eastward_Water_current']
}

dict_to_AllData = open_data_input(**input_args)

#TODO: make a different PAR file for each year ?
dict_to_AllData['PAR'] = {
            # 'file_name': '/media/share/PAR/PAR_ref_OC_2020.nc',
            'file_name': '/media/global/PAR_ref_OC_2020.nc',
            'variable_name': 'par',
            'latitude_name': 'lat',
            'longitude_name': 'lon',
            'time_name': 'time',
            'depth_name': 'depth',
            'unit_conversion': 11.574,
            'time_zero': datetime.datetime(model_properties.year, 1, 1),
            'time_step': datetime.timedelta(days=1)
        }

### Initialize the netcdf reading interface
algaeData = AllData(dict_to_AllData)

### get the copernicus grid and mask

parms_run = list(model_properties.parameters['run'].values())[0]['parameters']
parms_farm = list(model_properties.parameters['farm'].values())[0]['parameters']
parms_harvest = list(model_properties.parameters['harvest'].values())[0]['parameters']
harvest_type = list(model_properties.parameters['harvest'].keys())[0]

sim_area = {
    # 'longitude': (-4, -3),
    # 'latitude': (48.5, 49),
    'longitude': (parms_run['min_lon'], parms_run['max_lon']),
    'latitude': (parms_run['min_lat'], parms_run['max_lat']),
    'time_index': 0,
    'depth': parms_farm['z']
}

longitudes, _ = algaeData.parameterData['Temperature'].getVariable('longitude', **sim_area)
latitudes, _ = algaeData.parameterData['Temperature'].getVariable('latitude', **sim_area)

mask1 = algaeData.parameterData['Temperature'].getVariable(**sim_area)[0].mask
#mask2 = algaeData.parameterData['eastward_Water_current'].getVariable(**sim_area)[0].mask
#mask3 = algaeData.parameterData['northward_Water_current'].getVariable(**sim_area)[0].mask
#mask = np.logical_or(mask1, np.logical_or(mask2, mask3))

n_slices = parms_run['n_cores']
lon_split = np.array_split(longitudes, n_slices)
mask_split = np.array_split(mask1, n_slices, 1)

model = MA_model_scipy(model_properties.parameters)

if harvest_type == "Summer_growth":
    time_axis = np.array(range(parms_harvest['deployment_month'], parms_harvest['harvesting_month'] + 1))
elif harvest_type == "Winter_growth":
    time_axis = np.array(list(range(parms_harvest['deployment_month'], 13)) +
                         list(range(1, parms_harvest['harvesting_month'] + 1)))

### Create datasets for monthly sim
for i, (lon_arr, mask_arr) in enumerate(zip(lon_split, mask_split)):
    initialize_result(
        model_properties.getMonthlySimulationsPath(i),
        time_axis, latitudes, lon_arr, model.names, mask_arr)

pool = mp.Pool(parms_run['n_cores'])
y0 = np.array([0, 0, 0, parms_harvest['deployment_Nf'], 0], dtype=np.float64)
t0 = time.time()

n_cells = pool.starmap_async(run_scenario_a_monthly, [(
        model_properties.getMonthlySimulationsPath(i),
        model_properties.parameters, 
        y0, 
        dict_to_AllData, 
        model_properties.year, 
        False, 
        True
    ) for i in range(n_slices)]).get()

with open(model_properties.destination_path + '/stats.log', 'w') as file:
    file.write(str(n_cells) + '\n')
    file.write(str((time.time() - t0) / sum(n_cells)) + '\n')


def main():
    print(1)
