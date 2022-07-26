import datetime
import os

import numpy as np
import pandas as pd
from numpy import ma
import netCDF4 as nc

from advectionPrototype.climatology_ellipses import degrees_to_meters
from advectionPrototype.quickestHadley import quickest, sortData, giveResol, givecoor, v2d_cgrid_cur, u2d_cgrid_cur, \
    giveCFL
from dataread.launch_model import MA_model_scipy
from dataread.make_runs import open_data_input
from dataread.models.ModelProperties import ModelProperties
from dataread.read_netcdf import AllData

dataRef: pd.DataFrame = pd.read_csv('/media/global/dataCmd.csv', delimiter=';')

model_properties = ModelProperties(os.getenv('DATASET_ID'), os.getenv('TASK_ID'))
try:
    model_properties.parse(os.getenv('PARAMETERS_JSON'))
except:
    raise RuntimeError('Cannot parse the value of the parameters_json environment variable')

if not model_properties.isDataDownloadTaskCompleted():
    raise RuntimeError('Data not downloaded')

parms_run = list(model_properties.parameters['run'].values())[0]['parameters']
parms_farm = list(model_properties.parameters['farm'].values())[0]['parameters']

paramNames = ['Nitrate', 'northward_Water_current', 'Ammonium', 'eastward_Water_current',
                      'Temperature']

dataCmdpath = './../global/dataCmd.csv'

input_args = {
        'zone': model_properties.attrs['metadata']['zone'],
        'file_adress': model_properties.file_template,
        'dataRef': dataRef,
        'paramNames': paramNames,
        'frequency': 'daily'
    }

dict_to_AllData = open_data_input(**input_args)

dict_to_AllData['PAR'] = {
    'file_name': model_properties.file_template,
    'variable_name': 'par',
    'latitude_name': 'lat',
    'longitude_name': 'lon',
    'time_name': 'time',
    'depth_name': 'depth',
    'unit_conversion': 11.574,
    'time_zero': datetime.datetime(model_properties.year, 1, 1),
    'time_step': datetime.timedelta(days=1)
}

sim_area = {
    'longitude': (parms_run['min_lon'], parms_run['max_lon']),
    'latitude': (parms_run['min_lat'], parms_run['max_lat']),
    'depth': (0, parms_farm['z']*1.4),
    'averagingDims': ('depth',),
    'weighted': False
}

### Initialize the netcdf reading interface
algaeData = AllData(dict_to_AllData)

dataNO3 = algaeData.parameterData['Nitrate'].getVariable(**sim_area)[0]
dataNH4 = algaeData.parameterData['Ammonium'].getVariable(**sim_area)[0]
dataTemp = algaeData.parameterData['Temperature'].getVariable(**sim_area)[0]
dataNwc = algaeData.parameterData['northward_Water_current'].getVariable(**sim_area)[0]
dataEwc = algaeData.parameterData['eastward_Water_current'].getVariable(**sim_area)[0]
dataPAR = algaeData.parameterData['PAR'].getVariable(**sim_area)[0]

sortedList = sortData(parms_run['dateBeginning'], parms_run['dateEnd'], len(dataNO3))

dataNH4 = ma.masked_outside(dataNH4, -1e4, 1e4)
dataNO3 = ma.masked_outside(dataNO3, -1e4, 1e4)
dataTemp = ma.masked_outside(dataTemp, -1e4, 1e4)
dataPAR = ma.masked_outside(dataPAR, -1e-2, 1e4)
dataPAR = dataPAR.filled(fill_value=8)

dataFin = pd.read_csv('./../global/dataCmd.csv', ';')

nwcDataLine = dataFin.loc[(dataFin["Parameter"] == 'Nitrate') & (dataFin["Place"] == model_properties.attrs['metadata']['zone'])]
nwcdataName = nwcDataLine.iloc[-1]["variable"]  # we find the data name in the dataset
resx, resy, km = giveResol(nwcDataLine)

longitudes, _ = algaeData.parameterData['Temperature'].getVariable('longitude', **sim_area)
latitudes, _ = algaeData.parameterData['Temperature'].getVariable('latitude', **sim_area)

longitudeMin, latitudeMin = givecoor(longitudes,latitudes, parms_run['min_lon'], parms_run['min_lat'])  # we get the indices of the wanted position
longitudeMax, latitudeMax = givecoor(longitudes,latitudes, parms_run['max_lon'], parms_run['max_lat'])  # we get the indices of the wanted position

firstday = datetime.datetime.strptime(parms_run['dateBeginning'], '%Y-%m-%d %H:%M:%S')

latRef = np.ones((np.shape(dataEwc[0])[1], np.shape(dataEwc[0])[0])) * latitudes[latitudeMin:latitudeMax]
decenturedEwc = u2d_cgrid_cur(dataEwc)
decenturedNwc = v2d_cgrid_cur(dataNwc)

if km:
    dxlist, dyMeter = resx * np.ones(np.shape(dataNwc[0])), resy  # baltic
else:
    dxlist, dyMeter = degrees_to_meters(resx, resy, latRef)
    dxlist = dxlist.T

dylist = dyMeter * np.ones(np.shape(dxlist))


maxCFL = 0
(nbrx, nbry) = np.shape(dataNO3[0])
for i in range(len(dataEwc)):
    CFL, cx, cy = giveCFL(dxlist, dyMeter, parms_run['dt'], decenturedEwc[i], decenturedNwc[i], nbrx, nbry)
    if np.max(CFL) > maxCFL:
        maxCFL = np.max(CFL)
if maxCFL>1:
    raise RuntimeError('CFL > 1')
model = MA_model_scipy(model_properties.parameters)
NO3field, NH4field, D, N_f, N_s, totalNH4deficit, totalNO3deficit = quickest(dyMeter, dxlist, parms_run['dt'],
                                                                             decenturedEwc, decenturedNwc, dataEwc,
                                                                             dataNwc, latRef.T, dataNO3, dataNH4,
                                                                             dataTemp, dataPAR, parms_run['Ks'], firstday, model,
                                                                             parms_farm['z']*1.4, parms_run['scenC'],sortedList)