### Shellfish model for EMFF shellfish and algae project

## Full details of scientific rationale can be found in associated ATBD
# Authors
# F. Coulibaly (ARGANS, France)  fcoulibaly@argans.eu to pass from R to Python language based M. Johnson (Bantry Marine Research Station, Ireland) mjohnson@bmrs.ie,


from dis import dis
import os
import sys
import copy
import datetime
from unittest import result
from dateutil.relativedelta import *
from fileinput import filename
from tkinter import N
from utils import import_json
import numpy as np
import numpy.ma as ma
import pandas as pd
from read_netcdf import *
from launch_model_Shellfish import *
import time
from scipy.integrate import solve_ivp
from scipy.integrate import OdeSolver
from scipy.interpolate import interp1d
from scipy.sparse import spdiags
import multiprocessing as mp

def initialize_result(fileName:str, times, latitudes, longitudes, 
                      variableNames:list, mask:np.array):
    """Initializes a netcdf file at fileName with a time, latitude, and
    longitude dimension. Corresponding variables are created with the values
    in times, latitudes, longitudes. These also defin the size of the
    dimensions.
    Variables are also created from the names stored in variableNames with
    dimensions (time, latitude, longitude).
    """

    ds = nc.Dataset(fileName, 'w', format='NETCDF4')

    days = ds.createDimension('time', len(times))
    lat = ds.createDimension('latitude', len(latitudes))
    lon = ds.createDimension('longitude', len(longitudes))

    tims = ds.createVariable('time', 'f4', ('time',))
    lats = ds.createVariable('latitude', 'f4', ('latitude',))
    lons = ds.createVariable('longitude', 'f4', ('longitude',))

    tims[:] = times
    lats[:] = latitudes
    lons[:] = longitudes

    full_mask = np.repeat(np.expand_dims(mask, 0), len(times), axis=0)
    print("full mask shape", np.shape(full_mask))

    for name in variableNames:
        var = ds.createVariable(name, 'f4', ('time', 'latitude', 'longitude',))
        var[:,:,:] = np.ma.masked_array(-1*np.ones(full_mask.shape), full_mask)
    ds.close()

def open_data_input(file_adress:str, zone:str, paramNames:list, dataRef: pd.DataFrame):

    fileNames = [file_adress.format(zone=zone, param=param) for param in paramNames]

    dataRows = [dataRef.index[(dataRef['Parameter']==param) & (dataRef['Place']==zone) &
                              (dataRef['type']=='model') & (dataRef['daily']=='monthly')][0]
                    for param in paramNames]

    variableNames = dataRef.iloc[dataRows]['variable'].tolist()
    latitudeNames = dataRef.iloc[dataRows]['latName'].fillna('latitude').tolist()
    longitudeNames = dataRef.iloc[dataRows]['longName'].fillna('longitude').tolist()
    timeNames = dataRef.iloc[dataRows]['timeName'].fillna('time').tolist()
    depthNames = dataRef.iloc[dataRows]['depthName'].fillna('depth').tolist()
    unitConversions = dataRef.iloc[dataRows]['unitFactor'].fillna(1).tolist()
    timeUnits = dataRef.iloc[dataRows]['time_units'].tolist()

    data = AllData(fileNameList=fileNames,
                   parameterNameList=paramNames,
                   variableNameList=variableNames,
                   latitudeNameList=latitudeNames,
                   longitudeNameList=longitudeNames,
                   timeNameList=timeNames,
                   depthNameList=depthNames,
                   unitConversionList=unitConversions,
                   timeUnitsList=timeUnits
    )
    return data

def giveEpsilon(temp, CHL, POP, spawnday, sd2, cCHL, SHE, STE,
                                    Nwc, Ewc, model):    
    data_in = {
        'SST': temp,
        'CHL_ext': CHL,
        'F_in': np.sqrt(Nwc ** 2 + Ewc ** 2),
    }

    y = dict(zip(["CHL", "SHE", "STE", "spawnday", "sd2", "POP"], [cCHL, SHE, STE, spawnday,sd2, POP]))
    return model.SF_derivative(y, data_in, model)

def run_scenario_a_monthly(fileName, year, firstday, sortedList, dt, dataEwc, dataNwc, dataTemp, dataCHL, model):
    """Runs simulation on all the grid points in fileName, reading data from
    inputData, and initializing at y0. The simulations are ran on monthly
    averaged data.
    Writes the results in fileName.
    input_args are passed to open_data_input() to open an AllData object.

    farm_at_gridSize: If True, the farm size parameters in the model are
        overruled and instead, a square farm with an area equal to the grid size
        is used. The grid size is obtained from latitude and longitude axes,
        they must have a constant step in degrees for this option to work
        properly.
    data_is_monthly: If True, the data opened by passing input_args to
        open_data_input() is assumed to already be averaged monthly. The data
        is then accessed by searching for the value that is nearest to the
        15th of each month.
    """
    ds = nc.Dataset(fileName, 'r')
    mask = ds['CHL'][0,:,:].mask
    times = ds['time'][:]
    ds.close()

    nbrStep = int(len(sortedList)*dt)
    print('nbrStep', nbrStep)
    daylist = [firstday]

    SHE = np.zeros(np.shape(dataCHL[0])) + 100
    STE = np.ones(np.shape(dataCHL[0]))*1000
    spawnday = np.zeros(np.shape(dataCHL[0])) + 365
    sd2 = np.zeros(np.shape(dataCHL[0])) + 365
    cCHL = np.zeros(np.shape(dataCHL[0]))
    POP = np.zeros(np.shape(dataCHL[0])) + 3000
    mask = dataCHL.mask[0,:,:]

    for k in range (nbrStep) :
        daylist += [firstday + relativedelta(minutes=int(k * 24 * 60 / dt))]

        dayNbr = sortedList[k // int(dt)]
        hNbr = dayNbr #if we don't use hourly speeds, hNbr = dayNbr
        month = dayNbr // int(30.5*dt)
        month = int(month)
        #startTime = datetime.datetime(year, month, 1, 0)
        if month == 12:
            endTime = datetime.datetime(year+1, 1, 1, 0)
        else:
            endTime = datetime.datetime(year, month+1, 1, 0)

        #values = np.ma.masked_array(-1*np.ones(y0_array.shape), full_mask)
        #days_start = (startTime - initTime).days
        #days_end = (endTime - initTime).days

        cCHL[mask] = 1e-5
        POP[mask] = 1e-5
        SHE[mask] = 1e-5
        STE[mask] = 1e-5
        spawnday[mask] = 1e-5
        sd2[mask] = 1e-5
        days = (daylist[-1] - datetime.datetime(daylist[-1].year, 1, 1,
                                                0)).days  # we get the number of the day in the year

        derivArray = giveEpsilon(dataTemp[dayNbr], dataCHL[dayNbr], POP, spawnday,sd2, cCHL, SHE, STE, dataNwc[dayNbr], dataEwc[dayNbr], model)
        #like np.array([dSHE, dSTE, dspawnday, dsd2, dPOP, dCHL_farm, dPOC_farm, dPHYC_farm, FW, DWW, SHL, NH4_production, CO2_production])
        cCHL += derivArray[5]*dt
        POP += derivArray[4]*dt
        SHE += derivArray[0] * dt
        STE += derivArray[1]*dt
        spawnday += derivArray[2]*dt
        sd2 += derivArray[3]*dt
        #maybe do here the spawing
        #Saving results to the netcdf
        if k%30 == 0:
            print('fileName', fileName)
            ds = nc.Dataset(fileName,'a')
            for j, name in enumerate(model.names):
                print(j, name)
                ds[name][month+1,:,:] = np.ma.masked_array(derivArray[j])
            ds.close()
        else: 
            continue

    return derivArray

def sortData(dateBeginning, dateEnd,lenDataArr):
    datetimeBeginning = datetime.datetime.strptime(dateBeginning, '%Y-%m-%d %H:%M:%S')
    datetimeEnd = datetime.datetime.strptime(dateEnd, '%Y-%m-%d %H:%M:%S')
    #firstdayNbr = (datetimeBeginning - datetime.datetime(datetimeBeginning.year, 1, 1, 0)).days
    firstdayNbr = (datetimeBeginning - datetime.datetime(datetimeBeginning.year, datetimeBeginning.month, datetimeBeginning.day, 0)).days
    #lastdayNbr = (datetimeEnd - datetime.datetime(datetimeBeginning.year, 1, 1, 0)).days #La date de debut doit etre datetimebeginning
    lastdayNbr = (datetimeEnd - datetime.datetime(datetimeBeginning.year, datetimeBeginning.month, datetimeBeginning.day, 0)).days

    if firstdayNbr<lastdayNbr:
        return np.arange(firstdayNbr,lastdayNbr)
    elif firstdayNbr==lastdayNbr:
        return np.arange(0,lenDataArr)
    else:
        return np.concatenate([np.arange(firstdayNbr,lenDataArr),np.arange(0,lastdayNbr)])

if __name__=="__main__":

    t0 = time.time()
    print("Start")
    zone = 'IBI'
    year = 2020
    input_args = {
        'zone' : zone,
        'file_adress' : '..\Data_merg\{zone}\{param}\{param}{zone}modelNetCDF2020-09to2021-06.nc',
        'dataRef' : pd.read_csv('.\dataCmd.csv', delimiter=';'),
        'paramNames' : ['Temperature', 'Chlorophyll-a','northward_Water_current', 'eastward_Water_current']
    }
    ShellfishData = open_data_input(**input_args)
    #print(ShellfishData)

    ### get the copernicus grid and mask
    sim_area = {
        #'longitude': (None, None),
        #'latitude': (None,None),
        'longitude': (-4, -3),
        'latitude': (48.5, 49),
        #'time_index': 0,
        'depth': 0
    }
    dateBeginning = '2020-09-01 00:00:00'
    dateEnd = '2021-06-30 00:00:00'

    dt = 1
    firstday = datetime.datetime.strptime(dateBeginning, '%Y-%m-%d %H:%M:%S')

    longitudes, _ = ShellfishData.parameterData['Temperature'].getVariable('longitude', **sim_area)
    latitudes, _ = ShellfishData.parameterData['Temperature'].getVariable('latitude', **sim_area)

    dataTemp = ShellfishData.parameterData['Temperature'].getVariable( **sim_area)[0]
    dataCHL = ShellfishData.parameterData['Chlorophyll-a'].getVariable( **sim_area)[0]
    dataNwc = ShellfishData.parameterData['northward_Water_current'].getVariable(**sim_area)[0]
    dataEwc = ShellfishData.parameterData['eastward_Water_current'].getVariable(**sim_area)[0]

    sortedList = sortData(dateBeginning, dateEnd, len(dataCHL))

    mask1 = ShellfishData.parameterData['Temperature'].getVariable(**sim_area)[0].mask
    mask2 = ShellfishData.parameterData['eastward_Water_current'].getVariable(**sim_area)[0].mask
    mask3 = ShellfishData.parameterData['northward_Water_current'].getVariable(**sim_area)[0].mask
    mask = np.logical_or(mask1, np.logical_or(mask2, mask3))[0] #just to take one time and get a matrix xith lat and long

    #Import Shellfish parameters from Json
    Shellfish_params = "user1_M_edulis_IBI_12-08-2022.json"
    json_data = import_json(Shellfish_params)

    #Initialize the shellfish model
    model = SF_model_scipy(json_data["parameters"])
    fileName = os.path.join('..','Results','simulations','monthly',"monthly_simulations.nc")
    initialize_result(fileName, np.array(range(1,11)), latitudes, longitudes, model.names, mask) #10 MOIS POUR CE TEST

    run_scenario_a_monthly(fileName, year, firstday, sortedList, dt, dataEwc, dataNwc, dataTemp, dataCHL, model)

    print("Time of processing", (time.time()-t0))
