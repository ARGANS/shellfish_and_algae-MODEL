import copy
import os
import time
import netCDF4 as nc
import numpy as np
import pandas as pd
import numpy.ma as ma
from scipy import interpolate
from scipy.sparse import dia_matrix
from scipy.sparse import spdiags
import datetime
from dateutil.relativedelta import *

# extract the data value at depth in the merged files (all the daily data merged in one file)
from advectionPrototype.saveAsTiff import saveAsTiff, giveMetadata
from dataread.launch_model import MA_model_scipy
from dataread.make_runs import open_data_input
from dataread.read_netcdf import extractVarSubset , AllData
from dataread.utils import import_json






#this function decentralize the u speeds
def u2d_cgrid_cur(var): 
    # var: 3D array ( time, lat, lon)
    varcg = np.zeros((var.shape[0], var.shape[1], var.shape[2] + 1))
    for i in range(var.shape[1]):
        for j in range(var.shape[2] - 1):
            if ((not ma.is_masked(var[0, i, j])) & (not ma.is_masked(var[0, i, j + 1]))):
                varcg[:, i, j] = (var[:, i, j] + var[:, i, j + 1]) / 2
    return varcg

#this function decentralize the v speeds
def v2d_cgrid_cur(var):
    varcg = np.zeros((var.shape[0], var.shape[1] + 1, var.shape[2]))
    for i in range(var.shape[1] - 1):
        for j in range(var.shape[2]):
            if ((not ma.is_masked(var[0, i, j])) & (not ma.is_masked(var[0, i + 1, j]))):
                varcg[:, i, j] = (var[:, i, j] + var[:, i + 1, j]) / 2
    return varcg

def degrees_to_meters(lonDist, latDist, refLat):
    """Converts degrees of latDist,lonDist to meters assuming that the latitude
    stays near refLat.
    """
    lat_degree = 111000 # conversion of latitude degrees to meters

    Dx = lonDist * lat_degree * np.cos(np.deg2rad(refLat))
    Dy = latDist * lat_degree

    return Dx, Dy

#return the CFL number, cx and cy
def giveCFL(dx, dy, dt, Ewc, Nwc):
    Cx = np.abs(Ewc[:, 1:] * dt / dx).flatten()
    Cy = np.abs(Nwc[1:, :] * dt / dy).flatten()
    return np.maximum(Cx, Cy), Cx, Cy

#take as input a mask (tuple), and return the position of the masked values in a vector (dim x * dim y)
def createListNan(maskPosition, nbry):
    listNan = []
    for k in range(len(maskPosition[0])):
        i, j = maskPosition[0][k], maskPosition[1][k]
        n = j + nbry * i
        listNan.append(n)
    return listNan

#return the nan or masked data position
def findNan(mask):
    # mask: 2D array (lat, lon) of booleans, typically the mask of a maskedArray

    nbry, nbrx = mask.shape

    nanPositions = np.empty((5,5), dtype=object)
    # examples:
    #   - nanPositions[0,0] = mask
    #   - nanPositions[0,1] = eastNan
    #   - nanPositions[0,-1] = westNan
    #   - nanPositions[2,0] = up2Nan / north2Nan
    #   - nanPositions[1,-1] = upWestNan / northWestNan

    for iRel in [-2, -1, 0, 1, 2]: # along latitude
        for jRel in [-2, -1, 0, 1, 2]: # along longitude
            # Matrix to shift mask along E-W by -jRel
            shiftMatE = spdiags([1]*nbrx, -jRel, nbrx, nbrx).todense()
            # Matrix to shift mask along N-S by -iRel
            shiftMatN = spdiags([1]*nbry, iRel, nbry, nbry).todense()
            nanPositions[iRel, jRel] = np.flatnonzero(shiftMatN.dot(mask).dot(shiftMatE))

    return nanPositions


#creates the matrix to compute the flux in u direction
def createMatE(decenteredEwc, nanLists, alpha1, alpha2):
    nbrx, nbry = decenteredEwc[:, 1:].shape

    offset = np.array([0, -1, 1, -2, 2])
    uGreater0 = ((decenteredEwc[:, 1:] > 0) * 1).flatten()
    termA = (1 - 2 * alpha1 + alpha2) * uGreater0 - (1 - uGreater0) * (2 * alpha1 + alpha2 - 1)
    termB = (alpha1 - 2 * alpha2 - 1) * uGreater0 - alpha1 * (1 - uGreater0)
    termC = alpha1 * uGreater0 + (1 - alpha1 - 2 * alpha2) * (1 - uGreater0)
    termD = uGreater0 * alpha2
    termE = alpha2 * (1 - uGreater0)

    termB[::nbry] = 0
    termC[nbry - 1::nbry] = 0

    termA[nanLists[0, -1]] = 0
    termB[nanLists[0, -1]] = 0

    termA[nanLists[0, 1]] = 0
    termC[nanLists[0, 1]] = 0

    termA[nanLists[0, -2]] = 0
    termB[nanLists[0, -2]] = 0

    termA[nanLists[0, 2]] = 0
    termC[nanLists[0, 2]] = 0

    data = np.zeros((5, nbrx * nbry))
    data[0, :] = termA
    data[1, :-1] = termB[1:]
    data[2, 1:] = termC[:-1]
    data[3, :-2] = termD[2:]
    data[4, 2:] = termE[:-2]
    Mat = dia_matrix((data, offset), shape=(nbrx * nbry, nbrx * nbry))

    return Mat

#creates the matrix to compute the flux in v direction
def createMatN(decenteredNwc, nanLists, alpha1, alpha2):
    nbrx, nbry = decenteredNwc[1:, :].shape

    offset = np.array([0, -nbry, nbry, -2 * nbry, 2 * nbry])
    vGreater0 = ((decenteredNwc[1:, :] > 0) * 1).flatten()
    termA = (1 - 2 * alpha1 + alpha2) * vGreater0 - (1 - vGreater0) * (2 * alpha1 + alpha2 - 1)
    termB = (alpha1 - 2 * alpha2 - 1) * vGreater0 - alpha1 * (1 - vGreater0)
    termC = alpha1 * vGreater0 + (1 - alpha1 - 2 * alpha2) * (1 - vGreater0)
    termD = vGreater0 * alpha2
    termE = alpha2 * (1 - vGreater0)

    termA[nanLists[-1, 0]] = 0
    termB[nanLists[-1, 0]] = 0

    termA[nanLists[1, 0]] = 0
    termC[nanLists[1, 0]] = 0

    termA[nanLists[-2, 0]] = 0
    termB[nanLists[-2, 0]] = 0

    termA[nanLists[2, 0]] = 0
    termC[nanLists[2, 0]] = 0

    data = np.zeros((5, nbrx * nbry))
    data[0] = termA
    data[1, :-nbry] = termB[nbry:]
    data[2, nbry:] = termC[:-nbry]
    data[3, :-2 * nbry] = termD[2 * nbry:]
    data[4, 2 * nbry:] = termE[:-2 * nbry]
    Mat = dia_matrix((data, offset), shape=(nbrx * nbry, nbrx * nbry))

    return Mat

#returns the space step
def giveDxDy(latitudes, longitudes):
    Dx = np.zeros((len(latitudes), len(longitudes)))
    Dy = np.zeros((len(latitudes), len(longitudes)))
    latRef = np.zeros((len(latitudes), len(longitudes)))

    Dx[:-1, :] = np.diff(latitudes)[np.newaxis].T
    Dx[-1, :] = Dx[-2, :]

    Dy[:, :-1] = np.diff(longitudes)
    Dy[:, -1] = Dy[:, -2]

    latRef[:,:] = latitudes[np.newaxis].T

    DxMeter, DyMeter = degrees_to_meters(Dx, Dy, latRef)
    return DxMeter, DyMeter

#return the variation of each quantities in the algae model
def giveEpsilon(day, temp, NH4, NO3, N_s, N_f, D, PAR, latRef, model, dt, Zmix):
    data_in = {
        'SST': temp,
        'PAR': PAR,
        'PO4_ext': 50,
        'K_d': 0.1,
        't_z': Zmix
    }
    y = dict(zip(["NH4", "NO3", "N_s", "N_f", "D"], [NH4, NO3, N_s, N_f, D]))
    return model.hadley_advection(day, y, data_in, dt, latRef)

#get in the dataLine comming from dataCmd.csv the resolution of the data
def giveResol(dataLine):
    resolString = dataLine.iloc[-1]["resolution"]
    splitedString = resolString.split('*')
    if len(splitedString[0].split('km'))==2:
        return float(splitedString[0].split('km')[0])*1e3, float(splitedString[1].split('km')[0])*1e3, True
    else:
        return float(splitedString[0])*1e-2, float(splitedString[1])*1e-2, False

def prepareDerivArrayScenC(derivArray,nanLists):
    for i in range(len(derivArray)):
        for nanL in nanLists.flatten():
            derivArrayLine = derivArray[i].flatten()
            derivArrayLine[nanL] = 1e-10
            derivArray[i] = derivArrayLine.reshape(grid_shape)
    return derivArray

def sortPAR(dateBeginning, dataArr):
    datetimeBeginning = datetime.datetime.strptime(dateBeginning, '%Y-%m-%d %H:%M:%S')
    firstmonthNbr = datetimeBeginning.month
    dataArr=dataArr.filled(fill_value=-999)
    dataArrShape = np.shape(dataArr)
    newdataArr = np.zeros(dataArrShape)
    newdataArr[:(dataArrShape[0]-firstmonthNbr)] = dataArr[firstmonthNbr:]
    newdataArr[(dataArrShape[0]-firstmonthNbr):] = dataArr[:firstmonthNbr]
    return ma.masked_outside(newdataArr, -1e-2, 1e4)

def sortData(dateBeginning, dateEnd,lenDataArr):
    datetimeBeginning = datetime.datetime.strptime(dateBeginning, '%Y-%m-%d %H:%M:%S')
    datetimeEnd = datetime.datetime.strptime(dateEnd, '%Y-%m-%d %H:%M:%S')
    firstdayNbr = (datetimeBeginning - datetime.datetime(datetimeBeginning.year, 1, 1, 0)).days
    lastdayNbr = (datetimeEnd - datetime.datetime(datetimeEnd.year, 1, 1, 0)).days
    if firstdayNbr<lastdayNbr:
        return np.arange(firstdayNbr,lastdayNbr)
    elif firstdayNbr==lastdayNbr:
        return np.arange(0,lenDataArr)
    else:
        return np.concatenate([np.arange(firstdayNbr,lenDataArr),np.arange(0,lastdayNbr)])

def run_simulation(out_file_name: str, model_json:dict, input_data: AllData):

    # Parse the input json info
    parms_run = list(model_json['parameters']['run'].values())[0]['parameters']
    parms_farm = list(model_json['parameters']['farm'].values())[0]['parameters']
    parms_harvest = list(model_json['parameters']['harvest'].values())[0]['parameters']
    harvest_type = list(model_json['parameters']['harvest'].keys())[0]

    year = model_json['dataset_parameters']['year']

    # Define beginning and end times of the simulation
    startDate = datetime.datetime(year, parms_harvest['deployment_month'], 1)
    if harvest_type == "Winter_growth":
        endDate = datetime.datetime(year + 1, parms_harvest['harvesting_month'], 1) + relativedelta(months=1) #first day of the next month at midnight
    else:
        endDate = datetime.datetime(year, parms_harvest['harvesting_month'], 1) + relativedelta(months=1) #first day of the next month at midnight

    # Data import information, except for the time
    data_kwargs = {
                'longitude': (-4, -2),
                'latitude': (49, 50), #TODO: get from json
                "depth": (0, (1 + parms_run['Von_Karman']) * parms_farm["z"]),
                "averagingDims": ("depth",)
                }

    # Time axis of each variable, stored in a dict
    time_axes, _ = input_data.getData(variable='time')

    # Initialize the working data at the start date, store the index used to detect when it will change.
    nearest_time_i = {}
    working_data = {}
    for par_name, par_data in input_data.parameterData.items():
        nearest_time_i[par_name] = iNearest(startDate, time_axes[par_name])
        working_data[par_name], _ = par_data.getVariable(time_index=nearest_time_i[par_name], **data_kwargs) #TODO: are dims really laways lat,lon ?

    grid_shape = working_data['Nitrate'].shape #the shape should be the same for all parameters

    # Initialize the model variables
    state_vars = {
        'cNO3': np.zeros(grid_shape),
        'cNH4': np.zeros(grid_shape),
        'N_s': np.zeros(grid_shape),
        'N_f': np.ones(grid_shape) * parms_harvest['deployment_Nf'],
        'D': np.zeros(grid_shape),
        'totalNH4deficit': np.zeros(grid_shape),
        'totalNO3deficit': np.zeros(grid_shape)
    }

    dt = 1/72 # TODO: make into parameter in json
    discr = 1/dt

    longitudes, _ = algaeData.parameterData['Nitrate'].getVariable('longitude', **data_kwargs)
    latitudes, _ = algaeData.parameterData['Nitrate'].getVariable('latitude', **data_kwargs)
    dxMeter, dyMeter = giveDxDy(latitudes, longitudes)

    # Simulation loop
    sim_date = startDate
    while sim_date < endDate:

        # Alter the date after new year in case of winter growth
        if harvest_type == "Winter_growth":
            data_date = sim_date.replace(year = year)
        else:
            data_date = sim_date

        # For each dataset, if the nearest i has changed, update the working data
        for par_name, par_data in algaeData.parameterData.items():
            new_i = iNearest(data_date, time_axes[par_name])
            if new_i != nearest_time_i[par_name]:
                nearest_time_i[par_name] = new_i
                working_data[par_name], _ = par_data.getVariable(time_index=new_i, **data_kwargs)

        update_model(state_vars=state_vars, working_data=working_data, dt=dt)

        sim_date += datetime.timedelta(days = dt)


def update_model(state_vars: dict, working_data: dict, dt):

    pass


#apply the quickest scheme
def quickest(dyMeter, dxlist, dt, Ewc, Nwc, latRef, dataNO3, dataNH4, dataTemp, dataPAR, Ks, firstday,
             model, Zmix, scenC, sortedList):
    #we initate the variables
    discr = 1/dt
    nbrStep = int(len(sortedList) * discr)
    cNO3 = np.zeros(np.shape(dataNO3[0]))
    cNH4 = np.zeros(np.shape(dataNH4[0]))
    N_s = np.ones(np.shape(dataNH4[0])) * 1000
    N_f = np.ones(np.shape(dataNH4[0])) * 1000
    D = np.ones(np.shape(dataNH4[0])) * 0.1

    totalNH4deficit = np.zeros(np.shape(dataNH4[0]))
    totalNO3deficit = np.zeros(np.shape(dataNO3[0]))
    grid_shape = np.shape(dataNO3[0])
    mask = dataNO3.mask[0, :, :]
    init = time.time()
    nanLists = findNan(mask)

    daylist = [firstday]

    maxCFL = 0
    dx = dxlist.flatten()
    dy = dyMeter.flatten()
    #for each time step
    for k in range(nbrStep):
        daylist += [firstday + relativedelta(minutes=int(k * 24 * 60 / discr))]
        #we compute the day, hour and month number
        dayNbr = sortedList[k // int(discr)]
        hNbr = dayNbr #if we don't use hourly speeds, hNbr = dayNbr
        month = dayNbr // int(30.5*discr)
        cNO3[mask] = 1e-5
        cNH4[mask] = 1e-5
        D[mask] = 1e-5
        N_f[mask] = 1e-5
        N_s[mask] = 1e-5
        #we compute the CFL
        CFL, cx, cy = giveCFL(dxlist, dyMeter, dt, Ewc[hNbr], Nwc[hNbr])
        alpha1 = (1 / 6) * (1 - CFL) * (2 - CFL)
        alpha2 = (1 / 6) * (1 - CFL) * (1 + CFL)
        if np.max(CFL)>maxCFL:
            maxCFL = np.max(CFL)
        matE = createMatE(Ewc[hNbr], nanLists, alpha1, alpha2)
        matN = createMatN(Nwc[hNbr], nanLists, alpha1, alpha2)

        days = (daylist[-1] - datetime.datetime(daylist[-1].year, 1, 1,
                                                0)).days  # we get the number of the day in the year

        oldNO3Cline = cNO3.flatten()
        oldNH4Cline = cNH4.flatten()
        oldDline = D.flatten()
        no3Hatu = matE.dot(oldNO3Cline)
        no3Hatv = matN.dot(oldNO3Cline)
        nh4Hatu = matE.dot(oldNH4Cline)
        nh4Hatv = matN.dot(oldNH4Cline)
        dHatu = matE.dot(oldDline)
        dHatv = matN.dot(oldDline)
        #we compute the advection terms
        advNO3 = - cy * no3Hatv - cx * no3Hatu + (dt / dy) * Ks * matN.dot(no3Hatv) + (
                    dt / dx) * Ks * matE.dot(no3Hatu)
        advNH4 = - cy * nh4Hatv - cx * nh4Hatu + (dt / dy) * Ks * matN.dot(nh4Hatv) + (
                    dt / dx) * Ks * matE.dot(nh4Hatu)
        advD = - cy * dHatv - cx * dHatu + (dt / dy) * Ks * matN.dot(dHatv) + (
                    dt / dx) * Ks * matE.dot(dHatu)

        CNO3line = oldNO3Cline + advNO3
        #cNO3 = np.minimum(CNO3line.reshape(grid_shape), 0)
        cNO3 = CNO3line.reshape(grid_shape)

        CNH4line = oldNH4Cline + advNH4
        #cNH4 = np.minimum(CNH4line.reshape(grid_shape), 0)
        cNH4 = CNH4line.reshape(grid_shape)

        CDline = oldDline + advD
        D = np.maximum(CDline.reshape(grid_shape),1e-4)

        derivArray = giveEpsilon(days, dataTemp[dayNbr], cNH4 + dataNH4[dayNbr], cNO3 + dataNO3[dayNbr], N_s, N_f, D,
                                 dataPAR[month], latRef, model, dt, Zmix)
        if scenC:
            derivArray = prepareDerivArrayScenC(derivArray, nanLists)

        CNO3line = CNO3line + derivArray[1].flatten() * dt
        #cNO3 = np.minimum(CNO3line.reshape(grid_shape), 0)
        cNO3 = CNO3line.reshape(grid_shape)

        totalNO3deficit += derivArray[1]

        CNH4line = CNH4line + derivArray[0].flatten() * dt
        #cNH4 = np.minimum(CNH4line.reshape(grid_shape), 0)
        cNH4 = CNH4line.reshape(grid_shape)

        totalNH4deficit += derivArray[0]

        CDline = CDline + derivArray[4].flatten() * dt
        D = np.maximum(CDline.reshape(grid_shape),1e-4)

        N_s += derivArray[2] * dt

        N_f += derivArray[3] * dt

        N_f = np.maximum(N_f,-1e-6)
        N_s = np.maximum(N_s, -1e-6)

        if k % int(30 * discr) == 0:
            print(k / discr)
            print(time.time() - init)
            print('Max CFL',maxCFL)

    NO3field, NH4field, D, N_f, N_s = cNO3 + dataNO3[dayNbr - 1], cNH4 + dataNH4[dayNbr - 1], D, N_f, N_s
    NO3field[mask] = np.nan
    NH4field[mask] = np.nan
    D[mask] = np.nan
    N_f[mask] = np.nan
    N_s[mask] = np.nan
    totalNH4deficit[mask] = np.nan
    totalNO3deficit[mask] = np.nan
    return NO3field, NH4field, D, N_f, N_s, totalNH4deficit, totalNO3deficit


if __name__ == "__main__":
    zone = "NWS"
    depth = 0
    PAR_year = 2020
    getAmmonium = False

    dateBeginning = '2020-09-01 00:00:00'
    dateEnd = '2020-04-30 00:00:00'

    scenC = False

    # discr = 144
    discr = 72  # Baltic, NWS, IBI
    #discr = 48 #BS
    dt = 1 / discr

    Ks = 1e-3
    Ks *= 60 * 60 * 24

    CPlat, CPlon = 153, 56

    model_params = "p:/Aquaculture/shellfish_and_algae-MODEL/macroalgae/macroalgae_model_parameters_input.json"
    json_data = import_json(model_params)

    parms_run = list(json_data['parameters']['run'].values())[0]['parameters']
    parms_farm = list(json_data['parameters']['farm'].values())[0]['parameters']
    parms_harvest = list(json_data['parameters']['harvest'].values())[0]['parameters']
    harvest_type = list(json_data['parameters']['harvest'].keys())[0]

    model = MA_model_scipy(json_data['parameters'])

    if getAmmonium:
        paramNames = ['Nitrate', 'northward_Water_current', 'Ammonium', 'eastward_Water_current',
                      'Temperature']  # IBI, Baltic and MED
    else:
        paramNames = ['Nitrate', 'northward_Water_current', 'eastward_Water_current', 'Temperature']  # NWS and BS

    firstday = datetime.datetime.strptime(dateBeginning, '%Y-%m-%d %H:%M:%S')

    dataCmdpath = './../global/dataCmd.csv'
    mergedFilepath = 'D:/Profils/mjaouen/Documents/alternance/EASME/data/{zone}/'.format(zone=zone)

    input_args = {
        'zone': zone,
        'file_adress': 'D:/Profils/mjaouen/Documents/alternance/EASME/data/{zone}/merged_{param}_{zone}.nc',
        'dataRef': pd.read_csv(dataCmdpath, delimiter=';'),
        'paramNames': paramNames,
        'frequency': 'daily'
        # 'with_PAR': PAR_year,
        # 'PAR_file': 'D:/Profils/mjaouen/Documents/alternance/EASME/data/{zone}/PAR_{zone}NewGrid.nc'.format(zone = zone)
    }

    dict_to_AllData = open_data_input(**input_args)

    dict_to_AllData['PAR'] = {
        'file_name': 'D:/Profils/mjaouen/Documents/alternance/EASME/data/{zone}/PAR_{zone}NewGrid.nc'.format(zone=zone),
        'variable_name': 'par',
        'latitude_name': 'lat',
        'longitude_name': 'lon',
        'time_name': 'time',
        'depth_name': 'depth',
        'unit_conversion': 11.574,
        'time_zero': datetime.datetime(PAR_year, 1, 1),
        'time_step': datetime.timedelta(days=1)
    }

    sim_area = {
        'longitude': (-4, -2),
        'latitude': (49, 50),
        #'depth': 0
        'depth': (0, (1 + parms_run['Von_Karman']) * parms_farm["z"]),
        'averagingDims': ('depth',),
        #'weighted': False
    }

    Zmix = (1 + parms_run['Von_Karman']) * parms_farm["z"]

    ### Initialize the netcdf reading interface
    algaeData = AllData(dict_to_AllData)

    dataNO3 = algaeData.parameterData['Nitrate'].getVariable(**sim_area)[0]
    if getAmmonium:
        dataNH4 = algaeData.parameterData['Ammonium'].getVariable(**sim_area)[0]
    else:
        dataNH4 = dataNO3 * 0.1
    dataTemp = algaeData.parameterData['Temperature'].getVariable(**sim_area)[0]
    dataNwc = algaeData.parameterData['northward_Water_current'].getVariable(**sim_area)[0]
    dataEwc = algaeData.parameterData['eastward_Water_current'].getVariable(**sim_area)[0]
    dataPAR = algaeData.parameterData['PAR'].getVariable(**sim_area)[0]

    print(np.shape(dataNO3))
    print(np.shape(dataNH4))
    print(np.shape(dataTemp))
    print(np.shape(dataNwc))
    print(np.shape(dataPAR))

    sortedList = sortData(dateBeginning, dateEnd, len(dataNO3))

    dataNH4 = ma.masked_outside(dataNH4, -1e4, 1e4)
    dataNO3 = ma.masked_outside(dataNO3, -1e4, 1e4)
    dataTemp = ma.masked_outside(dataTemp, -1e4, 1e4)
    dataPAR = ma.masked_outside(dataPAR, -1e-2, 1e4)
    dataPAR = dataPAR.filled(fill_value=8)

    fileV = 'D:/Profils/mjaouen/Documents/alternance/EASME/data/{zone}/merged_northward_Water_current_{zone}.nc'.format(zone=zone)
    dataBaseNwc = nc.Dataset(fileV)

    dataFin = pd.read_csv('./../global/dataCmd.csv', ';')

    nwcDataLine = dataFin.loc[(dataFin["Parameter"] == 'Nitrate') & (dataFin["Place"] == zone)]
    nwcdataName = nwcDataLine.iloc[-1]["variable"]  # we find the data name in the dataset
    resx, resy, km = giveResol(nwcDataLine)
    print(resx, resy, km)

    longitudes, _ = algaeData.parameterData['Temperature'].getVariable('longitude', **sim_area)
    latitudes, _ = algaeData.parameterData['Temperature'].getVariable('latitude', **sim_area)

    '''    print(type(dataEwc))
    dataEwc[np.where(np.abs(dataEwc) >1e8)]=0
    dataNwc[np.where(np.abs(dataNwc) >1e8)] = 0'''

    latRef = np.ones((np.shape(dataEwc[0])[1], np.shape(dataEwc[0])[0])) * latitudes

    decenteredEwc = u2d_cgrid_cur(dataEwc)
    decenteredNwc = v2d_cgrid_cur(dataNwc)

    dxlist, dyMeter = giveDxDy(latitudes,longitudes)

    dylist = dyMeter * np.ones(np.shape(dxlist))

    maxCFL = 0
    grid_shape = np.shape(dataNO3[0])
    for i in range(len(dataEwc)):
        CFL, cx, cy = giveCFL(dxlist, dyMeter, dt, decenteredEwc[i], decenteredNwc[i])
        if np.max(CFL) > maxCFL:
            maxCFL = np.max(CFL)
    print('maxCFL: '+ str(maxCFL))
    #xsize, ysize, ulx, uly, xres, yres = giveMetadata(latitudes, longitudes)
    #saveAsTiff(dataNO3[0], xsize, ysize, ulx, uly, xres, yres, "I:/work-he/apps/safi/data/IBI/test.tiff")
    NO3field, NH4field, D, N_f, N_s, totalNH4deficit, totalNO3deficit = quickest(dyMeter, dxlist, dt,
                                                                                 decenteredEwc, decenteredNwc,
                                                                                 latRef.T, dataNO3, dataNH4,
                                                                                 dataTemp, dataPAR, Ks, firstday, model,
                                                                                 Zmix,scenC,sortedList)


    print(NO3field.shape)
