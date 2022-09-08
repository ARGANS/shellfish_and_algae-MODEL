import copy
import os
import stat
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
import sys

# extract the data value at depth in the merged files (all the daily data merged in one file)
sys.path.append('p:/Aquaculture/shellfish_and_algae-MODEL/dataread/src/')
#from saveAsTiff import saveAsTiff, giveMetadata
try:
    from launch_model import MA_model_scipy
    from make_runs import open_data_input, initialize_result
    from read_netcdf import AllData, iNearest
    from utils import import_json
except ImportError:
    from dataread.launch_model import MA_model_scipy
    from dataread.make_runs import open_data_input, initialize_result
    from dataread.read_netcdf import AllData, iNearest
    from dataread.utils import import_json


class Resampler:
    def __init__(self, dxRatio, dyRatio, init_grid_shape):
        self.dxRatio = dxRatio
        self.dyRatio = dyRatio
        nbry, nbrx = init_grid_shape
        self.grid_shape = (int(nbry/dyRatio),int(nbrx/dxRatio))

    def findElmt(self, i, j):
        newi = np.floor(i * self.dyRatio).astype(int)
        newj = np.floor(j * self.dxRatio).astype(int)
        return newi, newj

    def giveNewMatCoor(self):
        rowCoor = np.zeros(self.grid_shape[0] * self.grid_shape[1])
        columnCoor = np.zeros(self.grid_shape[0] * self.grid_shape[1])
        ival, jval = self.findElmt(np.arange(self.grid_shape[0]),
                                   np.arange(self.grid_shape[1]))
        #for each row value
        for k in range(self.grid_shape[0]):
            rowCoor[k * self.grid_shape[1]:(k + 1) * self.grid_shape[1]] = ival[k]
            columnCoor[k * self.grid_shape[1]:(k + 1) * self.grid_shape[1]] = jval
        return rowCoor.astype(int), columnCoor.astype(int)

    def resampleLonLat(self,lon,lat):
        lat_id, lon_id = self.findElmt(np.arange(self.grid_shape[0]), np.arange(self.grid_shape[1]))
        return lon[lon_id], lat[lat_id]

    def resampleData(self, dataArray):
        rowCoor, columnCoor = self.giveNewMatCoor()
        resampledArray = dataArray[rowCoor, columnCoor].reshape(self.grid_shape)
        return resampledArray

#this function decentralize the u speeds
def u2d_cgrid_cur(var):
    # var: 3D array ( time, lat, lon)
    varcg = np.zeros((var.shape[0], var.shape[1], var.shape[2] + 1))
    for i in range(var.shape[1]):
        for j in range(var.shape[2] - 1):
            if ((not ma.is_masked(var[0, i, j])) & (not ma.is_masked(var[0, i, j + 1]))):
                varcg[:, i, j+1] = (var[:, i, j] + var[:, i, j + 1]) / 2
    return varcg

#this function decentralize the v speeds
def v2d_cgrid_cur(var):
    varcg = np.zeros((var.shape[0], var.shape[1] + 1, var.shape[2]))
    for i in range(var.shape[1] - 1):
        for j in range(var.shape[2]):
            if ((not ma.is_masked(var[0, i, j])) & (not ma.is_masked(var[0, i + 1, j]))):
                varcg[:, i+1, j] = (var[:, i, j] + var[:, i + 1, j]) / 2
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

    Dx[:, :-1] = np.diff(longitudes)
    Dx[:, -1] = Dx[:, -2]

    Dy[:-1, :] = np.diff(latitudes)[np.newaxis].T
    Dy[-1, :] = Dy[-2, :]

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

def prepareScenC(nitrogenArray,nanLists, grid_shape):
    for iRel in [-1, 0, 1]:  # along latitude
        for jRel in [-1, 0, 1]:  # along longitude
            nanL = nanLists[iRel,jRel]
            nitArrayLine = nitrogenArray.flatten()
            nitArrayLine[nanL] = 0
            nitrogenArray = nitArrayLine.reshape(grid_shape)
    return nitrogenArray

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

    t_init = time.time()

    # Parse the input json info
    parms_run = list(model_json['parameters']['run'].values())[0]['parameters']
    parms_farm = list(model_json['parameters']['farm'].values())[0]['parameters']
    parms_harvest = list(model_json['parameters']['harvest'].values())[0]['parameters']
    harvest_type = list(model_json['parameters']['harvest'].keys())[0]

    scenC = (model_json['metadata']['scenario']=="C")

    year = int(model_json['dataset_parameters']['year'])

    model = MA_model_scipy(model_json['parameters'])

    # Define beginning and end times of the simulation
    startDate = datetime.datetime(year, int(parms_harvest['deployment_month']), 1)
    if harvest_type == "Winter_growth":
        endDate = datetime.datetime(year + 1, int(parms_harvest['harvesting_month']), 1) + relativedelta(months=1) #first day of the next month at midnight
    else:
        endDate = datetime.datetime(year, int(parms_harvest['harvesting_month']), 1) + relativedelta(months=1) #first day of the next month at midnight

    # Data import information, except for the time
    data_kwargs = {
                'longitude': (parms_run["min_lon"], parms_run["max_lon"]),
                'latitude': (parms_run["min_lat"], parms_run["max_lat"]),
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

    for par_name, par_data in input_data.parameterData.items():
        print(par_name)
        print(working_data[par_name].shape)

    init_grid_shape = working_data['Nitrate'].shape #the shape should be the same for all parameters
    longitudes, _ = algaeData.parameterData['Nitrate'].getVariable('longitude', **data_kwargs)
    latitudes, _ = algaeData.parameterData['Nitrate'].getVariable('latitude', **data_kwargs)
    dxMeter, dyMeter = giveDxDy(latitudes, longitudes)
    latRef = np.zeros((len(latitudes), len(longitudes)))
    latRef[:, :] = latitudes[np.newaxis].T

    grid_shape = init_grid_shape
    if scenC:
        dxRatio = 1852 / np.mean(dxMeter) #1852 meters = 1 nautical mile
        dyRatio = 1852 / np.mean(dyMeter)
        resa = Resampler(dxRatio,dyRatio,init_grid_shape)
        grid_shape = resa.grid_shape
        for par_name, par_data in input_data.parameterData.items():
            working_data[par_name] = resa.resampleData(working_data[par_name])
        latRef = resa.resampleData(latRef)
        dyMeter = resa.resampleData(dyMeter)
        dxMeter = resa.resampleData(dxMeter)

    for par_name, par_data in input_data.parameterData.items():
        print(par_name)
        print(working_data[par_name].shape)

    mask = working_data['Nitrate'].mask
    nanLists = findNan(mask)

    # Initialize the model variables
    state_vars = {
        'cNO3': np.ma.masked_array(np.zeros(grid_shape), mask),
        'cNH4': np.ma.masked_array(np.zeros(grid_shape), mask),
        'N_s': np.ma.masked_array(np.zeros(grid_shape), mask),
        'N_f': np.ma.masked_array(np.ones(grid_shape) * parms_harvest['deployment_Nf'], mask),
        'D': np.ma.masked_array(np.ones(grid_shape) * parms_run["Detritus"], mask)
    }
    if scenC:
        state_vars['N_s'] = prepareScenC(state_vars['N_s'], nanLists, grid_shape)
        state_vars['N_f'] = prepareScenC(state_vars['N_f'], nanLists, grid_shape)

    working_data["decentered_U"] = np.ma.masked_array(np.zeros((grid_shape[0], grid_shape[1] + 1)))
    working_data["decentered_U"][:, 1:-1] = (working_data['eastward_Water_current'][:, 1:] + working_data['eastward_Water_current'][:, :-1]) / 2
    working_data["decentered_U"][working_data["decentered_U"].mask] = 0

    working_data["decentered_V"] = np.ma.masked_array(np.zeros((grid_shape[0] + 1, grid_shape[1])))
    working_data["decentered_V"][1:-1, :] = (working_data['northward_Water_current'][1:, :] + working_data['northward_Water_current'][:-1, :]) / 2
    working_data["decentered_V"][working_data["decentered_V"].mask] = 0

    dt = 1/72 # days # TODO: make into parameter in json

    Ks = 1e-3 * 60 * 60 * 24 # m2/s

    # Simulation loop
    sim_date = startDate
    while sim_date < endDate:
        print(f'{sim_date}')

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

                if scenC:
                    working_data[par_name] = resa.resampleData(working_data[par_name])
                # Update the centered currents as well
                if par_name == "eastward_Water_current":
                    working_data["decentered_U"][:, 1:-1] = (working_data['eastward_Water_current'][:, 1:] + working_data['eastward_Water_current'][:, :-1]) / 2
                    working_data["decentered_U"][working_data["decentered_U"].mask] = 0
                if par_name == "northward_Water_current":
                    working_data["decentered_V"][1:-1, :] = (working_data['northward_Water_current'][1:, :] + working_data['northward_Water_current'][:-1, :]) / 2
                    working_data["decentered_V"][working_data["decentered_V"].mask] = 0

        if (model_json['metadata']['scenario'] == "A"):
            advection_terms = advection_modelA(state_vars=state_vars, working_data=working_data,
                                               dt=dt, dxMeter=dxMeter, dyMeter=dyMeter,
                                               paramDetritus=parms_run["Detritus"])
        else:
            # Compute the advection terms
            advection_terms = advection_model(state_vars=state_vars, working_data=working_data,
                                              dt=dt, dxMeter=dxMeter, dyMeter=dyMeter, nanLists=nanLists,
                                              Ks=Ks)

        # Apply the advection
        for var_name in state_vars.keys():
            state_vars[var_name] += advection_terms[var_name] * dt


        days = (data_date - datetime.datetime(year, 1, 1)).days # Note: returns an integer only, that is sufficient precision for this

        bgc_terms = bgc_model(state_vars=state_vars, working_data=working_data, dt=dt,
                              model=model, parms_run=parms_run, days=days, latRef=latRef)


        # Testing: value of dt for an adapting time step
        '''
        print(f'    Value of dt_max (minutes): {dt_max*(60*24)}')
        dt_list = [dt_max]
        for var_name in ['N_s', 'N_f']:#state_vars.keys():
            value_array = state_vars[var_name].copy()
            if var_name == 'cNO3':
                value_array += working_data['Nitrate']
            if var_name == 'cNH4':
                value_array += working_data['Ammonium']

            dt_array = - value_array/bgc_terms[var_name]
            dt_array[dt_array <= 0] = np.inf
            dt_pos = np.amin(dt_array)
            print(f'    Value of dt for {var_name} > 0 (minutes): {dt_pos*(60*24)}')

            dt_list.append(dt_pos * 0.8)

        dt = min(dt_list)
        print(f'Value of dt used: {dt*(60*24)} minutes')
        '''

        # Apply the bgc terms
        if scenC:
            for var_name in state_vars.keys():
                bgc_terms[var_name] = prepareScenC(bgc_terms[var_name], nanLists, grid_shape)
        for var_name in state_vars.keys():
            state_vars[var_name] += bgc_terms[var_name] * dt


        # Safety maximum for now
        state_vars['N_s'] = np.maximum(state_vars['N_s'], 1e-6)
        state_vars['N_f'] = np.maximum(state_vars['N_f'], 1e-6)
        state_vars['D'] = np.maximum(state_vars['D'], 1e-4)

        # Debugging
        '''
        do_break = False
        for var_name in ['N_s', 'N_f', 'D']:
            if np.any(state_vars[var_name] < 0):
                print(f"Boke the simulation because {var_name} is negative")
                do_break = True
        if do_break:
            break

        if np.any(state_vars['cNO3'] + working_data['Nitrate'] < 0) or np.any(state_vars['cNH4'] + working_data['Ammonium'] < 0) :
            print("Boke the simulation because the total NO3 or NH4 is negative")
            break
        '''

        # Time dissipation of the signal
        #dissip_t = 3 #days
        #state_vars['cNH4'] = state_vars['cNH4'] * (1 - dt/dissip_t)
        #state_vars['cNO3'] = state_vars['cNO3'] * (1 - dt/dissip_t)
        #state_vars['D'] = state_vars['cNO3'] * (1 - dt/dissip_t) + parms_run["Detritus"] * dt/dissip_t

        sim_date += datetime.timedelta(days = dt)

    if scenC:
        latStep = (latitudes[-1]-latitudes[0])/(grid_shape[0]-1)
        lonStep = (longitudes[-1]-longitudes[0])/(grid_shape[1]-1)

        latitudes = latStep*np.arange(grid_shape[0])+latitudes[0]
        longitudes = lonStep * np.arange(grid_shape[1]) + longitudes[0]

    state_vars['NH4'] = working_data['Ammonium'] + state_vars['cNH4']
    state_vars['NO3'] = working_data['Nitrate'] + state_vars['cNO3']

    # Create output file
    initialize_result(out_file_name, times=[0], latitudes=latitudes, longitudes=longitudes,
                      variableNames=['NH4', 'NO3', 'N_s', 'N_f', 'D'], mask=mask)
    # Write values to file
    ds = nc.Dataset(out_file_name, 'a')
    for name in model.names:
        ds[name][0,:,:] = np.ma.masked_array(state_vars[name], mask)
    ds.close()

    return time.time() - t_init


#return the variation of each quantities in the algae model
def bgc_model(state_vars: dict, working_data: dict, dt, model, parms_run, days, latRef):

    data_in = {
        'SST': working_data['Temperature'],
        'PAR': working_data['par'],
        'PO4_ext': working_data['Phosphate'], #TODO: from working data
        'K_d': parms_run["K_d490"],
        't_z': (1 + parms_run['Von_Karman']) * model._parameters["z"]
    }

    y = {
        'NH4': working_data['Ammonium'] + state_vars['cNH4'],
        'NO3': working_data['Nitrate'] + state_vars['cNO3'],
        'N_s': state_vars['N_s'],
        'N_f': state_vars['N_f'],
        'D': state_vars['D'],
    }

    terms_list = model.hadley_advection(days, y, data_in, dt, latRef)
    all_terms = dict(zip(["cNH4", "cNO3", "N_s", "N_f", "D"], terms_list))

    return all_terms

def advection_modelA(state_vars: dict, working_data: dict, dt, dxMeter: np.array, dyMeter: np.array, paramDetritus):

    mask = working_data['Nitrate'].mask
    grid_shape = working_data['Nitrate'].shape

    uGreater0 = ((working_data["decentered_U"][:, 1:] > 0) * working_data["decentered_U"][:, 1:]).flatten()
    uLower0 = ((working_data["decentered_U"][:, :-1] < 0) * working_data["decentered_U"][:, :-1]).flatten()

    vGreater0 = ((working_data["decentered_V"][1:] > 0) * working_data["decentered_V"][1:]).flatten()
    vLower0 = ((working_data["decentered_V"][:-1] < 0) * working_data["decentered_V"][:-1]).flatten()

    cNO3_line = state_vars['cNO3'].flatten()
    cNH4_line = state_vars['cNH4'].flatten()
    D_line_eps = paramDetritus-state_vars['D'].flatten()
    dx = dxMeter.flatten()
    dy = dyMeter.flatten()

    #we compute the advection terms
    advNO3 = ((dt / dx) * (-uGreater0 + uLower0) + (dt / dy) * (-vGreater0 + vLower0))*cNO3_line
    advNH4 = ((dt / dx) * (-uGreater0 + uLower0) + (dt / dy) * (-vGreater0 + vLower0))*cNH4_line
    advD = -((dt / dx) * (-uGreater0 + uLower0) + (dt / dy) * (-vGreater0 + vLower0))*D_line_eps

    print(np.mean(advNO3),np.mean(advNH4), np.mean(advD))
    # reshape to the grid
    advNO3 = advNO3.reshape(grid_shape)
    advNH4 = advNH4.reshape(grid_shape)
    advD = advD.reshape(grid_shape)

    # reapply masks that were altered
    advNO3.mask = mask
    advNH4.mask = mask
    advD.mask = mask

    all_terms = {
        'cNO3': advNO3 / dt,
        'cNH4': advNH4 / dt,
        'N_s': 0,
        'N_f': 0,
        'D': advD / dt
    }

    return all_terms

def advection_model(state_vars: dict, working_data: dict, dt, dxMeter: np.array, dyMeter: np.array,
                    nanLists: np.array, Ks):

    mask = working_data['Nitrate'].mask
    grid_shape = working_data['Nitrate'].shape

    CFL, cx, cy = giveCFL(dxMeter, dyMeter, dt,
                          working_data["decentered_U"], working_data["decentered_V"])

    alpha1 = (1 / 6) * (1 - CFL) * (2 - CFL)
    alpha2 = (1 / 6) * (1 - CFL) * (1 + CFL)

    matE = createMatE(working_data["decentered_U"], nanLists, alpha1, alpha2)
    matN = createMatN(working_data["decentered_V"], nanLists, alpha1, alpha2)

    cNO3_line = state_vars['cNO3'].flatten()
    cNH4_line = state_vars['cNH4'].flatten()
    D_line = state_vars['D'].flatten()
    dx = dxMeter.flatten()
    dy = dyMeter.flatten()

    no3Hatu = matE.dot(cNO3_line)
    no3Hatv = matN.dot(cNO3_line)
    nh4Hatu = matE.dot(cNH4_line)
    nh4Hatv = matN.dot(cNH4_line)
    dHatu = matE.dot(D_line)
    dHatv = matN.dot(D_line)
    #we compute the advection terms
    advNO3 = - cy*no3Hatv - cx*no3Hatu + (dt / dy) * Ks * matN.dot(no3Hatv) + \
                                         (dt / dx) * Ks * matE.dot(no3Hatu)
    advNH4 = - cy*nh4Hatv - cx*nh4Hatu + (dt / dy) * Ks * matN.dot(nh4Hatv) + \
                                         (dt / dx) * Ks * matE.dot(nh4Hatu)
    advD = - cy*dHatv - cx*dHatu + (dt / dy) * Ks * matN.dot(dHatv) + \
                                   (dt / dx) * Ks * matE.dot(dHatu)

    # reshape to the grid
    advNO3 = advNO3.reshape(grid_shape)
    advNH4 = advNH4.reshape(grid_shape)
    advD = advD.reshape(grid_shape)

    # reapply masks that were altered
    advNO3.mask = mask
    advNH4.mask = mask
    advD.mask = mask

    all_terms = {
        'cNO3': advNO3 / dt,
        'cNH4': advNH4 / dt,
        'N_s': 0,
        'N_f': 0,
        'D': advD / dt
    }

    return all_terms


#apply the quickest scheme
def quickest(dyMeter, dxlist, dt, Ewc, Nwc, latRef, dataNO3, dataNH4, dataTemp, dataPAR, Ks, firstday,
             model, Zmix, scenC, sortedList):
    #we initate the variables
    discr = 1/dt
    nbrStep = int(len(sortedList) * discr)

    grid_shape = np.shape(dataNO3[0])
    mask = dataNO3.mask[0, :, :]
    init = time.time()
    nanLists = findNan(mask)

    print(mask)

    cNO3 = np.ma.masked_array(np.zeros(grid_shape), mask)
    cNH4 = np.ma.masked_array(np.zeros(grid_shape), mask)
    N_s = np.ma.masked_array(np.ones(grid_shape) * 1000, mask)
    N_f = np.ma.masked_array(np.ones(grid_shape) * 1000, mask)
    D = np.ma.masked_array(np.ones(grid_shape) * 0.1, mask)
    totalNH4deficit = np.ma.masked_array(np.zeros(grid_shape), mask)
    totalNO3deficit = np.ma.masked_array(np.zeros(grid_shape), mask)

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
        cNO3 = CNO3line.reshape(grid_shape)

        CNH4line = oldNH4Cline + advNH4
        cNH4 = CNH4line.reshape(grid_shape)

        CDline = oldDline + advD
        D = np.maximum(CDline.reshape(grid_shape), 1e-4)

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

        cNO3.mask = mask
        cNH4.mask = mask
        D.mask = mask
        N_s.mask = mask
        N_f.mask = mask
        totalNO3deficit.mask = mask
        totalNH4deficit.mask = mask

        if k % int(30 * discr) == 0:
            print(k / discr)
            print(time.time() - init)
            print('Max CFL',maxCFL)

    NO3field, NH4field, D, N_f, N_s = cNO3 + dataNO3[dayNbr - 1], cNH4 + dataNH4[dayNbr - 1], D, N_f, N_s

    return NO3field, NH4field, D, N_f, N_s, totalNH4deficit, totalNO3deficit


if __name__ == "__main__":
    zone = "NWS"
    depth = 0
    PAR_year = 2020
    getAmmonium = False

    dateBeginning = '2020-09-01 00:00:00'
    dateEnd = '2020-04-30 00:00:00'

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

    dataCmdpath = 'p:/Aquaculture/shellfish_and_algae-MODEL/global/dataCmd.csv'

    input_args = {
        'zone': zone,
        'file_adress': 'D:/data_scenario_B/{zone}/merged_{param}_{zone}.nc',
        'dataRef': pd.read_csv(dataCmdpath, delimiter=';'),
        'paramNames': paramNames,
        'frequency': 'daily'
    }

    #dict_to_AllData = open_data_input(**input_args)
    #print(dict_to_AllData)
    dict_to_AllData = {
        "Temperature": {
            'file_name': f'D:/data_scenario_B/{zone}/merged_Temperature_{zone}.nc',
            'variable_name': 'thetao',
            'unit_conversion': 1,
            'time_zero': datetime.datetime(1970, 1, 1),
            'time_step': datetime.timedelta(seconds=1)
        },
        "Nitrate": {
            'file_name': f'D:/data_scenario_B/{zone}/merged_Nitrate_{zone}.nc',
            'variable_name': 'no3',
            'unit_conversion': 14,
            'time_zero': datetime.datetime(1970, 1, 1),
            'time_step': datetime.timedelta(seconds=1)
        },
        "eastward_Water_current": {
            'file_name': f'D:/data_scenario_B/{zone}/merged_eastward_Water_current_{zone}.nc',
            'variable_name': 'uo',
            'unit_conversion': 86400,
            'time_zero': datetime.datetime(1970, 1, 1),
            'time_step': datetime.timedelta(seconds=1)
        },
        "northward_Water_current": {
            'file_name': f'D:/data_scenario_B/{zone}/merged_northward_Water_current_{zone}.nc',
            'variable_name': 'vo',
            'unit_conversion': 86400,
            'time_zero': datetime.datetime(1970, 1, 1),
            'time_step': datetime.timedelta(seconds=1)
        }
    }

    dict_to_AllData['Ammonium'] = dict_to_AllData['Nitrate'].copy()
    dict_to_AllData['Ammonium']['unit_conversion'] /= 10

    dict_to_AllData['par'] = {
        'file_name': f'D:/data_scenario_B/{zone}/PAR_NewGrid.nc',
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

    '''
    dataNO3 = algaeData.parameterData['Nitrate'].getVariable(**sim_area)[0]
    if getAmmonium:
        dataNH4 = algaeData.parameterData['Ammonium'].getVariable(**sim_area)[0]
    else:
        dataNH4 = dataNO3 * 0.1
    dataTemp = algaeData.parameterData['Temperature'].getVariable(**sim_area)[0]
    dataNwc = algaeData.parameterData['northward_Water_current'].getVariable(**sim_area)[0]
    dataEwc = algaeData.parameterData['eastward_Water_current'].getVariable(**sim_area)[0]
    dataPAR = algaeData.parameterData['par'].getVariable(**sim_area)[0]

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


    dataFin = pd.read_csv('p:/Aquaculture/shellfish_and_algae-MODEL/global/dataCmd.csv', ';')

    nwcDataLine = dataFin.loc[(dataFin["Parameter"] == 'Nitrate') & (dataFin["Place"] == zone)]
    nwcdataName = nwcDataLine.iloc[-1]["variable"]  # we find the data name in the dataset
    resx, resy, km = giveResol(nwcDataLine)
    print(resx, resy, km)

    longitudes, _ = algaeData.parameterData['Temperature'].getVariable('longitude', **sim_area)
    latitudes, _ = algaeData.parameterData['Temperature'].getVariable('latitude', **sim_area)

    #print(type(dataEwc))
    #dataEwc[np.where(np.abs(dataEwc) >1e8)]=0
    #dataNwc[np.where(np.abs(dataNwc) >1e8)] = 0

    latRef = np.ones((np.shape(dataEwc[0])[1], np.shape(dataEwc[0])[0])) * latitudes

    decenteredEwc = u2d_cgrid_cur(dataEwc)
    decenteredNwc = v2d_cgrid_cur(dataNwc)

    dxlist, dyMeter = giveDxDy(latitudes, longitudes)

    maxCFL = 0
    grid_shape = np.shape(dataNO3[0])
    for i in range(len(dataEwc)):
        CFL, cx, cy = giveCFL(dxlist, dyMeter, dt, decenteredEwc[i], decenteredNwc[i])
        if np.max(CFL) > maxCFL:
            maxCFL = np.max(CFL)
    print('maxCFL: '+ str(maxCFL))
    #xsize, ysize, ulx, uly, xres, yres = giveMetadata(latitudes, longitudes)
    #saveAsTiff(dataNO3[0], xsize, ysize, ulx, uly, xres, yres, "I:/work-he/apps/safi/data/IBI/test.tiff")
    """
    NO3field, NH4field, D, N_f, N_s, totalNH4deficit, totalNO3deficit = quickest(dyMeter, dxlist, dt,
                                                                                 decenteredEwc, decenteredNwc,
                                                                                 latRef.T, dataNO3, dataNH4,
                                                                                 dataTemp, dataPAR, Ks, firstday, model,
                                                                                 Zmix,scenC,sortedList)
    """

    '''
    result = run_simulation(out_file_name=f"D:/data_scenario_B/{zone}/test_2.nc", input_data=algaeData,
                            model_json=json_data)

    print(f"Computation time (seconds): {result}")