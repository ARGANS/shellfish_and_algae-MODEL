import datetime
import numpy as np
import pandas as pd
import time
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy import interpolate
from scipy.sparse import dia_matrix
from scipy.sparse import spdiags
import multiprocessing as mp
import datetime
from dateutil.relativedelta import *
from .read_netcdf import *
from .launch_model import *
from .utils import import_json

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

def degrees_to_meters(lonDist, latDist, refLat):
    """Converts degrees of latDist,lonDist to meters assuming that the latitude
    stays near refLat.
    """
    lat_degree = 111000 # conversion of latitude degrees to meters

    Dx = lonDist * lat_degree * np.cos(np.deg2rad(refLat))
    Dy = latDist * lat_degree

    return Dx, Dy


def meters_to_degrees(Dx, Dy, refLat):
    """Converts longitudinal and meridional distances Dx and Dy to degrees of
    longitude and latitude assuming that the latitude stays near refLat"""

    lat_degree = 111000 # conversion of latitude degrees to meters

    lonDist = Dx / (lat_degree * np.cos(np.deg2rad(refLat)))
    latDist = Dy / lat_degree

    return lonDist, latDist

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
def createMatE_available(decenteredEwc, nanLists):
    nbrx, nbry = decenteredEwc[:, 1:].shape

    offset = np.array([-1, 1])

    uGreater0 = ((decenteredEwc[:, :-1] > 0) * decenteredEwc[:, :-1]).flatten()
    uLower0 = ((decenteredEwc[:, 1:] < 0) * decenteredEwc[:, 1:]).flatten()

    termA = -uGreater0
    termB = uLower0


    termA[::nbry] = 0
    termB[nbry - 1::nbry] = 0

    termA[nanLists[0, -1]] = 0
    termB[nanLists[0, 1]] = 0

    data = np.zeros((2, nbrx * nbry))
    data[0, :-1] = termA[1:]
    data[1, 1:] = termB[:-1]
    Mat = dia_matrix((data, offset), shape=(nbrx * nbry, nbrx * nbry))

    return Mat

#creates the matrix to compute the flux in v direction
def createMatN_available(decenteredNwc, nanLists):
    nbrx, nbry = decenteredNwc[1:, :].shape

    offset = np.array([-nbry, nbry])

    vGreater0 = ((decenteredNwc[:-1] > 0) * decenteredNwc[:-1]).flatten()
    vLower0 = ((decenteredNwc[1:] < 0) * decenteredNwc[1:]).flatten()

    termA = -vGreater0
    termB = vLower0

    termA[nanLists[-1, 0]] = 0

    termB[nanLists[1, 0]] = 0

    data = np.zeros((2, nbrx * nbry))
    data[0, :-nbry] = termA[nbry:]
    data[1, nbry:] = termB[:-nbry]
    Mat = dia_matrix((data, offset), shape=(nbrx * nbry, nbrx * nbry))

    return Mat

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

def prepareScenC(nitrogenArray,nanLists, grid_shape):
    for iRel in [-1, 0, 1]:  # along latitude
        for jRel in [-1, 0, 1]:  # along longitude
            nanL = nanLists[iRel,jRel]
            nitArrayLine = nitrogenArray.flatten()
            nitArrayLine[nanL] = 0
            nitrogenArray = nitArrayLine.reshape(grid_shape)
    return nitrogenArray

def run_simulation(out_file_name: str, model_json:dict, input_data: AllData):

    t_init = time.time()

    dataBounds = {
        'northward_Water_current': [-1e7,1e7],
        'eastward_Water_current': [-1e7,1e7],
        'Nitrate': [-1e-2,1e4],
        'Ammonium': [-1e-2,1e4],
        'Temperature': [-1e4, 1e4],
        'Phosphate' : [-1e-2,1e4],
        'par' : [-1e-2,1e4]
    }

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
    longitudes, _ = input_data.parameterData['Nitrate'].getVariable('longitude', **data_kwargs)
    latitudes, _ = input_data.parameterData['Nitrate'].getVariable('latitude', **data_kwargs)
    dxMeter, dyMeter = giveDxDy(latitudes, longitudes)
    latRef = np.zeros((len(latitudes), len(longitudes)))
    latRef[:, :] = latitudes[np.newaxis].T

    mask = working_data['northward_Water_current'].mask
    for par_name, par_data in input_data.parameterData.items():
        working_data[par_name] = np.ma.masked_outside(working_data[par_name], dataBounds[par_name][0],
                                                      dataBounds[par_name][1])
        if par_name != 'par':
            working_data[par_name].filled(fill_value=0)
        else:
            working_data[par_name].filled(fill_value=8)
        #mask = np.ma.mask_or(mask, working_data[par_name].mask)

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
        mask = resa.resampleData(mask)
    nanLists = findNan(mask)

    for par_name, par_data in input_data.parameterData.items():
        print(par_name)
        print(working_data[par_name].shape)


    # Initialize the model variables
    state_vars = {
        'cNO3': np.ma.masked_array(np.zeros(grid_shape), mask),
        'cNH4': np.ma.masked_array(np.zeros(grid_shape), mask),
        'N_s': np.ma.masked_array(np.zeros(grid_shape), mask),
        'N_f': np.ma.masked_array(np.ones(grid_shape) * parms_harvest['deployment_Nf'], mask),
        'D': np.ma.masked_array(np.ones(grid_shape) * parms_run["Detritus"], mask)
    }

    availableNut = {
        'avNO3': np.ma.masked_array(np.zeros(grid_shape), mask),
        'avNH4': np.ma.masked_array(np.zeros(grid_shape), mask),
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
        for par_name, par_data in input_data.parameterData.items():
            new_i = iNearest(data_date, time_axes[par_name])
            if new_i != nearest_time_i[par_name]:
                nearest_time_i[par_name] = new_i
                working_data[par_name], _ = par_data.getVariable(time_index=new_i, **data_kwargs)

                if scenC:
                    working_data[par_name] = resa.resampleData(working_data[par_name])
                # Update the centered currents as well
                if par_name == "eastward_Water_current":
                    working_data["decentered_U"][:, 1:-1] = (working_data['eastward_Water_current'][:, 1:] + working_data['eastward_Water_current'][:, :-1]) / 2
                if par_name == "northward_Water_current":
                    working_data["decentered_V"][1:-1, :] = (working_data['northward_Water_current'][1:, :] + working_data['northward_Water_current'][:-1, :]) / 2
                if par_name != 'par':
                    working_data[par_name].filled(fill_value=0)
                else:
                    working_data[par_name].filled(fill_value=8)

        availableNut_term = give_availableNut(working_data=working_data,
                                              dt=dt, dxMeter=dxMeter, dyMeter=dyMeter, nanLists=nanLists)
        for var_name in ['avNO3', 'avNH4']:
            availableNut[var_name] += availableNut_term[var_name]

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
                      variableNames=['NH4', 'NO3', 'N_s', 'N_f', 'D', 'avNH4', 'avNO3', 'cNO3', 'cNH4'], mask=mask)

    # Write values to file
    ds = nc.Dataset(out_file_name, 'a')
    for name in model.names:
        ds[name][0,:,:] = np.ma.masked_array(state_vars[name], mask)
    for name in ['cNO3', 'cNH4']:
        ds[name][0,:,:] = np.ma.masked_array(state_vars[name], mask)
    for name in ['avNH4', 'avNO3']:
        ds[name][0, :, :] = np.ma.masked_array(availableNut[name], mask)
        availableNut[name] = np.ma.masked_array(availableNut[name], mask)
        availableNut[name][availableNut[name].mask] = np.nan
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

def give_availableNut(working_data: dict, dt, dxMeter: np.array, dyMeter: np.array,nanLists: np.array):
    mask = working_data['Nitrate'].mask
    grid_shape = working_data['Nitrate'].shape

    NO3_line = working_data['Nitrate'].flatten()
    NH4_line = working_data['Ammonium'].flatten()
    dx = dxMeter.flatten()
    dy = dyMeter.flatten()

    matE = createMatE_available(working_data["decentered_U"], nanLists)
    matN = createMatN_available(working_data["decentered_V"], nanLists)

    # we compute the advection terms
    advNO3 = -(dt / dx) * matE.dot(NO3_line)-(dt / dy) * matN.dot(NO3_line)
    advNH4 = -(dt / dx) * matE.dot(NH4_line)-(dt / dy) * matN.dot(NH4_line)

    # reshape to the grid
    advNO3 = advNO3.reshape(grid_shape)
    advNH4 = advNH4.reshape(grid_shape)

    # reapply masks that were altered
    advNO3 = np.ma.masked_array(advNO3, mask)
    advNH4 = np.ma.masked_array(advNH4, mask)

    all_terms = {
        'avNO3': advNO3,
        'avNH4': advNH4,
    }

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

    for name in variableNames:
        var = ds.createVariable(name, 'f4', ('time', 'latitude', 'longitude',))
        var[:,:,:] = np.ma.masked_array(np.nan*np.ones(full_mask.shape), full_mask)

    ds.close()


def open_data_input(file_adress:str, zone:str, paramNames:list, dataRef: pd.DataFrame, frequency='monthly', dataType='model'):

    dataRows = [dataRef.index[(dataRef['Parameter']==param) & (dataRef['Place']==zone) &
                              (dataRef['frequency']==frequency) & (dataRef['type']==dataType)][0]
                    for param in paramNames]

    # Gives the argument to ParamData() corresponding to a column in dataRef
    columns_arguments = {
        'variable': 'variable_name',
        'latName': 'latitude_name',
        'longName': 'longitude_name',
        'timeName': 'time_name',
        'depthName': 'depth_name',
        'unitFactor': 'unit_conversion'
    }

    # Gives the default value corresponding to a column in dataRef
    fill_na = {
        'variable': 'N/A',
        'latName': 'latitude',
        'longName': 'longitude',
        'timeName': 'time',
        'depthName': 'depth',
        'unitFactor': 1
    }

    parameter_dict = {parName: {} for parName in paramNames}
    for parName, iRow in zip(paramNames, dataRows):

        parameter_dict[parName]['file_name'] = file_adress.format(zone=zone, param=parName)

        # Fill all argNames except file_name, time_zero, and time_step
        for colName, argName in columns_arguments.items():
            parameter_dict[parName][argName] = dataRef[colName].fillna(fill_na[colName])[iRow]

        # timeOrigin must be expressed in ISO 8601"
        parameter_dict[parName]['time_zero'] = datetime.datetime.strptime(dataRef['timeOrigin'][iRow], "%Y-%m-%dT%H:%M:%SZ")

        # convert dataRef['timeUnit'][iRow] to a float if possible
        try:
            timeUnit = float(dataRef['timeUnit'][iRow])
        except:
            timeUnit = dataRef['timeUnit'][iRow]

        if type(timeUnit) is str:
            parameter_dict[parName]['time_step'] = datetime.timedelta(**{timeUnit: 1})
        else: # should be a number
            parameter_dict[parName]['time_step'] = datetime.timedelta(seconds=timeUnit)

    return parameter_dict


def dataCmd_to_AllData(dataCmdDict: dict, adress_format:str):
    """
    Converts a dictionary containing items from lines in dataCmd.csv to a
    dictionary of inputs that can be passed to the AllData() constructor.
    """

    # Gives the argument to ParamData() corresponding to a column in dataRef
    columns_arguments = {
        'variable': 'variable_name',
        'latName': 'latitude_name',
        'longName': 'longitude_name',
        'timeName': 'time_name',
        'depthName': 'depth_name',
        'unitFactor': 'unit_conversion'
    }

    parameter_dict = {parName: {} for parName in dataCmdDict.keys()}
    for parName, line in dataCmdDict.items():

        # adress_format can use any value from a dataCmd column (Place, Parameter, or more)
        parameter_dict[parName]['file_name'] = adress_format.format(**line)

        # Fill all argNames except file_name, time_zero, and time_step
        for colName, argName in columns_arguments.items():
            if line[colName] != "":
                parameter_dict[parName][argName] = line[colName]

        # timeOrigin must be expressed in ISO 8601"
        parameter_dict[parName]['time_zero'] = datetime.datetime.strptime(line['timeOrigin'], "%Y-%m-%dT%H:%M:%SZ")

        # convert line['timeUnit'] to a float if possible
        try:
            timeUnit = float(line['timeUnit'])
        except:
            timeUnit = line['timeUnit']

        if type(timeUnit) is str:
            parameter_dict[parName]['time_step'] = datetime.timedelta(**{timeUnit: 1})
        else: # should be a number
            parameter_dict[parName]['time_step'] = datetime.timedelta(seconds=timeUnit)

    return parameter_dict


if __name__=="__main__":

    input_args = {
        'zone' : "IBI",
        'file_adress' : '/media/share/data/{zone}/{param}/{param}{zone}modelNetCDF2021-01to2022-01.nc',
        'dataRef' : pd.read_csv('/media/global/dataCmd.csv', delimiter=';'),
        'paramNames' : ['Ammonium', 'Nitrate', 'Temperature', 'northward_Water_current', 'eastward_Water_current']
    }
    ### Initialize the netcdf reading interface
    algaeData = open_data_input(**input_args)


    ### get the copernicus grid and mask

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

    ###


    model_params = "macroalgae_model_parameters_input.json"
    json_data = import_json(model_params)
    model = MA_model_scipy(json_data['parameters'])

    n_slices = 10

    lon_split = np.array_split(longitudes, n_slices)
    mask_split = np.array_split(mask, n_slices, 1)


    ### Create datasets for monthly sim
    for i, (lon_arr, mask_arr) in enumerate(zip(lon_split, mask_split)):
        initialize_result(f"/media/share/results/simulations/monthly/monthly_simulations_{i:03d}.nc",
        np.array(range(1,13)), latitudes, lon_arr, model.names, mask_arr)

    pool = mp.Pool(10)

    y0 = np.array([0, 0, 0, 1000, 0], dtype=np.float64)


    t0 = time.time()
    n_cells = pool.starmap_async(run_scenario_a_monthly,
        [(f"/media/share/results/simulations/monthly/monthly_simulations_{i:03d}.nc", 
            json_data['parameters'], y0, input_args, 2021, True, True) for i in range(n_slices)]).get()
    #n_cells = run_scenario_a_monthly("/media/share/results/complete_simulations_monthly_test_0.nc", model, y0, input_args, 2021)


    print(n_cells)
    print((time.time()-t0)/sum(n_cells))
