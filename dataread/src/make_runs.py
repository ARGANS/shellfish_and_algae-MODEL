import datetime
import numpy as np
import pandas as pd
import time
import scipy.sparse as sp
from skimage.morphology import label
import datetime
from dateutil.relativedelta import *
from .read_netcdf import *
from .launch_model import *
from .utils import import_json
import os

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


class Decenterer:
    def __init__(self, mask, dx, dy):

        self._mask = mask.copy()

        m, n = self._mask.shape

        # gradl(L(U, V, l)) = A * X + B

        # Blocks of A
        A_uu = sp.eye(m*(n-1)) * 2
        A_vv = sp.eye((m-1)*n) * 2
        A_uv = None
        A_vu = None

        diags_ul = np.ones((2, m*n))
        diags_ul[1, :] = -1
        ul_block = sp.dia_matrix((diags_ul, [0, 1]), shape=(n-1, n)) * dy
        A_ul = sp.block_diag([ul_block]*m)

        A_vl = sp.lil_matrix(((m-1)*n, m*n))
        for k in range((m-1)*n):
            A_vl[k, k] = dx
            A_vl[k, k + n] = -dx

        A_lu = A_ul.T
        A_lv = A_vl.T
        A_ll = None

        A = sp.bmat([[A_uu, A_uv, A_ul],
                     [A_vu, A_vv, A_vl],
                     [A_lu, A_lv, A_ll]])


        # masks for where currents should be 0, obtained from diff.
        fake_UVc = np.ma.masked_array(np.zeros((m,n)), self._mask)
        Uc_mask = np.ma.getmaskarray(np.diff(fake_UVc, axis=1).flatten())
        Vc_mask = np.ma.getmaskarray(np.diff(fake_UVc, axis=0).flatten())

        # Find all independent areas (SLOW)
        labs, n_labs = label(fake_UVc.filled(1), background=1, return_num=True, connectivity=1)

        # mask where the lambda constraints are unnecessary/redundant i.e.:
        #   - where Uc/Vc is masked
        #   - one cell per independent area
        l_mask = self._mask.flatten()
        for i in range(1, n_labs+1):
            first_lab_index = np.where(labs.flatten() == i)[0][0]
            l_mask[first_lab_index] = True

        # Define the rows/columns of the matrix that will form a non singular matrix
        self._where_not_zero = np.where(np.logical_not(np.concatenate((Uc_mask, Vc_mask, l_mask))))[0]

        # Remove rows and columns where the values are set to 0
        A = A.tocsc()[:, self._where_not_zero]
        A = A[self._where_not_zero, :]

        self._A = A

        # invert the matrix
        #self._Ainv = sp.linalg.inv(A)

    def apply(self, Uc, Vc):

        m, n = self._mask.shape

        # definition of B
        B_u = - (Uc[:, 1:].filled(0) + Uc[:, :-1].filled(0)).flatten()
        B_v = - (Vc[1:, :].filled(0) + Vc[:-1, :].filled(0)).flatten()
        B_l = np.zeros((m, n)).flatten()

        B = np.concatenate((B_u, B_v, B_l))[np.newaxis].T

        # Initialize result array with zeros and mask
        X_full = np.ma.masked_array(np.zeros(B.shape), np.ones(B.shape))

        B = B[self._where_not_zero]

        #X = - self._Ainv.dot(B)
        X = sp.linalg.spsolve(self._A, -B)

        # Apply the result where relevant in X_full
        X_full[self._where_not_zero] = X[np.newaxis].T


        U = np.ma.masked_array(np.zeros((m, n+1)))
        V = np.ma.masked_array(np.zeros((m+1, n)))
        U[:, 1:-1] = X_full[0:(m*(n-1))].reshape((m, n-1))
        V[1:-1, :] = X_full[(m*(n-1)):(m*(n-1) + (m-1)*n)].reshape((m-1, n))

        return U.filled(0), V.filled(0)


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

#return dt number in function of the current speed and time steps
def give_dt(dx, dy, Ewc, Nwc):
    Cx = np.abs(Ewc[:, 1:] / dx).flatten()
    Cy = np.abs(Nwc[1:, :] / dy).flatten()
    return 0.9/np.max([Cx, Cy])

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
            shiftMatE = sp.spdiags([1]*nbrx, -jRel, nbrx, nbrx).todense()
            # Matrix to shift mask along N-S by -iRel
            shiftMatN = sp.spdiags([1]*nbry, iRel, nbry, nbry).todense()
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
    Mat = sp.dia_matrix((data, offset), shape=(nbrx * nbry, nbrx * nbry))

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
    Mat = sp.dia_matrix((data, offset), shape=(nbrx * nbry, nbrx * nbry))

    return Mat


def createMatEupwind(decenteredEwc):

    decenteredEwcPlus = decenteredEwc[:, 1:]
    decenteredEwcMinus = decenteredEwc[:, :-1]

    nbrx, nbry = decenteredEwcPlus.shape

    uPlusGreater0 = ((decenteredEwcPlus > 0) * 1).flatten()
    uMinusGreater0 = ((decenteredEwcMinus > 0) * 1).flatten()
    termAPlus = uPlusGreater0
    termBPlus = 1 - uPlusGreater0
    termAMinus = 1 - uMinusGreater0
    termBMinus = uMinusGreater0

    # Plus half matrix
    offset = np.array([0, 1])
    diags = np.zeros((2, nbrx * nbry))
    diags[0, :] = termAPlus
    diags[1, 1:] = termBPlus[:-1]
    matEplus = sp.dia_matrix((diags, offset), shape=(nbrx * nbry, nbrx * nbry))

    # Minus half matrix
    offset = np.array([0, -1])
    diags = np.zeros((2, nbrx * nbry))
    diags[0, :] = termAMinus
    diags[1, :-1] = termBMinus[1:]
    matEminus = sp.dia_matrix((diags, offset), shape=(nbrx * nbry, nbrx * nbry))

    return matEplus, matEminus


def createMatNupwind(decenteredNwc):

    decenteredNwcPlus = decenteredNwc[1:, :]
    decenteredNwcMinus = decenteredNwc[:-1, :]

    nbrx, nbry = decenteredNwcPlus.shape

    uPlusGreater0 = ((decenteredNwcPlus > 0) * 1).flatten()
    uMinusGreater0 = ((decenteredNwcMinus > 0) * 1).flatten()
    termAPlus = uPlusGreater0
    termBPlus = 1 - uPlusGreater0
    termAMinus = 1 - uMinusGreater0
    termBMinus = uMinusGreater0

    # Plus half matrix
    offset = np.array([0, nbry])
    diags = np.zeros((2, nbrx * nbry))
    diags[0, :] = termAPlus
    diags[1, nbry:] = termBPlus[:-nbry]
    matNplus = sp.dia_matrix((diags, offset), shape=(nbrx * nbry, nbrx * nbry))

    # Minus half matrix
    offset = np.array([0, -nbry])
    diags = np.zeros((2, nbrx * nbry))
    diags[0, :] = termAMinus
    diags[1, :-nbry] = termBMinus[nbry:]
    matNminus = sp.dia_matrix((diags, offset), shape=(nbrx * nbry, nbrx * nbry))

    return matNplus, matNminus


#returns the space step
def giveDxDy(latitudes, longitudes):

    lon_step_average = (longitudes[-1] - longitudes[0]) / (len(longitudes) - 1)
    lat_step_average = (latitudes[-1] - latitudes[0]) / (len(latitudes) - 1)

    latRef_average = (latitudes[-1] + latitudes[0]) / 2

    dxMeter, dyMeter = degrees_to_meters(lon_step_average, lat_step_average, latRef_average)

    return dxMeter, dyMeter


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

def latLon_to_xy(lat,lon, longitudes, latitudes, nanLists):
    lonNear, latNear = iNearest(lon, longitudes), iNearest(lat, latitudes)
    mat_mask = np.zeros(len(latitudes)*len(longitudes))
    mat_mask[nanLists[0,0]] = 1
    mat_mask = mat_mask.reshape(len(latitudes),len(longitudes))
    if mat_mask[latNear,lonNear]==0:
        return latNear,lonNear
    else:
        print('---------------------------------------')
        latRef_average = (latitudes[-1] + latitudes[0]) / 2
        lonListMeter, latListMeter = degrees_to_meters(np.abs(longitudes[lonNear-1:lonNear+2]-longitudes[lonNear]), np.abs(latitudes[latNear-1:latNear+2]-latitudes[latNear]), latRef_average)
        listPosition = []
        for i in [-1,0,1]:
            for j in [-1,0,1]:
                if (mat_mask[latNear+i,lonNear+j] == 0):
                    print(i,j)
                    dist = np.sqrt(lonListMeter[j]**2+latListMeter[i]**2)
                    listPosition.append((latNear+i,lonNear+j,dist))
        finalArray = np.array(listPosition,dtype=[('poslat',float),('poslon',float),('dist',float)])
        if len(finalArray)>0:
            bestCoor = np.sort(finalArray,order='dist')[0]
        else:
            bestCoor = [latNear,lonNear]
        print(bestCoor)
        return int(bestCoor[0]),int(bestCoor[1])
        

def giveFarmPos(farmList, longitudes, latitudes, nanLists):
    xList = []
    yList = []
    for i in range(len(farmList)):
        if (farmList[i][0]<latitudes[-1]) and (farmList[i][0]>latitudes[0]) and (farmList[i][1]<longitudes[-1]) and (farmList[i][1]>longitudes[0]):
            lat,lon = farmList[i][0], farmList[i][1]
            yi, xi= latLon_to_xy(lat,lon, longitudes, latitudes, nanLists)
            xList.append(xi)
            yList.append(yi)
    return (np.array(yList),np.array(xList))

def run_simulation(out_file_name: str, model_json:dict, input_data: AllData, farm_pos_file=None):

    t_init = time.time()

    dataBounds = {
        'northward_Water_current': [-1e6, 1e6],
        'eastward_Water_current': [-1e6, 1e6],
        'Nitrate': [-1e-2, 1e4],
        'Ammonium': [-1e-2, 1e4],
        'Temperature': [-1e4, 1e4],
        'Phosphate': [-1e-2, 1e4],
        'par': [-1e-2, 1e4]
    }

    #dt = 1 / 72  # days # TODO: make into parameter in json

    # Parse the input json info
    parms_run = list(model_json['parameters']['run'].values())[0]['parameters']
    parms_farm = list(model_json['parameters']['farm'].values())[0]['parameters']
    parms_harvest = list(model_json['parameters']['harvest'].values())[0]['parameters']
    harvest_type = list(model_json['parameters']['harvest'].keys())[0]

    scenC = (model_json['metadata']['scenario']=="C" and (not farm_pos_file))

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
        working_data[par_name], data_dims = par_data.getVariable(time_index=nearest_time_i[par_name], **data_kwargs) #TODO: are dims really laways lat,lon ?
        print(par_name)
        print(working_data[par_name].shape)
        print(data_dims)

    init_grid_shape = working_data['Nitrate'].shape #the shape should be the same for all parameters
    longitudes, _ = input_data.parameterData['Nitrate'].getVariable('longitude', **data_kwargs)
    latitudes, _ = input_data.parameterData['Nitrate'].getVariable('latitude', **data_kwargs)
    dxMeter, dyMeter = giveDxDy(latitudes, longitudes)

    latRef = np.zeros((len(latitudes), len(longitudes)))
    latRef[:, :] = latitudes[np.newaxis].T

    for par_name, par_data in input_data.parameterData.items():
        working_data[par_name] = np.ma.masked_outside(working_data[par_name], dataBounds[par_name][0],dataBounds[par_name][1])
        if (par_name == "Nitrate") or (par_name == "Ammonium"):
            working_data[par_name] = np.maximum(working_data[par_name], 0)

    # Iniitializing the mask to be used, based on the first time step.
    for par_name, par_data in input_data.parameterData.items():
        working_data[par_name] = np.ma.masked_outside(working_data[par_name], dataBounds[par_name][0],dataBounds[par_name][1])
        if (par_name == "Nitrate") or (par_name == "Ammonium"):
            working_data[par_name] = np.maximum(working_data[par_name], 0)

    mask = working_data['Temperature'].mask.copy()
    for par_name in ['Nitrate', 'Ammonium', 'Phosphate', 'northward_Water_current', 'eastward_Water_current']:
        mask = np.ma.mask_or(mask, working_data[par_name].mask)

    # In the case of scenario C, resample the data to near 1 mile.
    grid_shape = init_grid_shape
    if scenC:
        dxRatio = 1852 / np.mean(dxMeter) #1852 meters = 1 nautical mile
        dyRatio = 1852 / np.mean(dyMeter)
        resa = Resampler(dxRatio,dyRatio,init_grid_shape)
        grid_shape = resa.grid_shape
        for par_name, par_data in input_data.parameterData.items():
            working_data[par_name] = resa.resampleData(working_data[par_name])
        latRef = resa.resampleData(latRef)
        '''dyMeter = resa.resampleData(dyMeter)
        dxMeter = resa.resampleData(dxMeter)'''
        mask = resa.resampleData(mask)
    nanLists = findNan(mask)

    #if we only put farms in the location given in farm_pos_file
    if farm_pos_file:
        farmList = np.loadtxt(farm_pos_file, dtype=float, delimiter=';', skiprows=1, usecols=(1, 2))
        mask_farm = giveFarmPos(farmList, longitudes, latitudes, nanLists)
    else:
        mask_farm = None

    for par_name, par_data in input_data.parameterData.items():
        print(par_name)
        print(working_data[par_name].shape)

    # Initialize the model variables
    if mask_farm:
        state_vars = {
            'cNO3': np.ma.masked_array(np.zeros(grid_shape), mask.copy()),
            'cNH4': np.ma.masked_array(np.zeros(grid_shape), mask.copy()),
            'N_s': np.ma.masked_array(np.zeros(grid_shape), mask.copy()),
            'N_f': np.ma.masked_array(np.zeros(grid_shape), mask.copy()),
            'D': np.ma.masked_array(np.zeros(grid_shape), mask.copy())
        }
        state_vars['N_f'][mask_farm] = parms_harvest['deployment_Nf']
    else:
        state_vars = {
            'cNO3': np.ma.masked_array(np.zeros(grid_shape), mask.copy()),
            'cNH4': np.ma.masked_array(np.zeros(grid_shape), mask.copy()),
            'N_s': np.ma.masked_array(np.zeros(grid_shape), mask.copy()),
            'N_f': np.ma.masked_array(np.ones(grid_shape)*parms_harvest['deployment_Nf'], mask.copy()),
            'D': np.ma.masked_array(np.zeros(grid_shape), mask.copy())
        }

    availableNut = {
        'avNO3': np.ma.masked_array(np.zeros(grid_shape), mask.copy()),
        'avNH4': np.ma.masked_array(np.zeros(grid_shape), mask.copy()),
    }

    if scenC:
        state_vars['N_s'] = prepareScenC(state_vars['N_s'], nanLists, grid_shape)
        state_vars['N_f'] = prepareScenC(state_vars['N_f'], nanLists, grid_shape)

    # Intitialize the decenterer and apply it to the first U/V data
    t_init_decenterer = time.time()
    print('Starting initialization of the decenterer')
    decenterer = Decenterer(mask, dxMeter, dyMeter)
    print(f'End of decenterer initialization, time taken: {(time.time() - t_init_decenterer)} seconds')

    t_init_decenterer = time.time()
    print('Starting applying the decenterer ofr the first time')

    working_data['decentered_U'], working_data['decentered_V'] = decenterer.apply(working_data['eastward_Water_current'],
                                                                                  working_data['northward_Water_current'])

    dt_phys = give_dt(dxMeter, dyMeter, working_data["decentered_U"],working_data["decentered_V"])

    print(f'End of first decenterer application, time taken: {(time.time() - t_init_decenterer)} seconds')

    # Simulation loop
    date_month = 0
    deficitNbr = 0
    sim_date = startDate
    while sim_date < endDate:
        print(f'{sim_date}')

        # Alter the date after new year in case of winter growth
        if harvest_type == "Winter_growth":
            data_date = sim_date.replace(year = year)
        else:
            data_date = sim_date

        # For each dataset, if the nearest i has changed, update the working data
        reapply_decenterer = False
        for par_name, par_data in input_data.parameterData.items():
            new_i = iNearest(data_date, time_axes[par_name])
            if new_i != nearest_time_i[par_name]:
                nearest_time_i[par_name] = new_i
                working_data[par_name], _ = par_data.getVariable(time_index=new_i, **data_kwargs)
                working_data[par_name] = np.ma.masked_outside(working_data[par_name], dataBounds[par_name][0],
                                                              dataBounds[par_name][1])
                if scenC:
                    working_data[par_name] = resa.resampleData(working_data[par_name])
                # Update the centered currents as well
                if par_name == "eastward_Water_current" or par_name == "northward_Water_current":
                    reapply_decenterer = True
                elif (par_name == "Nitrate") or (par_name == "Ammonium"):
                    working_data[par_name] = np.maximum(working_data[par_name], 0)
        if reapply_decenterer:
            working_data['decentered_U'], working_data['decentered_V'] = decenterer.apply(working_data['eastward_Water_current'],
                                                                                          working_data['northward_Water_current'])
            dt_phys = give_dt(dxMeter, dyMeter, working_data["decentered_U"],
                                  working_data["decentered_V"])


        # Compute the BGC terms
        days = (data_date - datetime.datetime(year, 1, 1)).days # Note: returns an integer only, that is sufficient precision for this

        bgc_terms = bgc_model(state_vars=state_vars, working_data=working_data,
                              model=model, parms_run=parms_run, days=days, latRef=latRef)

        # Value of dt for Ns and Nf > 0
        dt_list = []
        for var_name in ['N_s', 'N_f']:
            dt_array = - state_vars[var_name]/bgc_terms[var_name]
            dt_array[dt_array <= 0] = np.inf
            dt_positive = np.amin(dt_array)

            dt_list.append(dt_positive * 0.9)
        dt_bgc = min(dt_list)

        # Update dt to respect the CFL and the BGC terms staying positive
        print(f'dt for phys: {dt_phys * 24 * 60} min ; dt for bgc: {dt_bgc * 24 * 60} min')
        dt = min(dt_phys, dt_bgc)

        # Apply the bgc terms
        if scenC:
            for var_name in state_vars.keys():
                bgc_terms[var_name] = prepareScenC(bgc_terms[var_name], nanLists, grid_shape)
        for var_name in state_vars.keys():
            state_vars[var_name][mask_farm] += bgc_terms[var_name][mask_farm] * dt

        # "Kill" functionally dead plants
        dead_cells = (state_vars['N_f'] < 1e-2)
        state_vars['N_f'][dead_cells] = 0
        state_vars['N_s'][dead_cells] = 0

        # Compute the maximum available nutrients
        availableNut_term = give_availableNut(working_data=working_data,
                                              dt=dt, dxMeter=dxMeter, dyMeter=dyMeter, nanLists=nanLists)
        for var_name in ['avNO3', 'avNH4']:
            availableNut[var_name] += availableNut_term[var_name]

        # Compute the advection terms
        if ((model_json['metadata']['scenario'] == "A") and (not farm_pos_file)):
            advection_terms = advection_modelA(state_vars=state_vars, working_data=working_data,
                                               dxMeter=dxMeter, dyMeter=dyMeter)
        else:
            advection_terms = advection_model(state_vars=state_vars, working_data=working_data,
                                              dxMeter=dxMeter, dyMeter=dyMeter)

        # Apply the advection
        for var_name in state_vars.keys():
            state_vars[var_name] += advection_terms[var_name] * dt


        # Time dissipation of the signal
        dissip_t = 10 #days
        state_vars['cNH4'] = state_vars['cNH4'] * (1 - dt/dissip_t)
        state_vars['cNO3'] = state_vars['cNO3'] * (1 - dt/dissip_t)
        state_vars['D'] = state_vars['D'] * (1 - dt/dissip_t)

        if date_month != sim_date.month:
            date_month = sim_date.month
            unitsDict = {
                'cNO3': 'mg N/m^3',
                'cNH4': 'mg N/m^3',
                'CMEMS_NO3': 'mg N/m^3',
                'CMEMS_NH4': 'mg N/m^3'}
            tempFileName = out_file_name[:-3]+f'{deficitNbr:02d}'+'.nc'
            print(tempFileName)
            # Create output file
            initialize_result(tempFileName, times=[0], latitudes=latitudes, longitudes=longitudes,
                            variableNames=['cNO3', 'cNH4', 'CMEMS_NO3', 'CMEMS_NH4'], unitsDict=unitsDict, mask=mask)
            # Write values to file
            ds = nc.Dataset(tempFileName, 'a')
            for name in ['cNO3', 'cNH4']:
                ds[name][0,:,:] = np.ma.masked_array(state_vars[name], mask.copy())
            ds['CMEMS_NO3'][0,:,:] = np.ma.masked_array(working_data['Nitrate'], mask.copy())
            ds['CMEMS_NH4'][0,:,:] = np.ma.masked_array(working_data['Ammonium'], mask.copy())
            ds.close()
            deficitNbr += 1

        sim_date += datetime.timedelta(days = dt)

    if scenC:
        latStep = (latitudes[-1] - latitudes[0]) / (grid_shape[0] - 1)
        lonStep = (longitudes[-1] - longitudes[0]) / (grid_shape[1] - 1)

        latitudes = latStep * np.arange(grid_shape[0]) + latitudes[0]
        longitudes = lonStep * np.arange(grid_shape[1]) + longitudes[0]

    state_vars['CMEMS_NH4'] = working_data['Ammonium']
    state_vars['CMEMS_NO3'] = working_data['Nitrate']

    unitsDict = {'CMEMS_NO3': 'mg/m^3',
                'CMEMS_NH4': 'mg/m^3',
                'cNO3': 'mg N/m^3',
                'cNH4': 'mg N/m^3',
                'D': 'mg N/m^3',
                'N_f': 'mg N/m^3',
                'N_s': 'mg N/m^3',
                'avNO3': 'mg N/m^3',
                'avNH4': 'mg N/m^3'}
    # Create output file
    initialize_result(out_file_name, times=[0], latitudes=latitudes, longitudes=longitudes,
                      variableNames=['CMEMS_NH4', 'CMEMS_NO3', 'N_s', 'N_f', 'D', 'avNH4', 'avNO3', 'cNO3', 'cNH4'], unitsDict=unitsDict, mask=mask)
                      
    # Write values to file
    ds = nc.Dataset(out_file_name, 'a')
    for name in model.names:
        ds[name][0,:,:] = np.ma.masked_array(state_vars[name], mask.copy())
    for name in ['cNO3', 'cNH4']:
        ds[name][0,:,:] = np.ma.masked_array(state_vars[name], mask.copy())
    for name in ['avNH4', 'avNO3']:
        ds[name][0, :, :] = np.ma.masked_array(availableNut[name], mask.copy())
        availableNut[name] = np.ma.masked_array(availableNut[name], mask.copy())
        availableNut[name][availableNut[name].mask] = np.nan
    ds.close()
    print(time.time() - t_init)
    return time.time() - t_init


#return the variation of each quantities in the algae model
def bgc_model(state_vars: dict, working_data: dict, model, parms_run, days, latRef):

    data_in = {
        'SST': working_data['Temperature'],
        'PAR': working_data['par'],
        'PO4_ext': working_data['Phosphate'],
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

    terms_list = model.hadley_advection(days, y, data_in, latRef, state_vars['cNH4'])
    all_terms = dict(zip(["cNH4", "cNO3", "N_s", "N_f", "D"], terms_list))

    return all_terms


def give_availableNut(working_data: dict, dt, dxMeter: np.array, dyMeter: np.array, nanLists: np.array):
    mask = working_data['Nitrate'].mask.copy()
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
    advNO3 = np.ma.masked_array(advNO3, mask.copy())
    advNH4 = np.ma.masked_array(advNH4, mask.copy())

    all_terms = {
        'avNO3': advNO3,
        'avNH4': advNH4,
    }

    return all_terms

def advection_modelA(state_vars: dict, working_data: dict, dxMeter: np.array, dyMeter: np.array):
    mask = working_data['Nitrate'].mask.copy()
    grid_shape = working_data['Nitrate'].shape

    uGreater0 = ((working_data["decentered_U"][:, 1:] > 0) * working_data["decentered_U"][:, 1:]).flatten()
    uLower0 = ((working_data["decentered_U"][:, :-1] < 0) * working_data["decentered_U"][:, :-1]).flatten()

    vGreater0 = ((working_data["decentered_V"][1:] > 0) * working_data["decentered_V"][1:]).flatten()
    vLower0 = ((working_data["decentered_V"][:-1] < 0) * working_data["decentered_V"][:-1]).flatten()

    cNO3_line = state_vars['cNO3'].flatten()
    cNH4_line = state_vars['cNH4'].flatten()
    D_line_eps = state_vars['D'].flatten()
    dx = dxMeter.flatten()
    dy = dyMeter.flatten()

    #we compute the advection terms
    advNO3 = ((1 / dx) * (-uGreater0 + uLower0) + (1 / dy) * (-vGreater0 + vLower0))*cNO3_line
    advNH4 = ((1 / dx) * (-uGreater0 + uLower0) + (1 / dy) * (-vGreater0 + vLower0))*cNH4_line
    advD = ((1 / dx) * (-uGreater0 + uLower0) + (1 / dy) * (-vGreater0 + vLower0))*D_line_eps

    # reshape to the grid
    advNO3 = advNO3.reshape(grid_shape)
    advNH4 = advNH4.reshape(grid_shape)
    advD = advD.reshape(grid_shape)

    # reapply masks that were altered # TODO: check if still necessary
    advNO3.mask = mask.copy()
    advNH4.mask = mask.copy()
    advD.mask = mask.copy()

    all_terms = {
        'cNO3': advNO3,
        'cNH4': advNH4,
        'N_s': 0,
        'N_f': 0,
        'D': advD
    }

    return all_terms

def advection_model(state_vars: dict, working_data: dict, dxMeter: float, dyMeter: float):

    grid_shape = working_data['Nitrate'].shape

    matE_plus_half, matE_less_half = createMatEupwind(working_data["decentered_U"])
    matN_plus_half, matN_less_half = createMatNupwind(working_data["decentered_V"])

    dx = dxMeter.flatten()
    dy = dyMeter.flatten()

    all_terms = {
        'N_s': 0,
        'N_f': 0
    }

    for var_name in ['cNO3', 'cNH4', 'D']:

        var_line = state_vars[var_name].flatten()

        # Obtain the left or right (up or down) value of var depending on the sign of currents
        var_hat_plus_half_u = matE_plus_half.dot(var_line)
        var_hat_plus_half_v = matN_plus_half.dot(var_line)

        var_hat_less_half_u = matE_less_half.dot(var_line)
        var_hat_less_half_v = matN_less_half.dot(var_line)

        # Compute the left and right (up and down) flow rates for each cell
        Fmat_plus_half = var_hat_plus_half_u * working_data["decentered_U"][:, 1:].flatten()
        Gmat_plus_half = var_hat_plus_half_v * working_data["decentered_V"][1:, :].flatten()

        Fmat_less_half = var_hat_less_half_u * working_data["decentered_U"][:, :-1].flatten()
        Gmat_less_half = var_hat_less_half_v * working_data["decentered_V"][:-1, :].flatten()

        # compute the advection terms
        advected = - (1 / dy) * (Gmat_plus_half - Gmat_less_half) - (1 / dx) * (Fmat_plus_half - Fmat_less_half)

        # reshape to the grid and add to output
        all_terms[var_name] = advected.reshape(grid_shape)

    return all_terms


def initialize_result(fileName:str, times, latitudes, longitudes,
                      variableNames:list, unitsDict:dict, mask:np.array):
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
        var.units = unitsDict[name]
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
    pass
