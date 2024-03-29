import datetime
import numpy as np
import pandas as pd
import time
import scipy.sparse as sp
import datetime
from dateutil.relativedelta import *
from .read_netcdf import *
from .launch_model import *
from .utils import import_json


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


#returns the space step
def giveDxDy(latitudes, longitudes):

    lon_step_average = (longitudes[-1] - longitudes[0]) / (len(longitudes) - 1)
    lat_step_average = (latitudes[-1] - latitudes[0]) / (len(latitudes) - 1)

    latRef_average = (latitudes[-1] + latitudes[0]) / 2

    dxMeter, dyMeter = degrees_to_meters(lon_step_average, lat_step_average, latRef_average)

    return dxMeter, dyMeter

def run_simulation(out_file_name: str, model_json:dict, input_data: AllData):

    t_init = time.time()

    dataBounds = {
        'northward_Water_current': [-1e7,1e7],
        'eastward_Water_current': [-1e7,1e7],
        'Chlorophyll-a': [-1e-2,1e4],
        'Temperature': [-1e4, 1e4],
    }

    # Parse the input json info
    parms_run = list(model_json['parameters']['run'].values())[0]['parameters']
    parms_farm = list(model_json['parameters']['farm'].values())[0]['parameters']

    *especes, = model_json.\
        get('parameters', {}).\
        get('species', {})

    print(especes[0])

    stocking_density = model_json.\
        get('parameters', {}).\
        get('species', {}).\
        get(especes[0], {}).\
        get('parameters', {}).\
        get('stocking_density', 1000)

    #parms_harvest = list(model_json['parameters']['harvest'].values())[0]['parameters']
    #harvest_type = list(model_json['parameters']['harvest'].keys())[0]

    year = int(model_json['dataset_parameters']['year'])

    model = SF_model_scipy(model_json['parameters'])

    # Define beginning and end times of the simulation #change the name of variable 
    #startDate = datetime.datetime(year, int(parms_harvest['deployment_month']), 1)
    #endDate = datetime.datetime(year + 2, int(parms_harvest['harvesting_month']), 1) + relativedelta(months=1) 
    startDate = datetime.datetime(year, 1, 1)
    endDate = datetime.datetime(year + 1, 12, 1) + relativedelta(months=1) 

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
    init_grid_shape = working_data['Chlorophyll-a'].shape #the shape should be the same for all parameters
    longitudes, _ = input_data.parameterData['Chlorophyll-a'].getVariable('longitude', **data_kwargs)
    latitudes, _ = input_data.parameterData['Chlorophyll-a'].getVariable('latitude', **data_kwargs)
    #dxMeter, dyMeter = giveDxDy(latitudes, longitudes)

    #latRef = np.zeros((len(latitudes), len(longitudes)))
    #latRef[:, :] = latitudes[np.newaxis].T

    for par_name, par_data in input_data.parameterData.items():
        working_data[par_name] = np.ma.masked_outside(working_data[par_name], dataBounds[par_name][0],
                                                      dataBounds[par_name][1])
        working_data[par_name].filled(fill_value=0)
    
    mask = working_data['Chlorophyll-a'].mask

    grid_shape = init_grid_shape

    #nanLists = findNan(mask)

    for par_name, par_data in input_data.parameterData.items():
        print(par_name)
        print(working_data[par_name].shape)


    # Initialize the model variables
    state_vars = {
        'CHL': np.ma.masked_array(np.zeros(grid_shape), mask),
        'SHE': np.ma.masked_array(np.ones(grid_shape)*100, mask),
        'STE': np.ma.masked_array(np.ones(grid_shape)*1000, mask),
        'POP': np.ma.masked_array(np.ones(grid_shape)*stocking_density, mask),
        'spawnday': np.ma.masked_array(np.ones(grid_shape)*365, mask),
        'POC': np.ma.masked_array(np.ones(grid_shape)*100, mask),
        'PHYC': np.ma.masked_array(np.ones(grid_shape), mask),
        'sd2' : np.ma.masked_array(np.ones(grid_shape)*365, mask)
    }

    dt = 1 # days # TODO: make into parameter in json

    # Simulation loop
    sim_date = startDate
    while sim_date < endDate:
        print(f'{sim_date}')

        # Alter the date after new year to replicate data 
        data_date = sim_date.replace(year = year)

        # For each dataset, if the nearest i has changed, update the working data
        for par_name, par_data in input_data.parameterData.items():
            new_i = iNearest(data_date, time_axes[par_name])
            if new_i != nearest_time_i[par_name]:
                nearest_time_i[par_name] = new_i
                working_data[par_name], _ = par_data.getVariable(time_index=new_i, **data_kwargs)
                working_data[par_name] = np.ma.masked_outside(working_data[par_name], dataBounds[par_name][0],
                                                              dataBounds[par_name][1])
                working_data[par_name].filled(fill_value=0)

        days = (data_date - datetime.datetime(year, 1, 1)).days # Note: returns an integer only, that is sufficient precision for this

        bgc_terms = bgc_model(state_vars=state_vars, working_data=working_data,
                              model=model)

        for var_name in state_vars.keys():
            state_vars[var_name] += bgc_terms[var_name] * dt

        #state_vars = spawn_evenfunc(state_vars=state_vars)

        sim_date += datetime.timedelta(days = dt)

    

    output_data = output_dict(state_vars=state_vars, working_data=working_data,model=model)

    unitsDict = {"DSTW": 'kg/m^2', "STE": 'MJ/m^2', "FW": 'kg/m^2', "DWW": 'kg/m^2', "SHL": 'cm',
                 "NH4_production": 'g(N)/day', "CO2_production": 'g(C)/day'}
    # Create output file
    initialize_result(out_file_name, times=[0], latitudes=latitudes, longitudes=longitudes,
                      variableNames=["DSTW", "STE", "FW", "DWW", "SHL", "NH4_production", "CO2_production"],
                      unitsDict=unitsDict, mask=mask)

    # Write values to file
    ds = nc.Dataset(out_file_name, 'a')
    for name in model.names:
        ds[name][0,:,:] = np.ma.masked_array(output_data[name], mask)
    ds.close()
    return time.time() - t_init


#return the variation of each quantities in the algae model
def bgc_model(state_vars: dict, working_data: dict, model):

    data_in = {
        'SST': working_data['Temperature'],
        'CHL_ext': working_data['Chlorophyll-a'],
        'F_in': np.sqrt(working_data['northward_Water_current'] ** 2 + working_data['eastward_Water_current'] ** 2),
        't_z':  1.4, #(1 + parms_run['Von_Karman']) * model._parameters["z"]
        'POC': 100,
        'PHYC': state_vars['PHYC'],
    }

    y = {
        'CHL': state_vars['CHL'],
        'SHE': state_vars['SHE'],
        'STE': state_vars['STE'],
        'spawnday': state_vars['spawnday'],
        'sd2': state_vars['sd2'],
        'POP': state_vars['POP'],
        'POC': 100,
    }

    terms_list = model.derivative(y, data_in, model)
    all_terms = dict(zip(["CHL", "SHE", "STE", "spawnday", "sd2", 'POP', 'POC', 'PHYC'], terms_list))

    return all_terms

def output_dict(state_vars: dict, working_data: dict, model):
    # a remplacer pour t_z = (1 + parms_run['Von_Karman']) * model._parameters["z"]
    data_in = {
        'SST': working_data['Temperature'],
        'CHL_ext': working_data['Chlorophyll-a'],
        'F_in': np.sqrt(working_data['northward_Water_current'] ** 2 + working_data['eastward_Water_current'] ** 2),
        't_z': 1.4 
    }

    y = {
        'CHL': state_vars['CHL'],
        'SHE': state_vars['SHE'],
        'STE': state_vars['STE'],
        'spawnday': state_vars['spawnday'],
        'sd2': state_vars['sd2'],
        'POP': state_vars['POP'],
    }
    output_data = model.get_output(y,data_in, model)
    return output_data

def spawn_evenfunc(state_vars: dict):
    y = {
        'CHL': state_vars['CHL'],
        'SHE': state_vars['SHE'],
        'STE': state_vars['STE'],
        'spawnday': state_vars['spawnday'],
        'sd2': state_vars['sd2'],
        'POP': state_vars['POP'],
    }
    spawn_data = model.spawn_evenfunc(y, model)
    return(spawn_data)

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
