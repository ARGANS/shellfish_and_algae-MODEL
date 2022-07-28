import datetime
import numpy as np
import pandas as pd
import time
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import multiprocessing as mp
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
        var[:,:,:] = np.ma.masked_array(-1*np.ones(full_mask.shape), full_mask)

    ds.close()


def run_scenario_a_monthly(fileName:str, model_params:dict, y0:list, input_args:dict, year:int,
                           farm_at_gridSize=False, data_is_monthly=False):
    """Runs simulation on all the grid points in fileName, reading data from
    inputData, and initializing at y0. The simulations are ran on monthly
    averaged data.
    Writes the results in fileName.
    input_args is passed to AllData() to open an AllData object.

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

    algaeData = AllData(input_args)

    n_cells = 0

    ds = nc.Dataset(fileName, 'r')
    mask = ds['NH4'][0,:,:].mask
    latitudes = ds['latitude'][:]
    longitudes = ds['longitude'][:]
    times = ds['time'][:]
    ds.close()

    y0_array = np.tile(np.expand_dims(y0,(1,2)), (1, len(latitudes), len(longitudes)))
    full_mask = np.repeat(np.expand_dims(mask, 0), len(y0), axis=0)

    initTime = datetime.datetime(year, 1, 1, 0) # still initiating at 1st january for the light scheme

    gridLat = latitudes[1] - latitudes[0]
    gridLon = longitudes[1] - longitudes[0]

    parms_run = list(model_params['run'].values())[0]['parameters']

    for i_month, month in enumerate(times): # month is integer 1-12
        month = int(month)
        #print(f"MONTH: {month}")
        startTime = datetime.datetime(year, month, 1, 0)
        if month == 12:
            endTime = datetime.datetime(year+1, 1, 1, 0)
        else:
            endTime = datetime.datetime(year, month+1, 1, 0)

        values = np.ma.masked_array(-1*np.ones(y0_array.shape), full_mask)

        days_start = (startTime - initTime).days
        days_end = (endTime - initTime).days

        for i, lat in enumerate(latitudes):
            #print(f"Latitude: {lat}")

            # The model is recreated at every latitude, not too slow
            model = MA_model_scipy(model_params)
            if farm_at_gridSize:
                lon_m, lat_m = degrees_to_meters(gridLon, gridLat, lat)
                square_side = np.sqrt(lon_m * lat_m)
                model._parameters.update({'y_farm': square_side,
                                          'x_farm': square_side
                                        })
                # Ths size does not change with longitude

            for j, lon in enumerate(longitudes):
                #print(f"    Longitude: {lon}")
                if mask[i, j]:
                    continue

                if i_month==0:
                    n_cells += 1

                data_kwargs = {
                    "longitude": lon,
                    "latitude": lat,
                    "depth": (0, (1 + parms_run['Von_Karman']) * model._parameters["z"])
                }
                if data_is_monthly:
                    data_kwargs["time"] = startTime + datetime.timedelta(days=14)
                    data_kwargs["averagingDims"] = ("depth",)
                else:
                    data_kwargs["time"] = (startTime, endTime)
                    data_kwargs["averagingDims"] = ("time", "depth")

                data, dims = algaeData.getData(**data_kwargs)

                # The nearest data for one parameter is masked, ignore this point from now on
                if np.ma.core.MaskedConstant() in data.values():
                    mask[i,j] = True
                    continue

                data_in = {
                    'SST': data['Temperature'],
                    'PAR': data['PAR'],
                    'NH4_ext': data['Ammonium'],
                    'NO3_ext': data['Nitrate'],
                    'PO4_ext': data['Phosphate'],
                    'K_d': parms_run['K_d490'],
                    'F_in': np.sqrt(data['northward_Water_current']**2 + data['eastward_Water_current']**2),
                    't_z': (1 + parms_run['Von_Karman']) * model._parameters["z"],
                    'D_ext': parms_run['Detritus']
                }

                result = solve_ivp(MA_model_scipy.derivative_fast, (days_start, days_end), y0_array[:,i,j], args=(data_in, lat, model),
                                jac=MA_model_scipy.jacobian_fast,
                                rtol=0.05, method='BDF')

                if not result.success:
                    # try again with strict tolerance
                    result = solve_ivp(MA_model_scipy.derivative_fast, (days_start, days_end), y0_array[:,i,j], args=(data_in, lat, model),
                                jac=MA_model_scipy.jacobian_fast,
                                method='BDF')

                    if not result.success:
                        # Still failing
                        # Do not attempt to write if it failed
                        print(result.message)
                        continue

                values[:,i,j] = np.squeeze(result.y[:,-1])

        # Write values to file
        ds = nc.Dataset(fileName, 'a')
        for k, name in enumerate(model.names):
            ds[name][i_month,:,:] = np.ma.masked_array(values[k,:,:], mask)
        ds.close()

        # pass result as y0 for next step
        y0_array[:,:,:] = values[:,:,:]

    return n_cells


def run_scenario_a(fileName:str, model, y0:list, input_args:dict):
    """Runs simulation on all the grid points in fileName, reading data from
    inputData, and initializing at y0. Writes the results in fileName.
    input_args are passed to open_data_input() to open an AllData object.
    """

    algaeData = open_data_input(**input_args)

    n_cells = 0

    ds = nc.Dataset(fileName, 'r')
    mask = ds['NH4'][0,:,:].mask
    latitudes = ds['latitude'][:]
    longitudes = ds['longitude'][:]
    ds.close()

    interpKind = "nearest"

    for i, lat in enumerate(latitudes):
        #print(f"Latitude: {lat}")
        for j, lon in enumerate(longitudes):
            #print(f"    Longitude: {lon}")
            if mask[i, j]:
                continue

            n_cells += 1

            input_data = algaeData.getTimeSeries(lat, lon, (startDate, endDate), 3)

            time_axis = [(date - input_data['date'][0]).days for date in input_data['date']]
            data_fun = {
                'SST': interp1d(time_axis, input_data['Temperature'], kind=interpKind, assume_sorted=True),
                'PAR': lambda _ : 500,
                'NH4_ext': interp1d(time_axis, input_data['Ammonium'], kind=interpKind, assume_sorted=True),
                'NO3_ext': interp1d(time_axis, input_data['Nitrate'], kind=interpKind, assume_sorted=True),
                'PO4_ext': lambda _ : 50,
                'K_d': lambda _ : 0.1,
                'F_in': interp1d(time_axis, np.sqrt(input_data['northward_Water_current']**2 + input_data['eastward_Water_current']**2), kind=interpKind, assume_sorted=True),
                'h_z_SML': lambda _ : 30,
                't_z': lambda _ : 10,
                'D_ext': lambda _ : 0.1
            }

            result = solve_ivp(MA_model_scipy.derivative, (0, time_axis[-1]), y0, args=(data_fun, lat, model),
                            jac=MA_model_scipy.jacobian, t_eval=time_axis,
                            #rtol=0.05, atol=[0.5, 0.5, 100, 100, 0.04], method='BDF')
                            #atol=[0.5, 0.5, 100, 100, 0.04], method='BDF')
                            rtol=0.05, method='BDF')
                            #method='BDF')

            if not result.success:
                # Do not attempt to write if it failed
                print(result.message)
                continue

            ds = nc.Dataset(fileName, 'a')
            for k, name in enumerate(model.names):
                ds[name][:,i,j] = result.y[k,:]
            ds.close()

    return n_cells


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
