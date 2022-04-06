import datetime
import numpy as np
from read_netcdf import *
from launch_model import *
import time
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import multiprocessing as mp


def yearly_runs(data, latitudes, longitude, startDate, endDate, mask, model):
    result = np.zeros(mask.shape)
    result = np.ma.masked_array(result, mask)

    n_cells = 0
    time_model = 0
    T0 = time.time()

    for i, lat in enumerate(latitudes):
        print(f"LATITUDE: {lat}")
        for j, lon in enumerate(longitudes):
            if mask[i, j] :
                continue

            print(f"lon: {lon}")

            n_cells += 1

            # Get data at 3m
            df = algaeData.getTimeSeries(lat, lon, (startDate, endDate), 3)

            # This translation should be done somewhere else ?
            dataToR = pd.DataFrame({
                'time': [(date - df['date'][0]).days + 1 for date in df['date']],
                'SST': df['Temperature'],
                'PAR': 500,
                'NH4_ext': df['Ammonium'],
                'NO3_ext': df['Nitrate'],
                'PO4_ext': 50,
                'K_d': 0.1,
                'F_in': np.sqrt(df['northward_Water_current']**2 + df['eastward_Water_current']**2),
                'h_z_SML': 30,
                't_z': 10,
                'D_ext': 0.1
                })

            t0 = time.time()

            #out = model.apply_on(dataToR, float(lat))
            y0 = pd.DataFrame({"NH4": [0],
                               "NO3": 0,
                               "N_s": 1,
                               "N_f": 1,
                               "D": 0,
                               "Yield": 0,
                               "Yield_per_m": 0})

            out_rows = []
            for t in dataToR.index:

                # twice the same row
                input_row = dataToR.iloc[[t,t]].reset_index()

                # start at 0, count the number of days
                #ndays = dataToR['time'][t] - dataToR['time'][t-1] if t!=0 else dataToR['time'][0]
                #input_row['time'] = [0, ndays]
                input_row.loc[1, ['time']] = dataToR['time'][t+1] if t+1 in dataToR.index else dataToR['time'][t] + 1

                row = model.apply_on(input_row, float(lat), y0).iloc[[1]]

                out_rows.append(row)

                # change y0 for next step
                y0 = row[["NH4", "NO3", "N_s", "N_f", "D", "Yield", "Yield_per_m"]]
            
            out = pd.concat(out_rows, ignore_index=True)

            time_model += time.time() - t0

            result[i,j] = np.sum(out['f_NO3'])

    print(result)
    print(f"Average model run: {time_model/n_cells}\n" +
          f"Total time: {time.time() - T0}"
        )


def monthly_simulations(data, latitudes, longitude, startDate, endDate, mask, model):
    ### Create the necdf output
    ds = nc.Dataset("/media/share/results/monthly_simulations_test.nc", 'w', format='NETCDF4')

    months = ds.createDimension('time', None)
    lat = ds.createDimension('latitude', len(latitudes))
    lon = ds.createDimension('longitude', len(longitudes))

    lats = ds.createVariable('latitude', 'f4', ('latitude',))
    lons = ds.createVariable('longitude', 'f4', ('longitude',))

    lats[:] = latitudes
    lons[:] = longitudes



    for name in ["NH4", "NO3", "N_s", "N_f", "D", "Yield", "Yield_per_m"]:
        var = ds.createVariable(name, 'f4', ('time', 'latitude', 'longitude',))

    ds.close()


    n_cells = 0
    T0 = time.time()

    initTime = datetime.datetime(2021, 1, 1, 0)

    results = {}
    results_names = ["NH4", "NO3", "N_s", "N_f", "D", "Yield", "Yield_per_m"]
    for param in results_names:
        zeros = np.zeros(mask.shape)
        results[param] = np.ma.masked_array(zeros, mask)

    y0 = pd.DataFrame({"NH4": [0],
                       "NO3": 0,
                       "N_s": 1,
                       "N_f": 1,
                       "D": 0,
                       "Yield": 0,
                       "Yield_per_m": 0})

    y0_array = np.empty(mask.shape, dtype=object)

    # !!! All the same reference, but should not matter
    for i in range(mask.shape[0]):
        for j in range(mask.shape[1]):
            y0_array[i,j] = y0

    for month in range(1,13):
        print(f"MONTH: {month}")
        startTime = datetime.datetime(2021, month, 1, 0)
        if month == 12:
            endTime = datetime.datetime(2022, 1, 1, 0)
        else:
            endTime = datetime.datetime(2021, month+1, 1, 0)
        data, dims = algaeData.getData(longitude = sim_area['longitude'], 
                                       latitude = sim_area['latitude'],
                                       depth = sim_area['depth'],
                                       time = (startTime, endTime),
                                       averagingDims = ('time',),
                                       weighted = False
                                       )

        for i, lat in enumerate(latitudes):
            print(f"    latitude: {lat}")
            for j, lon in enumerate(longitudes):
                if mask[i, j]:
                    continue
                n_cells += 1
                #print(f"        longitude: {lon}")

                # This translation should be done somewhere else ?
                dataToR = pd.DataFrame({
                    'time': [(startTime - initTime).days, (endTime - initTime).days],
                    'SST': data['Temperature'][i,j],
                    'PAR': 500,
                    'NH4_ext': data['Ammonium'][i,j],
                    'NO3_ext': data['Nitrate'][i,j],
                    'PO4_ext': 50,
                    'K_d': 0.1,
                    'F_in': np.sqrt(data['northward_Water_current'][i,j]**2 + data['eastward_Water_current'][i,j]**2),
                    'h_z_SML': 30,
                    't_z': 10,
                    'D_ext': 0.1
                    })

                # run the month
                model_res = model.apply_on(dataToR, float(lat), y0_array[i,j]).iloc[[1]]

                for param in results_names:
                    results[param][i,j] = model_res[param][0]

                y0_array[i,j] = model_res[["NH4", "NO3", "N_s", "N_f", "D", "Yield", "Yield_per_m"]]

        ds = nc.Dataset("/media/share/results/monthly_simulations_test.nc", 'a')
        for name, values in results.items():
            ds[name][month-1,:,:] = values
        ds.close()

    total_time = time.time() - T0
    print(f"Total time: {total_time}" +
          f"Time per cell per month: {total_time/n_cells}"
    )


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


def run_scenario_a_monthly(fileName:str, model, y0:list, input_args:dict, year:int):
    """Runs simulation on all the grid points in fileName, reading data from
    inputData, and initializing at y0. The simulations are ran on monthly
    averaged data.
    Writes the results in fileName.
    input_args are passed to open_data_input() to open an AllData object.
    """

    algaeData = open_data_input(**input_args)

    n_cells = 0

    ds = nc.Dataset(fileName, 'r')
    mask = ds['NH4'][0,:,:].mask
    latitudes = ds['latitude'][:]
    longitudes = ds['longitude'][:]
    ds.close()

    y0_array = np.tile(np.expand_dims(y0,(1,2)), (1, len(latitudes), len(longitudes)))
    full_mask = np.repeat(np.expand_dims(mask, 0), len(y0), axis=0)

    initTime = datetime.datetime(year, 1, 1, 0)

    for month in range(1,13):
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
            for j, lon in enumerate(longitudes):
                #print(f"    Longitude: {lon}")
                if mask[i, j]:
                    continue

                if month==1:
                    n_cells += 1

                data, dims = algaeData.getData(longitude = lon,
                                       latitude = lat,
                                       depth = 3,
                                       time = (startTime, endTime)
                                        )

                data_in = {
                    'SST': np.average(data['Temperature']),
                    'PAR': 500,
                    'NH4_ext': np.average(data['Ammonium']),
                    'NO3_ext': np.average(data['Nitrate']),
                    'PO4_ext': 50,
                    'K_d': 0.1,
                    'F_in': np.average(np.sqrt(data['northward_Water_current']**2 + data['eastward_Water_current']**2)),
                    'h_z_SML': 30,
                    't_z': 10,
                    'D_ext': 0.1
                }

                result = solve_ivp(MA_model_scipy.derivative_fast, (days_start, days_end), y0_array[:,i,j], args=(data_in, lat, model),
                                jac=MA_model_scipy.jacobian_fast, t_eval=[days_end],
                                rtol=0.05, method='BDF')

                if not result.success:
                    # Do not attempt to write if it failed
                    print(result.message)
                    continue

                values[:,i,j] = np.squeeze(result.y)

        # Write values to file
        ds = nc.Dataset(fileName, 'a')
        for k, name in enumerate(model.names):
            ds[name][month-1,:,:] = values[k,:,:]
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


def open_data_input(mainpath:str, zone:str, paramNames:list, fileRef:str):

    dataRef = pd.read_csv(fileRef, delimiter=';')

    fileNames = [mainpath + f"{zone}/{param}/{zone}_{param}_merged.nc" for param in paramNames]

    #dataRows = [dataRef.index[(dataRef['Parameter']==param) & (dataRef['Place']==zone)][0] for param in paramNames]
    dataRows = [dataRef.index[(dataRef['Parameter']==param) & (dataRef['Place']==zone)][-1] for param in paramNames]

    variableNames = dataRef.iloc[dataRows]['variable'].tolist()
    latitudeNames = dataRef.iloc[dataRows]['latName'].fillna('latitude').tolist()
    longitudeNames = dataRef.iloc[dataRows]['longName'].fillna('longitude').tolist()
    timeNames = dataRef.iloc[dataRows]['timeName'].fillna('time').tolist()
    depthNames = dataRef.iloc[dataRows]['depthName'].fillna('depth').tolist()
    unitConversions = dataRef.iloc[dataRows]['unitFactor'].fillna(1).tolist()

    data = AllData(fileNameList=fileNames,
                   parameterNameList=paramNames,
                   variableNameList=variableNames,
                   latitudeNameList=latitudeNames,
                   longitudeNameList=longitudeNames,
                   timeNameList=timeNames,
                   depthNameList=depthNames,
                   unitConversionList=unitConversions
    )

    return data


if __name__=="__main__":

    input_args = {
        'zone' : "IBI",
        'mainpath' : '/media/share/data_merged/',
        'fileRef' : './dataCmd.csv',
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

    startDate = datetime.datetime(2021, 1, 1, 12)
    endDate = datetime.datetime(2022, 1, 1, 12)

    longitudes, _ = algaeData.parameterData['Temperature'].getVariable('longitude', **sim_area)
    latitudes, _ = algaeData.parameterData['Temperature'].getVariable('latitude', **sim_area)
    times, _ = algaeData.parameterData['Temperature'].getVariable('time', time=(startDate, endDate), rawTime=True)

    mask1 = algaeData.parameterData['Temperature'].getVariable(**sim_area)[0].mask
    mask2 = algaeData.parameterData['eastward_Water_current'].getVariable(**sim_area)[0].mask
    mask3 = algaeData.parameterData['northward_Water_current'].getVariable(**sim_area)[0].mask

    mask = np.logical_or(mask1, np.logical_or(mask2, mask3))

    model = MA_model_scipy("macroalgae_model_parameters.json")

    n_slices = 20

    lon_split = np.array_split(longitudes, n_slices)
    mask_split = np.array_split(mask, n_slices, 1)

    ### Create datasets
    for i, (lon_arr, mask_arr) in enumerate(zip(lon_split, mask_split)):
        initialize_result(f"/media/share/results/complete_simulations_test_{i}.nc",
        times, latitudes, lon_arr, model.names, mask_arr)
    ### Create datasets for monthly sim
    for i, (lon_arr, mask_arr) in enumerate(zip(lon_split, mask_split)):
        initialize_result(f"/media/share/results/complete_simulations_monthly_test_{i}.nc",
        np.array(range(1,13)), latitudes, lon_arr, model.names, mask_arr)

    pool = mp.Pool(10)

    y0 = np.array([0, 0, 0, 1000, 0], dtype=np.float64)


    t0 = time.time()
    #n_cells = pool.starmap_async(run_scenario_a, 
    #    [(f"/media/share/results/complete_simulations_test_{i}.nc", model, y0, input_args) for i in range(n_slices)]).get()
    n_cells = pool.starmap_async(run_scenario_a_monthly, 
        [(f"/media/share/results/complete_simulations_monthly_test_{i}.nc", model, y0, input_args, 2021) for i in range(n_slices)]).get()
    #n_cells = run_scenario_a_monthly("/media/share/results/complete_simulations_monthly_test_0.nc", model, y0, input_args, 2021)


    print(n_cells)
    print((time.time()-t0)/sum(n_cells))