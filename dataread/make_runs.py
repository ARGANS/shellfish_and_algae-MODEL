import datetime
import numpy as np
from read_netcdf import *
from launch_model import *
import time
from scipy.integrate import odeint


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



if __name__=="__main__":

    zone = "IBI"

    mainpath = '/media/share/data_merged/'

    #dataRef = pd.read_csv('/profils/qjutard/shellfish_and_algae-MODEL/dataimport/src/dataCmd.csv', delimiter=';')
    dataRef = pd.read_csv('./dataCmd.csv', delimiter=';')

    ### Initialize the netcdf reading interface

    #paramNames = ['Ammonium', 'Nitrate', 'Temperature', 'northward_Water_current', 'eastward_Water_current', 'ocean_mixed_layer_thickness', 'par']
    paramNames = ['Ammonium', 'Nitrate', 'Temperature', 'northward_Water_current', 'eastward_Water_current']
    fileNames = [mainpath + f"{zone}/{param}/{zone}_{param}_merged.nc" for param in paramNames]

    #dataRows = [dataRef.index[(dataRef['Parameter']==param) & (dataRef['Place']==zone)][0] for param in paramNames]
    dataRows = [dataRef.index[(dataRef['Parameter']==param) & (dataRef['Place']==zone)][-1] for param in paramNames]

    variableNames = dataRef.iloc[dataRows]['variable'].tolist()
    latitudeNames = dataRef.iloc[dataRows]['latName'].fillna('latitude').tolist()
    longitudeNames = dataRef.iloc[dataRows]['longName'].fillna('longitude').tolist()
    timeNames = dataRef.iloc[dataRows]['timeName'].fillna('time').tolist()
    depthNames = dataRef.iloc[dataRows]['depthName'].fillna('depth').tolist()
    unitConversions = dataRef.iloc[dataRows]['unitFactor'].fillna(1).tolist()

    algaeData = AllData(fileNameList=fileNames,
                        parameterNameList=paramNames,
                        variableNameList=variableNames,
                        latitudeNameList=latitudeNames,
                        longitudeNameList=longitudeNames,
                        timeNameList=timeNames,
                        depthNameList=depthNames,
                        unitConversionList=unitConversions
    )


    ### get the copernicus grid and mask

    sim_area = {
        'longitude': (-4, -3.5),
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

    startDate = datetime.datetime(2021, 1, 1, 12)
    endDate = datetime.datetime(2022, 1, 1, 12)

    model = MA_model_scipy("macroalgae_model.R", "macroalgae_model_parameters.json")

    for i, lat in enumerate(latitudes):
        print(f"Latitude: {lat}")
        for j, lon in enumerate(longitudes):
            print(f"    Longitude: {lon}")
            if mask[i, j]:
                continue


