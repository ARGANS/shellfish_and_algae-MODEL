import datetime
import numpy as np
from read_netcdf import *
from launch_model import *
import time


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
        'longitude': (-4, -3),
        'latitude': (48, 50),
        'time_index': 0,
        'depth_index': 0
    }

    longitudes = algaeData.parameterData['Temperature'].getVariable('longitude', **sim_area)
    latitudes = algaeData.parameterData['Temperature'].getVariable('latitude', **sim_area)
    mask = algaeData.parameterData['Temperature'].getVariable(**sim_area).mask

    startDate = datetime.datetime(2021, 1, 1, 12)
    endDate = datetime.datetime(2022, 1, 1, 12)

    model = MA_model("macroalgae_model.R", "macroalgae_model_parameters.json")

    result = np.zeros(mask.shape)
    result = np.ma.masked_array(result, mask)

    t0 = time.time()
    for i, lat in enumerate(latitudes):
        print(f"LATITUDE: {lat}")
        for j, lon in enumerate(longitudes):
            print(f"lon: {lon}")
            if mask[i, j] :
                continue

            # Get data at 3m
            df = algaeData.getTimeSeries(lat, lon, (startDate, endDate), 3)

            # This translation should be done somewhere else ?
            dataToR = pd.DataFrame({
                'time': [(date - df['date'][0]).days + 1 for date in df['date']],
                'SST': df['Temperature'],
                'PAR': 500,
                'NH4_in': df['Ammonium'],
                'NO3_in': df['Nitrate'],
                'PO4_in': 50,
                'K_d': 0.1,
                'F_in': np.sqrt(df['northward_Water_current']**2 + df['eastward_Water_current']**2),
                'h_z_SML': 30,
                't_z': 0,
                'D_in': 0.1
                })
            print(dataToR)

            out = model.apply_on(dataToR, float(lat))

            result[i,j] = np.sum(out['f_NO3'])

    print(result)
    print(f"time taken: {t0-time.time()}")
