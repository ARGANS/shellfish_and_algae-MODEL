import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import STAP
from rpy2.robjects.conversion import localconverter
import pandas as pd
import numpy as np
import datetime
import json
import time
from scipy.integrate import odeint

from read_netcdf import *


class MA_model:

    def __init__(self, model_path, parms_path):

        # Opens model_path to read the run_MA_model function
        with open(model_path, 'r') as f:
            body = f.read()
            self._model = STAP(body, "run_MA_model")

        ### This section will need an update to use the input json rather than defaults json
        with open(parms_path, 'r') as f:
            parms = json.loads(f.read())

        param_lists = []
        for _,section in parms.items():
            # just take the first option for each section for now
            first_default = section['defaults'][next(iter(section['defaults']))]

            param_lists.append(robjects.ListVector(first_default['parameters']))
            param_lists.append(robjects.ListVector(first_default['options']))

        # get the R concatenator function
        self._conc = robjects.r('c')

        # concatenate all parameter lists
        self._parameters = self._conc(*param_lists)


    def apply_on(self, input_data, latitude, y0=robjects.r('c("NH4" = 0, \
                                                              "NO3" = 0, \
                                                              "N_s" = 1, \
                                                              "N_f" = 1, \
                                                              "D" = 0, \
                                                              "Yield" = 0,\
                                                              "Yield_per_m" = 0)')):
        # Applies the model on the input data (pd.DataFrame), latitude has to be specified

        parms_with_lat = self._conc(self._parameters, robjects.ListVector({"latitude": latitude}))

        # Using a converter between R data.frame and python pd.DataFrame
        with localconverter(robjects.default_converter + pandas2ri.converter):
            out = self._model.run_MA_model(input = input_data,
                                           parameters = parms_with_lat,
                                           y0 = y0)

        return out


class MA_model_scipy:
    def __init__(self, model_path, parms_path):

        # Opens model_path to read the run_MA_model function
        with open(model_path, 'r') as f:
            body = f.read()
            self._model = STAP(body, "MA_model_const")
            self._jacobian = STAP(body, "MA_model_jacobian")

        ### This section will need an update to use the input json rather than defaults json
        with open(parms_path, 'r') as f:
            parms = json.loads(f.read())

        param_lists = []
        for _,section in parms.items():
            # just take the first option for each section for now
            first_default = section['defaults'][next(iter(section['defaults']))]

            param_lists.append(robjects.ListVector(first_default['parameters']))
            param_lists.append(robjects.ListVector(first_default['options']))

        # get the R concatenator function
        self._conc = robjects.r('c')
        # concatenate all parameter lists
        self._parameters = self._conc(*param_lists)

        self._names = ["NH4", "NO3", "N_s", "N_f", "D"]

        self.time_reading = 0
        self.time_computing = 0
        self.time_reading_J = 0
        self.time_computing_J = 0


    @staticmethod
    def derivative(y: np.array, t: float, data: pd.DataFrame, latitude: float, 
                   longitude: float, startDate: datetime.datetime, self):

        t1 = time.time()
        yToR = pd.Series(y, index=self._names)

        iData = iNearest(t, data['time'], "closest")
        dataToR = data.iloc[[iData]]
        self.time_reading += time.time() - t1

        t1 = time.time()
        with localconverter(robjects.default_converter + pandas2ri.converter):
            out = self._model.MA_model_const(t = t,
                                             y = yToR,
                                             parms = self._parameters,
                                             data = dataToR,
                                             latitude = latitude
                                             )
        self.time_computing += time.time() - t1

        return np.squeeze(np.array(out))

    @staticmethod
    def jacobian(y: np.array, t: float, data: pd.DataFrame, latitude: float, 
                   longitude: float, startDate: datetime.datetime, self):
        t1 = time.time()
        yToR = pd.Series(y, index=self._names)

        iData = iNearest(t, data['time'], "closest")
        dataToR = data.iloc[[iData]]
        self.time_reading_J += time.time() - t1

        t1 = time.time()
        with localconverter(robjects.default_converter + pandas2ri.converter):
            out = self._jacobian.MA_model_jacobian(t = t,
                                                   y = yToR,
                                                   parms = self._parameters,
                                                   data = dataToR,
                                                   latitude = latitude
                                                   )
        self.time_computing_J += time.time() - t1
        return out



if __name__ =="__main__":

    # Read from CSV for now
    bantryData = pd.read_csv('bantry_data/bantry_3m.csv', delimiter=';')
    bantryData['date'] = [datetime.datetime.fromisoformat(str_date) for str_date in bantryData['date']]

    dataToR = pd.DataFrame({
        'time': [(date - bantryData['date'][0]).days + 1 for date in bantryData['date']],
        'SST': 15,
        'PAR': bantryData['par'],
        'NH4_ext': bantryData['Ammonium'],
        'NO3_ext': bantryData['Nitrate'],
        'PO4_ext': 50,
        'K_d': 0.1,
        'F_in': 100,
        'h_z_SML': 30,
        't_z': 10,
        'D_ext': 0.1
    })

    dataToR = pd.DataFrame({
        'time': [0, 1],
        'SST': 15,
        'PAR': bantryData['par'][100],
        'NH4_ext': bantryData['Ammonium'][100],
        'NO3_ext': bantryData['Nitrate'][100],
        'PO4_ext': 50,
        'K_d': 0.1,
        'F_in': 100,
        'h_z_SML': 30,
        't_z': 10,
        'D_ext': 0.1
    })

    y0 = pd.DataFrame({"NH4": [0],
                       "NO3": 0,
                       "N_s": 1,
                       "N_f": 1,
                       "D": 0,
                       "Yield": 0,
                       "Yield_per_m": 0})

    model = MA_model("macroalgae_model.R", "macroalgae_model_parameters.json")

    print(dataToR)

    out = model.apply_on(dataToR, 51., y0)

    print(out)


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

    startDate = datetime.datetime(2021, 1, 1, 12)
    endDate = datetime.datetime(2022, 1, 1, 12)

    model = MA_model_scipy("macroalgae_model.R", "macroalgae_model_parameters.json")

    input_data = algaeData.getTimeSeries(54, 0, (startDate, endDate), 3)
    data = pd.DataFrame({
                'time': [(date - input_data['date'][0]).days for date in input_data['date']],
                'SST': input_data['Temperature'],
                'PAR': 500,
                'NH4_ext': input_data['Ammonium'],
                'NO3_ext': input_data['Nitrate'],
                'PO4_ext': 50,
                'K_d': 0.1,
                'F_in': np.sqrt(input_data['northward_Water_current']**2 + input_data['eastward_Water_current']**2),
                'h_z_SML': 30,
                't_z': 10,
                'D_ext': 0.1
                    })

    y0 = np.array([0, 0, 1, 1, 0])

    dy0 = MA_model_scipy.derivative(y0, 0, data, 49, -3, startDate, model)
    J0 = MA_model_scipy.jacobian(y0, 0, data, 49, -3, startDate, model)
    print(dy0)
    print(J0)

    t0 = time.time()

    result = odeint(MA_model_scipy.derivative, y0, args=(data, 54, 0, startDate, model),
                    t=np.array([0, 365])

    print(f"One run: {time.time() - t0}\nReading: {model.time_reading}\nComputing (R): {model.time_computing}")

    t0 = time.time()

    result = odeint(MA_model_scipy.derivative, y0, args=(data, 54, 0, startDate, model),
                    Dfun = MA_model_scipy.jacobian, col_deriv=True,
                    t=np.array([0, 365])

    print(f"One run with jacobian: {time.time() - t0}\nReading: {model.time_reading}\nComputing (R): {model.time_computing}" +
          f"\nReading jacobian: {model.time_reading_J}\nComputing jacobian (R): {model.time_computing_J}"
           )

