import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import STAP
from rpy2.robjects.conversion import localconverter
import pandas as pd
import datetime
import json
import time
#from scipy.integrate import odeint

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






if __name__ =="__main__":

    # Read from CSV for now
    bantryData = pd.read_csv('bantry_data/bantry_3m.csv', delimiter=';')
    bantryData['date'] = [datetime.datetime.fromisoformat(str_date) for str_date in bantryData['date']]

    dataToR = pd.DataFrame({
        'time': [(date - bantryData['date'][0]).days + 1 for date in bantryData['date']],
        'SST': 15,
        'PAR': bantryData['par'],
        'NH4_in': bantryData['Ammonium'],
        'NO3_in': bantryData['Nitrate'],
        'PO4_in': 50,
        'K_d': 0.1,
        'F_in': 100,
        'h_z_SML': 30,
        't_z': 0,
        'D_in': 0.1
    })

    dataToR = pd.DataFrame({
        'time': [0, 1],
        'SST': 15,
        'PAR': bantryData['par'][100],
        'NH4_in': bantryData['Ammonium'][100],
        'NO3_in': bantryData['Nitrate'][100],
        'PO4_in': 50,
        'K_d': 0.1,
        'F_in': 100,
        'h_z_SML': 30,
        't_z': 0,
        'D_in': 0.1
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



