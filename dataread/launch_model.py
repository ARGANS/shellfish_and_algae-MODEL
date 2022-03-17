import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import STAP
from rpy2.robjects.conversion import localconverter
import pandas as pd
import datetime
import json

class MA_model:

    def __init__(self, model_path, parms_path):

        # Opens model_path to read the run_MA_model function
        with open(model_path, 'r') as f:
            body = f.read()
            self._model = STAP(body, "run_MA_model")

        self._y0 = robjects.r('c("NH4" = 0, \
                                 "NO3" = 0, \
                                 "N_s" = 1, \
                                 "N_f" = 1, \
                                 "D" = 0, \
                                 "Yield" = 0,\
                                 "Yield_per_m" = 0)')

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


    def apply_on(self, input_data, latitude):
        # Applies the model on the input data (pd.DataFrame), latitude has to be specified

        parms_with_lat = self._conc(self._parameters, robjects.ListVector({"latitude": latitude}))

        # Using a converter between R data.frame and python pd.DataFrame
        with localconverter(robjects.default_converter + pandas2ri.converter):
            out = self._model.run_MA_model(input = dataToR,
                                           parameters = parms_with_lat,
                                           y0 = self._y0)

        return out



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

    model = MA_model("macroalgae_model.R", "macroalgae_model_parameters.json")

    out = model.apply_on(dataToR, 51.)

    print(out)



