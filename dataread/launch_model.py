import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import STAP
from rpy2.robjects.conversion import localconverter
import pandas as pd
import datetime
import json


if __name__ =="__main__":

    with open('P:/Aquaculture/shellfish_and_algae-MODEL/macroalgae/macroalgae_model.R', 'r') as f:
        body = f.read()
        model = STAP(body, "run_MA_model")

    # Read from CSV for now
    bantryData = pd.read_csv('P:/Aquaculture/shellfish_and_algae-MODEL/macroalgae/bantry_data/bantry_3m.csv', delimiter=';')
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

    with open('P:/Aquaculture/shellfish_and_algae-MODEL/macroalgae/model_parameters.json', 'r') as f:
        default_parms = json.loads(f.read())

    param_lists = [robjects.ListVector({"latitude": 55.})]
    for _,section in default_parms.items():
        # just take the first option for each section for now
        first_default = section['defaults'][next(iter(section['defaults']))]

        param_lists.append(robjects.ListVector(first_default['parameters']))
        param_lists.append(robjects.ListVector(first_default['options']))

    # get the R concatenator function
    conc = robjects.r('c')

    # concatenate all parameter lists
    parameters = conc(*param_lists)

    y0 = robjects.r('c("NH4" = 0, \
                       "NO3" = 0, \
                       "N_s" = 1, \
                       "N_f" = 1, \
                       "D" = 0, \
                       "Yield" = 0,\
                       "Yield_per_m" = 0)')


    # Using a converter between R data.frame and python pd.DataFrame
    with localconverter(robjects.default_converter + pandas2ri.converter):
       out = model.run_MA_model(input=dataToR, parameters=parameters, y0=y0)

    print(out)



