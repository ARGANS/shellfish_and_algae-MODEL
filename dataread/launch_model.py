import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import STAP
from rpy2.robjects.conversion import localconverter
import pandas as pd
import datetime
import json


if __name__ =="__main__":

    with open('P:/Aquaculture/shellfish_and_algae-MODEL/macroalgae/run_MA.R', 'r') as f:
        body = f.read()
        run_MA = STAP(body, "run_model")

    # Read from CSV for now
    bantryData = pd.read_csv('P:/Aquaculture/shellfish_and_algae-MODEL/macroalgae/bantry_data/bantry_MLDaveraged.csv', delimiter=';')
    bantryData['date'] = [datetime.datetime.fromisoformat(str_date) for str_date in bantryData['date']]

    dataToR = pd.DataFrame({
        'time': [(date - bantryData['date'][0]).days for date in bantryData['date']],
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

    print(dataToR)

    with open('P:/Aquaculture/shellfish_and_algae-MODEL/macroalgae/model_parameters.json', 'r') as f:
        default_parms = json.loads(f.read())


    parms_species = default_parms['algae']['species']['ulva']['parameters']
    parms_farm = default_parms['algae']['farm_parameters']
    parms_run = default_parms['algae']['run_parameters']

    parms_run['light_scheme'] = 2

    parms_species = robjects.ListVector(parms_species)
    parms_farm = robjects.ListVector(parms_farm)
    parms_run = robjects.ListVector(parms_run)
    parms_additional = robjects.ListVector({"latitude": 55.})

    # get the R concatenator function
    conc = robjects.r('c')

    parameters = conc(parms_species, parms_farm, parms_run, parms_additional)

    y0 = robjects.ListVector({"NH4": 0,
                              "NO3": 0,
                              "N_s": 1,
                              "N_f": 1,
                              "D": 0,
                              "Yield": 0
                            })

    #print(parms)

    # Using a converter between R data.frame and python pd.DataFrame
    with localconverter(robjects.default_converter + pandas2ri.converter):
       out = run_MA.run_model(input=dataToR, parameters=parameters, y0=y0)

    print(out)



