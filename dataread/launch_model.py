import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import STAP
from rpy2.robjects.conversion import localconverter
import pandas as pd
import datetime


if __name__=="__main__":

    with open('P:/Aquaculture/shellfish_and_algae-MODEL/macroalgae/run_MA.R', 'r') as f:
        body = f.read()
        #body = body.replace('\r\n', '\n')
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


    with localconverter(robjects.default_converter + pandas2ri.converter):
       out = run_MA.run_model(input = dataToR)

    print(out)



