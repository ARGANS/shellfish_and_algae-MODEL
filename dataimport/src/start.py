import os
import json
from general import giveDateslist, getData
import pandas as pd
from pprint import pprint

def cleanFinalSlash(value: str) -> str:
    return value[:-1] if value.endswith('/') else value

# Read dataset parameters
### parameters_json_value == {"zone":"Artic","depth_min":0,"depth_max": 0,"year": 2022}
parameters_json_value:str = os.getenv('INPUT_PARAMETERS') 
try:
    input_parameters:dict = json.loads(parameters_json_value)
except:
    raise RuntimeError('Cannot parse the value of the parameters_json environment variable')
		

year = int(input_parameters['year'])
zone = input_parameters['zone']
deepthmin = int(input_parameters['depth_min'])
deepthmax = int(input_parameters['depth_max'])


outputDirectory = cleanFinalSlash(os.getenv('INPUT_DESTINATION'))

wantedData = ['Temperature', 'Nitrate', 'Ammonium', 'Phosphate', 'eastward_Water_current', 
              'northward_Water_current', 'pCO2','disolved_inorganic_carbon','primary_producer_POC']

dateBeginning = f'{year}-01-01 00:00:00'
dateEnd = f'{year + 1}-01-01 00:00:00'
## Per month
# frequency = 1
## Per day
frequency = 'monthly'

dataFin = pd.read_csv('./dataCmd.csv',';')
datesList = giveDateslist(dateBeginning, dateEnd, frequency)

for dat in wantedData:
    dataOutputDirectory = outputDirectory + '/' + dat + '/'
    dataLine = dataFin.loc[dataFin["Parameter"] == dat]
    print(f'dataLine {dataOutputDirectory}')
    pprint(dataLine)
    if dataLine.iloc[0]["daily"] != 'permanent':
        if frequency == 'monthly':
            getData(dat, zone, dataFin, deepthmin, deepthmax,  dataOutputDirectory, datesList[0], datesList[1],frequency)
        elif frequency == 'daily':
            for (dateBeg, dateE) in zip(datesList[0], datesList[1]):
                getData(dat, zone, dataFin, deepthmin, deepthmax, dataOutputDirectory, dateBeg, dateE, frequency)
    else:
        getData(dat, zone, dataFin, deepthmin, deepthmax, dataOutputDirectory)
