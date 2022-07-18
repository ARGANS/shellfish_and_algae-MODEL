import os
import json
from general import giveDateslist, getData
import pandas as pd
from pprint import pprint

def cleanFinalSlash(value: str) -> str:
    return value[:-1] if value.endswith('/') else value

# Read dataset parameters
parameters_json_value:str = os.getenv('INPUT_PARAMETERS') 
try:
    input_parameters:dict = json.loads(parameters_json_value)
except:
    raise RuntimeError('Cannot parse the value of the parameters_json environment variable')
		

# DEPRECATED
# since the script gets the dataset properties from an environment variable, the code below is no longer needed:
year = int(input_parameters['first obs'][-4:])
zone = input_parameters['Place']
deepthmin = int(input_parameters['depth-min'])
deepthmax = int(input_parameters['depth-max'])
type = input_parameters['type']
frequency = input_parameters['daily']

# TODO update the getData() function and use the input_parameters object which contains the data directly from dataCmd.csv

outputDirectory = cleanFinalSlash(os.getenv('INPUT_DESTINATION'))

wantedData = ['Temperature', 'Nitrate', 'Ammonium', 'Phosphate', 'eastward_Water_current', 
              'northward_Water_current', 'pCO2','disolved_inorganic_carbon','primary_producer_POC']

dateBeginning = f'{year}-01-01 00:00:00'
dateEnd = f'{year + 1}-01-01 00:00:00'
dataFin = pd.read_csv('/media/global/dataCmd.csv',';')
datesList = giveDateslist(dateBeginning, dateEnd, frequency)

for dat in wantedData:
    dataOutputDirectory = outputDirectory + '/' + dat + '/'
    dataLine = dataFin.loc[dataFin["Parameter"] == dat]

    if dataLine.iloc[0]['daily'] != 'permanent':
        for (dateBeg, dateE) in zip(datesList[0], datesList[1]):
            getData(dat, zone, dataFin, deepthmin, deepthmax, dataOutputDirectory, dateBeg, dateE, frequency=frequency, type=type)
    else:
        getData(dat, zone, dataFin, deepthmin, deepthmax, dataOutputDirectory, type=type)
