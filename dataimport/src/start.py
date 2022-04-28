import os
import json
from general import giveDateslist, getData
import pandas as pd
from pprint import pprint

# Read dataset parameters
### parameters_json_value == {"zone":"Artic","depth_min":0,"depth_max": 0,"year": 2022}
parameters_json_value:str = os.getenv('parameters_json') 
try:
    input_parameters:dict = json.loads(parameters_json_value)
except:
    raise RuntimeError('Cannot parse the value of the parameters_json environment variable')
		

year = input_parameters['year']
zone = input_parameters['zone']
deepthmin = input_parameters['depth_min']
deepthmax = input_parameters['depth_max']


outputDirectory = os.getenv('AC_OUTPUT_DIR')
# TODO remove finalizing / in outputDirectory

wantedData = ['Temperature', 'Nitrate', 'Ammonium', 'eastward_Water_current', 'northward_Water_current']

dateBeginning = f'{year}-01-01 00:00:00'
dateEnd = f'{year + 1}-01-01 00:00:00'
frequency = 2

dataFin = pd.read_csv('./dataCmd.csv',';')
datesList = giveDateslist(dateBeginning, dateEnd, frequency)

# print('dataFin')
# pprint(dataFin)
# print('datesList')
# pprint(datesList)
for dat in wantedData:
    dataOutputDirectory = outputDirectory + '/' + dat + '/'
    dataLine = dataFin.loc[dataFin["Parameter"] == dat]
    print(f'dataLine {dataOutputDirectory}')
    pprint(dataLine)
    if dataLine.iloc[0]["daily"] > 0:
        getData(dat, zone, dataFin, deepthmin, deepthmax,  dataOutputDirectory, datesList[0], datesList[1],frequency)
    else:
        getData(dat, zone, dataFin, deepthmin, deepthmax, dataOutputDirectory)
