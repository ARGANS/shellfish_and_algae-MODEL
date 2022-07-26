import os
import json
from general import processCollectionOfProperties


def cleanFinalSlash(value: str) -> str:
    return value[:-1] if value.endswith('/') else value

MTDT_KEY = 'metadata'
DTSTS_KEY = 'datasets'


# Read dataset parameters
parameters_json_value:str = os.getenv('INPUT_PARAMETERS') 
try:
    input_parameters:dict = json.loads(parameters_json_value)
except:
    raise RuntimeError('Cannot parse the value of the INPUT_PARAMETERS environment variable')

outputDirectory = cleanFinalSlash(os.getenv('INPUT_DESTINATION'))		

year = int(input_parameters['year'])
deepthmin = int(input_parameters['depth-min'])
deepthmax = int(input_parameters['depth-max'])


datasets = input_parameters.get(DTSTS_KEY, None)

if datasets is None:
    raise RuntimeError(f'{DTSTS_KEY} not found')

processCollectionOfProperties(datasets, outputDirectory, year, deepthmin, deepthmax)

