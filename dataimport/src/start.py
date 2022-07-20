import os
import json
from general import processCollectionOfProperties


def cleanFinalSlash(value: str) -> str:
    return value[:-1] if value.endswith('/') else value

VAR_PROP = 'variables'


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


if not VAR_PROP in input_parameters:
    raise RuntimeError(f'{VAR_PROP} not found')

processCollectionOfProperties(input_parameters[VAR_PROP], outputDirectory, year, deepthmin, deepthmax)

