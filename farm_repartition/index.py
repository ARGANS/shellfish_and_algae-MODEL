import os

from code.optimal_farming import optimal_farming
from code.utils.tiff import read_tiff
from code.models.TiffImage import TiffImage
from code.utils.json import import_json, dump_json, loads

def parse_parameters(model_parameters: dict):
    '''
    Create a dictonary with the needed information for farm distribution
    This information come from the model_parameters dictioanry
    '''
    print(f'Model parameters: {model_parameters}')

    type = model_parameters.get('type', 'Algae')

    if type == 'Algae':
        params = ['DW_PUA','FW_PUA','protein_PUA','kcal_PUA','Biomass_CO2','CO2_uptake_PUA']
    elif type == 'Shellfish':
        params = ["NH4_production", "CO2_production", "DSTW", "SHL", "STE", "FW", "DWW"]

    *especes, = model_parameters.\
        get('parameters', {}).\
        get('species', {})
    zone = model_parameters.\
        get('metadata', {}).\
        get('zone', None)
        
    scenario = model_parameters.\
        get('metadata', {}).\
        get('scenario', None)

    farm_distr = model_parameters.\
        get('parameters', {}).\
        get('OptimalFarmsRepartition', {}).\
        get('default', {}).\
        get('parameters', {})
    mask_file_path = '/'.join(DATA_DIR.split('/')[:-3]) + '/' + farm_distr.get('external_mask_file', '')
    if type == 'Algae':
        return {
    	'especes': especes,
        'scenario': scenario,
        'zones': [zone],
        'params': params,
        'farm_distribution': {
            "production_required":  farm_distr.get("production_required", 10),
            "maximum_bathymetry": farm_distr.get("maximum_bathymetry", -30),
            "farms_separated_by": farm_distr.get("farms_separated_by", 10),
            "mask_name": farm_distr.get("external_mask_file", 'maskFile.tif'),
            "surface": farm_distr.get("surface", 1),
            "minimum_production": farm_distr.get("minimum_production", 1000),
            "external_mask_file": mask_file_path
        },
        'DATA_DIR': DATA_DIR,
        'OUT_DIR': OUT_DIR
        }, type
    elif type == 'Shellfish':
        return {
    	'especes': especes,
        'scenario': scenario,
        'zones': [zone],
        'params': params,
        'farm_distribution': {
            "production_required":  farm_distr.get("production_required", 10),
            "maximum_bathymetry": farm_distr.get("maximum_bathymetry", -30),
            "farms_separated_by": farm_distr.get("farms_separated_by", 10),
            "mask_name": farm_distr.get("external_mask_file", 'maskFile.tif'),
            "surface": farm_distr.get("surface", 1),
            "minimum_production": farm_distr.get("minimum_production", 1000),
            "external_mask_file": mask_file_path
        },
        'DATA_DIR': DATA_DIR,
        'OUT_DIR': OUT_DIR
        }, type
    else:
        return {}, type

DATA_DIR=os.getenv('INPUT_SOURCE')
MAPS_DIR='/media/global/maps'
OUT_DIR=os.getenv('INPUT_DESTINATION')
TMP_DIR='/tmp'

serialized_properties = os.getenv('INPUT_MODEL_PROPERTIES_JSON')

if serialized_properties is not None:
    model_parameters:dict = loads(serialized_properties)
    print('INPUT_MODEL_PROPERTIES_JSON')
    print(f'INPUT_MODEL_PROPERTIES_JSON: {serialized_properties} // {serialized_properties.__class__}')
else:
    # we read the json file
    model_parameters_path = DATA_DIR + '/parameters.json'
    try: # we try to read the file
        model_parameters:dict = import_json(model_parameters_path)
    except OSError as e:
        raise RuntimeError(f'File {model_parameters_path} does not exists: {e}')
    
conf, type = parse_parameters(model_parameters)
print(f'Used variables of the current execution {conf}')
dump_json(conf, OUT_DIR + '/conf.json')

especes = conf['especes']
scenario = conf['scenario']
zones = conf['zones']
params = conf['params']

ficzee = MAPS_DIR + '/zee_europe.tif' # eez file
ficzee_image:TiffImage = read_tiff(ficzee)

mask_name = conf.get('farm_distribution', {}).get("mask_name") #look if a mask is defined
if os.path.exists(DATA_DIR + '/../'+mask_name) and len(mask_name)>1:
    mask = DATA_DIR + '/../'+mask_name
else:
    mask = None

#; Farms optimal distribution
#; -------------------------------

prods = [conf.get('farm_distribution', {}).get("production_required")] # Mega Tons to be produced
depth = conf.get('farm_distribution', {}).get("maximum_bathymetry") #; maximum depth in m
surf = conf.get('farm_distribution', {}).get("surface") #; farm surface in km²
dist = conf.get('farm_distribution', {}).get("farms_separated_by") #; minimal distance between two farms (km)
prodsmin = [conf.get('farm_distribution', {}).get("minimum_production")] #; minimal production for a farm in kg/m²/an

ficbathy = MAPS_DIR + '/Bathy.TIF' #; Bathymetry file
ficzee = MAPS_DIR + '/zee_europe.tif' #; EEZ mask file
ficbathyout = OUT_DIR + '/bathy_europe.tif' #output file for the bathymetry
#we search after the fresh weight per unit area file
if type == 'Algae':
    varname = 'FW_PUA'
elif type == 'Shellfish':
    varname = 'FW'
if zones[0] != "Europe":
    ficin=DATA_DIR + f'/concat/params_1km/{varname}_1km.tif'
else:
    ficin=DATA_DIR + f'/{varname}.tif'

#launch the optimisation
for espece, minprod in zip(especes, prodsmin):
    for prod in prods:
        optimal_farming(
            espece=espece,
            scenario=scenario,
            prod=prod,
            depth=depth,
            surf=surf,
            dist=dist,
            minprod=minprod,
            ficbathy=ficbathy,
            ficbathyout=ficbathyout,
            ficzee=ficzee,
            ficin=ficin,
            ficout=OUT_DIR + f'/opt',
            outDir = OUT_DIR,
            mask = mask
        )
