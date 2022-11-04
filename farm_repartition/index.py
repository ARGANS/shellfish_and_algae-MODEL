import os

from code.optimal_farming import optimal_farming
from code.utils.tiff import read_tiff, fusionDesFichiers, fusion_zones2
from code.models.TiffImage import TiffImage
from code.utils.json import import_json, dump_json

def parse_parameters(path):
    try:
        model_parameters:dict = import_json(path)
    except OSError as e:
        raise RuntimeError(f'File {path} does not exists')

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

    return {
    	'especes': especes,
        'scenario': scenario,
        'zones': [zone],
        'params': params,
        'farm_distribution': {
            '_mask_file_path_exists': os.path.exists(mask_file_path) and os.path.isfile(mask_file_path),
            # TODO merge JSONs
            "production_required":  farm_distr.get("production_required", 10),
            "maximum_bathymetry": farm_distr.get("maximum_bathymetry", -30),
            "farms_separated_by": farm_distr.get("farms_separated_by", 10),
            "surface": farm_distr.get("surface", 1),
            "minimum_production": farm_distr.get("minimum_production", 10000),
            "external_mask_file": mask_file_path
        },
        'DATA_DIR': DATA_DIR,
        'OUT_DIR': OUT_DIR
    }


DATA_DIR=os.getenv('INPUT_SOURCE')
MAPS_DIR='/media/global/maps'
OUT_DIR=os.getenv('INPUT_DESTINATION')
TMP_DIR='/tmp'

model_parameters_path = DATA_DIR + '/parameters.json'
conf = parse_parameters(model_parameters_path)
print(f'Used variables of the current execution {conf}')
dump_json(conf, OUT_DIR + '/conf.json')

# ==========================================================================
#  Scénario A ou B
# ==========================================================================
#  Zones => 1km sur grille ZEE
#  ---------------------------

especes = conf['especes']
scenario = conf['scenario']
zones = conf['zones']
params = conf['params']

ficzee = MAPS_DIR + '/zee_europe.tif' # Fichier des zee
ficzee_image:TiffImage = read_tiff(ficzee)


#; Farms optimal distribution
#; -------------------------------

prods = [conf.get('farm_distribution', {}).get("production_required")] # Mega Tons à produire
depth = conf.get('farm_distribution', {}).get("maximum_bathymetry") #; Profondeur maximum à considérer en m
surf = conf.get('farm_distribution', {}).get("surface") #; surface de la ferme en km²
dist = conf.get('farm_distribution', {}).get("farms_separated_by") #; pas de fermes à moins 2 pixels (km)
prodsmin = [conf.get('farm_distribution', {}).get("minimum_production")] #; Poduction minimmum d'une ferme en kg/m²/an

# TODO mask file

ficbathy = MAPS_DIR + '/Bathy.TIF' #; Fichier bathymetie globale
ficzee = MAPS_DIR + '/zee_europe.tif' #; Fichier des zee
ficbathyout = OUT_DIR + '/bathy_europe.tif'
if zones != "Europe":
    ficin=DATA_DIR + '/concat/params_1km/FW_PUA_1km.tif'
else:
    ficin=DATA_DIR + '/FW_PUA_1km.tif'

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
            ficin=DATA_DIR + '/concat/params_1km/FW_PUA_1km.tif',
            ficout=OUT_DIR + f'/ZEE_{espece}_sc{scenario}farms-for-{prod}MT-FW_surf{surf}_depth{depth}_dist{dist}',
        )
