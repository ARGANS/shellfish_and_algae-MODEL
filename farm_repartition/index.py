import os
from code.zones_2_zee import zones_2_zee
from code.fusion_zones import fusion_zones
from code.optimal_farming import optimal_farming
from code.create_cdf import create_cdf
from code.utils.tiff import read_tiff, fusionDesFichiers, fusion_zones2
from code.models.TiffImage import TiffImage
from code.utils.cmd import run
from code.utils.json import import_json, dump_json


DATA_DIR=os.getenv('INPUT_SOURCE', '/media/share/data/66/_dataread/b6f5a33f4e7ab4808f936f02c0a67c92/posttreatment')
MAPS_DIR='/media/global/maps'
OUT_DIR=os.getenv('INPUT_DESTINATION', '/media/share/out')
TMP_DIR='/tmp'


model_parameters_path = DATA_DIR + '/parameters.json'

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


for espece in especes:
    for param in params:
        for zone in zones:
            ficin = DATA_DIR + f'/{param}.tif'
            ficout = OUT_DIR + f'/out1/{zone}_{espece}_{param}_{scenario}_1km.tif'
            fictmp = TMP_DIR + f'/{zone}_{espece}_{param}_{scenario}_1km.tif'
            zones_2_zee(ficzee=ficzee_image, ficin=ficin, ficout=ficout, fictmp=fictmp,  scenario=scenario)

"""
# DISABLED
# Fusion de 2  zones 
# ------------------
zones=['NWS','IBI'] # 2 zones à fusionner
params=['DW_PUA','FW_PUA','protein_PUA','kcal_PUA','Biomass_CO2','CO2_uptake_PUA']

for espece in especes:
    for param in params:
        ficins = []
        fictmps = []
        for zone in zones:
            ficin = OUT_DIR + '/out1' + f'/{zone}_{espece}_{param}_{scenario}_1km.tif'
            fictmp = TMP_DIR + f'/{zone}_{espece}_{param}_{scenario}_1km.tif'
            ficins.append(ficin)
            fictmps.append(fictmp)
            cmd = f'gdal_translate -co compress=none {ficin} {fictmp}'
            run(cmd)

        fictmp=f'./tmp/tmp_{espece}_{param}_{scenario}_1km.tif'
        fusionDesFichiers(*fictmps, fictmp) 

        # ; Ecriture du résultat
        zone_range = '-'.join(zones)
        ficout=OUT_DIR + f'/out2/{zone_range}_{espece}_{param}_{scenario}_1km.tif'
        cmd=f'gdal_translate -co compress=deflate {fictmp} {ficout}'
        run(cmd)
        # run('rm -f ' + fictmp)
"""

# ; Fusion des zones => EU
# ; ----------------------
# zones=['ARC','BKS','BAL' ,'NWS-IBI','MED'] #	 ;zones à fusionner
# zones=['Baltic','BS','NWS','IBI','MED','Artic']
# especes=['alaria'] #,'saccharina','ulva']
# params=['DW_PUA','FW_PUA','protein_PUA','kcal_PUA','Biomass_CO2','CO2_uptake_PUA']
ficzee_europe = MAPS_DIR + '/zee_europe.tif'	#; Fichier des zee 
ficzee_europe_image:TiffImage = read_tiff(ficzee_europe)

for espece in especes:
    for param in params:
        fichiers_a_fusioner = ' '.join([f'{OUT_DIR}/out1/{zone}_{espece}_{param}_{scenario}_1km.tif' for zone in zones]) 
        ficout = OUT_DIR + f'/out3/EU_{espece}_{param}_{scenario}_1km.tif'
        fictmp = TMP_DIR + f'/EU_{espece}_{param}_{scenario}_1km.tif'
        fictmp_zee=TMP_DIR + f'/ZEE_{espece}_{param}_{scenario}_1km.tif'
        ficout_zee=OUT_DIR + f'/out3/ZEE_{espece}_{param}_{scenario}_1km.tif'
        fusion_zones(ficzee_europe_image, fichiers_a_fusioner, fictmp, ficout)
        fusion_zones2(ficzee_europe, fictmp, fictmp_zee)
        run('rm -f ' + fictmp)
        run(f'gdal_translate -co compress=deflate {fictmp_zee} {ficout_zee}')
        run('rm -f ' + fictmp_zee)


"""
# DISABLED
# ; Création du fichier netcdf
# ; --------------------------


# especes=['saccharina','alaria','ulva']
especes=['alaria']
libespeces=['Saccharina latissima','Alaria esculenta','Ulva lactuca']
params=['DW_PUA','FW_PUA','protein_PUA','kcal_PUA','Biomass_CO2','CO2_uptake_PUA']
libparams=['Dry Weight (kg.m-2.an-1)','Fresh Weight (kg.m-2.an-1)','Energy (kcal.m-2.an-1)', 'CO2 (kcal.m-2.an-1)','Protein production (kg.m-2.an-1)','Equivalent emissions CO2 (kg.m-2.an-1']
unites=['kg.m-2.an-1','kg.m-2.an-1','kg.m-2.an-1', 'kcal.m-2.an-1', 'kg.m-2.an-1', 'kg.m-2.an-1']
spatial_resolution='1km'
libzone='European EEZ'
source='CMEMS models, EMODNET, NASA/OCENACOLOR'


for espece, libespece in zip(especes, libespeces):
    for param, libparam, unite in zip(params, libparams, unites):
        create_cdf(
            ficin=OUT_DIR + f'/out3/ZEE_{espece}_{param}_{scenario}_1km.tif',
            fictmp=f'./tmp/ZEE_{espece}_{param}_{scenario}_1km.tif',
            ficout=OUT_DIR + f'/out3/ZEE_{espece}_{param}_{scenario}_1km.nc',
            espece=espece, 
            libespece=libespece,
            nomparam=param, 
            libparam=libparam,
            scenario=scenario, 
            unite=unite,
            libzone=libzone, 
            annee=2021, 
            title='Maps of '+ libparam,
			spatial_resolution=spatial_resolution, 
            source=source
        )
"""
#; Répartition optimale des fermes
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
            ficin=OUT_DIR + f'/out3/ZEE_{espece}_FW_PUA_{scenario}_1km.tif',
            ficout=OUT_DIR + f'/out3/ZEE_{espece}_sc{scenario}farms-for-{prod}MT-FW_surf{surf}_depth{depth}_dist{dist}',
        )
