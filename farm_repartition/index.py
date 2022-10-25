from code.zones_2_zee import zones_2_zee
from code.fusion_zones import fusion_zones
from code.optimal_farming import optimal_farming
from code.create_cdf import create_cdf
from code.utils.tiff import read_tiff, fusionDesFichiers, fusion_zones2
from code.models.TiffImage import TiffImage
from code.utils.cmd import run

DATA_DIR='/media/share/data/66/_dataread/b6f5a33f4e7ab4808f936f02c0a67c92/posttreatment' # TODO from env variable
MAPS_DIR='/media/global/maps'
OUT_DIR='/media/share/out' # TODO from env variables
TMP_DIR='/tmp'



# ==========================================================================
#  Scénario A ou B
# ==========================================================================
#  Zones => 1km sur grille ZEE
#  ---------------------------

# especes=['alaria','ulva','saccharina']
especes=['saccharina'] # TODO from json
scenario='B' # TODO from json
zones=["IBI"] # TODO from json
params=['DW_PUA','FW_PUA','protein_PUA','kcal_PUA','Biomass_CO2','CO2_uptake_PUA']
#zones=["Arctic", "Baltic", "NWS", "IBI", "MED", "BS"] # zones à traite

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

# especes=['alaria', 'ulva']
prods = [1.,2.,5.,10.] #	    		; Mega Tons à produire
depth = -30. # 			    ; Profondeur maximum à considérer en m
surf = 1. #   			    ; surface de la ferme en km² 
dist = 3. #     			; pas de fermes à moins 2 pixels (km)
prodsmin = [10.,3.,1.] #    		    ; Poduction minimmum d'une ferme en kg/m²/an
ficbathy = MAPS_DIR + '/Bathy.TIF' #  	; Fichier bathymetie globale
ficzee = MAPS_DIR + '/zee_europe.tif' #	; Fichier des zee 
ficbathyout = MAPS_DIR + '/bathy_europe.tif'


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
