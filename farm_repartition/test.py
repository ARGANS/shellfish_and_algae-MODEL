from code.zones_2_zee import zones_2_zee
from code.utils.tiff import read_tiff, fusionDesFichiers
from code.models.TiffImage import TiffImage
from code.utils.cmd import run

test_path1 = './tmp/IBI_alaria_Biomass_CO2_B_1km.tif'
test_path2 = './tmp/NWS_alaria_Biomass_CO2_B_1km.tif'
fusionDesFichiers(test_path1, test_path2)
