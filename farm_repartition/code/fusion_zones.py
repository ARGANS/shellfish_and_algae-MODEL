from code.utils.cmd import run
from code.models.TiffImage import TiffImage

COMMAND1 = 'gdalwarp -overwrite -srcnodata "-999. nan" {lst} {tmp} -tr {pasx} {pasy} -r near -te {lonmin} {latmin} {lonmax} {latmax} -dstnodata "-999."  -t_srs EPSG:4326'

def fusion_zones(ficzee:TiffImage, fichiers_a_fusioner:str, fictmp, ficout):
    cmd1 = COMMAND1.format(
        lst=fichiers_a_fusioner, 
        tmp=fictmp,
        pasx=ficzee.pasx,
        pasy=ficzee.pasy,
        lonmin=ficzee.lonmin,
        latmin=ficzee.latmin,
        lonmax=ficzee.lonmax,
        latmax=ficzee.latmax
    )
    cmd2 = f'gdal_translate -co compress=deflate {fictmp} {ficout}'
    run(cmd1)
    run(cmd2)
    # fw=rotate(read_tiff(fictmp, geo=geo), 7)
    
