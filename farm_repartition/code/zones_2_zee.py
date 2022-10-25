from code.utils.cmd import run
from code.models.TiffImage import TiffImage

__all__ = ['zones_2_zee']

COMMAND = 'gdalwarp -overwrite -srcnodata {coef} {ficin} {fictmp} -tr {pasx} {pasy} -r cubicspline -te {lonmin} {latmin} {lonmax} {latmax} -dstnodata "-999.".  -t_srs EPSG:4326'

def zones_2_zee(ficzee:TiffImage, ficin:str, ficout:str, fictmp, scenario:str):
    cmds = [
        COMMAND.format(
            coef=('"9.96921e+36"' if scenario == 'A' else '"nan"'),
            ficin=ficin,
            fictmp=fictmp,
            pasx=ficzee.pasx,
            pasy=ficzee.pasy,
            lonmin=ficzee.lonmin,
            lonmax=ficzee.lonmax,
            latmin=ficzee.latmin,
            latmax=ficzee.latmax
        ),
        f'gdal_translate -co compress=deflate {fictmp} {ficout}',
        f'rm -f {fictmp}'
    ]

    for cmd in cmds:
        run(cmd)
