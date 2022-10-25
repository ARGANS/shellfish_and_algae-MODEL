from code.utils.cmd import run
from code.utils.tiff import read, read_tiff, printBand, createRaster
import numpy as np
import math
import os    

def optimal_farming(espece, scenario, prod, depth, surf, dist, minprod, ficbathy, ficbathyout, ficzee, ficin, ficout):
    print('optimal_farming')

    # ; Fichier ZEE
    # ; -----------
    tiffImage = read_tiff(ficzee)
    dataset_zee, band_zee, array_zee = read(ficzee)
    nb_col=dataset_zee.RasterXSize
    nb_lig=dataset_zee.RasterYSize
    lon = [(x*tiffImage.pasx + tiffImage.lonmin + tiffImage.pasx / 2) for x in range(nb_col)]
    lat = [(x*tiffImage.pasy + tiffImage.latmin + tiffImage.pasy / 2) for x in range(nb_lig)]
    
    #------------------------------------------------------------------------------------
    # zee=rotate(read_tiff(ficzee, geo=geo), 7)
    # sz=size(zee, /dim)
    # nbcol=sz(0)
    # nblig=sz(1)
    # pasx=geo.MODELPIXELSCALETAG(0)
    # pasy=geo.MODELPIXELSCALETAG(1)
    # lonmin=geo.MODELTIEPOINTTAG(3)
    # latmax=geo.MODELTIEPOINTTAG(4)
    # lonmax=lonmin+(nbcol*pasx)
    # latmin=latmax-(nblig*pasy)
    # lats=findgen(5031)*pasy+latmin+pasy
    # lons=findgen(8005)*pasx+lonmin+pasx
    #------------------------------------------------------------------------------------

    # ; Extraction de la bathymétrie sur la grille zee
    # ; ----------------------------------------------
    cmd = 'gdalwarp -overwrite {ficbathy} -tr {pasx} {pasy} -r cubicspline -te {lonmin} {latmin} {lonmax} {latmax} {bathyout} -t_srs EPSG:4326'
    run(cmd.format(
        ficbathy=ficbathy,
        pasx=tiffImage.pasx,
        pasy=tiffImage.pasy,
        lonmin=tiffImage.lonmin,
        latmin=tiffImage.latmin,
        lonmax=tiffImage.lonmax,
        latmax=tiffImage.latmax,
        bathyout=ficbathyout
    ))

    # bath=rotate(read_tiff(bathyout, geo=geob), 7)
    dataset_bath, band_bath, array_bath = read(ficbathyout)

    # bath[where(~finite(bath), /null)] = 0
    for y in range(band_bath.YSize):
        for x in range(band_bath.XSize):
            if not math.isfinite(array_bath[y][x]):
                array_bath[y][x] = 0

    # ; Lecture du fichier production
    # ; -----------------------------
    fictmp = './tmp/' + os.path.basename(ficin)
    run('gdal_translate -co compress=none {ficin} {fictmp}'.format(ficin=ficin, fictmp=fictmp))
    dataset_prod, band_prod, array_prod = read(fictmp)
    run('rm -f ' + fictmp)

    for y in range(band_prod.YSize):
        for x in range(band_prod.XSize):
            if array_prod[y, x] > 1000:
                array_prod[y, x] = 0
            elif not math.isfinite(array_prod[y][x]):
                array_prod[y, x] = 0
            else:
                array_prod[y, x] *= 1000
            # ind=where(r GT 1000,ct)
            # IF (ct GT 0) THEN r(ind)=0
            # r[where(~finite(r), /null)] = 0
            # r=r*1000. 						; pour passer en T/km²

            # ; Masquage par Bathy et Zee        
            if array_bath[y, x] <= depth or array_zee[y, x] <=0:
                array_prod[y, x] = 0
            # ind=where(bath LE depth, ct)
            # IF (ct GT 0) THEN r(ind)=0
            # ind=where(zee LE 0, ct)
            # IF (ct GT 0) THEN r(ind)=0

    ncol=dataset_prod.RasterXSize
    nlig=dataset_prod.RasterYSize

    # ; Boucle sur les pixels
    # ; ---------------------
    farm =  (np.arange(ncol * nlig) * 0).reshape(ncol, nlig)
    # farm=bytarr(ncol, nlig)
    farmprod =  (np.arange(ncol * nlig) * 0).reshape(ncol, nlig)
    # farmprod=fltarr(ncol, nlig)
    prodfin=0
    nbfarms=0

    ulfic = open(ficout+'_posfarms.txt', 'w')

    wprod=prod*1000000 # ; pour passer en T
    _, max_prod = band_prod.ComputeRasterMinMax()
    while prodfin < wprod:
        for y in range(band_prod.YSize):
            for x in range(band_prod.XSize):
                if array_prod[y, x] == max_prod and y > 0 and x > 0:
                    farm[y, x] = 1
                    farmprod[y, x] = surf * array_prod[y, x]
                    prodfin += farmprod[y, x]

                fprod = array_prod[y, x]
                flat = lat[y]
                flon = lon[x]

                x1=int(x-dist)
                if not (x1 > 0):
                    x1 = 0 
                y1=int(y-dist) 
                if not (y1 > 0):
                    y1 = 0
                x2=int(x+dist)
                if not (x2 < ncol):
                    x2 = ncol
                y2=int(y+dist)
                if not (y2 < nlig):
                    y2=nlig
                array_prod[x1:x2, y1:y2] = 0
                nbfarms=nbfarms+1
                ulfic.write(';'.join([str(item) for item in [nbfarms, flat, flon, fprod]]) + '\n')
    ulfic.close()

    # ; Ecriture des résultats
    # ; -----------------------
    createRaster(
        farm, 
        ficout+'.tif',
        band_prod.DataType, 
        dataset_prod.GetGeoTransform(), 
        band_prod.GetNoDataValue()
    )
    createRaster(
        farmprod, 
        ficout+'_prod.tif',
        band_prod.DataType, 
        dataset_prod.GetGeoTransform(), 
        band_prod.GetNoDataValue()
    )

    minprod=np.min(farmprod[farmprod > 0])
    maxprod=np.max(farmprod)
    with open(ficout + '.txt', 'w') as f:
        f.write('readme')
        f.write('Species: '+ espece + '\n')
        f.write('Scenario: ' + scenario + '\n')
        f.write('Production to be reached: ' + prod + 'MT/year' + '\n')
        f.write('Maximal depth: ' + depth + 'm' + '\n')
        f.write('Surface of each farm: ' + surf + 'km²' + '\n')
        f.write('Distance between farms: ' + dist + 'km' + '\n')
        f.write('Number of farms: ' + nbfarms + '\n')
        f.write('Minimal production of a farm: ' + minprod + 'kg/year' + '\n')
        f.write('Maximal production of a farm: ' + maxprod + 'kg/year' + '\n')
