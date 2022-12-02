from code.utils.cmd import run
from code.utils.tiff import read, read_tiff, printBand, createRaster
import numpy as np
import math
import os

def optimal_farming(espece: str, scenario: str, prod: float, depth: float, surf: float, dist: float, minprod: float, ficbathy: str, ficbathyout: str, ficzee: str, ficin: str, ficout: str, outDir, mask:str):
    '''
    espece(str): the species studied
    scenario(str): the scenario considered
    prod(str): wanted production in MT
    depth(flaot): maximum depth in m
    surf(flaot): the farm surface in km2
    dist(float): the minimum distance between two farms in km
    minprod(float): the minimum production for a farm to be considered in kg/m2/year
    ficbathy(str): bathymetry file address
    ficbathyout(str): bathymetry output file address
    ficzee(str): european exclusive economic zone file address
    ficin(str): fresh weight in kg/m2/year PUA file address
    ficout(str): output file address
    outDir(str): output directory
    mask(str): the mask address
    '''
    print('optimal_farming')

    # ; EEZ file
    # ; -----------
    tiffImage = read_tiff(ficzee)
    dataset_zee, band_zee, array_zee = read(ficzee)
    nb_col = dataset_zee.RasterXSize
    nb_lig = dataset_zee.RasterYSize
    lon = [(x * tiffImage.pasx + tiffImage.lonmin) for x in range(nb_col)]
    lat = [(x * tiffImage.pasy + tiffImage.latmin) for x in range(nb_lig)]


    # ; Bathymetry extraction on the european eez
    # ; ----------------------------------------------
    cmd = 'gdalwarp -overwrite {ficbathy} -tr {pasx} {pasy} -dstnodata "9999" -r cubicspline -te {lonmin} {latmin} {lonmax} {latmax} {bathyout} -t_srs EPSG:4326'
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

    dataset_bath, band_bath, array_bath = read(ficbathyout)

    # we put at 0 the masked pixels
    array_bath[np.where(~np.isfinite(array_bath))] = 0

    # ; File reading
    # ; -----------------------------
    fictmp = '/tmp/' + os.path.basename(ficin)
    run('gdal_translate -co compress=none {ficin} {fictmp}'.format(ficin=ficin, fictmp=fictmp))
    dataset_prod, band_prod, array_prod = read(fictmp)
    run('rm -f ' + fictmp)
    # we put at 0 the masked pixels
    array_prod[np.where(array_prod > 100)] = 0 #we put to 0 the unrealistic productions
    array_prod[np.where(~np.isfinite(array_prod))] = 0
    array_prod = array_prod * 1000  # convertion between kg/m2 and T/km2

    array_prod[np.where((array_bath > 0) + (array_bath < depth))] = 0 #we mask the value with a too hight or too low bathimetry
    array_prod[np.where(array_zee <= 0)] = 0 # we mask the data out of the zee area

    if mask: # if we use a mask
        cmd = 'gdalwarp -overwrite {mask} -tr {pasx} {pasy} -dstnodata "-999" -r near -te {lonmin} {latmin} {lonmax} {latmax} {maskout} -t_srs EPSG:4326'
        run(cmd.format(
            mask=mask,
            pasx=tiffImage.pasx,
            pasy=tiffImage.pasy,
            lonmin=tiffImage.lonmin,
            latmin=tiffImage.latmin,
            lonmax=tiffImage.lonmax,
            latmax=tiffImage.latmax,
            maskout=outDir + '/reshaped_mask.tif'
        )) # we put the mask on the same grid as the others data

        dataset_mask, band_mask, array_mask = read(outDir + '/reshaped_mask.tif')
        array_prod[np.where(array_mask <= 0)] = 0


    ncol, nlig = np.shape(array_prod)

    farm = np.zeros((ncol, nlig)) # creating the array containing the selected farms production
    farmprod = np.zeros((ncol, nlig)) # creating the array containing the selected farms production
    prodfin = 0
    nbfarms = 0

    ulfic = open(ficout + '_posfarms.txt', 'w')

    wprod = prod * 1e6  # convertion from T/km2 to T
    _, max_prod = band_prod.ComputeRasterMinMax()
    lineNbr = 0
    print('minprod')
    print(minprod)
    ulfic.write('farm;lat;lon;fw \n')
    while (prodfin < wprod) and (lineNbr < 2e4):
        lineNbr += 1
        ind = np.where(array_prod == np.max(array_prod))  # we search the maximal production coordinates
        indi = ind[0][0]
        indj = ind[1][0]
        fprod = array_prod[indi, indj]  # we get the maximal production value
        print('fprod')
        print(fprod)
        if fprod < minprod:  # if the production is lower than the minimal production we stop searching after farms
            break
        farm[indi, indj] = 1  # we place the farm on the optimal farms map
        farmprod[indi, indj] = surf * array_prod[
            indi, indj]  # we place the farm production on the optimal farms production map
        prodfin += farmprod[indi, indj]
        flat = lat[::-1][indi]
        flon = lon[indj]

        i1 = max(indi - dist, 0)
        j1 = max(indj - dist, 0)
        i2 = min(indi + dist, ncol - 1)
        j2 = min(indj + dist, nlig - 1)
        array_prod[i1:i2,j1:j2] = 0  # we put to 0 the farms arround the selected farm to make sure to don't select farms to close to each others
        nbfarms = nbfarms + 1
        ulfic.write(';'.join([str(item) for item in [nbfarms, flat, flon,
                                                     fprod]]) + '\n')  # we save the selected farm coordinates and production value
    ulfic.close()
    print('nbr of farms:')
    print(str(nbfarms))
    print('final production: ')
    print(prodfin)

    # ; Writing the results
    # ; -----------------------
    createRaster(
        farm,
        ficout + '.tif',
        band_prod.DataType,
        dataset_prod.GetGeoTransform(),
        band_prod.GetNoDataValue()
    )
    createRaster(
        farmprod,
        ficout + '_prod.tif',
        band_prod.DataType,
        dataset_prod.GetGeoTransform(),
        band_prod.GetNoDataValue()
    )

    if np.sum((farmprod > 0) * 1) == 0: # if we have not selected any farm
        print('The minimal production is too hight')
    else: # otherwise we create a file which sum up the farm repartition
        minprod = np.min(farmprod[farmprod > 0])
        maxprod = np.max(farmprod)
        with open(ficout + '.txt', 'w') as f:
            f.write('readme')
            f.write('Species: ' + espece + '\n')
            f.write('Scenario: ' + scenario + '\n')
            f.write('Production to be reached: ' + str(prod) + 'MT/year' + '\n')
            f.write('Maximal depth: ' + str(depth) + 'm' + '\n')
            f.write('Surface of each farm: ' + str(surf) + 'km²' + '\n')
            f.write('Distance between farms: ' + str(dist) + 'km' + '\n')
            f.write('Number of farms: ' + str(nbfarms) + '\n')
            f.write('Minimal production of a farm: ' + str(minprod) + 'kg/year' + '\n')
            f.write('Maximal production of a farm: ' + str(maxprod) + 'kg/year' + '\n')
