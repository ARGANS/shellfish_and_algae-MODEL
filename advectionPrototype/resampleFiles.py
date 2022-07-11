import os

import pandas as pd

import netCDF4 as nc

from advectionPrototype.saveAsTiff import getMetadata
from dataread.read_netcdf import extractVarSubset


def getwantedMergeData(data, depth, dataCmdpath, zone, mergedFilepath='D:/Profils/mjaouen/Documents/alternance/EASME/data/MED/'):
    csvFile = pd.read_csv(dataCmdpath, ';')
    DataLine = csvFile.loc[(csvFile["Parameter"] == data) & (csvFile["type"] == 'model') & (csvFile["Place"] == zone)]
    variable = DataLine.iloc[0]["variable"]
    # we search after the file containing the wanted data
    for r, d, f in os.walk(mergedFilepath):
        for i in range(len(f)):
            filename = f[i]
            # if it is the wanted file
            if filename[:10] == 'merged_' + data[:3]:
                fn = mergedFilepath + f[i]
                ncDataset = nc.Dataset(fn)
                break
    return extractVarSubset(ncDataset, variable, depth=depth)[0], ncDataset

if __name__ == "__main__":
    zone = 'Arctic'
    variableName = 'Temperature'
    depth = 0
    dataCmdpath = 'D:/Profils/mjaouen/Documents/alternance/EASME/gitcode/shellfish_and_algae-MODEL/dataimport/src/dataCmd.csv'
    mergedFilepath = 'D:/Profils/mjaouen/Documents/alternance/EASME/data/{zone}/'.format(zone=zone)

    _, ds = getwantedMergeData('Nitrate', depth, dataCmdpath, zone, mergedFilepath)
    dataFin = pd.read_csv('./../dataimport/src/dataCmd.csv', ';')
    ewcDataLine = dataFin.loc[(dataFin["Parameter"] == 'Nitrate') & (dataFin["Place"] == zone)]
    ewcdataName = ewcDataLine.iloc[-1]["variable"]  # we find the data name in the dataset
    xsize, ysize, ulx, uly, xres, yres = getMetadata(ds, ewcDataLine.iloc[-1]["latName"],
                                                     ewcDataLine.iloc[-1]["longName"])
    longMini = ulx
    latiMax = uly
    longMax = ulx + xsize * xres
    latiMini = uly - ysize * yres
    variable = (dataFin.loc[(dataFin["Parameter"] == variableName) & (dataFin["Place"] == zone)]).iloc[-1]["variable"]
    ficInitial = f"D:/Profils/mjaouen/Documents/alternance/EASME/data/Arctic/merged_{variableName}_ArcticStere.nc".format(variableName = variableName)
    #ficInitial = '"D:/Profils/mjaouen/Documents/alternance/EASME/data/reference_monthly_PAR_2021_filled_4_max0.nc"'
    ficin = "I:/work-he/apps/safi/data/{zone}/{variableName}_Stereo.tiff".format(zone= zone, variableName = variableName)
    ficout = "I:/work-he/apps/safi/data/{zone}/{variableName}_EPSG.tiff".format(zone= zone, variableName = variableName)
    ngtiffFile = "I:/work-he/apps/safi/data/{zone}/{variableName}_NewGrid.tiff".format(zone= zone, variableName = variableName)
    ficFinal = "I:/work-he/apps/safi/data/{zone}/{variableName}_NewGrid.nc".format(zone= zone, variableName = variableName)

    fillval = str(nc.Dataset(ficInitial).variables[variable]._FillValue)

    print(fillval, type(fillval))

    os.system("gdal_translate -a_nodata "+'"'+fillval+'"'+" NETCDF:"+ ficInitial +":"+variable+" " + ficin)
    #os.system("gdal_translate -of GTiff -a_nodata "+'"-999"'+" "+ficInitial +" "+ ficin)

    '''os.system('gdalwarp -s_srs ' + ficInitial + ' -tr ' + str(xres) + ' ' + str(
        yres) + ' -r cubicspline -srcnodata ' + '"-999"' + ' -te ' + str(longMini) + ' ' + str(latiMini) + ' ' + str(
        longMax) + ' ' + str(latiMax) + ' ' + ficout + ' -t_srs EPSG:4326')'''
    # To get the data in EPSG:4326 .tiff from Stereographic .nc
    os.system('gdalwarp -s_srs "+proj=stere +lat_0=90 +lat_ts=90 +lon_0=-45 +k=0.00001 +x_0=0 +y_0=0 +a=6378273 +b=6378273  +units=m +no_defs" '+ficin+' -srcnodata ' +'"'+fillval+'"'+ ' -t_srs EPSG:4326  "'+ ficout+'" -overwrite')
    os.system('gdalwarp '+ ficout+' -tr '+str(xres)+' '+str(yres)+' -r cubicspline -srcnodata '+'"'+fillval+'"'+' -te '+str(longMini)+' '+str(latiMini)+' '+str(longMax)+' '+str(latiMax)+' '+ngtiffFile+' -t_srs EPSG:4326')

    os.system("gdal_translate -a_nodata "+'"'+fillval+'"'+" -of NetCDF "+ ngtiffFile + " " + ficFinal)