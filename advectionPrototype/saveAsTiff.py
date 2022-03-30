from osgeo import gdal
from osgeo import osr
import numpy as np
import netCDF4 as nc
import os
import pandas as pd
import datetime

# from advectionPrototype.advection import getAverageAlongDepth
from dataread.read_netcdf import extractWithAverage

#return ds metadatas
def getMetadata(ds):
    lats = ds['latitude'][:]
    lons = ds['longitude'][:]

    xres = lons[1] - lons[0]
    yres = lats[1] - lats[0]

    ysize = len(lats)
    xsize = len(lons)

    ulx = lons[0]
    uly = lats[-1]

    return xsize, ysize, ulx, uly, xres, yres

#save dataArray in a geotiff file
def saveAsTiff(dataArray,xsize, ysize, ulx, uly, xres, yres, filepath ):
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(filepath, xsize, ysize, 1, gdal.GDT_Float32)

    # this assumes the projection is Geographic lat/lon WGS 84
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    ds.SetProjection(srs.ExportToWkt())

    gt = [ulx, xres, 0, uly, 0, -yres]
    ds.SetGeoTransform(gt)

    outband = ds.GetRasterBand(1)
    outband.WriteArray(np.flip(dataArray, 0))

    del ds


if __name__ == "__main__":
    dataCmdpath = 'I:/work-he/apps/safi/data/IBI/dataCmd.csv'
    mainpath = 'I:/work-he/apps/safi/data/IBI/'
    datanorth = 'northward_Water_current'
    data = 'Nitrate'

    path = mainpath + data + '/'
    csvFile = pd.read_csv(dataCmdpath, ';')
    DataLine = csvFile.loc[csvFile["Parameter"] == data]
    variable = DataLine.iloc[0]["variable"]
    date = datetime.datetime(2020, 4, 16)

    for r, d, f in os.walk(path):
        for i in range(len(f)):
            filename = f[i]
            # we read the date corresponding to the file
            stringDate = filename[-13:-3]
            fileDate = datetime.datetime(int(stringDate[0:4]), int(stringDate[5:7]), int(stringDate[8:10]))
            # if it is the wanted date
            if date == fileDate:
                fn = path + f[i]
                ds = nc.Dataset(fn)

    #nitrateAverage = getAverageAlongDepth('Nitrate', (0, 10), date, mainpath)
    #xsize, ysize, ulx, uly, xres, yres = getMetadata(ds)
    #saveAsTiff(nitrateAverage, xsize, ysize, ulx, uly, xres, yres, "I:/work-he/apps/safi/data/IBI/NO3Average16042020.tiff")

