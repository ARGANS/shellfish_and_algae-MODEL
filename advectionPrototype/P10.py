import time
import datetime
import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import pandas as pd
from dataread.read_netcdf import extractWithAverage, extractVarSubset
from saveAsTiff import getMetadata, saveAsTiff

data = 'Nitrate'
mergedFilePath = 'I:/work-he/apps/safi/data/IBI/merged_files/'
dataCmdpath = 'I:/work-he/apps/safi/data/IBI/dataCmd.csv'
depth = 3


#extract the data value at depth in the merged files (all the daily data merged in one file)
def getwantedMergeData(data,depth,dataCmdpath,mergedFilepath = 'D:/Profils/mjaouen/Documents/alternance/EASME/data/'):
    csvFile = pd.read_csv(dataCmdpath, ';')
    DataLine = csvFile.loc[csvFile["Parameter"] == data]
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
    initialTime = time.time()
    averageData, ds = getwantedMergeData(data,depth,dataCmdpath)
    print(np.shape(averageData))
    mask = averageData.mask
    maskposition = np.where(mask == True)  # we get the position of the masked data
    averageData[maskposition] = np.nan
    P10nut = np.percentile(averageData, 10, axis=0)#we compute the 10th percentile
    xsize, ysize, ulx, uly, xres, yres = getMetadata(ds)

    finaltime = time.time() - initialTime
    print('running time : ', finaltime, ' seconds')
    saveAsTiff(P10nut, xsize, ysize, ulx, uly, xres, yres,
               "I:/work-he/apps/safi/data/IBI/P102020.tiff")
    fig1, ax1 = plt.subplots()
    plt.imshow(P10nut)
    ax1.invert_yaxis()
    ax1.set_title('P10 NO3 IBI 2020')
    plt.clim(0, 8)
    plt.colorbar()

    fig, ax2 = plt.subplots()
    plt.imshow(averageData[0])
    ax2.invert_yaxis()
    ax2.set_title('NO3 IBI 18/01/2020')
    plt.clim(0, 8)
    plt.colorbar()
    plt.show()