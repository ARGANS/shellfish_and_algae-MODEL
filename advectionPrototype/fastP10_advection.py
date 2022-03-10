import copy
import time
import datetime
import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import pandas as pd
from saveAsTiff import getMetadata, saveAsTiff

from dataread.read_netcdf import extractVarSubset

#mainpath = 'I:/work-he/apps/safi/data/IBI/'
mainpath = '/mount/internal/work-he/apps/safi/data/IBI'
datanorth = 'northward_Water_current'
date = datetime.datetime(2020,4,16)

#extract the data value at depth in the merged files (all the daily data merged in one file)
def getwantedMergeData(data,depth,dataCmdpath,mergedFilepath = 'I:/work-he/apps/safi/data/IBI/merged_files/'):
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

# this function put 0 to the stream that goes outside the studied area
def correctBorder(streamMat):
    issueBorderUp = np.where(streamMat[0]==-1)
    streamMat[0][issueBorderUp]=0
    issueBorderDown = np.where(streamMat[len(streamMat)-1] == 1)
    streamMat[len(streamMat)-1][issueBorderDown] = 0
    return streamMat

# give the stream depending of the northward and eastward current
def giveStream(northcurrentAverage, eastcurrentAverage):
    #we compute the cos and the sin of the angle between the eastward vector and the stream
    norm2 = northcurrentAverage**2 + eastcurrentAverage**2
    sin = northcurrentAverage/np.sqrt(norm2)
    cos = eastcurrentAverage/np.sqrt(norm2)

    #depending of the value of the cos and sin, we indentify the direction taken by the stream
    eastStream = np.zeros(np.shape(eastcurrentAverage))
    eastStream += (cos>0.5)*1
    eastStream += (cos < -0.5) * -1

    northStream = np.zeros(np.shape(eastcurrentAverage))
    northStream += (sin > 0.5) * 1
    northStream += (sin < -0.5) * -1

    #we correct the stream direction to avoid stream to go out of the studied area
    northStream = correctBorder(northStream)
    eastStream = correctBorder(eastStream.T).T
    return northStream, eastStream

#compute the number of upstream cell for each cell
def giveNbrOfUpstreamCell(northStream, eastStream):
    nbrUpStrm=np.zeros(np.shape(northStream))
    for i in range(len(northStream)):
        for j in range(len(northStream[0])):
            nbrUpStrm[i+int(northStream[i,j]),int(j+eastStream[i,j])] +=1
    return nbrUpStrm

#compute the nutriment at day plus one, taking into account a consumption of deficitPct% of the day before nutriment
def calcDeficit(nitrateAverageDayPlus1, nitrateAverage, northStream, eastStream,deficitPct):
    nitrateAverageDayPlus1def = copy.deepcopy(nitrateAverageDayPlus1)
    for i in range(len(northStream)):
        for j in range(len(northStream[0])):
            nitrateAverageDayPlus1def[i+int(northStream[i,j]),int(j+eastStream[i,j])] -= nitrateAverage[i,j]*deficitPct
    return nitrateAverageDayPlus1def

if __name__ == "__main__":
    #we get the average on the 10 first meters
    nwc = 'northward_Water_current'
    ewc = 'eastward_Water_current'
    nut = 'Nitrate'
    depth = 3
    dataCmdpath = 'I:/work-he/apps/safi/data/IBI/dataCmd.csv'
    averageDataNwc, _ = getwantedMergeData(nwc, depth, dataCmdpath)
    averageDataEwc, _ = getwantedMergeData(ewc, depth, dataCmdpath)
    averageDataNut, ds = getwantedMergeData(nut, depth, dataCmdpath)
    Nmax = len(averageDataEwc)
    mask = averageDataNut.mask
    maskposition = np.where(mask == True)  # we get the position of the masked data
    generalList = np.zeros((Nmax,len(averageDataNut[0]),len(averageDataNut[0][0])))
    i = 0
    while i<Nmax-1:
        print(i)
        northStream, eastStream = giveStream(averageDataNwc[i], averageDataEwc[i])
        # we compute the nutriment at day D+1, with a deficit of 10% of day D nutriment
        nitrateAverageDayPlus1deficit = calcDeficit(averageDataNut[i+1], averageDataNut[i], northStream, eastStream, 0.2)
        generalList[i] = nitrateAverageDayPlus1deficit
        i+=1
    generalList[maskposition] = np.nan
    P10nut = np.percentile(generalList, 10, axis=0)#we compute the 10th percentile
    xsize, ysize, ulx, uly, xres, yres = getMetadata(ds)

    saveAsTiff(P10nut, xsize, ysize, ulx, uly, xres, yres,"I:/work-he/apps/safi/data/IBI/P10with20pctDeficit3m.tiff")

    fig1, ax1 = plt.subplots()
    plt.imshow(P10nut)
    ax1.invert_yaxis()
    ax1.set_title('P10 NO3 IBI 2020')
    plt.clim(0, 8)
    plt.colorbar()
    plt.show()