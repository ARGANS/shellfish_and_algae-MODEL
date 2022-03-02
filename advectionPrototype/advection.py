import time
import datetime
import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import pandas as pd

#import sys
#sys.path.append('../../dataImport/dataread/')
from dataread.read_netcdf import extractWithAverage

mainpath = 'I:/work-he/apps/safi/data/IBI/'
datanorth = 'northward_Water_current'
date = datetime.datetime(2020,4,16)
dateplus1 = date + datetime.timedelta(days=1)

# give the average along depth between the two depth written in the depth tuple
def getAverageAlongDepth(data:str, depth:tuple, date:datetime.datetime, mainpath:str, dataCmdpath = 'I:/work-he/apps/safi/data/IBI/dataCmd.csv' ):
    path = mainpath + data + '/'
    csvFile=pd.read_csv(dataCmdpath,';')
    DataLine = csvFile.loc[csvFile["Parameter"] == data]
    variable = DataLine.iloc[0]["variable"]
    #we search after the file containing the wanted data at the right date
    for r, d, f in os.walk(path):
        for i in range(len(f)):
            filename = f[i]
            #we read the date corresponding to the file
            stringDate = filename[-13:-3]
            fileDate = datetime.datetime(int(stringDate[0:4]), int(stringDate[5:7]), int(stringDate[8:10]))
            #if it is the wanted date
            if date == fileDate:
                fn=path+f[i]
                ncDataset = nc.Dataset(fn)
                break
    return extractWithAverage(ncDataset, variable, ('depth',), depth=depth, time_index=0)[0]

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

#compute the nutriment at day plus one, taking into account a consumption of 10% of the day before nutriment
def calcDeficit(nitrateAverageDayPlus1, nitrateAverage, northStream, eastStream):
    for i in range(len(northStream)):
        for j in range(len(northStream[0])):
            nitrateAverageDayPlus1[i+int(northStream[i,j]),int(j+eastStream[i,j])] -= nitrateAverage[i,j]*0.1
    return nitrateAverageDayPlus1

if __name__ == "__main__":
    initialTime = time.time()
    #we get the average on the 10 first meters
    northcurrentAverage = getAverageAlongDepth('northward_Water_current', (0,10), date, mainpath)
    eastcurrentAverage = getAverageAlongDepth('eastward_Water_current', (0,10), date, mainpath)
    mask = eastcurrentAverage.mask
    maskposition = np.where(mask == True)#we get the position of the masked data
    #we compute the average of nutriment at day D and day D+1
    nitrateAverage = getAverageAlongDepth('Nitrate', (0,10), date, mainpath)
    nitrateAverageDayPlus1 = getAverageAlongDepth('Nitrate', (0, 10), dateplus1, mainpath)
    northStream, eastStream = giveStream(northcurrentAverage, eastcurrentAverage)
    nbrUpStrm = giveNbrOfUpstreamCell(northStream, eastStream)
    finaltime = time.time()-initialTime
    print('running time : ',finaltime,' seconds')
    # we compute the nutriment at day D+1, with a deficit of 10% of day D nutriment
    nitrateAverageDayPlus1 = calcDeficit(nitrateAverageDayPlus1, nitrateAverage, northStream, eastStream)

    northStream[maskposition] = np.nan
    eastStream[maskposition] = np.nan
    nbrUpStrm[maskposition] = np.nan
    nitrateAverageDayPlus1[maskposition] = np.nan
    fig, ax = plt.subplots()
    im = ax.imshow((np.abs(northStream) + np.abs(eastStream) - 1) * (np.abs(northStream) + np.abs(eastStream) - 2) / 2)
    ax.invert_yaxis()
    ax.set_title('grid cells that points to themself')
    fig.colorbar(im)
    plt.show()

    fig, ax = plt.subplots()
    plt.imshow(nitrateAverage)
    ax.set_title('nitrate Average')
    ax.invert_yaxis()
    plt.clim(0, 8)
    plt.colorbar()
    plt.show()

    fig, ax = plt.subplots()
    plt.imshow(nitrateAverageDayPlus1)
    ax.set_title('nitrate Average the day after')
    ax.invert_yaxis()
    plt.clim(0, 8)
    plt.colorbar()
    plt.show()

    fig, ax = plt.subplots()
    plt.imshow((nitrateAverageDayPlus1<0)*1.)
    ax.set_title('nitrate Average the day after, neagtives values')
    ax.invert_yaxis()
    plt.colorbar()
    plt.show()

    fig, ax = plt.subplots()
    im = ax.imshow(nbrUpStrm)
    ax.invert_yaxis()
    ax.set_title('Number of upstream cells')
    fig.colorbar(im)
    plt.show()

    anyUpcell = (nbrUpStrm==0)*1.
    anyUpcell[maskposition] = np.nan
    fig, ax = plt.subplots()
    im = ax.imshow(anyUpcell)
    ax.invert_yaxis()
    ax.set_title('0 upstream cells')
    fig.colorbar(im)
    plt.show()

    morethan1 = (nbrUpStrm > 1) * 1.
    morethan1[maskposition] = np.nan
    fig, ax = plt.subplots()
    im = ax.imshow(morethan1)
    ax.invert_yaxis()
    ax.set_title('Multiple upstream cells')
    fig.colorbar(im)
    plt.show()

    fig, ax = plt.subplots()
    im = ax.imshow(eastStream)
    ax.set_title('east stream')
    ax.invert_yaxis()
    fig.colorbar(im)
    plt.show()

    fig, ax = plt.subplots()
    im = ax.imshow(northStream)
    ax.set_title('north stream')
    ax.invert_yaxis()
    fig.colorbar(im)
    plt.show()