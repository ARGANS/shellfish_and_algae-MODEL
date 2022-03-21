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

#compute the nutriment at day plus one, taking into account a consumption of deficitPct% of the day before nutriment
def calcDeficit(nitrateAverageDayPlus1, nitrateAverage, northStream, eastStream,deficitPct):
    nitrateAverageDayPlus1def = copy.deepcopy(nitrateAverageDayPlus1)
    for i in range(len(northStream)):
        for j in range(len(northStream[0])):
            nitrateAverageDayPlus1def[i+int(northStream[i,j]),int(j+eastStream[i,j])] -= nitrateAverage[i,j]*deficitPct
    return nitrateAverageDayPlus1def

#compute the nutriment at day plus one, taking into account a consumption of deficitPct% of the day before nutriment
def calcDeficit(nitrateAverageDayPlus1, nitrateAverage, meanEast,meanNorth, stdPar, stdPer, angle, deficitPct):
    nitrateAverageDayPlus1def = copy.deepcopy(nitrateAverageDayPlus1)
    yCoordMat = np.ones((len(nitrateAverage),len(nitrateAverage[0])))*np.arange(len(nitrateAverage[0]))
    xcoordMat = (np.ones((len(nitrateAverage[0]),len(nitrateAverage)))*np.arange(len(nitrateAverage))).T
    for i in range(len(nitrateAverage)):
        for j in range(len(nitrateAverage[0])):
            xB2, yB2 = giveCoorNewBasis(angle[i,j], xcoordMat, yCoordMat, i+meanEast, j+meanNorth)
            inEllipsemat = inEllipse(stdPar, stdPer, xB2, yB2)
            # we substract the deficit to the cells in the ellipse
            nitrateAverageDayPlus1def[inEllipsemat] -= nitrateAverage[i,j]*deficitPct*(1/len(inEllipsemat[0]))
    return nitrateAverageDayPlus1def

#give the coordiantes of a point in the basis centured at x1, y1 with angleb between the new and the original basis axis
def giveCoorNewBasis(angleb, x, y, x1, y1):
    cos = np.cos(angleb)
    sin = np.sin(angleb)
    distx = x-x1
    disty = y-y1
    return distx*cos+disty*sin, -distx*sin+disty*cos

def inEllipse(a,b,x,y):
    r = (x/a)**2+(y/b)**2
    return r<=1

if __name__ == "__main__":
    #we get the average on the 10 first meters
    nwc = 'northward_Water_current'
    ewc = 'eastward_Water_current'
    nut = 'Nitrate'
    depth = 3
    dataCmdpath = 'I:/work-he/apps/safi/data/IBI/dataCmd.csv'
    dataNwc, _ = getwantedMergeData(nwc, depth, dataCmdpath)
    dataEwc, _ = getwantedMergeData(ewc, depth, dataCmdpath)
    dataNut, ds = getwantedMergeData(nut, depth, dataCmdpath)
    Nmax = len(dataEwc)
    mask = dataEwc.mask
    maskposition = np.where(mask == True) # we get the position of the masked data

    meanEast = np.mean(dataEwc, axis=0)
    meanNorth = np.mean(dataNwc, axis=0)

    angleCos = meanEast/np.sqrt(meanEast**2+meanNorth**2)
    angle = np.arccos(angleCos)

    stdE = np.std(dataEwc, axis=0)
    stdN = np.std(dataNwc, axis=0)
    stdPar = np.cos(angle)*stdE - np.sin(angle)*stdN #stdPar parallel to the mean vector
    stdPer = np.sin(angle)*stdE + np.cos(angle)*stdN #perpendicular to the mean vector

    maskposition = np.where(mask == True) # we get the position of the masked data
    generalList = np.zeros((Nmax, len(dataNwc[0]), len(dataNwc[0][0])))
    i = 0
    while i < Nmax - 1:
        print(i)
        # we compute the nutriment at day D+1, with a deficit of 10% of day D nutriment
        nitrateAverageDayPlus1deficit = calcDeficit(dataNut[i + 1], dataNut[i], meanEast, meanNorth, stdPar, stdPer, angle, 0.2)
        generalList[i] = nitrateAverageDayPlus1deficit
        i += 1
    generalList[maskposition] = np.nan
    P10nut = np.percentile(generalList, 10, axis=0)  # we compute the 10th percentile

    xsize, ysize, ulx, uly, xres, yres = getMetadata(ds)

    saveAsTiff(P10nut, xsize, ysize, ulx, uly, xres, yres,"I:/work-he/apps/safi/data/IBI/P10with20pctEllipseDeficit3m.tiff")

    fig1, ax1 = plt.subplots()
    plt.imshow(P10nut)
    ax1.invert_yaxis()
    ax1.set_title('P10 NO3 IBI 2020')
    plt.clim(0, 8)
    plt.colorbar()
    plt.show()