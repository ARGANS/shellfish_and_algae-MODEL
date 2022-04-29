import copy
import os
import time
import netCDF4 as nc
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import numpy.ma as ma
from scipy.sparse import dia_matrix

# extract the data value at depth in the merged files (all the daily data merged in one file)
from advectionPrototype.climatology_ellipses import degrees_to_meters
from dataread.read_netcdf import extractVarSubset

def getwantedMergeData(data, depth, dataCmdpath, mergedFilepath='D:/Profils/mjaouen/Documents/alternance/EASME/data/'):
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

#give the indices coresponding to lonval, and latval in the list of coordinates
def givecoor(ds,lonval,latval,dataName,dataFin):
    DataLine = dataFin.loc[dataFin["Parameter"] == dataName]
    # we get the longitude and latitude list
    lonList = ds[DataLine.iloc[0]["longName"]][:]
    latList = ds[DataLine.iloc[0]["latName"]][:]
    i=0
    loni = lonList[i]
    # we browse the data until we find a coordiate bigger than the wanted coordiante
    while i+1<len(lonList) and lonval>loni:
        i+=1
        loni = lonList[i]
    j = 0
    lati = latList[j]
    while j+1 < len(latList) and latval > lati:
        j += 1
        lati = latList[j]
    return i, j

def u2d_cgrid_cur(var):
    varcg = np.zeros((var.shape[0], var.shape[1], var.shape[2] + 1))
    for i in range(var.shape[1]):
        for j in range(var.shape[2] - 2,-1,-1):
            if ((not ma.is_masked(var[0,i, j])) & (not ma.is_masked(var[0,i, j + 1]))):
                #varcg[:, i , j+ 1] = 2 * var[:, i , j+ 1] - varcg[:, i, j + 2]
                varcg[:,i, j ] = (var[:,i, j] + var[:,i, j+1])/2
                '''            elif ((not ma.is_masked(var[0,i, j+1])) & (ma.is_masked(var[0,i, j + 2]))):
                                varcg[:, i, j + 1] = var[:, i, j + 1]'''
    return varcg

def v2d_cgrid_cur(var):
    varcg = np.zeros((var.shape[0],var.shape[1] + 1, var.shape[2]))
    for i in range(var.shape[1] - 2,-1,-1):
        for j in range(var.shape[2]):
            if ((not ma.is_masked(var[0,i, j])) & (not ma.is_masked(var[0,i + 1, j]))):
                varcg[:, i, j] =  (var[:, i, j] + var[:, i + 1, j])/2
                #varcg[:,i + 1, j] = 2 * var[:,i+1, j] - varcg[:,i+2, j]
    return varcg

def giveCFL(dx,dy,dt,Ewc,Nwc,nbrx,nbry):
    Cx = np.abs(Ewc[:,1:]*dt/dx)
    Cy = np.abs(Nwc[1:]*dt/dy)
    return np.maximum(Cx,Cy).reshape(nbrx * nbry), Cx.reshape(nbrx * nbry), Cy.reshape(nbrx * nbry)

def createListNan(maskPosition,nbry):
    listeNan = []
    for k in range(len(maskPosition[0])):
        i, j = maskPosition[0][k], maskPosition[1][k]
        n = j + nbry * i
        listeNan += [n]
    return listeNan

def findNan(dataNO3,nbry):
    mask = dataNO3.mask
    maskpos2D = np.where(mask[0] == True)
    dataNO3[0][maskpos2D] = np.nan
    westNanij_1 = np.where(np.isnan(dataNO3[0][:, :-1]))
    westNan = (westNanij_1[0], westNanij_1[1] + 1) #we have a Nan at the west

    eastNanij_1 = np.where(np.isnan(dataNO3[0][:, 1:]))
    eastNan = (eastNanij_1[0], eastNanij_1[1]) #we have a Nan at the east

    downNani_1j = np.where(np.isnan(dataNO3[0][:-1]))
    downNan = (downNani_1j[0] + 1, downNani_1j[1])

    upNani_1j = np.where(np.isnan(dataNO3[0][1:]))
    upNan = (upNani_1j[0], upNani_1j[1])

    up2Nani_1j = np.where(np.isnan(dataNO3[0][2:]))
    up2Nan = (up2Nani_1j[0], up2Nani_1j[1])

    down2Nani_1j = np.where(np.isnan(dataNO3[0][:-2]))
    down2Nan = (down2Nani_1j[0] + 2, down2Nani_1j[1])

    east2Nani_1j = np.where(np.isnan(dataNO3[0][:,2:]))
    east2Nan = (east2Nani_1j[0], east2Nani_1j[1])

    west2Nani_1j = np.where(np.isnan(dataNO3[0][:,:-2]))
    west2Nan = (west2Nani_1j[0], west2Nani_1j[1]+2)

    upNanList = createListNan(upNan, nbry)
    downNanList = createListNan(downNan, nbry)
    eastNanList = createListNan(eastNan, nbry)
    westNanList = createListNan(westNan, nbry)

    up2NanList = createListNan(up2Nan, nbry)
    down2NanList = createListNan(down2Nan, nbry)
    east2NanList = createListNan(east2Nan, nbry)
    west2NanList = createListNan(west2Nan, nbry)
    return westNanList, eastNanList, upNanList, downNanList, west2NanList, east2NanList, up2NanList, down2NanList

def createMatE(nbrx, nbry,centuredEwc,eastNanList, westNanList, alpha1, alpha2):
    offset = np.array([0, -1, 1, -2, 2])
    uGreater0 = ((centuredEwc[:, 1:] > 0) * 1).reshape(nbrx * nbry)
    termA = (1 - 2 * alpha1 + alpha2) * uGreater0 - (1 - uGreater0) * (2 * alpha1 + alpha2 - 1)
    termB = (alpha1 - 2 * alpha2 - 1) * uGreater0 - alpha1 * (1 - uGreater0)
    termC = alpha1 * uGreater0 + (1 - alpha1 - 2 * alpha2) * (1 - uGreater0)
    termD = uGreater0 * alpha2
    termE = alpha2 * (1 - uGreater0)

    termB[::nbry] = 0
    termC[nbry - 1::nbry] = 0

    termA[westNanList] = 0
    termB[westNanList] = 0

    termA[eastNanList] = 0
    termC[eastNanList] = 0
    data = np.zeros((5, nbrx * nbry))
    data[0] = termA
    data[1][:-1] = termB[1:]
    data[2][1:] = termC[:-1]
    data[3][:-2] = termD[2:]
    data[4][2:] = termE[:-2]
    Mat = dia_matrix((data, offset), shape=(nbrx * nbry, nbrx * nbry))

    return Mat


def createMatN(nbrx, nbry, centuredNwc, downNanList, upNanList, alpha1, alpha2):
    offset = np.array([0, -nbry, nbry, -2*nbry, 2*nbry])
    vGreater0 = ((centuredNwc[1:] > 0) * 1).reshape(nbrx * nbry)
    termA = (1-2*alpha1+alpha2) * vGreater0 - (1 - vGreater0) * (2*alpha1+alpha2-1)
    termB = (alpha1-2*alpha2-1) * vGreater0 - alpha1 * (1 - vGreater0)
    termC = alpha1 * vGreater0 + (1-alpha1-2*alpha2) * (1 - vGreater0)
    termD = vGreater0*alpha2
    termE = alpha2*(1 - vGreater0)

    termA[downNanList] = 0
    termB[downNanList] = 0

    termA[upNanList] = 0
    termC[upNanList] = 0
    data = np.zeros((5, nbrx * nbry))
    data[0] = termA
    data[1][:-nbry] = termB[nbry:]
    data[2][nbry:] = termC[:-nbry]
    data[3][:-2*nbry] = termD[2*nbry:]
    data[4][2*nbry:] = termE[:-2*nbry]
    Mat = dia_matrix((data, offset), shape=(nbrx * nbry, nbrx * nbry))

    return Mat

def quickest(dyMeter, dxlist, dt,discr, centuredEwc, centuredNwc, dataNO3):
    nbrStep = len(centuredNwc) * discr
    C = np.zeros(np.shape(dataNO3[0]))
    (nbrx, nbry) = np.shape(dataNO3[0])
    mask = dataNO3.mask
    maskpos2D = np.where(mask[0] == True)
    init = time.time()
    westNanList, eastNanList, upNanList, downNanList, west2NanList, east2NanList, up2NanList, down2NanList = findNan(dataNO3,nbry)
    for k in range(nbrStep - 1):
        dayNbr =  k// (discr)
        hNbr = k//int(discr/24)
        C[maskpos2D] = 0
        CFL, cx, cy = giveCFL(dyMeter, dxlist, dt, centuredEwc[hNbr], centuredNwc[hNbr],nbrx, nbry)
        #print('CFL max: ', max(CFL))
        alpha1 = (1 / 6) * (1 - CFL) * (2 - CFL)
        alpha2 = (1 / 6) * (1 - CFL) * (1 + CFL)
        matE = createMatE(nbrx, nbry, centuredEwc[hNbr], eastNanList, westNanList, alpha1, alpha2)
        matN = createMatN(nbrx, nbry, centuredNwc[hNbr], downNanList, upNanList, alpha1, alpha2)

        #B = np.zeros((nbrx , nbry))
        B = -dt*0.1*(C +dataNO3[dayNbr]).reshape(nbrx * nbry)
        #B = -B.reshape(nbrx * nbry)

        Cline = C.reshape(nbrx * nbry) - cy*matN.dot(C.reshape(nbrx * nbry)) - cx*matE.dot(C.reshape(nbrx * nbry)) + B
        C = Cline.reshape(nbrx, nbry)
        if k % int(discr) == 0:
            print(time.time() - init)
            print(k / discr)
            fig, ax = plt.subplots()
            CpnIm = copy.deepcopy(C)
            CpnIm[maskpos2D] = np.nan
            plt.imshow(CpnIm)
            plt.colorbar()
            ax.invert_yaxis()
            ax.set_title('c')
            plt.show()
            fig, ax = plt.subplots()
            CpnIm = copy.deepcopy(C +dataNO3[dayNbr])
            CpnIm[maskpos2D] = np.nan
            plt.imshow(CpnIm)
            plt.colorbar()
            ax.invert_yaxis()
            ax.set_title('C prim')
            plt.show()

if __name__ == "__main__":

    zone = 'IBI'
    depth = 0

    latmin = 43
    latmax = 48

    lonmin = -4
    lonmax = -0.5

    dxdeg = 0.028
    dydeg = 0.028

    discr = 50 * 80
    dt = 1 / discr

    dataCmdpath = 'I:/work-he/apps/safi/data/IBI/dataCmd.csv'
    dataNO3, ds = getwantedMergeData('Nitrate', depth, dataCmdpath)

    fileU = 'D:/Profils/mjaouen/Documents/alternance/EASME/data/cmems_mod_ibi_phy_u0.nc'
    dataBaseEwc = nc.Dataset(fileU)

    fileV = 'D:/Profils/mjaouen/Documents/alternance/EASME/data/cmems_mod_ibi_phy_v0.nc'
    dataBaseNwc = nc.Dataset(fileV)

    dataFin = pd.read_csv('./../dataimport/src/dataCmd.csv', ';')

    ewcDataLine = dataFin.loc[(dataFin["Parameter"] == 'eastward_Water_current') & (dataFin["Place"] == zone)]
    ewcdataName = ewcDataLine.iloc[-1]["variable"]  # we find the data name in the dataset

    nwcDataLine = dataFin.loc[(dataFin["Parameter"] == 'northward_Water_current') & (dataFin["Place"] == zone)]
    nwcdataName = nwcDataLine.iloc[-1]["variable"]  # we find the data name in the dataset

    longitudeMin, latitudeMin = givecoor(dataBaseEwc, lonmin, latmin, 'eastward_Water_current', dataFin)  # we get the indices of the wanted position
    longitudeMax, latitudeMax = givecoor(dataBaseEwc, lonmax, latmax, 'eastward_Water_current', dataFin)  # we get the indices of the wanted position

    dataNO3 = dataNO3[:, latitudeMin:latitudeMax, longitudeMin:longitudeMax]
    dataNwc = dataBaseNwc[nwcdataName][:, latitudeMin:latitudeMax, longitudeMin:longitudeMax]
    dataEwc = dataBaseEwc[ewcdataName][:, latitudeMin:latitudeMax, longitudeMin:longitudeMax]

    dataNwc *= 60 * 60 * 24
    dataEwc *= 60 * 60 * 24

    latRef = np.ones((np.shape(dataEwc[0])[1],np.shape(dataEwc[0])[0]))*np.array(dataBaseEwc[ewcDataLine.iloc[-1]["latName"]][latitudeMin:latitudeMax])

    centuredEwc = u2d_cgrid_cur(dataEwc)
    centuredNwc = v2d_cgrid_cur(dataNwc)

    dxlist, dyMeter = degrees_to_meters(dxdeg, dydeg, latRef)

    dxlist = dxlist.T
    dylist = dyMeter*np.ones(np.shape(dxlist))

    quickest(dyMeter, dxlist, dt, discr, centuredEwc, centuredNwc, dataNO3)

    fig, ax = plt.subplots()
    plt.imshow(dataNwc[1])
    plt.colorbar()
    ax.invert_yaxis()
    # plt.clim(-10, 10)
    ax.set_title('Northward current')

    fig1, ax1 = plt.subplots()
    plt.imshow(centuredNwc[1])
    plt.colorbar()
    ax1.invert_yaxis()
    # plt.clim(-10, 10)
    ax1.set_title('centured Northward current')
    plt.show()

    fig, ax = plt.subplots()
    plt.imshow(dataEwc[1])
    plt.colorbar()
    ax.invert_yaxis()
    # plt.clim(-10, 10)
    ax.set_title('Eastward current')

    fig1, ax1 = plt.subplots()
    plt.imshow(centuredEwc[1])
    plt.colorbar()
    ax1.invert_yaxis()
    # plt.clim(-10, 10)
    ax1.set_title('centured Eastward current')
    plt.show()