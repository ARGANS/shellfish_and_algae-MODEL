import copy
import os
import time
import netCDF4 as nc
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import numpy.ma as ma
from scipy import interpolate
from scipy.sparse import dia_matrix
import datetime
from dateutil.relativedelta import *

# extract the data value at depth in the merged files (all the daily data merged in one file)
from advectionPrototype.climatology_ellipses import degrees_to_meters
from advectionPrototype.saveAsTiff import saveAsTiff, getMetadata
from dataread.launch_model import MA_model_scipy
from dataread.make_runs import open_data_input
from dataread.read_netcdf import extractVarSubset
from dataread.utils import import_json


def getwantedMergeData(data, depth, dataCmdpath, mergedFilepath='D:/Profils/mjaouen/Documents/alternance/EASME/data/'):
    csvFile = pd.read_csv(dataCmdpath, ';')
    DataLine = csvFile.loc[(csvFile["Parameter"] == data) & (csvFile["type"] == 'model')]
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


# sort the list of data from the older to the newer
def sortDateList(listValue, ldate):
    sortLval = []  # we define the sorted data list
    sortldate = []  # we define the sorted date list
    # we read the list we have to sort
    for k in range(len(ldate)):
        i = 0
        # while we didn't red all the sorted list
        while i < len(sortLval):
            # if the list element is graeter than the element we want to place in the list
            if ldate[k] < sortldate[i]:
                # we place this element before the element who is greater than it
                sortldate = sortldate[:i] + [ldate[k]] + sortldate[i:]
                sortLval = sortLval[:i] + [listValue[k]] + sortLval[i:]
                i = len(sortLval) + 1
            i += 1
        if i == len(sortLval) - 1 or i == len(sortLval):
            sortldate += [ldate[k]]
            sortLval += [listValue[k]]
    return np.array(sortLval), np.array(sortldate)


def getData(path, data, zone, depth, latitudeMinwc, latitudeMaxwc, longitudeMinwc, longitudeMaxwc):
    dataFin = pd.read_csv('./../dataimport/src/dataCmd.csv', ';')
    wantedDataLine = dataFin.loc[(dataFin["Parameter"] == data) & (dataFin["Place"] == zone)]
    data = wantedDataLine.iloc[-1]["variable"]  # we find the data name in the dataset
    listValue = []
    for r, d, f in os.walk(path):
        listValue = np.zeros(((len(f) - 10) * 361, latitudeMaxwc - latitudeMinwc, longitudeMaxwc - longitudeMinwc))
        for i in range(len(f) - 10):
            fn = path + f[i]
            print(fn)
            # we read the file
            ds = nc.Dataset(fn)
            # we get the date
            # ldate += [ds['time'][0]]
            # we get the data
            listValue[i * 361:(i + 1) * 361] = ds[data][:, latitudeMinwc:latitudeMaxwc, longitudeMinwc:longitudeMaxwc]
    return listValue


# give the indices coresponding to lonval, and latval in the list of coordinates
def givecoor(ds, lonval, latval, dataName, dataFin):
    DataLine = dataFin.loc[dataFin["Parameter"] == dataName]
    # we get the longitude and latitude list
    lonList = ds[DataLine.iloc[0]["longName"]][:]
    latList = ds[DataLine.iloc[0]["latName"]][:]
    i = 0
    loni = lonList[i]
    # we browse the data until we find a coordiate bigger than the wanted coordiante
    if lonval == None:
        i= None
    else:
        while i + 1 < len(lonList) and lonval > loni:
            i += 1
            loni = lonList[i]
    j = 0
    lati = latList[j]
    if latval == None:
        j = None
    else:
        while j + 1 < len(latList) and latval > lati:
            j += 1
            lati = latList[j]
    return i, j


def u2d_cgrid_cur(var):
    varcg = np.zeros((var.shape[0], var.shape[1], var.shape[2] + 1))
    for i in range(var.shape[1]):
        for j in range(var.shape[2] - 2, -1, -1):
            if ((not ma.is_masked(var[0, i, j])) & (not ma.is_masked(var[0, i, j + 1]))):
                varcg[:, i, j] = (var[:, i, j] + var[:, i, j + 1]) / 2
    return varcg


def v2d_cgrid_cur(var):
    varcg = np.zeros((var.shape[0], var.shape[1] + 1, var.shape[2]))
    for i in range(var.shape[1] - 2, -1, -1):
        for j in range(var.shape[2]):
            if ((not ma.is_masked(var[0, i, j])) & (not ma.is_masked(var[0, i + 1, j]))):
                varcg[:, i, j] = (var[:, i, j] + var[:, i + 1, j]) / 2
    return varcg


def giveCFL(dx, dy, dt, Ewc, Nwc, nbrx, nbry):
    Cx = np.abs(Ewc[:, 1:] * dt / dx)
    Cy = np.abs(Nwc[1:] * dt / dy)
    return np.maximum(Cx, Cy).reshape(nbrx * nbry), Cx.reshape(nbrx * nbry), Cy.reshape(nbrx * nbry)


def createListNan(maskPosition, nbry):
    listeNan = []
    for k in range(len(maskPosition[0])):
        i, j = maskPosition[0][k], maskPosition[1][k]
        n = j + nbry * i
        listeNan += [n]
    return listeNan


def findNan(dataNO3, nbry):
    dataNO3Copy = copy.deepcopy(dataNO3)
    mask = dataNO3Copy.mask
    maskpos2D = np.where(mask[0] == True)
    dataNO3Copy[0][maskpos2D] = np.nan

    westNanij_1 = np.where(np.isnan(dataNO3Copy[0][:, :-1]))
    westNan = (westNanij_1[0], westNanij_1[1] + 1)  # we have a Nan at the west

    eastNan = np.where(np.isnan(dataNO3Copy[0][:, 1:]))  # we have a Nan at the east

    downNani_1j = np.where(np.isnan(dataNO3Copy[0][:-1]))
    downNan = (downNani_1j[0] + 1, downNani_1j[1])

    upNan = np.where(np.isnan(dataNO3Copy[0][1:]))

    up2Nan = np.where(np.isnan(dataNO3Copy[0][2:]))

    down2Nani_1j = np.where(np.isnan(dataNO3Copy[0][:-2]))
    down2Nan = (down2Nani_1j[0] + 2, down2Nani_1j[1])

    east2Nan = np.where(np.isnan(dataNO3Copy[0][:, 2:]))

    west2Nani_1j = np.where(np.isnan(dataNO3Copy[0][:, :-2]))
    west2Nan = (west2Nani_1j[0], west2Nani_1j[1] + 2)

    upEastNan = np.where(np.isnan(dataNO3Copy[0][1:, 1:]))

    downWestNani_1j_1 = np.where(np.isnan(dataNO3Copy[0][:-1, :-1]))
    downWestNan = (downWestNani_1j_1[0] + 1, downWestNani_1j_1[1] + 1)

    upWestNanj_1 = np.where(np.isnan(dataNO3Copy[0][1:, :-1]))
    upWestNan = (upWestNanj_1[0], upWestNanj_1[1] + 1)

    downEastNani_1 = np.where(np.isnan(dataNO3Copy[0][:-1, 1:]))
    downEastNan = (downEastNani_1[0] + 1, downEastNani_1[1])

    upNanList = createListNan(upNan, nbry)
    downNanList = createListNan(downNan, nbry)
    eastNanList = createListNan(eastNan, nbry)
    westNanList = createListNan(westNan, nbry)

    up2NanList = createListNan(up2Nan, nbry)
    down2NanList = createListNan(down2Nan, nbry)
    east2NanList = createListNan(east2Nan, nbry)
    west2NanList = createListNan(west2Nan, nbry)

    downEastNanList = createListNan(downEastNan, nbry)
    upWestNanList = createListNan(upWestNan, nbry)
    downWestNanList = createListNan(downWestNan, nbry)
    upEastNanList = createListNan(upEastNan, nbry)

    return westNanList, eastNanList, upNanList, downNanList, west2NanList, east2NanList, up2NanList, down2NanList, \
           downEastNanList, upWestNanList, downWestNanList, upEastNanList


def createMatE(nbrx, nbry, decenturedEwc, eastNanList, westNanList, east2NanList, west2NanList, alpha1, alpha2):
    offset = np.array([0, -1, 1, -2, 2])
    uGreater0 = ((decenturedEwc[:, 1:] > 0) * 1).reshape(nbrx * nbry)
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

    termA[west2NanList] = 0
    termB[west2NanList] = 0

    termA[east2NanList] = 0
    termC[east2NanList] = 0
    data = np.zeros((5, nbrx * nbry))
    data[0] = termA
    data[1][:-1] = termB[1:]
    data[2][1:] = termC[:-1]
    data[3][:-2] = termD[2:]
    data[4][2:] = termE[:-2]
    Mat = dia_matrix((data, offset), shape=(nbrx * nbry, nbrx * nbry))

    return Mat


def createMatN(nbrx, nbry, decenturedNwc, downNanList, upNanList, down2NanList, up2NanList, alpha1, alpha2):
    offset = np.array([0, -nbry, nbry, -2 * nbry, 2 * nbry])
    vGreater0 = ((decenturedNwc[1:] > 0) * 1).reshape(nbrx * nbry)
    termA = (1 - 2 * alpha1 + alpha2) * vGreater0 - (1 - vGreater0) * (2 * alpha1 + alpha2 - 1)
    termB = (alpha1 - 2 * alpha2 - 1) * vGreater0 - alpha1 * (1 - vGreater0)
    termC = alpha1 * vGreater0 + (1 - alpha1 - 2 * alpha2) * (1 - vGreater0)
    termD = vGreater0 * alpha2
    termE = alpha2 * (1 - vGreater0)

    termA[downNanList] = 0
    termB[downNanList] = 0

    termA[upNanList] = 0
    termC[upNanList] = 0

    termA[down2NanList] = 0
    termB[down2NanList] = 0

    termA[up2NanList] = 0
    termC[up2NanList] = 0
    data = np.zeros((5, nbrx * nbry))
    data[0] = termA
    data[1][:-nbry] = termB[nbry:]
    data[2][nbry:] = termC[:-nbry]
    data[3][:-2 * nbry] = termD[2 * nbry:]
    data[4][2 * nbry:] = termE[:-2 * nbry]
    Mat = dia_matrix((data, offset), shape=(nbrx * nbry, nbrx * nbry))

    return Mat


def giveEpsilon(day, temp, NH4, NO3,cNO3, cNH4, advNO3, advNH4, N_s, N_f, D, Nwc, Ewc, latRef, model, nbrx, nbry, dt):
    data_in = {
        'SST': temp,
        'PAR': 500,
        'NH4_ext': NH4,
        'NO3_ext': NO3,
        'PO4_ext': 50,
        'K_d': 0.1,
        'F_in': np.sqrt(Nwc ** 2 + Ewc ** 2),
        'h_z_SML': 30,
        't_z': 10,
        'D_ext': 0.1
    }
    y = dict(zip(["NH4", "NO3", "N_s", "N_f", "D"], [NH4, NO3, N_s, N_f, D]))
    return model.derivative_fast_advection(day, y, data_in,cNO3, cNH4, advNO3, advNH4, nbrx, nbry, dt,
                                           latRef, model)


def quickest(dyMeter, dxlist, dt, discr, Ewc, Nwc, centEwc, centNwc, latRef, dataNO3, dataNH4, dataTemp, Ks, firstday,
             model):
    nbrStep = int(len(Nwc) * discr)
    cNO3 = np.zeros(np.shape(dataNO3[0]))
    cNH4 = np.zeros(np.shape(dataNH4[0]))
    N_s = np.ones(np.shape(dataNH4[0])) * 1000
    N_f = np.ones(np.shape(dataNH4[0])) * 1000
    D = np.ones(np.shape(dataNH4[0])) * 0.1
    (nbrx, nbry) = np.shape(dataNO3[0])
    mask = dataNO3.mask
    maskpos2D = np.where(mask[0] == True)
    init = time.time()
    westNanList, eastNanList, upNanList, downNanList, west2NanList, east2NanList, up2NanList, down2NanList, \
    downEastNanList, upWestNanList, downWestNanList, upEastNanList = findNan(dataNO3, nbry)

    daylist = [firstday]

    listNO3prim = [dataNO3[0][CPlat, CPlon]]
    listNO3 = [dataNO3[0][CPlat, CPlon]]

    listNH4prim = [dataNH4[0][CPlat, CPlon]]
    listNH4 = [dataNH4[0][CPlat, CPlon]]

    listDprim = [D[CPlat, CPlon]]

    listN_fprim = [N_f[CPlat, CPlon]]

    listN_sprim = [N_s[CPlat, CPlon]]

    zeroNO3 = np.zeros((nbrx, nbry))
    zeroNH4 = np.zeros((nbrx, nbry))
    zeroD = np.zeros((nbrx, nbry))
    zeroN_f = np.zeros((nbrx, nbry))
    zeroN_s = np.zeros((nbrx, nbry))
    maxCFL = 0
    for k in range(nbrStep):
        dayNbr = k // int(discr)
        hNbr = k // int(discr)
        cNO3[maskpos2D] = 0
        cNH4[maskpos2D] = 0
        D[maskpos2D] = 0
        N_f[maskpos2D] = 0
        N_s[maskpos2D] = 0
        CFL, cx, cy = giveCFL(dxlist, dyMeter, dt, Ewc[hNbr], Nwc[hNbr], nbrx, nbry)
        # print('CFL max: ', max(CFL))
        alpha1 = (1 / 6) * (1 - CFL) * (2 - CFL)
        alpha2 = (1 / 6) * (1 - CFL) * (1 + CFL)
        if np.max(CFL)>maxCFL:
            maxCFL = np.max(CFL)
        matE = createMatE(nbrx, nbry, Ewc[hNbr], eastNanList, westNanList, east2NanList, west2NanList, alpha1, alpha2)
        matN = createMatN(nbrx, nbry, Nwc[hNbr], downNanList, upNanList, down2NanList, up2NanList, alpha1, alpha2)

        days = (daylist[-1] - datetime.datetime(daylist[-1].year, 1, 1,
                                                0)).days  # we get the number of the day in the year

        oldNO3Cline = cNO3.reshape(nbrx * nbry)
        oldNH4Cline = cNH4.reshape(nbrx * nbry)
        advNO3 = - cy * matN.dot(oldNO3Cline) - cx * matE.dot(oldNO3Cline) + (dt / dyMeter) * Ks * matN.dot(
            matN.dot(oldNO3Cline)) + (dt / dxlist.reshape(nbrx * nbry)) * Ks * matE.dot(matN.dot(oldNO3Cline))
        advNH4 = - cy * matN.dot(oldNH4Cline) - cx * matE.dot(oldNH4Cline) + (
                dt / dyMeter) * Ks * matN.dot(matN.dot(oldNH4Cline)) + (
                         dt / dxlist.reshape(nbrx * nbry)) * Ks * matE.dot(matN.dot(oldNH4Cline))

        derivArray = giveEpsilon(days, dataTemp[dayNbr], cNH4 + dataNH4[dayNbr], cNO3 + dataNO3[dayNbr], cNO3, cNH4,
                                 advNO3, advNH4, N_s, N_f, D,
                                 centNwc[hNbr], centEwc[hNbr], latRef, model, nbrx, nbry, dt)

        CNO3line = oldNO3Cline +advNO3 + derivArray[1].reshape(nbrx * nbry) * dt
        cNO3 = np.minimum(CNO3line.reshape(nbrx, nbry), 0)
        #cNO3 = CNO3line.reshape(nbrx, nbry)

        CNH4line = oldNH4Cline +advNH4 + derivArray[0].reshape(nbrx * nbry) * dt
        cNH4 = np.minimum(CNH4line.reshape(nbrx, nbry), 0)
        #cNH4 = CNH4line.reshape(nbrx, nbry)

        oldDCline = D.reshape(nbrx * nbry)
        CDline = oldDCline - cy * matN.dot(oldDCline) - cx * matE.dot(oldDCline) + (
                dt / dyMeter) * Ks * matN.dot(matN.dot(oldDCline)) + (
                         dt / dxlist.reshape(nbrx * nbry)) * Ks * matE.dot(matN.dot(oldDCline)) + derivArray[4].reshape(
            nbrx * nbry) * dt
        D = np.maximum(CDline.reshape(nbrx, nbry),1e-4)

        N_s += derivArray[2] * dt

        N_f += derivArray[3] * dt

        N_f = np.maximum(N_f,-1e-6)
        listNO3 += [dataNO3[dayNbr][CPlat, CPlon]]

        zeroNO3[np.where((cNO3 + dataNO3[dayNbr]) < 0)] += 1
        zeroNH4[np.where((cNH4 + dataNH4[dayNbr]) < 0)] += 1
        zeroD[np.where(D < 0)] += 1
        zeroN_f[np.where(N_f < 0)] += 1
        zeroN_s[np.where(N_s < 0)] += 1

        listNH4 += [dataNH4[dayNbr][CPlat, CPlon]]
        if k % int(30 * discr) == 0:
            daylist += [firstday + relativedelta(minutes=int(k * 24 * 60 / discr))]
            print(k / discr)
            print(time.time() - init)
            print('Max CFL',maxCFL)
            listNH4prim += [cNH4[CPlat, CPlon] + dataNH4[dayNbr][CPlat, CPlon]]
            listNO3prim += [cNO3[CPlat, CPlon] + dataNO3[dayNbr][CPlat, CPlon]]
            listDprim += [D[CPlat, CPlon]]
            listN_fprim += [N_f[CPlat, CPlon]]
            listN_sprim += [N_s[CPlat, CPlon]]


    fig, ax = plt.subplots()
    CpnIm = copy.deepcopy(zeroNO3)
    CpnIm[maskpos2D] = np.nan
    plt.imshow(CpnIm)
    plt.colorbar()
    ax.invert_yaxis()
    ax.set_title("zeros NO3")

    fig, ax = plt.subplots()
    CpnIm = copy.deepcopy(zeroNH4)
    CpnIm[maskpos2D] = np.nan
    plt.imshow(CpnIm)
    plt.colorbar()
    ax.invert_yaxis()
    ax.set_title("zeros NH4")

    fig, ax = plt.subplots()
    CpnIm = copy.deepcopy(zeroN_f)
    CpnIm[maskpos2D] = np.nan
    plt.imshow(CpnIm)
    plt.colorbar()
    ax.invert_yaxis()
    ax.set_title("zeros N_f")

    fig, ax = plt.subplots()
    CpnIm = copy.deepcopy(zeroN_s)
    CpnIm[maskpos2D] = np.nan
    plt.imshow(CpnIm)
    plt.colorbar()
    ax.invert_yaxis()
    ax.set_title("zeros N_s")

    fig, ax = plt.subplots()
    CpnIm = copy.deepcopy(zeroD)
    CpnIm[maskpos2D] = np.nan
    plt.imshow(CpnIm)
    plt.colorbar()
    ax.invert_yaxis()
    ax.set_title("zeros D")

    plt.show()

    fig, ax = plt.subplots()
    CpnIm = copy.deepcopy(cNO3 + dataNO3[dayNbr - 1])
    plt.imshow(CpnIm)
    plt.colorbar()
    plt.clim(0, 280)
    ax.invert_yaxis()
    ax.set_title("NO3'")

    fig, ax = plt.subplots()
    CpnIm = copy.deepcopy(cNH4 + dataNH4[dayNbr - 1])
    plt.imshow(CpnIm)
    plt.colorbar()
    plt.clim(0, 140)
    ax.invert_yaxis()
    ax.set_title("NH4'")

    fig, ax = plt.subplots()
    ImD = copy.deepcopy(D)
    ImD[maskpos2D] = np.nan
    plt.imshow(ImD)
    plt.colorbar()
    ax.invert_yaxis()
    ax.set_title("D")

    fig, ax = plt.subplots()
    CpnIm = copy.deepcopy(N_s)
    plt.imshow(CpnIm)
    plt.colorbar()
    ax.invert_yaxis()
    ax.set_title("N_s")

    fig, ax = plt.subplots()
    CpnIm = copy.deepcopy(N_f)
    plt.imshow(CpnIm)
    plt.colorbar()
    ax.invert_yaxis()
    ax.set_title("N_f")

    fig, ax = plt.subplots()
    CpnIm = copy.deepcopy(derivArray[1] * dt)
    plt.imshow(CpnIm)
    plt.colorbar()
    plt.clim(0, 280)
    ax.invert_yaxis()
    ax.set_title("dNO3'")

    fig, ax = plt.subplots()
    CpnIm = copy.deepcopy(derivArray[0] * dt)
    CpnIm[maskpos2D] = np.nan
    plt.imshow(CpnIm)
    plt.colorbar()
    plt.clim(0, 114)
    ax.invert_yaxis()
    ax.set_title("dNH4'")

    fig, ax = plt.subplots()
    CpnIm = copy.deepcopy(derivArray[4] * dt)
    CpnIm[maskpos2D] = np.nan
    plt.imshow(CpnIm)
    plt.colorbar()
    ax.invert_yaxis()
    ax.set_title("dD")

    fig, ax = plt.subplots()
    CpnIm = copy.deepcopy(derivArray[2] * dt)
    CpnIm[maskpos2D] = np.nan
    plt.imshow(CpnIm)
    plt.colorbar()
    ax.invert_yaxis()
    ax.set_title("dN_s")

    fig, ax = plt.subplots()
    CpnIm = copy.deepcopy(derivArray[3] * dt)
    CpnIm[maskpos2D] = np.nan
    plt.imshow(CpnIm)
    plt.colorbar()
    ax.invert_yaxis()
    ax.set_title("dN_f")

    plt.show()

    '''fig, ax = plt.subplots()
    plt.plot(daylist, listNO3, label='NO3')
    ax.plot(daylist, listNO3prim, label="NO3' daily")
    ax.legend()
    ax.set_xlabel('date')
    ax.set_ylabel('Nitrate')

    fig, ax = plt.subplots()
    plt.plot(daylist, listNH4, label='NH4')
    ax.plot(daylist, listNH4prim, label="NH4' daily")
    ax.legend()
    ax.set_xlabel('date')
    ax.set_ylabel('Ammonium')

    fig, ax = plt.subplots()
    ax.plot(daylist, listDprim, label="D daily")
    ax.legend()
    ax.set_xlabel('date')
    ax.set_ylabel('Detritus')

    fig, ax = plt.subplots()
    ax.plot(daylist, listN_fprim, label="N_f' daily")
    ax.legend()
    ax.set_xlabel('date')
    ax.set_ylabel('fixed nitrate')

    fig, ax = plt.subplots()
    ax.plot(daylist, listN_sprim, label="N_s' daily")
    ax.legend()
    ax.set_xlabel('date')
    ax.set_ylabel('stored nitrate')
    plt.show()'''

    fig, ax = plt.subplots()
    CpnIm = copy.deepcopy(cNO3)
    CpnIm[maskpos2D] = np.nan
    plt.imshow(CpnIm)
    plt.colorbar()
    plt.clim(-140, 0)
    ax.invert_yaxis()
    ax.set_title('cNO3')
    plt.show()

    fig, ax = plt.subplots()
    CpnIm = copy.deepcopy(cNH4)
    CpnIm[maskpos2D] = np.nan
    plt.imshow(CpnIm)
    plt.colorbar()
    plt.clim(-80, 0)
    ax.invert_yaxis()
    ax.set_title('cNH4')
    plt.show()
    NO3field, NH4field, D, N_f, N_s = cNO3 + dataNO3[dayNbr - 1], cNH4 + dataNH4[dayNbr - 1], D, N_f, N_s
    NO3field[maskpos2D] = np.nan
    NH4field[maskpos2D] = np.nan
    D[maskpos2D] = np.nan
    N_f[maskpos2D] = np.nan
    N_s[maskpos2D] = np.nan
    print(listNH4prim)
    print(listNO3prim)
    print(listDprim)
    print(listN_fprim)
    print(listN_sprim)
    print(daylist)
    return NO3field, NH4field, D, N_f, N_s




if __name__ == "__main__":
    zone = 'IBI'
    depth = 0

    latmin = None
    latmax = None

    lonmin = None
    lonmax = None

    dxdeg = 0.028
    dydeg = 0.028

    discr = 72
    dt = 1 / discr

    Ks = 1e-3
    Ks *= 60 * 60 * 24

    CPlat, CPlon = 153, 56

    firstday = datetime.datetime.strptime('2020-01-18', '%Y-%m-%d')

    dataCmdpath = 'D:/Profils/mjaouen/Documents/alternance/EASME/gitcode/shellfish_and_algae-MODEL/dataimport/src/dataCmd.csv'
    _, ds = getwantedMergeData('Nitrate', depth, dataCmdpath)
    '''dataNH4, ds = getwantedMergeData('Ammonium', depth, dataCmdpath)
    dataTemp, ds = getwantedMergeData('Temperature', depth, dataCmdpath)'''

    input_args = {
        'zone': "IBI",
        'file_adress': 'D:/Profils/mjaouen/Documents/alternance/EASME/data/merged_{param}_{zone}.nc',
        'dataRef': pd.read_csv(dataCmdpath, delimiter=';'),
        'paramNames': ['Ammonium', 'Nitrate', 'Temperature', 'northward_Water_current', 'eastward_Water_current']
    }

    ### Initialize the netcdf reading interface
    algaeData = open_data_input(**input_args)

    sim_area = {
        'longitude': (lonmin, lonmax),
        'latitude': (latmin, latmax),
        'depth': depth
    }

    dataNO3 = algaeData.parameterData['Nitrate'].getVariable(**sim_area)[0]
    dataNH4 = algaeData.parameterData['Ammonium'].getVariable(**sim_area)[0]
    dataTemp = algaeData.parameterData['Temperature'].getVariable(**sim_area)[0]
    dataNwc = algaeData.parameterData['northward_Water_current'].getVariable(**sim_area)[0]
    dataEwc = algaeData.parameterData['eastward_Water_current'].getVariable(**sim_area)[0]

    fileU = 'D:/Profils/mjaouen/Documents/alternance/EASME/data/cmems_mod_ibi_phy_u0.nc'
    dataBaseEwc = nc.Dataset(fileU)

    fileV = 'D:/Profils/mjaouen/Documents/alternance/EASME/data/cmems_mod_ibi_phy_v0.nc'
    dataBaseNwc = nc.Dataset(fileV)

    dataFin = pd.read_csv('./../dataimport/src/dataCmd.csv', ';')

    ewcDataLine = dataFin.loc[(dataFin["Parameter"] == 'eastward_Water_current') & (dataFin["Place"] == zone)]
    ewcdataName = ewcDataLine.iloc[-1]["variable"]  # we find the data name in the dataset

    nwcDataLine = dataFin.loc[(dataFin["Parameter"] == 'northward_Water_current') & (dataFin["Place"] == zone)]
    nwcdataName = nwcDataLine.iloc[-1]["variable"]  # we find the data name in the dataset

    longitudeMin, latitudeMin = givecoor(dataBaseNwc, lonmin, latmin, 'northward_Water_current',
                                         dataFin)  # we get the indices of the wanted position
    longitudeMax, latitudeMax = givecoor(dataBaseNwc, lonmax, latmax, 'northward_Water_current',
                                         dataFin)  # we get the indices of the wanted position

    # we get daily data
    '''dataEwc, ds = getwantedMergeData('eastward_Water_current', depth, dataCmdpath)
    dataNwc, ds = getwantedMergeData('northward_Water_current', depth, dataCmdpath)
    dataNwc = dataNwc[:, latitudeMin:latitudeMax, longitudeMin:longitudeMax]
    dataEwc = dataEwc[:, latitudeMin:latitudeMax, longitudeMin:longitudeMax]
    dataNO3 = dataNO3[:, latitudeMin:latitudeMax, longitudeMin:longitudeMax]*14
    dataNH4 = dataNH4[:, latitudeMin:latitudeMax, longitudeMin:longitudeMax]*14
    dataTemp = dataTemp[:, latitudeMin:latitudeMax, longitudeMin:longitudeMax]

    dataNwc *= 60 * 60 * 24
    dataEwc *= 60 * 60 * 24'''

    '''    print(type(dataEwc))
    dataEwc[np.where(np.abs(dataEwc) >1e8)]=0
    dataNwc[np.where(np.abs(dataNwc) >1e8)] = 0'''

    latRef = np.ones((np.shape(dataEwc[0])[1], np.shape(dataEwc[0])[0])) * np.array(
        dataBaseEwc[ewcDataLine.iloc[-1]["latName"]][latitudeMin:latitudeMax])

    decenturedEwc = u2d_cgrid_cur(dataEwc)
    decenturedNwc = v2d_cgrid_cur(dataNwc)

    dxlist, dyMeter = degrees_to_meters(dxdeg, dydeg, latRef)

    dxlist = dxlist.T
    dylist = dyMeter * np.ones(np.shape(dxlist))

    model_params = "D:/Profils/mjaouen/Documents/alternance/EASME/gitcode/shellfish_and_algae-MODEL/macroalgae/macroalgae_model_parameters_input.json"
    json_data = import_json(model_params)

    model = MA_model_scipy(json_data['parameters'])

    maxCFL = 0
    (nbrx, nbry) = np.shape(dataNO3[0])
    for i in range(len(dataEwc)):
        CFL, cx, cy = giveCFL(dxlist, dyMeter, dt, decenturedEwc[i], decenturedNwc[i], nbrx, nbry)
        if np.max(CFL)>maxCFL:
            maxCFL = np.max(CFL)
            print(maxCFL)

    NO3field, NH4field,D, N_f, N_s =  quickest(dyMeter, dxlist, dt, discr, decenturedEwc, decenturedNwc, dataEwc, dataNwc, latRef.T, dataNO3, dataNH4,
             dataTemp, Ks, firstday, model)
    xsize, ysize, ulx, uly, xres, yres = getMetadata(ds,latitudeMin,latitudeMax,longitudeMin,longitudeMax)
    saveAsTiff(NO3field, xsize, ysize, ulx, uly, xres, yres, "I:/work-he/apps/safi/data/IBI/NO3field.tiff")
    saveAsTiff(NH4field, xsize, ysize, ulx, uly, xres, yres,
               "I:/work-he/apps/safi/data/IBI/NH4field.tiff")
    saveAsTiff(D, xsize, ysize, ulx, uly, xres, yres,
               "I:/work-he/apps/safi/data/IBI/D.tiff")
    saveAsTiff(N_f, xsize, ysize, ulx, uly, xres, yres,
               "I:/work-he/apps/safi/data/IBI/N_f.tiff")
    saveAsTiff(N_s, xsize, ysize, ulx, uly, xres, yres,
               "I:/work-he/apps/safi/data/IBI/N_s.tiff")
    saveAsTiff(N_f/10, xsize, ysize, ulx, uly, xres, yres,
               "I:/work-he/apps/safi/data/IBI/biomass.tiff")