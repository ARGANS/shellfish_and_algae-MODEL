import copy
import os
import time
import netCDF4 as nc
import numpy as np
import pandas as pd
import numpy.ma as ma
from scipy import interpolate
from scipy.sparse import dia_matrix
from scipy.sparse import spdiags
import datetime
from dateutil.relativedelta import *

# extract the data value at depth in the merged files (all the daily data merged in one file)
from advectionPrototype.saveAsTiff import saveAsTiff, giveMetadata
from dataread.launch_model import MA_model_scipy
from dataread.make_runs import open_data_input
from dataread.read_netcdf import extractVarSubset, AllData
from dataread.utils import import_json


class resample:
    def __init__(self, dxRatio, dyRatio):
        self.dxRatio = dxRatio
        self.dyRatio = dyRatio

    def findElmt(self, i, j):
        newi = np.round(i * self.dxRatio).astype(int)
        newj = np.round(j * self.dyRatio).astype(int)
        return newi, newj

    def giveNewMatCoor(self, matshape, imin, jmin):
        rowCoor = np.zeros(matshape[0] * matshape[1])
        columnCoor = np.zeros(matshape[0] * matshape[1])
        ival, jval = self.findElmt(np.arange(matshape[0]),
                                   np.arange(matshape[1]))  # comment retrouver ou on était dans new
        ival, jval = ival + imin, jval + jmin
        for k in range(matshape[0]):
            rowCoor[k * matshape[1]:(k + 1) * matshape[1]] = ival[k]
            columnCoor[k * matshape[1]:(k + 1) * matshape[1]] = jval
        return rowCoor.astype(int), columnCoor.astype(int)


def resampleData(resa, dxlist, latRef, dataNH4, dataPAR, dataNO3, dataEwc, dataNwc, dataTemp, newArrayShape, i, j):
    rowCoor, columnCoor = resa.giveNewMatCoor(newArrayShape, i, j)
    dxlistarray = dxlist[rowCoor, columnCoor].reshape(newArrayShape)
    latRefarray = latRef.T[rowCoor, columnCoor].reshape(newArrayShape)

    dataNH4array = dataNH4[:, rowCoor, columnCoor].reshape((len(dataNH4),) + newArrayShape)
    dataPARarray = dataPAR[:, rowCoor, columnCoor].reshape((len(dataPAR),) + newArrayShape)
    dataNO3array = dataNO3[:, rowCoor, columnCoor].reshape((len(dataNO3),) + newArrayShape)
    dataEwcarray = dataEwc[:, rowCoor, columnCoor].reshape((len(dataEwc),) + newArrayShape)
    dataNwcarray = dataNwc[:, rowCoor, columnCoor].reshape((len(dataNwc),) + newArrayShape)
    dataTemparray = dataTemp[:, rowCoor, columnCoor].reshape((len(dataTemp),) + newArrayShape)
    decenturedEwcarray = u2d_cgrid_cur(dataEwcarray)
    decenturedNwcarray = v2d_cgrid_cur(dataNwcarray)
    return dxlistarray, latRefarray, dataNH4array, dataPARarray, dataNO3array, dataTemparray, dataEwcarray, dataNwcarray, decenturedEwcarray, decenturedNwcarray


# this function decentralize the u speeds
def u2d_cgrid_cur(var):
    # var: 3D array ( time, lat, lon)
    varcg = np.zeros((var.shape[0], var.shape[1], var.shape[2] + 1))
    for i in range(var.shape[1]):
        for j in range(var.shape[2] - 1):
            if ((not ma.is_masked(var[0, i, j])) & (not ma.is_masked(var[0, i, j + 1]))):
                varcg[:, i, j] = (var[:, i, j] + var[:, i, j + 1]) / 2
    return varcg


# this function decentralize the v speeds
def v2d_cgrid_cur(var):
    varcg = np.zeros((var.shape[0], var.shape[1] + 1, var.shape[2]))
    for i in range(var.shape[1] - 1):
        for j in range(var.shape[2]):
            if ((not ma.is_masked(var[0, i, j])) & (not ma.is_masked(var[0, i + 1, j]))):
                varcg[:, i, j] = (var[:, i, j] + var[:, i + 1, j]) / 2
    return varcg


def degrees_to_meters(lonDist, latDist, refLat):
    """Converts degrees of latDist,lonDist to meters assuming that the latitude
    stays near refLat.
    """
    lat_degree = 111000  # conversion of latitude degrees to meters

    Dx = lonDist * lat_degree * np.cos(np.deg2rad(refLat))
    Dy = latDist * lat_degree

    return Dx, Dy


# return the CFL number, cx and cy
def giveCFL(dx, dy, dt, Ewc, Nwc, nbrx, nbry):
    Cx = np.abs(Ewc[:, 1:] * dt / dx).reshape(nbrx * nbry)
    Cy = np.abs(Nwc[1:, :] * dt / dy).reshape(nbrx * nbry)
    return np.maximum(Cx, Cy), Cx, Cy


# take as input a mask (tuple), and return the position of the masked values in a vector (dim x * dim y)
def createListNan(maskPosition, nbry):
    listNan = []
    for k in range(len(maskPosition[0])):
        i, j = maskPosition[0][k], maskPosition[1][k]
        n = j + nbry * i
        listNan.append(n)
    return listNan


# return the nan or masked data position
def findNan(dataNO3):
    # dataNO3: 3D array (time, lat, lon)

    # Why the copy ?
    dataNO3Copy = copy.deepcopy(dataNO3[0])
    mask = dataNO3Copy.mask
    nbrx, nbry = mask.shape

    nanPositions = np.empty((5, 5), dtype=np.object)
    # examples:
    #   - nanPositions[0,0] = mask
    #   - nanPositions[0,1] = eastNan
    #   - nanPositions[0,-1] = westNan
    #   - nanPositions[2,0] = up2Nan / north2Nan
    #   - nanPositions[1,-1] = upWestNan / northWestNan

    for iRel in [-2, -1, 0, 1, 2]:  # along latitude
        for jRel in [-2, -1, 0, 1, 2]:  # along longitude
            # Matrix to shift mask along E-W by -jRel
            shiftMatE = spdiags([1] * nbry, -jRel, nbry, nbry).todense()
            # Matrix to shift mask along N-S by -iRel
            shiftMatN = spdiags([1] * nbrx, iRel, nbrx, nbrx).todense()
            nanPositions[iRel, jRel] = np.flatnonzero(shiftMatN.dot(mask).dot(shiftMatE))

    return nanPositions


# creates the matrix to compute the flux in u direction
def createMatE(nbrx, nbry, decenturedEwc, nanLists, alpha1, alpha2):
    offset = np.array([0, -1, 1, -2, 2])
    uGreater0 = ((decenturedEwc[:, 1:] > 0) * 1).reshape(nbrx * nbry)
    termA = (1 - 2 * alpha1 + alpha2) * uGreater0 - (1 - uGreater0) * (2 * alpha1 + alpha2 - 1)
    termB = (alpha1 - 2 * alpha2 - 1) * uGreater0 - alpha1 * (1 - uGreater0)
    termC = alpha1 * uGreater0 + (1 - alpha1 - 2 * alpha2) * (1 - uGreater0)
    termD = uGreater0 * alpha2
    termE = alpha2 * (1 - uGreater0)

    termB[::nbry] = 0
    termC[nbry - 1::nbry] = 0

    termA[nanLists[0, -1]] = 0
    termB[nanLists[0, -1]] = 0

    termA[nanLists[0, 1]] = 0
    termC[nanLists[0, 1]] = 0

    termA[nanLists[0, -2]] = 0
    termB[nanLists[0, -2]] = 0

    termA[nanLists[0, 2]] = 0
    termC[nanLists[0, 2]] = 0

    data = np.zeros((5, nbrx * nbry))
    data[0] = termA
    data[1][:-1] = termB[1:]
    data[2][1:] = termC[:-1]
    data[3][:-2] = termD[2:]
    data[4][2:] = termE[:-2]
    Mat = dia_matrix((data, offset), shape=(nbrx * nbry, nbrx * nbry))

    return Mat


# creates the matrix to compute the flux in v direction
def createMatN(nbrx, nbry, decenturedNwc, nanLists, alpha1, alpha2):
    offset = np.array([0, -nbry, nbry, -2 * nbry, 2 * nbry])
    vGreater0 = ((decenturedNwc[1:] > 0) * 1).reshape(nbrx * nbry)
    termA = (1 - 2 * alpha1 + alpha2) * vGreater0 - (1 - vGreater0) * (2 * alpha1 + alpha2 - 1)
    termB = (alpha1 - 2 * alpha2 - 1) * vGreater0 - alpha1 * (1 - vGreater0)
    termC = alpha1 * vGreater0 + (1 - alpha1 - 2 * alpha2) * (1 - vGreater0)
    termD = vGreater0 * alpha2
    termE = alpha2 * (1 - vGreater0)

    termA[nanLists[-1, 0]] = 0
    termB[nanLists[-1, 0]] = 0

    termA[nanLists[1, 0]] = 0
    termC[nanLists[1, 0]] = 0

    termA[nanLists[-2, 0]] = 0
    termB[nanLists[-2, 0]] = 0

    termA[nanLists[2, 0]] = 0
    termC[nanLists[2, 0]] = 0

    data = np.zeros((5, nbrx * nbry))
    data[0] = termA
    data[1][:-nbry] = termB[nbry:]
    data[2][nbry:] = termC[:-nbry]
    data[3][:-2 * nbry] = termD[2 * nbry:]
    data[4][2 * nbry:] = termE[:-2 * nbry]
    Mat = dia_matrix((data, offset), shape=(nbrx * nbry, nbrx * nbry))

    return Mat


# returns the space step
def giveDxDy(latitudes, longitudes, latRef):
    Dx = np.zeros((len(longitudes), len(latitudes)))
    Dy = np.zeros((len(latitudes), len(longitudes)))
    Dy[:, :-1] = longitudes[1:] - longitudes[:-1]
    Dy[:, -1] = Dy[:, -2]
    Dx[:, :-1] = latitudes[1:] - latitudes[:-1]
    Dx[:, -1] = Dx[:, -2]
    DxMeter, DyMeter = degrees_to_meters(Dx, Dy, latRef)
    return DxMeter.T, DyMeter


# return the variation of each quantities in the algae model
def giveEpsilon(day, temp, NH4, NO3, cNO3, cNH4, advNO3, advNH4, N_s, N_f, D, Nwc, Ewc, PAR, latRef, model, nbrx, nbry,
                dt, Zmix):
    data_in = {
        'SST': temp,
        'PAR': PAR,
        'NH4_ext': NH4,
        'NO3_ext': NO3,
        'PO4_ext': 50,
        'K_d': 0.1,
        'F_in': np.sqrt(Nwc ** 2 + Ewc ** 2),
        'h_z_SML': 30,
        't_z': Zmix,
        'D_ext': 0.1
    }
    y = dict(zip(["NH4", "NO3", "N_s", "N_f", "D"], [NH4, NO3, N_s, N_f, D]))
    return model.derivative_fast_advection(day, y, data_in, cNO3, cNH4, advNO3, advNH4, nbrx, nbry, dt,
                                           latRef, model)


# get in the dataLine comming from dataCmd.csv the resolution of the data
def giveResol(dataLine):
    resolString = dataLine.iloc[-1]["resolution"]
    splitedString = resolString.split('*')
    if len(splitedString[0].split('km')) == 2:
        return float(splitedString[0].split('km')[0]) * 1e3, float(splitedString[1].split('km')[0]) * 1e3, True
    else:
        return float(splitedString[0]) * 1e-2, float(splitedString[1]) * 1e-2, False


def prepareDerivArrayScenC(derivArray, nanLists):
    for i in range(len(derivArray)):
        for nanL in nanLists.flatten():
            derivArrayLine = derivArray[i].reshape(nbrx * nbry)
            derivArrayLine[nanL] = 1e-10
            derivArray[i] = derivArrayLine.reshape(nbrx, nbry)
    return derivArray


def sortPAR(dateBeginning, dataArr):
    datetimeBeginning = datetime.datetime.strptime(dateBeginning, '%Y-%m-%d %H:%M:%S')
    firstmonthNbr = datetimeBeginning.month
    dataArr = dataArr.filled(fill_value=-999)
    dataArrShape = np.shape(dataArr)
    newdataArr = np.zeros(dataArrShape)
    newdataArr[:(dataArrShape[0] - firstmonthNbr)] = dataArr[firstmonthNbr:]
    newdataArr[(dataArrShape[0] - firstmonthNbr):] = dataArr[:firstmonthNbr]
    return ma.masked_outside(newdataArr, -1e-2, 1e4)


def sortData(dateBeginning, dateEnd, lenDataArr):
    datetimeBeginning = datetime.datetime.strptime(dateBeginning, '%Y-%m-%d %H:%M:%S')
    datetimeEnd = datetime.datetime.strptime(dateEnd, '%Y-%m-%d %H:%M:%S')
    firstdayNbr = (datetimeBeginning - datetime.datetime(datetimeBeginning.year, 1, 1, 0)).days
    lastdayNbr = (datetimeEnd - datetime.datetime(datetimeEnd.year, 1, 1, 0)).days
    if firstdayNbr < lastdayNbr:
        return np.arange(firstdayNbr, lastdayNbr)
    elif firstdayNbr == lastdayNbr:
        return np.arange(0, lenDataArr)
    else:
        return np.concatenate([np.arange(firstdayNbr, lenDataArr), np.arange(0, lastdayNbr)])


# apply the quickest scheme
def quickest(dyMeter, dxlist, dt, Ewc, Nwc, centEwc, centNwc, latRef, dataNO3, dataNH4, dataTemp, dataPAR, Ks, firstday,
             model, Zmix, scenC, sortedList):
    # we initate the variables
    discr = 1 / dt
    nbrStep = int(len(sortedList) * discr)
    cNO3 = np.zeros(np.shape(dataNO3[0]))
    cNH4 = np.zeros(np.shape(dataNH4[0]))
    N_s = np.ones(np.shape(dataNH4[0])) * 1000
    N_f = np.ones(np.shape(dataNH4[0])) * 1000
    D = np.ones(np.shape(dataNH4[0])) * 0.1

    totalNH4deficit = np.zeros(np.shape(dataNH4[0]))
    totalNO3deficit = np.zeros(np.shape(dataNO3[0]))
    (nbrx, nbry) = np.shape(dataNO3[0])
    mask = dataNO3.mask
    maskpos2D = np.where(mask[0] == True)
    init = time.time()
    nanLists = findNan(dataNO3)

    daylist = [firstday]

    maxCFL = 0
    dx = dxlist.reshape(nbrx * nbry)
    dy = dyMeter.reshape(nbrx * nbry)
    # for each time step
    for k in range(nbrStep):
        daylist += [firstday + relativedelta(minutes=int(k * 24 * 60 / discr))]
        # we compute the day, hour and month number
        dayNbr = sortedList[k // int(discr)]
        hNbr = dayNbr  # if we don't use hourly speeds, hNbr = dayNbr
        month = dayNbr // int(30.5 * discr)
        cNO3[maskpos2D] = 1e-5
        cNH4[maskpos2D] = 1e-5
        D[maskpos2D] = 1e-5
        N_f[maskpos2D] = 1e-5
        N_s[maskpos2D] = 1e-5
        # we compute the CFL
        CFL, cx, cy = giveCFL(dxlist, dyMeter, dt, Ewc[hNbr], Nwc[hNbr], nbrx, nbry)
        alpha1 = (1 / 6) * (1 - CFL) * (2 - CFL)
        alpha2 = (1 / 6) * (1 - CFL) * (1 + CFL)
        if np.max(CFL) > maxCFL:
            maxCFL = np.max(CFL)
        matE = createMatE(nbrx, nbry, Ewc[hNbr], nanLists, alpha1, alpha2)
        matN = createMatN(nbrx, nbry, Nwc[hNbr], nanLists, alpha1, alpha2)

        days = (daylist[-1] - datetime.datetime(daylist[-1].year, 1, 1,
                                                0)).days  # we get the number of the day in the year

        oldNO3Cline = cNO3.reshape(nbrx * nbry)
        oldNH4Cline = cNH4.reshape(nbrx * nbry)
        no3Hatu = matE.dot(oldNO3Cline)
        no3Hatv = matN.dot(oldNO3Cline)
        nh4Hatu = matE.dot(oldNH4Cline)
        nh4Hatv = matN.dot(oldNH4Cline)
        # we compute the advection terms
        advNO3 = - cy * no3Hatv - cx * no3Hatu + (dt / dy) * Ks * matN.dot(no3Hatv) + (
                dt / dx) * Ks * matE.dot(no3Hatu)
        advNH4 = - cy * nh4Hatv - cx * nh4Hatu + (dt / dy) * Ks * matN.dot(nh4Hatv) + (
                dt / dx) * Ks * matE.dot(nh4Hatu)

        derivArray = giveEpsilon(days, dataTemp[dayNbr], cNH4 + dataNH4[dayNbr], cNO3 + dataNO3[dayNbr], cNO3, cNH4,
                                 advNO3, advNH4, N_s, N_f, D,
                                 centNwc[hNbr], centEwc[hNbr], dataPAR[month], latRef, model, nbrx, nbry, dt, Zmix)
        if scenC:
            derivArray = prepareDerivArrayScenC(derivArray, nanLists)

        CNO3line = oldNO3Cline + advNO3 + derivArray[1].reshape(nbrx * nbry) * dt
        cNO3 = np.minimum(CNO3line.reshape(nbrx, nbry), 0)
        # cNO3 = CNO3line.reshape(nbrx, nbry)

        totalNO3deficit += derivArray[1]

        CNH4line = oldNH4Cline + advNH4 + derivArray[0].reshape(nbrx * nbry) * dt
        cNH4 = np.minimum(CNH4line.reshape(nbrx, nbry), 0)
        # cNH4 = CNH4line.reshape(nbrx, nbry)

        totalNH4deficit += derivArray[0]

        oldDCline = D.reshape(nbrx * nbry)
        CDline = oldDCline - cy * matN.dot(oldDCline) - cx * matE.dot(oldDCline) + (
                dt / dy) * Ks * matN.dot(matN.dot(oldDCline)) + (
                         dt / dx) * Ks * matE.dot(matN.dot(oldDCline)) + derivArray[4].reshape(
            nbrx * nbry) * dt
        D = np.maximum(CDline.reshape(nbrx, nbry), 1e-4)

        N_s += derivArray[2] * dt

        N_f += derivArray[3] * dt

        N_f = np.maximum(N_f, -1e-6)
        N_s = np.maximum(N_s, -1e-6)

        if k % int(30 * discr) == 0:
            print(k / discr)
            print(time.time() - init)
            print('Max CFL', maxCFL)

    NO3field, NH4field, D, N_f, N_s = cNO3 + dataNO3[dayNbr - 1], cNH4 + dataNH4[dayNbr - 1], D, N_f, N_s
    NO3field[maskpos2D] = np.nan
    NH4field[maskpos2D] = np.nan
    D[maskpos2D] = np.nan
    N_f[maskpos2D] = np.nan
    N_s[maskpos2D] = np.nan
    totalNH4deficit[maskpos2D] = np.nan
    totalNO3deficit[maskpos2D] = np.nan
    return NO3field, NH4field, D, N_f, N_s, totalNH4deficit, totalNO3deficit


if __name__ == "__main__":
    zone = "NWS"
    depth = 0
    PAR_year = 2020
    getAmmonium = False

    dateBeginning = '2020-09-01 00:00:00'
    dateEnd = '2020-04-30 00:00:00'

    latmin = None
    latmax = None

    lonmin = None
    lonmax = None

    scenC = True

    # discr = 144
    discr = 72  # Baltic, NWS, IBI
    # discr = 48 #BS
    dt = 1 / discr

    Ks = 1e-3
    Ks *= 60 * 60 * 24

    CPlat, CPlon = 153, 56

    model_params = "./../macroalgae/macroalgae_model_parameters_input.json"
    json_data = import_json(model_params)

    paramSacch = json_data['parameters']['species']['saccharina']['parameters']
    model = MA_model_scipy(json_data['parameters'])

    Q_min = paramSacch['Q_min']
    density_MA = 0.4
    h_MA = paramSacch['h_MA']
    w_MA = paramSacch['w_MA']
    DF_MA = 0.113
    kcal_MA = 2.29
    prot_MA = 0.08
    CN_MA = 21

    Zmix = h_MA * 1.4

    if getAmmonium:
        paramNames = ['Nitrate', 'northward_Water_current', 'Ammonium', 'eastward_Water_current',
                      'Temperature']  # IBI, Baltic and MED
    else:
        paramNames = ['Nitrate', 'northward_Water_current', 'eastward_Water_current', 'Temperature']  # NWS and BS

    firstday = datetime.datetime.strptime(dateBeginning, '%Y-%m-%d %H:%M:%S')

    dataCmdpath = './../global/dataCmd.csv'
    mergedFilepath = 'D:/Profils/mjaouen/Documents/alternance/EASME/data/{zone}/'.format(zone=zone)

    input_args = {
        'zone': zone,
        'file_adress': 'D:/Profils/mjaouen/Documents/alternance/EASME/data/{zone}/merged_{param}_{zone}.nc',
        'dataRef': pd.read_csv(dataCmdpath, delimiter=';'),
        'paramNames': paramNames,
        'frequency': 'daily'
        # 'with_PAR': PAR_year,
        # 'PAR_file': 'D:/Profils/mjaouen/Documents/alternance/EASME/data/{zone}/PAR_{zone}NewGrid.nc'.format(zone = zone)
    }

    dict_to_AllData = open_data_input(**input_args)

    dict_to_AllData['PAR'] = {
        'file_name': 'D:/Profils/mjaouen/Documents/alternance/EASME/data/{zone}/PAR_{zone}NewGrid.nc'.format(zone=zone),
        'variable_name': 'par',
        'latitude_name': 'lat',
        'longitude_name': 'lon',
        'time_name': 'time',
        'depth_name': 'depth',
        'unit_conversion': 11.574,
        'time_zero': datetime.datetime(PAR_year, 1, 1),
        'time_step': datetime.timedelta(days=1)
    }

    sim_area = {
        'longitude': (lonmin, lonmax),
        'latitude': (latmin, latmax),
        'depth': 0
        # 'depth': (0, Zmix),
        # 'averagingDims': ('depth',),
        # 'weighted': False
    }

    ### Initialize the netcdf reading interface
    algaeData = AllData(dict_to_AllData)

    dataNO3 = algaeData.parameterData['Nitrate'].getVariable(**sim_area)[0]
    if getAmmonium:
        dataNH4 = algaeData.parameterData['Ammonium'].getVariable(**sim_area)[0]
    else:
        dataNH4 = dataNO3 * 0.1
    dataTemp = algaeData.parameterData['Temperature'].getVariable(**sim_area)[0]
    dataNwc = algaeData.parameterData['northward_Water_current'].getVariable(**sim_area)[0]
    dataEwc = algaeData.parameterData['eastward_Water_current'].getVariable(**sim_area)[0]
    dataPAR = algaeData.parameterData['PAR'].getVariable(**sim_area)[0]

    print(np.shape(dataNO3))
    print(np.shape(dataNH4))
    print(np.shape(dataTemp))
    print(np.shape(dataNwc))
    print(np.shape(dataPAR))

    sortedList = sortData(dateBeginning, dateEnd, len(dataNO3))

    dataNH4 = ma.masked_outside(dataNH4, -1e4, 1e4)
    dataNO3 = ma.masked_outside(dataNO3, -1e4, 1e4)
    dataTemp = ma.masked_outside(dataTemp, -1e4, 1e4)
    dataPAR = ma.masked_outside(dataPAR, -1e-2, 1e4)
    dataPAR = dataPAR.filled(fill_value=8)

    fileV = 'D:/Profils/mjaouen/Documents/alternance/EASME/data/{zone}/merged_northward_Water_current_{zone}.nc'.format(
        zone=zone)
    dataBaseNwc = nc.Dataset(fileV)

    dataFin = pd.read_csv('./../global/dataCmd.csv', ';')

    nwcDataLine = dataFin.loc[(dataFin["Parameter"] == 'Nitrate') & (dataFin["Place"] == zone)]
    nwcdataName = nwcDataLine.iloc[-1]["variable"]  # we find the data name in the dataset
    resx, resy, km = giveResol(nwcDataLine)
    print(resx, resy, km)

    longitudes, _ = algaeData.parameterData['Temperature'].getVariable('longitude', **sim_area)
    latitudes, _ = algaeData.parameterData['Temperature'].getVariable('latitude', **sim_area)

    '''    print(type(dataEwc))
    dataEwc[np.where(np.abs(dataEwc) >1e8)]=0
    dataNwc[np.where(np.abs(dataNwc) >1e8)] = 0'''

    latRef = np.ones((np.shape(dataEwc[0])[1], np.shape(dataEwc[0])[0])) * latitudes

    decenturedEwc = u2d_cgrid_cur(dataEwc)
    decenturedNwc = v2d_cgrid_cur(dataNwc)

    dxlist, dyMeter = giveDxDy(latitudes, longitudes, latRef)

    dxRatio = 1852 / np.mean(dxlist)
    dyRatio = 1852 / np.mean(dyMeter)

    print(Q_min, Zmix)

    maxCFL = 0
    (nbrx, nbry) = np.shape(dataNO3[0])
    for i in range(len(dataEwc)):
        CFL, cx, cy = giveCFL(dxlist, dyMeter, dt, decenturedEwc[i], decenturedNwc[i], nbrx, nbry)
        if np.max(CFL) > maxCFL:
            maxCFL = np.max(CFL)
    print('maxCFL: ' + str(maxCFL))
    xsize, ysize, ulx, uly, xres, yres = giveMetadata(latitudes, longitudes)
    saveAsTiff(dataNO3[0], xsize, ysize, ulx, uly, xres, yres, "I:/work-he/apps/safi/data/IBI/test.tiff")
    resa = resample(dxRatio, dyRatio)
    newArrayShape = (int(50000 / 1852), int(50000 / 1852))
    print(int(50000 / resx))
    print(int(50000 / resy))
    # i, j are the coordinates in the top left corner of the studied area
    for i in range(0, nbrx, int(50000 / resx)):
        for j in range(0, nbry, int(50000 / resy)):
            print(i, j)
            NO3Oldarray = dataNO3[:, i:i + int(50000 / resx), j:j + int(50000 / resy)].filled(np.nan)
            if (np.sum(np.isnan(NO3Oldarray))) and (np.sum(~np.isnan(NO3Oldarray)) > 0):
                '''plt.imshow(NO3Oldarray[0])
                plt.show()'''
                dxlistarray, latRefarray, dataNH4array, dataPARarray, dataNO3array, \
                dataTemparray, dataEwcarray, dataNwcarray, decenturedEwcarray, \
                decenturedNwcarray = resa.resampleData(dxlist, latRef, dataNH4, dataPAR, dataNO3, dataEwc, dataNwc,
                                                       dataTemp, newArrayShape,
                                                       i, j)
                NO3field, NH4field, D, N_f, N_s, totalNH4deficit, totalNO3deficit = quickest(dyMeter, dxlistarray, dt,
                                                                                             decenturedEwcarray,
                                                                                             decenturedNwcarray,
                                                                                             dataEwcarray,
                                                                                             dataNwcarray, latRefarray,
                                                                                             dataNO3array, dataNH4array,
                                                                                             dataTemparray,
                                                                                             dataPARarray, Ks, firstday,
                                                                                             model,
                                                                                             Zmix, scenC, sortedList)
    if int(50000 / resx) != 50000 / resx:
        i = np.shape(dataNO3)[1] - 1 - int(50000 / resx)
        for j in range(0, nbry, int(50000 / resy)):
            print(i, j)
            NO3Oldarray = dataNO3[:, i:i + int(50000 / resx), j:j + int(50000 / resy)].filled(np.nan)
            if (np.sum(np.isnan(NO3Oldarray))) and (np.sum(~np.isnan(NO3Oldarray)) > 0):
                dxlistarray, latRefarray, dataNH4array, dataPARarray, dataNO3array, \
                dataTemparray, dataEwcarray, dataNwcarray, decenturedEwcarray, \
                decenturedNwcarray = resa.resampleData(dxlist, latRef, dataNH4, dataPAR, dataNO3, dataEwc, dataNwc,
                                                       dataTemp, newArrayShape,
                                                       i, j)
                NO3field, NH4field, D, N_f, N_s, totalNH4deficit, totalNO3deficit = quickest(dyMeter, dxlistarray, dt,
                                                                                             decenturedEwcarray,
                                                                                             decenturedNwcarray,
                                                                                             dataEwcarray,
                                                                                             dataNwcarray, latRefarray,
                                                                                             dataNO3array, dataNH4array,
                                                                                             dataTemparray,
                                                                                             dataPARarray, Ks, firstday,
                                                                                             model,
                                                                                             Zmix, scenC, sortedList)
    if int(50000 / resy) != 50000 / resy:
        j = np.shape(dataNO3)[2] - 1 - int(50000 / resy)
        for i in range(0, nbrx, int(50000 / resx)):
            print(i, j)
            NO3Oldarray = dataNO3[:, i:i + int(50000 / resx), j:j + int(50000 / resy)].filled(np.nan)
            if (np.sum(np.isnan(NO3Oldarray))) and (np.sum(~np.isnan(NO3Oldarray)) > 0):
                dxlistarray, latRefarray, dataNH4array, dataPARarray, dataNO3array, \
                dataTemparray, dataEwcarray, dataNwcarray, decenturedEwcarray, \
                decenturedNwcarray = resa.resampleData(dxlist, latRef, dataNH4, dataPAR, dataNO3, dataEwc, dataNwc,
                                                       dataTemp, newArrayShape,
                                                       i, j)
                NO3field, NH4field, D, N_f, N_s, totalNH4deficit, totalNO3deficit = quickest(dyMeter, dxlistarray,
                                                                                             dt,
                                                                                             decenturedEwcarray,
                                                                                             decenturedNwcarray,
                                                                                             dataEwcarray,
                                                                                             dataNwcarray,
                                                                                             latRefarray,
                                                                                             dataNO3array,
                                                                                             dataNH4array,
                                                                                             dataTemparray,
                                                                                             dataPARarray, Ks,
                                                                                             firstday, model,
                                                                                             Zmix, scenC,
                                                                                             sortedList)

    if (int(50000 / resx) != 50000 / resx) and (int(50000 / resy) != 50000 / resy):
        j = np.shape(dataNO3)[2] - 1 - int(50000 / resy)
        i = np.shape(dataNO3)[1] - 1 - int(50000 / resx)
        NO3Oldarray = dataNO3[:, i:i + int(50000 / resx), j:j + int(50000 / resy)].filled(np.nan)
        if (np.sum(np.isnan(NO3Oldarray))) and (np.sum(~np.isnan(NO3Oldarray)) > 0):
            dxlistarray, latRefarray, dataNH4array, dataPARarray, dataNO3array, \
            dataTemparray, dataEwcarray, dataNwcarray, decenturedEwcarray, \
            decenturedNwcarray = resa.resampleData(dxlist, latRef, dataNH4, dataPAR, dataNO3, dataEwc, dataNwc,
                                                   dataTemp, newArrayShape,
                                                   i, j)
            NO3field, NH4field, D, N_f, N_s, totalNH4deficit, totalNO3deficit = quickest(dyMeter, dxlistarray,
                                                                                         dt,
                                                                                         decenturedEwcarray,
                                                                                         decenturedNwcarray,
                                                                                         dataEwcarray,
                                                                                         dataNwcarray,
                                                                                         latRefarray,
                                                                                         dataNO3array,
                                                                                         dataNH4array,
                                                                                         dataTemparray,
                                                                                         dataPARarray, Ks,
                                                                                         firstday, model,
                                                                                         Zmix, scenC,
                                                                                         sortedList)

    DW = N_f / Q_min  # gDW m-3
    DW_line = DW * h_MA * w_MA / 1000  # kg/m (i.e. per m of rope)
    DW_PUA = DW * h_MA * density_MA / 1000  # kg/m^2 (multiply be density to account for unused space within farm)

    FW = DW / DF_MA  # gFW m-3
    FW_line = DW_line / DF_MA  # kg/m (i.e. per m of rope)
    FW_PUA = DW_PUA / DF_MA  # kg/m^2

    # Energy
    kcal_PUA = DW * h_MA * density_MA * kcal_MA  # kcal/m^2

    # protein
    protein_PUA = DW_PUA * prot_MA  # kg/m^2

    # CO2 uptake
    Biomass_CO2 = (N_f / 14) * CN_MA * 44 / 1000  # g (CO2) /m^3    (44 is rmm of CO2)
    CO2_uptake_PUA = Biomass_CO2 * h_MA * density_MA / 1000  # kg (CO2) / m^2

    saveAsTiff(NO3field, xsize, ysize, ulx, uly, xres, yres,
               "I:/work-he/apps/safi/data/{zone}/NO3field.tif".format(zone=zone))
    saveAsTiff(NH4field, xsize, ysize, ulx, uly, xres, yres,
               "I:/work-he/apps/safi/data/{zone}/NH4field.tif".format(zone=zone))
    saveAsTiff(D, xsize, ysize, ulx, uly, xres, yres,
               "I:/work-he/apps/safi/data/{zone}/D.tif".format(zone=zone))
    saveAsTiff(N_f, xsize, ysize, ulx, uly, xres, yres,
               "I:/work-he/apps/safi/data/{zone}/N_f.tif".format(zone=zone))
    saveAsTiff(N_s, xsize, ysize, ulx, uly, xres, yres,
               "I:/work-he/apps/safi/data/{zone}/N_s.tif".format(zone=zone))
    saveAsTiff(DW, xsize, ysize, ulx, uly, xres, yres,
               "I:/work-he/apps/safi/data/{zone}/DW.tif".format(zone=zone))
    saveAsTiff(totalNH4deficit, xsize, ysize, ulx, uly, xres, yres,
               "I:/work-he/apps/safi/data/{zone}/totalNH4deficit.tif".format(zone=zone))
    saveAsTiff(totalNO3deficit, xsize, ysize, ulx, uly, xres, yres,
               "I:/work-he/apps/safi/data/{zone}/totalNO3deficit.tif".format(zone=zone))
    saveAsTiff(DW_PUA, xsize, ysize, ulx, uly, xres, yres,
               "I:/work-he/apps/safi/data/{zone}/DW_PUA.tif".format(zone=zone))
    saveAsTiff(FW, xsize, ysize, ulx, uly, xres, yres,
               "I:/work-he/apps/safi/data/{zone}/FW.tif".format(zone=zone))
    saveAsTiff(FW_line, xsize, ysize, ulx, uly, xres, yres,
               "I:/work-he/apps/safi/data/{zone}/FW_line.tif".format(zone=zone))
    saveAsTiff(FW_PUA, xsize, ysize, ulx, uly, xres, yres,
               "I:/work-he/apps/safi/data/{zone}/FW_PUA.tif".format(zone=zone))
    saveAsTiff(kcal_PUA, xsize, ysize, ulx, uly, xres, yres,
               "I:/work-he/apps/safi/data/{zone}/kcal_PUA.tif".format(zone=zone))
    saveAsTiff(protein_PUA, xsize, ysize, ulx, uly, xres, yres,
               "I:/work-he/apps/safi/data/{zone}/protein_PUA.tif".format(zone=zone))
    saveAsTiff(Biomass_CO2, xsize, ysize, ulx, uly, xres, yres,
               "I:/work-he/apps/safi/data/{zone}/Biomass_CO2.tif".format(zone=zone))
    saveAsTiff(CO2_uptake_PUA, xsize, ysize, ulx, uly, xres, yres,
               "I:/work-he/apps/safi/data/{zone}/CO2_uptake_PUA.tif".format(zone=zone))
    saveAsTiff(DW_line, xsize, ysize, ulx, uly, xres, yres,
               "I:/work-he/apps/safi/data/{zone}/DW_line.tif".format(zone=zone))
