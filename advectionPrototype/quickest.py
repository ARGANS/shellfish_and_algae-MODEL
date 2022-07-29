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

#sort the list of data from the older to the newer
def sortDateList(listValue,ldate):
    sortLval = [] #we define the sorted data list
    sortldate = [] #we define the sorted date list
    # we read the list we have to sort
    for k in range(len(ldate)):
        i = 0
        # while we didn't red all the sorted list
        while i<len(sortLval):
            # if the list element is graeter than the element we want to place in the list
            if ldate[k]<sortldate[i]:
                # we place this element before the element who is greater than it
                sortldate = sortldate[:i]+[ldate[k]]+sortldate[i:]
                sortLval = sortLval[:i] + [listValue[k]] + sortLval[i:]
                i=len(sortLval)+1
            i+=1
        if i == len(sortLval)-1 or i == len(sortLval):
            sortldate += [ldate[k]]
            sortLval += [listValue[k]]
    return np.array(sortLval), np.array(sortldate)

def getData(path,data,zone, depth, latitudeMinwc,latitudeMaxwc, longitudeMinwc,longitudeMaxwc):
    dataFin = pd.read_csv('./../global/dataCmd.csv', ';')
    wantedDataLine = dataFin.loc[(dataFin["Parameter"] == data) & (dataFin["Place"] == zone)]
    data = wantedDataLine.iloc[-1]["variable"]  # we find the data name in the dataset
    listValue = []
    ldate = []
    dataTabl = []
    for r, d, f in os.walk(path):
        listValue = np.zeros(((len(f)-10)*361,latitudeMaxwc-latitudeMinwc,longitudeMaxwc-longitudeMinwc))
        for i in range(len(f)-10):
            fn = path + f[i]
            print(fn)
            # we read the file
            ds = nc.Dataset(fn)
            # we get the date
            #ldate += [ds['time'][0]]
            # we get the data
            listValue[i*361:(i+1)*361] = ds[data][:, latitudeMinwc:latitudeMaxwc, longitudeMinwc:longitudeMaxwc]
    #listValue, ldate = sortDateList(listValue, ldate)
    #dataTabl += [listValue.tolist()]
    return listValue

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

    eastNan = np.where(np.isnan(dataNO3[0][:, 1:])) #we have a Nan at the east

    downNani_1j = np.where(np.isnan(dataNO3[0][:-1]))
    downNan = (downNani_1j[0] + 1, downNani_1j[1])

    upNan = np.where(np.isnan(dataNO3[0][1:]))

    up2Nan = np.where(np.isnan(dataNO3[0][2:]))

    down2Nani_1j = np.where(np.isnan(dataNO3[0][:-2]))
    down2Nan = (down2Nani_1j[0] + 2, down2Nani_1j[1])

    east2Nan = np.where(np.isnan(dataNO3[0][:,2:]))

    west2Nani_1j = np.where(np.isnan(dataNO3[0][:,:-2]))
    west2Nan = (west2Nani_1j[0], west2Nani_1j[1]+2)

    upEastNan = np.where(np.isnan(dataNO3[0][1:,1:]))

    downWestNani_1j_1 = np.where(np.isnan(dataNO3[0][:-1,:-1]))
    downWestNan= (downWestNani_1j_1[0] + 1, downWestNani_1j_1[1]+1)

    upWestNanj_1 = np.where(np.isnan(dataNO3[0][1:,:-1]))
    upWestNan = (upWestNanj_1[0], upWestNanj_1[1] + 1)

    downEastNani_1 = np.where(np.isnan(dataNO3[0][:-1, 1:]))
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

    return westNanList, eastNanList, upNanList, downNanList, west2NanList, east2NanList, up2NanList, down2NanList,\
           downEastNanList, upWestNanList, downWestNanList, upEastNanList

def createMatE(nbrx, nbry,centuredEwc,eastNanList, westNanList, east2NanList, west2NanList, alpha1, alpha2):
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


def createMatN(nbrx, nbry, centuredNwc, downNanList, upNanList, down2NanList, up2NanList, alpha1, alpha2):
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

    termA[down2NanList] = 0
    termB[down2NanList] = 0

    termA[up2NanList] = 0
    termC[up2NanList] = 0
    data = np.zeros((5, nbrx * nbry))
    data[0] = termA
    data[1][:-nbry] = termB[nbry:]
    data[2][nbry:] = termC[:-nbry]
    data[3][:-2*nbry] = termD[2*nbry:]
    data[4][2*nbry:] = termE[:-2*nbry]
    Mat = dia_matrix((data, offset), shape=(nbrx * nbry, nbrx * nbry))

    return Mat

def quickest(dyMeter, dxlist, dt,discr, centuredEwc, centuredNwc, dataNO3, Ks, firstday):
    nbrStep = int(len(centuredNwc) * discr)
    C = np.zeros(np.shape(dataNO3[0]))
    (nbrx, nbry) = np.shape(dataNO3[0])
    mask = dataNO3.mask
    maskpos2D = np.where(mask[0] == True)
    init = time.time()
    westNanList, eastNanList, upNanList, downNanList, west2NanList, east2NanList, up2NanList, down2NanList, \
    downEastNanList, upWestNanList, downWestNanList, upEastNanList = findNan(dataNO3,nbry)
    listCprim = [dataNO3[0][CPlat,CPlon]]
    listC = [dataNO3[0][CPlat,CPlon]]
    daylist = [firstday]
    clist = [0]

    indx = np.random.choice(range(nbrx),int(nbrx*nbry/2))
    indy = np.random.choice(range(nbry),int(nbrx*nbry/2))
    posFarm = np.ones((nbrx, nbry))
    posFarm[indx,indy] = 0

    fig, ax = plt.subplots()
    CpnIm = copy.deepcopy(posFarm)
    CpnIm[maskpos2D] = np.nan
    plt.imshow(CpnIm)
    plt.colorbar()
    ax.invert_yaxis()
    ax.set_title('Farm position')
    plt.show()

    for k in range(nbrStep):
        dayNbr =  k// int(discr)
        daylist += [firstday+ relativedelta(minutes=int(k*24*60/discr))]
        hNbr = k // int(discr)
        if hNbr>len(centuredEwc)-1:
            break
        '''if k == 0:
            fig, ax = plt.subplots()
            B=np.ones(nbrx * nbry)
            B[downNanList] = 0
            B[upNanList] = 0
            B[eastNanList] = 0
            B[westNanList] = 0
            B[downEastNanList] = 0
            B[upWestNanList] = 0
            B[downWestNanList] = 0
            B[upEastNanList] = 0
            B = B.reshape(nbrx,nbry)
            B[maskpos2D] = np.nan
            plt.scatter([CPlon], [CPlat], c='red')
            plt.imshow(-B)
            plt.colorbar()
            ax.invert_yaxis()
            plt.show()'''
        C[maskpos2D] = 0
        CFL, cx, cy = giveCFL(dxlist, dyMeter, dt, centuredEwc[hNbr], centuredNwc[hNbr],nbrx, nbry)
        #print('CFL max: ', max(CFL))
        alpha1 = (1 / 6) * (1 - CFL) * (2 - CFL)
        alpha2 = (1 / 6) * (1 - CFL) * (1 + CFL)
        matE = createMatE(nbrx, nbry, centuredEwc[hNbr], eastNanList, westNanList, east2NanList, west2NanList, alpha1, alpha2)
        matN = createMatN(nbrx, nbry, centuredNwc[hNbr], downNanList, upNanList, down2NanList, up2NanList, alpha1, alpha2)

        '''B = np.zeros((nbrx , nbry))
        B[158,9] = 1*dt
        B[145,34] = 1*dt
        B[130,50] = 1*dt
        B[110,81] = 1*dt
        B[72,90] = 1*dt
        B[44,86] = 1*dt
        B[27,71] = 1*dt
        B[31,16] = 1*dt
        B = B.reshape(nbrx * nbry)'''
        B = -dt * 0.5 * (C + dataNO3[dayNbr])
        B[indx,indy] = 0
        B = B.reshape(nbrx * nbry)
        '''B[downNanList] = 0
        B[upNanList] = 0
        B[eastNanList] = 0
        B[westNanList] = 0

        B[downEastNanList] = 0
        B[upWestNanList] = 0
        B[downWestNanList] = 0
        B[upEastNanList] = 0'''

        oldCline = C.reshape(nbrx * nbry)
        Cline = oldCline - cy*matN.dot(oldCline) - cx*matE.dot(oldCline)+(dt/dyMeter)*Ks*matN.dot(matN.dot(oldCline))+(dt/dxlist.reshape(nbrx * nbry))*Ks*matE.dot(matN.dot(oldCline)) + B

        '''Cline[downNanList] = oldCline[downNanList]  + B[downNanList]
        Cline[upNanList] = oldCline[upNanList] + B[upNanList]
        Cline[eastNanList] = oldCline[eastNanList] + B[eastNanList]
        Cline[westNanList] = oldCline[westNanList] + B[westNanList]

        Cline[downEastNanList] = oldCline[downEastNanList] + B[downEastNanList]
        Cline[upWestNanList] = oldCline[upWestNanList] + B[upWestNanList]
        Cline[downWestNanList] = oldCline[downWestNanList] + B[downWestNanList]
        Cline[upEastNanList] = oldCline[upEastNanList] + B[upEastNanList]'''
        C = np.minimum(Cline.reshape(nbrx, nbry),0)
        #C = Cline.reshape(nbrx, nbry)

        listCprim += [C[CPlat,CPlon]+dataNO3[dayNbr][CPlat,CPlon]]
        listC += [dataNO3[dayNbr][CPlat,CPlon]]
        clist +=[C[CPlat,CPlon]]
        if k % int(30 * discr) == 0:
            print(k / discr)
            print(time.time() - init)

    fig, ax = plt.subplots()
    CpnIm = copy.deepcopy(C)
    CpnIm[maskpos2D] = np.nan
    plt.imshow(CpnIm)
    plt.colorbar()
    plt.clim(-10, 0)
    ax.invert_yaxis()
    ax.set_title('c')
    plt.show()

    fig, ax = plt.subplots()
    CpnIm = copy.deepcopy(C + dataNO3[dayNbr-1])
    CpnIm[maskpos2D] = np.nan
    plt.imshow(CpnIm)
    plt.colorbar()
    plt.clim(0, 20)
    # plt.scatter([9,34,50,81,90,86,71,16],[158,145,130,110,72,44,27,31], c='red')
    ax.invert_yaxis()
    ax.set_title("C'")
    plt.show()

    #print(listCprim)
    fig, ax = plt.subplots()
    plt.plot(daylist,listC, label='C')
    #plt.plot(daylist[:len(listCprimHourly)],listCprimHourly, label='C prime Hourly')
    ax.plot(daylist, listCprim, label="C' daily")
    ax.legend()
    ax.set_xlabel('date')
    ax.set_ylabel('Nitrate')
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

    discr = 288
    dt = 1 / discr

    Ks = 1e-3
    Ks *= 60 * 60 * 24

    CPlat,CPlon = 153, 56

    dataCmdpath = 'I:/work-he/apps/safi/data/IBI/dataCmd.csv'
    dataNO3, ds = getwantedMergeData('Nitrate', depth, dataCmdpath)

    fileU = 'D:/Profils/mjaouen/Documents/alternance/EASME/data/cmems_mod_ibi_phy_u0.nc'
    dataBaseEwc = nc.Dataset(fileU)

    fileV = 'D:/Profils/mjaouen/Documents/alternance/EASME/data/cmems_mod_ibi_phy_v0.nc'
    dataBaseNwc = nc.Dataset(fileV)

    fileV = 'D:/Profils/mjaouen/Documents/alternance/EASME/data/cmems_mod_ibi_phy_v0.nc'
    olddataBaseNwc = nc.Dataset(fileV)

    dataFin = pd.read_csv('./../global/dataCmd.csv', ';')

    ewcDataLine = dataFin.loc[(dataFin["Parameter"] == 'eastward_Water_current') & (dataFin["Place"] == zone)]
    ewcdataName = ewcDataLine.iloc[-1]["variable"]  # we find the data name in the dataset

    nwcDataLine = dataFin.loc[(dataFin["Parameter"] == 'northward_Water_current') & (dataFin["Place"] == zone)]
    nwcdataName = nwcDataLine.iloc[-1]["variable"]  # we find the data name in the dataset

    longitudeMin, latitudeMin = givecoor(olddataBaseNwc, lonmin, latmin, 'northward_Water_current', dataFin)  # we get the indices of the wanted position
    longitudeMax, latitudeMax = givecoor(olddataBaseNwc, lonmax, latmax, 'northward_Water_current', dataFin)  # we get the indices of the wanted position

    '''longitudeMinwc, latitudeMinwc = givecoor(dataBaseNwc, lonmin, latmin, 'northward_Water_current',
                                         dataFin)  # we get the indices of the wanted position
    longitudeMaxwc, latitudeMaxwc = givecoor(dataBaseNwc, lonmax, latmax, 'northward_Water_current',
                                         dataFin)  # we get the indices of the wanted position'''
    #we get the hourly data
    '''dataEwc = getData('D:/Profils/mjaouen/Documents/alternance/EASME/data/eastward_Water_current/', 'eastward_Water_current', zone, depth, latitudeMin, latitudeMax, longitudeMin, longitudeMax)
    dataNwc = getData('D:/Profils/mjaouen/Documents/alternance/EASME/data/northward_Water_current/',
                      'northward_Water_current', zone, depth, latitudeMin, latitudeMax, longitudeMin,
                      longitudeMax)
    print(np.shape(dataEwc))'''

    #print(np.shape(dataNwc))
    dataNO3 = dataNO3[:, latitudeMin:latitudeMax, longitudeMin:longitudeMax]

    #we get daily data
    dataEwc, ds = getwantedMergeData('eastward_Water_current', depth, dataCmdpath)
    dataNwc, ds = getwantedMergeData('northward_Water_current', depth, dataCmdpath)

    dataNwc = dataNwc[:, latitudeMin:latitudeMax, longitudeMin:longitudeMax]
    dataEwc = dataEwc[:, latitudeMin:latitudeMax, longitudeMin:longitudeMax]

    dataNwc *= 60 * 60 * 24
    dataEwc *= 60 * 60 * 24

    dataEwc[np.where(np.abs(dataEwc) >1e8)]=0
    dataNwc[np.where(np.abs(dataNwc) >1e8)] = 0

    print(np.shape(np.array(dataBaseEwc[ewcDataLine.iloc[-1]["latName"]][latitudeMin:latitudeMax])))
    print(np.shape(dataEwc[0])[1],np.shape(dataEwc[0])[0])

    latRef = np.ones((np.shape(dataEwc[0])[1],np.shape(dataEwc[0])[0]))*np.array(dataBaseEwc[ewcDataLine.iloc[-1]["latName"]][latitudeMin:latitudeMax])

    centuredEwc = u2d_cgrid_cur(dataEwc)
    centuredNwc = v2d_cgrid_cur(dataNwc)

    dxlist, dyMeter = degrees_to_meters(dxdeg, dydeg, latRef)

    dxlist = dxlist.T
    dylist = dyMeter*np.ones(np.shape(dxlist))

    firstday = datetime.datetime.strptime('2020-01-18', '%Y-%m-%d')

    ''' interdataNO3 = np.zeros((np.shape(dataNO3)[0]*24,np.shape(dataNO3)[1],np.shape(dataNO3)[2]))
    for i in range(np.shape(dataNO3)[1]):
        for j in range(np.shape(dataNO3)[2]):
            inter = interpolate.splrep(np.linspace(0,1,len(dataNO3)), dataNO3[:,i,j])
            interdataNO3[:,i,j] = interpolate.splev(np.linspace(0,1,len(dataNO3)*24), inter)
    maskNO3 = np.zeros(np.shape(interdataNO3))
    maskNO3[np.where(interdataNO3>1e8)] = 1
    interdataNO3 = ma.masked_array(interdataNO3, mask=maskNO3)'''

    quickest(dyMeter, dxlist, dt, discr, centuredEwc, centuredNwc, dataNO3, Ks, firstday)

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
