import datetime
import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import pandas as pd
from statsmodels.nonparametric.kernel_regression import KernelReg
import h5py

path = 'I:/work-he/apps/safi/data/IBI/'

#read the csv where we listed the location of the different data
def readcsv(file='./dataCmd.csv'):
    data = pd.read_csv(file, delimiter=',', header=0)
    numpData = data.to_numpy()
    for k in range(len(numpData)):
        numpData[k] = [numpData[k][0].split(';')]
    dataFin = []
    for k in range(len(numpData)):
        dataFin += [numpData[k][0]]
    dataFin = np.array(dataFin)
    return dataFin

def givecoor(path,lonval,latval,dataName):
    for r, d, f in os.walk(path):
        fn = path + f[0]
        ds = nc.Dataset(fn)
        if dataName == 'par':
            lonList = ds['lon'][:]
            latList = ds['lat'][:]
        else:
            lonList = ds['longitude'][:]
            latList = ds['latitude'][:]
        i=0
        loni = lonList[i]
        while i+1<len(lonList) and lonval>loni:
            i+=1
            loni = lonList[i]
        j = 0
        lati = latList[j]
        while j+1 < len(latList) and latval > lati:
            j += 1
            lati = latList[j]
    return i, j

def sortDateList(listValue,ldate):
    sortLval = []
    sortldate = []
    for k in range(len(ldate)):
        i = 0
        while i<len(sortLval):
            if ldate[k]<sortldate[i]:
                sortldate = sortldate[:i]+[ldate[k]]+sortldate[i:]
                sortLval = sortLval[:i] + [listValue[k]] + sortLval[i:]
                i=len(sortLval)+1
            i+=1
        if i == len(sortLval)-1 or i == len(sortLval):
            sortldate += [ldate[k]]
            sortLval += [listValue[k]]
    return np.array(sortLval), np.array(sortldate)

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

lat = 51.583332
lon = -9.916667

depth = 1
dataName = 'Temperature'
zone = 'IBI'
path=path+dataName+"/"
print(path)
dataFin=readcsv()
wantedDataLine = np.where((dataFin[:, 1] == dataName) & (dataFin[:, 2] == zone))
imgNb = wantedDataLine[0][0]
data = dataFin[imgNb][15]
longitude, latitude = givecoor(path,lon,lat,dataName)

listValue = []
ldate = []
for r, d, f in os.walk(path):
    for i in range(len(f)):
        print(f[i])
        fn=path+f[i]
        ds = nc.Dataset(fn,'r', format='NETCDF4')
        if dataName == 'par':
            date = fn[-46:-38]
            ldate += [datetime.datetime(int(date[0:4]), int(date[4:6]), int(date[6:8]))]
            listValue += [ds[data][latitude, longitude]]
        else:
            date = fn[-13:-3]
            ldate += [datetime.datetime(int(date[0:4]), int(date[5:7]), int(date[8:10]))]
            listValue += [ds[data][0, depth, latitude, longitude]]


listValue,ldate = sortDateList(listValue,ldate)

listValue = np.array(listValue)

for i in range(1,len(listValue)):
    if listValue[i]==np.nan and listValue[i-1]!=np.nan and listValue[i+1]!=np.nan:
        listValue = 0.5*(listValue[i-1] + listValue[i+1])

ldateNoNan = ldate[np.logical_not(np.isnan(listValue))]
listnoNAN=listValue[np.logical_not(np.isnan(listValue))]

print(listnoNAN)



valList = np.linspace(1,366,len(listnoNAN))



plt.plot(ldateNoNan,smooth(listnoNAN,10), label='lissage',color='red')
plt.scatter(ldate,listValue, label ='valeurs originales')
plt.legend()

plt.title('{}, lat: {}, lon : {}'.format(dataName,lon,lat))
plt.show()