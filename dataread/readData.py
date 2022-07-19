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
def readcsv(file='./../global/dataCmd.csv'):
    data = pd.read_csv(file, delimiter=',', header=0)
    numpData = data.to_numpy()
    for k in range(len(numpData)):
        numpData[k] = [numpData[k][0].split(';')]
    dataFin = []
    for k in range(len(numpData)):
        dataFin += [numpData[k][0]]
    dataFin = np.array(dataFin)
    return dataFin

#give the indices coresponding to lonval, and latval in the list of coordinates
def givecoor(path,lonval,latval,dataName,dataFin):
    for r, d, f in os.walk(path):
        #we take a data file
        fn = path + f[0]
        ds = nc.Dataset(fn)
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

#sort the list of data from the older to the newer
def sortDateList(listValue,ldate):
    sortLval = [] #we define the sorted data list
    sortldate = [] #we define the sorted date list
    #we read the list we have to sort
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

#smooth the data, with a box_pts size convolutional matrix
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

lat = 51.583332
lon = -9.916667

depth = 1
dataName = 'par'
zone = 'IBI'

#we construct the path to the data
path=path+dataName+"/"
print(path)

#we search after the name of the data we want in the dataset

dataFin=pd.read_csv('./../global/dataCmd.csv',';')
wantedDataLine = dataFin.loc[(dataFin["Parameter"] == dataName) & (dataFin["Place"] == zone)]
data = wantedDataLine.iloc[-1]["variable"]
print(data)
longitude, latitude = givecoor(path,lon,lat,dataName, dataFin)

listValue = []
ldate = []
#we read each file in the folder
for r, d, f in os.walk(path):
    for i in range(len(f)):
        print(f[i])
        fn=path+f[i]
        #we read the file
        ds = nc.Dataset(fn,'r', format='NETCDF4')
        # we look where does the data came from
        if dataName == 'par':
            date = fn[-46:-38]
            ldate += [datetime.datetime(int(date[0:4]), int(date[4:6]), int(date[6:8]))] #we get the date
            listValue += [ds[data][latitude, longitude]] #we get the data
        else:
            date = fn[-13:-3]
            ldate += [datetime.datetime(int(date[0:4]), int(date[5:7]), int(date[8:10]))] #we get the date
            listValue += [ds[data][0, depth, latitude, longitude]] #we get the data


listValue,ldate = sortDateList(listValue,ldate)

listValue = np.array(listValue)

#we attribute the mean of its neighbors the points without data
for i in range(1,len(listValue)):
    if listValue[i]==np.nan and listValue[i-1]!=np.nan and listValue[i+1]!=np.nan:
        listValue = 0.5*(listValue[i-1] + listValue[i+1])

ldateNoNan = ldate[np.logical_not(np.isnan(listValue))]
listnoNAN=listValue[np.logical_not(np.isnan(listValue))]



valList = np.linspace(1,366,len(listnoNAN))



plt.plot(ldateNoNan,smooth(listnoNAN,10), label='lissage',color='red')
plt.scatter(ldate,listValue, label ='valeurs originales')
plt.legend()

plt.title('{}, lat: {}, lon : {}'.format(dataName,lon,lat))
plt.show()