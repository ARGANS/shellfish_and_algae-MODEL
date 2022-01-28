import datetime
import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import pandas as pd

#read the csv where we listed the location of the different data
def readcsv(file='I:/work-he/apps/safi/data/IBI/dataCmd.csv'):
    data = pd.read_csv(file, delimiter=',', header=0)
    numpData = data.to_numpy()
    for k in range(len(numpData)):
        numpData[k] = [numpData[k][0].split(';')]
    dataFin = []
    for k in range(len(numpData)):
        dataFin += [numpData[k][0]]
    dataFin = np.array(dataFin)
    return dataFin

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

#compute the richardson number
def giveRichardsonNbr(dS,dT,dU,dV,dz,alpha=7.5e-5,beta=7.7e-4):
    M2 = dU**2 + dV**2
    g=9,80665
    N2 = (g/dz)*(alpha *dT -beta*dS)
    return M2/N2

#smooth the data, with a box_pts size convolutional matrix
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

#give the indices coresponding to lonval, and latval in the list of coordinates
def givecoor(path,lonval,latval,dataName):
    for r, d, f in os.walk(path):
        # we take a data file
        fn = path + f[0]
        ds = nc.Dataset(fn)
        # we look where does the data came from
        if dataName == 'par':
            # we get the longitude and latitude list
            lonList = ds['lon'][:]
            latList = ds['lat'][:]
        else:
            # we get the longitude and latitude list
            lonList = ds['longitude'][:]
            latList = ds['latitude'][:]
        i=0
        loni = lonList[i]
        # we browse the data until we find a coordiate bigger than the wanted coordiante
        while i<len(lonList) and lonval>loni:
            i+=1
            loni = lonList[i]
        j = 0
        lati = latList[j]
        while j+1 < len(latList) and latval > lati:
            j += 1
            lati = latList[j]
    return i, j

'''latitude = 771
longitude = 561'''

lat = 51.587433
lon = -9.897116


mainpath = 'I:/work-he/apps/safi/data/IBI/'
dataName = 'Temperature'
zone = 'IBI'
dataTabl=[]

path=mainpath+dataName+"/"
dataFin=readcsv()
wantedDataLine = np.where((dataFin[:, 1] == dataName) & (dataFin[:, 2] == zone))
imgNb = wantedDataLine[0][0]
data = dataFin[imgNb][15]
ldate = []

#we get the
for r, d, f in os.walk(path):
        fn=path+f[0]
        ds = nc.Dataset(fn)
        z = ds['depth'][:]
print(z)

# for each wanted data
for dat in ['Salinity', 'Nitrate']:
    path=mainpath+dat+"/"
    print(path)
    dataFin=readcsv()
    wantedDataLine = np.where((dataFin[:, 1] == dat) & (dataFin[:, 2] == zone))
    imgNb = wantedDataLine[0][0]
    data = dataFin[imgNb][15] #we find the data name in the dataset
    longitude, latitude = givecoor(path, lon, lat, dat) #we get the indices of the wanted position
    listValue = []
    ldate = []
    for r, d, f in os.walk(path):
        for i in range(len(f)):
            fn=path+f[i]
            print(fn)
            # we read the file
            ds = nc.Dataset(fn)
            #we get the date
            date = fn[-13:-3]
            ldate += [datetime.datetime(int(date[0:4]), int(date[5:7]), int(date[8:10]))]
            #we get the data
            listValue += [[ds[data][0,0,latitude,longitude],ds[data][0,4,latitude,longitude],ds[data][0,8,latitude,longitude],ds[data][0,11,latitude,longitude]]]
    listValue,ldate = sortDateList(listValue,ldate)
    dataTabl+=[listValue.tolist()]

dataTabl = np.array(dataTabl)

#nR=giveRichardsonNbr(dataTabl[1],dataTabl[0],dataTabl[2],dataTabl[3],dz)

#we plot the results
fig , ax  = plt.subplots()
ax.plot(ldate, smooth(dataTabl[0][:,0],10),  label='0.5 m')
ax.plot( ldate, smooth(dataTabl[0][:,1],10), label='5.1 m')
ax.plot( ldate, smooth(dataTabl[0][:,2],10), label='11.4 m')
ax.plot( ldate, smooth(dataTabl[0][:,3],10), label='18.5 m')
fig2 , ax2  = plt.subplots()
ax2.plot( ldate, smooth(dataTabl[1][:,0],10), label='0.5 m')
ax2.plot( ldate, smooth(dataTabl[1][:,1],10), label='5.1 m')
ax2.plot( ldate, smooth(dataTabl[1][:,2],10), label='11.4 m')
ax2.plot(ldate, smooth(dataTabl[1][:,3],10), label='18.5 m')


ax.legend()
ax.set_xlabel('date')
ax.set_ylabel('PSU')
ax.set_title("{} lat: {}  lon: {}".format(lat,lon))
ax2.legend()
ax2.set_xlabel('date')
ax2.set_ylabel('mmol m-3')
ax2.set_title("{} lat: {}  lon: {}".format(lat,lon))
plt.show()