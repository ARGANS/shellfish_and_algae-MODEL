import ftplib

from owslib.wcs import WebCoverageService
from owslib.wfs import WebFeatureService
import cdsapi
import os
import numpy as np
import pandas as pd
import datetime

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

#exctract the data from marine copernicus
def getdataFromMarineCopernicus(dataInfo, dateBeginning, dateEnd, outputDirectory, outputFile, deepthmin=10, deepthmax=15):
    if (dataInfo[18] == '') or (dataInfo[18] == 'NULL'):
        os.system('python -m motuclient --motu ' + dataInfo[3] +
                  ' --service-id ' + dataInfo[4] + ' --product-id ' + dataInfo[5] +
                  ' --longitude-min ' + dataInfo[6] + ' --longitude-max ' + dataInfo[
                      7] + ' --latitude-min ' + dataInfo[8] + ' --latitude-max ' + dataInfo[
                      9] + ' --date-min ' + dateBeginning +
                  ' --date-max ' + dateEnd + ' --variable ' + dataInfo[15] +
                  ' --out-dir ' + outputDirectory + ' --out-name ' + outputFile + ' --user mjaouen --pwd Azerty123456 ')
    else:

        os.system('python -m motuclient --motu ' + dataInfo[3] +
                  ' --service-id ' + dataInfo[4] + ' --product-id ' + dataInfo[5] +
                  ' --longitude-min ' + dataInfo[6] + ' --longitude-max ' + dataInfo[
                      7] + ' --latitude-min ' +
                  dataInfo[8] + ' --latitude-max ' + dataInfo[9] + ' --date-min ' + dateBeginning +
                  ' --date-max ' + dateEnd + ' --depth-min ' + str(deepthmin) + '  --depth-max ' +
                  str(deepthmax) + '  --variable ' +
                  dataInfo[15] +
                  ' --out-dir ' + outputDirectory + ' --out-name ' + outputFile + ' --user mjaouen --pwd Azerty123456 ')

#give the file complete name depending of the filetype
def giveFile(filename,filetype):
    if filetype == 'GeoTiFF':
        return filename+'.tiff'
    elif filetype == 'NetCDF':
        return filename+'.nc'
    else:
        return filename+'.'+filetype

def getdataFromFtp(dataFin):
    HOSTNAME = dataFin[16]
    USERNAME = dataFin[4]
    PASSWORD = dataFin[5]

    # Connect FTP Server
    ftp_server = ftplib.FTP(HOSTNAME, USERNAME, PASSWORD)

    # force UTF-8 encoding
    ftp_server.encoding = "utf-8"
    ftp_server.cwd('/'+dataFin[6])
    filelist = ftp_server.nlst()

    ftp_server.dir()
    path = 'I:/work-he/apps/safi/data/IBI/par/'
    # os.mkdir(path)
    #we read all the file in the ftp directory
    for filename in filelist:
        print(filename)
        fileDirect = path + filename
        #we download the file
        with open(fileDirect, "wb") as file:
            ftp_server.retrbinary(f"RETR {filename}", file.write)

def getdataWCS(url,layer,requestbbox,file, version,format='GeoTIFF'):
    wcs = WebCoverageService(url, version=version, timeout=320)
    sed = wcs[layer]  # this is necessary to get essential metadata from the layers advertised
    cx, cy = map(int, sed.grid.highlimits)
    bbox = sed.boundingboxes[0]['bbox']
    lx, ly, hx, hy = map(float, bbox)
    resx, resy = (hx - lx) / cx, (hy - ly) / cy

    gc = wcs.getCoverage(identifier=layer,
                         bbox=requestbbox,
                         coverage=sed,
                         format=format,
                         crs=sed.boundingboxes[0]['nativeSrs'], resx=resx, resy=resy)
    f = open(file, 'wb')
    f.write(gc.read())
    f.close()

#for a date this function gives the year the month and the day numbers
def splitDate(date):
    datelist=date.split("-")
    year=datelist[0].split('"')[1]
    month=datelist[1]
    day=datelist[2].split(" ")[0]
    return [year,month,day]

def givedatesForClimatCoper(begDate,endDate):
    years = [ str(yr) for yr in range(int(begDate[0]),int(endDate[0])+1)]
    #if we stay in the same year
    if int(begDate[0])==int(endDate[0]):
        months = [ str(mnt) for mnt in range(int(begDate[1]),int(endDate[1])+1)]
    #if we begin a year and we end the year after
    elif int(begDate[0])+1==int(endDate[0]):
        months = [str(mnt) for mnt in range(int(begDate[1]), 13)]
        months += [str(mnt) for mnt in range(1, int(endDate[1]) + 1)]
    else:
        months = [str(mnt) for mnt in range(1,13)]
        days = [str(day) for day in range(1, 32)]
        return years, months, days
    #if we begin a month and we end the same month (of the same year)
    if int(begDate[1])==int(endDate[1]) and int(begDate[0])==int(endDate[0]):
        days = [ str(mnt) for mnt in range(int(begDate[2]),int(endDate[2])+1)]
    #if we begin a month and we finish the month after (of the same year or the year after if we begin in december)
    elif (int(begDate[1])+1==int(endDate[1]) and int(begDate[0])==int(endDate[0])) \
            or (int(begDate[1])==12 and int(endDate[1])==1 and int(begDate[0])==int(endDate[0])-1) :
        days = [str(day) for day in range(int(begDate[2]), 32)]
        days += [str(day) for day in range(1, int(endDate[2]) + 1)]
    else:
        days = [str(days) for days in range(1,32)]
    return years, months, days


def getData(wantedData, zone, dataFin, deepthmin, deepthmax, dateBeginning, dateEnd,outputDirectory):
    begDate = splitDate(dateBeginning)
    endDate = splitDate(dateEnd)
    # we select the lines that contains the data on the right zone
    wantedDataLine = np.where((dataFin[:, 1] == wantedData) & (dataFin[:, 2] == zone))
    servicetypelist = dataFin[wantedDataLine, 14]
    for j in range(len(wantedDataLine[0])):
        servicetype = servicetypelist[0][j]
        imgNb = wantedDataLine[0][j]
        filename = wantedData + zone + dataFin[imgNb][0] + dateBeginning[1:11]
        outputFile = giveFile(filename, dataFin[imgNb, 20])
        if servicetype == 'marineCopernicus':
            getdataFromMarineCopernicus(dataFin[imgNb], dateBeginning, dateEnd, outputDirectory, outputFile, deepthmin,
                                        deepthmax)

        elif servicetype == 'WCS':
            # define the connection
            url = dataFin[imgNb, 16]

            # define variables
            # requestbbox = (lonOuest, latSud, lonEst, latNord)
            requestbbox = (-4.7, 48.55, -4.50, 48.65)
            # requestbbox = (2.0, 51.5, 5.0, 54.0)
            layer = dataFin[imgNb, 15]
            # get the data
            getdataWCS(url, layer, requestbbox, outputFile, dataFin[imgNb, 4], format=dataFin[imgNb, 20])

        elif servicetype == 'cdsapi':
            c = cdsapi.Client()
            variable = dataFin[imgNb, 6]
            fileformat = dataFin[imgNb, 20]
            prodtype = dataFin[imgNb, 3]
            prodname = dataFin[imgNb, 4]
            time = dataFin[imgNb, 9]
            years, months, days = givedatesForClimatCoper(begDate, endDate)
            c.retrieve(
                prodname,
                {
                    'product_type': prodtype,
                    'format': fileformat,
                    'variable': variable,
                    'year': years,
                    'month': months,
                    'day': days,
                    'time': time,
                },
                outputFile)

        elif servicetype == 'ftp':
            getdataFromFtp(dataFin[imgNb])

def giveDateslist(dateBeginning,dateEnd):
    begList = [dateBeginning]
    endList = []
    gDate = datetime.datetime(int(dateBeginning[1:5]), int(dateBeginning[6:8]), int(dateBeginning[9:11]))
    date = gDate + datetime.timedelta(days = 1)
    datestr = date.strftime('"%Y-%m-%d %H:%M:%S"')
    while datestr != dateEnd:
        begList += [datestr]
        endList += [datestr]
        date = date + datetime.timedelta(days=1)
        datestr = date.strftime('"%Y-%m-%d %H:%M:%S"')
    endList+=[dateEnd]
    return begList, endList


wantedData=['Temperature','Nitrate','Ammonium','Phosphate','Salinity','northward_Water_current','eastward_Water_current']
dateBeginning = '"2020-11-15 00:00:00"'
dateEnd = '"2020-11-22 00:00:00"'
zone='IBI'
deepthmin=0
deepthmax=20
outputDirectory = 'I:/work-he/apps/safi/data/IBI/'

dataFin=readcsv()
datesList=giveDateslist(dateBeginning,dateEnd)
for i in range(len(datesList[0])):
    dateBeg=datesList[0][i]
    dateE=datesList[1][i]
    for dat in wantedData:
        DataOutputDirectory=outputDirectory+dat+'/'
        getData(dat, zone, dataFin, deepthmin, deepthmax, dateBeg, dateE,DataOutputDirectory)