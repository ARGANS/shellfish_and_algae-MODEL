import ftplib
import math
from pprint import pprint

from owslib.wcs import WebCoverageService
from owslib.wfs import WebFeatureService
import cdsapi
import os
import numpy as np
import pandas as pd
import datetime
from dateutil.relativedelta import *


# read the csv where we listed the location of the different data
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


# exctract the data from marine copernicus
def getdataFromMarineCopernicus(dataInfo, dateBeginning, dateEnd, outputDirectory, outputFile, deepthmin=10,
                                deepthmax=15):
    """
    Since version 1.8.0:
        motuclient -h
    Before version 1.8.0:
        python -m motu-client -h
    Newest version (?):
        python -m motuclient (part of the code prosposed by CMEMS)
    """
    if (np.isnan(dataInfo["depth-min"])):
        os.system(f'python -m motuclient'
                  f' --motu {dataInfo["motu"]}'
                  f' --service-id {dataInfo["service-id"]}'
                  f' --product-id {dataInfo["product-id"]}'
                  f' --longitude-min {dataInfo["longitude-min"]}'
                  f' --longitude-max {dataInfo["longitude-max"]}'
                  f' --latitude-min {dataInfo["latitude-min"]}'
                  f' --latitude-max {dataInfo["latitude-max"]}'
                  f' --date-min {dateBeginning}'
                  f' --date-max {dateEnd}'
                  f' --variable {dataInfo["variable"]}'
                  f' --out-dir {outputDirectory}'
                  f' --out-name {outputFile}'
                  f' --user mjaouen --pwd Azerty123456'
                  )
    else:
        deepthmin = str(deepthmin)
        deepthmax = str(deepthmax)
        os.system(f'python -m motuclient'
                  f' --motu {dataInfo["motu"]}'
                  f' --service-id {dataInfo["service-id"]}'
                  f' --product-id {dataInfo["product-id"]}'
                  f' --longitude-min {dataInfo["longitude-min"]}'
                  f' --longitude-max {dataInfo["longitude-max"]}'
                  f' --latitude-min {dataInfo["latitude-min"]}'
                  f' --latitude-max {dataInfo["latitude-max"]}'
                  f' --date-min {dateBeginning}'
                  f' --date-max {dateEnd}'
                  f' --depth-min {deepthmin}'
                  f' --depth-max {deepthmax}'
                  f' --variable {dataInfo["variable"]}'
                  f' --out-dir {outputDirectory}'
                  f' --out-name {outputFile}'
                  f' --user "mjaouen" --pwd "Azerty123456"'
                  )


# give the file complete name depending of the filetype
def giveFile(filename, filetype):
    if filetype == 'GeoTiFF':
        return filename + '.tiff'
    elif filetype == 'NetCDF':
        return filename + '.nc'
    else:
        return filename + '.' + filetype


def getdataFromFtp(dataFin, outputDirectory):
    HOSTNAME = dataFin["lien"]
    USERNAME = dataFin["motu"]
    PASSWORD = dataFin["service-id"]
    # Connect FTP Server
    ftp_server = ftplib.FTP(HOSTNAME, USERNAME, PASSWORD)

    # force UTF-8 encoding
    ftp_server.encoding = "utf-8"
    try:
        print('ok')
        print(f'FTP {"/" + dataFin["product-id"]}')
        ftp_server.cwd('/' + dataFin["product-id"])
        filelist = ftp_server.nlst()

        ftp_server.dir()
        # os.mkdir(path)
        # we read all the file in the ftp directory
        for filename in filelist:
            print(filename)
            fileDirect = outputDirectory + filename
            # we download the file
            with open(fileDirect, "wb") as file:
                ftp_server.retrbinary(f"RETR {filename}", file.write)
    except Exception as exc:
        print(f'Exception {exc} // {exc.__class__}')


def getdataWCS(url, layer, requestbbox, file, version, format='GeoTIFF'):
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


# for a date this function gives the year the month and the day numbers
def splitDate(date):
    return [date.year, date.month, date.day]


def givedatesForClimatCoper(begDate, endDate):
    years = [str(yr) for yr in range(int(begDate[0]), int(endDate[0]) + 1)]
    # if we stay in the same year
    if int(begDate[0]) == int(endDate[0]):
        months = [str(mnt) for mnt in range(int(begDate[1]), int(endDate[1]) + 1)]
    # if we begin a year and we end the year after
    elif int(begDate[0]) + 1 == int(endDate[0]):
        months = [str(mnt) for mnt in range(int(begDate[1]), 13)]
        months += [str(mnt) for mnt in range(1, int(endDate[1]) + 1)]
    else:
        months = [str(mnt) for mnt in range(1, 13)]
        days = [str(day) for day in range(1, 32)]
        return years, months, days
    # if we begin a month and we end the same month (of the same year)
    if int(begDate[1]) == int(endDate[1]) and int(begDate[0]) == int(endDate[0]):
        days = [str(mnt) for mnt in range(int(begDate[2]), int(endDate[2]) + 1)]
    # if we begin a month and we finish the month after (of the same year or the year after if we begin in december)
    elif (int(begDate[1]) + 1 == int(endDate[1]) and int(begDate[0]) == int(endDate[0])) \
            or (int(begDate[1]) == 12 and int(endDate[1]) == 1 and int(begDate[0]) == int(endDate[0]) - 1):
        days = [str(day) for day in range(int(begDate[2]), 32)]
        days += [str(day) for day in range(1, int(endDate[2]) + 1)]
    else:
        days = [str(days) for days in range(1, 32)]
    return years, months, days


def validate(value, cls):
    if isinstance(value, cls):
        return value
    return None


def getData(wantedData, zone, dataFin, deepthmin, deepthmax, outputDirectory, dateBeginning=None, dateEnd=None,
            frequency='daily', type='model'):
    # we select the lines that contains the data on the right zone
    wantedDataLine = dataFin.loc[
        (dataFin["Parameter"] == wantedData) & (dataFin["Place"] == zone) & (dataFin["daily"] == frequency) & (dataFin["type"] == type)]
    for j in wantedDataLine.index.values:
        servicetype = dataFin.iloc[j]["source"]
        if wantedDataLine.iloc[0]["daily"] == 'daily':
            begDate = splitDate(dateBeginning)
            endDate = splitDate(dateEnd)
            filename = wantedData + zone + (validate(dataFin.iloc[j].get("type"), str) or '') + dataFin.iloc[j][
                "fileType"] + dateBeginning.strftime("%Y-%m-%d")+ 'to' + dateEnd.strftime("%Y-%m-%d")
        elif wantedDataLine.iloc[0]["daily"] == 'monthly':
            begDate = splitDate(dateBeginning)
            endDate = splitDate(dateEnd)
            filename = wantedData + zone + (validate(dataFin.iloc[j].get("type"), str) or '') + dataFin.iloc[j][
                "fileType"] + dateBeginning.strftime("%Y-%m") + 'to' + dateEnd.strftime("%Y-%m")
        elif wantedDataLine.iloc[0]["daily"] == 'hourly':
            begDate = splitDate(dateBeginning)
            endDate = splitDate(dateEnd)
            filename = wantedData + zone + (validate(dataFin.iloc[j].get("type"), str) or '') + dataFin.iloc[j][
                "fileType"] + dateBeginning.strftime("%Y-%m-%d%H%M%S") + 'to' + dateEnd.strftime("%Y-%m-%d%H%M%S")
        else:
            filename = wantedData + zone + dataFin.iloc[j]["fileType"]
        outputFile = giveFile(filename, dataFin.iloc[j]["fileType"])
        print(servicetype)
        if servicetype == 'marineCopernicus':
            getdataFromMarineCopernicus(dataFin.iloc[j], dateBeginning.strftime('"%Y-%m-%d %H:%M:%S"'),
                                        dateEnd.strftime('"%Y-%m-%d %H:%M:%S"'), outputDirectory,
                                        outputFile, deepthmin, deepthmax)

        elif servicetype == 'WCS':
            # define the connection
            url = dataFin.iloc[j]["lien"]
            # define variables
            # requestbbox = (lonOuest, latSud, lonEst, latNord)
            requestbbox = (-4.7, 48.55, -4.50, 48.65)
            # requestbbox = (2.0, 51.5, 5.0, 54.0)
            layer = dataFin.iloc[j]["variable"]
            outputFileAdress = outputDirectory + outputFile
            # get the data
            getdataWCS(url, layer, requestbbox, outputFileAdress, dataFin.iloc[j]["service-id"],
                       dataFin.iloc[j]["fileType"])

        elif servicetype == 'cdsapi':
            c = cdsapi.Client()
            variable = dataFin.iloc[j]["product-id"]
            fileformat = dataFin.iloc[j]["fileType"]
            prodtype = dataFin.iloc[j]["type"]
            prodname = dataFin.iloc[j]["service-id"]
            time = dataFin.iloc[j]["time"]
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
            getdataFromFtp(dataFin.iloc[j], outputDirectory)


def giveDateslist(dateBeginning, dateEnd, frequency, timestep = None):
    if timestep:
        ndays = (dateEnd - dateBeginning).days
        begList = [(dateBeginning + datetime.timedelta(days=timestep)) for i in range(int(ndays/timestep))]
        endList = [(dateBeginning + datetime.timedelta(days=timestep)) for i in range(1, int(ndays/timestep) + 1)]
    elif frequency == 'monthly':
        begList = [datetime.datetime.strptime(dateBeginning, '%Y-%m-%d %H:%M:%S')]
        endList = [datetime.datetime.strptime(dateEnd, '%Y-%m-%d %H:%M:%S')]
    elif frequency == 'hourly':
        datetimeBeginning = datetime.datetime.strptime(dateBeginning, '%Y-%m-%d %H:%M:%S')
        datetimeEnd = datetime.datetime.strptime(dateEnd, '%Y-%m-%d %H:%M:%S')

        nmonth = ((datetimeEnd.year - datetimeBeginning.year) * 12) + datetimeEnd.month - datetimeBeginning.month

        begList = [(datetimeBeginning + relativedelta(days=15*i)) for i in range(int(nmonth*2))]
        endList = [(datetimeBeginning + relativedelta(days=15*i)) for i in range(1, int(nmonth*2) + 1)]
    else:
        datetimeBeginning = datetime.datetime.strptime(dateBeginning, '%Y-%m-%d %H:%M:%S')
        datetimeEnd = datetime.datetime.strptime(dateEnd, '%Y-%m-%d %H:%M:%S')

        nmonth = ((datetimeEnd.year - datetimeBeginning.year) * 12) + datetimeEnd.month - datetimeBeginning.month

        begList = [(datetimeBeginning + relativedelta(months=i)) for i in range(nmonth)]
        endList = [(datetimeBeginning + relativedelta(months=i)) for i in range(1, nmonth + 1)]

    return begList, endList