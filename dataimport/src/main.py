from general import readcsv, giveDateslist, getData
import pandas as pd
import os

wantedData=['Temperature']
dateBeginning = '2021-01-01 00:00:00'
dateEnd = '2022-01-01 00:00:00'
zone='IBI'
deepthmin=0
deepthmax=20
outputDirectory = os.getenv('OUTPUT_DIR')
frequency = 2

#outputDirectory = 'I:/work-he/apps/safi/data/IBI/'
outputDirectory = '/media/share/data/IBI/'

dataFin=pd.read_csv('./dataCmd.csv',';')

datesList=giveDateslist(dateBeginning,dateEnd,frequency)

for dat in wantedData:
    DataLine = dataFin.loc[dataFin["Parameter"] == dat]
    if DataLine.iloc[0]["daily"] > 0:
        DataOutputDirectory = outputDirectory + dat + '/'
        getData(dat, zone, dataFin, deepthmin, deepthmax,  DataOutputDirectory, datesList[0], datesList[1],frequency)
    else:
        DataOutputDirectory = outputDirectory + dat + '/'
        getData(dat, zone, dataFin, deepthmin, deepthmax, DataOutputDirectory)
