from general import readcsv, giveDateslist, getData
import pandas as pd
import os
import pprint

wantedData=['par']
# TODO use a environment variables
dateBeginning = '2020-11-15 00:00:00'
dateEnd = '2020-11-22 00:00:00'
zone='IBI'
deepthmin=0
deepthmax=20
outputDirectory = os.getenv('OUTPUT_DIR')

dataFin=pd.read_csv('./dataCmd.csv',';')
# print('dataFin')
# pprint.pprint(dataFin)

datesList=giveDateslist(dateBeginning,dateEnd)
# print('datesList')
# pprint.pprint(datesList)

for dat in wantedData:
    DataLine = dataFin.loc[dataFin["Parameter"] == dat]
    if DataLine.iloc[0]["daily"] == 1:
        for (dateBeg, dateE) in zip(datesList[0], datesList[1]):
            DataOutputDirectory = outputDirectory + dat + '/'
            getData(dat, zone, dataFin, deepthmin, deepthmax,  DataOutputDirectory, dateBeg, dateE)
    else:
        DataOutputDirectory = outputDirectory + dat + '/'
        getData(dat, zone, dataFin, deepthmin, deepthmax, DataOutputDirectory)
