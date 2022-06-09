from general import readcsv, giveDateslist, getData
import pandas as pd
import os

wantedData=['Nitrate', 'Ammonium', 'Temperature', 'northward_Water_current', 'eastward_Water_current']
dateBeginning = '2021-01-01 00:00:00'
dateEnd = '2022-01-01 00:00:00'
zone='BS'
deepthmin=0
deepthmax=20
outputDirectory = os.getenv('OUTPUT_DIR')
frequency = 'daily'

#outputDirectory = 'I:/work-he/apps/safi/data/BS/'
outputDirectory = '/media/share/data/IBI/'

dataFin=pd.read_csv('./dataCmd.csv',';')

datesList=giveDateslist(dateBeginning,dateEnd,frequency)


for dat in wantedData:
    DataLine = dataFin.loc[dataFin["Parameter"] == dat]
    if DataLine.iloc[0]["daily"] != 'permanent':
        for (dateBeg, dateE) in zip(datesList[0], datesList[1]):
            DataOutputDirectory = outputDirectory + dat + '/'
            getData(dat, zone, dataFin, deepthmin, deepthmax, DataOutputDirectory, dateBeg, dateE, frequency)
    else:
        DataOutputDirectory = outputDirectory + dat + '/'
        getData(dat, zone, dataFin, deepthmin, deepthmax, DataOutputDirectory)
