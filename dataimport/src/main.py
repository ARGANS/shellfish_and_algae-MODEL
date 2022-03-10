from general import readcsv, giveDateslist, getData
import pandas as pd

wantedData=['Temperature','Nitrate','Ammonium','Phosphate','Salinity','northward_Water_current','eastward_Water_current']
dateBeginning = '2020-11-15 00:00:00'
dateEnd = '2020-11-22 00:00:00'
zone='IBI'
deepthmin=0
deepthmax=20

outputDirectory = 'I:/work-he/apps/safi/data/IBI/'
#outputDirectory = '/media/share/data/IBI/'

dataFin=pd.read_csv('./dataCmd.csv',';')

datesList=giveDateslist(dateBeginning,dateEnd)

for dat in wantedData:
    DataLine = dataFin.loc[dataFin["Parameter"] == dat]
    if DataLine.iloc[0]["daily"] == 1:
        for (dateBeg, dateE) in zip(datesList[0], datesList[1]):
            DataOutputDirectory = outputDirectory + dat + '/'
            getData(dat, zone, dataFin, deepthmin, deepthmax,  DataOutputDirectory, dateBeg, dateE)
    else:
        DataOutputDirectory = outputDirectory + dat + '/'
        getData(dat, zone, dataFin, deepthmin, deepthmax, DataOutputDirectory)
