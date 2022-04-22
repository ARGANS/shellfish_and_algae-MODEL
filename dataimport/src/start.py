import os
from general import readcsv, giveDateslist, getData
import pandas as pd
from pprint import pprint

year = int(os.getenv('AC_YEAR'))
zone = os.getenv('AC_ZONE')
outputDirectory = os.getenv('AC_OUTPUT_DIR')
deepthmin = int(os.getenv('AC_DEPTHMIN'))
deepthmax = int(os.getenv('AC_DEPTHMAX'))

wantedData = ['Temperature', 'Nitrate', 'Ammonium', 'eastward_Water_current', 'northward_Water_current']

dateBeginning = f'{year}-11-01 00:00:00'
dateEnd = f'{year}-12-31 23:59:59'


print(f'{year} {zone} {deepthmin} {deepthmax}')
dataFin = pd.read_csv('./dataCmd.csv',';')
datesList = giveDateslist(dateBeginning, dateEnd)

print('dataFin')
pprint(dataFin)
print('datesList')
pprint(datesList)
for dat in wantedData:
    dataOutputDirectory = outputDirectory + dat + '/'
    dataLine = dataFin.loc[dataFin["Parameter"] == dat]
    print(f'dataLine {dataOutputDirectory}')
    pprint(dataLine)
    if dataLine.iloc[0]["daily"] == 1:
        for (dateBeg, dateE) in zip(datesList[0], datesList[1]):
            getData(dat, zone, dataFin, deepthmin, deepthmax,  dataOutputDirectory, dateBeg, dateE)
    else:
        getData(dat, zone, dataFin, deepthmin, deepthmax, dataOutputDirectory)
