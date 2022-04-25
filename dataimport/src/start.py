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

dateBeginning = f'{year}-01-01 00:00:00'
dateEnd = f'{year + 1}-01-01 00:00:00'
frequency = 2

print(f'{year} {zone} {deepthmin} {deepthmax}')
dataFin = pd.read_csv('./dataCmd.csv',';')
datesList = giveDateslist(dateBeginning, dateEnd, frequency)

print('dataFin')
pprint(dataFin)
print('datesList')
pprint(datesList)
for dat in wantedData:
    dataOutputDirectory = outputDirectory + dat + '/'
    dataLine = dataFin.loc[dataFin["Parameter"] == dat]
    print(f'dataLine {dataOutputDirectory}')
    pprint(dataLine)
    if dataLine.iloc[0]["daily"] > 0:
        getData(dat, zone, dataFin, deepthmin, deepthmax,  dataOutputDirectory, datesList[0], datesList[1],frequency)
    else:
        getData(dat, zone, dataFin, deepthmin, deepthmax, dataOutputDirectory)
