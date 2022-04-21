import os
from general import readcsv, giveDateslist, getData
import pandas as pd
import pprint

year = os.getenv('AC_YEAR')
zone = os.getenv('AC_ZONE')
outputDirectory = os.getenv('AC_OUTPUT_DIR')
deepthmin = int(os.getenv('AC_DEPTHMIN'))
deepthmax = int(os.getenv('AC_DEPTHMAX'))

#  TODO
wantedData=['par']

dateBeginning = f'{year}-01-01 00:00:00'
dateEnd = f'{year}-01-01 00:00:00'


print(f'{year} {zone} {deepthmin} {deepthmax}')
dataFin = pd.read_csv('./dataCmd.csv',';')

print('dataFin')
pprint.pprint(dataFin)

datesList = giveDateslist(dateBeginning, dateEnd)
print('datesList')
pprint.pprint(datesList)
