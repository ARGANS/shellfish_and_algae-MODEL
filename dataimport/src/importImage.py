import os

import numpy as np
import pandas as pd

outputFile='pic.nc'
outputDirectory='./'

file='./dataCmd.csv'
data = pd.read_csv(file, delimiter=',', header = 0)
numpData=data.to_numpy()
for k in range(len(numpData)):
    numpData[k]=[numpData[k][0].split(';')]
#print(numpData)

dataFin=[]
for k in range(len(numpData)):
    dataFin+=[numpData[k][0]]
dataFin=np.array(dataFin)

imgNb=6
datebeginning='"2021-10-22 00:00:00"'
dateEnd='"2021-11-22 00:00:00"'
#print(dataFin[imgNb])
for imgNb in [1,6,14,16,20,22,23,25,26,27,28,29,32]:
    print('--------------------------------------------------------')
    print(imgNb)
    print('--------------------------------------------------------')
    if dataFin[imgNb][18]=='NULL':
        os.system('python -m motuclient --motu '+dataFin[imgNb][3]+
              ' --service-id '+dataFin[imgNb][4]+' --product-id '+dataFin[imgNb][5]+
              ' --longitude-min '+dataFin[imgNb][6]+' --longitude-max '+dataFin[imgNb][7]+' --latitude-min '+dataFin[imgNb][8]+' --latitude-max '+dataFin[imgNb][9]+' --date-min '+datebeginning+
              ' --date-max '+dateEnd+' --variable '+dataFin[imgNb][15]+
              ' --out-dir '+outputDirectory+' --out-name '+outputFile+' --user mjaouen --pwd Azerty123456 ')
    else:
        os.system('python -m motuclient --motu ' + dataFin[imgNb][3] +
                  ' --service-id ' + dataFin[imgNb][4] + ' --product-id ' + dataFin[imgNb][5] +
                  ' --longitude-min ' + dataFin[imgNb][6] + ' --longitude-max ' + dataFin[imgNb][7] + ' --latitude-min ' +
                  dataFin[imgNb][8] + ' --latitude-max ' + dataFin[imgNb][9] +' --date-min '+datebeginning+
              ' --date-max '+dateEnd+' --depth-min ' + dataFin[imgNb][18] + '  --depth-max ' + dataFin[imgNb][19] + '  --variable ' +
                  dataFin[imgNb][15] +
                  ' --out-dir ' + outputDirectory + ' --out-name ' + outputFile + ' --user mjaouen --pwd Azerty123456 ')
