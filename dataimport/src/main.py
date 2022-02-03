from general import readcsv, giveDateslist, getData

wantedData=['Temperature','Nitrate','Ammonium','Phosphate','Salinity','northward_Water_current','eastward_Water_current']
dateBeginning = '"2020-11-15 00:00:00"'
dateEnd = '"2020-11-22 00:00:00"'
zone='IBI'
deepthmin=0
deepthmax=20

# outputDirectory = 'I:/work-he/apps/safi/data/IBI/'
outputDirectory = '/media/share/data/IBI/'

dataFin=readcsv()
datesList=giveDateslist(dateBeginning,dateEnd)

for (dateBeg, dateE) in zip(datesList[0], datesList[1]):
    for dat in wantedData:
        DataOutputDirectory = outputDirectory + dat + '/'
        getData(dat, zone, dataFin, deepthmin, deepthmax, dateBeg, dateE, DataOutputDirectory)
