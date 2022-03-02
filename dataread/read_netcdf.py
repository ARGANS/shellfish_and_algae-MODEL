import datetime
import numpy as np
import netCDF4 as nc
import pandas as pd
import math


def iNearest(value, sortedArray, mode="closest"):
    # finds the index i in sortedArray such that sortedArray[i] is closest
    # to value. sortedArray should be 1D with no duplicates.
    # tip: if mode=="right" and value is in sortedArray then the function
    # outputs the position of value in sortedArray

    n = len(sortedArray)

    iInsert = np.searchsorted(sortedArray, value)

    # Out of bound cases
    if iInsert == 0:
        return 0
    if iInsert == n:
        return n - 1

    if mode=="closest":
        if (sortedArray[iInsert] - value) < (value - sortedArray[iInsert-1]):
            #value is nearer iInsert
            return iInsert
        else:
            return iInsert - 1
    elif mode=="left":
        return iInsert - 1
    elif mode=="right":
        return iInsert

def extractVarSubset(ncDataset, variable, **kwargs):
    # Extract a subset of ncDataset[variable] along any number of dimensions
    # to be specified in **kwargs. Returns the subset and a tuple containing
    # the names of the dimensions remaining in the subset.
    #
    # The argument/value pairs can be:
    #   - argument: either the name of a dimension of variable or one such
    #               name followed by the '_index' suffix.
    #   - value: either a single value or a pair of values in a tuple
    #
    # For example: 
    #   - 'longitude=54.3' will subset ncDataset[variable] along the value of
    #     longitude that is nearest 54.3.
    #   - 'longitude=(50.2, 56.4)' will get the data along all longitudes
    #     between 50.2 and 56.4, including the former and excluding the latter.
    #
    # The indices can also be specified directly: 
    #   - 'longitude_index=5' will get the data along the value of longitude at
    #     index 5.
    #   - 'longitude_index=(3, 14)' will get the data where the index of the
    #     longitude dimension is between 3 and 14, including the former and
    #     excluding the latter.
    #   - 'longitude_index=(3, 14, 2)' will get the data where the index of the
    #     longitude dimension is between 3 and 14, including the former and
    #     excluding the latter, with a step of 2.
    #
    # Note: any dimension name that is passed that is not a dimension of 
    #       variable will be silently ignored.
    # Note: If both 'dimension=...' and 'dimension_index=...' are used then
    #       the 'dimension_index=...' call will be ignored.
    # Note: tuple slices can contain None at one end or the other to extract
    #       all data to/from a given value/index. When subsetting by indices,
    #       negative indices can be used to count from the end, a negative step
    #       will also be accepted.
    # WARNING: when slicing a dimension by its value(s), this function
    #          assumes that the values that the dimnesion takes are in ascending
    #          order.

    dimensions = ncDataset[variable].get_dims()
    #sliceList = [slice(None) for _ in range(len(dimensions))]
    sliceList = [slice(None)] * len(dimensions)
    for iDim, dim in enumerate(dimensions):
        if dim.name in kwargs:
            sliceVal = kwargs[dim.name]
            if type(kwargs[dim.name]) is tuple:
                # find the indices of the lower and upper bounds, that are within
                # the demanded values.
                iLower = iNearest((sliceVal[0] or -math.inf), ncDataset[dim.name][:], mode="right")
                iUpper = iNearest((sliceVal[1] or math.inf), ncDataset[dim.name][:], mode="left")
                # include the lower, exclude the latter (in case of equality)
                sliceList[iDim] = slice(iLower, iUpper+1, None)
            else:
                # get the index where the dimension is nearest
                iSlice = iNearest(kwargs[dim.name], ncDataset[dim.name])
                sliceList[iDim] = iSlice

        elif dim.name+'_index' in kwargs:
            sliceVal = kwargs[dim.name+'_index']
            if type(sliceVal) is tuple:
                byArg = sliceVal[2] if (len(sliceVal) >= 3) else None
                sliceList[iDim] = slice(sliceVal[0], sliceVal[1], byArg)
            else:
                sliceList[iDim] = sliceVal

    # Find the dims that will remain after slicing (all except the dims where
    # the slice is an int).
    dimRemains = [type(a) is slice for a in sliceList]
    remainingDims = tuple(dim.name for i,dim in enumerate(dimensions) if dimRemains[i])

    return ncDataset[variable][sliceList], remainingDims


def weightsForAverage(sortedDim, leftLimit, rightLimit):
    # Computes the weights to be used for computed the weighed average between
    # leftLimit and rightLimit of a variable that takes its values at the
    # points in sortedDim. sortedDim must be in ascending order, leftLim and
    # rightLim must be values outside or at the extremes of sortedDim.

    n = len(sortedDim)

    fullDim = np.insert(sortedDim, 0, leftLimit)
    fullDim = np.append(fullDim, rightLimit)

    weights = ( fullDim[2:] - fullDim[:-2] ) / 2
    weights[0] = (sortedDim[0] + sortedDim[1]) / 2 - leftLimit
    weights[-1] = rightLimit - (sortedDim[-2] + sortedDim[-1]) / 2

    # normalize the result
    return weights / np.sum(weights)


def extractWithAverage(ncDataset, variable, averagingDims, weighted=True, **kwargs):
    # Extract a subset of ncDataset[variable] along any number of dimensions
    # to be specified in **kwargs, then average the dataset along the dimensions
    # specified in averagingDims (tuple of strings).
    # See extractVarSubset() for reference on the format of **kwargs.
    #
    # If weighed==True then the average is weighed. For dimension 'dim',
    # variable['dim'] is assumed to be a step function with bins centered around
    # the values that 'dim' takes. Only the values that 'dim' takes that are
    # within the bounds specified in **kwargs are considered.

    dataArray, dataDims = extractVarSubset(ncDataset, variable, **kwargs)

    weights = np.ones(dataArray.shape)
    if weighted:
        for nameAxis in averagingDims:

            dimVector, _ = extractVarSubset(ncDataset, nameAxis, **kwargs)

            # If axis bounds were specified by value, use these for weights.
            # This will fail if only one value of nameAxis is provided.
            leftLimit = kwargs[nameAxis][0] if (nameAxis in kwargs) else min(dimVector)
            rightLimit = kwargs[nameAxis][1] if (nameAxis in kwargs) else min(dimVector)

            dimWeights = weightsForAverage(dimVector, leftLimit, rightLimit)

            # Expand on all dataDims, with dimWeights aligned to the correct dim
            expandAxes = [i for i in range(len(dataDims)) if dataDims[i]!=nameAxis]
            dimWeights = np.expand_dims(dimWeights, expandAxes)

            weights = dimWeights * weights

    # indices of the axes in dataArray over which to average.
    averagingAxes = tuple(i for i,dim in enumerate(dataDims) if dim in averagingDims)

    averagedArray = np.average(dataArray, averagingAxes, weights=weights)

    # Dimensions that have not been averaged
    remainingDims = tuple(dim for dim in dataDims if dim not in averagingDims)

    return averagedArray, remainingDims

def averageOverMLD(array2D, mld, depthAxis):
    # array2D must be along (time, depth)
    # mld must be along (time)
    # depthAxis must be along (depth) (time dependent depth not handled right now)

    # Build weights to average over time dependent MLD
    weights = np.zeros(array2D.shape)
    for iTime in range(array2D.shape[0]):
        inMLD = (depthAxis <= mld[iTime])
        weights[iTime, inMLD] = weightsForAverage(depthAxis[inMLD], 0, mld[iTime])

    return np.average(array2D, 1, weights=weights)

def dateToNum(date, mode='Copernicus'):
    # Returns the numeric value associated to date in netCDF files
    if mode == "Copernicus":
        diff = date - datetime.datetime(1950, 1, 1)
        num = diff.days * 24 + diff.seconds/3600

    return num

def numToDate(num, mode='Copernicus'):
    # Returns the date associated to num that is found in netCDF files
    if mode == "Copernicus":
        date = datetime.datetime(1950, 1, 1) + datetime.timedelta(seconds = num*3600)

    return date



class ParamData:
    # Class to wrap all functions related the access to the data in one netCDF
    # file.
    # The objective is that all files can be accessed in the same manner from
    # higher level functions.

    def __init__(self, file_name, variable_name, latitude_name='latitude',
                 longitude_name='longitude', time_name='time', depth_name='depth',
                 date2num=dateToNum, num2date=numToDate):

        self.ds = nc.Dataset(file_name)
        self._variableName = variable_name
        self._dimNames = {
            'latitude': latitude_name,
            'longitude': longitude_name,
            'time': time_name,
            'depth': depth_name
        }

        self.date2num = date2num
        self.num2date = num2date


    def getVariable(self, variable=None, **kwargs):
        # **kwargs may contain latitude, longitude, time, and depth arguments.
        # The value of the arguments can be anything that is accepted by
        # extractVarSubset()
        # The time argument should be specified with datetime.datetime()
        # object(s).
        # variable can be None to extract the interest variable from self.ds or
        # any variable name. If variable is one of latitude/longitude/time/depth
        # then the name is transparently altered to fit ds.

        if variable is None:
            variableName = self._variableName
        elif variable in self._dimNames:
            variableName = self._dimNames[variable]

        # Change the time input(s) for datetime() to numeric values
        if 'time' in kwargs:
            if type(kwargs['time']) is tuple:
                kwargs['time'] = tuple(self.date2num(date) for date in kwargs['time'])
            else:
                kwargs['time'] = self.date2num(kwargs['time'])

        newKwargs = {self._dimNames[dim]: kwargs[dim] for dim in kwargs.keys()}

        output, _ = extractVarSubset(self.ds, variableName, **newKwargs)

        # Change the output from numeric values to datetime() if we output time
        if variable == 'time':
            output = np.ma.masked_array([self.num2date(a) for a in output])

        # TODO: ensure that the remaining dimensions are always output in the
        # same order, notably (time, depth)

        return output


    def __del__(self):
        self.ds.close()



class AllData:

    def __init__(self, fileNameList, parameterNameList, variableNameList, latitudeNameList, 
                 longitudeNameList, timeNameList, depthNameList):

        self.parameterData = {}
        for fileName, parameterName, variableName, latitudeName, longitudeName, timeName, \
            depthName in zip(fileNameList, parameterNameList, variableNameList, latitudeNameList, \
            longitudeNameList, timeNameList, depthNameList):

            self.parameterData[parameterName] = ParamData(file_name=fileName,
                    variable_name=variableName, latitude_name=latitudeName,
                    longitude_name=longitudeName, time_name=timeName,
                    depth_name=depthName)


    def getTimeSeries(self, latitude, longitude, dateRange, depth, parameters=None):
        # Gets the time series of parameters at the given coordinates and
        # within the dateRange.

        if parameters is None:
            parameters = self.parameterData.keys()

        df = pd.DataFrame()

        nDays = (dateRange[1] - dateRange[0]).days
        df['date'] = [dateRange[0] + datetime.timedelta(days=days) for days in range(nDays)]

        for param in parameters:
            data = self.parameterData[param].getVariable(latitude=latitude, longitude=longitude,
                                                time=dateRange, depth=depth)
            timeAxis = self.parameterData[param].getVariable(variable='time', latitude=latitude, 
                                                longitude=longitude, time=dateRange, depth=depth)

            addData = pd.DataFrame({'date': timeAxis, param:data})
            # change masked values to None
            addData.loc[np.ma.getmaskarray(data), [param]] = None

            df = pd.merge_ordered(df, addData, how='left')

        return(df)


    def getTimeSeriesInMLD(self, latitude, longitude, dateRange, parameters=None,
                           mldName='ocean_mixed_layer_thickness'):
        # Gets the time series of parameters at the given coordinates and
        # within the dateRange. Then averages the values that are within the
        # MLD, weighing in consideration of the depth axis.

        if parameters is None:
            parameters = list(self.parameterData.keys())

        df = pd.DataFrame()

        nDays = (dateRange[1] - dateRange[0]).days
        df['date'] = [(dateRange[0] + datetime.timedelta(days=days)) for days in range(nDays)]

        # Ensure that the MLD gets imported first
        paramNames = parameters.copy()
        if mldName in paramNames:
            paramNames.insert(0, paramNames.pop(paramNames.index(mldName)))
        else:
            paramNames.insert(0, mldName)

        for param in paramNames:
            print(param)
            data = self.parameterData[param].getVariable(latitude=latitude, longitude=longitude,
                                                time=dateRange)
            timeAxis = self.parameterData[param].getVariable(variable='time', latitude=latitude,
                                                longitude=longitude, time=dateRange)

            addData = pd.DataFrame({'date': timeAxis})
            print(addData.duplicated(['date']).sum())

            if data.ndim == 2: # data has depth
                # Get the MLD at the times in paramTime
                paramMLD = pd.merge_ordered(df, addData, on='date', how='right')[mldName]
                notnaMLD = np.array(paramMLD.notna())

                # Get the depth axis of param
                depthAxis = self.parameterData[param].getVariable(variable='depth',
                                            latitude=latitude, longitude=longitude, time=dateRange)

                # Compute the weighted average
                addParam = np.ma.masked_array([None] * len(timeAxis))
                addParam[notnaMLD] = averageOverMLD(data[notnaMLD,:], np.array(paramMLD[notnaMLD]), depthAxis)

            else: # data has no depth
                addParam = data
 
            addData[param] = addParam
            # change masked values to None
            addData.loc[np.ma.getmaskarray(addParam), [param]] = None

            df = pd.merge_ordered(df, addData, on='date', how='left')
            print(df.duplicated(['date']).sum())

        return df


    def __del__(self):
        for _, parData in enumerate(self.parameterData):
            del parData



if __name__ == "__main__":
    lat = 51.587433
    lon = -9.897116
    zone = 'IBI'

    startDate = datetime.datetime(2020, 1, 1, 12)
    endDate = datetime.datetime(2021, 1, 31, 12)

    mainpath = 'I:/work-he/apps/safi/data/IBI/'

    dataRef = pd.read_csv('I:/work-he/apps/safi/data/IBI/dataCmd.csv', delimiter=';')

    # The 'Names' informations should probably be in dataRef, which would make it easier to
    # define and read.
    paramNames = ['Ammonium', 'Nitrate', 'Temperature', 'northward_Water_current', 'eastward_Water_current', 'ocean_mixed_layer_thickness', 'par']
    fileNames = [mainpath + f"merged_files/merged_{param}_{zone}.nc" for param in paramNames]
    variableNames = [dataRef.loc[(dataRef['Parameter']==param) & (dataRef['Place']==zone)].reset_index()['variable'][0] for param in paramNames]

    nParams = len(paramNames)
    latitudeNames = ['latitude'] * nParams
    longitudeNames = ['longitude'] * nParams
    timeNames = ['time'] * nParams
    depthNames = ['depth'] * nParams
    # For par
    latitudeNames[-1] = 'lat'
    longitudeNames[-1] = 'lon'

    algaeData = AllData(fileNameList=fileNames,
                        parameterNameList=paramNames,
                        variableNameList=variableNames,
                        latitudeNameList=latitudeNames,
                        longitudeNameList=longitudeNames,
                        timeNameList=timeNames,
                        depthNameList=depthNames
    )

    #df = algaeData.getTimeSeriesInMLD(lat, lon, (startDate, endDate), parameters=['Ammonium', 'par'])
    #df = algaeData.getTimeSeriesInMLD(lat, lon, (startDate, endDate), parameters=['Ammonium', 'Nitrate', 'Temperature', 'ocean_mixed_layer_thickness', 'par'])
    df = algaeData.getTimeSeriesInMLD(lat, lon, (startDate, endDate))
    print(df)
    df.to_csv(mainpath+'Bantry_data/bantry_MLDaveraged.csv',
              index=False, sep=';')

    del algaeData