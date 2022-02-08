import datetime
import os
import numpy as np
import netCDF4 as nc
import pandas as pd
import time
import math


def iNearest(value, sortedArray):
    # finds the index i in sortedArray such that sortedArray[i] is closest
    # to value. sortedArray should be 1D.

    n = len(sortedArray)

    iInsert = np.searchsorted(sortedArray, value)

    # Out of bound cases
    if iInsert == 0:
        return 0
    if iInsert == n:
        return n-1

    if (sortedArray[iInsert] - value) < (value - sortedArray[iInsert-1]):
        #value is nearer iInsert
        return iInsert
    else:
        return iInsert - 1


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
    # WARNING: when slicing a dimension by its exact value, this function
    #          assumes that the values that the dimesion takes are in ascending
    #          order.

    dimensions = ncDataset[variable].get_dims()
    sliceList = [slice(None) for _ in range(len(dimensions))]
    for iDim, dim in enumerate(dimensions):
        if dim.name in kwargs:
            sliceVal = kwargs[dim.name]
            if type(kwargs[dim.name]) is tuple:
                lowerBound = ncDataset[dim.name][:] >= (sliceVal[0] or -math.inf)
                upperBound = ncDataset[dim.name][:] < (sliceVal[1] or math.inf)
                sliceList[iDim] = np.squeeze(np.where(np.logical_and(lowerBound, upperBound)))
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


def extractWithAverage(ncDataset, variable, averagingDims, weighed=True, **kwargs):
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
    if weighed:
        for nameAxis in averagingDims:

            dimVector, _, _ = extractVarSubset(ncDataset, nameAxis, **kwargs)

            dimWeights = weightsForAverage(dimVector, min(dimVector), max(dimVector))

            # Expand on all dataDims, with dimWeights aligned to the correct dim
            expandAxes = [i for i in range(len(dataDims)) if dataDims[i]!=nameAxis]
            dimWeights = np.expand_dims(dimWeights, expandAxes)

            weights = dimWeights * weights

    # indices of the axes in dataArray over which to average.
    averagingAxes = tuple(i for i,dim in enumerate(dataDims) if dim in averagingDims)

    averagedArray = np.average(dataArray, averagingAxes, weights=weights)

    return averagedArray


def dateToNum(date, mode):
    # Returns the numeric value associated to date in netCDF files
    if mode == "Copernicus":
        diff = date - datetime.datetime(1950, 1, 1)
        num = diff.days * 24 + diff.seconds/3600

    return num

def numToDate(num, mode):
    # Returns the date associated to num that is found in netCDF files
    if mode == "Copernicus":
        date = datetime.datetime(1950, 1, 1) + datetime.timedelta(seconds = num*3600)

    return date


if __name__ == "__main__":
    lat = 51.587433
    lon = -9.897116
    zone = 'IBI'

    startDate = datetime.datetime(2020, 1, 25, 12)
    endDate = datetime.datetime(2021, 1, 1, 12)

    mainpath = 'I:/work-he/apps/safi/data/IBI/'

    dataRef = pd.read_csv('I:/work-he/apps/safi/data/IBI/dataCmd.csv', delimiter=';')

    fn = mainpath + "merged_files/merged_Ammonium_IBI.nc"

    dsParam = nc.Dataset(fn)
    dsMLD = nc.Dataset(mainpath + "merged_files/merged_ocean_mixed_layer_thickness_IBI.nc")

    startTime = dateToNum(startDate, "Copernicus")
    endTime = dateToNum(endDate, "Copernicus")

    extractArgs = {'latitude': lat,
                   'longitude': lon,
                   'time': (startTime, endTime)
                   }

    mld, mldDims = extractVarSubset(dsMLD, 'mlotst', **extractArgs)
    mldTime, _ = extractVarSubset(dsMLD, 'time', **extractArgs)

    param, paramDims = extractVarSubset(dsParam, 'nh4', **extractArgs)
    paramTime, _ = extractVarSubset(dsParam, 'time', **extractArgs)
    paramDepth, _ = extractVarSubset(dsParam, 'depth', **extractArgs)

    # TODO: ensure paramDims is (time, depth)
    # TODO: ensure time axes match

    # Build weights to average over time dependent MLD
    weights = np.zeros(param.shape)
    for iTime in range(param.shape[0]):
        inMLD = (paramDepth <= mld[iTime])
        weights[iTime, inMLD] = weightsForAverage(paramDepth[inMLD], 0, mld[iTime])


    dsParam.close()
    dsMLD.close()

    dataDict = {}
    dataDict['date'] = [numToDate(a, "Copernicus") for a in paramTime]

    #paramNames = ['Ammonium','Nitrate','Temperature','northward_Water_current','eastward_Water_current', 'ocean_mixed_layer_thickness']
    datNames = ['Ammonium','Nitrate']
    for dat in datNames:

        fileName = mainpath + f"merged_files/merged_{dat}_{zone}.nc"
        paramName = dataRef.loc[(dataRef['Parameter']==dat) & (dataRef['Place']==zone)].reset_index()['variable'][0]

        dsParam = nc.Dataset(fileName)

        param, paramDims = extractVarSubset(dsParam, paramName, **extractArgs)

        paramAveraged = np.average(param, 1, weights=weights)

        dataDict[paramName] = paramAveraged

        dsParam.close()

    df = pd.DataFrame(dataDict)

    print(df)
    df.to_csv(mainpath+'Bantry_data/bantry_MLDaveraged_2020.csv',
                index=False, sep=';')