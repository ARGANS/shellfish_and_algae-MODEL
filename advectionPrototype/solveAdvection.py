import copy
import os
import time
import netCDF4 as nc
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from advectionPrototype.saveAsTiff import saveAsTiff, getMetadata
from dataread.read_netcdf import extractVarSubset

# extract the data value at depth in the merged files (all the daily data merged in one file)
def getwantedMergeData(data, depth, dataCmdpath, mergedFilepath='I:/work-he/apps/safi/data/IBI/merged_files/'):
    csvFile = pd.read_csv(dataCmdpath, ';')
    DataLine = csvFile.loc[csvFile["Parameter"] == data]
    variable = DataLine.iloc[0]["variable"]
    # we search after the file containing the wanted data
    for r, d, f in os.walk(mergedFilepath):
        for i in range(len(f)):
            filename = f[i]
            # if it is the wanted file
            if filename[:10] == 'merged_' + data[:3]:
                fn = mergedFilepath + f[i]
                ncDataset = nc.Dataset(fn)
                break
    return extractVarSubset(ncDataset, variable, depth=depth)[0], ncDataset


if __name__ == "__main__":
    # we get the data at 3 meters
    nwc = 'northward_Water_current'
    ewc = 'eastward_Water_current'
    nitrate = 'Nitrate'
    depth = 3

    dataCmdpath = 'I:/work-he/apps/safi/data/IBI/dataCmd.csv'
    dataNwc, _ = getwantedMergeData(nwc, depth, dataCmdpath)
    dataEwc, _ = getwantedMergeData(ewc, depth, dataCmdpath)
    dataNO3, ds = getwantedMergeData(nitrate, depth, dataCmdpath)

    mask = dataNwc.mask
    maskposition = np.where(mask == True)  # we get the position of the masked data
    maskpos2D = np.where(mask[0] == True)
    dataNwc[0][maskpos2D] = np.nan
    diffCelli_1j_1 = np.where(
        ~np.isnan(dataNwc[0][:-1, 1:]) & ~np.isnan(dataNwc[0][1:, 1:]) & ~np.isnan(dataNwc[0][1:, :-1]))#we get the position of the usable cells
    northNotDiffCelli_1j_1 = np.where(np.isnan(dataNwc[0][:-1, 1:]) | np.isnan(dataNwc[0][1:, 1:]))#we get the position of the unusable cells
    eastNotDiffCelli_1j_1 = np.where(np.isnan(dataNwc[0][1:, :-1]) | np.isnan(dataNwc[0][1:, 1:]))#we get the position of the unusable cells

    diffCell = (diffCelli_1j_1[0] + 1, diffCelli_1j_1[1] + 1)
    diffCelli_1j = (diffCelli_1j_1[0], diffCelli_1j_1[1] + 1)
    diffCellij_1 = (diffCelli_1j_1[0] + 1, diffCelli_1j_1[1])

    northNotDiffCellij = (northNotDiffCelli_1j_1[0] + 1, northNotDiffCelli_1j_1[1] + 1)
    northNotDiffCelli_1j = (northNotDiffCelli_1j_1[0], northNotDiffCelli_1j_1[1] + 1)

    eastNotDiffCellij = (eastNotDiffCelli_1j_1[0] + 1, eastNotDiffCelli_1j_1[1] + 1)
    eastNotDiffCellij_1 = (eastNotDiffCelli_1j_1[0] + 1, eastNotDiffCelli_1j_1[1])

    init = time.time()
    DC = np.zeros(np.shape(dataNwc[0]))
    #we define the time and space discretisation
    discr = 160
    dt = 1 / discr
    dx = 3
    dy = 3
    nbrStep = len(dataNwc) * discr
    randTime = 0
    dataNwc = dataNwc * 60 * 60 * 24 /(1000*discr)
    dataEwc = dataEwc * 60 * 60 * 24 /(1000*discr)
    #for each time step
    for k in range(nbrStep):
        dayNbr = k // discr #we find which day we are working
        DC[diffCell] = DC[diffCell] * (
                    1 - dt * (dataNwc[dayNbr][diffCell] / dx) - dt * (dataEwc[dayNbr][diffCell] / dy)) + DC[
                           diffCelli_1j] * dataNwc[dayNbr][diffCell] * dt / dx \
                       + DC[diffCellij_1] * dataEwc[dayNbr][diffCell] * dt / dx \
                       - np.random.rand(len(DC[diffCell])) * 0.5 * dt * (dataNO3[dayNbr][diffCell] + DC[diffCell]) #we update DC value
        DC[northNotDiffCellij] = DC[northNotDiffCelli_1j] #we impose neumann condition as boundary condition
        DC[eastNotDiffCellij] = DC[eastNotDiffCellij_1]
        DC[0] = DC[1]
        DC[:, 0] = DC[:, 1]
    print(time.time() - init - randTime)
    print(randTime)

    xsize, ysize, ulx, uly, xres, yres = getMetadata(ds)

    DC[maskpos2D] = np.nan

    Cfin = dataNO3[dayNbr] + DC
    Cfin[maskpos2D] = np.nan

    saveAsTiff(Cfin, xsize, ysize, ulx, uly, xres, yres,
               "I:/work-he/apps/safi/data/IBI/lastDayValueDeficitAdvection4.tiff")
    saveAsTiff(DC, xsize, ysize, ulx, uly, xres, yres,
               "I:/work-he/apps/safi/data/IBI/lastDayDifferenceDeficitAdvection4.tiff")
    fig, ax = plt.subplots()
    plt.imshow(Cfin)
    plt.colorbar()
    ax.invert_yaxis()
    plt.clim(0, 10)
    ax.set_title('final')
    plt.show()

    fig, ax = plt.subplots()
    plt.imshow(DC)
    plt.colorbar()
    ax.invert_yaxis()
    plt.clim(-10, 0)
    ax.set_title('final')
    plt.show()
