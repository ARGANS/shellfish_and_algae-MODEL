import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('p:/Aquaculture/shellfish_and_algae-MODEL/dataread/src/')
from read_netcdf import *


def plot_map(ds, variable, colorRange=(None, None), **kwargs):
    """Plots a 2D map of variable in ds, using kwargs to subset the data with
    extractVarSubset(), return a figure object.
    """
    data, dims = extractVarSubset(ds, variable, **kwargs)

    # There should remain exactly 2 dimensions i dims
    # TODO: deal with the case where the dimension has no corresponding
    # variable
    x, _ = extractVarSubset(ds, dims[1], **kwargs)
    y, _ = extractVarSubset(ds, dims[0], **kwargs)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    im = ax.pcolormesh(x, y, data, vmin=colorRange[0], vmax=colorRange[1])

    ax.set_xlabel(dims[1])
    ax.set_ylabel(dims[0])

    cbar = plt.colorbar(im)
    # TODO: include units and long name in netcdf to use here
    cbar.set_label(variable)


    return fig


if __name__=="__main__":

        ds = nc.Dataset("D:/data_scenario_B/NWS/test_2.nc")

        print(ds.variables)


        #fig = plot_map(ds, 'N_f', time_index=-1, longitude=(-5, 5), latitude=(48, 53))
        fig = plot_map(ds, "N_f", time_index=0)#, colorRange=(-1000,0))
        #fig = plot_map(ds, "par", time_index=11, colorRange=(0,30))
        plt.show()

        """
        variables = {
            "NH4": {max: 100},
            "NO3": {max: 1000},
            "N_s": {max: 120000},
            "N_f": {max: 100000},
            "D": {max: 1},
        }
        for var,params in variables.items():
            for month in range(1,13):
                im = plot_map(ds, var, colorRange=(-1, params[max]), time=month)
                im.set_size_inches(10, 10)
                im.savefig(f'D:/IBI/maps/monthly_simulations_{var}_month_{month:02d}.png')
        """