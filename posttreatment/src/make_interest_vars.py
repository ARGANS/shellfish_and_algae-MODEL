import argparse
import json
import numpy as np
import netCDF4 as nc

def import_json(path: str) -> dict:
    with open(path, 'r') as f:
        data = json.loads(f.read())
    return data 


def getParameters(json_data):
    specie = list(json_data['parameters']['species'].keys())[0]
    paramAlgae = json_data['parameters']['species'][specie]['parameters']
    density_MA = json_data['parameters']['farm']['default']['parameters']['density_MA']
    return paramAlgae, density_MA


def main():
    '''
    This script aligns the variable "variable" in the file at "dir/input" to a
    regular lat/lon grid based on the "latName"/"lonName" axes of the file at
    "inputGrid".
    The result is stored in "dir/output".
    '''

    # Create argument parser
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, help="The path to the input file")
    # ap.add_argument("-o", "--outputDir", required=True, help="The path to the output directory")
    ap.add_argument("-j", "--jsonFile", required=True, help="The path to the parameters json")

    args = vars(ap.parse_args())

    # Parse arguments
    ncFile = args["input"]
    # outDir = args["outputDir"]
    jsonFile = args["jsonFile"]

    unitsDict = {'CMEMS_NO3': 'mg/m^3',
                 'CMEMS_NH4': 'mg/m^3',
                 'DW': 'gDW m-3',
                 'DW_line': 'kg/m',
                 'DW_PUA': 'kg/m^2',
                 'FW': 'gFW m-3',
                 'FW_line': 'kg/m',
                 'FW_PUA': 'kg/m^2',
                 'kcal_PUA': 'kcal/m^2',
                 'protein_PUA': 'kg/m^2',
                 'Biomass_CO2': 'g (CO2) /m^3',
                 'CO2_uptake_PUA': 'kg (CO2) / m^2',
                 'NO3field': 'mg N/m^3',
                 'NH4field': 'mg N/m^3',
                 'D': 'mg N/m^3',
                 'N_f': 'mg N/m^3',
                 'N_s': 'mg N/m^3',
                 'avNO3': 'mg N/m^3',
                 'avNH4': 'mg N/m^3'}

    json_data = import_json(jsonFile)
    paramAlgae, density_MA = getParameters(json_data)
    dataDict = {}
    ds_read = nc.Dataset(ncFile, 'r')
    for dataName in ['CMEMS_NH4', 'CMEMS_NO3', 'N_s', 'N_f', 'D', 'avNH4', 'avNO3']:
        dataDict[dataName] = ds_read[dataName][-1,:,:]
        #dataDict[dataName].set_fill_value(np.nan)
        dataDict[dataName][dataDict[dataName].mask]=np.nan

    dataDict['DW'] = dataDict['N_f'] / paramAlgae['Q_min']  # gDW m-3
    dataDict['DW_line'] = dataDict['DW'] * paramAlgae['h_MA'] * paramAlgae['w_MA'] / 1000  # kg/m (i.e. per m of rope)
    dataDict['DW_PUA'] = dataDict['DW'] * paramAlgae[
        'h_MA'] * density_MA / 1000  # kg/m^2 (multiply be density to account for unused space within farm)

    dataDict['FW'] = dataDict['DW'] / paramAlgae['DF_MA']  # gFW m-3
    dataDict['FW_line'] = dataDict['DW_line'] / paramAlgae['DF_MA']  # kg/m (i.e. per m of rope)
    dataDict['FW_PUA'] = dataDict['DW_PUA'] / paramAlgae['DF_MA']  # kg/m^2

    # Energy
    dataDict['kcal_PUA'] = dataDict['DW'] * paramAlgae['h_MA'] * density_MA * paramAlgae['kcal_MA']  # kcal/m^2

    # protein
    dataDict['protein_PUA'] = dataDict['DW_PUA'] * paramAlgae['prot_MA']  # kg/m^2

    # CO2 uptake
    dataDict['Biomass_CO2'] = (dataDict['N_f'] / 14) * paramAlgae[
        'CN_MA'] * 44 / 1000  # g (CO2) /m^3    (44 is rmm of CO2)
    dataDict['CO2_uptake_PUA'] = dataDict['Biomass_CO2'] * paramAlgae['h_MA'] * density_MA / 1000  # kg (CO2) / m^2

    newOutDataNames = ['DW', 'DW_line', 'DW_PUA', 'FW', 'FW_line', 'FW_PUA', 'kcal_PUA', 'protein_PUA', 'Biomass_CO2',
                    'CO2_uptake_PUA']
    ds_read.close()

    ds_append = nc.Dataset(ncFile, 'a')
    print(ds_append.variables)
    for name in newOutDataNames:
        var = ds_append.createVariable(name, 'f4', ('latitude', 'longitude',))
        var[:, :] = dataDict[name]
        var.units = unitsDict[name]

    ds_append.close()

    '''for dataName in outDataNames:
        dataDict[dataName] = nc.Dataset(ncFile)[dataName][0]
        dataDict[dataName][np.where(dataDict[dataName] > 1e12)] = np.nan
        saveAsTiff(dataDict[dataName], xsize, ysize, ulx, uly, xres, yres, f"{outDir}/{dataName}.tif", unitsDict[dataName])'''


if __name__ == "__main__":
    main()