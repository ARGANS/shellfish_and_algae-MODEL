import argparse
import os
from osgeo import gdal, osr
import numpy as np

from dataread.src.make_runs import dataCmd_to_AllData
from dataread.src.utils import import_json
import rasterio
import os
import multiprocessing as mp
from dataread.src.make_runs import *
from dataread.src.launch_model import MA_model_scipy
from dataread.src.models.ModelProperties import ModelProperties
from dataread.src.read_netcdf import *
import netCDF4 as nc


# save dataArray in a geotiff file
def saveAsTiff(dataArray, xsize, ysize, ulx, uly, xres, yres, filepath, units):
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(filepath, xsize, ysize, 1, gdal.GDT_Float32)

    # this assumes the projection is Geographic lat/lon WGS 84
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    ds.SetProjection(srs.ExportToWkt())

    gt = [ulx, xres, 0, uly, 0, -yres]
    ds.SetGeoTransform(gt)

    outband = ds.GetRasterBand(1)
    outband.WriteArray(np.flip(dataArray, 0))
    # ds.units = 'meters'
    del ds
    meta_data_dict = {"units": units}
    with rasterio.open(filepath, 'r+') as dataset:
        dataset.update_tags(**meta_data_dict)

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

    args = vars(ap.parse_args())

    # Parse arguments
    ncFile = args["input"]

    model_properties = ModelProperties(os.getenv('INPUT_SOURCE'), os.getenv('INPUT_DESTINATION'))
    try:
        model_properties.parse(os.getenv('INPUT_MODEL_PROPERTIES_JSON'))
    except:
        raise RuntimeError('Cannot parse the value of the parameters_json environment variable')

    full_json = model_properties.attrs

    dict_dataCmd = full_json['dataset_parameters']['datasets']
    dict_to_AllData = dataCmd_to_AllData(dict_dataCmd, model_properties.file_template)

    json_data = AllData(dict_to_AllData)

    zone = "MED"

    unitsDict = {'NO3': 'mg/m^3',
                 'NH4': 'mg/m^3',
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
    working_path = f"I:/work-he/apps/safi/data/{zone}"

    latitudes = nc.Dataset(ncFile)['latitude'][:]
    longitudes = nc.Dataset(ncFile)['longitude'][:]
    xsize, ysize, ulx, uly, xres, yres = len(longitudes), len(latitudes), longitudes[0], latitudes[-1], longitudes[1] - \
                                         longitudes[0], latitudes[1] - latitudes[0]

    dataDict = {}
    for dataName in ['NH4', 'NO3', 'N_s', 'N_f', 'D', 'cNH4', 'cNO3']:
        dataDict[dataName] = nc.Dataset(ncFile)[dataName][0]
        dataDict[dataName][np.where(dataDict[dataName] > 1e12)] = np.nan

    paramSacch = json_data['parameters']['species']['alaria']['parameters']

    density_MA = json_data['parameters']['farm']['default']['parameters']['density_MA']

    dataDict['DW'] = dataDict['N_f'] / paramSacch['Q_min']  # gDW m-3
    dataDict['DW_line'] = dataDict['DW'] * paramSacch['h_MA'] * paramSacch['w_MA'] / 1000  # kg/m (i.e. per m of rope)
    dataDict['DW_PUA'] = dataDict['DW'] * paramSacch[
        'h_MA'] * density_MA / 1000  # kg/m^2 (multiply be density to account for unused space within farm)

    dataDict['FW'] = dataDict['DW'] / paramSacch['DF_MA']  # gFW m-3
    dataDict['FW_line'] = dataDict['DW_line'] / paramSacch['DF_MA']  # kg/m (i.e. per m of rope)
    dataDict['FW_PUA'] = dataDict['DW_PUA'] / paramSacch['DF_MA']  # kg/m^2

    # Energy
    dataDict['kcal_PUA'] = dataDict['DW'] * paramSacch['h_MA'] * density_MA * paramSacch['kcal_MA']  # kcal/m^2

    # protein
    dataDict['protein_PUA'] = dataDict['DW_PUA'] * paramSacch['prot_MA']  # kg/m^2

    # CO2 uptake
    dataDict['Biomass_CO2'] = (dataDict['N_f'] / 14) * paramSacch['CN_MA'] * 44 / 1000  # g (CO2) /m^3    (44 is rmm of CO2)
    dataDict['CO2_uptake_PUA'] = dataDict['Biomass_CO2'] * paramSacch['h_MA'] * density_MA / 1000  # kg (CO2) / m^2

    for dataName in ['NO3','NH4','D','N_f','N_s','DW','DW_PUA','FW','FW_line','FW_PUA','kcal_PUA','protein_PUA','Biomass_CO2','CO2_uptake_PUA','DW_line']:
        saveAsTiff(dataDict[dataName], xsize, ysize, ulx, uly, xres, yres, f"{working_path}/{dataName}.tif", unitsDict[dataName])

if __name__ == "__main__":
    main()

#python D:/Profils/mjaouen/Documents/alternance/EASME/gitcode/shellfish_and_algae-MODEL/pretreatment/src/resampleGrid.py -i D:/Profils/mjaouen/Documents/alternance/EASME/data/NWS/merged_Temperature_NWS.nc -d I:/work-he/apps/safi/data/NWS -o I:/work-he/apps/safi/data/NWS/output.nc -v thetao -x longitude -y latitude