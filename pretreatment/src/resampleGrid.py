import argparse
import os

import netCDF4 as nc



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
    ap.add_argument("-g", "--inputGrid", required=True, help="The path to the input grid file")
    ap.add_argument("-d", "--dir", required=True, help="The path to the working directory")
    ap.add_argument("-o", "--output", required=True, help="The path to the output file")
    ap.add_argument("-v", "--variableInput", required=True, help="The variable name in the input file")
    ap.add_argument("-w", "--variableGrid", required=True, help="The variable name in the grid file")
    ap.add_argument("-x", "--lonName", required=True, help="The longitude name")
    ap.add_argument("-y", "--latName", required=True, help="The latitude name")
    args = vars(ap.parse_args())

    # Parse arguments
    inputFile = args["input"]
    inputGrid = args["inputGrid"]
    outputFile = args["output"]
    directory = args["dir"]
    variableInput = args["variableInput"]
    variableGrid = args["variableGrid"]
    lonName = args["lonName"]
    latName = args["latName"]

    inputData = nc.Dataset(inputGrid)

    fillval = str(inputData.variables[variableGrid]._FillValue)

    lon_try_names = [lonName, 'longitude', 'lon', 'long']
    for try_name in lon_try_names:
        if try_name in inputData.dimensions.keys():
            longitudes = inputData.variables[try_name]
            break

    lat_try_names = [latName, 'latitude', 'lat']
    for try_name in lat_try_names:
        if try_name in inputData.dimensions.keys():
            latitudes = inputData.variables[try_name]
            break

    xres = longitudes[1] - longitudes[0]
    yres = latitudes[1] - latitudes[0]

    xmin, xmax = longitudes[0], longitudes[-1]
    ymin, ymax = latitudes[0], latitudes[-1]

    inputData.close()

    os.system(f"./resampleGrid.sh {directory} {inputFile} {outputFile} {variableInput} {xres} {yres} {fillval} {xmin} {ymin} {xmax} {ymax}")

if __name__ == "__main__":
    main()

#python D:/Profils/mjaouen/Documents/alternance/EASME/gitcode/shellfish_and_algae-MODEL/pretreatment/src/resampleGrid.py -i D:/Profils/mjaouen/Documents/alternance/EASME/data/NWS/merged_Temperature_NWS.nc -d I:/work-he/apps/safi/data/NWS -o I:/work-he/apps/safi/data/NWS/output.nc -v thetao -x longitude -y latitude