import argparse
import os

import netCDF4 as nc


def main():
    # Create argument parser
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, help="The path to the input file")
    ap.add_argument("-d", "--dir", required=True, help="The path to the working directory")
    ap.add_argument("-o", "--output", required=True, help="The path to the output file")
    ap.add_argument("-v", "--variable", required=True, help="The variable")
    ap.add_argument("-x", "--lonName", required=True, help="The longitude name")
    ap.add_argument("-y", "--latName", required=True, help="The latitude name")
    args = vars(ap.parse_args())

    # Parse arguments
    print("Parsing arguments...")
    inputFile = args["input"]
    outputFile = args["output"]
    directory = args["dir"]
    variable = args["variable"]
    lonName = args["lonName"]
    latName = args["latName"]

    inputData = nc.Dataset(inputFile)

    fillval = str(inputData.variables[variable]._FillValue)
    longitudes = inputData.variables[lonName]
    latitudes = inputData.variables[latName]
    xres = longitudes[1] - longitudes[0]
    yres = latitudes[1] - latitudes[0]

    os.system(f"./resampleGrid.sh {directory} {inputFile} {outputFile} {variable} {xres} {yres} {fillval} {longitudes[0]} {latitudes[0]} {longitudes[-1]} {latitudes[-1]}")

if __name__ == "__main__":
    main()

#python D:/Profils/mjaouen/Documents/alternance/EASME/gitcode/shellfish_and_algae-MODEL/pretreatment/src/resampleGrid.py -i D:/Profils/mjaouen/Documents/alternance/EASME/data/NWS/merged_Temperature_NWS.nc -d I:/work-he/apps/safi/data/NWS -o I:/work-he/apps/safi/data/NWS/output.nc -v thetao -x longitude -y latitude