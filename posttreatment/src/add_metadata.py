import argparse
import json
import netCDF4 as nc
from datetime import datetime
import numpy as np
from osgeo import gdal

def import_json(path: str) -> dict:
    with open(path, 'r') as f:
        data = json.loads(f.read())
    return data 

def tiffToArray(tiffPath):
    dataset = gdal.Open(tiffPath)
    im = dataset.ReadAsArray()
    return np.array(im)

def main():
    '''
    This script add metadata to the input file
    '''

    # Create argument parser
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, help="The path to the input file")
    ap.add_argument("-j", "--jsonFile", required=True, help="The path to the parameters json")
    ap.add_argument("-n", "--name", required=True, help="Data name")
    ap.add_argument("-z", "--zeeFile", required=True, help="The path to the zee mask file")

    args = vars(ap.parse_args())

    ncFile = args["input"]
    jsonFile = args["jsonFile"]
    name = args["name"]
    zeeFile = args["zeeFile"]
    json_data = import_json(jsonFile)

    #print(ncFile)
    

    ds_append = nc.Dataset(ncFile, 'a')

    if len(np.shape(ds_append[name]))==3:
        data_array = ds_append[name][0, :, :]
        zeeMask_array = tiffToArray(zeeFile)
        data_array=np.ma.masked_where(zeeMask_array[::-1]<1,data_array)
        ds_append[name][0, :, :] = data_array
    else:
        data_array = ds_append[name][:, :]
        zeeMask_array = tiffToArray(zeeFile)
        data_array=np.ma.masked_where(zeeMask_array[::-1]<1,data_array)
        ds_append[name][:, :] = data_array

    ds_append.project = 'Studies to support the European Green Deal - Lot 1 Shellfish and algae'
    ds_append.institution = 'ARGANS-FR, Bantry Marine Research Station (BMRS), Cofrepeche'
    ds_append.production = 'ARGANS-FR E-mails: contact@argans.eu'
    ds_append.Author_email = 'contact@argans.eu'
    ds_append.creation_time = datetime.now().strftime("%d/%m/%Y, %H:%M:%S")
    ds_append.spatial_resolution = '1km'
    ds_append.source = 'CMEMS models, EMODNET, NASA/OCEANCOLOR'
    ds_append.title = str(name) + ' ' + str(json_data["metadata"]["name"])
    ds_append.image_type = 'composite'
    ds_append.image_reference_date = str(json_data["dataset_parameters"]["year"])
    ds_append.area = json_data["metadata"]["zone"]

    ds_append.close()

if __name__ == "__main__":
    main()