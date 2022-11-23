import argparse
import json
import netCDF4 as nc
from datetime import datetime

def import_json(path: str) -> dict:
    with open(path, 'r') as f:
        data = json.loads(f.read())
    return data 

def main():
    '''
    This script add metadata to the input file
    '''

    # Create argument parser
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, help="The path to the input file")
    ap.add_argument("-j", "--jsonFile", required=True, help="The path to the parameters json")
    ap.add_argument("-n", "--name", required=True, help="Data name")

    args = vars(ap.parse_args())

    ncFile = args["input"]
    jsonFile = args["jsonFile"]
    name = args["name"]
    json_data = import_json(jsonFile)

    #print(ncFile)

    ds_append = nc.Dataset(ncFile, 'a')

    '''print(ds_append)
    for name in rootgrp.ncattrs():
        print("Global attr {} = {}".format(name, getattr(rootgrp, name)))'''

    ds_append.project = 'Studies to support the European Grenn Deal - Lot 1 Shellfish and algae: http://www.??????.com'
    ds_append.institution = 'ARGANS-FR, Bantry Marine Research Station (BMRS), Cofrepeche'
    ds_append.production = 'ARGANS-FR E-mails: contact@argans.eu'
    ds_append.Author_email = 'contact@argans.eu'
    ds_append.creation_time = datetime.now().strftime("%m/%d/%Y, %H:%M:%S")
    ds_append.spatial_resolution = '1km'
    ds_append.source = 'CMEMS models, EMODNET, NASA/OCEANCOLOR'
    ds_append.title = str(name) + ' ' + str(json_data["metadata"]["name"])
    ds_append.image_type = 'composite'
    ds_append.image_reference_date = str(json_data["dataset_parameters"]["year"])
    ds_append.area = json_data["metadata"]["zone"]

    print(ds_append)

    '''for name in rootgrp.ncattrs():
        print("Global attr {} = {}".format(name, getattr(rootgrp, name)))'''
    ds_append.close()

if __name__ == "__main__":
    main()