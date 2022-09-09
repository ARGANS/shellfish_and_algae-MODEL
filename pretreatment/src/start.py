import os
import json

def cleanFinalSlash(value: str) -> str:
    return value[:-1] if value.endswith('/') else value

source_path:str = cleanFinalSlash(os.getenv('INPUT_SOURCE')) 
destination_path:str = cleanFinalSlash(os.getenv('INPUT_DESTINATION')) 

try:
    with open(source_path + '/parameters.json') as file:
        input_parameters:dict = json.load(file)
except:
    raise RuntimeError('Cannot parse input parameters')

dict_dataCmd = input_parameters.get('datasets', None)

if dict_dataCmd is None:
    raise RuntimeError('The datasets key does not exist in the source manifest file')

year = int(input_parameters.get('year'))

for param, dataset_properties in dict_dataCmd.items():
    dir_data = source_path + '/' + param
    dir_data_pretreated = destination_path + '/' + param
    zone = dataset_properties.get('Place')
    file_name = f'{param}{zone}modelNetCDF{year}-01to{year + 1}-01.nc'
    method = dataset_properties["pretreatment"]
    frequency = dataset_properties["frequency"]

    if method == "Copernicus":
        if frequency == "daily":
            os.system(f"./concatenate_copernicus.sh {dir_data} {dir_data_pretreated} {file_name}")
        else: # monthly
            os.system(f"ln -s {dir_data}/* {dir_data_pretreated}/{file_name}")


    elif method == "Copernicus_Arctic":
        if frequency == "daily":
            os.system(f"./concatenate_copernicus.sh {dir_data} {dir_data_pretreated} conc.nc")
        else: # monthly
            os.system(f"ln -s {dir_data}/* {dir_data_pretreated}/conc.nc")

        os.system(f"./resample_Arctic.sh {dir_data_pretreated} conc.nc {file_name} {dataset_properties['variable']}")
        os.system(f"rm {dir_data_pretreated}/conc.nc")

    elif method == "Reference":
        # Link to a reference file that is already stored
        os.system(f"ln -s {dataset_properties['product-id']} {dir_data_pretreated}/{file_name}")

    # For when we add this step for Hermes or NASA data
    #elif line["pretreatment"] == "Hermes":
    #    ...

    else: # for the data that requires no pretreatment
        os.system(f"ln -s {dir_data}/* {dir_data_pretreated}/{file_name}")
