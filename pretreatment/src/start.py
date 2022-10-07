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


################################################################################
### Each dataset source is treated with the relevant script to get them to a
### uniform format.
### Format expectations:
###     - One file per variable containing data over lat/lon/time/(depth)
###     - latitude and longitude are 1D axes
###     - latitude is sorted in ascending order (as are other axes)
################################################################################



for param, dataset_properties in dict_dataCmd.items():
    dir_data = source_path + '/' + param
    dir_data_pretreated = destination_path + '/' + param
    zone = dataset_properties.get('Place')
    file_name = f'{param}{zone}modelNetCDF{year}-01to{int(year) + 1}-01.nc'
    method = dataset_properties["pretreatment"]
    frequency = dataset_properties["frequency"]
    print(param)

    if method == "Reference":
        print('ok')
        # Copy a reference file that is already stored.
        # Should not use symlink because it we change the file later.
        os.system(f"cp {dataset_properties['product-id']} {dir_data_pretreated}/{file_name}")

    else: # for the data that requires no pretreatment
        os.system(f"ln -s {dir_data}/* {dir_data_pretreated}/{file_name}")


################################################################################
### All datasets are resampled to the grid of one of the variables
################################################################################

reference_param = 'Nitrate' #TODO: make into a parameter
reference_lonName = dict_dataCmd[reference_param]['longName']
reference_latName = dict_dataCmd[reference_param]['latName']
reference_varName = dict_dataCmd[reference_param]['variable']
reference_zone = dict_dataCmd[reference_param]['Place'] #should be the same for all params anyway
file_reference = f'{reference_param}{reference_zone}modelNetCDF{year}-01to{int(year) + 1}-01.nc'
file_reference_full = f'{destination_path}/{reference_param}/{file_reference}'

for param, dataset_properties in dict_dataCmd.items():

    # wait to do the reference param last
    if param == reference_param:
        continue

    dir_data_pretreated = destination_path + '/' + param
    zone = dataset_properties.get('Place')
    file_name = f'{param}{zone}modelNetCDF{year}-01to{int(year) + 1}-01.nc'

    os.system(f'python resampleGrid.py --dir {dir_data_pretreated}'
                                    f' --input {file_name}'
                                    f' --output {file_name}'
                                    f' --variableInput {dataset_properties.get("variable")}'
                                    f' --inputGrid {file_reference_full}'
                                    f' --variableGrid {reference_varName}'
                                    f' --lonName {reference_lonName}'
                                    f' --latName {reference_latName}')

# Now do the reference. TODO: there probably is a better way to do these in the correct order
os.system(f'python resampleGrid.py --dir {destination_path + "/" + reference_param}'
                                f' --input {file_reference}'
                                f' --output {file_reference}'
                                f' --variableInput {reference_varName}'
                                f' --inputGrid {file_reference_full}'
                                f' --variableGrid {reference_varName}'
                                f' --lonName {reference_lonName}'
                                f' --latName {reference_latName}')
