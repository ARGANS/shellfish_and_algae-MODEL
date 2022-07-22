import os
from models.ModelProperties import ModelProperties

model_properties = ModelProperties(os.getenv('DATASET_ID'), os.getenv('TASK_ID'))
try:
    model_properties.parse(os.getenv('PARAMETERS_JSON'))
except:
    raise RuntimeError('Cannot parse the value of the parameters_json environment variable')

if not model_properties.isDataDownloadTaskCompleted():
    raise RuntimeError('Data not downloaded')

# TODO: how to get this ?
dict_dataCmd = model_properties.<something to get the dict>

for param, line in dict_dataCmd.items():
    # TODO: get these dir names (no last /)
    dir_data = 
    dir_data_pretreated =
    # TODO: I think the template contains the /media/share/... part, this should just be the file name
    file_name = model_properties.file_template.format(**line)

    if line["pretreatment"] == "Copernicus":
        if line["frequency"] == "daily":
            os.system(f"./concatenate_copernicus.sh {dir_data} {dir_data_pretreated} {file_name}")
        else: # monthly
            os.system(f"ln -s {dir_data}/* {dir_data_pretreated}/{file_name}")


    elif line["pretreatment"] == "Copernicus_Arctic":
        if line["frequency"] == "daily":
            os.system(f"./concatenate_copernicus.sh {dir_data} {dir_data_pretreated} conc.nc")
        else: # monthly
            os.system(f"ln -s {dir_data}/* {dir_data_pretreated}/conc.nc")

        os.system(f"./resample_Arctic.sh {dir_data_pretreated} conc.nc {file_name} {line['variable']}")
        os.system(f"rm {dir_data_pretreated}/conc.nc")

    # For when we add this step for Hermes or NASA data
    #elif line["pretreatment"] == "Hermes":
    #    ...

    else: # for the data that requires no pretreatment
        os.system(f"ln -s {dir_data}/* {dir_data_pretreated}/{file_name}")