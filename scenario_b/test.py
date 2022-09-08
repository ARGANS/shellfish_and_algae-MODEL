import datetime
from dataread.read_netcdf import ParamData

param_metadata ={ 
    'Nitrate': {
        'depth_name': 'depth',
        'file_name': '/media/share/data/00c539642f66a1856e637249744b6b32/_pretreated/Nitrate/NitrateNWSmodelNetCDF2021-01to2022-01.nc',
        'latitude_name': 'latitude',
        'longitude_name': 'longitude',
        'time_name': 'time',
        'time_step': datetime.timedelta(seconds=3600),
        'time_zero': datetime.datetime(1950, 1, 1, 0, 0),
        'unit_conversion': 14.0,
        'variable_name': 'no3'},
    'PAR': {
        'depth_name': 'depth',
        'file_name': '/media/global/PAR_ref_OC_2020.nc',
        'latitude_name': 'lat',
        'longitude_name': 'lon',
        'time_name': 'time',
        'time_step': datetime.timedelta(days=1),
        'time_zero': datetime.datetime(2021, 1, 1, 0, 0),
        'unit_conversion': 11.574,
        'variable_name': 'par'},
    'Temperature': {
        'depth_name': 'depth',
        'file_name': '/media/share/data/00c539642f66a1856e637249744b6b32/_pretreated/Temperature/TemperatureNWSmodelNetCDF2021-01to2022-01.nc',
        'latitude_name': 'latitude',
        'longitude_name': 'longitude',
        'time_name': 'time',
        'time_step': datetime.timedelta(seconds=3600),
        'time_zero': datetime.datetime(1950, 1, 1, 0, 0),
        'unit_conversion': 1.0,
        'variable_name': 'thetao'},
    'eastward_Water_current': {
        'depth_name': 'depth',
        'file_name': '/media/share/data/00c539642f66a1856e637249744b6b32/_pretreated/eastward_Water_current/eastward_Water_currentNWSmodelNetCDF2021-01to2022-01.nc',
        'latitude_name': 'latitude',
        'longitude_name': 'longitude',
        'time_name': 'time',
        'time_step': datetime.timedelta(seconds=3600),
        'time_zero': datetime.datetime(1950, 1, 1, 0, 0),
        'unit_conversion': 86400.0,
        'variable_name': 'uo'},
    'northward_Water_current': {
        'depth_name': 'depth',
        'file_name': '/media/share/data/00c539642f66a1856e637249744b6b32/_pretreated/northward_Water_current/northward_Water_currentNWSmodelNetCDF2021-01to2022-01.nc',
        'latitude_name': 'latitude',
        'longitude_name': 'longitude',
        'time_name': 'time',
        'time_step': datetime.timedelta(seconds=3600),
        'time_zero': datetime.datetime(1950, 1, 1, 0, 0),
        'unit_conversion': 86400.0,
        'variable_name': 'vo'
    }
}

sim_area = {'averagingDims': ('depth',),
 'depth': (0, 2.8),
 'latitude': (-90, 90),
 'longitude': (-180, 180),
 'weighted': False}

paramData = ParamData(**param_metadata['eastward_Water_current'])
variableName = paramData._variableName

print('paramData')
print(paramData)

print(f'\nDS {paramData.ds.__class__} variableName: {variableName}')
print(paramData.ds)

kwargs = {'lon': (-180, 180), 'lat': (-90, 90), 'depth': (0, 2.8)} 
variable = paramData.ds[variableName]
dimensions = variable.get_dims()
print(f'\nDimensions: {dimensions}')




if False:
    print('\nStart extraction')
    if False:
        # causes an OOM exception
        sliceList = [slice(None)] * len(dimensions)
        extract = variable[tuple(sliceList)]
    else:
        ## sliceList ==(slice(None, None, None), slice(0, 1, None), slice(0, 1240, None), slice(0, 958, None))
        extract = variable[(slice(0, 1, None), slice(0, 1240, None), slice(0, 958, None), slice(None, None, None))]
    print('Extract')
    print(extract)


variable = paramData.getVariable(**sim_area)
# print('variable')
# print(variable)
print('end')
