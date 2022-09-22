import netCDF4
import numpy
import xarray
import os

ds = xarray.open_mfdataset(os.path.join('..','Data','IBI','northward_water_current','northward_Water_currentIBImodelNetCDF*.nc'),combine = 'nested', concat_dim="time")
ds.to_netcdf(os.path.join('..','Data_merg','IBI','northward_water_current','northward_water_currentIBImodelNetCDF2020-09to2021-06.nc'))

#ds = xarray.open_mfdataset('..\Data\IBI\Chlorophyll-a\Chlorophyll-aIBImodelNetCDF*.nc',combine = 'nested', concat_dim="time")
#ds.to_netcdf('..\Data_merg\IBI\Chlorophyll-a\Chlorophyll-aIBImodelNetCDF2020-09to2021-06.nc')

#ds = xarray.open_mfdataset('..\Data\IBI\Temperature\TemperatureIBImodelNetCDF*.nc',combine = 'nested', concat_dim="time")
#ds.to_netcdf('..\Data_merg\IBI\Temperature\TemperatureIBImodelNetCDF2020-09to2021-06.nc')

#ds = xarray.open_mfdataset('..\Data\IBI\eastward_water_current\eastward_water_currentIBImodelNetCDF*.nc',combine = 'nested', concat_dim="time")
#ds.to_netcdf('..\Data_merg\IBI\eastward_water_current\eastward_water_currentIBImodelNetCDF2020-09to2021-06.nc')
