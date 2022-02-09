# Commands used to concatenate netcdf input files into single files.
# To be made into a proper script ultimately

### For Copernicus data:
# On first file, make 'time' into the record dimension:
ncks -O --mk_rec_dim time in1.nc in1.nc
# Then concatenate all files:
ncrcat -h in*.nc out.nc
# Finally change 'valid_max' attribute to the adequate value (TODO: determining this value):
ncatted -O -a valid_max,time,o,d,622788. out.nc


### For PAR:
# Create a 'time' dimension and make PAR_mean and associates depend on it,
# do so for every file:
ncecat -O -u time in1.nc in1.nc
ncwa -O -a time in1.nc in1.nc
ncecat -O -u time in1.nc in1.nc
# Note: I don't know how this works ?
# Then concatenate all files:
ncrcat -h in*.nc out.nc
# TODO: 'time' is a dimension with no variable, it needs to store time information
# Assuming we have data every day this could work, the starting time needs to be determined:
ncap2 -O -s 'time=array(621234.,24,$time)' out.nc out.nc
# Note: This also assumes that the alphabetical order of files corresponds to their order in time
# 'lat' needs to be reversed (ordered) to find adequate coordinates optimally
ncpdq -a '-lat' out.nc out.nc
# Maybe rename lat into latitude and lon into longitude ?