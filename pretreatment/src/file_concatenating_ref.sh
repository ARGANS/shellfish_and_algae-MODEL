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


#### Averages for masks
# Done in two steps for memory issues, could be done with just ncwa
# use '-v PAR_mean' to extract only the interest variable
# average over time first
ncra -d depth,0.,10. merged.nc averaged.nc
# Then average over depth
ncwa -O -a depth averaged.nc averaged.nc
# TODO: weights ?

### For PAR ocean color
# In combination with some commands from the hermes PAR process

# Add a time dimension like for hermes PAR
# clean filename
file_clean=${filename##*/}
# get the starting and ending day from the file name
starting=${file_clean:5:3}
ending=${file_clean:12:3}
# Compute the time at the middle of the month, in "days since 2020(or other)-01-01 at 00:00:00"
# it requires "bc" to be installed (with apt-get)
day=$(echo "($starting+$ending-1)/2" | bc -l)
# add the computed day
ncap2 -O -s 'time[time]='$day test.nc test.nc
# concatenate like hermes
# reverse lat like hermes
# add units attribute to time
ncatted -O -a units,time,c,c,"days since 2020-01-01 00:00:00" test.nc
# Cut data over 63°N in January, over 68° in February
ncap2 -O -s '*par_jan=par(0,:,:); where(lat>63) par_jan=par@_FillValue; par(0,:,:)=par_jan; ram_delete(par_jan)' test.nc test.nc
ncap2 -O -s '*par_feb=par(1,:,:); where(lat>68) par_feb=par@_FillValue; par(1,:,:)=par_feb; ram_delete(par_feb)' test.nc test.nc


### Arctic resampling
# First concatenate files if necessary (conc.nc in the following)
# Translate the desired variable to GeoTIFF
gdal_translate -ot Float64 NETCDF:conc.nc:temperature -unscale file_Stereo.tiff
# resampling
gdalwarp -s_srs "+proj=stere +lat_0=90 +lat_ts=90 +lon_0=-45 +k=0.00001 +x_0=0 +y_0=0 +a=6378273 +b=6378273  +units=m +no_defs" file_Stereo.tiff -t_srs EPSG:4326  "file_EPSG.tiff" -overwrite
##Function to get the value of a variable attribute
#function ncattget { ncks --trd -M -m ${3} | grep -E -i "^${2} attribute [0-9]+: ${1}" | cut -f 11- -d ' ' | sort ; }
##get the original fill value of the variable
#fillval=$(ncattget _FillValue temperature conc.nc)
# Translate back to NetCDF
gdal_translate -of NetCDF file_EPSG.tiff final.nc