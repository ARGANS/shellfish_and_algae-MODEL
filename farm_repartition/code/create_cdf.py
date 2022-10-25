from code.utils.cmd import run
from code.utils.tiff import read, read_tiff, printBand
from code.utils.ncdf_date import ncdf_date, hours_since_1900
import os
import netCDF4 as nc


def create_cdf(ficin, fictmp, ficout, espece, libespece, nomparam, libparam, scenario, unite, libzone, annee, title, spatial_resolution, source):
    print('create_cdf', '\n')
    # fictmp = './tmp/' + os.path.basename(ficin)
    run('gdal_translate -co compress=none {ficin} {fictmp}'.format(ficin=ficin, fictmp=fictmp))
    dataset1, band1, array1 = read(fictmp)
    tiffImage = read_tiff(fictmp)
    nb_col=dataset1.RasterXSize
    nb_lig=dataset1.RasterYSize
    lon = [(x*tiffImage.pasx + tiffImage.lonmin + tiffImage.pasx / 2) for x in range(nb_col)]
    lat = [(x*tiffImage.pasy + tiffImage.latmin + tiffImage.pasy / 2) for x in range(nb_lig)]
    n_hours = hours_since_1900(annee)
    # ; Date de creation
    # ; ----------------
    datejour=ncdf_date()
    # ; Valeur du masque
    # ; ----------------
    mask_value = band1.GetNoDataValue()
    # ; ==========================
    # ; Creation du fichier NETCDF
    # ; ==========================
    # spawn, '\rm -f '+ficout+'*'
    # dirout=file_dirname(ficout)
    # file_mkdir, dirout
    ds = nc.Dataset(ficout, 'w', format='NETCDF4', clobber=True)

    # Write the global attributes
    ds.setncattr('Conventions', 'CF-1.0')
    ds.setncattr('netcdf_version_id', '3.4')
    ds.setncattr('project', 'Studies to support the European Green Deal - Lot 1 Shellfish and algae: http://www.??????.com/')
    ds.setncattr('institution', 'ARGANS-FR, Bantry Marine reseach Station (BMRS), Cofrepeche')
    ds.setncattr('production', 'ARGANS-FR E-mails: contact@argans.eu')
    ds.setncattr('WEB_visualisation','http://www.????.com"')
    ds.setncattr('Author_e-mail', 'contact@argans')
    ds.setncattr('creation_time', datejour)
    ds.setncattr('title', title)
    ds.setncattr('file_name', ficout)
    ds.setncattr('spatial_resolution', spatial_resolution)
    ds.setncattr('source', source)
    ds.setncattr('image_type', 'composite')
    ds.setncattr('image_reference_date', annee)
    ds.setncattr('southernmost_latitude', tiffImage.latmin)
    ds.setncattr('northernmost_latitude', tiffImage.latmax)
    ds.setncattr('westernmost_longitude', tiffImage.lonmin)
    ds.setncattr('easternmost_longitude', tiffImage.lonmax)
    ds.setncattr('area', libzone)
    ds.setncattr('product_version', '1.0')

    # ; Time
    # ; ----
    time_id = ds.createDimension('time', 1)
    tab_time = ds.createVariable('time', np.int32, ('time',))
    tab_time = n_hours
    tab_time.setncattr('long_name','time')
    tab_time.setncattr('units', 'hours since 1900-1-1 0:0:0')

    # ; Latitude
    # ; --------
    latitude_id = ds.createDimension('lat', nb_lig)
    tab_lat = ds.createVariable('time', np.float32, ('lat',))
    tab_lat[:] = lat
    tab_lat.setncattr('long_name', 'latitude')
    tab_lat.setncattr('units', 'degrees_north')
    tab_lat.setncattr('standard_name', 'latitude')

    # ; Longitude
    # ; ---------
    longitude_id = ds.createDimension('lon', nb_col)
    tab_lon = ds.createVariable('time', np.float32, ('lon',))
    tab_lon[:] = lon
    tab_lon.setncattr('long_name', 'longitude')
    tab_lon.setncattr('units', 'degrees_east')
    tab_lon.setncattr('standard_name', 'longitude')
    #-------------------------------------------------------
    # ; Parametre
    # ; ---------
    long_name=libparam
    valid_min, valid_max = band1.ComputeRasterMinMax()
    #-------------------------------------------------------
    param_cdf = array1

    tab_param = ds.createVariable(nomparam, np.float64, ('lon', 'lat', 'time'))
    tab_param[:] = param_cdf
    tab_lon.setncattr('long_name', long_name)
    tab_lon.setncattr('units', unite)
    tab_lon.setncattr('valid_min', valid_min) # /FLOAT
    tab_lon.setncattr('valid_max', valid_max) # /FLOAT
    tab_lon.setncattr('_FillValue', mask_value) # /FLOAT
    tab_lon.setncattr('missing_value', mask_value) # /FLOAT


    #-------------------------------------------------------
    # ; Fermeture du fichier NETCDF
    # ; ---------------------------
    ds.close()
    run('rm -f ' + fictmp)
