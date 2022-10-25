from osgeo import gdal
from osgeo import osr
from code.models.TiffImage import TiffImage
import numpy as np

def read_tiff(path:str) -> 'TiffImage':
    image = gdal.Open(path, gdal.GA_ReadOnly)
    assert image.RasterCount > 0
    
    nbcol=image.RasterXSize
    nblig=image.RasterYSize
    
    band = image.GetRasterBand(1)
    geotransform = image.GetGeoTransform()
    assert geotransform is not None
    
    # TODO
    # zee=rotate(read_tiff(ficzee, geo=geo), 7)
    # 270 deg.
    # https://www.l3harrisgeospatial.com/docs/rotate.html

    lonmin = geotransform[0]
    latmax = geotransform[3]
    pasx = geotransform[1]
    pasy = abs(geotransform[5])
    lonmax = lonmin + (nbcol * pasx)
    latmin = latmax - (nblig * pasy)

    return TiffImage(
        lonmin=lonmin,
        lonmax=lonmax,
        latmin=latmin,
        latmax=latmax,
        pasx=pasx,
        pasy=pasy
    )

def read(path:str):
    dataset = gdal.Open(path, gdal.GA_ReadOnly)
    assert dataset.RasterCount > 0
    band = dataset.GetRasterBand(1)
    return dataset, band, band.ReadAsArray()

def printBand(band):
    # ds.GetProjection()
    print(
        'Xsize ', band.XSize,
        'Ysize ', band.YSize, 
        'DataType ', band.DataType,
        'NoDataValue ', band.GetNoDataValue(),
        'Max ', band.GetMaximum(),
        'Min ', band.GetMinimum(),
        'RasterMinMax ', band.ComputeRasterMinMax()
    )

def createRaster(array, output_path, gdal_data_type, geotransform, no_data):
    driver = gdal.GetDriverByName(r'GTiff')
    rows, columns = array.shape
    print('ROWS: %s\n COLUMNS: %s\n' % (rows, columns))
    raster = driver.Create(output_path, int(columns), int(rows), 1, eType = gdal_data_type)

    spatial_reference = osr.SpatialReference()
    spatial_reference.ImportFromEPSG(4326)
    raster.SetProjection(spatial_reference.ExportToWkt())
    raster.SetGeoTransform(geotransform)
    output_band = raster.GetRasterBand(1)
    output_band.SetNoDataValue(no_data)
    output_band.WriteArray(array)          
    output_band.FlushCache()
    output_band.ComputeStatistics(False)


def fusionDesFichiers(path1:str, path2: str, outpath:str):
    # ; Fusion des fichiers
    dataset1, band1, array1 = read(path1)
    dataset2, band2, array2 = read(path2)
    
    print('Array')
    print(array1, len(array1))
    array = array1
    for y in range(band1.YSize):
        # if y % 10 == 0:
            # print(f'Y {y}')
        for x in range(band1.XSize):
            if array2[y][x] > 0:
                array[y][x] = array2[y][x]

            if array2[y][x] > 0 and array1[y][x] > 0:
                array[y][x] = (array2[y][x] + array1[y][x]) / 2
    
    print('Completed')    
    createRaster(
        array, 
        outpath,
        band1.DataType, 
        dataset1.GetGeoTransform(), 
        band1.GetNoDataValue()
    )
    
def fusion_zones2(zee_path:str, tmp_path:str, zee_transformed_path:str):
    dataset_zee, band_zee, array_zee = read(zee_path)
    dataset_tmp, band_tmp, array_tmp = read(tmp_path)

    for y in range(band_tmp.YSize):
        # if y % 10 == 0:
        #     print(f'Y {y}')
        for x in range(band_tmp.XSize):
            if array_zee[y][x] == 0:
                array_tmp[y][x] = band_tmp.GetNoDataValue()

    createRaster(
        array_tmp, 
        zee_transformed_path,
        band_tmp.DataType, 
        dataset_tmp.GetGeoTransform(), 
        band_tmp.GetNoDataValue()
    )


"""
zee=rotate(read_tiff(ficzee, geo=geo), 7)
fw=rotate(read_tiff(fictmp, geo=geo), 7)

ind=where(zee EQ 0, ct)
IF (ct GT 0) THEN fw(ind) = -999.

fictmp='../TMP/ZEE_'+espece+'_'+param+'_'+scenario+'_1km.tif'
write_tiff, fictmp, rotate(fw, 7), geo=geo, /float

ficout=dirout+'ZEE_'+espece+'_'+param+'_'+scenario+'_1km.tif'
cmd='gdal_translate -co compress=deflate ' +fictmp+' '+ficout
spawn, cmd
spawn,'rm -f '+fictmp


"""

# ; Fusion des fichiers
# ; -------------------
# r=r1
# ind=where(r2 GT 0, ct)
# IF (ct GT 0) THEN r(ind)=r2(ind)

# ind=where(r1 GT 0 and R2 GT 0, ct)
# IF (ct GT 0) THEN BEGIN
# r(ind)=(r1(ind)+r2(ind))/2.
# ENDIF
