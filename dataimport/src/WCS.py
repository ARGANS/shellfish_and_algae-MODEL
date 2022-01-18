# import packages
import os
from owslib.wcs import WebCoverageService
from osgeo import gdal

# define the connection
url = 'http://ows.emodnet-bathymetry.eu/wcs?'
wcs = WebCoverageService(url, version='1.0.0', timeout=320)
#wcs.identification.type
#wcs.identification.title

# define variables
clipfile = r'temp.tif'
#requestbbox = (lonOuest, latSud, lonEst, latNord)
requestbbox = (-4.7, 48.55, -4.50, 48.65)
#requestbbox = (2.0, 51.5, 5.0, 54.0)
layer = 'emodnet:mean'

# get the data
bathyml = 'emodnet:mean'
sed = wcs[layer]  # this is necessary to get essential metadata from the layers advertised
cx, cy = map(int, sed.grid.highlimits)
bbox = sed.boundingboxes[0]['bbox']
lx, ly, hx, hy = map(float, bbox)
resx, resy = (hx - lx) / cx, (hy - ly) / cy
# width = cx/1000
# height = cy/1000

gc = wcs.getCoverage(identifier=bathyml,
                     bbox=requestbbox,
                     coverage=sed,
                     format='GeoTIFF',
                     crs=sed.boundingboxes[0]['nativeSrs'], resx=resx, resy=resy)


fn = clipfile
f = open(fn, 'wb')
f.write(gc.read())
f.close()

#filetiff = gdal.Open(clipfile)
#theImage = filetiff.GetRasterBand(1).ReadAsArray()
#os.unlink(clipfile)
