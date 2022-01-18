# WFS = vector
# import the packages
import pandas
import sys
import os
from owslib.fes import *
from owslib.etree import etree
from owslib.wfs import WebFeatureService

fn = r'bathymetry.xml'
wfs = 'https://ows.emodnet-seabedhabitats.eu/geoserver/emodnet_open/wfs'

# define the link to the service and carry out the request
wfs = WebFeatureService(url=wfs, version='1.1.0')#, timeout=320)
print('layer list:')
print(list(wfs.contents))
layername = 'emodnet_bathymetry'

print('To choose srs')
print(wfs[layername].crsOptions)
print('To choose the styles')
print(wfs[layername].styles)
print('To choose bbox')
print(wfs[layername].boundingBoxWGS84)
print('To choose format')
print(wfs.getOperationByName('GetFeature').formatOptions)
#requestbbox = (lonWest, latSouth, lonEast, latNorth)
requestbbox = (-17.1586666674327, 32.7916669909389, 34.1019547693972, 63.9904999997752)


response = wfs.getfeature(typename=layername,bbox=requestbbox, srsname='urn:x-ogc:def:crs:EPSG:4326')
if os.path.isfile(fn):
    os.unlink(fn)
out = open(fn, 'wb')
out.write(response.read())
out.close()


# define pandas dataframe for the output
#df = pandas.read_csv(fn, header=0, sep=',')
#print(df)