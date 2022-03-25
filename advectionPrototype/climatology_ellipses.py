from time import time
import numpy as np
from numbers import Number
from numpy.lib.function_base import extract
from sklearn.linear_model import LinearRegression
import netCDF4 as nc
import pickle

import time

import sys
sys.path.insert(1, 'P:\\Aquaculture\\shellfish_and_algae-MODEL\\')
from dataread.read_netcdf import *

def degrees_to_meters(lonDist, latDist, refLat):
    """Converts degrees of latDist,lonDist to meters assuming that the latitude
    stays near refLat.
    """
    lat_degree = 111000 # conversion of latitude degrees to meters

    Dx = lonDist * lat_degree * np.cos(np.deg2rad(refLat))
    Dy = latDist * lat_degree

    return Dx, Dy

def meters_to_degrees(Dx, Dy, refLat):
    """Converts longitudinal and meridional distances Dx and Dy to degrees of
    longitude and latitude assuming that the latitude stays near refLat"""

    lat_degree = 111000 # conversion of latitude degrees to meters

    lonDist = Dx / (lat_degree * np.cos(np.deg2rad(refLat)))
    latDist = Dy / lat_degree

    return lonDist, latDist


class Ellipse:

    def __init__(self, eastward_current: np.array, northward_current: np.array, 
                 lon: Number, lat: Number, time_factor: Number):

        Vx, Vy = eastward_current*time_factor, northward_current*time_factor

        self._origin = np.array((lon, lat))

        # center defined by averages
        self._center = np.array((Vx.mean(), Vy.mean()))
        #print(self._center)

        # angle defined by linear fit
        self._angle = np.arctan(LinearRegression().fit(Vx.reshape((-1, 1)), Vy).coef_[0])
        #print(self._angle)

        # rotation matrix to the base of the ellipse
        self._rotation = np.array([[np.cos(self._angle), np.sin(self._angle)],
                                   [-np.sin(self._angle), np.cos(self._angle)]])
        #print(self._rotation)

        currents_input = np.array([Vx, Vy])
        #print(currents_input)
        currents_rotated = np.matmul(self._rotation, currents_input)
        #print(currents_rotated)
        # radii defined by 2 sigma or around the 95th percentile in a gaussian
        self._radius = 2 * np.std(currents_rotated, 1)
        #print(self._radius)

    def __str__(self):
        sentence = f"Ellipse object with parameters:\n" + \
                   f"    .Origin: {self._origin} degrees lon,lat\n" + \
                   f"    .Center: {self._center} meters\n" + \
                   f"    .Angle: {self._angle} radians\n" + \
                   f"    .Radii: {self._radius} meters\n"
        return sentence

    def contains(self, lon: Number, lat: Number):
        """Checks wether the point at lon,lat is in the ellipse
        """

        # Get the position relative to the origin in meters
        position = degrees_to_meters(lonDist = lon - self._origin[0],
                                     latDist = lat - self._origin[1],
                                     refLat = self._origin[1])
        
        # distance in directions a and b from the center of the ellipse
        pos_ab = np.matmul(self._rotation, np.array(position) - self._center)

        return np.sum((pos_ab/self._radius)**2) <= 1
    


if __name__=="__main__":

    path_east = "D:\\IBI\\merged_eastward_Water_current_IBI.nc"
    path_north = "D:\\IBI\\merged_northward_Water_current_IBI.nc"
    path_ellipses = "D:\\IBI\\ellipses.pkl"

    ds_east = nc.Dataset(path_east)
    ds_north = nc.Dataset(path_north)

    current_map, dims = extractVarSubset(ds_east, "uo", time_index=0, depth=3)

    print(dims)
    print(current_map.shape)

    Vx = np.array([1,2,3,4,4,5])
    Vy = np.array([1,2,3,4,-4,5])
    #Vy = np.array([-1,2,-3,-4,4,5])

    test = Ellipse(Vx, Vy, 0, 0, 1)

    print(test)

    
    t0 = time.time()

    ellipses = np.empty(current_map.shape, dtype=object)

    lat = extractVarSubset(ds_east, "latitude", depth=3)[0]
    lon = extractVarSubset(ds_east, "longitude", depth=3)[0]

    """

    for i_lat in range(current_map.shape[0]):
        print(i_lat)
        for j_lon in range(current_map.shape[1]):
            if current_map.mask[i_lat, j_lon]:
                continue

            V_east, _ = extractVarSubset(ds_east, "uo", depth=3, latitude_index=i_lat, longitude_index=j_lon)
            V_north, _ = extractVarSubset(ds_north, "vo", depth=3, latitude_index=i_lat, longitude_index=j_lon)

            ellipses[i_lat, j_lon] = Ellipse(eastward_current=V_east, northward_current=V_north, 
                                             lon=lon[j_lon], lat=lat[i_lat], time_factor=24*3600)
    
    print(time.time() - t0)

    with open(path_ellipses,"wb") as f:
        pickle.dump(ellipses, f)
    
    """
    with open(path_ellipses, "rb") as f:
        ellipses = pickle.load(f)


    excludes_self = np.zeros(ellipses.shape, dtype=bool)

    for i in range(ellipses.shape[0]):
        for j in range(ellipses.shape[1]):
            if ellipses[i,j] is None:
                continue
            excludes_self[i,j] = not ellipses[i,j].contains(lon=lon[j], lat=lat[i])

    print(time.time() -t0)
    print(excludes_self)
    print(np.sum(excludes_self))
