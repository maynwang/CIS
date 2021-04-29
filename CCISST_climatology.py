'''
    Make CCISST climatology
'''

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import mpl_toolkits.basemap as bm
from netCDF4 import Dataset
from scipy import io
import os
import ecoliver as ecj

# File location
dataheader = '/home/oliver/data/sst/CoralTemp/'
path = '/home/mwang/results/'

# Load in coordinates
nc = Dataset(dataheader + '1985/coraltemp_v1.0_19850101.nc')
lon = nc.variables['lon'][:]
lat = nc.variables['lat'][:]
nc.close()

# Domain of interest
lon1 = -67.
lon2 = -51.
lat1 = 49.
lat2 = 63.
ii = (lon>=lon1) * (lon<=lon2)
jj = (lat>=lat1) * (lat<=lat2)
lon = lon[ii]
lat = lat[jj]
llon, llat = np.meshgrid(lon, lat)
X = len(lon)
Y = len(lat)
M = 12 # Number of months

#
# Time series
#

# Time and dates
t, dates, T, year, month, day, doy = ecj.timevector([1985,1,1], [2019,12,31])

# Initialise fields
sst = np.zeros((Y, X, T))
c = np.zeros((Y, X, T))

# Load in each year and contribute to average
for tt in range(T):
    print(dates[tt])
    nc = Dataset(dataheader + str(year[tt]) + '/coraltemp_v1.0_' + str(year[tt]) + str(month[tt]).zfill(2) + str(day[tt]).zfill(2) + '.nc')
    SST = nc.variables['analysed_sst'][0,:,ii][jj,:]
    C = nc.variables['sea_ice_fraction'][0,:,ii][jj,:]
    nc.close()
    sst[:,:,tt] = SST
    c[:,:,tt] = C
    #m = month[tt] - 1
    #sst[:,:,m] += SST
    #c[:,:,m] += C
    #N[:,:,m] += 1.

# Mask
sst[sst == sst.min()] = np.nan
c[c<0] = np.nan

# Save data
np.savez(path + 'CCI_Labrador.npz', lon=lon, lat=lat, sst=sst, c=c, t=t, dates=dates)

# Load data -- add time
#data = np.load(dataheader + 'CCI_Labrador.npz')  
#sst = data['sst']
#c = data['c']

#
# Make climatology and trends
#

# Climatology
sst_clim = np.zeros((Y, X, M))
c_clim = np.zeros((Y, X, M))
for m in range(12):
    tt = month == m+1
    sst_clim[:,:,m] = np.nanmean(sst[:,:,tt], axis=2)
    c_clim[:,:,m] = np.nanmean(c[:,:,tt], axis=2)

# Trends
sst_tr = np.zeros((Y, X, M))
c_tr = np.zeros((Y, X, M))
for j in range(Y):
    print(j+1, Y)
    for i in range(X):
        for m in range(12):
            tt = month == m+1
            sst_tr[j,i,m] = ecj.trend(t[tt], sst[j,i,tt])[1]*365.25*10 # deg C / decade
            c_tr[j,i,m] = ecj.trend(t[tt], c[j,i,tt])[1]*365.25*10 # fraction / decade

# Differences
tt1 = (year >= 1985) * (year <= 1999)
tt2 = (year >= 2005) * (year <= 2019)
sst_diff = np.zeros((Y, X, M))
c_diff = np.zeros((Y, X, M))
for m in range(12):
    ttm = month == m+1
    tt1m = tt1 * ttm
    tt2m = tt2 * ttm
    sst_diff[:,:,m] = np.nanmean(sst[:,:,tt2m], axis=2) - np.nanmean(sst[:,:,tt1m], axis=2)
    c_diff[:,:,m] = np.nanmean(c[:,:,tt2m], axis=2) - np.nanmean(c[:,:,tt1m], axis=2)

# Apply mask to ice concentration
for m in range(12):
    tmp = c_clim[:,:,m]
    tmp[np.isnan(sst_clim[:,:,m])] = np.nan
    c_clim[:,:,m] = tmp
    tmp = c_tr[:,:,m]
    tmp[np.isnan(sst_tr[:,:,m])] = np.nan
    c_tr[:,:,m] = tmp
    tmp = c_diff[:,:,m]
    tmp[np.isnan(sst_diff[:,:,m])] = np.nan
    c_diff[:,:,m] = tmp

#c[c == c.min()] = np.nan
#c[c < 0] = np.nan

np.savez(path + 'CCI_Climatology_1985_2019.npz', lon=lon, lat=lat, sst_clim=sst_clim, c_clim=c_clim, sst_tr=sst_tr, c_tr=c_tr, sst_diff=sst_diff, c_diff=c_diff)

