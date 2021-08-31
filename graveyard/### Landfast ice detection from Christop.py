### Landfast ice detection from Christoph

# import modules that we need
from pathlib import Path  # makes handling of paths and file names easier
import matplotlib.pyplot as plt
import xarray as xr       # convenient for reading NetCDF files and preserving metadata

import cartopy
import cartopy.crs as ccrs
import cmocean.cm as cmo
import cv2
import numpy as np

def landfast_ice_detection(date_str):

    '''Calculates distance and area of landfast ice. Inputs are date_str, which is the date we wish to analyze'''

    # input data directory
    path = Path('/extra-space1/data/tikoralukupload/cis-weekly/nc/')

    # region of interest
    region = 'HB'

    # Projection for mapping
    rot = ccrs.RotatedPole(pole_latitude=55,pole_longitude=150)

    ### Construct the full file name:
    # Since path is a `pathlib.Path` object, we can create the full path of the input file by 
    # appending the name with a forward slash "/". We also use the `.format()` method to use
    # the variable `region` in the string of the file name.
    fname = path / '{}.nc'.format(region)

    ### Load data into a `xarray.Dataset`:
    # Note that the data is loaded lazily which means that no actual data are actually written
    # into memory.
    ds = xr.open_dataset(fname)

    # subset `xarray.Dataset` by index for HB
    ds = ds.isel(x=slice(1700, None), y=slice(600, 1450))

    # Note that this subsets the entire dataset and therefore, you only have to do it once for all variables.

    # variable name
    vname = 'E_CT'

    # create shortcuts to the variables
    lon = ds.longitude
    lat = ds.latitude
    # juld = ds.juld
    land = ds.land
    E_CT = ds[vname]

    # Note that I chose the bracket notation for the last variable and used a variable `vname`
    # to specify which field we want to look at. This makes the code more flexible as we could
    # define a list of names that we can loop through to run the same analysis systematically
    # for multiple variables.

    ### Use OpenCV to detect contours

    ret, thresh = cv2.threshold(land.values,       # src
                                0.5,               # thresh
                                1,                 # maxval
                                cv2.THRESH_BINARY) # type

    # Find the contours in the array
    contours, hierarchy = cv2.findContours(thresh.astype(np.uint8), # image
                                        cv2.RETR_EXTERNAL,       # mode
                                        cv2.CHAIN_APPROX_NONE)   # method                           

    # pick the largest contour and assume it is surrounding the land
    cland = np.squeeze(max(contours, key=cv2.contourArea))

    # x- and y-indices of the land contour
    xinds = xr.DataArray(cland[:,0], dims=['coastline'])
    yinds = xr.DataArray(cland[:,1], dims=['coastline'])

    # longitude and latitude of the land contour
    lon_cland = lon.isel(x=xinds, y=yinds)
    lat_cland = lat.isel(x=xinds, y=yinds)

    # we manually find the minimum and maximum indexes of the coastline of interest
    indmin = 2087
    indmax = None

    # get longitude and latitude of coastline of interest
    lon_coast = lon_cland.isel(coastline=slice(indmin, indmax))
    lat_coast = lat_cland.isel(coastline=slice(indmin, indmax))

    # quick plot to check selected coastline
    fig, ax = plt.subplots(1, 1, figsize=(7.5, 9), subplot_kw={'projection': rot})
    ax.set_extent([-61.5, -56.5, 53.8, 61.5])
    ax.pcolormesh(lon, lat, land, transform=ccrs.PlateCarree())
    ax.plot(lon_coast, lat_coast, '-' ,transform=ccrs.PlateCarree(),zorder=4)

    # select a day 
    CT = E_CT.sel(juld=date_str)

    # - fill the NaN values with 11, i.e., pretend they are part of the landfast sea ice
    # - select grid points that either have CT = 11 or are land as defined in the
    #   land/sea mask `land`, i.e., not the open water area.
    # - set all other values to zero
    CT_sina = CT.fillna(11).where((CT==11) | (land==1), other=0)



    # create a binary image where all grid points that are not land or have CT = 11 are
    # set to zero. (Technically, we could skip this step, because our input "image" is
    # already binary.)
    ret, thresh = cv2.threshold(CT_sina.values,    # src
                                10,                # thresh
                                11,                # maxval
                                cv2.THRESH_BINARY) # type\

    # find the contours
    contours, hierarchy = cv2.findContours(thresh.astype(np.uint8), # image
                                        cv2.RETR_EXTERNAL,       # mode
                                        cv2.CHAIN_APPROX_NONE)   # method

    # Since the land grid points are part of this contour, we assume again that the
    # contour of interest is the larges one.
    tmp = np.squeeze(max(contours, key=cv2.contourArea))

    # x- and y-indices of the land+sina contour
    xinds = xr.DataArray(tmp[:,0], dims=['coastline'])
    yinds = xr.DataArray(tmp[:,1], dims=['coastline'])

    # longitude and latitude of the land+sina contour
    lon_edge = lon.isel(x=xinds, y=yinds)
    lat_edge = lat.isel(x=xinds, y=yinds)

    # The exclusive interesection between the contour `tmp` and the land contour `cland`
    # is the landfast sea ice "sina"

    # Generate masks of both contours
    msk_tmp = cv2.drawContours(np.zeros(CT.shape, np.uint8), [tmp], 0, 1, cv2.FILLED)
    msk_cland = cv2.drawContours(np.zeros(CT.shape, np.uint8), [cland], 0, 1, cv2.FILLED)

    # Now we find all the grid points that are contained in both masks
    sina = cv2.bitwise_xor(msk_tmp, msk_cland)

    # Create a xr.DataArray with the same dimensions and coordinates as the original data.
    # We are also masking out the land and area where we didnt have any data.
    sina = xr.DataArray(data=sina, dims=CT.dims, coords=CT.coords).where(land==0)

    # In this array all grid points with landfast sea ice are equal to one. Since each grid
    # cell represents an area of 1 km^2, we can get the total area of sina by computing
    # the sum:
    sina_area = sina.sum()

    print('The total area of landfast sea ice is {} km^2'.format(sina_area.values))

    # create plot
    fig, ax = plt.subplots(1, 1, figsize=(7.5, 9), subplot_kw={'projection': rot})
    ax.set_extent([-61.5, -56.5, 53.8, 61.5])
    ax.contourf(lon, lat, sina, transform=ccrs.PlateCarree())

    ### Distance from coast to edge of LF sea ice

    # x-and y coordinates of coastline and landfast ice edge as complex numbers
    a = lon_coast.x + 1j * lon_coast.y 
    b = lon_edge.x + 1j * lon_edge.y

    # a is a "vector" of length N, where N is the number of grid points along the coast.
    # b is a "vector" of length M, where M is the number of grid points along the ice edge.

    # create a MxN matrix where each row is a copy of a
    A = np.array([a,] * len(b))

    # create a NxM matrix where each row is a copy of b
    B = np.array([b,] * len(a)) # NxM matrix 

    # Flatten the arrays which means that we append all rows an array  into one "vector".
    # Note that we take the transpose of B before we flatten the array.
    cst = A.flatten()
    edge = B.T.flatten()

    # Both "vectors" now have the same number (NxM) of elements that are all possible
    # combinations of grid points along the coast and ice edge.

    # Since the coordinates are complex numbers it is straightforward to calculate the distances
    # between all grid points along the coast and along the ice edge. By reshaping the "vector"
    # back into a MxN 
    dst = np.abs(cst - edge).reshape((len(b), -1))

    # the shortest distance between the coast and the edge of the landfast sea ice
    # is the minimum along each column (axis=0).
    edge_dist = np.min(dst, axis=0)

    # plot the ice edge distance as a function of the coastline
    fig, ax = plt.subplots(1, 1, figsize=(9, 4.5))
    ax.plot(edge_dist / 1000.)
    ax.set(xlabel='Coast Grid Points',
        ylabel='Distance [km]',
        title='Distance Between Coast and Landfast Ice Edge')



    # find unique indices of grid points along the ice edge where the distance to land
    # is shortest
    edge_ind = np.unique(np.argmin(dst, axis=0))

    # select the part of the ice edge that is between the first and last of those indices
    lf_lon = lon_edge.isel(coastline=slice(edge_ind[0], edge_ind[-1]))
    lf_lat = lat_edge.isel(coastline=slice(edge_ind[0], edge_ind[-1]))

    # select only those points along the ice edge that are closest to land
    shortest_lon = lon_edge.isel(coastline=edge_ind)
    shortest_lat = lat_edge.isel(coastline=edge_ind)

    # create plot
    fig, ax = plt.subplots(1, 1, figsize=(20, 20), subplot_kw={'projection': rot})
    ax.set_extent([-61.5, -56.5, 53.8, 61.5])

    #Declare the land and ocean parameters
    LAND_highres = cartopy.feature.NaturalEarthFeature('physical', 'land', '10m',
    edgecolor='black',
    facecolor=('silver'),
    linewidth=1)
    OCEAN_highres = cartopy.feature.NaturalEarthFeature('physical', 'ocean', '10m',
    facecolor='dimgrey')

    ax.add_feature(LAND_highres,zorder=2)
    ax.add_feature(OCEAN_highres,zorder=3)
    ax.coastlines(resolution='10m',linewidth=0.35,zorder=3)

    cs = ax.contourf(lon, lat, CT,
                    vmin=0,
                    vmax=12,
                    cmap=cmo.ice,
                    transform=ccrs.PlateCarree(),
                    zorder=4)
    ax.plot(shortest_lon, shortest_lat,
            '.', color='C1',
            transform=ccrs.PlateCarree(),
            zorder=4)
    ax.plot(lf_lon, lf_lat,
            color='k',
            linewidth=2,
            transform=ccrs.PlateCarree(),
            zorder=4)

    plt.colorbar(cs, ax=ax)

    return(lf_lon, lf_lat, shortest_lon, shortest_lat)