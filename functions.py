# Decompose time series, calculate trends and plot

import numpy as np 
import cmocean
from netCDF4 import Dataset
import glob
import matplotlib.pyplot as plt
import datetime 
from datetime import date 
import mpl_toolkits.basemap as bm
import matplotlib.animation as animation
import numpy.ma as ma
import matplotlib.cm as cm
from matplotlib.patches import Polygon
from scipy import signal, stats
import time
from mpl_toolkits.axes_grid1 import make_axes_locatable
import warnings; warnings.simplefilter('ignore')

def trend(x, y, alpha=0.05):
    '''
    Calculates the trend of y given the linear
    independent variable x. Outputs the mean,
    trend, and alpha-level (e.g., 0.05 for 95%)
    confidence limit on the trend.
    returns mean, trend, dtrend_95
    '''
    valid = ~np.isnan(y)
    if valid.sum() <= 1:
        return np.nan, np.nan, np.nan
    else:
        X = np.array([np.ones(len(x)), x-x.mean()])
        beta = np.linalg.lstsq(X[:,valid].T, y[valid])[0]
        yhat = np.sum(beta*X.T, axis=1)
        t_stat = stats.t.isf(alpha/2, len(x[valid])-2)
        s = np.sqrt(np.sum((y[valid] - yhat[valid])**2) / (len(x[valid])-2))
        Sxx = np.sum(X[1,valid]**2) - (np.sum(X[1,valid])**2)/len(x[valid]) # np.var(X, axis=1)[1]
        return beta[0], beta[1], t_stat * s / np.sqrt(Sxx)

def load_ice(filepath, region):
    nc = Dataset(filepath + region + '.nc', 'r')
    lon = nc.variables['longitude'][600:1450,1700:]
    lat = nc.variables['latitude'][600:1450,1700:]
    land = nc.variables['land'][600:1450,1700:]
    juld = nc.variables['juld'][:]
    E_CT = nc.variables['E_CT'][:,600:1450,1700:]
    E_FA = nc.variables['E_FA'][:,600:1450,1700:]
    E_SA = nc.variables['E_SA'][:,600:1450,1700:]
    nc.close()

    # Convert concentration egg codes to actual concentrations
    CT = E_CT 
    CT[E_CT==1] = 10
    CT[E_CT==2] = 20
    CT[E_CT==3] = 30
    CT[E_CT==4] = 40
    CT[E_CT==5] = 50
    CT[E_CT==6] = 60
    CT[E_CT==7] = 70
    CT[E_CT==8] = 80
    CT[E_CT==9] = 90
    CT[E_CT==10] = 95
    CT[E_CT==11] = 100
    return lon, lat, juld, CT, E_FA, E_SA

def find_trends(filepath, region, t_domain, lon, lat, juld, array):
    '''
    Calculates montly means and trends. t_domain is a string that specifies what domain to analyze (monthly or total)
    array is the data to be analyzed (ie. concentration, thickness)
    '''
    # lon, lat, juld, CT, E_FA, E_SA = load_ice(filepath,region)

    if input == 'monthly':
        # Convert dates
        d0ord = date(1950,1,1).toordinal()
        dt_ordinal = d0ord + juld
        dates = [date.fromordinal(dt_ordinal[tt]) for tt in range(len(juld))]
        months = [dates[tt].month for tt in range(len(juld))]
        months_unique = np.unique(months)

        # Calculate monthly means and subtract them from time series
        dCT_monthly=[]
        month_occurences=[]
        for month in months_unique:
            inds = np.array(months)==month
            # Monthly means
            # ybar_monthly.append(np.nanmean(CT[inds,:,:],axis=0))
            # Subtract mean from time series
            dCT_monthly.append(CT[inds,:,:] - np.nanmean(CT[inds,:,:],axis=0))
            # Indices of dates within each month
            month_occurences.append(juld[inds])

        n = np.shape(lon)[0]
        m = np.shape(lon)[1]

        ytrend_monthly = np.zeros((12,n,m))
        dtrend_95_monthly = np.zeros((12,n,m))
        ybar_monthly = np.zeros((12,n,m))

        # Calculate trends along time axis
        for month in range(len(ytrend_monthly)):
            for i in range(n):
                for j in range(m):
                    mean, tr, dt95 = trend(month_occurences[month], dCT_monthly[month][:,i,j])
                    ytrend_monthly[month,i,j] = tr
                    dtrend_95_monthly[month,i,j] = dt95
                    ybar_monthly[month,i,j] = mean

        ##### Anomolies #####

        # # Create 3d juld array, subtract mean from the array and multiply it by dCT
        # juld3d = np.repeat(juld[:, np.newaxis], n, axis=1)
        # juld3d = np.repeat(juld3d[:,:, np.newaxis], m, axis=2)

        # ybarCT = np.nanmean(CT, axis=0)
        # # Subtract the mean from the time series
        # dCT = CT - ybarCT
        # dCT_norm = dCT*(juld3d-juld3d.mean())

        # dCT_monthly_norm=[]
        # month_occurences3d=[]
        # for month in months_unique:
        #     inds = np.array(months)==month
        #     dCT_monthly_norm.append(dCT_norm[inds])
        #     month_occurences3d.append(juld3d[inds])

        # # Multiply ytrend by normalized juld 
        # yvar_monthly = []
        # for i in range(len(ytrend_monthly2)):
        #     days_norm = month_occurences3d[i] - month_occurences3d[i].mean()
        #     ytrend_norm = ytrend_monthly[i,:,:]*days_norm
        #     yvar_monthly.append(dCT_monthly[i] - ytrend_norm[i,:,:])

        
        return ybar_monthly, ytrend_monthly, dtrend_95_monthly, lon, lat, juld
    
    if input == 'total':
        # Mean of entire time series
        ybarCT = np.nanmean(CT, axis=0)

        # Subtract the mean from the time series
        dCT = CT - ybarCT

        n = np.shape(lon)[0]
        m = np.shape(lon)[1]
        ytrend=np.zeros((n,m))
        dtrend_95=np.zeros((n,m))
        ybar=np.zeros((n,m))

        # Calculate trends along time axis
        for i in range(n):
            for j in range(m):
                mean, tr, dt95 = trend(juld, dCT[:,i,j])
                ytrend[i,j] = tr
                dtrend_95[i,j] = dt95
                ybar[i,j] = mean

        return ybar, ytrend, dtrend_95, lon, lat, juld



