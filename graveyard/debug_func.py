
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
# import warnings; warnings.simplefilter('ignore')

path = '/extra-space1/data/tikoralukupload/cis-weekly/nc/'
region = 'HB'

nc = Dataset(path + region + '.nc', 'r')
x = nc.variables['juld'][:]
CT = nc.variables['E_CT'][:,600:1450,1700:]



n=10
m=10

for i in range(n):
    for j in range(m):
        y = CT[:,i,j]
        valid = ~np.isnan(y)
        if valid.sum() <= 1:
            print('stop')
        else:
            X = np.array([x-x.mean(), np.ones(len(x))])
            beta = np.linalg.lstsq(X[:,valid].T, y[valid])[0]
            yhat = np.sum(beta*X.T, axis=1)
            t_stat = stats.t.isf(alpha/2, len(x[valid])-2)
            s = np.sqrt(np.sum((y[valid] - yhat[valid])**2) / (len(x[valid])-2))
            Sxx = np.sum(X[1,valid]**2) - (np.sum(X[1,valid])**2)/len(x[valid]) # np.var(X, axis=1)[1]
            # return beta[0], beta[1], t_stat * s / np.sqrt(Sxx)
            print(beta[0])