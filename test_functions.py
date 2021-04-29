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

from functions1 import load_ice
from functions1 import find_trends

path = '/extra-space1/data/tikoralukupload/cis-weekly/nc/'
region = 'HB'

ybar, ytrend, dtrend_95 = find_trends(path, region, 'total')