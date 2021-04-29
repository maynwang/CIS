import os
import gdal
import rasterio as rio
import numpy as np
import matplotlib.pyplot as plt
import numpy as np

path = '/home/mwang/data/cis-weekly/tif/n-ct/'
dataset = gdal.Open(path + 'WA_2019-01-07_n-ct.tif')
# dataset = gdal.Open('GeoTiff_Image.tif', gdal.GA_ReadOnly) 
# Note GetRasterBand() takes band no. starting from 1 not 0
band = dataset.GetRasterBand(1)
arr = band.ReadAsArray()
plt.imshow(arr)
plt.show()


