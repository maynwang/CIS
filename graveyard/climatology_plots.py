'''
    Plot maps of GLORYS climatology for use in CONOC Atlas Mapping
'''

import numpy as np
import matplotlib as mpl
mpl.interactive(True)
from matplotlib import pyplot as plt
import mpl_toolkits.basemap as bm
import cmocean

# Load in data

data = np.load('/data/home/data/GLORYS_Labrador/CMEMS_Climatology_1993_2017.npz')
lon_G = data['lon']
lat_G = data['lat']
depth = data['depth']
#temp = data['temp']
#salt = data['salt']
u = data['u']
v = data['v']
#c = data['c']
#h = data['h']

data = np.load('/home/mwang/results/CCI_Climatology_1985_2019.npz')
lon = data['lon']
lat = data['lat']
sst = data['sst_clim']
c = 100.*data['c_clim']
sst_tr = data['sst_tr']
c_tr = 100.*data['c_tr']
sst_diff = data['sst_tr']
c_diff = 100.*data['c_tr']

# Apply Lake Melville mask to CCI SST
i0 = 174
j1 = 80
j2 = 110
sst[j1:j2+1,:i0,:] = np.nan
c[j1:j2+1,:i0,:] = np.nan
sst_tr[j1:j2+1,:i0,:] = np.nan
c_tr[j1:j2+1,:i0,:] = np.nan
sst_diff[j1:j2+1,:i0,:] = np.nan
c_diff[j1:j2+1,:i0,:] = np.nan

# Plots

domain = [53.00, -63.00, 59.00, -55.50]
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='h')
llon, llat = np.meshgrid(lon, lat)
llon_G, llat_G = np.meshgrid(lon_G, lat_G)
lonproj, latproj = proj(llon, llat)
lonproj_G, latproj_G = proj(llon_G, llat_G)
lonproj_uv, latproj_uv = proj(-62., 55.)

#month = 9 # 1=Jan, 2=Feb, etc.
#m = month - 1

sstlevels = [-2,-1,0,1,2,3,4,5,6,7,8]
#clevels = [0,0.1,0.3,0.5,0.7,0.9,1.0]
#clevels = [0,10,30,50,70,90,100]
clevels = [0,0.1,25,50,75,100]
d = 3
sc = 8

plt.clf()
plt.subplot(1,4,1)
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.75)
plt.contourf(lonproj, latproj, sst[:,:,9-1], levels=sstlevels, cmap=cmocean.cm.thermal)
plt.subplot(1,4,2)
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.75)
plt.contourf(lonproj, latproj, sst[:,:,12-1], levels=sstlevels, cmap=cmocean.cm.thermal)
plt.subplot(1,4,3)
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.75)
plt.contourf(lonproj, latproj, sst[:,:,3-1], levels=sstlevels, cmap=cmocean.cm.thermal)
plt.subplot(1,4,4)
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.75)
plt.contourf(lonproj, latproj, sst[:,:,6-1], levels=sstlevels, cmap=cmocean.cm.thermal)


plt.clf()
plt.subplot(1,4,1)
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.75)
plt.contourf(lonproj, latproj, c[:,:,9-1], levels=clevels, cmap=cmocean.cm.ice)
plt.subplot(1,4,2)
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.75)
plt.contourf(lonproj, latproj, c[:,:,12-1], levels=clevels, cmap=cmocean.cm.ice)
plt.subplot(1,4,3)
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.75)
plt.contourf(lonproj, latproj, c[:,:,3-1], levels=clevels, cmap=cmocean.cm.ice)
plt.subplot(1,4,4)
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.75)
plt.contourf(lonproj, latproj, c[:,:,6-1], levels=clevels, cmap=cmocean.cm.ice)



plt.clf()
plt.subplot(1,4,1)
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.75)
plt.quiver(lonproj_G[::d,::d], latproj_G[::d,::d], u[0,::d,::d,9-1], v[0,::d,::d,9-1])
plt.subplot(1,4,2)
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.75)
plt.quiver(lonproj_G[::d,::d], latproj_G[::d,::d], u[0,::d,::d,12-1], v[0,::d,::d,12-1])
plt.subplot(1,4,3)
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.75)
plt.quiver(lonproj_G[::d,::d], latproj_G[::d,::d], u[0,::d,::d,3-1], v[0,::d,::d,3-1])
plt.subplot(1,4,4)
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.75)
plt.quiver(lonproj_G[::d,::d], latproj_G[::d,::d], u[0,::d,::d,6-1], v[0,::d,::d,6-1])

plt.show()


plt.streamplot(lonproj_G, latproj_G, u[0,:,:,9-1], v[0,:,:,9-1]) #, linewidth=s[0,:,:,9-1], density=10, color='k')
plt.subplot(1,4,2)
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.75)
plt.contourf(lonproj, latproj, u[0,:,:,12-1], levels=clevels, cmap=cmocean.cm.ice)
plt.subplot(1,4,3)
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.75)
plt.contourf(lonproj, latproj, u[0,:,:,3-1], levels=clevels, cmap=cmocean.cm.ice)
plt.subplot(1,4,4)
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.75)
plt.contourf(lonproj, latproj, u[0,:,:,6-1], levels=clevels, cmap=cmocean.cm.ice)

H = plt.colorbar()
H.set_label('SST')
#H.set_label('Depth-averaged currents from GLORYS (m s$^{-1}$)')

#
# Figures for CONOC Atlas
#

# Mean states

plt.figure(1)
plt.clf()
plt.subplot(2,3,1, facecolor='0.85'); m=9-1
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.5)
CS = plt.contourf(lonproj, latproj, sst[:,:,m], levels=sstlevels, cmap=cmocean.cm.thermal)
plt.clim(-4,8)
CI = plt.contourf(lonproj, latproj, c[:,:,m], levels=clevels[2:], cmap=cmocean.cm.ice)
plt.clim(0,100)
plt.contour(lonproj, latproj, c[:,:,m], levels=[clevels[2]], colors='k', linewidths=2, linestyles='solid')
plt.quiver(lonproj_G[::d,::d], latproj_G[::d,::d], u[0,::d,::d,m], v[0,::d,::d,m], scale=sc, width=0.005, pivot='middle', color='0.15')
plt.title('September')

plt.subplot(2,3,2, facecolor='0.85'); m=12-1
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.5)
plt.contourf(lonproj, latproj, sst[:,:,m], levels=sstlevels, cmap=cmocean.cm.thermal)
plt.clim(-4,8)
plt.contourf(lonproj, latproj, c[:,:,m], levels=clevels[2:], cmap=cmocean.cm.ice)
plt.clim(0,100)
plt.contour(lonproj, latproj, c[:,:,m], levels=[clevels[2]], colors='k', linewidths=2, linestyles='solid')
plt.quiver(lonproj_G[::d,::d], latproj_G[::d,::d], u[0,::d,::d,m], v[0,::d,::d,m], scale=sc, width=0.005, pivot='middle', color='0.15')
plt.title('December')

plt.subplot(2,3,4, facecolor='0.85'); m=3-1
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.5)
plt.contourf(lonproj, latproj, sst[:,:,m], levels=sstlevels, cmap=cmocean.cm.thermal)
plt.clim(-4,8)
plt.contourf(lonproj, latproj, c[:,:,m], levels=clevels[2:], cmap=cmocean.cm.ice)
plt.clim(0,100)
plt.contour(lonproj, latproj, c[:,:,m], levels=[clevels[2]], colors='k', linewidths=2, linestyles='solid')
plt.quiver(lonproj_G[::d,::d], latproj_G[::d,::d], u[0,::d,::d,m], v[0,::d,::d,m], scale=sc, width=0.005, pivot='middle', color='0.15')
plt.title('March')

plt.subplot(2,3,5, facecolor='0.85'); m=6-1
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.5)
plt.contourf(lonproj, latproj, sst[:,:,m], levels=sstlevels, cmap=cmocean.cm.thermal)
plt.clim(-4,8)
plt.contourf(lonproj, latproj, c[:,:,m], levels=clevels[2:], cmap=cmocean.cm.ice)
plt.clim(0,100)
plt.contour(lonproj, latproj, c[:,:,m], levels=[clevels[2]], colors='k', linewidths=2, linestyles='solid')
plt.quiver(lonproj_G[::d,::d], latproj_G[::d,::d], u[0,::d,::d,m], v[0,::d,::d,m], scale=sc, width=0.005, pivot='middle', color='0.15')
plt.title('June')

plt.subplot(3,6,5)
HS = plt.colorbar(CS)
HS.set_label('Sea surface temperature ($^\circ$C)')
plt.subplot(3,6,6); m=6-1
HI = plt.colorbar(CI)
HI.set_label('Sea ice concentration (% area)')

# plt.savefig('Labrador_SST_ice_uv_mean_orig.pdf', bbox_inches='tight')

# Trends

ssttrlevels = [-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5]
ctrlevels = [-25,-20,-15,-10,-5,0,5]
sst_tr_tmp = sst_tr.copy()
sst_tr_tmp[sst_tr_tmp>0.5] = 0.5
sst_tr_tmp[sst_tr_tmp<-0.5] = -0.5
c_tr_tmp = c_tr.copy()
c_tr_tmp[c<clevels[2]] = np.nan;
c_tr_tmp[c_tr_tmp<-25] = -25

plt.figure(2)
plt.clf()
plt.subplot(2,3,1, facecolor='0.85'); m=9-1
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.5)
CS = plt.contourf(lonproj, latproj, sst_tr_tmp[:,:,m], levels=ssttrlevels, cmap=cmocean.cm.balance)
plt.clim(-0.7, 0.7)
CI = plt.contourf(lonproj, latproj, c_tr_tmp[:,:,m], levels=ctrlevels, cmap=cmocean.cm.diff)
plt.clim(-30, 30)
plt.contour(lonproj, latproj, c[:,:,m], levels=[clevels[2]], colors='k', linewidths=2, linestyles='solid')
plt.title('September')

plt.subplot(2,3,2, facecolor='0.85'); m=12-1
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.5)
CS = plt.contourf(lonproj, latproj, sst_tr_tmp[:,:,m], levels=ssttrlevels, cmap=cmocean.cm.balance)
plt.clim(-0.7, 0.7)
CI = plt.contourf(lonproj, latproj, c_tr_tmp[:,:,m], levels=ctrlevels, cmap=cmocean.cm.diff)
plt.clim(-30, 30)
plt.contour(lonproj, latproj, c[:,:,m], levels=[clevels[2]], colors='k', linewidths=2, linestyles='solid')
plt.title('December')

plt.subplot(2,3,4, facecolor='0.85'); m=3-1
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.5)
CS = plt.contourf(lonproj, latproj, sst_tr_tmp[:,:,m], levels=ssttrlevels, cmap=cmocean.cm.balance)
plt.clim(-0.7, 0.7)
CI = plt.contourf(lonproj, latproj, c_tr_tmp[:,:,m], levels=ctrlevels, cmap=cmocean.cm.diff)
plt.clim(-30, 30)
plt.contour(lonproj, latproj, c[:,:,m], levels=[clevels[2]], colors='k', linewidths=2, linestyles='solid')
plt.title('March')

plt.subplot(2,3,5, facecolor='0.85'); m=6-1
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.5)
CS = plt.contourf(lonproj, latproj, sst_tr_tmp[:,:,m], levels=ssttrlevels, cmap=cmocean.cm.balance)
plt.clim(-0.7, 0.7)
CI = plt.contourf(lonproj, latproj, c_tr_tmp[:,:,m], levels=ctrlevels, cmap=cmocean.cm.diff)
plt.clim(-30, 30)
plt.contour(lonproj, latproj, c[:,:,m], levels=[clevels[2]], colors='k', linewidths=2, linestyles='solid')
plt.title('June')

plt.subplot(3,6,5)
HS = plt.colorbar(CS)
HS.set_label('Sea surface temperature trend ($^\circ$C per decade)')
plt.subplot(3,6,6); m=6-1
HI = plt.colorbar(CI)
HI.set_label('Sea ice concentration trend (% per decade)')

# plt.savefig('Labrador_SST_ice_uv_trends_orig.pdf', bbox_inches='tight')

# Change from 1985-1999 to 2005-2019

sstdflevels = [0, 0.2, 0.4, 0.6, 0.8, 1]
cdflevels = [-30,-25,-20,-15,-10,-5,0]
sst_diff_tmp = sst_diff.copy()
#sst_diff_tmp[sst_diff_tmp>0.5] = 0.5
sst_diff_tmp[sst_diff_tmp<0.] = 0.
c_diff_tmp = c_diff.copy()
c_diff_tmp[c<clevels[2]] = np.nan;
c_diff_tmp[c_diff_tmp>0.] = 0.

plt.figure(3)
plt.clf()
plt.subplot(2,3,1, facecolor='0.85'); m=9-1
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.5)
CS = plt.contourf(lonproj, latproj, sst_diff_tmp[:,:,m], levels=sstdflevels, cmap=cmocean.cm.balance)
plt.clim(-1, 1)
CI = plt.contourf(lonproj, latproj, c_diff_tmp[:,:,m], levels=cdflevels, cmap=cmocean.cm.diff)
plt.clim(-30, 30)
plt.contour(lonproj, latproj, c[:,:,m], levels=[clevels[2]], colors='k', linewidths=2, linestyles='solid')
plt.title('September')

plt.subplot(2,3,2, facecolor='0.85'); m=12-1
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.5)
CS = plt.contourf(lonproj, latproj, sst_diff_tmp[:,:,m], levels=sstdflevels, cmap=cmocean.cm.balance)
plt.clim(-1, 1)
CI = plt.contourf(lonproj, latproj, c_diff_tmp[:,:,m], levels=cdflevels, cmap=cmocean.cm.diff)
plt.clim(-30, 30)
plt.contour(lonproj, latproj, c[:,:,m], levels=[clevels[2]], colors='k', linewidths=2, linestyles='solid')
plt.title('December')

plt.subplot(2,3,4, facecolor='0.85'); m=3-1
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.5)
CS = plt.contourf(lonproj, latproj, sst_diff_tmp[:,:,m], levels=sstdflevels, cmap=cmocean.cm.balance)
plt.clim(-1, 1)
CI = plt.contourf(lonproj, latproj, c_diff_tmp[:,:,m], levels=cdflevels, cmap=cmocean.cm.diff)
plt.clim(-30, 30)
plt.contour(lonproj, latproj, c[:,:,m], levels=[clevels[2]], colors='k', linewidths=2, linestyles='solid')
plt.title('March')

plt.subplot(2,3,5, facecolor='0.85'); m=6-1
proj.fillcontinents(color='0.85', lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines(color='k', linewidth=0.5)
CS = plt.contourf(lonproj, latproj, sst_diff_tmp[:,:,m], levels=sstdflevels, cmap=cmocean.cm.balance)
plt.clim(-1, 1)
CI = plt.contourf(lonproj, latproj, c_diff_tmp[:,:,m], levels=cdflevels, cmap=cmocean.cm.diff)
plt.clim(-30, 30)
plt.contour(lonproj, latproj, c[:,:,m], levels=[clevels[2]], colors='k', linewidths=2, linestyles='solid')
plt.title('June')

plt.subplot(3,6,5)
HS = plt.colorbar(CS)
HS.set_label('Sea surface temperature change ($^\circ$C)')
plt.subplot(3,6,6); m=6-1
HI = plt.colorbar(CI)
HI.set_label('Sea ice concentration change (% area)')

# plt.savefig('Labrador_SST_ice_uv_change_orig.pdf', bbox_inches='tight')

