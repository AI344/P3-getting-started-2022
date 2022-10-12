# %% 
import netCDF4 as nc
import pylab as plt
import numpy as np
# %%
data_to_plot=nc.Dataset('u-be647_ann_mean_surface_level_o3.nc')
# %%
data_to_plot.variables['o3']
# %%
# %% naive mean without surface area weighting
meano3 = np.mean(data_to_plot.variables['o3'][:]*28.8/48*1e9*surface(lat,lon)/sum(surface), axis=(0,2,3))

# %%
# test code to calculate surface area weighting - AI to check and fix
for ilat in lats:
    for ilon in lons:
        wtdo3 = o3[lat,lon]*surfacearea(lat,lon)
        
areawtmeano3 = np.sum(wtdo3)/np.sum(surfacearea)

###---### using xarray module

# %%
# plot O3 as a colormap
plt.pcolormesh(data_to_plot.variables['o3'][6,0,:,:]*28.8/48*1e9, vmin=0, vmax=50)
plt.colorbar()
# %%
import cartopy.crs as ccrs

# %%
ax = plt.axes(projection=ccrs.PlateCarree())
plt.pcolormesh(data_to_plot.variables['longitude'][:],
               data_to_plot.variables['latitude'][:],
               data_to_plot.variables['o3'][0,0,:,:]*28.8/48*1e9, 
               transform=ccrs.PlateCarree(),
               vmin=0,
               vmax=50.
               )
#ax.set_global()
ax.coastlines()
#ax.stock_img()
plt.show()
#plt.colorbar()
# %%
