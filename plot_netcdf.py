# %% 
import netCDF4 as nc
import pylab as plt
# %%
data_to_plot=nc.Dataset('u-be647_ann_mean_surface_level_o3.nc')
# %%
data_to_plot.variables['o3']
# %%

# plot mean O3 timeseries

# plot surface area weighted mean O3

# plot O3 as a colormap
plt.pcolormesh(data_to_plot.variables['o3'][0,0,:,:]*28.8/48*1e9, vmin=0, vmax=50)
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
# %%
