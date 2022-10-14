# %% 
import netCDF4 as nc
import pylab as plt
import numpy as np
from datetime import datetime
# %%
data_to_plot=nc.Dataset('u-be647_ann_mean_surface_level_o3.nc')
# %%
data_to_plot.variables['o3']

# %% naive mean without surface area weighting and with saw 
def area_grid(lat,lon):
    '''
    Function to calculate surface area per gridbox
    Units: m2
    S = R^2*(lon2-lon1)*(sin lat2 - sin lat1)
    lon in radians, R = 6371 km
    '''
    import numpy as np
    Pi           = np.float64(3.141592653589793238462643383279)
    Earth_Radius = np.float64(6371.0*1.0E3)#equator radius:6378.1*1E3
    lat_bound    = np.float64(89.999999999999999999999999)
    lon          = np.float64(lon)
    lat          = np.float64(lat)
    rlon         = (lon[:]/np.float64(180.0))*Pi
    rlat         = (lat[:]/np.float64(180.0))*Pi
    dlat         = (rlat[1] - rlat[0])/2.0
    dlon         = (rlon[1] - rlon[0])/2.0
    #
    area = np.zeros((len(rlat),len(rlon)),np.float64)
    j=0
    while j < len(rlat):
        if (lat[j] >= lat_bound):
            lat1 = rlat[j]
            lat2 = rlat[j] - dlat/2.0
        elif (lat[j] <= -1.0*lat_bound):
            lat1 = rlat[j] + dlat/2.0
            lat2 = rlat[j]
        else:
            lat1 = rlat[j] + dlat
            lat2 = rlat[j] - dlat
        i=0
        while i < len(rlon):
            lon1 = rlon[i] - dlon
            lon2 = rlon[i] + dlon
            area[j,i] = (Earth_Radius**2)*(abs(np.sin(lat1)-np.sin(lat2))*abs(lon1-lon2))
            i += 1
        j += 1
    return area

grid_areas = area_grid(data_to_plot.variables['latitude'][:], data_to_plot.variables['longitude'][:])
tot_surface = np.sum(grid_areas)

o3_mean_list = []
o3_mean_adj = []
for timestamp in data_to_plot.variables['o3'][:,:,:,:]:
    o3_mean_list.append(np.mean(timestamp)*28.8/48*1e9)
    adj_o3 = np.multiply(timestamp,grid_areas)
    o3_mean_adj.append(np.sum(adj_o3)*28.8*1e9/(48*tot_surface))

raw_time_nums = data_to_plot.variables['time'][:]
dates = nc.num2date(raw_time_nums, data_to_plot.variables['time'].units, data_to_plot.variables['time'].calendar)
list_of_dates = []
for i in dates:
    j = datetime.strptime(str(i)[:-9], '%Y-%m-%d')
    list_of_dates.append(j)

plt.plot(list_of_dates, o3_mean_list, label='Unweighted')
plt.plot(list_of_dates, o3_mean_adj, label='Surface weighted')
plt.legend()
plt.title('Surface weighted vs unweighted global O3 concs')
plt.xlabel('Date')
plt.ylabel('Mean Global O3 Conc (ppb)')
plt.show()

# %%
# rewritten using xarray module
import xarray as xr
import cftime

ds=xr.open_dataset('u-be647_ann_mean_surface_level_o3.nc')
# times = ds.time.values
# o3_data = ds.o3.values
# lats = ds.latitude.values
# longs = ds.longitude.values

grid_areas = area_grid(ds.latitude.values, ds.longitude.values)
tot_surface = np.sum(grid_areas)

o3_mean_list = []
o3_mean_adj = []
for timestamp in ds.o3.values[:,:,:,:]:
    o3_mean_list.append(np.mean(timestamp)*28.8/48*1e9)
    adj_o3 = np.multiply(timestamp,grid_areas)
    o3_mean_adj.append(np.sum(adj_o3)*28.8*1e9/(48*tot_surface))

times = ds.time.values
conv_time = xr.times.convert_calendar('standard')


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
