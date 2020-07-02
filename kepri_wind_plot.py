import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from metpy.units import units

# Open the netCDF file as a xarray Dataset
data = xr.open_dataset('CERSAT-GLO-BLENDED_WIND_L4-V6-OBS_FULL_TIME_SERIE_201810-201901_KEPRI.nc')

# Resample and parse the dataset into monthly
data = data.resample(time="1M").mean().metpy.parse_cf()

# To rename variables, supply a dictionary between old and new names to the rename method
data = data.rename({
    'eastward_wind': 'u',
    'northward_wind': 'v',
})

# Getting the cartopy coordinate reference system (CRS) of the projection of a DataArray
data_crs = data['wind_speed'].metpy.cartopy_crs

# Get multiple coordinates (for example, in just the x and y direction)
x, y = data['wind_speed'].metpy.coordinates('x', 'y')

# Or, we can just get a coordinate from the property
time = data['wind_speed'].metpy.time

# Select the data for this time
time_idx = 3
data_month = data.isel(time=time_idx)

# Create the matplotlib figure and axis
fig, ax = plt.subplots(1, 1, figsize=(12, 8), subplot_kw={'projection': data_crs})

# Add geographic features
ax.add_feature(cfeature.LAND.with_scale('10m'), facecolor=cfeature.COLORS['land'], zorder=1)
ax.add_feature(cfeature.COASTLINE.with_scale('10m'), edgecolor='black', zorder=1)

# Plot wind speed as filled contours
levels = np.arange(0,11)
c = ax.contourf(x, y, data_month['wind_speed'], cmap='jet')
fig.colorbar(c)

# Normalise the data for uniform arrow size
u = data_month['u']
v = data_month['v']
u_norm = u / np.sqrt(u ** 2.0 + v ** 2.0)
v_norm = v / np.sqrt(u ** 2.0 + v ** 2.0)

# Plot wind quiver
q = ax.quiver(x, y, u_norm, v_norm)

# Set extent
ax.set_extent([103, 107, -0.35, 2.5])

# Show gridlines
gl = ax.gridlines(draw_labels=True)
gl.xlabels_top = False
gl.ylabels_right = False

# Set a title and show the plot
ax.set_title('Wind speed (m/s) and direction at '
             + time[time_idx].dt.strftime('%Y-%m').item())
plt.show()

