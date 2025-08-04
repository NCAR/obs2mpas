import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as tck
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import NaturalEarthFeature
import h5py as h5
import netCDF4

def main():

    # MPAS init file
    init_path = "/glade/campaign/mmm/parc/taosun/pandac/15kmMeshGFS/MPAS_IC/2018041500/x1.2621442.init.2018-04-15_00.00.00.nc"
    init_file = netCDF4.Dataset(init_path, 'r')

    # @ model cells
    lats = np.array( init_file.variables['latCell'][:] ) * 180.0 / np.pi
    lons = np.array( init_file.variables['lonCell'][:] ) * 180.0 / np.pi
    lons = (lons + 180) % 360 - 180
    # Check ranges
    print("Latitude range:", lats.min(), lats.max())
    print("Longitude range:", lons.min(), lons.max())

    fname = 'newfile.nc'
    nc = netCDF4.Dataset(fname, 'r')
    variables = list(nc.variables.keys())
    
    for varname in variables:
   
        data = nc.variables[varname][0, :]  # time step zero
        fill = nc.variables[varname].getncattr('_FillValue') if '_FillValue' in nc.variables[varname].ncattrs() else -999.0
        mask = (data != fill) & np.isfinite(data)
        d = data[mask]
        lat = np.ravel(lats[mask])
        lon = np.ravel(lons[mask])
        
        print('Making plot for: '+varname)
        scatter(lons,lats,d,varname,'jet',varname)  


def scatter(lon,lat,data,title,colormap,savename):
    proj = ccrs.PlateCarree(central_longitude=180)
    extent = [-180,180,-90,90]
    fig = plt.figure(figsize=(8,8))
    ax = plt.axes(projection=proj)
    background(ax, extent)
    vmin = np.nanmin(data)
    vmax = np.nanmax(data)
    cmap = colormap

    cntr = ax.scatter(lon,lat, c=data, s=0.5, vmin=vmin, vmax=vmax, cmap=cmap, transform=ccrs.PlateCarree())
    plt.title( title+', nlocs: '+str(len(data)) )
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom",size="3%", pad=0.2,axes_class=plt.Axes)
    cbar = plt.colorbar(cntr,cax=cax,orientation='horizontal',extend='both')

    outputfolder = './scatter'
    if not os.path.exists(outputfolder):
      os.makedirs(outputfolder)

    plt.savefig(outputfolder+'/'+savename+'.png',dpi=300,bbox_inches='tight')
    plt.close()


def background(ax,extent):
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='black', alpha=0.5, linestyle='dotted')
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = True
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
    ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}
    return ax


if __name__ == '__main__':
    main()
