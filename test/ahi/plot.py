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

def main():

    fname = 'mpas_iodav1.nc'
    ncf = h5.File(fname, 'r')
    lat = ncf['latitude@MetaData'][:]
    lon = ncf['longitude@MetaData'][:]
    bt = ncf['brightness_temperature_13@ObsValue'][:]
    cld = ncf['cloudAmount@MetaData'][:]

    scatter(lon,lat,bt,'Superobbed BT ch13','jet','btch13_so_ahi')    
    scatter(lon,lat,cld,'Superobbed CloudAmount','Greys','cldamount_so_ahi')

def scatter(lon,lat,data,title,colormap,savename):
    proj = ccrs.PlateCarree(central_longitude=180)
    lon = (lon + 180) % 360 - 180
    #extent = [-140,-10,-65,65]
    extent = [-180,180,-90,90]
    fig = plt.figure(figsize=(8,8))
    ax = plt.axes(projection=proj)
    background(ax, extent)
    if colormap == 'Greys':
       vmin = 0
       vmax = 1
    else:
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
