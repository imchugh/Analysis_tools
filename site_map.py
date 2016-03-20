# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 18:04:27 2016

@author: imchugh
"""

import os
import numpy as np
import linecache
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap
import gdal

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

path = '/media/Data/Dropbox/Data_sites non flux/Site data plots and non-flux/Sites/Whroo/Data/DEM/DEM_1sec'
f = 'whroo_dem.asc'

bound_box_llcrnrlon = 144.99344
bound_box_llcrnrlat = -36.70387
bound_box_urcrnrlon = 145.06777
bound_box_urcrnrlat = -36.65189
centre_lat = (bound_box_llcrnrlat + bound_box_urcrnrlat) / 2
centre_lon = (bound_box_llcrnrlon + bound_box_urcrnrlon) / 2
tower_lat = -36.673215
tower_lon = 145.029247

target = os.path.join(path, f)

hdr = [linecache.getline(target, i) for i in range(1,7)]
values = [float(h.split(" ")[-1].strip()) for h in hdr]
cols, rows, lx, ly, cell, nd = values

# DEM Raster datum is WGS84; EPSG code 4326
# Is the epsg kwarg in the Basemap class instance specifying a projection?!
map = Basemap(projection = 'merc',
              lat_0 = centre_lat, lon_0 = centre_lon,
              llcrnrlat = bound_box_llcrnrlat,
              llcrnrlon = bound_box_llcrnrlon,
              urcrnrlat = bound_box_urcrnrlat,
              urcrnrlon = bound_box_urcrnrlon)
#              epsg = 4326)

#data = np.loadtxt(target, skiprows=6)
#x = np.linspace(map.llcrnrx, map.urcrnrx, data.shape[1])
#y = np.linspace(map.llcrnry, map.urcrnry, data.shape[0])

ds = gdal.Open(target)
data = ds.ReadAsArray()
x = np.linspace(map.llcrnrx, map.urcrnrx, data.shape[1])
y = np.linspace(map.llcrnry, map.urcrnry, data.shape[0])
if map.llcrnrlat < 0:
    y = y [::-1]

xx, yy = np.meshgrid(x, y)

fig, ax = plt.subplots(1, 1, figsize = (12, 9))
fig.patch.set_facecolor('white')
#ax.set_title('Whroo tower site', y = 1.03, fontsize = 20)

cmap = plt.get_cmap('gray')
new_cmap = truncate_colormap(cmap, 0.3, 1)
color_min = int(data.min() / 10.0) * 10
color_max = int(data.max() / 10.0 + 1) * 10
colormesh = map.pcolormesh(xx, yy, data, vmin = color_min, 
                           vmax = color_max, cmap=new_cmap)
contour = map.contour(xx, yy, data, colors = 'k')
plt.clabel(contour, inline=True, fmt='%1.0f', fontsize=12, colors='k')
cb = map.colorbar(colormesh, pad = '5%')
cb.set_label('Elevation (m)', labelpad = 10)


point_x, point_y = map(tower_lon, tower_lat)
map.plot(point_x, point_y, marker = 's', color = 'black', markersize = 8)
plt.text(point_x + 250, point_y - 80, 'Tower', fontsize = 18)
         
clr_ax = cb.ax
text = clr_ax.yaxis.label
font = matplotlib.font_manager.FontProperties(size = 18)
text.set_font_properties(font)

map.drawmapscale(145.005, -36.7, tower_lon, tower_lat, 2000, units = 'm', barstyle = 'fancy')

plt.show()
