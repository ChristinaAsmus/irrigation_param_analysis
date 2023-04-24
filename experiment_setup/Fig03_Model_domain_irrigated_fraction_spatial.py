#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
creates Fig03 
plottin gof the model domain together with its analysis regions and example grid cell used in 4.2. 

"""


import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os 

import sys
from  analysis_functions.functions_reading_files import *
from  analysis_functions.functions_correcting_time import * 
from  analysis_functions.functions_plotting import * 
from  analysis_functions.functions_calculations import *


dir_working=os.getcwd()
# creates dir in parent directory 
dir_out=os.path.join(os.pardir,'Figures') 
if not os.path.exists(dir_out):
   os.makedirs(dir_out)
print('Output directory is: ', dir_out)

# background map: irrifrac
remo_dir = '/work/ch0636/g300099/SIMULATIONS/GAR11/remo_results/067016/2017/xt/'
remo_files = 'e067016t2017060100.nc'
remo_tfile = xr.open_dataset(remo_dir+remo_files)
irrifrac=remo_tfile.IRRIFRAC[0] 

# In[]: Irrifrac plot with regions

# add regions to a background map
regions=True
background=irrifrac
x_Italy=[background.rlon[60],background.rlon[110],background.rlon[110],background.rlon[60],background.rlon[60]]
y_Italy=[background.rlat[70],background.rlat[70],background.rlat[50],background.rlat[50],background.rlat[70]]

x_Spain=[background.rlon[1],background.rlon[40],background.rlon[40],background.rlon[1],background.rlon[1]]
y_Spain=[background.rlat[70],background.rlat[70],background.rlat[30],background.rlat[30],background.rlat[70]]

x_France=[background.rlon[30],background.rlon[50],background.rlon[50],background.rlon[30],background.rlon[30]]
y_France=[background.rlat[105],background.rlat[105],background.rlat[90],background.rlat[90],background.rlat[105]]


x_gridbox=background.rlon[85]
y_gridbox=background.rlat[63]

irrifrac=irrifrac.where(irrifrac>0)



fig = plt.figure(figsize=(18, 4))

ax1 = fig.add_subplot(1, 1, 1,  projection=rotated_pole)
cax1= fig.add_axes([0.62, 0.12, 0.01, 0.76]) 
levels1 = np.arange(0,1.1,0.1).round(1)
ticks=levels1[::2]
rotplot=plot_rotvar(fig,irrifrac, ax1, cax1, '[-]', 'irrigated fraction  ','viridis_r',\
            levels1, ticks, 'neither' ,'vertical' )


if regions==True :
    ax1.plot(x_Italy, y_Italy, transform=rotated_pole, color='k',zorder=4, linewidth=0.8)
    ax1.plot(x_Spain, y_Spain, transform=rotated_pole, color='k',zorder=4, linewidth=0.8)
    ax1.plot(x_France, y_France, transform=rotated_pole, color='k',zorder=4, linewidth=0.8)
ax1.set_title(' ')
ax1.annotate('IT', (-8,-4.3))
ax1.annotate('CF', (-11.7,-0.5))
ax1.annotate('SF', (-14.7,-4.3))
plt.scatter(x=x_gridbox.rlon,y=y_gridbox.rlat,marker='x', s=40, facecolors='firebrick', alpha=1,   zorder=4)
ax1.annotate("example\n grid cell", xy=(-5.5,-5), xytext=(-3.5,-4.1), arrowprops=dict(arrowstyle="->",color='firebrick', lw=2), \
             color='firebrick', fontsize=10, zorder=4)
plt.savefig(str(dir_out)+'/Fig03.png',dpi=300, bbox_inches='tight')
plt.show()
