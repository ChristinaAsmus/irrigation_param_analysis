#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
creates Fig A2:
plotting of soil moisture (wsechirr, wsmx), drainage (DRAINIRR) and runoff (RUNOFFIR) 
for the irrigated fraction as time series for test month June for all water application schemes

"""

import xarray as xr
import numpy as np
import pandas as pd
from datetime import datetime
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import os

import sys
sys.path.append('/home/g/g300099/pyprograms/Paper1_for_publishing_final/plot_figures/functions/') 
from functions_reading_files import *
from functions_correcting_time import * 
from functions_plotting import * 

# In[]: appendix for schemes 

# read data 

year=2017
month=6

exp_number_noirri='067015'
exp_number_irri_prescribed='067020'
exp_number_irri_flextime='067019'
exp_number_irri_adapt='067017'

# background map: irrifrac
remo_dir = '/work/ch0636/g300099/SIMULATIONS/GAR11/remo_results/067016/2017/xt/'
remo_files = 'e067016t2017060100.nc'
remo_tfile = xr.open_dataset(remo_dir+remo_files)
irrifrac=remo_tfile.IRRIFRAC[0]


varlist=[ 'WSECHIRR','WSMXIRR',  'RUNOFFIR', 'DRAINIRR']
var_num_list=['701', '728',  '713', '700']
for var, var_num in zip(varlist, var_num_list): 
    single_var_data_adapt=read_efiles(var, var_num, exp_number_irri_adapt, year, month)
    single_var_data_prescribed=read_efiles(var, var_num, exp_number_irri_prescribed, year, month)
    single_var_data_flextime=read_efiles(var, var_num, exp_number_irri_flextime, year, month)
    single_var_data_noirri=read_efiles(var, var_num, exp_number_noirri, year, month)
    if var==varlist[0]: 
        ds_var_irri_adapt=single_var_data_adapt
        ds_var_irri_prescribed=single_var_data_prescribed
        ds_var_irri_flextime=single_var_data_flextime
        ds_var_noirri=single_var_data_noirri
    else:
        ds_var_irri_adapt=xr.merge([ds_var_irri_adapt, single_var_data_adapt])
        ds_var_irri_prescribed=xr.merge([ds_var_irri_prescribed, single_var_data_prescribed])
        ds_var_irri_flextime=xr.merge([ds_var_irri_flextime, single_var_data_flextime])
        ds_var_noirri=xr.merge([ds_var_noirri, single_var_data_noirri])
ds_var_noirri=xr.merge([ds_var_noirri, irrifrac])
ds_var_irri_flextime=xr.merge([ds_var_irri_flextime, irrifrac])
ds_var_irri_adapt=xr.merge([ds_var_irri_adapt, irrifrac])
ds_var_irri_prescribed=xr.merge([ds_var_irri_prescribed, irrifrac])

dsirr_adapt = correct_timedim(ds_var_irri_adapt)
dsirr_flextime=correct_timedim(ds_var_irri_flextime)
dsirr_prescribed=correct_timedim(ds_var_irri_prescribed)
dsnoirr=correct_timedim(ds_var_noirri)

# In[]: define plot directory

dir_working=os.getcwd()
# creates dir in parent directory 
dir_out=os.path.join(os.pardir,'Figures') 
if not os.path.exists(dir_out):
   os.makedirs(dir_out)
print('Output directory is: ', dir_out)

# In[]: plot DRAINIRR & RUNOFFIR

# gridcell to analyse 
x=63
y=82 

drain_irri_prescribed=dsirr_prescribed['DRAINIRR'][:,x,y].to_series()
drain_irri_flextime=dsirr_flextime['DRAINIRR'][:,x,y].to_series()
drain_irri_adapt=dsirr_adapt['DRAINIRR'][:,x,y].to_series()
drain_noirri=dsnoirr['DRAINIRR'][:,x,y].to_series()

runoff_irri_prescribed=dsirr_prescribed['RUNOFFIR'][:,x,y].to_series()
runoff_irri_flextime=dsirr_flextime['RUNOFFIR'][:,x,y].to_series()
runoff_irri_adapt=dsirr_adapt['RUNOFFIR'][:,x,y].to_series()
runoff_noirri=dsnoirr['RUNOFFIR'][:,x,y].to_series()

ws_irri_prescribed=(dsirr_prescribed['WSECHIRR'][:,x,y]/dsirr_prescribed['WSMXIRR'][:,x,y]).to_series()
ws_irri_flextime=(dsirr_flextime['WSECHIRR'][:,x,y]/dsirr_flextime['WSMXIRR'][:,x,y]).to_series()
ws_irri_adapt=(dsirr_adapt['WSECHIRR'][:,x,y]/dsirr_adapt['WSMXIRR'][:,x,y]).to_series()
ws_noirri=(dsnoirr['WSECHIRR'][:,x,y]/dsnoirr['WSMXIRR'][:,x,y]).to_series()



# figure for appendix 
legend_labels=['prescribed','flextime','adaptive']

fig = plt.figure(figsize=(18,10))  
ax1 = fig.add_subplot(311)
ws_irri_prescribed.plot.line(ax=ax1)
ws_irri_flextime.plot.line(ax=ax1)
ws_irri_adapt.plot.line(ax=ax1)
#ax1.legend(legend_labels)
ax1.set_ylabel('fraction of wsmx [-]')
ax1.grid()

ax2 = fig.add_subplot(312, sharex = ax1)
runoff_irri_prescribed.plot.line(ax=ax2)
runoff_irri_flextime.plot.line(ax=ax2)
runoff_irri_adapt.plot.line(ax=ax2)
#ax2.legend(legend_labels)
ax2.set_ylim(0,1)
ax2.set_ylabel('runoff [m]')
ax2.grid()

ax3 = fig.add_subplot(313, sharex = ax1)
drain_irri_prescribed.plot.line(ax=ax3)
drain_irri_flextime.plot.line(ax=ax3)
drain_irri_adapt.plot.line(ax=ax3)
ax3.legend(legend_labels, loc='upper right')
ax3.set_ylabel('drainage [m]')
ax3.set_xlabel(' ')
ax3.grid()
plt.savefig(str(dir_out)+'/FigA2.png',dpi=300, bbox_inches='tight')
plt.show()
