#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 11:28:27 2023

@author: g300099
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

# In[]:experiment for analyzing 

year=2017
month=6

exp_number_noirri='067015'
exp_number_irri_prescribed='067020'
exp_number_irri_flextime='067019'
exp_number_irri_adapt='067017'

# In[]: define plot directory

dir_out='/work/ch0636/g300099/SIMULATIONS/GAR11/plot/plot_for_paper1/'


# In[]: read data 

# background map: irrifrac
remo_dir = '/work/ch0636/g300099/SIMULATIONS/GAR11/remo_results/067016/2017/xt/'
remo_files = 'e067016t2017060100.nc'
remo_tfile = xr.open_dataset(remo_dir+remo_files)
irrifrac=remo_tfile.IRRIFRAC[0]


varlist=['IRRWR']
var_num_list=['796']
for var, var_num in zip(varlist, var_num_list): 
    single_var_data_adapt=read_efiles(var, var_num, exp_number_irri_adapt, year, month)
    single_var_data_prescribed=read_efiles(var, var_num, exp_number_irri_prescribed, year, month)
    if var==varlist[0]: 
        ds_var_irri_adapt=single_var_data_adapt
        ds_var_irri_prescribed=single_var_data_prescribed
    else:
        ds_var_irri_adapt=xr.merge([ds_var_irri_adapt, single_var_data_adapt])
        ds_var_irri_prescribed=xr.merge([ds_var_irri_prescribed, single_var_data_prescribed])
ds_var_irri_adapt=xr.merge([ds_var_irri_adapt, irrifrac])
ds_var_irri_prescribed=xr.merge([ds_var_irri_prescribed, irrifrac])


varlist_flextime=['IRRWR', 'IRRDUR']
var_num_list=['796', '795']
for var, var_num in zip(varlist_flextime, var_num_list): 
    single_var_data=read_efiles(var, var_num, exp_number_irri_flextime, year, month)
    if var==varlist[0]: 
        ds_var_irri=single_var_data
    else:
        ds_var_irri=xr.merge([ds_var_irri, single_var_data])
ds_var_irri_flextime=xr.merge([ds_var_irri, irrifrac])


dsirr_adapt = correct_timedim(ds_var_irri_adapt)
dsirr_flextime=correct_timedim(ds_var_irri_flextime)
dsirr_prescribed=correct_timedim(ds_var_irri_prescribed)

   
# irrwr
irrwr_prescribed=dsirr_prescribed.IRRWR.groupby('time.hour')[18].sum('time')*1000
irrwr_prescribed=irrwr_prescribed.where(irrwr_prescribed>0)

irrwr_flextime=dsirr_flextime.IRRWR.groupby('time.hour')[23].sum('time')*1000
irrwr_flextime=irrwr_flextime.where(irrwr_flextime>0)
irrwr_adapt=dsirr_adapt.IRRWR.where(dsirr_adapt.IRRWR != 0).resample(time='1D').max().sum('time')*1000
irrwr_adapt=irrwr_adapt.where(irrwr_adapt>0)


# spatial plot 

fig = plt.figure(figsize=(18, 4))
ax1 = fig.add_subplot(1, 3, 1, projection=rotated_pole)
cax1 = fig.add_axes([0.148, 0.01, 0.19, 0.04])
levels1=[0,50, 100, 150, 200, 250,300, 350, 400]
ticks1=levels1
#irrwr prescribed
plot_rotvar(fig, irrwr_prescribed, ax1, cax1, '[mm]', 'irrigation water [mm] ','viridis_r',\
            levels1, ticks1,'max', 'horizontal' )

ax2 = fig.add_subplot(1, 3, 2, projection=rotated_pole)
cax2 = fig.add_axes([0.42, 0.01, 0.19, 0.04])
#irrwr flextime 
plot_rotvar(fig, irrwr_flextime, ax2, cax2, '[mm]', 'irrigation water [mm] ','viridis_r',\
            levels1, ticks1,'max', 'horizontal' )
    
ax3 = fig.add_subplot(1, 3, 3, projection=rotated_pole)
cax3 = fig.add_axes([0.69, 0.01, 0.19, 0.04])
#irrwr adapt 
plot_rotvar(fig, irrwr_adapt, ax3, cax3, '[mm]', 'irrigation water [mm] ','viridis_r',\
            levels1, ticks1,'max', 'horizontal' )
#plt.show()
#plt.savefig(str(dir_out)+'/schemes_irrwr_'+str(month)+'spatial_new.png',dpi=300, bbox_inches='tight')
