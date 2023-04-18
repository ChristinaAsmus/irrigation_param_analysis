#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 12:35:52 2023

@author: g300099
"""

import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as mtick 

import sys
sys.path.append('/home/g/g300099/pyprograms/Paper1_for_publishing_final/plot_figures/functions/') 
from functions_reading_files import *
from functions_correcting_time import * 
from functions_plotting import * 
from functions_calculations import *


# In[]: experiment number and year
  
year=2017
# for all months month=0
month=0

exp_number_irri='067016'
exp_number_noirri='067015'


# In[]: read efiles 

# background map: irrifrac
remo_dir = '/work/ch0636/g300099/SIMULATIONS/GAR11/remo_results/067016/2017/xt/'
remo_files = 'e067016t2017060100.nc'
remo_tfile = xr.open_dataset(remo_dir+remo_files)
irrifrac=remo_tfile.IRRIFRAC[0]
landfrac=remo_tfile.BLA[0]

# analysis variables 
varlist=['WSECHIRR','APRL','APRC','TEMP2','T2MAX', 'T2MIN']
var_num_list=['701','142','143','167','201','202']

for var, var_num in zip(varlist, var_num_list): 
    single_var_data=read_efiles(var, var_num, exp_number_noirri, year, month)
    if var==varlist[0]: 
        ds_var_noirri=single_var_data
    else:
        ds_var_noirri=xr.merge([ds_var_noirri, single_var_data])
ds_var_noirri=xr.merge([ds_var_noirri, irrifrac, landfrac])
dsnoirr_newtime = correct_timedim(ds_var_noirri)

# In[]: define plot directory

dir_out='/work/ch0636/g300099/SIMULATIONS/GAR11/plot/plot_for_paper1/'
    
# In[]: Define regions 

rlat_list=[(50,70),(30,70),(90,105)]
rlon_list=[(60,110),(1,40),(30,50)]
title_list=['IT','SF','CF']


# In[]: soil moisture precipitation series plot for regions

# select months of interest AMJJA
dsnoirr_newtime=dsnoirr_newtime.sel(time=dsnoirr_newtime.time.dt.month.isin([4, 5, 6, 7, 8]))
ws_var=dsnoirr_newtime['WSECHIRR']


fig = plt.figure(figsize=(18, 4))

for i in range(len(rlat_list)):
     ws_series=((ws_var.isel(rlat=slice(rlat_list[i][0],rlat_list[i][1]),rlon=slice(rlon_list[i][0],rlon_list[i][1])))\
     .resample(time='D').mean(skipna=True).mean(dim=['rlon','rlat'], skipna=True)).to_series()
     precip_series=calculate_sum_precip('APRL','APRC',dsnoirr_newtime, rlat_list, rlon_list, i).to_series()
     
     df_combined=pd.DataFrame() 
     df_combined['ws']=ws_series
     df_combined['precip']=precip_series
     df_combined.index=df_combined.index.strftime("%d.%m.")

         
     ax1 = fig.add_subplot(1, 3, i+1)
     df_combined.precip.plot.bar(ax=ax1, secondary_y=True, color='royalblue',  width=0.6, zorder=-1)
     df_combined.ws.plot.line(ax=ax1, color='k', fontsize=16, use_index=False, zorder=10)
     ax1.set_ylim(0,0.4)
     ax1.set_ylabel('soil moisture [m]', fontsize=16)
     ax1.right_ax.set_ylabel('precipitation [mm]', color='royalblue', fontsize=16)
     ax1.right_ax.set_ylim(0, 20)
     ax1.right_ax.set_yticks(np.linspace(0,20,9))
     ax1.right_ax.set_yticklabels(labels=np.linspace(0,20,9), color='royalblue', fontsize=16)
     tick_loc=mtick.FixedLocator([0,14,30,44,61,75,91,106,121,137,152])
     ax1.xaxis.set_major_locator(tick_loc)
     for label in ax1.right_ax.get_yticklabels()[1::2]:
         label.set_visible(False)
     for label in ax1.get_yticklabels()[1::2]:
         label.set_visible(False)
     for tick in ax1.get_xticklabels():
         tick.set_rotation(45)     
     ax1.set_title(title_list[i], fontsize=16)
     ax1.yaxis.grid() # horizontal lines
wspace = 0.7   
fig.subplots_adjust(wspace=wspace)
#plt.show()
#plt.savefig(str(dir_out)+'/'+str(exp_number_irri)+'_abs_soil_moisture_precip_AMJJA_new_new.png',dpi=300, bbox_inches='tight')


# In[]: temperature series plot for regions  

# select months AMJ 
dsnoirr_newtime=dsnoirr_newtime.sel(time=dsnoirr_newtime.time.dt.month.isin([4, 5, 6, 7, 8]))
title_list=[' ', '','']

fig = plt.figure(figsize=(18, 4))
for i in range(len(rlat_list)):
     temp2=(dsnoirr_newtime.isel(rlat=slice(rlat_list[i][0],rlat_list[i][1]),rlon=slice(rlon_list[i][0],rlon_list[i][1]))\
     ['TEMP2'][:,0,:,:]).resample(time='D').mean(skipna=True).mean(dim=['rlon','rlat'], skipna=True)-273.15
     t2max=(dsnoirr_newtime.isel(rlat=slice(rlat_list[i][0],rlat_list[i][1]),rlon=slice(rlon_list[i][0],rlon_list[i][1]))\
     ['T2MAX'][:,0,:,:]).resample(time='D').max(skipna=True).mean(dim=['rlon','rlat'], skipna=True)-273.15
     t2min=(dsnoirr_newtime.isel(rlat=slice(rlat_list[i][0],rlat_list[i][1]),rlon=slice(rlon_list[i][0],rlon_list[i][1]))\
     ['T2MIN'][:,0,:,:]).resample(time='D').min(skipna=True).mean(dim=['rlon','rlat'], skipna=True)-273.15    
         
     ax1 = fig.add_subplot(1, 3, i+1)
     temp2.plot.line(ax=ax1, color='blue', label='2m temperature')
     t2max.plot.line(ax=ax1, color='orange', label='2m max. temperature')
     t2min.plot.line(ax=ax1, color='green', label='2m min. temperature')
     ax1.set_ylim(0,40)
     ax1.set_ylabel('temperature [°C]', fontsize=16)
     ax1.xaxis.set_major_formatter(mdates.DateFormatter('%d.%m.'))
     ax1.xaxis.set_major_locator(mdates.MonthLocator(bymonthday=[1]))
     plt.tick_params(axis='both', labelsize=16) 
     ax1.set_xlabel('2017', fontsize=16) 
     ax1.set_title(title_list[i], fontsize=16)
     ax1.grid() 
     if i==0: 
         ax1.legend(labels=['T2Mean','T2Max','T2Min'],loc='lower right', fontsize=14)

wspace = 0.7   # the amount of width reserved for blank space between subplots
fig.subplots_adjust(wspace=wspace)
#plt.show()
#plt.savefig(str(dir_out)+'/'+str(exp_number_irri)+'_temperature_AMJJA_2.png',dpi=300, bbox_inches='tight')

