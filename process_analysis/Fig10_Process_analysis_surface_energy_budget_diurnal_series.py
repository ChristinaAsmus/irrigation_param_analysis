#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 14:04:53 2023

@author: g300099
"""


import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sys
sys.path.append('/home/g/g300099/pyprograms/Paper1_for_publishing_final/plot_figures/functions/') 
from functions_reading_files import *
from functions_correcting_time import * 
from functions_plotting import * 
from functions_calculations import *

# In[]: select experiment, year and month

year=2017
month=0
exp_number_irri='067016'
exp_number_noirri='067015'


# In[]: read files 

# background map: irrifrac
remo_dir = '/work/ch0636/g300099/SIMULATIONS/GAR11/remo_results/067016/2017/xt/'
remo_files = 'e067016t2017060100.nc'
remo_tfile = xr.open_dataset(remo_dir+remo_files)
irrifrac=remo_tfile.IRRIFRAC[0]

varlist=['SRADS','TRADS','AHFLIRR','AHFSIRR']
var_num_list=['176','177', '736','734']

for var, var_num in zip(varlist, var_num_list): 
    single_var_data_irri=read_efiles(var, var_num, exp_number_irri, year, month)
    single_var_data_noirri=read_efiles(var, var_num, exp_number_noirri, year, month)
    if var==varlist[0]: 
        ds_var_irri=single_var_data_irri
        ds_var_noirri=single_var_data_noirri
    else:
        ds_var_irri=xr.merge([ds_var_irri, single_var_data_irri])
        ds_var_noirri=xr.merge([ds_var_noirri, single_var_data_noirri])
ds_var_irri=xr.merge([ds_var_irri, irrifrac])
ds_var_noirri=xr.merge([ds_var_noirri, irrifrac])

dsirr_newtime = correct_timedim(ds_var_irri)
dsnoirr_newtime = correct_timedim(ds_var_noirri)

# In[]: read mfiles
mds_irri=read_mfiles(exp_number_irri, 2017,0)
mds_noirri=read_mfiles(exp_number_noirri, 2017,0)

mdsirri_extended = correct_timedim_mfiles(xr.merge([mds_irri, irrifrac]))
mdsnoirri_extended = correct_timedim_mfiles(xr.merge([mds_noirri, irrifrac]))

# In[]: define plot directory

dir_out='/work/ch0636/g300099/SIMULATIONS/GAR11/plot/plot_for_paper1/'
    
# In[]: Define regions 

rlat_list=[(50,70),(30,70),(90,105)]
rlon_list=[(60,110),(1,40),(30,50)]
title_list=['IT','SF','CF']

# In[]: surface energy balance 
irrilimit=0.0
#use only fully irrigated months AMJ
dsirr_newtime_AMJ=dsirr_newtime.sel(time=dsirr_newtime.time.dt.month.isin([4, 5, 6]))
dsnoirr_newtime_AMJ=dsnoirr_newtime.sel(time=dsnoirr_newtime.time.dt.month.isin([4, 5, 6]))

sbalance_varlist=['SRADS','TRADS', 'AHFLIRR', 'AHFSIRR']

labels= {'SRADS':'shortwave radiation','TRADS':'longwave radiation','RN':'net radiation','AHFSIRR':'sensible heat flux',\
         'AHFLIRR':'latent heat flux','GHFL':'ground heat flux*'}
colors= {'SRADS': 'gold', 'TRADS':'red','AHFSIRR':'green','AHFLIRR':'blue','GHFL':'brown', 'RN':'black'}


# plot over all regions absolut and diff subplots
fig = plt.figure(figsize=(18, 4))

params = {'legend.fontsize':18,
         'axes.labelsize': 18,
         'axes.titlesize':18,
         'xtick.labelsize':16,
         'ytick.labelsize':16}
plt.rcParams.update(params)

for i in range(len(rlat_list)):
    irri_df=pd.DataFrame()
    noirri_df=pd.DataFrame()
    for var in sbalance_varlist: 
        # we have to change the sign for heat fluxes because it is negative for upward in the model
        if var=='AHFLIRR' or var=='AHFSIRR': 
             irri_df[var]=-(prepare_var(var,i,irrilimit, rlat_list, rlon_list, dsirr_newtime_AMJ, dsnoirr_newtime_AMJ)[0])
             noirri_df[var]=-(prepare_var(var,i,irrilimit, rlat_list, rlon_list, dsirr_newtime_AMJ, dsnoirr_newtime_AMJ)[1])
        else:    
            irri_df[var]=prepare_var(var,i,irrilimit,  rlat_list, rlon_list, dsirr_newtime_AMJ, dsnoirr_newtime_AMJ)[0]
            noirri_df[var]=prepare_var(var,i,irrilimit,  rlat_list, rlon_list, dsirr_newtime_AMJ, dsnoirr_newtime_AMJ)[1]
    irri_df['RN']=(irri_df.SRADS+irri_df.TRADS)
    noirri_df['RN']=(noirri_df.SRADS+noirri_df.TRADS)
    #calculation of ground heat flux GFHL as residuum
    irri_df['GHFL']=(irri_df.SRADS+irri_df.TRADS)-(irri_df.AHFSIRR+irri_df.AHFLIRR)
    noirri_df['GHFL']=(noirri_df.SRADS+noirri_df.TRADS)-(noirri_df.AHFSIRR+noirri_df.AHFLIRR)
    #
    df_diff_surface=(irri_df-noirri_df)
    #
    ax3 = fig.add_subplot(1, 3, (i+1))
    df_diff_surface[['RN','AHFLIRR','AHFSIRR','GHFL']].plot.line(ax=ax3,legend=False, color=[colors.get(x) for x in df_diff_surface[['RN','AHFLIRR','AHFSIRR','GHFL']].columns],\
       marker='.')
    ax3.set_title(title_list[i])
    ax3.set_xlabel('hour') 
    ax3.set_ylabel('[ΔWm$^{-2}$]')
    ax3.set_ylim(-200,250)
    ax3.grid(True)
    ax3.tick_params(axis='both')
ax3.legend([labels.get(x) for x in df_diff_surface[['RN','AHFLIRR','AHFSIRR','GHFL']].columns], loc='upper right',  bbox_to_anchor=(2.0, 1.05))
plt.tight_layout()
#plt.savefig(str(dir_out)+'/surface_energy_balance_hourly_regions_new.png',dpi=300, bbox_inches='tight')


