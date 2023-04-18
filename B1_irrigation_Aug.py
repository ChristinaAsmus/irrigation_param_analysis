#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 16:51:49 2023

@author: g300099
"""

import xarray as xr
import matplotlib.pyplot as plt

from functions_reading_files import *
from functions_correcting_time import * 
from functions_plotting import * 

# In[]: appendix irrimask for August

year=2017
exp_number_irri='067016'


varlist=['IRRIMASK']
var_num_list=['797']
month=8



for var, var_num in zip(varlist, var_num_list): 
    single_var_data=read_efiles(var, var_num, exp_number_irri, year, month)
    if var==varlist[0]: 
        ds_var_irri=single_var_data
    else:
        ds_var_irri=xr.merge([ds_var_irri, single_var_data])
ds_var_irri=xr.merge([ds_var_irri, irrifrac])
dsirr_newtime = correct_timedim(ds_var_irri)

month_name='June'
nirrihours = 9  # da in den Output Daten IRRIMASK nur für 5h aktiv ist
imasksum = (dsirr_newtime.IRRIMASK.sum('time')/nirrihours)
imaskvalues = imasksum.where(imasksum > 0, np.nan)


params = {'legend.fontsize':14,
         'axes.labelsize': 16,
         'axes.titlesize':18,
         'xtick.labelsize':14,
         'ytick.labelsize':14}
plt.rcParams.update(params)

fig = plt.figure(figsize=(12, 7))
ax = fig.add_subplot(1, 1, 1, projection=rotated_pole)
cax = fig.add_axes([0.265, 0.05, 0.5, 0.04])
levels=np.arange(0,31,2)
ticks=levels
rotplot = plot_rotvar(fig, imaskvalues, ax, cax, '[n]', 'number of irrigated days', 'GnBu', levels, ticks, 'max','horizontal')
                       
#plt.savefig(str(dir_out)+'/app_irrimask_'+str(month_name)+'.png',dpi=300, bbox_inches='tight')