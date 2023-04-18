#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 14:44:36 2023

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

# In[]: read mfiles
mds_irri=read_mfiles(exp_number_irri, 2017,0)
mds_noirri=read_mfiles(exp_number_noirri, 2017,0)

mdsirri_extended = correct_timedim_mfiles(xr.merge([mds_irri, irrifrac]))
mdsnoirri_extended = correct_timedim_mfiles(xr.merge([mds_noirri, irrifrac]))

# In[]: define plot directory

dir_out='/work/ch0636/g300099/SIMULATIONS/GAR11/plot/plot_for_paper1/'
    
# In[]: 
# LAI 

irrilimit=0.0

var='ALAI_PFI'


# monthly for irrigated fraction
data_irri=(mdsirri_extended.where(mdsirri_extended.IRRIFRAC>irrilimit)[var][:,0,:,:].mean(dim=['rlon','rlat'], skipna=True))
data_noirri=(mdsnoirri_extended.where(mdsnoirri_extended.IRRIFRAC>irrilimit)[var][:,0,:,:].mean(dim=['rlon','rlat'], skipna=True))
laimonthdiff=data_irri-data_noirri

#
month_list=['April','May','June']
mdsirri_AMJ=mdsirri_extended.sel(time=mdsirri_extended.time.dt.month.isin([4,5,6]))
mdsnoirri_AMJ=mdsnoirri_extended.sel(time=mdsnoirri_extended.time.dt.month.isin([4,5,6]))
laispatialdiff=mdsirri_AMJ.ALAI_PFI[:,0,:,:]-mdsnoirri_AMJ.ALAI_PFI[:,0,:,:]

#NPP
irrilimit=0.0
var='ANPP_ACI'

# monthly for irrigated fraction

data_irri=(mdsirri_extended.where(mdsirri_extended.IRRIFRAC>irrilimit)[var][:,0,:,:].mean(dim=['rlon','rlat'], skipna=True))*(30*24)
data_noirri=(mdsnoirri_extended.where(mdsnoirri_extended.IRRIFRAC>irrilimit)[var][:,0,:,:].mean(dim=['rlon','rlat'], skipna=True))*(30*24)
nppmonthdiff=data_irri-data_noirri
# we still have to change the unit! 
# ANPP_ACI is in gCm2/s ich brauch gCm2/month, Richtwert: Indien max. 100 in monsoon time


month_list=['April','May','June']
mdsirri_AMJ=(mdsirri_extended.sel(time=mdsirri_extended.time.dt.month.isin([4,5,6])))
mdsnoirri_AMJ=(mdsnoirri_extended.sel(time=mdsnoirri_extended.time.dt.month.isin([4,5,6])))
nppspatialdiff=((mdsirri_AMJ[var][:,0,:,:])*(30*24))-((mdsnoirri_AMJ[var][:,0,:,:])*(30*24))
# line plot  

fig = plt.figure(figsize=(16, 12))

params = {'legend.fontsize':20,
         'axes.labelsize': 20,
         'axes.titlesize':18,
         'xtick.labelsize':20,
         'ytick.labelsize':20}
plt.rcParams.update(params)

ax1 = fig.add_subplot(2, 2, 1)
laimonthdiff.plot.line(ax=ax1, marker='.')
ax1.set_xticks(ticks=laimonthdiff.time.values, labels=laimonthdiff.time.dt.strftime('%b').values)
data_irri.plot.line(ax=ax1,label='irrigated', add_legend=False)
data_noirri.plot.line(ax=ax1, label='not irrigated', add_legend=False)
ax1.set_xlabel('months',) 
ax1.set_ylabel('LAI [m$^2$m$^{-2}$]')
ax1.grid(True)
ax1.set_title(' ')
ax1.legend(loc='upper right')
ax1.tick_params(axis='both')

ax2 = fig.add_subplot(2, 2, 2)
ax2.set_xticks(ticks=laimonthdiff.time.values, labels=laimonthdiff.time.dt.strftime('%b').values)
laimonthdiff.plot.line(ax=ax2, marker='.')
ax2.set_ylim(-0.2,0.2)
ax2.set_xlabel('months') 
ax2.set_ylabel('LAI [Δm$^2$m$^{-2}$]')
ax2.grid(True)
ax2.set_title(' ')
ax2.tick_params(axis='both')
plt.tight_layout()

ax3 = fig.add_subplot(2, 2, 3)
ax3.sharex(ax1)
ax3.set_xticks(ticks=nppmonthdiff.time.values, labels=nppmonthdiff.time.dt.strftime('%b').values)
data_irri.plot.line(ax=ax3,label='irrigated', add_legend=False)
data_noirri.plot.line(ax=ax3, label='not irrigated', add_legend=False)
ax3.set_xlabel('months') 
ax3.set_ylabel('NPP [gCm$^{-2}$month$^{-1}$]')
ax3.grid(True)
ax3.set_title(' ')
ax3.legend(loc='upper right')
ax3.tick_params(axis='both')

ax4 = fig.add_subplot(2, 2, 4)
ax4.set_xticks(ticks=nppmonthdiff.time.values, labels=nppmonthdiff.time.dt.strftime('%b').values)
nppmonthdiff.plot.line(ax=ax4, marker='.')
ax4.set_xlabel('months') 
ax4.set_ylabel('NPP [ΔgCm$^{-2}$month$^{-1}$]')
ax4.grid(True)
ax4.set_title(' ')
ax4.tick_params(axis='both')
plt.tight_layout()
#plt.savefig(str(dir_out)+'/NPP_lines_months_diff.png',dpi=300, bbox_inches='tight')

