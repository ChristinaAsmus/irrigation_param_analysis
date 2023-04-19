#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
creates Fig 8 
creates combined subplot of spatial overview, diurnal cycle and annual cycle of irrigation effects 
for soil moisture, surface&soil temperature and evapotranspiration 
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

varlist=['WSECHIRR','TSECHIRR','TD3IRR','TD4IRR','TD5IRR','TDIRR','TDCLIRR','EVAP',\
         'EVAPIRR','ESKIN','ETRANS','EBSOIL']
var_num_list=['701','732','721','722','723','715','716','182',\
              '739','302','012','013']

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
    

# In[]: soil and surface variables
irrilimit=0.0
# soil moisture
#use only fully irrigated months AMJ
dsirr_newtime_AMJ=dsirr_newtime.sel(time=dsirr_newtime.time.dt.month.isin([4, 5, 6]))
dsnoirr_newtime_AMJ=dsnoirr_newtime.sel(time=dsnoirr_newtime.time.dt.month.isin([4, 5, 6]))

#spatial plot for mean of AJM 
wsspatialdiff=calculate_meandiff('WSECHIRR',dsirr_newtime_AMJ, dsnoirr_newtime_AMJ)
wshourdiff=calculate_hourmeans_diff('WSECHIRR',irrilimit, dsirr_newtime_AMJ, dsnoirr_newtime_AMJ)
wsmonthdiff=calculate_monthmeans_diff('WSECHIRR',irrilimit, mdsirri_extended, mdsnoirri_extended)


# soil temperatures 
#spatial plot for mean of AJM 
stempspatialdiff=calculate_meandiff('TSECHIRR',dsirr_newtime_AMJ, dsnoirr_newtime_AMJ)

# consider all levels of soil temperature for line plot
depths_str=['0.0 m', '-0.065 m', '-0.319 m', '-1.232 m', '-4.134 m']
stemp_list=['TD3IRR','TD4IRR','TD5IRR','TDIRR','TDCLIRR']

month_index=np.arange(dsirr_newtime.time.dt.month.values.min(),dsirr_newtime.time.dt.month.values.max()+1)
hour_index=np.arange(dsirr_newtime_AMJ.time.dt.hour.values.min(),dsirr_newtime_AMJ.time.dt.hour.values.max()+1)

stemp_month_df=pd.DataFrame()
stemp_month_df['month']=month_index
stemp_month_df.index=stemp_month_df.month
stemp_month_df.drop(['month'],axis=1, inplace=True)

stemp_day_df=pd.DataFrame()
stemp_day_df['hour']=hour_index
stemp_day_df.index=stemp_day_df.hour
stemp_day_df.drop(['hour'],axis=1, inplace=True)


for var in stemp_list: 
    stemp_day_df[var]=calculate_hourmeans_diff(var,irrilimit, dsirr_newtime_AMJ, dsnoirr_newtime_AMJ)
    stemp_month_df[var]=calculate_monthmeans_diff(var,irrilimit, mdsirri_extended, mdsnoirri_extended)

# evapotranspiration 
# spatial plot for mean AMJ
evapospatialdiff=-(calculate_monthlysum_mean_evapos('EVAPIRR',irrilimit, \
                                      dsirr_newtime_AMJ, dsnoirr_newtime_AMJ))
# consider all all evapotranspiration fractions
evap_list=['EVAP','EVAPIRR','ESKIN','ETRANS','EBSOIL']

# combine evapotranspiration fractions in 2 dataframes (diurnal and annual cycle)
month_index=np.arange(dsirr_newtime.time.dt.month.values.min(),dsirr_newtime.time.dt.month.values.max()+1)
hour_index=np.arange(dsirr_newtime_AMJ.time.dt.hour.values.min(),dsirr_newtime_AMJ.time.dt.hour.values.max()+1)

evap_month_df=pd.DataFrame()
evap_month_df['month']=month_index
evap_month_df.index=evap_month_df.month

evap_day_df=pd.DataFrame()
evap_day_df['hour']=hour_index
evap_day_df.index=evap_day_df.hour

for var in evap_list: 
    evap_month_df[var]=-(calculate_monthlysum_evapos(var, irrilimit, dsirr_newtime,\
                                                dsnoirr_newtime))
    evap_day_df[var]=-(calculate_hourlymean_evapos(var, irrilimit, \
                dsirr_newtime_AMJ, dsnoirr_newtime_AMJ))

evap_month_df['EVAP_all']=evap_month_df['EVAP']+evap_month_df['EVAPIRR']
evap_month_df.drop(['EVAP','EVAPIRR','month'], axis=1, inplace=True)

evap_day_df['EVAP_all']=evap_day_df['EVAP']+evap_day_df['EVAPIRR']
evap_day_df.drop(['EVAP','EVAPIRR','hour'], axis=1, inplace=True)
 

labels=['evaporation of \n skin reservoir','transpiration \n of vegetation',\
        'evaporation \n of bare soil', \
        'evapotranspi- \n ration summed']

# plot 3x3 subplots  

xticklabels=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

fig=plt.figure(figsize=(18,18))

ax1 = fig.add_subplot(3, 3, 1, projection=rotated_pole)
cax1 = fig.add_axes([0.138, 0.655, 0.18, 0.015])
levels=[-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6]
ticks=[-0.6,-0.4,-0.2,0,0.2,0.4,0.6]
rotplot = plot_rotvar(fig,wsspatialdiff, ax1,cax1, label='soil moisture [Δm]', unit='[Δm]',
    cmap='RdBu_r', levels=levels,  ticks=ticks,extend_scale='both', cbar_orient='horizontal')

ax2 = fig.add_subplot(3, 3, 2)
wshourdiff.plot.line(ax=ax2, marker='.')
ax2.set_xlabel('hours') 
ax2.set_ylabel('[Δm]')
ax2.grid(True)

ax3 = fig.add_subplot(3, 3, 3)
wsmonthdiff.plot.line(ax=ax3, marker='.')
ax3.set_xticks(ticks=wsmonthdiff.time.values, labels=wsmonthdiff.time.dt.strftime('%b').values)
ax3.set_xticks(ax3.get_xticks()[::2])
ax3.set_xlabel('months') 
ax3.set_ylabel('[Δm]')
ax3.grid(True)


ax4 = fig.add_subplot(3, 3, 4, projection=rotated_pole)
cax4 = fig.add_axes([0.138, 0.376, 0.18, 0.015])
levels=[-4,-3,-2,-1,1,2,3,4]
ticks=[-4,-2,0,2,4]
rotplot = plot_rotvar(fig,stempspatialdiff, ax4, cax4, label='surface temperature [ΔK]', unit='[ΔK]',
    cmap='RdBu_r', levels=levels,ticks=ticks, extend_scale='both', cbar_orient='horizontal' )

ax5 = fig.add_subplot(3, 3, 5)
stemp_day_df.plot(ax=ax5, legend=False, marker='.')
ax5.set_xticks(stemp_day_df.index[0::5], minor=False)
ax5.set_xlabel('hours') 
ax5.set_ylabel('[ΔK]')
ax5.grid(True)

ax6 = fig.add_subplot(3, 3, 6)
stemp_month_df.plot(ax=ax6, legend=False, marker='.')
ax6.set_xticks(stemp_month_df.index, minor=False)
ax6.set_xticklabels(xticklabels, rotation=45)
ax6.set_xticks(ax6.get_xticks()[::2])
ax6.set_xlabel('month') 
ax6.set_ylabel('[ΔK]')
ax6.legend(depths_str, bbox_to_anchor=(1.1, 1))
ax6.grid(True)


ax7 = fig.add_subplot(3, 3, 7, projection=rotated_pole)
cax7 = fig.add_axes([0.138, 0.095, 0.18, 0.015])
levels=[-150, -125, -100, -75, -50, -25, 25, 50, 75, 100,  125, 150]
ticks=[-150, -100, -50, 0,50, 100, 150]
rotplot = plot_rotvar(fig, evapospatialdiff, ax7, cax7, label='evapotranspiration [Δmm]',\
   unit='[Δmm]', cmap='RdBu_r', levels=levels,  ticks=ticks,extend_scale='both', cbar_orient='horizontal')

ax8 = fig.add_subplot(3, 3, 8)
evap_day_df.plot(ax=ax8, legend=False, marker='.')
ax8.set_xticks(evap_day_df.index[0::5], minor=False)

ax8.set_xlabel('hours') 
ax8.set_ylabel('[Δmm]')
ax8.grid(True)

ax9 = fig.add_subplot(3, 3, 9)
evap_month_df.plot(ax=ax9, legend=False, marker='.')
ax9.set_xticks(evap_month_df.index, minor=False)
ax9.set_xticklabels(xticklabels, rotation=45)
ax9.set_xticks(ax9.get_xticks()[::2])
ax9.set_xlabel('month') 
ax9.set_ylabel('[Δmm]')
ax9.legend(labels, bbox_to_anchor=(1.1, 1))
ax9.grid(True)
#
wspace = 0.35  
hspace = 0.45  
fig.subplots_adjust(wspace=wspace, hspace=hspace)
#plt.savefig(str(dir_out)+'/'+str(exp_number_irri)+'soil_processes_combined_new.png',dpi=300, bbox_inches='tight')
plt.show()
