#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 16:53:56 2023

@author: g300099

appendix station location for comparison with observation data
"""


import xarray as xr 
import numpy as np
import pandas as pd
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import glob

from functions_reading_files import *
from functions_correcting_time import * 
from functions_plotting import * 
from functions_calculations import *
from functions_rotation import *
from functions_idw import *

varlist=['T2M_MEAN','T2M_MAXMEAN','T2M_MINMEAN']
varremolist=['TEMP2','T2MAX','T2MIN']

for var, varremo in zip(varlist, varremolist): 
    print('Calculating variable:', varremo)
    
    
    # read station data
    
    if varremo=='TEMP2':
        filetype='mfiles'
    elif varremo=='T2MAX' or varremo== 'T2MIN': 
        filetype='efiles'
    else: 
        print('ATTENTION: No filetype set! Check var.')
    
    
    station_dir='/work/ch0636/g300099/data/SCIA/'
    station_file='stations.csv'
    data_dir='/work/ch0636/g300099/data/SCIA/'
    if varremo=='TEMP2' and var=='T2M_MEAN': 
        temp_dir='T2M_MEAN/'
        temp_files=glob.glob(data_dir+temp_dir+'tempmean*.csv')
    elif varremo=='T2MAX' and var=='T2M_MAXMEAN': 
        temp_dir='T2M_MAXMEAN/'
        temp_files=glob.glob(data_dir+temp_dir+'tempmaxmean*.csv')
    elif varremo=='T2MIN' and var=='T2M_MINMEAN': 
        temp_dir='T2M_MINMEAN/'
        temp_files=glob.glob(data_dir+temp_dir+'tempminmean*.csv')
    else:
        print('Variable combination wrong!')
    
    
    station=pd.read_csv(station_dir+station_file)
    
    # create dataframe with station's position and variable
    data_files=temp_files
    
    data=[]
    for filename in data_files: 
        df=pd.read_csv(filename, index_col=None)
        data.append(df)
    
    df_data=pd.concat(data,axis=0,ignore_index=True)
    
    # station data | extract month data, merge lon & lat to data, renaming
    
    #month=8
    monthlist=[4,5,6,7,8]
    year=2017
    
    for month in monthlist:
        print('calculating month:',month)
        var_month=df_data[(df_data.mese==month) & (df_data.anno==year)]
        station.rename(columns={'nomestazione':'anagrafica','nomerete':'rete'}, inplace=True)
        
        station_data=pd.merge(var_month,station) 
        if varremo=='TEMP2' and var=='T2M_MEAN': 
            station_data.rename(columns={'rete':'network','anagrafica':'station_name','anno':'year',\
                                        'mese':'month','Temperatura media':'t2mean','2 dev standard':'std',\
                                        'ndati':'ndata','codice':'code','longitudine':'longitude',\
                                        'latitudine':'latitude'}, inplace=True)
        elif varremo=='T2MAX' and var=='T2M_MAXMEAN': 
                station_data.rename(columns={'rete':'network','anagrafica':'station_name','anno':'year',\
                                        'mese':'month','Temperatura massima (media)':'t2max','2 dev standard':'std',\
                                        'ndati':'ndata','codice':'code','longitudine':'longitude',\
                                        'latitudine':'latitude'}, inplace=True)
        elif varremo=='T2MIN' and var=='T2M_MINMEAN': 
                station_data.rename(columns={'rete':'network','anagrafica':'station_name','anno':'year',\
                                        'mese':'month','Temperatura minima (media)':'t2min','2 dev standard':'std',\
                                        'ndati':'ndata','codice':'code','longitudine':'longitude',\
                                    'latitudine':'latitude'}, inplace=True)
        else: 
            print('ATTENTION: var is wrong set!')
        station_data.drop_duplicates(ignore_index=True)
        
        #  remo data | read files
        
        exp_number_irri='067016' 
        exp_number_noirri='067015'
        
        # background map: irrifrac
        remo_dir='/work/ch0636/g300099/SIMULATIONS/GAR11/remo_results/067016/2017/xt/'
        remo_files='e067016t20170'+str(month)+'0100.nc'
        remo_irrifrac = xr.open_dataset(remo_dir+remo_files)
        irrifrac=remo_irrifrac.IRRIFRAC[0]
        
        #remo data 
        if varremo=='T2MIN' and var=='T2M_MINMEAN': 
            remo_irri_var= ((read_efiles('T2MIN', 202, exp_number_irri, year, month)[varremo][:,0,:,:])-273.15).drop('height2m')
            remo_irri_var_newtime = correct_timedim(remo_irri_var.to_dataset())
            remo_irri_daily_ext=remo_irri_var_newtime.groupby("time.day").min().mean(dim='day')[varremo]
            remo_irri_data=xr.merge([irrifrac,remo_irri_daily_ext],compat='override')
           
            remo_noirri_var= ((read_efiles('T2MIN', 202, exp_number_noirri, year, month)[varremo][:,0,:,:])-273.15).drop('height2m') 
            remo_noirri_var_newtime = correct_timedim(remo_noirri_var.to_dataset())
            remo_noirri_daily_ext=remo_noirri_var_newtime.groupby("time.day").min().mean(dim='day')[varremo]
            remo_noirri_data=xr.merge([irrifrac,remo_noirri_daily_ext],compat='override')
        
        elif varremo=='T2MAX' and var=='T2M_MAXMEAN': 
            remo_irri_var= ((read_efiles('T2MAX', 201, exp_number_irri, year, month)[varremo][:,0,:,:])-273.15).drop('height2m') 
            remo_irri_var_newtime = correct_timedim(remo_irri_var.to_dataset())
            remo_irri_daily_ext=remo_irri_var_newtime.groupby("time.day").max().mean(dim='day')[varremo]
            remo_irri_data=xr.merge([irrifrac,remo_irri_daily_ext],compat='override')
            
            remo_noirri_var= ((read_efiles('T2MAX', 201, exp_number_noirri, year, month)[varremo][:,0,:,:])-273.15).drop('height2m') 
            remo_noirri_var_newtime = correct_timedim(remo_noirri_var.to_dataset())
            remo_noirri_daily_ext=remo_noirri_var_newtime.groupby("time.day").max().mean(dim='day')[varremo]
            remo_noirri_data=xr.merge([irrifrac,remo_noirri_daily_ext],compat='override')
            
        elif varremo=='TEMP2' and var=='T2M_MEAN': 
            remo_irri_var=((read_mfiles(exp_number_irri, year, month)[varremo][0,0,:,:])-273.15).drop('height2m')
            remo_noirri_var=((read_mfiles(exp_number_noirri, year, month)[varremo][0,0,:,:])-273.15).drop('height2m')
            remo_irri_data=xr.merge([irrifrac,remo_irri_var],compat='override')
            remo_noirri_data=xr.merge([irrifrac,remo_noirri_var],compat='override')
        
        else:
            print('Your variable is wrong! Var:',varremo)
            
        
        # rotate station coordinates 
        # delete duplicates, one station hats2 different entries as mean values (t2mean)--> had to be deleted
        if var=='T2M_MEAN':
           station_data_sel=station_data[['t2mean','latitude','longitude']].drop_duplicates(ignore_index=True)
        elif var=='T2M_MINMEAN': 
           station_data_sel=station_data[['t2min','latitude','longitude']].drop_duplicates(ignore_index=True)
        elif var=='T2M_MAXMEAN': 
           station_data_sel=station_data[['t2max','latitude','longitude']].drop_duplicates(ignore_index=True)
        else: 
            print('Station data is not cleaned from duplicates.')
        
        np_lat=39.25 
        np_lon=-162             
        
        stations_rot=[]
        rot_longitude=[]
        rot_latitude=[]
        for i in range(len(station_data_sel.longitude)): 
            stations_rot.append(rotated_coord_transform(station_data_sel.longitude.values[i], \
                            station_data_sel.latitude.values[i],np_lon,np_lat, direction='geo2rot'))
            rot_longitude.append(round(stations_rot[i][0],5))
            rot_latitude.append(round(stations_rot[i][1],5))
        station_data_sel['rlon']=rot_longitude
        station_data_sel['rlat']=rot_latitude
        
        # clean from duplicates in rotated grid 
        station_data_corr=station_data_sel[(station_data_sel.rlon!=station_data_sel.rlon.value_counts().loc[lambda x: x>1].index[0])\
                                &(station_data_sel.rlat!=station_data_sel.rlat.value_counts().loc[lambda x: x>1].index[0])].reset_index()
        df_stations=station_data_corr.set_index(['rlat','rlon'])
        da_stations=df_stations.to_xarray().astype('float64')
        
        
        # filter data with irrigated fraction 
        
        remo_irri_data_cut=remo_irri_data.isel(rlat=slice(50,70),rlon=slice(60,110))
        remo_noirri_data_cut=remo_noirri_data.isel(rlat=slice(50,70),rlon=slice(60,110))
        irrifrac=remo_irri_data.IRRIFRAC[50:70,60:110]
        
        # find good threshold:
        # 1. filter all grid cells with irrigated fraction > 0 
        # 2. filter all grid cells with irrigated fraction > mean(irrifrac)
        remo_irri_data_irrifrac=remo_irri_data_cut.where(remo_irri_data_cut.IRRIFRAC>0)
        remo_noirri_data_irrifrac=remo_noirri_data_cut.where(remo_noirri_data_cut.IRRIFRAC>0)
        #print('Mean irrifrac:', remo_irri_data_irrifrac.IRRIFRAC.mean().values)
        
        # remo_irri and remo_no_irri are the same in irrifrac --> come up with same result for filtering gridcells & stations
        remo_irri_data_filtered=remo_irri_data_cut.where(remo_irri_data_cut.IRRIFRAC>remo_irri_data_irrifrac.IRRIFRAC.mean().values)
        remo_noirri_data_filtered=remo_noirri_data_cut.where(remo_noirri_data_cut.IRRIFRAC>remo_noirri_data_irrifrac.IRRIFRAC.mean().values)
        
        # IDW interpolation for model data to station location 
        
        # station positions as target positions
        rlat_target=station_data_corr.rlat.values
        rlon_target=station_data_corr.rlon.values
        
        # model grid cells as source 
        rlat_source=remo_irri_data_filtered.rlat.values 
        rlon_source=remo_irri_data_filtered.rlon.values
        data_irri=remo_irri_data_filtered[varremo].values
        data_noirri=remo_noirri_data_filtered[varremo].values
        
        # applying IDW 
        df_inttemp_irri=idw_for_model(rlat_source, rlon_source, rlat_target, rlon_target, data_irri)
        df_inttemp_noirri=idw_for_model(rlat_source, rlon_source, rlat_target, rlon_target, data_noirri)
        
        df_inttemp_irri_clean=df_inttemp_irri.dropna(axis=0).reset_index(drop=True)
        
        # count used values
        station_numb=np.sum(~np.isnan(df_inttemp_irri.temp2m.values))
        print('station number:',station_numb)

    irrifrac=irrifrac.where(irrifrac>0)
    
    fig1=plt.figure(figsize=(10,3))
    ax1=fig1.add_subplot(1,1,1,projection=rotated_pole) 
    ax1.add_feature(land,zorder=1)
    ax1.add_feature(coastline,edgecolor='black')
    ax1.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax1.set_xmargin(0)
    ax1.set_ymargin(0) 
    levels=np.arange(0,1.1,0.1)
    ticks=levels[::2].round(2)
    rotplot=irrifrac.plot.pcolormesh(ax=ax1,cmap='YlGnBu',levels=levels,add_colorbar=False)
    for k in range(len(df_inttemp_irri_clean)):
       plt.scatter(x=df_inttemp_irri_clean.rlon[k],y=df_inttemp_irri_clean.rlat[k], s=20, facecolors='darkorange', alpha=1, edgecolor='darkorange',\
                  transform=rotated_pole)
    cbar=plt.colorbar(rotplot, ticks=ticks)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(ticks)
    cbar.ax.set_title('fraction [-]', fontsize=15, pad=17)
    cbar.ax.tick_params(labelsize=15)
    ax1.gridlines()
    plt.title('SCIA station locations selected \n'+str(varremo),fontsize=15)
    plt.tight_layout()
    #plt.savefig('/work/ch0636/g300099/EVALUATION/plots/Eval_Scia_stations_location_selected_'+str(var)+'_new.png',dpi=300, bbox_inches='tight')
     
    
       