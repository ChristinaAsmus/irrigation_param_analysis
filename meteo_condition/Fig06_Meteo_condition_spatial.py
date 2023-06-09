#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 11:39:17 2023

@author: g300099
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from analysis_functions.functions_correcting_time import correct_timedim
from analysis_functions.functions_plotting import plot_rotvar, rotated_pole
from analysis_functions.functions_reading_files import read_efiles

# In[]: experiment number and year

year = 2017
# for all months month=0
month = 0

exp_number_irri = "067016"
exp_number_noirri = "067015"


# In[]: read efiles

# paths to the data
data_path = "../data"

# background map: irrifrac
remo_irrifrac_dir = str(exp_number_noirri) + "/irrifrac/"
remo_irrifrac_files = "e" + str(exp_number_noirri) + "e_c743_201706.nc"
remo_irrifrac = xr.open_dataset(
    str(data_path) + "/" + str(remo_irrifrac_dir) + str(remo_irrifrac_files)
)
irrifrac = remo_irrifrac.IRRIFRAC[0]
# land-sea-mask
remo_landfrac_dir = str(exp_number_noirri) + "/bla/"
remo_landfrac_files = "e" + str(exp_number_noirri) + "e_c172_201706.nc"
remo_landfrac = xr.open_dataset(
    str(data_path) + "/" + str(remo_landfrac_dir) + str(remo_landfrac_files)
)
landfrac = remo_landfrac.BLA[0]

# analysis variables
varlist = [
    "WSECHIRR",
    "APRL",
    "APRC",
    "TEMP2",
    "T2MAX",
    "T2MIN",
]
var_num_list = ["701", "142", "143", "167", "201", "202"]

for var, var_num in zip(varlist, var_num_list):
    single_var_data = read_efiles(
        data_path, var, var_num, exp_number_noirri, year, month
    )
    if var == varlist[0]:
        ds_var_noirri = single_var_data
    else:
        ds_var_noirri = xr.merge([ds_var_noirri, single_var_data])
ds_var_noirri = xr.merge([ds_var_noirri, irrifrac, landfrac])
dsnoirr_newtime = correct_timedim(ds_var_noirri)

# In[]: define plot directory


dir_working = os.getcwd()
# creates dir in parent directory
dir_out = os.path.join(os.pardir, "Figures")
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
print("Output directory is: ", dir_out)
# In[]: Define regions

rlat_list = [(50, 70), (30, 70), (90, 105)]
rlon_list = [(60, 110), (1, 40), (30, 50)]
title_list = ["IT", "SF", "CF"]


# In[]: spatial plot for overview AMJ

# select months of interest AMJ
dsnoirr_newtime_AMJ = dsnoirr_newtime.sel(
    time=dsnoirr_newtime.time.dt.month.isin([4, 5, 6])
)
lsm = dsnoirr_newtime_AMJ.BLA

ws_mean_AMJ = (dsnoirr_newtime_AMJ.WSECHIRR).mean("time")
precip_mean_AMJ = (
    (dsnoirr_newtime_AMJ.APRL + dsnoirr_newtime_AMJ.APRC)
    .groupby("time.month")
    .sum()
    .mean("month")
)
temp2_mean_AMJ = dsnoirr_newtime_AMJ.TEMP2.mean("time")[0, :, :] - 273.15

ws_mean_AMJ = ws_mean_AMJ.where(irrifrac > 0)
temp2_mean_AMJ = temp2_mean_AMJ.where(lsm > 0)
precip_mean_AMJ = precip_mean_AMJ.where(lsm > 0)


# Heatwave (3-5Aug)
dsnoirr_newtime_Aug = dsnoirr_newtime.sel(
    time=(
        dsnoirr_newtime.time.dt.month.isin([8])
        & dsnoirr_newtime.time.dt.day.isin(np.arange(3, 6, 1))
    )
)

temp2_mean_Aug = dsnoirr_newtime_Aug.TEMP2.mean("time")[0, :, :] - 273.15
ws_mean_Aug = (dsnoirr_newtime_Aug.WSECHIRR).mean("time")
precip_mean_Aug = (dsnoirr_newtime_Aug.APRL + dsnoirr_newtime_Aug.APRC).sum("time")

ws_mean_Aug = ws_mean_Aug.where(irrifrac > 0)
temp2_mean_Aug = temp2_mean_Aug.where(lsm > 0)
precip_mean_Aug = precip_mean_Aug.where(lsm > 0)

# In[]: plot
fig = plt.figure(figsize=(18, 12))


params = {
    "legend.fontsize": 18,
    "legend.markerscale": 15,
    "axes.labelsize": 18,
    "axes.titlesize": 18,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
}
plt.rcParams.update(params)

# temperature
ax1 = fig.add_subplot(2, 3, 1, projection=rotated_pole)
cax1 = fig.add_axes([0.127, 0.53, 0.22, 0.02])
levels1 = np.arange(0, 45, 2.5)
ticks1 = [round(x) for x in levels1][::2]
plot_rotvar(
    fig,
    temp2_mean_AMJ,
    ax1,
    cax1,
    "[°C]",
    "2 m temperature [°C] ",
    "jet",
    levels1,
    ticks1,
    "max",
    "horizontal",
)
ax1.set_title(" ")

# soil moisture
ax2 = fig.add_subplot(2, 3, 2, projection=rotated_pole)
cax2 = fig.add_axes([0.4, 0.53, 0.22, 0.02])
levels2 = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
ticks2 = levels2[::2]
plot_rotvar(
    fig,
    ws_mean_AMJ,
    ax2,
    cax2,
    "[m]",
    "soil moisture [m] ",
    "YlGnBu",
    levels2,
    ticks2,
    "max",
    "horizontal",
)
ax2.set_title(" ")

# precipitation
ax3 = fig.add_subplot(2, 3, 3, projection=rotated_pole)
cax3 = fig.add_axes([0.677, 0.53, 0.215, 0.02])
levels3 = np.arange(0, 400, 25)
ticks3 = levels3[::3]
plot_rotvar(
    fig,
    precip_mean_AMJ,
    ax3,
    cax3,
    "[mm]",
    "precipitation [mm] ",
    "YlGnBu",
    levels3,
    ticks3,
    "max",
    "horizontal",
)


# temperature
ax4 = fig.add_subplot(2, 3, 4, projection=rotated_pole)
cax4 = fig.add_axes([0.127, 0.08, 0.22, 0.02])
levels4 = np.arange(0, 45, 2.5).round()
ticks4 = [round(x) for x in levels4][::2]
plot_rotvar(
    fig,
    temp2_mean_Aug,
    ax4,
    cax4,
    "[°C]",
    "2 m temperature [°C] ",
    "jet",
    levels4,
    ticks4,
    "max",
    "horizontal",
)
ax4.set_title(" ")

# soil moisture
ax5 = fig.add_subplot(2, 3, 5, projection=rotated_pole)
cax5 = fig.add_axes([0.4, 0.08, 0.22, 0.02])
levels5 = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
ticks5 = levels5[::2]
plot_rotvar(
    fig,
    ws_mean_Aug,
    ax5,
    cax5,
    "[m]",
    "soil moisture [m] ",
    "YlGnBu",
    levels5,
    ticks5,
    "max",
    "horizontal",
)
ax5.set_title(" ")

# precipitation
ax6 = fig.add_subplot(2, 3, 6, projection=rotated_pole)
cax6 = fig.add_axes([0.678, 0.08, 0.22, 0.02])
levels6 = np.arange(0, 400, 25)
ticks6 = levels6[::3]
plot_rotvar(
    fig,
    precip_mean_Aug,
    ax6,
    cax6,
    "[mm]",
    "precipitation [mm] ",
    "YlGnBu",
    levels6,
    ticks6,
    "max",
    "horizontal",
)
hspace = 0.5
fig.subplots_adjust(hspace=hspace)
plt.savefig(str(dir_out) + "/Fig06.png", dpi=300, bbox_inches="tight")
plt.show()
