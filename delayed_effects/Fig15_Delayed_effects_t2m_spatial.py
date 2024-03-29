#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
creates Fig 15
creates spatial plot of temperature effects during the heatwave
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from analysis_functions.functions_correcting_time import correct_timedim
from analysis_functions.functions_plotting import (
    draw_rectangel,
    plot_rotvar,
    rotated_pole,
)
from analysis_functions.functions_reading_files import read_efiles

# In[]: important definitions

year = 2017
exp_number_irri = "067016"
exp_number_noirri = "067015"

# In[]: read files


# paths to the data
data_path = "../data"
# background map: irrifrac
remo_dir = str(exp_number_noirri) + "/irrifrac/"

remo_files = "e" + str(exp_number_noirri) + "e_c743_201706.nc"
remo_tfile = xr.open_dataset(str(data_path) + "/" + str(remo_dir) + str(remo_files))
irrifrac = remo_tfile.IRRIFRAC[0]

varlist = ["TEMP2", "T2MAX", "T2MIN"]
var_num_list = ["167", "201", "202"]
month = 8


for var, var_num in zip(varlist, var_num_list):
    single_var_data_irri = read_efiles(
        data_path, var, var_num, exp_number_irri, year, month
    )
    single_var_data_noirri = read_efiles(
        data_path, var, var_num, exp_number_noirri, year, month
    )
    if var == varlist[0]:
        ds_var_irri = single_var_data_irri
        ds_var_noirri = single_var_data_noirri
    else:
        ds_var_irri = xr.merge([ds_var_irri, single_var_data_irri])
        ds_var_noirri = xr.merge([ds_var_noirri, single_var_data_noirri])
ds_var_irri = xr.merge([ds_var_irri, irrifrac])
ds_var_noirri = xr.merge([ds_var_noirri, irrifrac])

dsirr_newtime = correct_timedim(ds_var_irri)
dsnoirr_newtime = correct_timedim(ds_var_noirri)


# In[]: define plot directory

dir_working = os.getcwd()
# creates dir in parent directory
dir_out = os.path.join(os.pardir, "Figures_corr")
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
print("Output directory is: ", dir_out)
# In[]: spatial plot for heat wave temperature


# 2m temp max (3-5Aug)
dsirr_newtime = correct_timedim(ds_var_irri)
dsirr_newtime = dsirr_newtime.sel(
    time=(
        dsirr_newtime.time.dt.month.isin([8])
        & dsirr_newtime.time.dt.day.isin(np.arange(3, 5, 1))
    )
)


dsnoirr_newtime = correct_timedim(ds_var_noirri)
dsnoirr_newtime = dsnoirr_newtime.sel(
    time=(
        dsnoirr_newtime.time.dt.month.isin([8])
        & dsnoirr_newtime.time.dt.day.isin(np.arange(3, 5, 1))
    )
)


tempspatialdiff = (
    (dsirr_newtime["TEMP2"][:, 0, :, :])
    .resample(time="D")
    .mean(skipna=True)
    .mean(dim="time", skipna=True)
    - 273.15
) - (
    (dsnoirr_newtime["TEMP2"][:, 0, :, :])
    .resample(time="D")
    .mean(skipna=True)
    .mean(dim="time", skipna=True)
    - 273.15
)
t2minspatialdiff = (
    (dsirr_newtime["T2MIN"][:, 0, :, :])
    .resample(time="D")
    .min(skipna=True)
    .mean(dim="time", skipna=True)
    - 273.15
) - (
    (dsnoirr_newtime["T2MIN"][:, 0, :, :])
    .resample(time="D")
    .min(skipna=True)
    .mean(dim="time", skipna=True)
    - 273.15
)

t2maxspatialdiff = (
    (dsirr_newtime["T2MAX"][:, 0, :, :])
    .resample(time="D")
    .max(skipna=True)
    .mean(dim="time", skipna=True)
    - 273.15
) - (
    (dsnoirr_newtime["T2MAX"][:, 0, :, :])
    .resample(time="D")
    .max(skipna=True)
    .mean(dim="time", skipna=True)
    - 273.15
)


levels = [-4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]
ticks = [-4, -3, -2, -1, 0, 1, 2, 3, 4]

fig = plt.figure(figsize=(18, 4))

params = {
    "legend.fontsize": 18,
    "axes.labelsize": 18,
    "axes.titlesize": 16,
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
}
plt.rcParams.update(params)


ax1 = fig.add_subplot(1, 3, 1, projection=rotated_pole)
cax1 = fig.add_axes([0.145, 0.01, 0.19, 0.04])
rotplot = plot_rotvar(
    fig,
    tempspatialdiff,
    ax1,
    cax1,
    unit="[K]",
    label="T2Mean [K]",
    cmap="RdBu_r",
    levels=levels,
    extend_scale="both",
    ticks=ticks,
    cbar_orient="horizontal",
)
ax1.set_title(" ")
region_IT = draw_rectangel(tempspatialdiff, 60, 110, 70, 50)
ax1.plot(
    region_IT[0],
    region_IT[1],
    transform=rotated_pole,
    color="k",
    zorder=4,
    linewidth=1.4,
    linestyle="--",
)


ax2 = fig.add_subplot(1, 3, 2, projection=rotated_pole)
cax2 = fig.add_axes([0.42, 0.01, 0.19, 0.04])
rotplot = plot_rotvar(
    fig,
    t2maxspatialdiff,
    ax2,
    cax2,
    unit="[K]",
    label="Δ T2Max [K]",
    cmap="RdBu_r",
    levels=levels,
    extend_scale="both",
    ticks=ticks,
    cbar_orient="horizontal",
)
ax2.set_title(" ")
region_IT = draw_rectangel(t2maxspatialdiff, 60, 110, 70, 50)
ax2.plot(
    region_IT[0],
    region_IT[1],
    transform=rotated_pole,
    color="k",
    zorder=4,
    linewidth=1.4,
    linestyle="--",
)


ax3 = fig.add_subplot(1, 3, 3, projection=rotated_pole)
cax3 = fig.add_axes([0.69, 0.01, 0.19, 0.04])
rotplot = plot_rotvar(
    fig,
    t2minspatialdiff,
    ax3,
    cax3,
    unit="[K]",
    label="Δ T2Min [K]",
    cmap="RdBu_r",
    levels=levels,
    extend_scale="both",
    ticks=ticks,
    cbar_orient="horizontal",
)
ax3.set_title(" ")
region_IT = draw_rectangel(t2maxspatialdiff, 60, 110, 70, 50)
ax3.plot(
    region_IT[0],
    region_IT[1],
    transform=rotated_pole,
    color="k",
    zorder=4,
    linewidth=1.4,
    linestyle="--",
)

plt.savefig(str(dir_out) + "/Fig15.png", dpi=300, bbox_inches="tight")
plt.show()
