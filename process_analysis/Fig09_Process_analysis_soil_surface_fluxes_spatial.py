#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
creates Fig 9
creates spatial plot of irrigation effects on sensible & latent heat flux and Bowen ratio
"""


import os
import sys

import matplotlib.pyplot as plt
import xarray as xr

sys.path.append(
    "/home/g/g300099/pyprograms/Paper1_for_publishing_final/plot_figures/functions/"
)
from functions_calculations import *
from functions_correcting_time import *
from functions_plotting import *
from functions_reading_files import *

# In[]: select experiment, year and month

year = 2017
month = 0
exp_number_irri = "067016"
exp_number_noirri = "067015"


# In[]: read files

# background map: irrifrac
remo_dir = (
    "/work/ch0636/g300099/SIMULATIONS/GAR11/remo_results/"
    + str(exp_number_noirri)
    + "/2017/var_series/IRRIFRAC/"
)
remo_files = "e" + str(exp_number_noirri) + "e_c743_201706.nc"
remo_tfile = xr.open_dataset(remo_dir + remo_files)
irrifrac = remo_tfile.IRRIFRAC[0]

varlist = ["AHFLIRR", "AHFSIRR"]
var_num_list = ["736", "734"]

for var, var_num in zip(varlist, var_num_list):
    single_var_data_irri = read_efiles(var, var_num, exp_number_irri, year, month)
    single_var_data_noirri = read_efiles(var, var_num, exp_number_noirri, year, month)
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

# In[]: read mfiles
mds_irri = read_mfiles(exp_number_irri, 2017, 0)
mds_noirri = read_mfiles(exp_number_noirri, 2017, 0)

mdsirri_extended = correct_timedim_mfiles(xr.merge([mds_irri, irrifrac]))
mdsnoirri_extended = correct_timedim_mfiles(xr.merge([mds_noirri, irrifrac]))

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


# In[]: fluxes

irrilimit = 0.0
mdsirri_AMJ = mdsirri_extended.sel(time=mdsirri_extended.time.dt.month.isin([4, 5, 6]))
mdsnoirri_AMJ = mdsnoirri_extended.sel(
    time=mdsnoirri_extended.time.dt.month.isin([4, 5, 6])
)
# diff fluxes
lhfl_diff = -(calculate_meandiff("AHFLIRR", mdsirri_AMJ, mdsnoirri_AMJ))
shfl_diff = -(calculate_meandiff("AHFSIRR", mdsirri_AMJ, mdsnoirri_AMJ))

# Bowen
bowen_irri_mean = (calculate_bowen_ratio(mdsirri_AMJ)).mean("time")
bowen_irri_mean_clean = bowen_irri_mean.where((bowen_irri_mean < 10 ^ 12))
bowen_noirri_mean = (calculate_bowen_ratio(mdsnoirri_AMJ)).mean("time")
bowen_noirri_mean_clean = bowen_noirri_mean.where((bowen_noirri_mean < 10 ^ 12))

bratio_diff = bowen_irri_mean_clean - bowen_noirri_mean_clean


fig = plt.figure(figsize=(18, 4))

params = {
    "legend.fontsize": 18,
    "axes.labelsize": 18,
    "axes.titlesize": 18,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
}
plt.rcParams.update(params)


ax1 = fig.add_subplot(1, 3, 1, projection=rotated_pole)
cax1 = fig.add_axes([0.145, 0.01, 0.19, 0.04])
levels = [-150, -125, -100, -75, -50, -25, 25, 50, 75, 100, 125, 150]
ticks = [-150, -100, -50, 50, 100, 150]
rotplot = plot_rotvar(
    fig,
    lhfl_diff,
    ax1,
    cax1,
    label="latent heat flux [ΔWm$^{-2}$]",
    unit="[ΔWm$^{-2}$]",
    cmap="RdBu_r",
    levels=levels,
    extend_scale="both",
    ticks=ticks,
    cbar_orient="horizontal",
)

ax2 = fig.add_subplot(1, 3, 2, projection=rotated_pole)
cax2 = fig.add_axes([0.42, 0.01, 0.19, 0.04])

rotplot = plot_rotvar(
    fig,
    shfl_diff,
    ax2,
    cax2,
    label="sensible heat flux [ΔWm$^{-2}$]",
    unit="[ΔWm$^{-2}$]",
    cmap="RdBu_r",
    levels=levels,
    extend_scale="both",
    ticks=ticks,
    cbar_orient="horizontal",
)


ax3 = fig.add_subplot(1, 3, 3, projection=rotated_pole)
cax3 = fig.add_axes([0.69, 0.01, 0.19, 0.04])
levels3 = [-1.0, -0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8, 1.0]
ticks = [
    -0.8,
    -0.4,
    0.4,
    0.8,
]
rotplot = plot_rotvar(
    fig,
    bratio_diff,
    ax3,
    cax3,
    label="Bowen ratio [ΔWm$^{-2}$]",
    unit="[ΔWm$^{-2}$]",
    cmap="RdBu_r",
    levels=levels3,
    extend_scale="both",
    ticks=ticks,
    cbar_orient="horizontal",
)
plt.savefig(str(dir_out) + "/Fig09.png", dpi=300, bbox_inches="tight")
plt.show()
