#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
creates Fig 05
creates plot of irrigation water used for different water application schemes
"""


import os

import matplotlib.pyplot as plt
import xarray as xr

from analysis_functions.functions_correcting_time import correct_timedim
from analysis_functions.functions_plotting import plot_rotvar, rotated_pole
from analysis_functions.functions_reading_files import read_efiles

# In[]:experiment for analyzing

year = 2017
month = 6

exp_number_noirri = "067015"
exp_number_irri_prescribed = "067020"
exp_number_irri_flextime = "067019"
exp_number_irri_adapt = "067017"

# In[]: define plot directory


dir_working = os.getcwd()
# creates dir in parent directory
dir_out = os.path.join(os.pardir, "Figures")
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
print("Output directory is: ", dir_out)

# In[]: read data
# paths to the data
data_path = "../data"

# background map: irrifrac
remo_dir = str(exp_number_noirri) + "/irrifrac/"
remo_files = "e" + str(exp_number_noirri) + "e_c743_201706.nc"
remo_tfile = xr.open_dataset(str(data_path) + "/" + str(remo_dir) + str(remo_files))
irrifrac = remo_tfile.IRRIFRAC[0]


varlist = ["IRRWR"]
var_num_list = ["796"]
for var, var_num in zip(varlist, var_num_list):
    single_var_data_adapt = read_efiles(
        data_path, var, var_num, exp_number_irri_adapt, year, month
    )
    single_var_data_prescribed = read_efiles(
        data_path, var, var_num, exp_number_irri_prescribed, year, month
    )
    if var == varlist[0]:
        ds_var_irri_adapt = single_var_data_adapt
        ds_var_irri_prescribed = single_var_data_prescribed
    else:
        ds_var_irri_adapt = xr.merge([ds_var_irri_adapt, single_var_data_adapt])
        ds_var_irri_prescribed = xr.merge(
            [ds_var_irri_prescribed, single_var_data_prescribed]
        )
ds_var_irri_adapt = xr.merge([ds_var_irri_adapt, irrifrac])
ds_var_irri_prescribed = xr.merge([ds_var_irri_prescribed, irrifrac])


varlist_flextime = ["IRRWR", "IRRDUR"]
var_num_list = ["796", "795"]
for var, var_num in zip(varlist_flextime, var_num_list):
    single_var_data = read_efiles(
        data_path, var, var_num, exp_number_irri_flextime, year, month
    )
    if var == varlist[0]:
        ds_var_irri = single_var_data
    else:
        ds_var_irri = xr.merge([ds_var_irri, single_var_data])
ds_var_irri_flextime = xr.merge([ds_var_irri, irrifrac])


dsirr_adapt = correct_timedim(ds_var_irri_adapt)
dsirr_flextime = correct_timedim(ds_var_irri_flextime)
dsirr_prescribed = correct_timedim(ds_var_irri_prescribed)


# irrwr
irrwr_prescribed = dsirr_prescribed.IRRWR.groupby("time.hour")[18].sum("time") * 1000
irrwr_prescribed = irrwr_prescribed.where(irrwr_prescribed > 0)

irrwr_flextime = dsirr_flextime.IRRWR.groupby("time.hour")[23].sum("time") * 1000
irrwr_flextime = irrwr_flextime.where(irrwr_flextime > 0)
irrwr_adapt = (
    dsirr_adapt.IRRWR.where(dsirr_adapt.IRRWR != 0)
    .resample(time="1D")
    .max()
    .sum("time")
    * 1000
)
irrwr_adapt = irrwr_adapt.where(irrwr_adapt > 0)


# spatial plot

fig = plt.figure(figsize=(18, 4))
ax1 = fig.add_subplot(1, 3, 1, projection=rotated_pole)
cax1 = fig.add_axes([0.148, 0.01, 0.19, 0.04])
levels1 = [0, 50, 100, 150, 200, 250, 300, 350, 400]
ticks1 = levels1
# irrwr prescribed
plot_rotvar(
    fig,
    irrwr_prescribed,
    ax1,
    cax1,
    "[mm]",
    "irrigation water [mm] ",
    "viridis_r",
    levels1,
    ticks1,
    "max",
    "horizontal",
)

ax2 = fig.add_subplot(1, 3, 2, projection=rotated_pole)
cax2 = fig.add_axes([0.42, 0.01, 0.19, 0.04])
# irrwr flextime
plot_rotvar(
    fig,
    irrwr_flextime,
    ax2,
    cax2,
    "[mm]",
    "irrigation water [mm] ",
    "viridis_r",
    levels1,
    ticks1,
    "max",
    "horizontal",
)

ax3 = fig.add_subplot(1, 3, 3, projection=rotated_pole)
cax3 = fig.add_axes([0.69, 0.01, 0.19, 0.04])
# irrwr adapt
plot_rotvar(
    fig,
    irrwr_adapt,
    ax3,
    cax3,
    "[mm]",
    "irrigation water [mm] ",
    "viridis_r",
    levels1,
    ticks1,
    "max",
    "horizontal",
)
plt.savefig(str(dir_out) + "/Fig05.png", dpi=300, bbox_inches="tight")
plt.show()
