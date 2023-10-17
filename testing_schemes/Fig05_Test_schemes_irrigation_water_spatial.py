#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
creates Fig 05
creates plot of irrigation water used for different water application schemes
"""


import os

import matplotlib.pyplot as plt
import numpy as np
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

# calculate the differences between T1 - T2 and T1 - T3
diff_prescribed_flextime = irrwr_prescribed - irrwr_flextime
# diff_flextime_adative    = irrwr_flextime   - irrwr_adapt
diff_prescribed_adaptive = irrwr_prescribed - irrwr_adapt


# spatial plot

fig = plt.figure(figsize=(18, 9))
ax1 = fig.add_subplot(2, 3, 1, projection=rotated_pole)
cax1 = fig.add_axes([0.148, 0.51, 0.19, 0.03])
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
ax1.text(
    0.02,
    0.976,
    "T1",
    transform=ax1.transAxes,
    fontsize=16,
    fontweight="bold",
    va="top",
    bbox=dict(facecolor="white", edgecolor="black", pad=5.0),
    zorder=+6,
)

ax2 = fig.add_subplot(2, 3, 2, projection=rotated_pole)
cax2 = fig.add_axes([0.42, 0.51, 0.19, 0.03])
levels2 = [-40, -35, -30, -25, -20, -15, -10, -5, 5, 10, 15, 20, 25, 30, 35, 40]
ticks2 = [-40, -30, -20, -10, 10, 20, 30, 40]
# irrwr prescribed
plot_rotvar(
    fig,
    diff_prescribed_flextime,
    ax2,
    cax2,
    "[mm]",
    "irrigation water [mm] ",
    "RdBu_r",
    levels2,
    ticks2,
    "both",
    "horizontal",
)
ax2.text(
    0.02,
    0.976,
    "T1-T2",
    transform=ax2.transAxes,
    fontsize=16,
    fontweight="bold",
    va="top",
    bbox=dict(facecolor="white", edgecolor="black", pad=5.0),
    zorder=+6,
)


ax3 = fig.add_subplot(2, 3, 3, projection=rotated_pole)
cax3 = fig.add_axes([0.69, 0.51, 0.19, 0.03])
# irrwr adapt
plot_rotvar(
    fig,
    diff_prescribed_adaptive,
    ax3,
    cax3,
    "[mm]",
    "irrigation water [mm] ",
    "RdBu_r",
    levels2,
    ticks2,
    "both",
    "horizontal",
)
ax3.text(
    0.02,
    0.976,
    "T1-T3",
    transform=ax3.transAxes,
    fontsize=16,
    fontweight="bold",
    va="top",
    bbox=dict(facecolor="white", edgecolor="black", pad=5.0),
    zorder=+6,
)
# plt.savefig(str(dir_out) + "/Fig05_diff_revised.png", dpi=300, bbox_inches="tight")
# plt.show()

# distribution

ax4 = fig.add_subplot(2, 3, 5)
xr.plot.hist(
    irrwr_prescribed,
    ax=ax4,
    alpha=0.5,
    bins=np.arange(0, 500, 20),
    density=False,
    label="prescribed",
)
xr.plot.hist(
    irrwr_flextime,
    ax=ax4,
    alpha=0.5,
    bins=np.arange(0, 500, 20),
    density=False,
    label="flextime",
)
xr.plot.hist(
    irrwr_adapt,
    ax=ax4,
    alpha=0.5,
    bins=np.arange(0, 500, 20),
    density=False,
    label="adaptive",
)
plt.legend(["prescribed (T1)", "flextime (T2)", "adaptve (T3)"])
ax4.set_xlabel("irrigation water [mm]")
ax4.set_ylabel("number of grid cells")
# plt.tight_layout()
fig.subplots_adjust(
    left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.35
)

plt.savefig(str(dir_out) + "/Fig05.png", dpi=300, bbox_inches="tight")

# In[]: distribution of irrigation water

irrwr_prescribed_1d = irrwr_prescribed.values.reshape(
    len(irrwr_prescribed.rlon) * len(irrwr_prescribed.rlat)
)
irrwr_flextime_1d = irrwr_flextime.values.reshape(
    len(irrwr_flextime.rlon) * len(irrwr_flextime.rlat)
)
irrwr_adapt_1d = irrwr_adapt.values.reshape(
    len(irrwr_adapt.rlon) * len(irrwr_adapt.rlat)
)


fig, ax1 = plt.subplots()
# only kde plot
xr.plot.hist(
    irrwr_prescribed,
    ax=ax1,
    alpha=0.5,
    bins=np.arange(0, 500, 20),
    density=False,
    label="prescribed",
)
xr.plot.hist(
    irrwr_flextime,
    ax=ax1,
    alpha=0.5,
    bins=np.arange(0, 500, 20),
    density=False,
    label="flextime",
)
xr.plot.hist(
    irrwr_adapt,
    ax=ax1,
    alpha=0.5,
    bins=np.arange(0, 500, 20),
    density=False,
    label="adaptive",
)
ax1.set_xlabel("irrigation water [mm]")
ax1.set_ylabel("number of grid cells")
plt.legend()
# plt.savefig(str(dir_out) + "/hist_irrwr.png", dpi=300, bbox_inches="tight")
