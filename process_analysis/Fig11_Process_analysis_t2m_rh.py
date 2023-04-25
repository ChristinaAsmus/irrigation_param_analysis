#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
creates Fig 11
creates spatial plot (mean), boxplot for diurnal cycle and time series for annaul cycle for 2m temperature and 2m relative humdity

"""


import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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

varlist = ["DEW2", "TEMP2"]
var_num_list = ["168", "167"]

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


# In[]: atmosphere
# t2m
irrilimit = 0.0

# use only fully irrigated months AMJ
dsirr_AMJ = dsirr_newtime.sel(time=dsirr_newtime.time.dt.month.isin([4, 5, 6]))
dsnoirr_AMJ = dsnoirr_newtime.sel(time=dsnoirr_newtime.time.dt.month.isin([4, 5, 6]))

varlist = ["TEMP2", "T2MAX", "T2MIN"]
temp_str = ["T2Mean", "T2Max", "T2Min"]


# spatial plot for mean of AJM
tempspatialdiff = calculate_meandiff("TEMP2", dsirr_AMJ, dsnoirr_AMJ)[0, :, :]

# hourly boxplot
tempdiff = (dsirr_newtime.TEMP2.where(dsirr_AMJ["IRRIFRAC"] > irrilimit))[
    :, 0, :, :
] - (dsnoirr_newtime.TEMP2.where(dsnoirr_AMJ["IRRIFRAC"] > irrilimit))[:, 0, :, :]
tempdiffhour = (
    tempdiff.groupby("time.hour")
    .mean("time", skipna=True)
    .stack(space=["rlon", "rlat"])
    .reset_index("space")
    .drop(["rlon", "rlat", "lon", "lat"])
)
df_tempdiffhour = tempdiffhour.to_dataframe()


# Relative Humidity

mdsirri_AMJ = mdsirri_extended.sel(time=mdsirri_extended.time.dt.month.isin([4, 5, 6]))
mdsnoirri_AMJ = mdsnoirri_extended.sel(
    time=mdsnoirri_extended.time.dt.month.isin([4, 5, 6])
)

# for Magnus inverse
tdew_irri = mdsirri_AMJ.DEW2[:, 0, :, :] - 273.15
t_irri = mdsirri_AMJ.TEMP2[:, 0, :, :] - 273.15
tdew_noirri = mdsnoirri_AMJ.DEW2[:, 0, :, :] - 273.15
t_noirri = mdsnoirri_AMJ.TEMP2[:, 0, :, :] - 273.15

# spatial plot
rh_irri_AMJ = (calculate_rh(tdew_irri, t_irri)).mean("time")
rh_noirri_AMJ = (calculate_rh(tdew_noirri, t_noirri)).mean("time")
rhdiff_AMJ = rh_irri_AMJ - rh_noirri_AMJ

# diurnal / hourly
dsirri_AMJ = dsirr_newtime.sel(time=dsirr_newtime.time.dt.month.isin([4, 5, 6]))
dsnoirri_AMJ = dsnoirr_newtime.sel(time=dsnoirr_newtime.time.dt.month.isin([4, 5, 6]))


tdew_irri_series_AMJ = dsirri_AMJ.DEW2[:, 0, :, :] - 273.15
t_irri_series_AMJ = dsirri_AMJ.TEMP2[:, 0, :, :] - 273.15
tdew_noirri_series_AMJ = dsnoirri_AMJ.DEW2[:, 0, :, :] - 273.15
t_noirri_series_AMJ = dsnoirri_AMJ.TEMP2[:, 0, :, :] - 273.15

rh_irri_series_AMJ = calculate_rh(tdew_irri_series_AMJ, t_irri_series_AMJ)
rh_noirri_series_AMJ = calculate_rh(tdew_noirri_series_AMJ, t_noirri_series_AMJ)

irrilimit = 0.0

rhhourdiff = (
    (rh_irri_series_AMJ.where(mdsirri_extended["IRRIFRAC"] > irrilimit))
    .mean(["rlat", "rlon"], skipna=True)
    .groupby("time.hour")
    .mean(skipna=True)
) - (
    (rh_noirri_series_AMJ.where(mdsnoirri_extended["IRRIFRAC"] > irrilimit))
    .mean(["rlat", "rlon"], skipna=True)
    .groupby("time.hour")
    .mean(skipna=True)
)


rhdiff = (rh_irri_series_AMJ.where(mdsirri_extended["IRRIFRAC"] > irrilimit)) - (
    rh_noirri_series_AMJ.where(mdsnoirri_extended["IRRIFRAC"] > irrilimit)
)
rhdiffhour = (
    rhdiff.groupby("time.hour")
    .mean("time", skipna=True)
    .stack(space=["rlon", "rlat"])
    .reset_index("space")
    .drop(["rlon", "rlat", "lon", "lat"])
)
s_rhdiffhour = rhdiffhour.to_series()
df_rhdiffhour = s_rhdiffhour.to_frame("RH2")


# annual / monthly
tdew_irri_months = mdsirri_extended.DEW2[:, 0, :, :] - 273.15
t_irri_months = mdsirri_extended.TEMP2[:, 0, :, :] - 273.15
tdew_noirri_months = mdsnoirri_extended.DEW2[:, 0, :, :] - 273.15
t_noirri_months = mdsnoirri_extended.TEMP2[:, 0, :, :] - 273.15

rh_irri_months = calculate_rh(tdew_irri_months, t_irri_months)
rh_irri_months = (rh_irri_months.where(mdsirri_extended["IRRIFRAC"] > irrilimit)).mean(
    ["rlat", "rlon"], skipna=True
)
rh_noirri_months = calculate_rh(tdew_noirri_months, t_noirri_months)
rh_noirri_months = (
    rh_noirri_months.where(mdsirri_extended["IRRIFRAC"] > irrilimit)
).mean(["rlat", "rlon"], skipna=True)

rhdiff_months = rh_irri_months - rh_noirri_months

# plot rh and temperature

fig = plt.figure(figsize=(18, 10))

params = {
    "legend.fontsize": 18,
    "axes.labelsize": 18,
    "axes.titlesize": 18,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
}
plt.rcParams.update(params)

ax1 = fig.add_subplot(2, 3, 1, projection=rotated_pole)
cax1 = fig.add_axes([0.12, 0.55, 0.18, 0.02])
levels = [-2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2]
ticks = [-2, -1, 0, 1, 2]
rotplot = plot_rotvar(
    fig,
    tempspatialdiff,
    ax1,
    cax1,
    label="T2Mean [ΔK]",
    unit="[ΔK]",
    cmap="RdBu_r",
    levels=levels,
    extend_scale="both",
    ticks=ticks,
    cbar_orient="horizontal",
)
ax1.set_title(" ")

ax2 = fig.add_subplot(2, 3, 2)
box = df_tempdiffhour.boxplot(
    column="TEMP2", by="hour", ax=ax2, whis=[5, 95]
)  # , showfliers=False)
ax2.set_xticks(ax2.get_xticks()[::10])
ax2.set_ylabel("[ΔK]")
ax2.set_ylim(-4.5, 1.0)
ax2.set_title(" ")
plt.suptitle("")
ax2.grid(axis="x")

ax3 = fig.add_subplot(2, 3, 3)
xticklabels = [
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec",
]
irrilimit = 0.0
for var in varlist:
    tempmonthdiff = calculate_monthmeans_diff(
        var, irrilimit, mdsirri_extended, mdsnoirri_extended
    )[:, 0]
    tempmonthdiff.plot.line(ax=ax3, marker=".")
    ax3.set_xticks(
        ticks=tempmonthdiff.time.values,
        labels=tempmonthdiff.time.dt.strftime("%b").values,
    )
    ax3.set_xticklabels(xticklabels, rotation=45)
    ax3.set_xticks(ax3.get_xticks()[::2])
    ax3.set_xlabel("month")
    ax3.set_ylabel("[ΔK]")
    ax3.legend(temp_str, bbox_to_anchor=(1.1, 1.05))
    ax3.set_title(" ")
    ax3.grid(True)

ax4 = fig.add_subplot(2, 3, 4, projection=rotated_pole)
cax4 = fig.add_axes([0.12, 0.1, 0.18, 0.02])
levels4 = [-16, -14, -12, -10, -8, -6, -4, -2, 2, 4, 6, 8, 10, 12, 14, 16]
ticks = [-16, -12, -8, -4, 4, 8, 12, 16]
rotplot = plot_rotvar(
    fig,
    rhdiff_AMJ,
    ax4,
    cax4,
    unit="[Δ%]",
    label="relative humidity [Δ%]",
    cmap="RdBu_r",
    extend_scale="both",
    levels=levels4,
    ticks=ticks,
    cbar_orient="horizontal",
)
ax4.set_title(" ")


ax5 = fig.add_subplot(2, 3, 5)
box = df_rhdiffhour.boxplot(
    column="RH2", by="hour", ax=ax5, whis=[5, 95]
)  # , showfliers=False)
ax5.set_xticks(ax5.get_xticks()[::10])
ax5.set_ylabel("[Δ%]")
# ax5.set_ylim(-4.5,1.)
ax5.set_title(" ")
plt.suptitle("")  # that's what you're after
ax5.grid(axis="x")


xticklabels = [
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec",
]
ax6 = fig.add_subplot(2, 3, 6)
rhdiff_months.plot(ax=ax6, marker=".")
ax6.set_xticks(
    rhdiff_months.time.values, labels=rhdiff_months.time.dt.strftime("%b").values
)
ax6.set_xticklabels(xticklabels, rotation=45)
ax6.set_xticks(ax6.get_xticks()[::2])
ax6.set_xlabel("month")
ax6.set_ylabel("[Δ%]")
ax6.set_title(" ")
ax6.set_ylim(0, 3.7)
ax6.grid(True)

wspace = 0.35
hspace = 0.5
fig.subplots_adjust(wspace=wspace, hspace=hspace)
plt.savefig(str(dir_out) + "/Fig11.png", dpi=300, bbox_inches="tight")
plt.show()
