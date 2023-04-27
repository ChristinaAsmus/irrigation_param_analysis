#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
creates Fig12
crates plot for irrigation effects as spatial plot, barplot (irri vs. noirri), diff plot on precipitation
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from analysis_functions.functions_calculations import (
    calculate_meandiff_precip,
    calculate_monthly_sum_precip,
)
from analysis_functions.functions_correcting_time import (
    correct_timedim,
    correct_timedim_mfiles,
)
from analysis_functions.functions_plotting import plot_rotvar, rotated_pole
from analysis_functions.functions_reading_files import read_efiles, read_mfiles

# In[]: select experiment, year and month

year = 2017
month = 0
exp_number_irri = "067016"
exp_number_noirri = "067015"


# In[]: read files


# paths on levante
levante_path_noirri = "/work/ch0636/g300099/Paper1/data/not_irrigated/sim/"
levante_path_irri_sim = "/work/ch0636/g300099/Paper1/data/irrigated/sim/"

# background map: irrifrac
remo_dir = str(exp_number_noirri) + "/irrifrac/"
remo_files = "e" + str(exp_number_noirri) + "e_c743_201706.nc"
remo_tfile = xr.open_dataset(levante_path_noirri + remo_dir + remo_files)
irrifrac = remo_tfile.IRRIFRAC[0]

varlist = ["APRL", "APRC"]
var_num_list = ["142", "143"]

for var, var_num in zip(varlist, var_num_list):
    single_var_data_irri = read_efiles(
        levante_path_irri_sim, var, var_num, exp_number_irri, year, month
    )
    single_var_data_noirri = read_efiles(
        levante_path_noirri, var, var_num, exp_number_noirri, year, month
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

# In[]: read mfiles
mds_irri = read_mfiles(levante_path_irri_sim, exp_number_irri, 2017, 0)
mds_noirri = read_mfiles(levante_path_noirri, exp_number_noirri, 2017, 0)

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

# In[]: precipitation

irrilimit = 0.0

precip_varlist = ["APRL", "APRC"]

# precipitation barplot for months
x = np.arange(
    dsirr_newtime.time.dt.month.values.min(),
    dsirr_newtime.time.dt.month.values.max() + 1,
)
x = np.arange(0, 12)
month_index = [
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sept",
    "Oct",
    "Nov",
    "Dec",
]

simlist = ["irri", "noirri"]

data_df = pd.DataFrame()
data_df["month"] = month_index
data_df.index = data_df.month
data_df.drop("month", axis=1, inplace=True)
for s, sim in enumerate(simlist):
    data_df[sim] = calculate_monthly_sum_precip(
        precip_varlist[0], precip_varlist[1], irrilimit, dsirr_newtime, dsnoirr_newtime
    )[s]


# use only fully irrigated months AMJ - spatial
dsirr_AMJ = dsirr_newtime.sel(time=dsirr_newtime.time.dt.month.isin([4, 5, 6]))
dsnoirr_AMJ = dsnoirr_newtime.sel(time=dsnoirr_newtime.time.dt.month.isin([4, 5, 6]))
precipdiff = calculate_meandiff_precip(
    dsirr_AMJ, dsnoirr_AMJ, precip_varlist[0], precip_varlist[1]
)

# plot
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
cax1 = fig.add_axes([0.135, 0.01, 0.19, 0.04])
levels = [-100, -80, -60, -40, -20, 20, 40, 60, 80, 100]
ticks = [-80, -40, 40, 80]
rotplot = plot_rotvar(
    fig,
    precipdiff,
    ax1,
    cax1,
    unit="[Δmm]",
    label="precipitation [Δmm]",
    cmap="RdBu_r",
    levels=levels,
    extend_scale="both",
    ticks=ticks,
    cbar_orient="horizontal",
)


ax2 = fig.add_subplot(1, 3, 2)
data_df.plot.bar(ax=ax2, rot=45, zorder=3)
ax2.grid(axis="y", linewidth=0.5)
ax2.legend(["irrigated", "not irrigated"])
ax2.set_ylabel("[mm]")

# monthly diff
diff_df = pd.DataFrame()
diff_df["month"] = month_index
diff_df.index = diff_df.month
diff_df.drop("month", axis=1, inplace=True)
diff_df["PRECIP"] = (
    calculate_monthly_sum_precip(
        precip_varlist[0], precip_varlist[1], irrilimit, dsirr_newtime, dsnoirr_newtime
    )
)[2]


ax3 = fig.add_subplot(1, 3, 3)
diff_df.plot.line(ax=ax3, rot=45, legend=False)
ax3.set_xticks(x, diff_df.index, minor=False)
ax3.set_xticks(ax3.get_xticks()[::2])
ax3.set_xlabel("month")
ax3.set_ylabel("[Δmm]")
ax3.set_ylim(-0.6, 1.0)
ax3.grid(True)
wspace = 0.35
fig.subplots_adjust(wspace=wspace)
plt.savefig(str(dir_out) + "/Fig12.png", dpi=300, bbox_inches="tight")
plt.show()
