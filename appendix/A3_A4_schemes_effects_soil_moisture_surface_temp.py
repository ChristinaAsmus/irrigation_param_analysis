#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
creates Fig A3 and A4
plotting of irrigation effects on soil moisture and surface temperature using the different water application schemes
"""


import os

import matplotlib.pyplot as plt
import xarray as xr

from analysis_functions.functions_correcting_time import correct_timedim
from analysis_functions.functions_plotting import plot_rotvar, rotated_pole

# sys.path.append('/home/g/g300099/pyprograms/Paper1_for_publishing_final/plot_figures/functions/')
from analysis_functions.functions_reading_files import read_efiles

# In[]: appendix for schemes

# read data

year = 2017
month = 6

exp_number_noirri = "067015"
exp_number_irri_prescribed = "067020"
exp_number_irri_flextime = "067019"
exp_number_irri_adapt = "067017"


# paths to the data
data_path = "../data"

# background map: irrifrac
remo_dir = str(exp_number_noirri) + "/irrifrac/"
remo_files = "e" + str(exp_number_noirri) + "e_c743_201706.nc"
remo_tfile = xr.open_dataset(str(data_path) + "/" + str(remo_dir) + str(remo_files))
irrifrac = remo_tfile.IRRIFRAC[0]

varlist = ["WSECHIRR", "WSMXIRR", "TSECHIRR"]
var_num_list = ["701", "728", "732"]
for var, var_num in zip(varlist, var_num_list):
    single_var_data_adapt = read_efiles(
        data_path, var, var_num, exp_number_irri_adapt, year, month
    )
    single_var_data_prescribed = read_efiles(
        data_path, var, var_num, exp_number_irri_prescribed, year, month
    )
    single_var_data_flextime = read_efiles(
        data_path, var, var_num, exp_number_irri_flextime, year, month
    )
    single_var_data_noirri = read_efiles(
        data_path, var, var_num, exp_number_noirri, year, month
    )
    if var == varlist[0]:
        ds_var_irri_adapt = single_var_data_adapt
        ds_var_irri_prescribed = single_var_data_prescribed
        ds_var_irri_flextime = single_var_data_flextime
        ds_var_noirri = single_var_data_noirri
    else:
        ds_var_irri_adapt = xr.merge([ds_var_irri_adapt, single_var_data_adapt])
        ds_var_irri_prescribed = xr.merge(
            [ds_var_irri_prescribed, single_var_data_prescribed]
        )
        ds_var_irri_flextime = xr.merge(
            [ds_var_irri_flextime, single_var_data_flextime]
        )
        ds_var_noirri = xr.merge([ds_var_noirri, single_var_data_noirri])
ds_var_noirri = xr.merge([ds_var_noirri, irrifrac])
ds_var_irri_flextime = xr.merge([ds_var_irri_flextime, irrifrac])
ds_var_irri_adapt = xr.merge([ds_var_irri_adapt, irrifrac])
ds_var_irri_prescribed = xr.merge([ds_var_irri_prescribed, irrifrac])

dsirr_adapt = correct_timedim(ds_var_irri_adapt)
dsirr_flextime = correct_timedim(ds_var_irri_flextime)
dsirr_prescribed = correct_timedim(ds_var_irri_prescribed)
dsnoirr = correct_timedim(ds_var_noirri)

# In[]: define plot directory

dir_working = os.getcwd()
# creates dir in parent directory
dir_out = os.path.join(os.pardir, "Figures")
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
print("Output directory is: ", dir_out)

# In[]:
####################### diff between schemes and not irrigated ###############################
# wsechirr
var_ws = "WSECHIRR"
var_ws_prescribed_diff = dsirr_prescribed[var_ws].mean("time") - dsnoirr[var_ws].mean(
    "time"
)
var_ws_flextime_diff = dsirr_flextime[var_ws].mean("time") - dsnoirr[var_ws].mean(
    "time"
)
var_ws_adapt_diff = dsirr_adapt[var_ws].mean("time") - dsnoirr[var_ws].mean("time")

# tsechirr
var_ts = "TSECHIRR"
var_ts_prescribed_diff = dsirr_prescribed[var_ts].mean("time") - dsnoirr[var_ts].mean(
    "time"
)
var_ts_flextime_diff = dsirr_flextime[var_ts].mean("time") - dsnoirr[var_ts].mean(
    "time"
)
var_ts_adapt_diff = dsirr_adapt[var_ts].mean("time") - dsnoirr[var_ts].mean("time")

# subplot
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
# wsechirr prescribed
ax1 = fig.add_subplot(2, 3, 1, projection=rotated_pole)
cax1 = fig.add_axes([0.127, 0.52, 0.22, 0.02])
levels1 = [
    -0.6,
    -0.5,
    -0.4,
    -0.3,
    -0.2,
    -0.1,
    0.1,
    0.2,
    0.3,
    0.4,
    0.5,
    0.6,
]
ticks1 = [-0.6, -0.4, -0.2, 0.2, 0.4, 0.6]
plot_rotvar(
    fig,
    var_ws_prescribed_diff,
    ax1,
    cax1,
    "[Δm]",
    "soil moisture [Δm] ",
    "RdBu_r",
    levels1,
    ticks1,
    "both",
    cbar_orient="horizontal",
)

# wsechirr flextime
ax2 = fig.add_subplot(2, 3, 2, projection=rotated_pole)
cax2 = fig.add_axes([0.4, 0.52, 0.22, 0.02])
plot_rotvar(
    fig,
    var_ws_flextime_diff,
    ax2,
    cax2,
    "[Δm]",
    "soil moisture [Δm] ",
    "RdBu_r",
    levels1,
    ticks1,
    "both",
    cbar_orient="horizontal",
)

# wsechirr adapt
ax3 = fig.add_subplot(2, 3, 3, projection=rotated_pole)
cax3 = fig.add_axes([0.677, 0.52, 0.215, 0.02])
plot_rotvar(
    fig,
    var_ws_adapt_diff,
    ax3,
    cax3,
    "[Δm]",
    "soil moisture [Δm] ",
    "RdBu_r",
    levels1,
    ticks1,
    "both",
    cbar_orient="horizontal",
)

# tsechirr
# tsechirr prescribed
ax4 = fig.add_subplot(2, 3, 4, projection=rotated_pole)
cax4 = fig.add_axes([0.127, 0.11, 0.22, 0.02])
levels4 = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]
ticks4 = levels4
plot_rotvar(
    fig,
    var_ts_prescribed_diff,
    ax4,
    cax4,
    "[ΔK]",
    "soil temperature [ΔK] ",
    "RdBu_r",
    levels4,
    ticks4,
    "both",
    cbar_orient="horizontal",
)

# tsechirr flextime
ax5 = fig.add_subplot(2, 3, 5, projection=rotated_pole)
cax5 = fig.add_axes([0.4, 0.11, 0.22, 0.02])
plot_rotvar(
    fig,
    var_ts_flextime_diff,
    ax5,
    cax5,
    "[ΔK]",
    "soil temperature [ΔK] ",
    "RdBu_r",
    levels4,
    ticks4,
    "both",
    cbar_orient="horizontal",
)

# tsechirr adapt
ax6 = fig.add_subplot(2, 3, 6, projection=rotated_pole)
cax6 = fig.add_axes([0.678, 0.11, 0.22, 0.02])
plot_rotvar(
    fig,
    var_ts_adapt_diff,
    ax6,
    cax6,
    "[ΔK]",
    "soil temperature [ΔK] ",
    "RdBu_r",
    levels4,
    ticks4,
    "both",
    cbar_orient="horizontal",
)
plt.savefig(str(dir_out) + "FigA3.png", dpi=300, bbox_inches="tight")
plt.show()

############################## diff between schemes #######################
# difference between schemes
var_ws_flextime_prescribed = dsirr_flextime[var_ws].mean("time") - dsirr_prescribed[
    var_ws
].mean("time")
var_ws_adapt_prescribed = dsirr_adapt[var_ws].mean("time") - dsirr_prescribed[
    var_ws
].mean("time")
# differences between schemes
var_ts_flextime_prescribed = dsirr_flextime[var_ts].mean("time") - dsirr_prescribed[
    var_ts
].mean("time")
var_ts_adapt_prescribed = dsirr_adapt[var_ts].mean("time") - dsirr_prescribed[
    var_ts
].mean("time")


fig = plt.figure(figsize=(12, 12))
# #wsechirr prescribed
# ax1 = fig.add_subplot(2, 2, 1, projection=rotated_pole)
# cax1 = fig.add_axes([0.148, 0.01, 0.19, 0.04])
# levels1=[ -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
# ticks1=[ -0.6,-0.4, -0.2,  0.2,  0.4, 0.6]
# plot_rotvar(fig, var_ws_prescribed_diff, ax1, cax1, '[Δm]', 'soil moisture [Δm] ','RdBu_r',\
#             levels1,ticks1, 'both' , cbar_orient='horizontal')
# wsechirr
# wsechirr flextime-prescribed
ax2 = fig.add_subplot(2, 2, 1, projection=rotated_pole)
cax2 = fig.add_axes([0.12, 0.51, 0.35, 0.02])
levels2 = [-0.1, -0.08, -0.06, -0.04, -0.02, 0.02, 0.04, 0.06, 0.08, 0.1]
ticks2 = [-0.08, -0.04, 0.04, 0.08]
plot_rotvar(
    fig,
    var_ws_flextime_prescribed,
    ax2,
    cax2,
    "[Δm]",
    "soil moisture [Δm] ",
    "RdBu_r",
    levels2,
    ticks2,
    "both",
    cbar_orient="horizontal",
)

# wsechirr adapt-prescribed
ax3 = fig.add_subplot(2, 2, 2, projection=rotated_pole)
cax3 = fig.add_axes([0.55, 0.51, 0.35, 0.02])
plot_rotvar(
    fig,
    var_ws_adapt_prescribed,
    ax3,
    cax3,
    "[Δm]",
    "soil moisture [Δm] ",
    "RdBu_r",
    levels2,
    ticks2,
    "both",
    cbar_orient="horizontal",
)

# tsechirr
# #tsechirr adapt
# ax1 = fig.add_subplot(1, 3, 1, projection=rotated_pole)
# cax1 = fig.add_axes([0.148, 0.01, 0.19, 0.04])
# levels1=[-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]
# ticks1=levels1
# plot_rotvar(fig, var_ts_prescribed_diff, ax1, cax1, '[ΔK]', 'soil temperature [ΔK] ','RdBu_r',\
#             levels1, ticks1,'both' , cbar_orient='horizontal')

# tsechirr flextime-prescribed
ax5 = fig.add_subplot(2, 2, 3, projection=rotated_pole)
cax5 = fig.add_axes([0.12, 0.1, 0.35, 0.02])
levels5 = [-0.5, -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.5]
ticks5 = [-0.4, -0.2, 0.2, 0.4]
plot_rotvar(
    fig,
    var_ts_flextime_prescribed,
    ax5,
    cax5,
    "[ΔK]",
    "soil temperature [ΔK] ",
    "RdBu_r",
    levels5,
    ticks5,
    "both",
    cbar_orient="horizontal",
)

# tsechirr adapt-prescribed
ax6 = fig.add_subplot(2, 2, 4, projection=rotated_pole)
cax6 = fig.add_axes([0.55, 0.1, 0.35, 0.02])
plot_rotvar(
    fig,
    var_ts_adapt_prescribed,
    ax6,
    cax6,
    "[ΔK]",
    "soil temperature [ΔK] ",
    "RdBu_r",
    levels5,
    ticks5,
    "both",
    cbar_orient="horizontal",
)
plt.savefig(str(dir_out) + "/FigA4.png", dpi=300, bbox_inches="tight")
plt.show()
