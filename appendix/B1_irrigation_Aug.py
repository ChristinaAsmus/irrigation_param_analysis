#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
creates Fig B1
plotting of irrigation days in August using the irrigation mask
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from analysis_functions.functions_correcting_time import correct_timedim
from analysis_functions.functions_plotting import plot_rotvar, rotated_pole
from analysis_functions.functions_reading_files import read_efiles

# In[]: appendix irrimask for August

year = 2017
exp_number_irri = "067016"


varlist = ["IRRIMASK"]
var_num_list = ["797"]
month = 8


# paths on levante
levante_path_irri_sim = "/work/ch0636/g300099/Paper1/data/irrigated/sim/"

for var, var_num in zip(varlist, var_num_list):
    single_var_data = read_efiles(
        levante_path_irri_sim, var, var_num, exp_number_irri, year, month
    )
    if var == varlist[0]:
        ds_var_irri = single_var_data
    else:
        ds_var_irri = xr.merge([ds_var_irri, single_var_data])
dsirr_newtime = correct_timedim(ds_var_irri)

dir_working = os.getcwd()
# creates dir in parent directory
dir_out = os.path.join(os.pardir, "Figures")
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
print("Output directory is: ", dir_out)

month_name = "August"
# IRRIMASK is summed up over the month. Therefore we have to divide it by the number of active irrigation hours in the output (10h-1h)
nirrihours = 9
imasksum = dsirr_newtime.IRRIMASK.sum("time") / nirrihours
imaskvalues = imasksum.where(imasksum > 0, np.nan)


params = {
    "legend.fontsize": 14,
    "axes.labelsize": 16,
    "axes.titlesize": 18,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
}
plt.rcParams.update(params)

fig = plt.figure(figsize=(12, 7))
ax = fig.add_subplot(1, 1, 1, projection=rotated_pole)
cax = fig.add_axes([0.265, 0.05, 0.5, 0.04])
levels = np.arange(0, 31, 2)
ticks = levels
rotplot = plot_rotvar(
    fig,
    imaskvalues,
    ax,
    cax,
    "[n]",
    "number of irrigated days",
    "GnBu",
    levels,
    ticks,
    "max",
    "horizontal",
)
plt.savefig(str(dir_out) + "/FigB1.png", dpi=300, bbox_inches="tight")
plt.show()
