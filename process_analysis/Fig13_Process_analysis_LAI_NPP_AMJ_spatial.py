#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
creates Fig 13
creates spatial plot for irrigation effects on LAI and NPP
Be careful with the unit of NPP, we use gC( of freshmatter)/m2month1
"""


import os

import matplotlib.pyplot as plt
import xarray as xr

from analysis_functions.functions_correcting_time import correct_timedim_mfiles
from analysis_functions.functions_plotting import plot_rotvar, rotated_pole
from analysis_functions.functions_reading_files import read_mfiles

# In[]: select experiment, year and month

year = 2017
month = 0
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

# In[]: read mfiles
mds_irri = read_mfiles(data_path, exp_number_irri, 2017, 0)
mds_noirri = read_mfiles(data_path, exp_number_noirri, 2017, 0)

mdsirri_extended = correct_timedim_mfiles(xr.merge([mds_irri, irrifrac]))
mdsnoirri_extended = correct_timedim_mfiles(xr.merge([mds_noirri, irrifrac]))

# In[]: define plot directory


dir_working = os.getcwd()
# creates dir in parent directory
dir_out = os.path.join(os.pardir, "Figures_corr")
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
print("Output directory is: ", dir_out)

# In[]: LAI and NPP
# LAI

irrilimit = 0.0

var = "ALAI_PFI"


# monthly for irrigated fraction
data_irri = mdsirri_extended.where(mdsirri_extended.IRRIFRAC > irrilimit)[var][
    :, 0, :, :
].mean(dim=["rlon", "rlat"], skipna=True)
data_noirri = mdsnoirri_extended.where(mdsnoirri_extended.IRRIFRAC > irrilimit)[var][
    :, 0, :, :
].mean(dim=["rlon", "rlat"], skipna=True)
laimonthdiff = data_irri - data_noirri

#
month_list = ["April", "May", "June"]
mdsirri_AMJ = mdsirri_extended.sel(time=mdsirri_extended.time.dt.month.isin([4, 5, 6]))
mdsnoirri_AMJ = mdsnoirri_extended.sel(
    time=mdsnoirri_extended.time.dt.month.isin([4, 5, 6])
)
laispatialdiff = mdsirri_AMJ.ALAI_PFI[:, 0, :, :] - mdsnoirri_AMJ.ALAI_PFI[:, 0, :, :]

# NPP
irrilimit = 0.0
var = "ANPP_ACI"

# monthly for irrigated fraction

data_irri = (
    mdsirri_extended.where(mdsirri_extended.IRRIFRAC > irrilimit)[var][:, 0, :, :].mean(
        dim=["rlon", "rlat"], skipna=True
    )
) * (30 * 24)
data_noirri = (
    mdsnoirri_extended.where(mdsnoirri_extended.IRRIFRAC > irrilimit)[var][
        :, 0, :, :
    ].mean(dim=["rlon", "rlat"], skipna=True)
) * (30 * 24)
nppmonthdiff = data_irri - data_noirri

month_list = ["April", "May", "June"]
mdsirri_AMJ = mdsirri_extended.sel(time=mdsirri_extended.time.dt.month.isin([4, 5, 6]))
mdsnoirri_AMJ = mdsnoirri_extended.sel(
    time=mdsnoirri_extended.time.dt.month.isin([4, 5, 6])
)
nppspatialdiff = ((mdsirri_AMJ[var][:, 0, :, :]) * (30 * 24)) - (
    (mdsnoirri_AMJ[var][:, 0, :, :]) * (30 * 24)
)

# spatial plot
fig = plt.figure(figsize=(14, 10))

params = {
    "legend.fontsize": 18,
    "axes.labelsize": 18,
    "axes.titlesize": 18,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
}
plt.rcParams.update(params)

cax1 = fig.add_axes([0.105, 0.53, 0.27, 0.02])
# cax2 = fig.add_axes([0.425, 0.54, 0.175, 0.02])
# cax3 = fig.add_axes([0.7, 0.54, 0.175, 0.02])
cbaraxes_list = [cax1, cax1, cax1]

for i in range(len(laispatialdiff)):
    ax1 = fig.add_subplot(2, 3, 1 + i, projection=rotated_pole)
    levels = [-1, -0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8, 1]
    ticks = [
        -0.8,
        -0.4,
        0,
        0.4,
        0.8,
    ]
    rotplot = plot_rotvar(
        fig,
        laispatialdiff[i],
        ax1,
        cbaraxes_list[i],
        label="Δ LAI [m$^2$m$^{-2}$]",
        unit="[-]",
        cmap="RdBu_r",
        levels=levels,
        extend_scale="both",
        ticks=ticks,
        cbar_orient="horizontal",
    )
    ax1.set_title(month_list[i])

# NPP
# diff spatial
cax4 = fig.add_axes([0.105, 0.11, 0.27, 0.02])
# cax5 = fig.add_axes([0.425, 0.07, 0.175, 0.02])
# cax6 = fig.add_axes([0.7, 0.07, 0.175, 0.02])
cbaraxes_list = [cax4, cax4, cax4]

for i in range(len(nppspatialdiff)):
    ax4 = fig.add_subplot(2, 3, 3 + (1 + i), projection=rotated_pole)
    levels = [
        -1200,
        -1100,
        -1000,
        -900,
        -800,
        -700,
        -600,
        -500,
        -400,
        -300,
        -200,
        -100,
        100,
        200,
        300,
        400,
        500,
        600,
        700,
        800,
        900,
        1000,
        1100,
        1200,
    ]
    ticks = [-1200, -400, 0, 400, 1200]
    rotplot = plot_rotvar(
        fig,
        nppspatialdiff[i],
        ax4,
        cbaraxes_list[i],
        label="Δ NPP [gCm$^{-2}$month$^{-1}$]",
        unit="[gCm$^{-2}$month$^{-1}$]",
        cmap="RdBu_r",
        levels=levels,
        extend_scale="both",
        ticks=ticks,
        cbar_orient="horizontal",
    )
    ax4.set_title(" ")
hspace = 0.2
fig.subplots_adjust(hspace=hspace)
plt.savefig(str(dir_out) + "/Fig13.png", dpi=300, bbox_inches="tight")
plt.show()
