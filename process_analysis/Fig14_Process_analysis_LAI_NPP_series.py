#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
creates Fig 14
creates time series plots of LAI and NPP
"""


import os

import matplotlib.pyplot as plt
import xarray as xr

from analysis_functions.functions_correcting_time import correct_timedim_mfiles
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

# In[]:
# LAI

irrilimit = 0.0

var = "ALAI_PFI"


# monthly for irrigated fraction
lai_irri = mdsirri_extended.where(mdsirri_extended.IRRIFRAC > irrilimit)[var][
    :, 0, :, :
].mean(dim=["rlon", "rlat"], skipna=True)
lai_noirri = mdsnoirri_extended.where(mdsnoirri_extended.IRRIFRAC > irrilimit)[var][
    :, 0, :, :
].mean(dim=["rlon", "rlat"], skipna=True)
laimonthdiff = lai_irri - lai_noirri

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

npp_irri = (
    mdsirri_extended.where(mdsirri_extended.IRRIFRAC > irrilimit)[var][:, 0, :, :].mean(
        dim=["rlon", "rlat"], skipna=True
    )
) * (30 * 24)
npp_noirri = (
    mdsnoirri_extended.where(mdsnoirri_extended.IRRIFRAC > irrilimit)[var][
        :, 0, :, :
    ].mean(dim=["rlon", "rlat"], skipna=True)
) * (30 * 24)
nppmonthdiff = npp_irri - npp_noirri


month_list = ["April", "May", "June"]
mdsirri_AMJ = mdsirri_extended.sel(time=mdsirri_extended.time.dt.month.isin([4, 5, 6]))
mdsnoirri_AMJ = mdsnoirri_extended.sel(
    time=mdsnoirri_extended.time.dt.month.isin([4, 5, 6])
)
nppspatialdiff = ((mdsirri_AMJ[var][:, 0, :, :]) * (30 * 24)) - (
    (mdsnoirri_AMJ[var][:, 0, :, :]) * (30 * 24)
)
# line plot

fig = plt.figure(figsize=(16, 12))

params = {
    "legend.fontsize": 20,
    "axes.labelsize": 20,
    "axes.titlesize": 18,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
}
plt.rcParams.update(params)

ax1 = fig.add_subplot(2, 2, 1)
# laimonthdiff.plot.line(ax=ax1, marker=".")
# ax1.set_xticks(
#    ticks=laimonthdiff.time.values, labels=laimonthdiff.time.dt.strftime("%b").values
# )
lai_irri.plot.line(ax=ax1, label="irrigated", add_legend=False)
lai_noirri.plot.line(ax=ax1, label="not irrigated", add_legend=False)
ax1.set_xlabel(
    "months",
)
ax1.set_ylabel("LAI [m$^2$m$^{-2}$]")
ax1.grid(True)
ax1.set_title(" ")
ax1.legend(loc="upper right")
ax1.tick_params(axis="both")

ax2 = fig.add_subplot(2, 2, 2)
ax2.set_xticks(
    ticks=laimonthdiff.time.values, labels=laimonthdiff.time.dt.strftime("%b").values
)
laimonthdiff.plot.line(ax=ax2, marker=".")
ax2.set_ylim(-0.2, 0.2)
ax2.set_xlabel("months")
ax2.set_ylabel("Δ LAI [m$^2$m$^{-2}$]")
ax2.grid(True)
ax2.set_title(" ")
ax2.tick_params(axis="both")
plt.tight_layout()

ax3 = fig.add_subplot(2, 2, 3)
ax3.sharex(ax1)
ax3.set_xticks(
    ticks=nppmonthdiff.time.values, labels=nppmonthdiff.time.dt.strftime("%b").values
)
npp_irri.plot.line(ax=ax3, label="irrigated", add_legend=False)
npp_noirri.plot.line(ax=ax3, label="not irrigated", add_legend=False)
ax3.set_xlabel("months")
ax3.set_ylabel("NPP [gCm$^{-2}$month$^{-1}$]")
ax3.grid(True)
ax3.set_title(" ")
ax3.legend(loc="upper right")
ax3.tick_params(axis="both")

ax4 = fig.add_subplot(2, 2, 4)
ax4.set_xticks(
    ticks=nppmonthdiff.time.values, labels=nppmonthdiff.time.dt.strftime("%b").values
)
nppmonthdiff.plot.line(ax=ax4, marker=".")
ax4.set_xlabel("months")
ax4.set_ylabel("Δ NPP [gCm$^{-2}$month$^{-1}$]")
ax4.grid(True)
ax4.set_title(" ")
ax4.tick_params(axis="both")
plt.tight_layout()
plt.savefig(str(dir_out) + "/Fig14.png", dpi=300, bbox_inches="tight")
plt.show()
