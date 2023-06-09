#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
creates Fig 04
creates subplot with increase of soil moisture development and hourly irrigation distribution of the different water application schemes
"""


import os
from datetime import datetime

import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr

# sys.path.append(
#    "/home/g/g300099/pyprograms/Paper1_for_publishing_final/plot_figures/analysis_functions/"
# )
from analysis_functions.functions_correcting_time import correct_timedim
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


varlist = ["IRRWR", "WSECHIRR", "WSMXIRR"]
var_num_list = ["796", "701", "728"]
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


varlist = ["IRRWR", "WSECHIRR", "WSMXIRR", "IRRDUR"]
var_num_list = ["796", "701", "728", "795"]
for var, var_num in zip(varlist, var_num_list):
    single_var_data = read_efiles(
        data_path, var, var_num, exp_number_irri_flextime, year, month
    )
    if var == varlist[0]:
        ds_var_irri = single_var_data
    else:
        ds_var_irri = xr.merge([ds_var_irri, single_var_data])
ds_var_irri_flextime = xr.merge([ds_var_irri, irrifrac])


varlist = ["WSECHIRR", "WSMXIRR"]
var_num_list = ["701", "728", "142", "143", "713", "700", "732"]
for var, var_num in zip(varlist, var_num_list):
    single_var_data = read_efiles(
        data_path, var, var_num, exp_number_noirri, year, month
    )
    if var == varlist[0]:
        ds_var_noirri = single_var_data
    else:
        ds_var_noirri = xr.merge([ds_var_noirri, single_var_data])
ds_var_noirri = xr.merge([ds_var_noirri, irrifrac])


dsirr_adapt = correct_timedim(ds_var_irri_adapt)
dsirr_flextime = correct_timedim(ds_var_irri_flextime)
dsirr_prescribed = correct_timedim(ds_var_irri_prescribed)
dsnoirr = correct_timedim(ds_var_noirri)


# In[]: for paper wsechirr time series 1. day

# gridcell to analyse
x = 63
y = 82

ws_adapt = (
    dsirr_adapt["WSECHIRR"][0:24, x, y] / dsirr_adapt["WSMXIRR"][0:24, x, y]
).to_series()
ws_flextime = (
    dsirr_flextime["WSECHIRR"][0:24, x, y] / dsirr_flextime["WSMXIRR"][0:24, x, y]
).to_series()
ws_prescribed = (
    dsirr_prescribed["WSECHIRR"][0:24, x, y] / dsirr_prescribed["WSMXIRR"][0:24, x, y]
).to_series()
ws_noirri = (
    dsnoirr["WSECHIRR"][0:24, x, y] / dsnoirr["WSMXIRR"][0:24, x, y]
).to_series()

ws_df = (pd.concat([ws_prescribed, ws_flextime, ws_adapt, ws_noirri], axis=1)).rename(
    columns={0: "prescribed", 1: "flextime", 2: "adaptive", 3: "not irrigated"}
)
ws_df.index = ws_df.index.strftime("%H:%M")
labels = ["prescribed", "flextime", "adative", "not irrigated"]


# In[]: for paper irrwr
# prescribed
irrwater_gridcell_prescribed = dsirr_prescribed.IRRWR[0:24, x, y].to_series() * 10**3
irrwater_gridcell_hour_prescribed = (
    irrwater_gridcell_prescribed.diff()
    .where(irrwater_gridcell_prescribed.diff() > 0)
    .dropna()
)
irrwater_gridcell_hour_prescribed.index = (
    irrwater_gridcell_hour_prescribed.index.strftime("%H:%M")
)


# flextime
# for flextime we have to calculate the amount per hour using irrdur because it is calculated only as sum over all time
irrwr_flextime = dsirr_flextime["IRRWR"][0:24, x, y]
irrdur_flextime = dsirr_flextime["IRRDUR"][0:24, x, y]

# calculate time
# irrstart is the first output timestep which shows irrwr (in the model irrigation starts with /:01)
irrstart = 8
irrend = irrdur_flextime.time[irrstart] + (
    pd.Timedelta(minutes=float((irrdur_flextime / (60)).max().values))
)

irrwr_df_minutes = pd.DataFrame(
    index=pd.date_range(
        datetime(2017, 6, 1, hour=irrstart, minute=0),
        periods=float((irrdur_flextime / (60)).max().values),
        freq="1min",
    )
)
irrwr_df_minutes["flextime"] = (irrwr_flextime.max().values * 1000) / (
    irrdur_flextime.max().values / 60
)
irrwr_df_hours = irrwr_df_minutes.resample("1H").sum()
irrwr_df_hours.index = irrwr_df_hours.index.strftime("%H:%M")
irrwater_gridcell_hour_flextime = irrwr_df_hours

# adaptive
irrwater_gridcell_adapt = dsirr_adapt.IRRWR[0:24, x, y].to_series() * 10**3
irrwater_gridcell_hour_adapt = (
    irrwater_gridcell_adapt.diff().where(irrwater_gridcell_adapt.diff() > 0).dropna()
)
irrwater_gridcell_hour_adapt.index = irrwater_gridcell_hour_adapt.index.strftime(
    "%H:%M"
)


irrwater_merge = pd.concat(
    [
        irrwater_gridcell_hour_prescribed,
        irrwater_gridcell_hour_flextime,
        irrwater_gridcell_hour_adapt,
    ],
    axis=1,
)


# plot all together

fig = plt.figure(figsize=(22, 7))

ax1 = fig.add_subplot(121)

ws_df.plot.line(ax=ax1)
ax1.legend(labels, fontsize=14, loc=(0.045, 0.17))
ax1.grid()
ax1.set_ylim(0.4, 1.05)
ax1.set_ylabel("fraction of max. water-holding capacity [-]", fontsize=16)
ax1.set_xlabel("01.06.2017", fontsize=16)
ax1.tick_params(axis="both", which="major", labelsize=14)
plt.xticks(rotation=45)  # Rotates X-Axis Ticks by 45-degrees
ax1.axhline(0.75, color="k", linestyle="--")
ax1.axhline(1.0, color="k", linestyle="--")
ax1.annotate("irrigation threshold", (0.045, 0.76), fontsize=14)
ax1.annotate("irrigation target", (0.045, 1.01), fontsize=14)


ax2 = fig.add_subplot(122)
irrwater_merge.plot.bar(ax=ax2, rot=45)
ax2.set_ylabel("irrigation water [mm]", fontsize=16)
ax2.tick_params(axis="both", which="major", labelsize=14)
ax2.set_xlabel("01.06.2017", fontsize=16)
ax2.legend(["prescribed", "flextime", "adaptive"], fontsize=14)
ax2.yaxis.grid(color="gray")
plt.savefig(str(dir_out) + "/Fig04.png", dpi=300, bbox_inches="tight")
plt.show()
