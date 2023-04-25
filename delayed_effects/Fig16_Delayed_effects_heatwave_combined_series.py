#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
creates Fig 16
creates combined time series of temperature, soil moisture, precipitation evpotranspiration development during the heat wave
"""
import os
import sys

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import xarray as xr

sys.path.append(
    "/home/g/g300099/pyprograms/Paper1_for_publishing_final/plot_figures/functions/"
)
from functions_calculations import *
from functions_correcting_time import *
from functions_plotting import *
from functions_reading_files import *

# In[]: important definitions

year = 2017
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

varlist = [
    "TEMP2",
    "T2MAX",
    "T2MIN",
    "WSECHIRR",
    "WSMXIRR",
    "APRL",
    "APRC",
    "EVAPIRR",
    "EVAP",
    "EBSOIL",
    "ETRANS",
    "ALAI_PFI",
    "ESKIN",
]
var_num_list = [
    "167",
    "201",
    "202",
    "701",
    "728",
    "142",
    "143",
    "739",
    "182",
    "013",
    "012",
    "746",
    "302",
]
month = 8


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


# In[]: define plot directory


dir_working = os.getcwd()
# creates dir in parent directory
dir_out = os.path.join(os.pardir, "Figures")
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
print("Output directory is: ", dir_out)
# In[]: all as series


# 2m temp max (1-7Aug, hottest days were 3.-5. Aug)
dsirr_newtime = correct_timedim(ds_var_irri)
dsirr_newtime = dsirr_newtime.sel(
    time=(
        dsirr_newtime.time.dt.month.isin([8])
        & dsirr_newtime.time.dt.day.isin(np.arange(1, 8, 1))
    )
)


dsnoirr_newtime = correct_timedim(ds_var_noirri)
dsnoirr_newtime = dsnoirr_newtime.sel(
    time=(
        dsnoirr_newtime.time.dt.month.isin([8])
        & dsnoirr_newtime.time.dt.day.isin(np.arange(1, 8, 1))
    )
)

# time series for regions

rlat_list = [(50, 70)]
rlon_list = [(60, 110)]
title_list = ["IT"]

fig = plt.figure(figsize=(15, 19))
params = {
    "legend.fontsize": 20,
    "axes.labelsize": 20,
    "axes.titlesize": 20,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
}
plt.rcParams.update(params)

# temperature
for i in range(len(rlat_list)):
    temp2_noirr = (
        dsnoirr_newtime.isel(
            rlat=slice(rlat_list[i][0], rlat_list[i][1]),
            rlon=slice(rlon_list[i][0], rlon_list[i][1]),
        )["TEMP2"][:, 0, :, :]
    ).resample(time="D").mean(skipna=True).mean(
        dim=["rlon", "rlat"], skipna=True
    ) - 273.15
    temp2_irr = (
        dsirr_newtime.isel(
            rlat=slice(rlat_list[i][0], rlat_list[i][1]),
            rlon=slice(rlon_list[i][0], rlon_list[i][1]),
        )["TEMP2"][:, 0, :, :]
    ).resample(time="D").mean(skipna=True).mean(
        dim=["rlon", "rlat"], skipna=True
    ) - 273.15
    # temp2_series=temp2.to_series()
    t2max_noirr = (
        dsnoirr_newtime.isel(
            rlat=slice(rlat_list[i][0], rlat_list[i][1]),
            rlon=slice(rlon_list[i][0], rlon_list[i][1]),
        )["T2MAX"][:, 0, :, :]
    ).resample(time="D").max(skipna=True).mean(
        dim=["rlon", "rlat"], skipna=True
    ) - 273.15
    t2max_irr = (
        dsirr_newtime.isel(
            rlat=slice(rlat_list[i][0], rlat_list[i][1]),
            rlon=slice(rlon_list[i][0], rlon_list[i][1]),
        )["T2MAX"][:, 0, :, :]
    ).resample(time="D").max(skipna=True).mean(
        dim=["rlon", "rlat"], skipna=True
    ) - 273.15
    # t2max_series=t2max.to_series()
    t2min_noirr = (
        dsnoirr_newtime.isel(
            rlat=slice(rlat_list[i][0], rlat_list[i][1]),
            rlon=slice(rlon_list[i][0], rlon_list[i][1]),
        )["T2MIN"][:, 0, :, :]
    ).resample(time="D").min(skipna=True).mean(
        dim=["rlon", "rlat"], skipna=True
    ) - 273.15
    t2min_irr = (
        dsirr_newtime.isel(
            rlat=slice(rlat_list[i][0], rlat_list[i][1]),
            rlon=slice(rlon_list[i][0], rlon_list[i][1]),
        )["T2MIN"][:, 0, :, :]
    ).resample(time="D").min(skipna=True).mean(
        dim=["rlon", "rlat"], skipna=True
    ) - 273.15
    # t2min_series=t2min.to_series()

    ax1 = fig.add_subplot(4, 1, i + 1)
    plt.rcParams.update(params)

    temp2_irr.plot.line(ax=ax1, color="blue", linestyle="-", label="T2Mean irrigated")
    temp2_noirr.plot.line(
        ax=ax1, color="blue", linestyle="--", label="T2Mean not irrigated"
    )
    t2max_irr.plot.line(ax=ax1, color="orange", linestyle="-", label="T2Max irrigated")
    t2max_noirr.plot.line(
        ax=ax1, color="orange", linestyle="--", label="T2Max not irrigated"
    )
    t2min_irr.plot.line(ax=ax1, color="green", linestyle="-", label="T2Min irrigated")
    t2min_noirr.plot.line(
        ax=ax1, color="green", linestyle="--", label="T2Min not irrigated"
    )
    ax1.set_ylim(20, 38)
    ax1.set_ylabel("temperature \n [°C]")

    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%d.%m."))
    ax1.xaxis.set_major_locator(mdates.MonthLocator(bymonthday=[1, 3, 5, 7]))

    ax1.set_xlabel("2017")
    ax1.set_title(title_list[i])
    ax1.grid()
ax1.legend(bbox_to_anchor=(1.0, 1.05))

# soil moisture


for i in range(len(rlat_list)):
    relws_noirr = dsnoirr_newtime["WSECHIRR"] / dsnoirr_newtime["WSMXIRR"]
    ws_noirr = (
        (
            relws_noirr.isel(
                rlat=slice(rlat_list[i][0], rlat_list[i][1]),
                rlon=slice(rlon_list[i][0], rlon_list[i][1]),
            )
        )
        .resample(time="D")
        .mean(skipna=True)
        .mean(dim=["rlon", "rlat"], skipna=True)
    )
    # ws_noirr_series=ws_noirr.to_series()

    relws_irr = dsirr_newtime["WSECHIRR"] / dsirr_newtime["WSMXIRR"]
    ws_irr = (
        (
            relws_irr.isel(
                rlat=slice(rlat_list[i][0], rlat_list[i][1]),
                rlon=slice(rlon_list[i][0], rlon_list[i][1]),
            )
        )
        .resample(time="D")
        .mean(skipna=True)
        .mean(dim=["rlon", "rlat"], skipna=True)
    )
# ws_irr_series=ws_irr.to_series()

ax2 = fig.add_subplot(4, 1, 2, sharex=ax1)
plt.rcParams.update(params)
ws_irr.plot.line(ax=ax2, color="blue", linestyle="-", label="soil moisture irrigated")
ws_noirr.plot.line(
    ax=ax2, color="blue", linestyle="--", label="soil moisture not irrigated"
)
ax2.set_ylim(0, 1)
ax2.set_ylabel("rel. soil moisture \n [-]")
ax2.grid()
ax2.legend(bbox_to_anchor=(1.0, 1.05))

# precipitation
ax3 = fig.add_subplot(4, 1, 3)
plt.rcParams.update(params)


def calculate_sum_precip(data, var1, var2, i):
    precip = (
        data.isel(
            rlat=slice(rlat_list[i][0], rlat_list[i][1]),
            rlon=slice(rlon_list[i][0], rlon_list[i][1]),
        )[var1]
    ) + (
        data.isel(
            rlat=slice(rlat_list[i][0], rlat_list[i][1]),
            rlon=slice(rlon_list[i][0], rlon_list[i][1]),
        )[var2]
    )
    precip_mean = (
        precip.resample(time="D").sum().mean(dim=["rlon", "rlat"], skipna=True)
    )
    return precip_mean


precip_noirr_series = calculate_sum_precip(
    dsnoirr_newtime, "APRL", "APRC", i
).to_series()
precip_irr_series = calculate_sum_precip(dsirr_newtime, "APRL", "APRC", i).to_series()

df_precip = pd.DataFrame()
df_precip["precip_noirr"] = precip_noirr_series
df_precip["precip_irr"] = precip_irr_series
df_precip.index = df_precip.index.strftime("%d.%m.")
labels = ["precipitation irrigated", "precipitation not irrigated"]


df_precip[["precip_irr", "precip_noirr"]].plot.bar(
    ax=ax3, legend=True, width=0.6, zorder=1
)
# plt.legend(labels,loc='upper right',  bbox_to_anchor=(1.0, 1.0))
ax3.set_ylim(0.0, 6)
ax3.set_ylabel("precipitation \n [mm]")
ax3.set_xlabel(" ")
# tick_loc=mtick.FixedLocator([0,1,2,3,4,5,6])
# ax3.xaxis.set_major_locator(tick_loc)

for tick in ax3.get_xticklabels():
    #    tick.set_rotation(45)
    tick.set_visible(False)
# ax3.set_xlabel('2017')
ax3.yaxis.grid()  # horizontal lines
ax3.legend(labels, bbox_to_anchor=(1.42, 1.05))
#

# evapotranspiration

colors = {
    "evapotranspiration": "black",
    "evaporation of soil": "brown",
    "evaporation of skin reservoir": "blue",
    "transpiration of vegetation": "green",
}


irrilimit = 0.0  # maybe I have to take this out


# calculate monthly sum
def calculate_daily_sum_evapos(var, i, irrilimit):
    dsirr_month = (
        (
            dsirr_newtime.isel(
                rlat=slice(rlat_list[i][0], rlat_list[i][1]),
                rlon=slice(rlon_list[i][0], rlon_list[i][1]),
            )[var]
        )
        .groupby("time.day")
        .sum()
    ).mean(dim=["rlon", "rlat"], skipna=True)

    dsnoirr_month = (
        (
            dsnoirr_newtime.isel(
                rlat=slice(rlat_list[i][0], rlat_list[i][1]),
                rlon=slice(rlon_list[i][0], rlon_list[i][1]),
            )[var]
        )
        .groupby("time.day")
        .sum()
    ).mean(dim=["rlon", "rlat"], skipna=True)
    diff_month = dsirr_month - dsnoirr_month
    return dsirr_month, dsnoirr_month, diff_month


datelist = pd.to_datetime(
    [
        "2017-08-01",
        "2017-08-02",
        "2017-08-03",
        "2017-08-04",
        "2017-08-05",
        "2017-08-06",
        "2017-08-07",
    ],
    format="%Y-%m-%d",
)

for i in range(len(rlat_list)):
    irri_df = pd.DataFrame()
    noirri_df = pd.DataFrame()
    varlist = ["EVAPIRR", "EVAP", "EBSOIL", "ETRANS", "ESKIN"]
    for var in varlist:
        irri_df[var] = -(calculate_daily_sum_evapos(var, i, irrilimit))[0]
        irri_df["date"] = datelist.strftime("%d.%m.")
        irri_df.index = irri_df.date
        noirri_df[var] = -(calculate_daily_sum_evapos(var, i, irrilimit))[1]
        noirri_df["date"] = datelist.strftime("%d.%m.")
        noirri_df.index = noirri_df.date
    irri_df["EVAP_all"] = irri_df["EVAP"] + irri_df["EVAPIRR"]
    # irri_df.drop(['EVAP','EVAPIRR','date'], axis=1, inplace=True)
    noirri_df["EVAP_all"] = noirri_df["EVAP"] + noirri_df["EVAPIRR"]
    noirri_df.drop(["EVAP", "EVAPIRR", "date"], axis=1, inplace=True)

    b_1 = irri_df.EBSOIL.values
    a_1 = irri_df.EVAP_all.values
    c_1 = irri_df.ETRANS.values
    d_1 = irri_df.ESKIN.values

    b_2 = noirri_df.EBSOIL.values
    a_2 = noirri_df.EVAP_all.values
    c_2 = noirri_df.ETRANS.values
    d_2 = noirri_df.ESKIN.values

    data = {
        "date": np.concatenate(
            (
                irri_df.index,
                noirri_df.index,
                irri_df.index,
                noirri_df.index,
                irri_df.index,
                noirri_df.index,
                irri_df.index,
                noirri_df.index,
            )
        ),
        "[mm]": np.concatenate((a_1, a_2, b_1, b_2, c_1, c_2, d_1, d_2)),
        "variable": ["evapotranspiration"] * len(a_1)
        + ["evapotranspiration"] * len(a_2)
        + ["evaporation of soil"] * len(b_1)
        + ["evaporation of soil"] * len(b_2)
        + ["transpiration of vegetation"] * len(c_1)
        + ["transpiration of vegetation"] * len(c_2)
        + ["evaporation of skin reservoir"] * len(d_1)
        + ["evaporation of skin reservoir"] * len(d_2),
        "simulation": ["irrigated"] * len(a_1)
        + ["not irrigated"] * len(a_2)
        + ["irrigated"] * len(b_1)
        + ["not irrigated"] * len(b_2)
        + ["irrigated"] * len(c_1)
        + ["not irrigated"] * len(c_2)
        + ["irrigated"] * len(d_1)
        + ["not irrigated"] * len(d_2),
    }

    ax4 = fig.add_subplot(4, 1, 4)
    plt.rcParams.update(params)
    p = sns.lineplot(
        data=data,
        x="date",
        y="[mm]",
        hue="variable",
        style="simulation",
        palette=colors,
        ax=ax4,
    )
    p.legend_.remove()
    p.xaxis.set_major_locator(mdates.MonthLocator(bymonthday=[1, 3, 5, 7]))
    datelist_ticks = datelist[::2]
    p.set_xticklabels(datelist_ticks.strftime("%d.%m"), rotation=45)
    ax4.grid(True)
    # ax4.set_title(title_list[i])
    ax4.set_ylim(-1, 4.2)
    ax4.set_ylabel("evapotranspiration fractions \n [mmd$^{-1}$]")
    p.legend(bbox_to_anchor=(1.0, 1.05))  # loc='upper right')
plt.savefig(str(dir_out) + "/Fig16.png", dpi=300, bbox_inches="tight")
plt.show()
