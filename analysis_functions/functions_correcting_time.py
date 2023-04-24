#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 08:56:35 2023

@author: g300099
"""
from datetime import datetime

import numpy as np
import pandas as pd


def correct_timedim(data):
    datnew = []
    timenew = []
    datlist = list(map(str, np.floor(data.time.values).astype(int)))
    timelist = [
        str(item)
        for item in [
            int(item)
            for item in (round(num) for num in [t % 1 * 24 for t in data.time.values])
        ]
    ]

    for i in range(len(datlist)):
        datnew.append(datetime.strptime(datlist[i], "%Y%m%d").date())
        timenew.append(datetime.strptime(timelist[i], "%H").time())
    date_time = pd.DataFrame([f"{d} {t}" for d, t in zip(datnew, timenew)])
    time_corr = pd.to_datetime(date_time[0][:])
    dsnew = data.assign_coords({"time": ("time", time_corr)})
    return dsnew


def correct_timedim_mfiles(data):
    datnew = []
    datlist = list(map(str, (np.floor(data.time.values) / 100).astype(int)))

    for i in range(len(datlist)):
        datnew.append(datetime.strptime(datlist[i], "%Y%m").date())
        time_corr = pd.to_datetime(datnew[:])

    dsnew = data.assign_coords(
        {"time": time_corr}
    )  # .swap_dims({'dim_0':'time'}).drop_vars('dim_0')
    return dsnew
