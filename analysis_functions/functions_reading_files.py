#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 09:05:55 2023

@author: g300099
"""

import os

import xarray as xr


def read_efiles(dir_path, var, var_num, exp_number, year, month):
    data_path = (
        str(dir_path)
        + "/"
        # "/work/ch0636/g300099/SIMULATIONS/GAR11/remo_results/"
        + str(exp_number)
        # + "/2017/"
    )
    if month >= 1 and month <= 12:
        efiles = (
            "hourly/var_series/"
            + str(var)
            + "/e"
            + str(exp_number)
            + "e_c"
            + str(var_num)
            + "_20170"
            + str(month)
            + ".nc"
        )
    else:
        efiles = (
            "hourly/var_series/"
            + str(var)
            + "/e"
            + str(exp_number)
            + "e_c"
            + str(var_num)
            + "_2017*.nc"
        )
    data = xr.open_mfdataset(os.path.join(data_path, efiles))  # , parallel=True)
    return data


def read_tfiles(exp_number, year, month):
    dir_path = (
        "/work/ch0636/g300099/SIMULATIONS/GAR11/remo_results/"
        + str(exp_number)
        + "/2017/"
    )
    if month >= 1 and month < 10:
        tfiles = (
            "xt/e" + str(exp_number) + "t" + str(year) + "0" + str(month) + "*.nc"
        )  #
    elif month >= 10 and month <= 12:
        tfiles = "xt/e" + str(exp_number) + "t" + str(year) + str(month) + "*.nc"  #
    else:
        tfiles = "xt/e" + str(exp_number) + "t" + str(year) + "*.nc"  #
    data = xr.open_mfdataset(os.path.join(dir_path, tfiles))  # , parallel=True)
    return data


def read_mfiles(dir_path, exp_number, year, month):
    data_path = str(dir_path) + "/" + str(exp_number) + "/monthly/"
    if month >= 1 and month <= 12:
        mfiles = "e" + str(exp_number) + "m20170" + str(month) + ".nc"
    else:
        mfiles = "e" + str(exp_number) + "m2017*.nc"
    data = xr.open_mfdataset(os.path.join(data_path, mfiles))  # , parallel=True)
    return data
