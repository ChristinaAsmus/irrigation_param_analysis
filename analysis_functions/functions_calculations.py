#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 14:35:00 2023

@author: g300099
"""
import numpy as np


# convective and large-scale rain has to be summed up to show total rain
def calculate_sum_precip(var1, var2, dsnoirr_newtime, rlat_list, rlon_list, i):
    precip = (
        dsnoirr_newtime.isel(
            rlat=slice(rlat_list[i][0], rlat_list[i][1]),
            rlon=slice(rlon_list[i][0], rlon_list[i][1]),
        )[var1]
    ) + (
        dsnoirr_newtime.isel(
            rlat=slice(rlat_list[i][0], rlat_list[i][1]),
            rlon=slice(rlon_list[i][0], rlon_list[i][1]),
        )[var2]
    )
    precip_mean = (
        precip.resample(time="D").sum().mean(dim=["rlon", "rlat"], skipna=True)
    )
    return precip_mean


def calculate_bowen_ratio(array):
    bowen = -(array["AHFSIRR"]) / (-array["AHFLIRR"])
    return bowen


def calculate_hourmeans_diff(var, irrilimit, irri_array, noirri_array):
    hourdiff = (
        (
            (irri_array.where(irri_array["IRRIFRAC"] > irrilimit)[var]).mean(
                dim=["rlon", "rlat"], skipna=True
            )
        )
        .groupby("time.hour")
        .mean(skipna=True)
    ) - (
        (
            (noirri_array.where(noirri_array["IRRIFRAC"] > irrilimit)[var]).mean(
                dim=["rlon", "rlat"], skipna=True
            )
        )
        .groupby("time.hour")
        .mean(skipna=True)
    )
    return hourdiff


# monthly
def calculate_monthmeans_diff(var, irrilimit, irri_array, noirri_array):
    monthdiff = (
        (irri_array[var].where(irri_array["IRRIFRAC"] > irrilimit))
        - (noirri_array[var].where(noirri_array["IRRIFRAC"] > irrilimit))
    ).mean(dim=["rlon", "rlat"], skipna=True)
    return monthdiff


def calculate_meandiff(var, ds_irri, ds_noirri):
    varirrmean = ds_irri[var].mean("time")
    varnoirrmean = ds_noirri[var].mean("time")
    diff = varirrmean - varnoirrmean
    return diff


def calculate_meanvar(var, data):
    varmean = data[var].mean("time")
    return varmean


def calculate_monthlysum_evapos(var, irrilimit, irri_array, noirri_array):
    diffmonth = (
        (irri_array[var].where(irri_array["IRRIFRAC"] > irrilimit))
        .groupby("time.month")
        .sum()
    ).mean(dim=["rlon", "rlat"], skipna=True) - (
        (noirri_array[var].where(noirri_array["IRRIFRAC"] > irrilimit))
        .groupby("time.month")
        .sum()
    ).mean(
        dim=["rlon", "rlat"], skipna=True
    )
    return diffmonth


def calculate_hourlymean_evapos(var, irrilimit, irri_array, noirri_array):
    diffhour = (
        (irri_array[var].where(irri_array["IRRIFRAC"] > irrilimit))
        .groupby("time.hour")
        .mean()
    ).mean(dim=["rlon", "rlat"], skipna=True) - (
        (noirri_array[var].where(noirri_array["IRRIFRAC"] > irrilimit))
        .groupby("time.hour")
        .mean()
    ).mean(
        dim=["rlon", "rlat"], skipna=True
    )
    return diffhour


def calculate_monthlysum_mean_evapos(var, irrilimit, irri_array, noirri_array):
    diffmonth = (irri_array[var].where(irri_array["IRRIFRAC"] > irrilimit)).groupby(
        "time.month"
    ).sum().mean("month", skipna=True) - (
        (noirri_array[var].where(noirri_array["IRRIFRAC"] > irrilimit))
        .groupby("time.month")
        .sum()
    ).mean(
        "month", skipna=True
    )
    return diffmonth


def calculate_abs_monthlysum_evapos(var, data_array):
    datamonth_sum_mean = (
        (-data_array[var]).groupby("time.month").sum().mean("month", skipna=True)
    )
    return datamonth_sum_mean


def prepare_var(var, i, irrilimit, rlat_list, rlon_list, irri_array, noirri_array):
    sel_irri = (
        (
            irri_array.isel(
                rlat=slice(rlat_list[i][0], rlat_list[i][1]),
                rlon=slice(rlon_list[i][0], rlon_list[i][1]),
            )[var].where(irri_array["IRRIFRAC"] > irrilimit)
        )
        .mean(dim=["rlon", "rlat"], skipna=True)
        .groupby("time.hour")
        .mean()
    )
    sel_noirri = (
        (
            noirri_array.isel(
                rlat=slice(rlat_list[i][0], rlat_list[i][1]),
                rlon=slice(rlon_list[i][0], rlon_list[i][1]),
            )[var].where(noirri_array["IRRIFRAC"] > irrilimit)
        )
        .mean(dim=["rlon", "rlat"], skipna=True)
        .groupby("time.hour")
        .mean()
    )
    return sel_irri, sel_noirri


def calculate_rh(tdew, t):
    rh = 100 * (
        (np.exp((17.625 * tdew) / (243.04 + tdew)))
        / (np.exp((17.625 * t) / (243.04 + t)))
    )
    return rh


# calculate monthly sum
def calculate_monthly_sum_precip(var1, var2, irrilimit, dsirr_newtime, dsnoirr_newtime):
    dsirr_precip = dsirr_newtime[var1] + dsirr_newtime[var2]
    dsirr_month = (
        (dsirr_precip.where(dsirr_newtime["IRRIFRAC"] > irrilimit))
        .groupby("time.month")
        .sum()
    ).mean(dim=["rlon", "rlat"], skipna=True)
    dsnoirr_precip = dsnoirr_newtime[var1] + dsnoirr_newtime[var2]
    dsnoirr_month = (
        (dsnoirr_precip.where(dsnoirr_newtime["IRRIFRAC"] > irrilimit))
        .groupby("time.month")
        .sum()
    ).mean(dim=["rlon", "rlat"], skipna=True)
    diff_month = dsirr_month - dsnoirr_month
    return dsirr_month, dsnoirr_month, diff_month


# spatial
def calculate_meandiff_precip(array_irri, array_noirri, var1, var2):
    dsirr_precip = array_irri[var1] + array_irri[var2]
    dsirr_month = dsirr_precip.groupby("time.month").sum()
    dsnoirr_precip = array_noirri[var1] + array_noirri[var2]
    dsnoirr_month = dsnoirr_precip.groupby("time.month").sum()
    diff_month = (dsirr_month - dsnoirr_month).mean("month")
    return diff_month
