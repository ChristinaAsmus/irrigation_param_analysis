#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 13:07:15 2023

@author: g300099
"""

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


def idw_for_model(rlat_source, rlon_source, rlat_target, rlon_target, data):
    # adapted from Dueymes, G. 2020: https://www.guillaumedueymes.com/post/netcdf_interpolation/

    rlon_source2d, rlat_source2d = np.meshgrid(rlon_source, rlat_source)
    rlon_target2d, rlat_target2d = np.meshgrid(rlon_target, rlat_target)
    temp_source2d = data.flatten()

    def lon_lat_to_cartesian(lon, lat):
        # WGS 84 reference coordinate system parameters
        A = 6378.137  # major axis [km]
        E2 = 6.69437999014e-3  # eccentricity squared

        lon_rad = np.radians(lon)
        lat_rad = np.radians(lat)
        # convert to cartesian coordinates
        r_n = A / (np.sqrt(1 - E2 * (np.sin(lat_rad) ** 2)))
        x = r_n * np.cos(lat_rad) * np.cos(lon_rad)
        y = r_n * np.cos(lat_rad) * np.sin(lon_rad)
        z = r_n * (1 - E2) * np.sin(lat_rad)
        return x, y, z

    xs, ys, zs = lon_lat_to_cartesian(rlon_source2d.flatten(), rlat_source2d.flatten())
    xt, yt, zt = lon_lat_to_cartesian(rlon_target, rlat_target)

    tree = cKDTree(np.column_stack((xs, ys, zs)))
    # idw
    d, inds = tree.query(np.column_stack((xt, yt, zt)), k=4)
    w = 1.0 / d**2
    temp_idw = np.sum(w * temp_source2d[inds], axis=1) / np.sum(w, axis=1)

    df_inttemp = pd.DataFrame()
    df_inttemp["rlon"] = rlon_target
    df_inttemp["rlat"] = rlat_target
    df_inttemp["temp2m"] = temp_idw
    # da_inttemp=df_inttemp.set_index(['rlat','rlon']).to_xarray()
    return df_inttemp
