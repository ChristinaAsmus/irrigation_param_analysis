#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 08:28:55 2022

@author: christina
"""


import math


def rotated_coord_transform(lon, lat, np_lon, np_lat, direction="rot2geo"):
    """Transforms a coordinate into a rotated grid coordinate and vice versa.

    The coordinates have to given in degree and will be returned in degree.

    **Arguments:**
        *lon:*
            Longitude coordinate.
        *lat:*
            Latitude coordinate.
        *np_lon:*
            Longitude coordinate of the rotated pole.
        *np_lat:*
            Latitude coordinate of the rotated pole.
        *direction:*
            Direction of the rotation.
            Options are: 'rot2geo' (default) for a transformation to regular
            coordinates from rotated. 'geo2rot' transforms regular coordinates
            to rotated.

    **Returns:**
        *lon_new:*
            New longitude coordinate.
        *lat_new:*
            New latitude coordinate.

    Written by Kevin Sieck
    """

    # Convert degrees to radians
    lon = (lon * math.pi) / 180.0
    lat = (lat * math.pi) / 180.0

    #    SP_lon = SP_coor(1)
    #    SP_lat = SP_coor(2)

    theta = 90.0 - np_lat  # Rotation around y-axis
    phi = np_lon + 180.0  # Rotation around z-axis

    # Convert degrees to radians
    phi = (phi * math.pi) / 180.0
    theta = (theta * math.pi) / 180.0

    # Convert from spherical to cartesian coordinates
    x = math.cos(lon) * math.cos(lat)
    y = math.sin(lon) * math.cos(lat)
    z = math.sin(lat)

    # Regular -> Rotated
    if direction == "geo2rot":
        x_new = (
            math.cos(theta) * math.cos(phi) * x
            + math.cos(theta) * math.sin(phi) * y
            + math.sin(theta) * z
        )
        y_new = -math.sin(phi) * x + math.cos(phi) * y
        z_new = (
            -math.sin(theta) * math.cos(phi) * x
            - math.sin(theta) * math.sin(phi) * y
            + math.cos(theta) * z
        )

    # Rotated -> Regular
    elif direction == "rot2geo":
        phi = -phi
        theta = -theta

        x_new = (
            math.cos(theta) * math.cos(phi) * x
            + math.sin(phi) * y
            + math.sin(theta) * math.cos(phi) * z
        )
        y_new = (
            -math.cos(theta) * math.sin(phi) * x
            + math.cos(phi) * y
            - math.sin(theta) * math.sin(phi) * z
        )
        z_new = -math.sin(theta) * x + math.cos(theta) * z

    # Convert cartesian back to spherical coordinates
    lon_new = math.atan2(y_new, x_new)
    lat_new = math.asin(z_new)

    # Convert radians back to degrees
    lon_new = (lon_new * 180.0) / math.pi
    lat_new = (lat_new * 180.0) / math.pi

    return (lon_new, lat_new)
