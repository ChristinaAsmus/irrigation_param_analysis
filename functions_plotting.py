#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 08:54:54 2023

@author: g300099
"""

import cartopy.feature as cfeature
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

def earth_feature(name_obj, res, col, cat):
    lines=cfeature.NaturalEarthFeature(category=cat, name=name_obj, scale=res, facecolor=col)
    return lines
def earth_feature_area(name_obj, res, cat,col):
    area=cfeature.NaturalEarthFeature(category=cat, name=name_obj, scale=res, facecolor=col)
    return area

states=earth_feature('admin_1_states_provinces_lines','10m','none','cultural')
land=earth_feature_area('land','10m','physical','none')
rivers=earth_feature('rivers_lake_centerlines','10m','none','physical')
coastline=earth_feature('coastline','10m','none','physical')
borders=earth_feature('admin_0_boundary_lines_land','10m','none','cultural')
ocean=earth_feature_area('ocean','10m','physical','none')

              
def plot_rotvar(fig, varvalue, axs, cbaxes, unit, label, cmap, levels, ticks, extend_scale, cbar_orient):
    axs.add_feature(land, zorder=1)
    axs.add_feature(coastline, edgecolor='black', linewidth=0.6)
    axs.add_feature(borders, edgecolor='black', linewidth=0.6)
    axs.set_xmargin(0)
    axs.set_ymargin(0)
    cmap = plt.get_cmap(cmap)
    rotplot = varvalue.plot.pcolormesh(
        ax=axs, cmap=cmap, levels=levels, extend=extend_scale, add_colorbar=False)
    #colorbar below the plot 
    cbar=fig.colorbar(rotplot, cax=cbaxes, orientation=cbar_orient, label=label, ticks=ticks)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(ticks)
    cbar.outline.set_visible(False)
    axs.gridlines(linewidth=0.7, color='gray', alpha=0.8, linestyle='--', zorder=3)
    return rotplot

rotated_pole = ccrs.RotatedPole(pole_latitude=39.25, pole_longitude=-162)


def draw_rectangel(background,x1, x2, y1,y2):
    x_Italy=[background.rlon[x1],background.rlon[x2],background.rlon[x2],background.rlon[x1],background.rlon[x1]]
    y_Italy=[background.rlat[y1],background.rlat[y1],background.rlat[y2],background.rlat[y2],background.rlat[y1]]
    return x_Italy, y_Italy