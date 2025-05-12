#!/usr/bin/env python

# this basically makes it so that the machine doesn't have to figure
# out x-windows stuff:
import matplotlib as mpl
mpl.use('Agg')

from datetime import datetime, timedelta
import numpy as np

import matplotlib.pyplot as plt
import glob
import os
import argparse

from gitm_routines import *
from vista_read import *

# -------------------------------------------------------------------
# This gives the axes for multi-row, but single column plots.
# Inputs:
#    - figIn: this is the figure (e.g., fig = plt.figure(figsize=(10, 10)))
#    - nPlots: number of plots in total
#    - yBottom: space desired below the bottom plot
#    - yTop: space desired above the top plot
#    - yBuffer: space desired between each plot
# Optional Inputs:
#    - xLeft: space desired to the left of the plots
#    - xRight: space desired to the right of the plots
# Outputs:
#    - ax: array of axes
# -------------------------------------------------------------------

def get_axes_one_column(figIn,
                        nPlots,
                        yBottom,
                        yTop,
                        yBuffer,
                        xLeft = 0.11,
                        xRight = 0.02):
    
    # plot size in y direction:
    ySize = (1.0 - yBottom - yTop) / nPlots - yBuffer * (nPlots - 1) / nPlots
    xSize = 1.0 - xLeft - xRight
    
    ax = []
    for iPlot in range(nPlots):
        # I want to reverse this, so that 0 is the top plot:
        i = nPlots - iPlot - 1
        ax.append(figIn.add_axes([xLeft,
                                  yBottom + i * (ySize + yBuffer),
                                  xSize,
                                  ySize]))
    
    return ax



#-----------------------------------------------------------------------------
# get arguments from the user
#-----------------------------------------------------------------------------

def get_args():
    parser = argparse.ArgumentParser(
        description = 'Compare GITM and VISTA results')
    parser.add_argument('filelist', nargs='+', \
                        help = 'list files to use for generating plots')
    args = parser.parse_args()
    return args

#-----------------------------------------------------------------------------
# vertically integrate the 3D data given the altitudes.
#-----------------------------------------------------------------------------

def vertically_integrate(value, alts, calc3D = False):
    [nLons, nLats, nAlts] = value.shape
    integrated = np.zeros((nLons, nLats, nAlts))
    descending = np.arange(nAlts-2, -1, -1)
    dz = alts[:,:,-1] - alts[:,:,-2]
    integrated[:,:,-1] = value[:,:,-1] * dz
    for i in descending:
        dz = alts[:,:,i+1] - alts[:,:,i]
        integrated[:,:,i] = integrated[:,:,i+1] + value[:,:,i] * dz
    if (not calc3D):
        integrated = integrated[:,:,0]
    return integrated
        
#-----------------------------------------------------------------------------
# 
#-----------------------------------------------------------------------------

def calc_tec_stats(time, lats, lons, tecmap):
    # calculate the average in different regions of the globe:
    #  - assume grid is in degrees

    ut = time.hour + time.minute/60.0 + time.second/3600.0
    lts = (lons/15.0 + ut) % 24
    
    lts2d, lats2d = np.meshgrid(lts, lats)

    #localtimes = lons2d/15.0 #(lons2d/15.0 + ut) % 24.0

    # assume a uniform grid to calculate the area of each cell:
    dlat = lats[1] - lats[0]
    dlon = lons[1] - lons[0]
    # take 120 km altitude, convert to m
    r = (6372.0 + 120.0) * 1000.0 
    area = r * r * dlat * dlon * np.cos(lats2d * np.pi / 180.0)

    stats = {}
    std = {}

    # total TEC:
    stats['global'] = np.sum(tecmap * area) / np.sum(area)
    std['global'] = np.std(tecmap)
    
    # north polar TEC:
    polar = 60.0
    stats['north_polar'] = np.sum(tecmap[lats2d >= polar] * \
                                  area[lats2d >= polar]) / \
                                  np.sum(area[lats2d >= polar])
    std['north_polar'] = np.std(tecmap[lats2d >= polar])
    
    stats['south_polar'] = np.sum(tecmap[lats2d <= -polar] * \
                                  area[lats2d <= -polar]) / \
                                  np.sum(area[lats2d <= -polar])
    std['south_polar'] = np.std(tecmap[lats2d <= -polar])

    mid = 30.0
    stats['equator'] = np.sum(tecmap[np.abs(lats2d) < mid] * \
                              area[np.abs(lats2d) < mid]) / \
                              np.sum(area[np.abs(lats2d) < mid])
    std['equator'] = np.std(tecmap[np.abs(lats2d) < mid])
    # make a normalized grid (north mid):
    dl = (polar - mid)/2
    c = (polar + mid)/2.0
    l2d = np.abs(lats2d - c) / dl
    areas = area[(l2d <= 1)]
    stats['north_mid'] = np.sum(tecmap[(l2d <= 1)] * \
                                  area[(l2d <= 1)]) / \
                                  np.sum(area[(l2d <= 1)])
    std['north_mid'] = np.std(tecmap[(l2d <= 1)])
    # make a normalized grid (south mid):
    dl = (polar - mid)/2
    c = -(polar + mid)/2.0
    l2d = np.abs(lats2d - c) / dl
    areas = area[(l2d <= 1)]
    stats['south_mid'] = np.sum(tecmap[(l2d <= 1)] * \
                                area[(l2d <= 1)]) / \
                                np.sum(area[(l2d <= 1)])
    std['south_mid'] = np.std(tecmap[(l2d <= 1)])

    lats2dmod = np.abs(lats2d) + 0.001
    tonoon = np.abs(lts2d - 12.0)
    # daytime is defined as tonoon < 5.0
    lats2dmod[tonoon > 5.0] = lats2dmod[tonoon > 5.0] + 99.0
    stats['day_eq'] = np.mean(tecmap[lats2dmod < mid])
    #stats['day_eq'] = np.sum(tecmap[lats2dmod < mid] * \
    #                          area[lats2dmod < mid]) / \
    #                          np.sum(area[lats2dmod < mid])
    std['day_eq'] = np.std(tecmap[lats2dmod < mid])

    lats2dmod = lats2d + 0.001
    tonoon = np.abs(lts2d - 12.0) 
    # nighttime is defined as tonoon > 7.0
    # regions to EXCLUDE:
    lats2dmod[tonoon < 7.0] = lats2dmod[tonoon < 7.0] + 99.0

    # make a normalized grid (north mid):
    polar = 50.0
    mid = 40.0
    dl = (polar - mid)/2
    c = (polar + mid)/2.0
    l2d = np.abs(lats2dmod - c) / dl
    areas = area[(l2d <= 1)]
    stats['n_mid_n'] = np.sum(tecmap[(l2d <= 1)] * \
                                  area[(l2d <= 1)]) / \
                                  np.sum(area[(l2d <= 1)])
    std['n_mid_n'] = np.std(tecmap[(l2d <= 1)])

    # make a normalized grid (north mid):
    dl = (polar - mid)/2
    c = -(polar + mid)/2.0
    l2d = np.abs(lats2dmod - c) / dl
    areas = area[(l2d <= 1)]
    stats['s_mid_n'] = np.sum(tecmap[(l2d <= 1)] * \
                                  area[(l2d <= 1)]) / \
                                  np.sum(area[(l2d <= 1)])
    std['s_mid_n'] = np.std(tecmap[(l2d <= 1)])

        
    return stats, std
    
#-----------------------------------------------------------------------------
# Main Code:
#-----------------------------------------------------------------------------

# define some constants
rtod =  180.0 / np.pi
dtor = 1.0 / rtod

dir = '/Users/ridley/Data/VISTA/'

file = dir + 'VISTA_110805.hdf5'

vista = read_vista_file(file)
vLats = vista['lats']
vLocalTimes = vista['localtimes']

# Get the input arguments
args = get_args()

headers = read_gitm_headers(files = args.filelist)

vars = remap_variable_names(headers["vars"])
ie_ = vars.index('[e-] (/m3)')

# get lons, lats, alts:
iVars_ = [0,1,2]
iFile = 0
data = read_gitm_one_file(headers["filename"][iFile], iVars_)
Alts3d = data[2];
Alts = data[2][0,0,:]/1000.0;
Lons = data[0][:,0,0]*rtod;
Lats = data[1][0,:,0]*rtod;
[nLons, nLats, nAlts] = data[0].shape

gTec = []
gTimes = []

for iFile, file in enumerate(headers["filename"]):

    iVars_ = [ie_]
    data = read_gitm_one_file(file, iVars_)
    gTimes.append(headers["time"][iFile])
    
    eDensity = np.array(data[ie_])
    tec = vertically_integrate(eDensity, Alts3d)/1e16
    gTec.append(tec)
    
time = vista['times'][0]
vTec = vista['tec'][:,:,0]

stats1 = calc_tec_stats(time, vLats, vLocalTimes*15.0, vTec)

def find_time(tArray, tGoal):
    tDiff = []
    for t in tArray:
        tDiff.append(np.abs((t - tGoal).total_seconds()))
    i_ = np.argmin(tDiff)
    return i_

keys = ['global', 'north_polar', 'north_mid', 'n_mid_n', \
        'equator', 'day_eq', \
        'south_mid', 's_mid_n', 'south_polar']

nStats = len(keys)
gStats = np.zeros((len(gTimes), nStats))
vStats = np.zeros((len(gTimes), nStats))
gStd = np.zeros((len(gTimes), nStats))
vStd = np.zeros((len(gTimes), nStats))

for iTime, gTime in enumerate(gTimes):
    
    iT_ = find_time(vista['times'], gTime)
    vTecLT = vista['tec'][:, :, iT_]

    vTime = vista['times'][iT_]
    UT = gTime.hour + gTime.minute/60.0 + gTime.second/3600.0
    #vLons = ((vLocalTimes - UT) * 15.0 + 360.0) % 360.0
    vLons = vLocalTimes * 15.0
    iShift = int(UT * 15.0)

    vTec = np.roll(vTecLT, shift = -iShift, axis = 1)

    fig = plt.figure(figsize = (10,10))
    plt.rcParams.update({'font.size': 14})

    ax1 = fig.add_axes([0.075, 0.535, 1.0, 0.42])
    ax2 = fig.add_axes([0.075, 0.06, 1.0, 0.42])

    cmap = mpl.cm.plasma

    all = [vTec, gTec[iTime]]
    maxi = np.max([np.max(vTec), np.max(gTec[iTime])])
    mini = 0.0
    
    cax1 = ax1.pcolor(vLons, vLats, vTec, cmap = cmap, vmin = mini, vmax = maxi)
    cbar1 = fig.colorbar(cax1, ax = ax1, shrink = 0.5, pad=0.01)
    #ax1.set_ylim(miniLat, maxiLat)
    #ax1.set_xlim(np.min(mTime), np.max(mTime))
    ax1.set_ylabel('Latitude (deg)')
    plotTime = vTime.strftime(' at %B %d, %Y %H%M UT')
    ax1.set_title('VISTA TEC' + plotTime)
    ax1.axhline(-60)
    ax1.axhline(-30)
    ax1.axhline(30)
    ax1.axhline(60)

    cax2 = ax2.pcolor(Lons, Lats, np.transpose(gTec[iTime]), \
                      cmap = cmap, shading = 'auto', vmin = mini, vmax = maxi)
    cbar2 = fig.colorbar(cax2, ax = ax2, shrink = 0.5, pad=0.01)
    ax2.set_ylabel('Latitude (deg)')
    ax2.set_xlabel('Longitude (deg)')
    plotTime = gTime.strftime(' at %B %d, %Y %H%M UT')
    ax2.set_title('GITM TEC' + plotTime)
    ax2.axhline(-60)
    ax2.axhline(-30)
    ax2.axhline(30)
    ax2.axhline(60)

    sTime = gTime.strftime('%y%m%d_%H%M%S')

    plotfile = 'vista_gitm_' + sTime + '.png'
    print('Writing plotfile : ', plotfile)
    fig.savefig(plotfile)
    plt.close()

    g, gs = calc_tec_stats(gTime, Lats, Lons, np.transpose(gTec[iTime]))
    v, vs = calc_tec_stats(vTime, vLats, vLons, vTec)
    for i, key in enumerate(keys):
        gStats[iTime, i] = g[key]
        gStd[iTime, i] = gs[key]
        vStats[iTime, i] = v[key]
        vStd[iTime, i] = vs[key]

fig = plt.figure(figsize = (10,10))
plt.rcParams.update({'font.size': 14})

yBottom = 0.06
yTop = 0.04
yBuffer = 0.005
ax = get_axes_one_column(fig,
                         len(keys),
                         yBottom,
                         yTop,
                         yBuffer,
                         xLeft = 0.11,
                         xRight = 0.02)

c = np.zeros((len(keys)))
r = np.zeros((len(keys)))

nT = len(gStats[:,i])
for i, key in enumerate(keys):
    ax[i].plot(gTimes, gStats[:,i], label = 'GITM', color = 'b')
    ax[i].plot(gTimes, vStats[:,i], label = 'VISTA', color = 'r')
    for iT in range(nT):
        x = [gTimes[iT], gTimes[iT]]
        y = gStats[iT, i] + [gStd[iT, i], -gStd[iT, i]]
        ax[i].plot(x, y, color = 'b', linewidth = 3, alpha = 0.2)
        y = vStats[iT, i] + [vStd[iT, i], -vStd[iT, i]]
        ax[i].plot(x, y, color = 'r', linewidth =2, alpha = 0.2)
                   
        
    ax[i].set_ylabel(key + '\nTEC')
    ax[i].legend()
    mini = np.min([gStats[:,i], vStats[:,i]])
    maxi = np.max([gStats[:,i], vStats[:,i]])
    c[i] = (maxi + mini)/2
    r[i] = (maxi - mini)

half = np.max([gStd, vStd])

range = np.max(r) * 1.1
if (half > range):
    range = half
    
for i, key in enumerate(keys):
    ax[i].set_ylim(c[i] - range/2, c[i] + range/2)
        
sTime = gTimes[0].strftime('%b %d, %Y %H:%M UT') + ' - ' + \
    gTimes[-1].strftime('%b %d, %Y %H:%M UT Hours')
    
ax[-1].set_xlabel(sTime)
ax[0].set_title('GITM Comparisons with VISTA')
    
sTime = gTimes[0].strftime('%y%m%d')
plotfile = 'vista_gitm_lines_' + sTime + '.png'
print('Writing plotfile : ', plotfile)
fig.savefig(plotfile)
plt.close()


package = {'times': gTimes,
           'alt': -1.0}
for i, key in enumerate(keys):
    package['gitm_' + key] = gStats[:, i]
    package['vista_' + key] = vStats[:, i]

write_log(package, fileHeader = 'vista_gitm',
          message = 'From thermo_vista.py')
