#!/usr/bin/env python

# this basically makes it so that the machine doesn't have to figure
# out x-windows stuff:
import matplotlib as mpl
mpl.use('Agg')

from datetime import datetime
from datetime import timedelta
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from gitm_routines import *
from pylab import cm
import argparse

#-----------------------------------------------------------------------------
# get arguments from the user
#-----------------------------------------------------------------------------

def get_args():
    parser = argparse.ArgumentParser(
        description = 'Calculate O/N2 from model results')
    parser.add_argument('-on2',  \
                        action='store_true', default = False, \
                        help = 'plot only on2')
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
# Make a plot with a color bar
#-----------------------------------------------------------------------------

def make_plot(ax, lons, lats, vals, var, doLabelY = True, minmax = []):

    if (len(minmax) == 2):
        mini = minmax[0]
        maxi = minmax[1]
    else:
        maxi = np.max(np.abs(vals))
        mini = np.min(vals)

    cmap = mpl.cm.plasma        
    cax = ax.pcolor(lons, lats, vals, vmin = mini, vmax = maxi, cmap = cmap)

    ax.set_aspect(1)
    ax.set_ylim(-90, 90)
    ax.set_xlim(0, 360)
    ax.set_title(var)
    ax.set_xlabel('Longitude (deg)')
    if (doLabelY):
        ax.set_ylabel('Latitude (deg)')
    else:
        ax.set_yticks([])
    
    cbar = fig.colorbar(cax, ax = ax, shrink = 0.5, pad=0.01)
    #cbar.set_label(var, rotation=90)

    return


# -------------------------------------------------------------------
# This gives the axes for multi-row, multi-column plots.
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

def get_axes(figIn,
             nPlots,
             yBottom,
             yTop,
             yBuffer,
             xLeft,
             xRight,
             xBuffer):
    
    n = np.sqrt(nPlots)
    nX = int(np.ceil(n))
    nY = int(np.ceil(nPlots / nX))

    # plot size in y direction:
    ySize = (1.0 - yBottom - yTop) / nY - yBuffer * (nY - 1) / nY
    xSize = (1.0 - xLeft - xRight) / nX - xBuffer * (nX - 1) / nX
    
    ax = []
    for iPlot in range(nPlots):
        iX = int(iPlot % nX)
        iY = int(iPlot / nX)
        # I want to reverse y, so that 0 is the top left plot:
        iY = nY - iY - 1
        ax.append(figIn.add_axes([xLeft + iX * (xSize + xBuffer),
                                  yBottom + iY * (ySize + yBuffer),
                                  xSize,
                                  ySize]))
    
    return ax

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Main Code!
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

# define some constants
rtod =  180.0 / np.pi
dtor = 1.0 / rtod

# Get the input arguments
args = get_args()

#headers = read_gitm_headers('3DALL')
headers = read_gitm_headers(files = args.filelist)

nGitmFiles = len(headers["time"])

iFile = 0

# get lons, lats, alts:
iVars_ = [0,1,2]
data = read_gitm_one_file(headers["filename"][iFile], iVars_)
Alts3d = data[2]
Alts = data[2][0,0,:]/1000.0;
Lons2d = data[0][:,:,0]*rtod;
Lats2d = data[1][:,:,0]*rtod;
[nLons, nLats, nAlts] = data[0].shape

vars = remap_variable_names(headers["vars"])

iO_ = vars.index('[O] (/m3)')
iN2_ = vars.index('[N2] (/m3)')
ie_ = vars.index('[e-] (/m3)')

on2all = []
tecall = []
nmf2all = []
hmf2all = []
times = []

for iFile, file in enumerate(headers["filename"]):

    # Read O and N2:
    iVars_ = [iO_, iN2_, ie_]
    data = read_gitm_one_file(file, iVars_)
    times.append(headers["time"][iFile])
    oDensity = data[iO_]
    n2Density = data[iN2_]
    eDensity = np.array(data[ie_])

    # 1e21/m2 is the integration boundary.  
    n2Int = vertically_integrate(n2Density, Alts3d, calc3D = True)
    oInt = vertically_integrate(oDensity, Alts3d, calc3D = True)
    tec = vertically_integrate(eDensity, Alts3d)/1e16

    on2 = np.zeros((nLons, nLats))
    nmf2 = np.zeros((nLons, nLats))
    hmf2 = np.zeros((nLons, nLats))

    iAlts = np.arange(nAlts)
    for iLat in range(nLats):
        for iLon in range(nLons):
            n21d = n2Int[iLon,iLat,:]/1e21
            o1d = oInt[iLon,iLat,:]
            e1d = eDensity[iLon,iLat,:]

            i = iAlts[n21d < 1.0][0]
            r = (1.0 - n21d[i]) / (n21d[i-1] - n21d[i])
            n2 = (r * n21d[i-1] + (1.0 - r) * n21d[i]) * 1e21
            o = r * o1d[i-1] + (1.0 - r) * o1d[i]
            on2[iLon, iLat] = o / n2
            ie = np.argmax(e1d)
            nmf2[iLon, iLat] = e1d[ie]
            hmf2[iLon, iLat] = Alts3d[iLon, iLat, ie]/1000.0

    on2all.append(on2)
    tecall.append(tec)
    nmf2all.append(nmf2)
    hmf2all.append(hmf2)


on2MM = [np.min(on2all), np.max(on2all)]

for iFile, file in enumerate(headers["filename"]):

    fig = plt.figure(figsize=(10, 6))

    ax = []
    cax = []
    cbar = []
    time = times[iFile]
    sTime = time.strftime('%y%m%d_%H%M%S')
    plotTime = time.strftime(' at %B %d, %Y %H%M UT')

    if (not args.on2):
        nPlots = 4
        yBot = 0.00
        yTop = 0.00
        yBuf = 0.00
        xLef = 0.065
        xRig = 0.00
        xBuf = 0.00
        lead = 'all_'
    else:
        nPlots = 1
        yBot = 0.00
        yTop = 0.00
        yBuf = 0.00
        xLef = 0.065
        xRig = 0.00
        xBuf = 0.00
        lead = 'on2_'
    ax = get_axes(fig, nPlots, yBot, yTop, yBuf, xLef, xRig, xBuf)

    make_plot(ax[0], Lons2d, Lats2d, on2all[iFile], \
              'O/N2 Ratio' + plotTime, doLabelY = True,
              minmax = on2MM)
    if (not args.on2):
        make_plot(ax[1], Lons2d, Lats2d, tec, 'TEC (1e16/m2)', doLabelY = False)
        make_plot(ax[2], Lons2d, Lats2d, nmf2/1e12, 'NMF2 (1e12/m3)')
        make_plot(ax[3], Lons2d, Lats2d, hmf2, 'HMF2 (km)', doLabelY = False)

    outfile = lead + sTime + '.png'
    print('--> Writing file : ' + outfile)
    plt.savefig(outfile)
    plt.close()

