#!/usr/bin/env python3

import re
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as dates
import argparse, os

from swmf_imf import *
from read_ae import *
from dst import *

#------------------------------------------------------------------------------
# Get closest day boundary
#------------------------------------------------------------------------------

def get_day_boundary(time):

    year = time.year
    month = time.month
    day = time.day
    tBoundary = dt.datetime(year, month, day)
    diff = (tBoundary - time).total_seconds()/3600.0
    if (diff < -12.0):
        tBoundary = tBoundary + dt.timedelta(days = 1)
        
    return tBoundary

#------------------------------------------------------------------------------
# 
#------------------------------------------------------------------------------

def plot_day_boundaries(ax, times):

    tStart = get_day_boundary(times[0])
    tEnd = get_day_boundary(times[-1])

    nDays = int((tEnd - tStart).total_seconds()/86400.0)

    for iDay in range(nDays):
        d = tStart + dt.timedelta(days=iDay)
        ax.axvline(d, color = 'k', ls = '--')

    ax.set_xlim(tStart, tEnd)
    if (nDays > 2):
        ax.xaxis.set_major_formatter(dates.DateFormatter('%m-%d'))
    else:
        ax.xaxis.set_major_formatter(dates.DateFormatter('%d-%H'))
    return


def main(imf=None, ae=None, dst=None):
    """
    Makes plots of IMF and/or AE

    inputs:
     
        imf: returned from swmf_imf.read_swmf_file()
        ae: returned from read_ae.readae()

    Outputs:
        nothing

    A plot is created called 'driversYYYYMMDD.png' (min date given)

    """

    fig = plt.figure(figsize = (7,10))
    nY = 4
    if imf is not None:
        zeros = np.array(imf["bz"]) * 0.0
    elif ae is not None:
        zeros = np.array(ae["ae"]) * 0.0
    else:
        raise ValueError("Neither AE nor IMF were provided!")

    if dst is not None:
        nY = nY + 1
    
    plt.rcParams.update({'font.size': 14})

    ax = []

    xStart = 0.15
    xSize = 0.81
    
    yGap = 0.04
    yEnd = 0.95

    ySize = (yEnd - 0.06 - yGap * (nY-1)) / nY 

    iY = 0

    if imf is not None:
        iY += 1
        ax.append(fig.add_subplot(nY*100 + 10 + iY))
        yStart = yEnd - ySize * iY - yGap * (iY-1)
        ax[-1].set_position([xStart, yStart, xSize, ySize])

        ax[-1].plot(imf["times"], imf["bz"])
        ax[-1].plot(imf["times"], zeros, 'k:')
        ax[-1].set_ylabel('(a) IMF Bz (nT)')
        ax[-1].set_title('IMF, Solar Wind, and AE Drivers for GITM')
        plot_day_boundaries(ax[-1], imf["times"])

        iY += 1
        ax.append(fig.add_subplot(nY*100 + 10 + iY))
        yStart = yEnd - ySize * iY - yGap * (iY-1)
        ax[-1].set_position([xStart, yStart, xSize, ySize])
        ax[-1].plot(imf["times"], imf["by"])
        ax[-1].plot(imf["times"], zeros, 'k:')
        ax[-1].set_ylabel('(b) IMF By (nT)')
        plot_day_boundaries(ax[-1], imf["times"])

        iY += 1
        ax.append(fig.add_subplot(nY*100 + 10 + iY))
        yStart = yEnd - ySize * iY - yGap * (iY-1)
        ax[-1].set_position([xStart, yStart, xSize, ySize])
        ax[-1].plot(imf["times"], imf["vx"])
        ax[-1].set_ylabel('(c) SW Vx (km/s)')
        plot_day_boundaries(ax[-1], imf["times"])
        mint = np.min(imf["times"])
        maxt = np.max(imf["times"])

    if ae is not None:
        iY += 1
        ax.append(fig.add_subplot(nY*100 + 10 + iY))
        yStart = yEnd - ySize * iY - yGap * (iY-1)
        ax[-1].set_position([xStart, yStart, xSize, ySize])
        ax[-1].plot(ae["time"], ae["ae"])
        ax[-1].set_ylabel('(d) AE (nT)')
        ax[-1].set_ylim(0.0, np.max(ae["ae"])*1.05)

        plot_day_boundaries(ax[-1], ae["time"])
        if (np.min(ae["time"]) > mint):
            mint = np.min(ae["time"])
        if (np.max(ae["time"]) < maxt):
            maxt = np.max(ae["time"])

    if dst is not None:
        iY += 1
        ax.append(fig.add_subplot(nY*100 + 10 + iY))
        yStart = yEnd - ySize * iY - yGap * (iY-1)
        ax[-1].set_position([xStart, yStart, xSize, ySize])
        ax[-1].plot(dst["times"], dst["dst"])
        ax[-1].set_ylabel('(e) Dst (nT)')
        maxi = np.max([np.max(dst["dst"]) * 1.05, 0.0])
        mini = np.floor(np.min(dst["dst"])/100.0) * 100.0
        ax[-1].set_ylim(mini, maxi)
        plot_day_boundaries(ax[-1], dst["times"])
        zeros = np.array(dst["dst"]) * 0.0
        ax[-1].plot(dst["times"], zeros, 'k:')
        
    sTime = mint.strftime('%b %d, %Y %H:%M UT') + ' - ' + \
        maxt.strftime('%b %d, %Y %H:%M UT')
    ax[-1].set_xlabel(sTime)
    for iY in range(nY):
        ax[iY].set_xlim(mint, maxt)

    plotfile = mint.strftime('drivers%Y%m%d.png')
    print('-> Writing: ', plotfile)
    fig.savefig(plotfile)
    plt.close()


def parse_args():
    parser = argparse.ArgumentParser(
        description = 'Plot AE and IMF data')

    parser.add_argument('-imf', type=str, default='imf.txt',
                        help='Path to IMF file')

    parser.add_argument('-ae', type=str, default='ae.txt',
                        help='Path to AE file')

    parser.add_argument('-dst', \
                        help='plot dst with drivers', \
                        action="store_true")
    
    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_args()

    if os.path.exists(args.imf):
        imf = read_swmf_file(args.imf)
    else:
        imf = None
    
    if os.path.exists(args.ae):
        ae = read_ae(args.ae)
    else:
        ae=None

    if (args.dst):
        sTime = imf["times"][0].strftime('%Y-%m-%d')
        print(sTime)
        dst = get_dst(times = [sTime])
    else:
        dst = None
        
    main(imf, ae, dst = dst)

