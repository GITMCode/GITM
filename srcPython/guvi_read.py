#!/usr/bin/env python

import matplotlib as mpl
mpl.use('Agg')

import datetime as dt
import numpy as np
from scipy.io import readsav
import matplotlib.pyplot as plt
import matplotlib.dates as dates

def convert_fractional_date_to_datetime(year, days):
    basetime = dt.datetime(year, 1, 1, 0, 0, 0)
    times = []
    for doy in days:
        seconds = (doy - 1) * 24.0 * 3600.0
        times.append(basetime + dt.timedelta(seconds = seconds))
    return times

def read_guvi_sav_file(savFile):

    print('Reading file : ', savFile)
    guvi = readsav(savFile)
    nOrbits = len(guvi['saved_data'])
    year = int(guvi['year'])

    allLats = []
    allLons = []
    allOn2 = []
    allTimes = []
    startTimes = []
    for iOrbit in range(nOrbits):
        nPts = guvi['saved_data'][iOrbit][5]
        lats = guvi['saved_data'][iOrbit][0][0:nPts]
        lons = guvi['saved_data'][iOrbit][1][0:nPts]
        on2 =  guvi['saved_data'][iOrbit][3][0:nPts]
        dates = guvi['saved_data'][iOrbit][4][0:nPts]
        times = convert_fractional_date_to_datetime(year, dates)
        allLats = np.concatenate((allLats, lats))
        allLons = np.concatenate((allLons, lons))
        allOn2 = np.concatenate((allOn2, on2))
        allTimes = np.concatenate((allTimes, times))
        startTimes.append(times[0])

    periods = []
    for i, sts in enumerate(startTimes[1:]):
        periods.append((sts - startTimes[i]).total_seconds())
    data = {'lats' : allLats,
            'lons' : allLons,
            'times' : allTimes,
            'on2' : allOn2,
            'period' : np.array(periods)}
    return data

def merge_dicts(dict1, dict2):
    combined = {}
    for key in dict1.keys():
        combined[key] = np.concatenate(( dict1[key], dict2[key]))
    return combined

def copy_dict(dict1):
    combined = {}
    for key in dict1.keys():
        combined[key] = dict1[key]
    return combined

def read_list_of_guvi_files(filelist):
    nFiles = len(filelist)
    if (nFiles == 1):
        alldata = read_guvi_sav_file(filelist[0])
    else:
        for iFile, file in enumerate(filelist):
            onefile = read_guvi_sav_file(file)
            if (iFile > 0):
                print(' -> merging')
                alldata = merge_dicts(alldata, onefile)
            else:
                print(' -> copying')
                alldata = copy_dict(onefile)
    return alldata


def plot_guvi(guvi):
    
    fig = plt.figure(figsize = (10,10))
    plt.rcParams.update({'font.size': 14})

    ax = fig.add_subplot(111)
    ax.set_position([0.1, 0.1, 0.8, 0.8])

    lats = guvi['lats']
    lons = guvi['lons']
    on2 = guvi['on2']
    times = guvi['times']
    cmap = mpl.cm.plasma
    cax = ax.scatter(lons, lats, c = on2, cmap = cmap)
    cbar = fig.colorbar(cax, ax = ax, shrink = 0.5, pad=0.01)

    sTime = times[0].strftime('GUVI O/N2 : %b %d, %Y %H:%M UT') + ' - ' + \
        times[-1].strftime('%b %d, %Y %H:%M UT')
    ax.set_title(sTime)

    plotfile = times[-1].strftime('guvi_on2%Y%m%d.png')

    print('Writing plotfile : ', plotfile)
    fig.savefig(plotfile)
    plt.close()
