#!/usr/bin/env python3

from datetime import datetime
from datetime import timedelta
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from glob import glob
import os
from gitm_routines import *
from scipy.stats import pearsonr
import argparse

# ------------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------------

def get_args():

    parser = argparse.ArgumentParser(
        description = 'Plot Aether / GITM model results')
    
    parser.add_argument('-list',  \
                        action='store_true', default = False, \
                        help = 'list variables in file')

    parser.add_argument('-var',  \
                        default = 3, type = int, \
                        help = 'variable to plot (number)')

    parser.add_argument('-alt', metavar = 'alt', default = 350.0, type = int, \
                        help = 'altitude :  alt in km (closest)')

    parser.add_argument('filelist', nargs='+', \
                        help = 'list files to use for generating plots')
    
    args = parser.parse_args()

    return args
    

def my_correlate(a, b):
    af = a.flatten()
    bf = b.flatten()
    a = (af - np.mean(af)) / (np.std(af) * len(af))
    b = (bf - np.mean(bf)) / (np.std(bf))
    c = pearsonr(a, b)
    return c[0]

# ------------------------------------------------------------------------------------
# Main code
# ------------------------------------------------------------------------------------

# Get the input arguments
args = get_args()

targetAlt = args.alt

files = args.filelist

if (args.list):
    data = read_gitm_one_file(files[0], [0])
    for iVar, var in enumerate(data['vars']):
        print('%2d. ' % iVar, var.decode('utf-8').replace(" ",""))
    exit()

nTimes = 12 * 3
nStep = 6

rtod =  180.0 / np.pi
dtor = 1.0 / rtod

iVar = args.var
               
nTimesTotal = len(files)
nSteps = int((nTimesTotal - nTimes) / nStep)

# get lons, lats, alts:
vars = [0,1,2]
data = read_gitm_one_file(files[0], vars)
Alts = data[2][0][0]/1000.0;
Lons = data[0][:,0,0]*rtod;
Lats = data[1][0,:,0]*rtod;
Lats2d = data[1][:,:,0]*rtod;

nLons = len(Lons)
nLats = len(Lats)
nAlts = len(Alts)

allData = np.zeros((nTimesTotal, nLons, nLats))

dLon = Lons[1]-Lons[0]
dLat = Lats[1]-Lats[0]

dAlt = np.abs(Alts - targetAlt)
iAlt = np.argmin(dAlt)


correlations = np.zeros((nSteps, nTimes))
correlationsNorth = np.zeros((nSteps, nTimes))
correlationsSouth = np.zeros((nSteps, nTimes))
dts = np.zeros((nSteps, nTimes))
allTimes = []

for iTime, file in enumerate(files):
    data4d = read_gitm_one_file(file, [iVar])
    allData[iTime, :, :] = data4d[iVar][:, :, iAlt]
    allTimes.append( data4d["time"])

for i in range(nSteps):
    iStart = i * nStep
    data0 = allData[iStart, :, :]
    for iTime in range(nTimes):
        data = allData[iTime + iStart, :, :]
        correlations[i, iTime] = my_correlate(data, data0)
        correlationsNorth[i, iTime] = my_correlate(data[Lats2d > 45.0], data0[Lats2d > 45.0])
        correlationsSouth[i, iTime] = my_correlate(data[Lats2d < -45.0], data0[Lats2d < -45.0])
    
        dts[i, iTime] = (allTimes[iTime + iStart] - allTimes[iStart]).total_seconds()

corrMeans = np.zeros(nTimes)
corrMeansNorth = np.zeros(nTimes)
corrMeansSouth = np.zeros(nTimes)
for iTime in range(nTimes):
    corrMeans[iTime] = np.mean(correlations[:, iTime])
    corrMeansNorth[iTime] = np.mean(correlationsNorth[:, iTime])
    corrMeansSouth[iTime] = np.mean(correlationsSouth[:, iTime])

fig = plt.figure(figsize=(10, 8.5))

xLeft = 0.07
xSize = 1.0 - xLeft - 0.02
yBot = 0.06
nY = 3
ySpace = 0.03
yS = (1.0 - ySpace * nY - yBot) / nY

i = 2
axG = fig.add_axes([xLeft, yBot + (yS+ySpace) * i, xSize, yS])
i = 1
axN = fig.add_axes([xLeft, yBot + (yS+ySpace) * i, xSize, yS])
i = 0
axS = fig.add_axes([xLeft, yBot + (yS+ySpace) * i, xSize, yS])

for iStep in range(nSteps):
    axG.plot(dts[iStep, :]/60.0, correlations[iStep, :], color = 'k', alpha = 0.1)
    axN.plot(dts[iStep, :]/60.0, correlationsNorth[iStep, :], color = 'b', alpha = 0.1)
    axS.plot(dts[iStep, :]/60.0, correlationsSouth[iStep, :], color = 'r', alpha = 0.1)

axG.plot(dts[1, :]/60.0, corrMeans, label = 'Global', color = 'k')
axN.plot(dts[1, :]/60.0, corrMeansNorth, label = 'North', color = 'b')
axS.plot(dts[1, :]/60.0, corrMeansSouth, label = 'South', color = 'r')

axG.set_ylabel('Correlation (Global)')
axN.set_ylabel('Correlation (North Polar)')
axS.set_ylabel('Correlation (South Polar)')

axS.set_xlabel('Delta-T (min)')

sTime = allTimes[0].strftime('%b %d, %Y %H:%M UT') + ' - ' + \
    allTimes[-1].strftime('%b %d, %Y %H:%M UT')

sTimeOut = allTimes[0].strftime('%Y%m%d%H') + '-' + \
    allTimes[-1].strftime('%Y%m%d%H')

#varName = remap_variable_names(data4d['vars'][iVar].decode('utf-8').replace(" ",""))
varName = data4d['vars'][iVar].decode('utf-8').replace(" ","")
print(varName)

axG.set_title(varName + '; ' + sTime)

axG.set_ylim(0.0, 1.02)
axN.set_ylim(0.0, 1.02)
axS.set_ylim(0.0, 1.02)

# Report out things at 90 minutes:
t0 = 90.0 # minutes

diff = np.abs(t0 - dts[0, :]/60.0)
i90 = np.argmin(diff)

axG.axvline(t0, linestyle = ':')
axG.text(t0 + 1, 0.9, 'Mean Corr : %4.2f' % corrMeans[i90])
mini = np.min(correlations[:, i90])
axG.scatter(t0, mini, color = 'k')
axG.text(t0 + 1, 0.05, 'Min Corr = %4.2f' % mini)

axN.axvline(t0, linestyle = ':')
axN.text(t0 + 1, 0.9, 'Mean Corr : %4.2f' % corrMeansNorth[i90], color = 'b')
mini = np.min(correlationsNorth[:, i90])
axN.scatter(t0, mini, color = 'b')
axN.text(t0 + 1, 0.05, 'Min Corr = %4.2f' % mini, color = 'b')

axS.axvline(t0, linestyle = ':')
axS.text(t0 + 1, 0.9, 'Mean Corr : %4.2f' % corrMeansSouth[i90], color = 'r')
mini = np.min(correlationsSouth[:, i90])
axS.scatter(t0, mini, color = 'r')
axS.text(t0 + 1, 0.05, 'Min Corr = %4.2f' % mini, color = 'r')

fileout = 'corr_' + sTimeOut + '_var%02d.png' % iVar
print('Writing file : ' + fileout)
fig.savefig(fileout)

