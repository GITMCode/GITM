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
from guvi_read import *
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
# make a time vs latitude map:
#-----------------------------------------------------------------------------

def make_map(ix, iy, values, lats, lsts, nX, nY):

    asc = 0.0
    ascN = 0
    nPts = len(ix)
    mapA = np.zeros([nX, nY+1])
    mapAc = np.zeros([nX, nY+1])
        
    for iPt in range(nPts-1):
        i = int(ix[iPt])
        j = int(iy[iPt])
        mapA[i][j] = mapA[i][j] + values[iPt]
        mapAc[i][j] = mapAc[i][j] + 1
        if (np.abs(lats[iPt]) < 45.0):
            asc = asc + lsts[iPt]
            ascN = ascN + 1
        
    if (ascN == 0):
        ascN = 1
    asc = asc / ascN
    mapA[mapAc > 0] = mapA[mapAc > 0] / mapAc[mapAc > 0]
    mapA[mapAc == 0] = np.nan

    return mapA, asc
    
#-----------------------------------------------------------------------------
# make a time vs latitude map:
#-----------------------------------------------------------------------------

def plot_map(mTime, mLats, gitmM, guviM, ltnode, node, \
             miniV, maxiV, sTimeL,
             miniLat, maxiLat, plotfile):

    fig = plt.figure(figsize = (10,10))
    plt.rcParams.update({'font.size': 14})

    ax1 = fig.add_axes([0.075, 0.535, 1.0, 0.42])
    ax2 = fig.add_axes([0.075, 0.06, 1.0, 0.42])

    cmap = mpl.cm.plasma

    cax1 = ax1.pcolor(mTime, mLats, guviM, \
                      vmin = miniV, vmax = maxiV, cmap = cmap)

    cbar1 = fig.colorbar(cax1, ax = ax1, shrink = 0.5, pad=0.01)
    ax1.set_ylim(miniLat, maxiLat)
    ax1.set_xlim(np.min(mTime), np.max(mTime))
    ax1.set_ylabel('Latitude (deg)')
    ax1.set_title('GUVI O/N2 ' + ltnode + node)

    cax2 = ax2.pcolor(mTime, mLats, gitmM, \
                      vmin = miniV, vmax = maxiV, cmap = cmap)
    cbar2 = fig.colorbar(cax2, ax = ax2, shrink = 0.5, pad=0.01)
    ax2.set_ylim(miniLat, maxiLat)
    ax2.set_xlim(np.min(mTime), np.max(mTime))
    ax2.set_xlabel(sTimeL)
    ax2.set_ylabel('Latitude (deg)')
    ax2.set_title('GITM O/N2')

    print('Writing plotfile : ', plotfile)
    fig.savefig(plotfile)
    plt.close()


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

# get lons, lats, alts:
iVars_ = [0,1,2]
iFile = 0
data = read_gitm_one_file(headers["filename"][iFile], iVars_)
Alts3d = data[2];
Alts = data[2][0,0,:]/1000.0;
Lons = data[0][:,0,0]*rtod;
Lats = data[1][0,:,0]*rtod;
[nLons, nLats, nAlts] = data[0].shape

dLon = Lons[1]-Lons[0]
dLat = Lats[1]-Lats[0]

vars = remap_variable_names(headers["vars"])

iO_ = vars.index('[O] (/m3)')
iN2_ = vars.index('[N2] (/m3)')

nGitmFiles = len(headers["time"])

on2all = []
times = []

for iFile, file in enumerate(headers["filename"]):

    # Read O and N2:
    iVars_ = [iO_, iN2_]
    data = read_gitm_one_file(file, iVars_)
    times.append(headers["time"][iFile])
    oDensity = data[iO_]
    n2Density = data[iN2_]

    intBnd = 1.0e21
    # 1e21/m2 is the integration boundary.  
    n2Int = vertically_integrate(n2Density, Alts3d, calc3D = True)
    oInt = vertically_integrate(oDensity, Alts3d, calc3D = True)
    on2 = np.zeros((nLons, nLats))

    iAlts = np.arange(nAlts)
    for iLat in range(nLats):
        for iLon in range(nLons):
            n21d = n2Int[iLon,iLat,:]/intBnd
            o1d = oInt[iLon,iLat,:]

            i = iAlts[n21d < 1.0][0]
            r = (1.0 - n21d[i]) / (n21d[i-1] - n21d[i])
            n2 = (r * n21d[i-1] + (1.0 - r) * n21d[i]) * intBnd
            o = r * o1d[i-1] + (1.0 - r) * o1d[i]
            on2[iLon, iLat] = o / n2

    on2all.append(on2)


#guvidir = '/Users/ridley/Data/Gitm/run.202405/guvi/'
guvidir = '/backup/Data/Guvi/'

# brute force!
oldday = -1
guvifilelist = []
for t in times:
    if (t.day != oldday):
        guvifile = guvidir + t.strftime('%Y/ON2_%Y_%jm.sav')
        guvifilelist.append(guvifile)
        oldday = t.day
        
guvi = read_list_of_guvi_files(guvifilelist)

period = np.mean(guvi['period'])

saveTimes = []
saveLats = []
saveLons = []
saveLts = []
guviOn2 = []
gitmOn2 = []

iAfter = 0
for iTime, time in enumerate(guvi['times']):

    if ((time >= times[0]) & (time <= times[-1])):

        while (time > times[iAfter]):
            iAfter = iAfter+1
            if (iAfter >= nGitmFiles-1):
                break
        
        iBefore = iAfter-1

        dt = (times[iAfter] - \
              times[iBefore]).total_seconds()
        xt = (time - times[iBefore]).total_seconds() / dt

        lon = guvi['lons'][iTime]
        lat = guvi['lats'][iTime]

        xLon = (lon-Lons[0])/dLon
        iLon = int(xLon)
        xLon = xLon - iLon
        
        yLat = (lat-Lats[0])/dLat
        jLat = int(yLat)
        yLat = yLat - jLat

        on2Before = \
            (1-xLon)*(1-yLat)*on2all[iBefore][iLon, jLat]+\
            (  xLon)*(1-yLat)*on2all[iBefore][iLon+1, jLat]+\
            (1-xLon)*(  yLat)*on2all[iBefore][iLon, jLat+1]+\
            (  xLon)*(  yLat)*on2all[iBefore][iLon+1, jLat+1]
        on2After = \
            (1-xLon)*(1-yLat)*on2all[iAfter][iLon, jLat]+\
            (  xLon)*(1-yLat)*on2all[iAfter][iLon+1, jLat]+\
            (1-xLon)*(  yLat)*on2all[iAfter][iLon, jLat+1]+\
            (  xLon)*(  yLat)*on2all[iAfter][iLon+1, jLat+1]
        
        guviOn2.append(guvi['on2'][iTime])
        gitmOn2.append((1-xt) * on2Before + xt * on2After)
        saveTimes.append(time)
        saveLats.append(lat)
        saveLons.append(lon)
        ut = time.hour + time.minute/60.0 + time.second/3600.0
        localtime = (lon/15.0 + ut) % 12.0
        saveLts.append(localtime)

nPts = len(saveTimes)
print(' -> Found ', nPts, ' GITM/GUVI points that match in time')

# now convert to maps:

StartTime = datetime(saveTimes[0].year, \
                     saveTimes[0].month, \
                     saveTimes[0].day, \
                     saveTimes[0].hour)

tInS = []
for t in saveTimes:
    tInS.append((t - StartTime).total_seconds())
tInS = np.array(tInS)

EndTimeS = tInS[-1]
dLat = 1.0
nX = int(np.round(EndTimeS/period))+1
nY = int(180.0/dLat)

mapTime = []
mapTime1D = []
mapLats = np.zeros([nX, nY+1])

for i in range(nX):
    mapLats[i,:] = np.arange(nY+1) * dLat - 90.0 #+ dLat/2.0
    dt = i * period - period/2
    temp = [ dt ] * (nY+1)
    mapTime.append(temp)
    mapTime1D.append(StartTime + timedelta(seconds = dt))
    
saveLats = np.array(saveLats)
mapTime = np.array(mapTime)/3600.0

ix = np.round(tInS/period)
iy = np.round((saveLats + 90.0)/dLat)

maxLat = 82.0

gitmOn2 = np.array(gitmOn2)
guviOn2 = np.array(guviOn2)
saveLats = np.array(saveLats)
saveLts = np.array(saveLts)
ix = np.array(ix)
iy = np.array(iy)
# combined ascending and descending nodes:
#gitmOn2Map, localtime = make_map(ix, iy, gitmOn2, saveLats, saveLts, nX, nY)
#guviOn2Map, localtime = make_map(ix, iy, guviOn2, saveLats, saveLts, nX, nY)

all = [gitmOn2, guviOn2]
maxi = np.max(all)
mini = np.min(all)

sTime = saveTimes[0].strftime('%b %d, %Y %H:%M UT') + ' - ' + \
    saveTimes[-1].strftime('%b %d, %Y %H:%M UT Hours')
minLat = np.min(saveLats)-5
maxLat = np.max(saveLats)+5

dlat = np.concatenate((saveLats[1:] - saveLats[0:-1],[0]))

# ascending only:
gitmOn2a = gitmOn2[dlat > 0]
guviOn2a = guviOn2[dlat > 0]
saveLatsa = saveLats[dlat > 0]
saveLtsa = saveLts[dlat > 0]
ixa = ix[dlat > 0]
iya = iy[dlat > 0]

gitmOn2MapA, localtimeA = \
    make_map(ixa, iya, gitmOn2a, saveLatsa, saveLtsa, nX, nY)
guviOn2MapA, localtimeA = \
    make_map(ixa, iya, guviOn2a, saveLatsa, saveLtsa, nX, nY)
sLt = '(Local Time = %4.1f, ' % localtimeA

outfile = saveTimes[-1].strftime('guvi_gitm_on2%Y%m%d_asc.png')
node = 'Ascending)'

plot_map(mapTime, mapLats, gitmOn2MapA, guviOn2MapA, sLt, node, \
             mini, maxi, sTime,
             minLat, maxLat, outfile)

# descending only:
gitmOn2d = gitmOn2[dlat < 0]
guviOn2d = guviOn2[dlat < 0]
saveLatsd = saveLats[dlat < 0]
saveLtsd = saveLts[dlat < 0]
ixd = ix[dlat < 0]
iyd = iy[dlat < 0]

gitmOn2MapD, localtimeD = \
    make_map(ixd, iyd, gitmOn2d, saveLatsd, saveLtsd, nX, nY)
guviOn2MapD, localtimeD = \
    make_map(ixd, iyd, guviOn2d, saveLatsd, saveLtsd, nX, nY)
sLt = '(Local Time = %4.1f, ' % localtimeD

outfile = saveTimes[-1].strftime('guvi_gitm_on2%Y%m%d_dec.png')
node = 'Descending)'

plot_map(mapTime, mapLats, gitmOn2MapD, guviOn2MapD, sLt, node, \
             mini, maxi, sTime,
             minLat, maxLat, outfile)

# output some statistics:

diff = gitmOn2 - guviOn2

globalMedianGuvi = np.median(guviOn2)
globalMedianGitm = np.median(gitmOn2)
globalMedianDiff = np.median(diff) / globalMedianGuvi * 100.0
print('Global Medians (GUVI, GITM) : ', globalMedianGuvi, globalMedianGitm, globalMedianDiff)

northMeanGuvi = np.mean(guviOn2[saveLats > 0.0])
northMeanGitm = np.mean(gitmOn2[saveLats > 0.0])
northRMS = np.sqrt(np.mean(diff**2)) / northMeanGuvi * 100.0
print('North Means (GUVI, GITM, RMS) : ', northMeanGuvi, northMeanGitm, northRMS)
northStdGuvi = np.std(guviOn2[saveLats > 0.0]) / northMeanGuvi * 100.0
northStdGitm = np.std(gitmOn2[saveLats > 0.0]) / northMeanGitm * 100.0
print('North Std (GUVI, GITM) : ', northStdGuvi, northStdGitm)
northPolarMeanGuvi = np.mean(guviOn2[saveLats > 40.0])
northPolarMeanGitm = np.mean(gitmOn2[saveLats > 40.0])
print('North Polar Mean (GUVI, GITM) : ', northPolarMeanGuvi, northPolarMeanGitm)
southPolarMeanGuvi = np.mean(guviOn2[saveLats < -40.0])
southPolarMeanGitm = np.mean(gitmOn2[saveLats < -40.0])
print('South Polar Mean (GUVI, GITM) : ', southPolarMeanGuvi, southPolarMeanGitm)
northPolarStdGuvi = np.std(guviOn2[saveLats > 40.0]) / northPolarMeanGuvi * 100.0
northPolarStdGitm = np.std(gitmOn2[saveLats > 40.0]) / northPolarMeanGitm * 100.0
print('North Polar Std (GUVI, GITM) : ', northPolarStdGuvi, northPolarStdGitm)
southPolarStdGuvi = np.std(guviOn2[saveLats < -40.0]) / southPolarMeanGuvi * 100.0
southPolarStdGitm = np.std(gitmOn2[saveLats < -40.0]) / southPolarMeanGitm * 100.0
print('South Polar Std (GUVI, GITM) : ', southPolarStdGuvi, southPolarStdGitm)
