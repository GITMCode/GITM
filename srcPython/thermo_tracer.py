#!/usr/bin/env python

from datetime import datetime
from datetime import timedelta
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from gitm_routines import *
from numpy import random

def thermo_get_points(lats_in, lons_in, alts_in, times_in, verr = 0.0):

    vars_in = [16, 17, 18]
    
    rtod =  180.0 / np.pi
    dtor = 1.0 / rtod

    print("Reading gitm headers")
    headers = read_gitm_headers('3DALL')

    nGitmFiles = len(headers["time"])

    # get lons, lats, alts:
    vars = [0,1,2]
    data = read_gitm_one_file(headers["filename"][0], vars)
    Alts = data[2][0][0]/1000.0;
    Lons = data[0][:,0,0]*rtod;
    Lats = data[1][0,:,0]*rtod;
    [nLons, nLats, nAlts] = data[0].shape

    dLon = Lons[1]-Lons[0]
    dLat = Lats[1]-Lats[0]
        
    print("Looping through times...")

    times_out = []
    values_out = []

    iBefore0 = -1
    iAfter0 = -1
    # here we just assume that we will start with the first set of files.
    # this shouldn't matter...
    iBefore = 0
    iAfter = 1
    nVars = len(vars_in)
    AfterVals = np.zeros(nVars)
    BeforeVals = np.zeros(nVars)

    dt = (times_in[1] - times_in[0]).total_seconds()
    Re = 6372.0 * 1000.0
    cDtoR = np.pi / 180.0

    nTimes = len(times_in)
    nPts = len(lons_in)

    latsOut = np.zeros((nTimes+1, nPts))
    lonsOut = np.zeros((nTimes+1, nPts))
    altsOut = np.zeros((nTimes+1, nPts))

    latsOut[0, : ] = lats_in
    lonsOut[0, : ] = lons_in
    altsOut[0, : ] = alts_in
    iT = 0
    
    for i, time in enumerate(times_in):

        while (time > headers["time"][iAfter]):
            iAfter = iAfter+1
            if (iAfter >= nGitmFiles-1):
                break
        
        if (iAfter == nGitmFiles):
            break

        iBefore = iAfter-1

        if (iBefore != iBefore0):
            file = headers["filename"][iBefore]
            BeforeData = read_gitm_one_file(file, vars_in)
        if (iAfter != iAfter0):
            file = headers["filename"][iAfter]
            AfterData = read_gitm_one_file(file, vars_in)

        if (time >= headers["time"][iBefore]):

            times_out.append(time)
            
            dtFiles = (headers["time"][iAfter] - \
                       headers["time"][iBefore]).total_seconds()
            xt = (time - headers["time"][iBefore]).total_seconds() / dtFiles

            iT += 1
            for iPt in range(nPts):
                lon = lonsOut[iT-1, iPt]
                lat = latsOut[iT-1, iPt]
                alt = altsOut[iT-1, iPt]
                R = Re + alt * 1000.0
            
                xLon = (lon-Lons[0])/dLon
                iLon = int(xLon)
                xLon = xLon - iLon
        
                yLat = (lat-Lats[0])/dLat
                jLat = int(yLat)
                yLat = yLat - jLat

                kAlt = 0
                zAlt = 0.0
                if ((alt > Alts[0]) and (nAlts > 1)):
                    if (alt > Alts[nAlts-1]):
                        # above domain:
                        kAlt = nAlts-2
                        zAlt = 1.0
                    else:
                        while (Alts[kAlt] < alt):
                            kAlt = kAlt + 1
                        kAlt = kAlt - 1
                        zAlt = (alt - Alts[kAlt]) / (Alts[kAlt+1] - Alts[kAlt])
                    kAltp1 = kAlt + 1
                else:
                    kAltp1 = kAlt

                for i, v in enumerate(vars_in):
                    BeforeVals[i] = \
                        (1-xLon)*(1-yLat)*(1-zAlt)*BeforeData[v][iLon][jLat][kAlt]+\
                        (  xLon)*(1-yLat)*(1-zAlt)*BeforeData[v][iLon+1][jLat][kAlt]+\
                        (1-xLon)*(  yLat)*(1-zAlt)*BeforeData[v][iLon][jLat+1][kAlt]+\
                        (  xLon)*(  yLat)*(1-zAlt)*BeforeData[v][iLon+1][jLat+1][kAlt]+\
                        (1-xLon)*(1-yLat)*(  zAlt)*BeforeData[v][iLon][jLat][kAltp1]+\
                        (  xLon)*(1-yLat)*(  zAlt)*BeforeData[v][iLon+1][jLat][kAltp1]+\
                        (1-xLon)*(  yLat)*(  zAlt)*BeforeData[v][iLon][jLat+1][kAltp1]+\
                        (  xLon)*(  yLat)*(  zAlt)*BeforeData[v][iLon+1][jLat+1][kAltp1]
                    AfterVals[i] = \
                        (1-xLon)*(1-yLat)*(1-zAlt)*AfterData[v][iLon][jLat][kAlt]+\
                        (  xLon)*(1-yLat)*(1-zAlt)*AfterData[v][iLon+1][jLat][kAlt]+\
                        (1-xLon)*(  yLat)*(1-zAlt)*AfterData[v][iLon][jLat+1][kAlt]+\
                        (  xLon)*(  yLat)*(1-zAlt)*AfterData[v][iLon+1][jLat+1][kAlt]+\
                        (1-xLon)*(1-yLat)*(  zAlt)*AfterData[v][iLon][jLat][kAltp1]+\
                        (  xLon)*(1-yLat)*(  zAlt)*AfterData[v][iLon+1][jLat][kAltp1]+\
                        (1-xLon)*(  yLat)*(  zAlt)*AfterData[v][iLon][jLat+1][kAltp1]+\
                        (  xLon)*(  yLat)*(  zAlt)*AfterData[v][iLon+1][jLat+1][kAltp1]
            
                vals = []
                for i, v in enumerate(vars_in):
                    vals.append( (1-xt) * BeforeVals[i] + xt * AfterVals[i] )

                errs = random.normal(loc=0, scale=verr, size=(3))
                    
                ve = vals[0] + errs[0]
                vn = vals[1] + errs[1]
                vr = vals[2] + errs[2]

                rNew = R + vr * dt
                altNew = (rNew - Re) / 1000.0
                rAve = (R + rNew) / 2
            
                latNew = lat + vn * dt / rAve / cDtoR
                doWrap = False
                if (latNew > 90.0):
                    latNew = 180.0 - latNew
                    doWrap = True
                if (latNew < -90.0):
                    latNew = -(180.0 + latNew)
                    doWrap = True
                latAve = (lat + latNew)/2.0
                cosLat = np.cos(latAve * dtor)
                lonNew = (lon + ve * dt / cosLat / rAve / cDtoR + 360.0) % 360.0
                if (doWrap):
                    lonNew = (lonNew + 180.0) % 360.0

                lonsOut[iT][iPt] = lonNew
                latsOut[iT][iPt] = latNew
                altsOut[iT][iPt] = altNew
            
        iBefore0 = iBefore
        iAfter0 = iAfter

    return lonsOut, latsOut, altsOut


# ------------------------------------------------------------------------------------
# Main code
# ------------------------------------------------------------------------------------

headers = read_gitm_headers('3DALL')

totalDt = (headers["time"][-1] - headers["time"][0]).total_seconds()

dt = 10.0
tSeconds = np.arange(0, totalDt, dt)
times = []
for t in tSeconds:
    times.append(headers["time"][0] + timedelta(seconds = t))

nPts = 3    
lats = random.rand(nPts) * 180.0 - 90.0
lons = random.rand(nPts) * 360.0
alts = np.zeros(nPts) + 350.0

lonsOut, latsOut, altsOut = thermo_get_points(lats, lons, alts, times, verr = 0.0)

fig = plt.figure(figsize=(10, 8.5))
ax = fig.add_subplot(111)
nPts = len(lats)
for iPt in range(nPts):
    ax.plot(lonsOut[:, iPt], latsOut[:, iPt])

ax.set_ylim(-90,90)
ax.set_xlim(0,360)

fileout = 'test.png'
print('Writing file : ' + fileout)
fig.savefig(fileout)
