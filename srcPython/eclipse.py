#!/usr/bin/env python

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import re
import argparse
import datetime as dt

# ----------------------------------------------------------------------
# 
# ----------------------------------------------------------------------

months = {'Jan' : 1,
          'Feb' : 2,
          'Mar' : 3,
          'Apr' : 4,
          'May' : 5,
          'Jun' : 6,
          'Jul' : 7,
          'Aug' : 8,
          'Sep' : 9,
          'Oct' : 10,
          'Nov' : 11,
          'Dec' : 12}

def read_file(file):

    fpin = open(file, 'r')

    lines = fpin.readlines()

    iStart = 0
    iEnd = len(lines)
    iLine = iStart

    inBetween = False

    lats = []
    lons = []
    times = []
    
    while (iLine < iEnd):
        line = lines[iLine]

        m = re.match(r'.*Path of Total Solar Eclipse of (\d\d\d\d) (...) (\d\d).*', line)
        if m:
            year = int(m.group(1))
            month = months[m.group(2)]
            day = int(m.group(3))
        
        m = re.match(r'.*>Limits<',line)
        if m:
            if (inBetween):
                iLine = iEnd
            else:
                inBetween = True
        else:

            if (inBetween):
                parts = line.split('</td><td>')
                # figure out time:
                hour = int(parts[0][-5:-3])
                minute = int(parts[0][-2:])
                times.append(dt.datetime(year, month, day, hour, minute, 0))
                # Grab Latitudes
                lats_deg = float(parts[5][0:2])
                lats_min = float(parts[5][3:7]) / 60.0
                if (parts[5][-1] == 'S'):
                    fac = -1.0
                else:
                    fac = 1.0
                lats.append(fac * (lats_deg + lats_min))
                # Grab Longitudes
                lons_deg = float(parts[6][0:3])
                lons_min = float(parts[6][4:8]) / 60.0
                if (parts[6][-1] == 'W'):
                    lon = 360.0 - (lons_deg + lons_min)
                else:
                    lon = lons_deg + lons_min
                lons.append(lon)
        iLine += 1

    data = {
        'lats': np.array(lats),
        'lons': np.array(lons),
        'times': times}

    return data
        
# ----------------------------------------------------------------------
# 
# ----------------------------------------------------------------------

def calc_declinations(times):
    declinations = []
    maxDec = 23.5
    year = times[0].year
    equinox = dt.datetime(year, 3, 21, 0, 0, 0)
    sInYear = 365.25 * 24.0 * 3600.0
    for t in times:
        rads = (t - equinox).total_seconds() / sInYear * np.pi * 2
        declinations.append(maxDec * np.sin(rads))
    return np.array(declinations)

# ----------------------------------------------------------------------
# 
# ----------------------------------------------------------------------

def calc_localtimes(lons, times):
    localtimes = []
    for i, t in enumerate(times):
        ut = t.hour + t.minute/60.0 + t.second/3600.0
        localtimes.append((lons[i]/15.0 + ut) % 24.0)
    return np.array(localtimes)


def calc_gse(latitudes, localtimes, declinations):

    r = 6371.0
    lts = localtimes * np.pi / 12.0
    decs = declinations * np.pi / 180.0
    lats = latitudes * np.pi / 180.0
    x = r * np.cos(lts - np.pi) * np.cos(lats)
    y = r * np.sin(lts - np.pi) * np.cos(lats)
    z = r * np.sin(lats)
    # rotate around y-axis for the Earth Tilt:
    xp = x * np.cos(-decs) - z * np.sin(-decs)
    yp = y
    zp = x * np.sin(-decs) + z * np.cos(-decs)

    gse = {'x' : xp,
           'y' : yp,
           'z' : zp}
    return gse

# ----------------------------------------------------------------------
# 
# ----------------------------------------------------------------------

def extend(var, n, dir):
    if (dir > 0):
        dv = var[-1] - var[-2]
        new = var[-1] + dv * n
    else:
        dv = var[0] - var[1]
        new = var[0] + dv * n
    return new

# ----------------------------------------------------------------------
# 
# ----------------------------------------------------------------------

def extend_time(var, n, dir):
    if (dir > 0):
        dv = (var[-1] - var[-2]).total_seconds()
        new = var[-1] + dt.timedelta(seconds = dv * n)
    else:
        dv = (var[0] - var[1]).total_seconds()
        new = var[0] + dt.timedelta(seconds = dv * n)
    return new

# -------------------------------------------------------------------
# Make a label that shows the whole time range
# -------------------------------------------------------------------

def get_label_time_range(times):
    startTime = times[0]
    endTime = times[-1]    
    label = startTime.strftime('%b %d, %Y %H:%M UT') + ' - ' + \
        endTime.strftime('%b %d, %Y %H:%M UT (Hours)')
    return label
# ----------------------------------------------------------------------
# 
# ----------------------------------------------------------------------

# Code tested with this website:
# https://eclipsewise.com/solar/SEpath/2001-2100/SE2024Apr08Tpath.html

file = 'eclipse_20240408.html'
data = read_file(file)
decs = calc_declinations(data['times'])
lts = calc_localtimes(data['lons'], data['times'])

gse = calc_gse(data['lats'], lts, decs)

fig = plt.figure(figsize = (10,9))
plt.rcParams.update({'font.size': 14})

ax = fig.add_axes([0.12, 0.08, 0.85, 0.85])

ax.scatter(gse['y'], gse['z'], color = 'blue')
nP = 20
yP = extend(gse['y'], nP, 1.0)
zP = extend(gse['z'], nP, 1.0)
yM = extend(gse['y'], nP, -1.0)
zM = extend(gse['z'], nP, -1.0)
tP = extend_time(data['times'], nP, 1.0)
tM = extend_time(data['times'], nP, -1.0)
ax.scatter([yM, yP], [zM, zP], color = 'red')
ax.set_ylim(-9000, 9000)
ax.set_xlim(-9000, 9000)

sTime = get_label_time_range([tM, tP])
ax.set_title('Eclipse from : ' + sTime)
ax.set_xlabel('GSE-Y (km)')
ax.set_ylabel('GSE-Z (km)')

theta = np.arange(0, 361.0, 1.0) * np.pi / 180.0
r = 6371.0
ax.plot(r * np.cos(theta), r * np.sin(theta), color = 'k')
ax.set_aspect(1)

ax.text(yM+100, zM, 'Start', color = 'red')
ax.text(yP+100, zP, 'End', color = 'red')


print('#ECLIPSE')
print(tM.strftime('%Y %m %d %H %M %S           Start Time'))
print(tP.strftime('%Y %m %d %H %M %S           End Time'))
print('%7.1f       Start GSE-Y' % yM)
print('%7.1f       Start GSE-Z' % zM)
print('%7.1f       End GSE-Y' % yP)
print('%7.1f       End GSE-Z' % zP)

peak = 0.9
width = 3000.0
minExp = 0.07
expWidth = 400.0

print('  %5.3f       Peak Depletion' % peak)
print('%7.1f       Width' % width)
print('  %5.3f       Min Exp Amplitude' % minExp)
print('%7.1f       Exp width' % expWidth)

plotfile = 'eclipse.png'
print('writing : ',plotfile)    
fig.savefig(plotfile)
plt.close()
