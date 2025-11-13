#!/usr/bin/env python

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import re
import argparse
import datetime as dt

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def get_args_timeline():

    parser = argparse.ArgumentParser(description =
                                     'Post process and move model results')
    parser.add_argument('-plotfile',
                        help = 'output file for plot',
                        default = 'timeline.png')

    parser.add_argument('-vars',
                        help = 'var(s) to plot (e.g. -vars 0 or -vars 0 1 2',
                        default = [0], nargs = '+', type = int)
    
    parser.add_argument('filelist', nargs='+', \
                        help = 'list files to use for generating plots')
    
    parser.add_argument('-gitm',  \
                        action='store_true', default = False, \
                        help = 'GITM log file')
    
    parser.add_argument('-start', default = '0000', \
                        help = 'start date as YYYYMMDD[.HHMM]')
    parser.add_argument('-end', default = '0000', \
                        help = 'end date as YYYYMMDD[.HHMM]')
    
    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------
# convert a yyyymmdd.hhmm string to date time
# ----------------------------------------------------------------------

def convert_time(sTime):
    
    yTime = sTime[0:4]
    mo = sTime[4:6]
    da = sTime[6:8]

    if (len(sTime) >= 11):
        hr = sTime[9:11]
        if (len(sTime) >= 13):
            mi = sTime[11:13]
        else:
            mi = '00'
    else:
        hr = '00'
        mi = '00'

    daTime = dt.datetime(int(yTime), int(mo), int(da), int(hr), int(mi), 0)
    diTime =  {
        0: [yTime, ' year'],
        1: [mo, ' month'],
        2: [da, ' day'],
        3: [hr, ' hour'],
        4: [mi, ' minute'],
        5: ['00', ' second']
    }
    return daTime, diTime

# ----------------------------------------------------------------------
# Find closest time, given an array of datetimes and a datetime to find
# ----------------------------------------------------------------------

def find_closest_time(timeArray, timeToFind):
    deltas = []
    for time in timeArray:
        deltas.append((timeToFind - time).total_seconds())
    iTime_ = np.argmin(np.abs(deltas))
    return iTime_

# ----------------------------------------------------------------------
# Read file
# ----------------------------------------------------------------------

def read_timeline_file(file):

    data = {"times" : [],
            'vars': [],
            "integral" : 'values',
            "file": file}

    fpin = open(file, 'r')

    lines = fpin.readlines()

    iStart = 0
    iEnd = len(lines)
    iLine = iStart

    while (iLine < iEnd):
        line = lines[iLine]

        m = re.match(r'#VAR',line)
        if m:
            # skip year, month, day, hour, minute, second:
            if (lines[iLine+1].strip() == 'Year'):
                iLine += 7
            else:
                iLine += 1
            while (len(lines[iLine].strip()) > 1):
                data["vars"].append(lines[iLine].strip())
                iLine += 1
            data['var'] = data["vars"][0]

        m = re.match(r'#INTEGRAL',line)
        if m:
            iLine += 1
            data["integral"] = lines[iLine].strip()

        m = re.match(r'#ALTITUDE',line)
        if m:
            iLine += 1
            data["alt"] = float(lines[iLine].strip())

        m = re.match(r'#DIRECTORY',line)
        if m:
            iLine += 1
            data["dir"] = lines[iLine].strip()

        m = re.match(r'#START',line)
        if m:
            iStart = iLine + 1
            break

        iLine += 1

    nVars = len(data['vars'])
    nTimes = iEnd - iStart
    allValues = np.zeros([nVars, nTimes])
    for i in range(iStart, iEnd):
        aline = lines[i].split()
        year = int(aline[0])
        month = int(aline[1])
        day = int(aline[2])
        hour = int(aline[3])
        minute = int(aline[4])
        second = int(aline[5])
        t = dt.datetime(year, month, day, hour, minute, second)
        
        data["times"].append(t)
        for iVar in range(nVars):
            allValues[iVar, i - iStart] = float(aline[6 + iVar])

    data['values'] = allValues
            
    return data

# ----------------------------------------------------------------------
# Read GITM Run log file
# GITM log files have a specific (old) format with a header that can
# basically be ignored, and you just have to look for the #START
# The variables are the line after that and then the data starts.
# time lines include:
# iStep yyyy mm dd hh mm ss  ms data1 data2 ....
# ----------------------------------------------------------------------

def read_gitm_log_file(file):

    data = {"times" : [],
            "vars": [],
            "integral" : 'values',
            "file": file,
            "alt": 0}

    fpin = open(file, 'r')

    lines = fpin.readlines()

    iStart = 0
    iEnd = len(lines)
    iLine = iStart

    while (iLine < iEnd):
        line = lines[iLine]

        m = re.match(r'#START',line)
        if m:
            iLine += 1
            aline = lines[iLine].split()
            data["vars"] = aline[8:]
            iStart = iLine + 1
            break

        iLine += 1

    nVars = len(data['vars'])
    nTimes = iEnd - iStart
    allValues = np.zeros([nVars, nTimes])
    for i in range(iStart, iEnd):
        aline = lines[i].split()
        year = int(aline[1])
        month = int(aline[2])
        day = int(aline[3])
        hour = int(aline[4])
        minute = int(aline[5])
        second = int(aline[6])
        t = dt.datetime(year, month, day, hour, minute, second)
        
        data["times"].append(t)
        for iVar in range(nVars):
            allValues[iVar, i - iStart] = float(aline[8 + iVar])

    data['values'] = allValues
            
    return data



# ----------------------------------------------------------------------
# match filename with colors and linestyles
# ----------------------------------------------------------------------

def assign_var_to_color(var):

    color = 'black'
    line = 'solid'
    label = file
    
    m = re.match('.*GOCE.*', var)
    if m:
        color = 'black'
        line = 'solid'
        label = 'GOCE'

    m = re.match('.*GITM.*', var)
    if m:
        color = 'blue'
        line = 'solid'
        label = 'GITM'
        
    m = re.match('.*Smooth.*', var)
    if m:
        line = 'dashed'
        label = label + ' (smoothed)'
        
    return color, line, label

# ----------------------------------------------------------------------
# match filename with colors and linestyles
# ----------------------------------------------------------------------

def assign_file_to_color(file):

    color = 'black'
    line = 'solid'
    label = file
    
    m = re.match(r'base', file)
    if m:
        color = 'black'
        line = 'solid'
        label = 'Idealized Model'

    m = re.match(r'fre', file)
    if m:
        color = 'grey'
        line = 'solid'
        label = 'FRE Model'

    m = re.match(r'ovation', file)
    if m:
        color = 'blue'
        line = 'dashed'
        label = 'Ovation'

    m = re.match(r'avee_050', file)
    if m:
        color = 'darkblue'
        line = 'dashed'
    m = re.match(r'avee_075', file)
    if m:
        color = 'blue'
        line = 'dashed'

    m = re.match(r'avee_150', file)
    if m:
        color = 'firebrick'
        line = 'dashdot'
        
    m = re.match(r'avee_200', file)
    if m:
        color = 'red'
        line = 'dashdot'
        
    m = re.match(r'fta', file)
    if m:
        color = 'forestgreen'
        line = 'dashed'
        label = 'FTA Model'

    m = re.match(r'swmf', file)
    if m:
        color = 'blue'
        line = 'dashdot'
        label = 'SWMF Driven'

    m = re.match(r'ae', file)
    if m:
        color = 'plum'
        line = 'dashed'
        label = 'FRE, AE Power'
        
    m = re.match(r'ions', file)
    if m:
        color = 'plum'
        line = 'dashdot'
        label = 'Ions'
        
    m = re.match(r'.*no', file)
    if m:
        color = 'cyan'
        label = 'NO Cooling'
    m = re.match(r'.*ch', file)
    if m:
        color = 'red'
        label = 'Chemical Heating'
    m = re.match(r'.*ht', file)
    if m:
        color = 'blue'
        label = 'Heat Transfer'
    m = re.match(r'.*jh', file)
    if m:
        label = 'Joule Heating'
        
    return color, line, label

# Needed to run main script as the default executable from the command line
if __name__ == '__main__':

    args = get_args_timeline()

    iVar = args.vars[0]

    nFiles = len(args.filelist)
    nVars = len(args.vars)

    useStartTime = False
    useEndTime = False
    if (args.start != '0000'):
        startTime, dummy = convert_time(args.start)
        useStartTime = True
        print('  -> startTime set to : ', startTime)
    if (args.end != '0000'):
        endTime, dummy = convert_time(args.end)
        useEndTime = True
        print('  -> endTime set to : ', endTime)

    fig = plt.figure(figsize = (10,10))
    ax = fig.add_subplot(111)
    
    for file in args.filelist:
        print("Reading file : ", file)
        if (args.gitm):
            data = read_gitm_log_file(file)
        else:
            data = read_timeline_file(file)
        if (not useStartTime):
            if (file == args.filelist[0]):
                startTime = data["times"][0]
            else:
                if (data["times"][0] < startTime):
                    startTime = data["times"][0]
        if (not useEndTime):
            if (file == args.filelist[0]):
                endTime = data["times"][-1]
            else:
                if (data["times"][-1] > endTime):
                    endTime = data["times"][-1]
                        
        for iVar in args.vars:
            if (nFiles > 1):
                color, line, label = assign_file_to_color(file.lower())
                if (nVars > 1):
                    label = data['vars'][iVar] + ' (' + label + ')'
                    if (iVar == args.vars[0]):
                        line = 'solid'
                    #else:
                    #    line = 'dashed'
            else:
                fileColor, fileLine, fileLabel = assign_file_to_color(file.lower())
                color, line, label = assign_var_to_color(data['vars'][iVar])
                print(fileLabel, ' -> ', label)
                if (label == 'GITM'):
                    label = label + '(' + fileLabel + ')'
                    print('Label Changed : ', label)
            print(label)
            ax.plot(data["times"], data["values"][iVar], label = label,
                    color = color, linestyle = line, linewidth = 2.0)

    iVar = args.vars[0]
    ytitle = data["integral"] + " of " + \
        data["vars"][iVar] + " at " + \
        "%5.1f" % data["alt"] + " km Alt"
    ax.set_ylabel(ytitle)

    ax.legend()
    ax.set_xlim(startTime, endTime)
    sTime = startTime.strftime('%b %d, %Y %H:%M UT') + ' - ' + \
        endTime.strftime('%b %d, %Y %H:%M UT (Hours)')
    ax.set_xlabel(sTime)
    
    if (nFiles == 1):
        iStart_ = find_closest_time(data["times"], startTime)   
        iEnd_ = find_closest_time(data["times"], endTime)

        print('Calculating Statistics between start and end indices :', iStart_, iEnd_)
        print('  --> Start Time : ', data["times"][iStart_])
        print('  --> End Time : ', data["times"][iEnd_])
        iRef = args.vars[0]
        var0 = data['vars'][iRef]
        for iVar in args.vars[1:]:
            var1 = data['vars'][iVar]
            print(' -> Comparing variables : ', var0, ' (ref) to ', var1)

            vals0 = data["values"][iRef][iStart_ : iEnd_]
            vals1 = data["values"][iVar][iStart_ : iEnd_]
            diff = vals1 - vals0
            
            ratio_ave = np.mean(data["values"][iVar][iStart_ : iEnd_]) / \
                np.mean(data["values"][iRef][iStart_ : iEnd_])
            min1 = np.min(data["values"][iVar][iStart_ : iEnd_]) 
            max1 = np.max(data["values"][iVar][iStart_ : iEnd_]) 
            min0 = np.min(data["values"][iRef][iStart_ : iEnd_]) 
            max0 = np.max(data["values"][iRef][iStart_ : iEnd_])
            range1 = (max1 - min1)
            range0 = (max0 - min0)
            ratio_max = max1 / max0
            ratio_range = range1 / range0

            mean_diff = np.mean(vals1 - vals0)
            nrmse = np.sqrt(np.mean(diff**2)) / (max0 - min0)
            pe = 1.0 - np.sqrt(np.mean(diff**2)) / np.sqrt(np.mean(vals0**2))

            sRatioMeans = 'Ratio of means : %4.2f' % ratio_ave
            sRatioMaxs = 'Ratio of maxes : %4.2f' % ratio_max
            sRatioRanges = 'Ratio of ranges : %4.2f' % ratio_range
            sMeanDiff =  'Mean Difference : %10.4e' % mean_diff
            sNrms = 'Normalized RMSE : %4.2f' % nrmse
            sPE = 'Prediction Eff. : %5.2f' % pe
            print('  -> Ratio of means : ', ratio_ave)
            print('  -> Ratio of maxs : ', ratio_max)
            print('  -> Ratio of ranges : ', ratio_range)
            print('  -> Mean Difference : ', mean_diff)
            print('  -> nrmse : ', nrmse)
            print('  -> pe : ', pe)

            if (len(args.vars) == 2):
                mini = np.min([min1, min0])
                maxi = np.max([max1, max0])
                r = maxi - mini
                mini = mini - 0.15 * r
                maxi = maxi + 0.02 * r
                deltaT = (data["times"][iEnd_] - data["times"][iStart_]).total_seconds()
                t10 = data["times"][iStart_] + dt.timedelta(seconds = deltaT * 0.1)
                t40 = data["times"][iStart_] + dt.timedelta(seconds = deltaT * 0.4)
                ax.set_ylim(mini, maxi)
                ax.text(t10, mini + 0.01 * r, sRatioMeans)
                ax.text(t10, mini + 0.05 * r, sRatioMaxs)
                ax.text(t10, mini + 0.09 * r, sRatioRanges)
                ax.text(t40, mini + 0.01 * r, sMeanDiff)
                ax.text(t40, mini + 0.05 * r, sNrms)
                ax.text(t40, mini + 0.09 * r, sPE)
                
    plotfile = args.plotfile
    print('writing : ',plotfile)    
    fig.savefig(plotfile)
    plt.close()

