#!/usr/bin/env python
""" Standard model visualization routines
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import re
from glob import glob
import argparse
import sys

from aetherpy.io import read_routines
from aetherpy.utils import inputs, time_conversion
from useful_functions import *

from scipy.stats import ks_2samp

def get_args():

    parser = argparse.ArgumentParser(
        description = 'Plot Aether / GITM model histograms')
    
    parser.add_argument('-list',  \
                        action='store_true', default = False, \
                        help = 'list variables in file')
    parser.add_argument('-var',  \
                        default = 3, type = int, \
                        help = 'variable to plot (number)')
    parser.add_argument('-cut', metavar = 'cut',  default ='alt', \
                        choices = ['alt', 'lat', 'lon'], 
                        help = 'alt,lat,lon : which cut you would like')

    parser.add_argument('-alt', metavar = 'alt', default =400.0, type = int, \
                        help = 'altitude :  alt in km (closest)')
    parser.add_argument('-lat', metavar = 'lat',  default =-100.0, \
                        help = 'latitude : latitude in degrees (closest)')
    parser.add_argument('-lon', metavar = 'lon',  default =-100.0,\
                        help = 'longitude in degrees (closest)')
    
    parser.add_argument('-winds', default = False,\
                        help='overplot winds (doesnt do anything!!!)', \
                        action="store_true")
    parser.add_argument('-mean',  \
                        action='store_true', default = False, \
                        help = '(ignored!)')
    parser.add_argument('-tec',  \
                        action='store_true', default = False, \
                        help = '(ignored!)')
    
    parser.add_argument('-filelist', nargs='+', \
                        help = 'list files to use for generating plots')

    parser.add_argument('-list2', nargs='+', \
                        help = '2nd list files to use for generating plots')
    
    parser.add_argument('-label', default = 'Series 1', \
                        help = 'Legend label for series 1')
    parser.add_argument('-label2', default = 'Series 2', \
                        help = 'Legend label for series 2')
    
    args = parser.parse_args()

    return args

# Get the input arguments
args = get_args()

# first set of files:

header = get_file_info(args)
if (args.list):
    list_vars(header)
    exit()
if (args.var >= len(header["vars"])):
    raise ValueError("requested variable doesn't exist: {:d}>{:d}". \
                     format(args.var, len(header["vars"])))
data = read_in_model_files(args, header)

doCompare = False
if (args.list2 != None):
    filelist = args.filelist
    args.filelist = args.list2
    header2 = get_file_info(args)
    data2 = read_in_model_files(args, header2)
    doCompare = True

times = data['times']
slices = data['slices']
maxi = np.max(np.abs(slices)) * 1.01
alpha = 1.0
if (doCompare):
    slices2 = data2['slices']
    maxi = np.max([maxi, np.max(np.abs(slices2))])
    alpha = 0.5
nBins = 25
if (np.min(slices) < 0):
    mini = -maxi
    nBins = nBins * 2 + 1
else:
    mini = np.min(slices) * 0.99

hist, binEdges = np.histogram(slices, bins = nBins, range = (mini, maxi))
hist = hist / np.sum(hist) * 100.0
binCenters = (binEdges[0:-1] + binEdges[1:])/2

fig = plt.figure(figsize = (10,10))
plt.rcParams.update({'font.size': 14})
ax = fig.add_subplot(111)

db = binEdges[1] - binEdges[0]
ax.bar(binCenters, hist, width = db * 0.9, color = 'blue', \
       label = args.label, alpha = alpha)

if (doCompare):
    hist2, be = np.histogram(slices2, bins = nBins, range = (mini, maxi))
    hist2 = hist2 / np.sum(hist2) * 100.0
    ax.bar(binCenters, hist2, width = db * 0.9, color = 'red', \
           label = args.label2, alpha = alpha)
    # Performs the two-sample Kolmogorov-Smirnov test for goodness of fit.
    ks_test = ks_2samp(slices.ravel(), slices2.ravel())
    pvalue = ks_test.pvalue
    sCompare = ' (Dists. are '
    if (pvalue < 0.05):
        sCompare = sCompare + 'different: '
    else:
        sCompare = sCompare + 'same: '
    sCompare = sCompare + '%5.3f)' % pvalue
    print('pvalue = ', pvalue)
    ax.legend()
else:
    sCompare = ''

print(header['vars'][args.var])

iVar = args.var
xtitle = header["vars"][iVar] + " at " + \
    "%5.1f" % data["z_val"] + " km Alt" + sCompare
ax.set_xlabel(xtitle)

ax.set_ylabel('Occurance Rate (%)')

sTime = times[0].strftime('%b %d, %Y %H:%M UT') + ' - ' + \
    times[-1].strftime('%b %d, %Y %H:%M UT')
title = header["vars"][iVar] + ' from ' + sTime
ax.set_title(title)

    
plotfile = 'test.png'
print('writing : ',plotfile)    
fig.savefig(plotfile)
plt.close()
