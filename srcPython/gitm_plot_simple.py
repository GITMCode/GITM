#!/usr/bin/env python3
""" Standard model visualization routines
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import argparse

from pylab import cm
from gitm_routines import *

def get_args():

    parser = argparse.ArgumentParser(
        description = 'Plot Aether / GITM model results')
    
    parser.add_argument('-list',  \
                        action='store_true', default = False, \
                        help = 'list variables in file')

    parser.add_argument('-var',  \
                        default = 3, type = int, \
                        help = 'variable to plot (number)')

    parser.add_argument('-cut', metavar = 'cut',  default ='alt', \
                        choices = ['alt', 'lat', 'lon'], 
                        help = 'alt,lat,lon : which cut you would like')

    parser.add_argument('-alt', metavar = 'alt', default = 250.0, \
                        help = 'altitude :  alt in km (closest)')
    parser.add_argument('-lat', metavar = 'lat',  default = 0.0, \
                        help = 'latitude : latitude in degrees (closest)')
    parser.add_argument('-lon', metavar = 'lon',  default = 0.0,\
                        help = 'longitude in degrees (closest)')

    parser.add_argument('-log',  default = False,
                        action="store_true",
                        help = 'plot the log of the variable')
    
    parser.add_argument('-mini',  default = 1e32, type = float, \
                        help = 'manually set the minimum value for the plots')
    parser.add_argument('-maxi',  default = -1e32, type = float, \
                        help = 'manually set the maxiumum value for the plots')

    parser.add_argument('filelist', nargs='+', \
                        help = 'list files to use for generating plots')
    
    args = parser.parse_args()

    return args
    
# ----------------------------------------------------------------------------
# get header info from a file
# ----------------------------------------------------------------------------

def get_file_info(args):

    header = read_gitm_header(args.filelist)
    header['vars'] = remap_variable_names(header['vars'])
    
    return header

# ----------------------------------------------------------------------------
# list the variables in the header
# ----------------------------------------------------------------------------

def list_vars(header):
    for k, v in header.items():
        if (k != 'vars'):
            print(k, '-> ', v)
        else:
            print('vars : ')
            for i, var in enumerate(v):
                print(i, var)
    return

# ----------------------------------------------------------------------------
# Read in data from the files:
# ----------------------------------------------------------------------------

def read_in_model_files(filelist, varlist):

    # first read in spatial information:
    vars = [0, 1, 2]
    spatialData = read_gitm_one_file(filelist[0], vars)

    lons = np.degrees(spatialData[0])  # Convert from rad to deg
    nLons = len(lons[:, 0, 0])
    lats = np.degrees(spatialData[1])  # Convert from rad to deg
    nLats = len(lats[0, :, 0])
    alts = spatialData[2] / 1000.0  # Convert from m to km
    nAlts = len(alts[0, 0, :])
    
    nTimes = len(filelist)
    nVars = len(varlist)
    if (nVars == 1):
        allData = np.zeros((nTimes, nLons, nLats, nAlts))
    else:
        allData = np.zeros((nTimes, nVars, nLons, nLats, nAlts))
    allTimes = []
    for iTime, filename in enumerate(filelist):
        data = read_gitm_one_file(filename, varlist)
        allTimes.append(data["time"])
        if (nVars == 1):
            allData[iTime, :, :, :] = data[varlist[0]][:, :, :]
        else:
            for iVar, var in enumerate(varlist):
                allData[iTime, iVar, :, :, :] = data[var][:, :, :]
    vars = []
    for var in varlist:
        vars.append(data['vars'][var])

    data = {'times': allTimes,
            'data': allData,
            'vars': vars,
            'lons': lons,
            'lats': lats,
            'alts': alts,
            'nTimes': nTimes,
            'nVars': nVars,
            'nLons' : nLons,
            'nLats': nLats,
            'nAlts': nAlts}
    
    return data

# ----------------------------------------------------------------------------
# take a dictionary containing all of the model data and
# return slices. The data should have the shape:
# [nTimes, nVars, nLons, nLats, nAlts]
# or
# [nTimes, nLons, nLats, nAlts]
# ----------------------------------------------------------------------------

def data_slice(allData3D, iLon = -1, iLat = -1, iAlt = -1):

    nTimes = allData3D['nTimes']
    nVars = allData3D['nVars']
    nLons = allData3D['nLons']
    nLats = allData3D['nLats']
    nAlts = allData3D['nAlts']

    if (nVars > 1):
        if (iAlt > -1):
            slices = np.zeros((nTimes, nVars, nLons, nLats))
            slices[:, :, :, :] = allData3D['data'][:, :, :, :, iAlt]
        elif (iLat > -1):
            slices = np.zeros((nTimes, nVars, nLons, nAlts))
            slices[:, :, :, :] = allData3D['data'][:, :, :, iLat, :]
        else:
            slices = np.zeros((nTimes, nVars, nLats, nAlts))
            slices[:, :, :, :] = allData3D['data'][:, :, iLon, :, :]
    else:
        if (iAlt > -1):
            slices = np.zeros((nTimes, nLons, nLats))
            slices[:, :, :] = allData3D['data'][:, :, :, iAlt]
        elif (iLat > -1):
            slices = np.zeros((nTimes, nLons, nAlts))
            slices[:, :, :] = allData3D['data'][:, :, iLat, :]
        else:
            slices = np.zeros((nTimes, nLats, nAlts))
            slices[:, :, :] = allData3D['data'][:, iLon, :, :]

    return slices

#-----------------------------------------------------------------------------
# find which cut direction to make, then which cut to take
#-----------------------------------------------------------------------------

def find_cut(args, allData):

    cutValue = -1e32
    cutString = ''
    cutShort = ''
    iLon = -1
    iLat = -1
    iAlt = -1
    xLabel = ''
    yLabel = ''
    lons1d = allData['lons'][:, 0, 0]
    lats1d = allData['lats'][0, :, 0]
    alts1d = allData['alts'][0, 0, :]
    xRange = [0,0]
    yRange = [0,0]

    if (args.cut == 'alt'):
        altGoal = float(args.alt)
        diff = np.abs(alts1d - altGoal)
        iAlt = np.argmin(diff)
        cutValue = alts1d[iAlt]
        cutString = 'Alt : %d km' % int(cutValue)
        cutShort = 'alt%04d_' % iAlt
        xLabel = 'Longitude (deg)'
        yLabel = 'Latitude (deg)'
        xPos1d = lons1d
        yPos1d = lats1d
        xRange = [0, 360]
        yRange = [-90.0, 90.0]
    if (args.cut == 'lon'):
        lonGoal = float(args.lon)
        diff = np.abs(lons1d - lonGoal)
        iLon = np.argmin(diff)
        cutValue = lons1d[iLon]
        cutString = 'Lon : %d deg' % int(cutValue)
        cutShort = 'lon%04d_' % iLon
        xLabel = 'Latitude (deg)'
        yLabel = 'Altitude (km)'
        xPos1d = lats1d
        yPos1d = alts1d
        xRange = [-90.0, 90.0]
        yRange = [alts1d[2], alts1d[-3]]
    if (args.cut == 'lat'):
        latGoal = float(args.lat)
        diff = np.abs(lats1d - latGoal)
        iLat = np.argmin(diff)
        cutValue = lats1d[iLat]
        cutString = 'Lat : %d deg' % int(cutValue)
        cutShort = 'lat%04d_' % iLat
        xLabel = 'Longitude (deg)'
        yLabel = 'Altitude (km)'
        xPos1d = lons1d
        yPos1d = alts1d
        xRange = [0.0, 360.0]
        yRange = [alts1d[2], alts1d[-3]]

    cut = {'iLon': iLon,
           'iLat': iLat,
           'iAlt': iAlt,
           'cutValue': cutValue,
           'cutString': cutString,
           'cutShort': cutShort,
           'xLabel': xLabel,
           'yLabel': yLabel,
           'xRange': xRange,
           'yRange': yRange,
           'xPos': xPos1d,
           'yPos': yPos1d}

    return cut

# ----------------------------------------------------------------------------
# This function calculates the edges of cells based on the centers of the cells
# it assumes a 1D array.
# ----------------------------------------------------------------------------

def move_centers_to_edges(pos):
    edges = (pos[1:] + pos[:-1])/2
    dpLeft = pos[1] - pos[0]
    dpRight = pos[-1] - pos[-2]
    edges = np.append(edges[0] - dpLeft, edges)
    edges = np.append(edges, dpRight + edges[-1])
    return edges

# ----------------------------------------------------------------------------
# Get max and min values
# model_data is [time, xpositions, ypositions]
# ----------------------------------------------------------------------------

def get_min_max_data(allSlices, yPos, \
                     yMin = -1e32, yMax = 1e32, \
                     color = 'default', \
                     minVal = 1e32, maxVal = -1e32,
                     isLog = False):
    
    symmetric = False
    if (color == 'default'):
        cmap = mpl.cm.plasma
    if (color == 'red'):
        cmap = mpl.cm.YlOrRd

    mask = ((yPos >= yMin) & (yPos <= yMax))
    if (mask.max()):    
        doPlot = True
        maxi = allSlices[:, :, mask].max() * 1.01
        mini = allSlices[:, :, mask].min() * 0.99
    
        if ((mini < 0.0) and (not isLog)):
            symmetric = True
            cmap = mpl.cm.bwr
            maxi = abs(allSlices[:, :, mask]).max() * 1.05
            mini = -maxi

    else:
        doPlot = False
        mini = -1.0
        maxi = 1.0

    if (minVal < 1e31):
        mini = minVal
    if (maxVal > -1e31):
        maxi = maxVal
        
    min_max_data = {'mini' : mini,
                    'maxi' : maxi,
                    'cmap' : cmap,
                    'symmetric' : symmetric,
                    'doPlot': doPlot,
                    'mask': mask}
    
    return min_max_data

# ----------------------------------------------------------------------------
# Plot a series of slices
#   allSlices: 3D - [time, xpositions, ypositions]
#   xPosEdges: 1D - [xpositions] (edges of cells, not centers)
#   yPosEdges: 1D - [ypositions] (edges of cells, not centers)
#   dataMinMax: dictionary that contains info on min, max, and color map
#   varName: puts this string on the colorbar
#   titleAddOn: puts this string after the time on the title
#   xLabel: string label for the x axis
#   yLabel: string label for the y axis
#   filenamePrefix: string to add before time for filename (like varX_)
#   xLimits: limits for the x axis - (array of 2 elements)
#   yLimits: limits for the y axis - (array of 2 elements)
#   dpi: dots per inch for file
# ----------------------------------------------------------------------------

def plot_series_of_slices(allSlices,
                          allTimes,
                          xPosEdges,
                          yPosEdges,
                          dataMinMax,
                          varName = '',
                          titleAddOn = '',
                          xLabel = '',
                          yLabel = '',
                          filenamePrefix = '',
                          xLimits = [0, 0],
                          yLimits = [0, 0],
                          dpi = 120):

    for iTime, uTime in enumerate(allTimes):

        fig = plt.figure(figsize=(10, 5.5), dpi = dpi)
        ax = fig.add_axes([0.07, 0.08, 0.97, 0.87])

        title = uTime.strftime("%d %b %Y %H:%M:%S UT") + titleAddOn
        value2d = allSlices[iTime, :, :].transpose()
        con = ax.pcolormesh(xPosEdges, yPosEdges, value2d, \
                            cmap = dataMinMax['cmap'], \
                            vmin = dataMinMax['mini'], \
                            vmax = dataMinMax['maxi'])
        if (xLimits[1] > xLimits[0]):
            ax.set_xlim(xLimits)
        if (yLimits[1] > yLimits[0]):
            ax.set_ylim(yLimits)
        #ax.set_aspect(1.0)
        ax.set_title(title)
        ax.set_xlabel(xLabel)
        ax.set_ylabel(yLabel)

        cbar = fig.colorbar(con, ax = ax, shrink = 0.75, pad = 0.02)
        cbar.set_label(varName, rotation=90)

        sTimeOut = uTime.strftime('%y%m%d_%H%M%S')
        outFile = filenamePrefix + sTimeOut + '.png'
        print(" ==> Writing file : ", outFile)
        fig.savefig(outFile, dpi = dpi)
        plt.close(fig)

#-----------------------------------------------------------------------------
# Main code
#-----------------------------------------------------------------------------

if __name__ == '__main__':

    # Get the input arguments
    args = get_args()

    header = get_file_info(args)
        
    if (args.list):
        list_vars(header)
        exit()

    allData = read_in_model_files(args.filelist, [args.var])
    cut = find_cut(args, allData)
    
    allSlices = data_slice(allData, \
                           iLon = cut['iLon'], \
                           iLat = cut['iLat'], \
                           iAlt = cut['iAlt'])
    allTimes = allData['times']

    # get min and max values, plus color table:
    dataMinMax = get_min_max_data(allSlices, cut['yPos'], \
                                  color = 'red', \
                                  minVal = args.mini, maxVal = args.maxi)
    
    varName = allData['vars'][0]
    if (args.log):
        allSlices = np.log10(allSlices)
        varName = 'log10(' + varName + ')'
        dataMinMax['mini'] = np.log10(dataMinMax['mini'])
        dataMinMax['maxi'] = np.log10(dataMinMax['maxi'])
    sVarNum = 'var%03d_' % args.var
    sFilePre = sVarNum + cut['cutShort']
    sTitleAdd = '; ' + cut['cutString']
    xEdges = move_centers_to_edges(cut['xPos'])
    yEdges = move_centers_to_edges(cut['yPos'])
    plot_series_of_slices(allSlices,
                          allTimes,
                          xEdges,
                          yEdges,
                          dataMinMax,
                          xLabel = cut['xLabel'],
                          yLabel = cut['yLabel'],
                          varName = varName,
                          titleAddOn = sTitleAdd,
                          filenamePrefix = sFilePre,
                          xLimits = cut['xRange'],
                          yLimits = cut['yRange'])
