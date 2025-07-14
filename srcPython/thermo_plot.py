#!/usr/bin/env python3
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
from aetherpy.plot import data_prep, movie_routines

from read_hiwinds import read_hiwind_file

def get_args():

    parser = argparse.ArgumentParser(
        description = 'Plot Aether / GITM model results')
    
    parser.add_argument('-list',  \
                        action='store_true', default = False, \
                        help = 'list variables in file')

    parser.add_argument('-timeplot',  \
                        action='store_true', default = False, \
                        help = 'Plot integrated (or mean) value vs. time')
    parser.add_argument('-mean',  \
                        action='store_true', default = False, \
                        help = 'Plot mean value instead of integrated value')
    parser.add_argument('-timefile',  default ='none', \
                        help = 'output filename for timeline file')
    
    parser.add_argument('-label',  \
                        action='store_true', default = False, \
                        help = 'Add label (e.g., (a), (b)..) to title')

    parser.add_argument('-color',  default ='default', \
                        choices = ['default', 'red'], 
                        help = 'set color bar')

    
    parser.add_argument('-var',  \
                        default = 3, type = int, \
                        help = 'variable to plot (number)')
    parser.add_argument('-cut', metavar = 'cut',  default ='alt', \
                        choices = ['alt', 'lat', 'lon'], 
                        help = 'alt,lat,lon : which cut you would like')
    parser.add_argument('-ext',  default ='png', \
                        choices = ['png', 'jpg', 'pdf'], 
                        help = 'plot type file extention')

    parser.add_argument('-hiwind',  default ='none', \
                        help = 'HIWIND file to plot location')
    
    parser.add_argument('-winds', default = False,\
                        help='overplot winds', \
                        action="store_true")
    parser.add_argument('-nstep', metavar = 'nstep',
                        default =4, type = int, \
                        help = 'number of steps between wind quivers')
    
    parser.add_argument('-nopole', default = False,\
                        help='dont plot polar regions', \
                        action="store_true")
    parser.add_argument('-north', default = False,\
                        help='only plot northern hemisphere results', \
                        action="store_true")
    parser.add_argument('-south', default = False,\
                        help='only plot southern hemisphere results', \
                        action="store_true")
    
    parser.add_argument('-alt', metavar = 'alt', default =400.0, type = int, \
                        help = 'altitude :  alt in km (closest)')
    parser.add_argument('-lat', metavar = 'lat',  default =-100.0, \
                        help = 'latitude : latitude in degrees (closest)')
    parser.add_argument('-lon', metavar = 'lon',  default =-100.0,\
                        help = 'longitude in degrees (closest)')
    parser.add_argument('-alog',  default = False,
                        action="store_true",
                        help = 'plot the log of the variable')
    parser.add_argument('-IsLog', default =False,
                        help='plot the log of the variable', 
                        action="store_true")    

    parser.add_argument('-mini',  default = 1e32, type = float, \
                        help = 'manually set the minimum value for the plots')
    parser.add_argument('-maxi',  default = -1e32, type = float, \
                        help = 'manually set the maxiumum value for the plots')

    parser.add_argument('-percent', default = False, 
                        action = 'store_true',
                        help = 'plot percentage difference of files')
    parser.add_argument('-diff', default = False, 
                        action = 'store_true',
                        help = 'plot difference of files (2 files needed)')
    parser.add_argument('-backdir', metavar = 'backdir',
                        default = 'none',
                        help = 'Subtract files in this directory')

    parser.add_argument('-mkv',
                        action = 'store_true',
                        default = False,
                        help = 'movie format = mkv')
    parser.add_argument('-mp4',
                        action = 'store_true',
                        default = True,
                        help = 'movie format = mp4')
    parser.add_argument('-gif',
                        action='store_true',
                        default = False, 
                        help = 'movie format = gif')
    parser.add_argument('-movie',  default = False,\
                        action='store_true',
                        help = 'Make a movie out of results')
    parser.add_argument('-tec',  default = False, \
                        action='store_true',
                        help = 'plot total electron content (TEC)')
    parser.add_argument('-on2',  default = False, \
                        action='store_true',
                        help = 'plot O/N2 ratio')
    parser.add_argument('-rate',  default =30,\
                        help = 'framerate for movie')
    parser.add_argument('filelist', nargs='+', \
                        help = 'list files to use for generating plots')
    
    args = parser.parse_args()

    return args


def determine_file_type(file):

    IsGitm = False
    HasHeader = False
    m = re.match(r'(.*)bin', file)
    if m:
        IsGitm = True
        # check for a header file:
        checkFile = glob(m.group(1)+"header")
        if (len(checkFile) > 0):
            if (len(checkFile[0]) > 1):
                HasHeader = True

    return IsGitm, HasHeader

#-----------------------------------------------------------------------------
# GITM's output variables are weird.  Fix them.
#-----------------------------------------------------------------------------

def fix_vars(vars):
    newvars = []
    for v in vars:
        nv = re.sub('!U', '', v)
        nv = re.sub('!N', '', nv)
        nv = re.sub('!D', '', nv)
        newvars.append(nv)

    return newvars


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


# ----------------------------------------------------------------------------
# Read in all of the model files:
# ----------------------------------------------------------------------------

def read_in_model_files(args, header):

    if (args.tec):
        ivar = 34
        plot_vars = [0, 1, 2, ivar]
        args.var = ivar
    elif (args.on2):
        iO_ = 4
        iN2_ = 6
        ivar = iO_
        plot_vars = [0, 1, 2, iO_, iN2_]
    else:
        plot_vars = [0, 1, 2, args.var]
        ivar = args.var

    # Update plotting variables to include the wind, if desired
    if args.winds:
        if (ivar > 18):
            plot_vars.append(37 if args.cut in ['alt', 'lat'] else 38)
            plot_vars.append(39 if args.cut in ['lat', 'lon'] else 38)
        else:
            plot_vars.append(16 if args.cut in ['alt', 'lat'] else 17)
            plot_vars.append(18 if args.cut in ['lat', 'lon'] else 17)

    all_winds_x = []
    all_winds_y = []

    # Prepare to load the desired file data
    all_2dim_data = []
    all_times = []
    all_int_data = []
    
    for j, filename in enumerate(header['filename']):
        # Read in the data file
        if header['IsGitm']:
            print('=> Reading GITM file : ', filename)
            data = read_routines.read_gitm_file(filename, plot_vars)
            ivar = args.var
        else:
            if j == 0:
                var_list = []
                for pvar in plot_vars:
                    var_list.append(header["vars"][pvar])
            if (header["HasHeader"]):
                data = read_routines.read_aether_one_binary_file(header, j,
                                                                 plot_vars)
                ivar = args.var
            else:
                data = read_routines.read_aether_file(filename, var_list)
                ivar = 3

        # For the first file, initialize the necessary plotting data
        if j == 0:
            # Get 1D arrays for the coordinates
            alts = data[2][0][0] / 1000.0  # Convert from m to km
            nAlts = len(alts)
            Alts3d = data[2]
            lons = np.degrees(data[0][:, 0, 0])  # Convert from rad to deg
            nLons = len(lons)
            lats = np.degrees(data[1][0, :, 0])  # Convert from rad to deg
            nLats = len(lats)
            # Find the desired index to cut along to get a 2D slice
            isgrid = False

            if (args.cut == 'lon'):
                pos = float(args.lon)
            if (args.cut == 'lat'):
                pos = float(args.lat)
            if (args.cut == 'alt'):
                pos = args.alt
                if (len(alts) == 1):
                    print("Only one alt found, setting alt pos = 0");
                    pos = 0
                    isgrid = True
                    print(alts)
                lat2d = data[1][:, :, 0]  # Convert from rad to deg
                dlon = data[0][1, 0, 0] - data[0][0, 0, 0]
                dlat = data[1][0, 1, 0] - data[1][0, 0, 0]
                area = np.cos(lat2d) * dlon * dlat * \
                    ((6372.0 + 100.0)*1000.0)**2
            icut, cut_data, x_pos, y_pos, z_val = data_prep.get_cut_index(
                lons, lats, alts, pos, isgrid, args.cut)
                
        if (args.cut == 'alt'):
            int_data = data[ivar][cut_data] * area
            if (args.mean):
                int_data = int_data / np.sum(area)
            all_int_data.append(np.sum(int_data))

        # Save the time data
        all_times.append(data["time"])

        # Save the z-axis data
        if args.tec:
            all_2dim_data.append(data_prep.calc_tec(alts, data[ivar], 2, -4))
        else:
            if (args.on2):
                oDensity = data[iO_]
                n2Density = data[iN2_]
                # 1e21/m2 is the integration boundary.  
                n2Int = vertically_integrate(n2Density, Alts3d, calc3D = True)
                oInt = vertically_integrate(oDensity, Alts3d, calc3D = True)
                on2 = np.zeros((nLons, nLats))
                iAlts = np.arange(nAlts)
                for iLat in range(nLats):
                    for iLon in range(nLons):
                        n21d = n2Int[iLon,iLat,:]/1e21
                        o1d = oInt[iLon,iLat,:]
                        i = iAlts[n21d < 1.0][0]
                        r = (1.0 - n21d[i]) / (n21d[i-1] - n21d[i])
                        n2 = (r * n21d[i-1] + (1.0 - r) * n21d[i]) * 1e21
                        o = r * o1d[i-1] + (1.0 - r) * o1d[i]
                        on2[iLon, iLat] = o / n2
                all_2dim_data.append(on2)
            else:
                all_2dim_data.append(data[ivar][cut_data])

        if (args.winds):
            all_winds_x.append(data[plot_vars[-2]][cut_data])
            all_winds_y.append(data[plot_vars[-1]][cut_data])

    # Convert data list to a numpy array
    all_2dim_data = np.array(all_2dim_data)
    
    if args.winds:
        all_winds_x = np.array(all_winds_x)
        all_winds_y = np.array(all_winds_y)

    data = {'winds_x' : np.array(all_winds_x),
            'winds_y' : np.array(all_winds_y),
            'times' : all_times,
            'slices' : np.array(all_2dim_data),
            'integrated_data': np.array(all_int_data),
            'icut' : icut,
            'x_pos' : x_pos,
            'y_pos' : y_pos,
            'z_val' : z_val}
    
    return data

# ----------------------------------------------------------------------------
# get header info from a file
# ----------------------------------------------------------------------------

def get_file_info(args):

    # determine what kind of files we are dealing with
    IsGitm, HasHeader = determine_file_type(args.filelist[0])
    
    if ((IsGitm) and (not HasHeader)):
        header = read_routines.read_gitm_headers(args.filelist, finds = 0)
    else:
        if (HasHeader):
            header = read_routines.read_aether_ascii_header(args.filelist)
            IsGitm = 0
        else:
            header = read_routines.read_aether_header(args.filelist)

    header['vars'] = fix_vars(header['vars'])
    header['IsGitm'] = IsGitm
    header['HasHeader'] = HasHeader
    
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
# Rescale the data if it is too large to too small
#   Want to eliminate axes with exponents in them
# ----------------------------------------------------------------------------

def rescale_data(original_data):
    factorString = ''
    new_data = original_data
    if (np.min(original_data) > 0.0):
        maxi = np.max(original_data)
        if ((np.log10(maxi) > 5.0) or (np.log10(maxi) < -5.0)):
            factor = 10**float(int(np.log10(maxi)))
            new_data = original_data / factor
            factorString = '(x%7.1e)' % factor
    
    return new_data, factorString


# ----------------------------------------------------------------------------
# Get max and min values
# model_data is [time, xpositions, ypositions]
# ----------------------------------------------------------------------------

def get_min_max_data(model_data, y_pos, args, minDomain, maxDomain):
    
    symmetric = False
    if (args.color == 'default'):
        cmap = mpl.cm.plasma
    if (args.color == 'red'):
        cmap = mpl.cm.YlOrRd
    mask = ((y_pos >= minDomain) & (y_pos <= maxDomain))

    if (mask.max()):    
        doPlot = True
        maxi = model_data[:, :, mask].max() * 1.01
        mini = model_data[:, :, mask].min() * 0.99
    
        if ((mini < 0.0) and (not args.alog)):
            symmetric = True
            cmap = mpl.cm.bwr
            maxi = abs(model_data[:, :, mask]).max() * 1.05
            mini = -maxi

    else:
        doPlot = False
        mini = -1.0
        maxi = 1.0

    if (args.mini < 1e31):
        mini = args.mini
    if (args.maxi > -1e31):
        maxi = args.maxi
        
    min_max_data = {'mini' : mini,
                    'maxi' : maxi,
                    'cmap' : cmap,
                    'symmetric' : symmetric,
                    'doPlot': doPlot,
                    'mask': mask}
    
    return min_max_data

# ----------------------------------------------------------------------------
# Get position at boundary with ghostcells
# ----------------------------------------------------------------------------

def get_boundaries(locations, nGhostCells):
    minLoc = (locations[nGhostCells-1] + locations[nGhostCells]) / 2.0
    maxLoc = (locations[-nGhostCells] + locations[-(nGhostCells+1)]) / 2.0
    return minLoc, maxLoc
    

# ----------------------------------------------------------------------------
# Plot value vs time
# ----------------------------------------------------------------------------

def plot_value_vs_time(times, values, var, loc, args, filename):

    fig = plt.figure(figsize=(10, 8.5))
    ax = fig.add_subplot(111)
    ax.plot(times, values)

    start = times[0].strftime('%b %d, %Y %H:%M')
    end = times[-1].strftime('%b %d, %Y %H:%M')
    ax.set_xlabel(start + ' to ' + end)
    if (args.mean):
        type = 'mean'
    else:
        type = 'integral'
    title = 'Global ' + type + ' of ' + var
    ax.set_ylabel(title)

    title = title + ' at {:.2f} km'.format(loc)
    ax.set_title(title)
        
    stime = times[0].strftime('%y%m%d')
    fileout = filename + '_' + stime + '.' + args.ext
    print(' ==> Writing file : ' + fileout)
    fig.savefig(fileout)

    return

# ----------------------------------------------------------------------------
# Write file with value vs time
# ----------------------------------------------------------------------------

def write_value_vs_time(times, values, var, loc, args, filename):

    if (args.mean):
        integral = 'global mean'
    else:
        integral = 'global integral'

    print(' ==> Writing file : ' + filename)

    fp = open(filename, "w")
    fp.write("\n")
    fp.write("#VAR\n")
    fp.write(var + "\n")

    fp.write("\n")
    fp.write("#INTEGRAL\n")
    fp.write(integral + "\n")

    fp.write("\n")
    fp.write("#ALTITUDE\n")
    ls = "%8.2f" % (loc)
    fp.write(ls + "\n")

    pwd = os.getcwd()
    fp.write("\n")
    fp.write("#DIRECTORY\n")
    fp.write(pwd + "\n")
    
    fp.write("\n")
    fp.write("#START\n")

    for i, t in enumerate(times):
        ts = t.strftime("%Y %m %d %H %M %S ")
        vs = "%e" % values[i]
        fp.write(ts + vs + "\n")

    fp.close()

    return


# ----------------------------------------------------------------------------
# Plot polar region
# ----------------------------------------------------------------------------

def plot_polar_region(fig, axplot, x_pos, y_pos, values,
                      utime, mask, whichPole,
                      miniPole, maxiPole, cmap, title, cbar_label):

    if (whichPole =='North'):
        ylabels = [r'80$^\circ$', r'70$^\circ$', r'60$^\circ$',
                   r'50$^\circ$']
        fac = 1.0
    else:
        ylabels = [r'-80$^\circ$', r'-70$^\circ$', r'-60$^\circ$',
                   r'-50$^\circ$']
        fac = -1.0
        
    # Find rotation for convertion to local time
    shift = time_conversion.calc_time_shift(utime)

    xlabels = []
    xlabelpos = []

    ylabelpos = [10.0, 20.0, 30.0, 40.0]
    xticks = np.arange(0, 2 * np.pi, np.pi / 2.0)
    yticks = np.arange(10, 50, 10)

    yp = 90.0 - fac * y_pos[mask]
    dy = (int(100.0*(yp[1]-yp[0]))/100.0)/2.0
    yp = np.append(yp - dy, yp[-1] + dy)
    xp = np.radians(x_pos + shift - 90.0)
    dx = (xp[1] - xp[0])/2
    xp = np.append(xp - dx, xp[-1] + dx)
    z = values[mask, :]
    axplot.grid(False)
    conn = axplot.pcolormesh(xp, yp, z,
                             shading = 'auto',
                             vmin = miniPole, vmax = maxiPole,
                             cmap = cmap)
    axplot.text(132.0/180. * np.pi, 72.5, title,
                verticalalignment = 'top',
                horizontalalignment = 'left',
                fontsize = 14)
    axplot.text(55.0/180.0 * np.pi, 46.0, cbar_label,
                verticalalignment = 'bottom',
                horizontalalignment = 'left',
                fontsize = 14)
    axplot.set_xticks(xlabelpos)
    axplot.set_xticklabels(xlabels)
    axplot.text(-np.pi/2, 46.0, '00 LT',
             verticalalignment='top', horizontalalignment='center')
    axplot.text(np.pi/2, 45.0, '12 LT',
             verticalalignment='bottom', horizontalalignment='center')
    axplot.text(-np.pi, 47.0, '18 LT',
             verticalalignment='center', horizontalalignment='center',
             rotation = 90)
    axplot.text(3*np.pi/4, 45.0, whichPole,
             verticalalignment='bottom',
             horizontalalignment='center', rotation = 45)
    axplot.set_yticks(ylabelpos)
    axplot.set_yticklabels(ylabels)
    axplot.grid(linestyle=':', color='black')
    axplot.set_xticks(xticks)
    axplot.set_yticks(yticks)
    axplot.set_ylim([0, 45])

    cbar = fig.colorbar(conn, ax = axplot, shrink=0.5, pad=0.005)

    return

# ----------------------------------------------------------------------------
# Get min index for time array
# ----------------------------------------------------------------------------

def get_min_time_index(times, utime):

    iTime = 0
    dtSave = np.abs((times[0] - utime).total_seconds())
    for i, t in enumerate(times):
        dt = np.abs((t - utime).total_seconds())
        if (dt < dtSave):
            iTime = i
            dtSave = dt
    return iTime
    
# ----------------------------------------------------------------------------
# Plot data locations on polar plots
# ----------------------------------------------------------------------------

def plot_data_locations_polar(axplot, DataLat, DataLon, DataTime,
                              utime, whichPole, data = [],
                              doPlotWind = False):

    # Find rotation for convertion to local time
    shift = time_conversion.calc_time_shift(utime)

    iTimeData = get_min_time_index(DataTime, utime)
    
    doPlotDataLoc = False
    if ((np.max(DataLat) > 45.0) and (whichPole == "North")):
        doPlotDataLoc = True
        r = 90.0 - DataLat
    if ((np.min(DataLat) < -45.0) and (whichPole == "Sorth")):
        doPlotDataLoc = True
        r = 90.0 + DataLat
    
    if (doPlotDataLoc):
        axplot.plot((DataLon + shift - 90.0) * np.pi/180.0,
                    r, color = 'yellow')

        if (doPlotWind):
            axplot.plot((DataLon[iTimeData] + shift - 90.0) * np.pi/180.0,
                        r[iTimeData], 'o', markersize = 4, color = 'yellow')
            timeM = data['merid_time']
            timeZ = data['zonal_time']
            zonal = data['zonal_wind']
            merid = data['merid_wind']
            iTz_ = get_min_time_index(timeZ, utime)
            iTm_ = get_min_time_index(timeM, utime)
            north = merid[iTm_]
            east = zonal[iTz_]
            x = (DataLon[iTimeData] + shift - 90.0) * np.pi/180.0
            y = r[iTimeData]
            fac = 1.0
            xwind = - fac * north * np.cos(x) - east * np.sin(x)
            ywind = - fac * north * np.sin(x) + east * np.cos(x)
            
            axplot.quiver(x, y, xwind, ywind, scale = 2500.0, color = 'grey')
        else:
            axplot.plot((DataLon[iTimeData] + shift - 90.0) * np.pi/180.0,
                        r[iTimeData], 'o', markersize = 10, color = 'yellow')
            
    return

# ----------------------------------------------------------------------------
# Plot winds in polar coordinates
# ----------------------------------------------------------------------------

def plot_winds_polar(axplot, x_pos, y_pos, xwinds, ywinds,
                     utime, mask, whichPole, nStep, scale):

    # Find rotation for convertion to local time
    shift = time_conversion.calc_time_shift(utime)

    if (whichPole =='North'):
        fac = 1.0
    else:
        fac = -1.0

    ewind = xwinds[mask, :]
    nwind = ywinds[mask, :]
    xp = np.array(np.radians(x_pos + shift - 90.0))
    yp = np.array(90.0 - fac * y_pos[mask])

    t, r = np.meshgrid(xp, yp)
    # north is towards pole in north and away in south:
    xwind = - fac * nwind * np.cos(t) - ewind * np.sin(t)
    ywind = - fac * nwind * np.sin(t) + ewind * np.cos(t)

    nx, ny = np.shape(xwind)
    for ix in range(0, nx, nStep):
        nStepL = int(nStep / np.sin(yp[ix]*np.pi/180.0))
        axplot.quiver(xp[::nStepL], yp[ix],
                      xwind[ix][::nStepL],
                      ywind[ix][::nStepL], scale = scale * 5.0)
    
    return

# ----------------------------------------------------------------------------
# 
# ----------------------------------------------------------------------------

def copy_dictionary(dict_in):

    dict_out = {}
    for key in dict_in.keys():
        dict_out[key] = dict_in[key]
    
    return dict_out


# ----------------------------------------------------------------------------
# 
# ----------------------------------------------------------------------------

def add_dir_to_files(header, dir):
    new_header = copy_dictionary(header)
    if (header['IsGitm']):
        ext = '.bin'
    else:
        ext = '.nc'
    DidWork = True
    for i, file in enumerate(new_header['filename']):
        m = re.match(r'(.*)(\d)'+ext, file)
        if m:
            file = glob(dir + '/' + m.group(1) + '?' + ext)
            if (len(file) > 0):
                new_header['filename'][i] = file[0]
            else:
                print('Can not file background file : ',file)
                DidWork = False
    new_header['DidWork'] = DidWork

    return new_header
            
# ----------------------------------------------------------------------------
# Define the main plotting routine

def plot_model_results():

    # Get the input arguments
    args = get_args()

    noPole = args.nopole
    if args.cut == 'lon':
        noPole = True
    if args.cut == 'lat':
        noPole = True

    header = get_file_info(args)
        
    if (args.list):
        list_vars(header)
        exit()

    if (args.var >= len(header["vars"])):
        raise ValueError("requested variable doesn't exist: {:d}>{:d}". \
                         format(args.var, len(header["vars"])))

    data = read_in_model_files(args, header)

    doSubtrackBack = False
    if (args.backdir != 'none'):
        header_back = add_dir_to_files(header, args.backdir)
        if (header_back['DidWork']):
            databack = read_in_model_files(args, header_back)
            doSubtrackBack = True
        else:
            print('Error! Need to Stop!!!')
            exit()
        all_winds_x = data['winds_x'] - databack['winds_x']
        all_winds_y = data['winds_y'] - databack['winds_y']
        if (args.percent):
            all_2dim_data = (data['slices'] - databack['slices']) / \
                databack['slices']
            all_int_data = (data['integrated_data'] - \
                            databack['integrated_data']) / \
                            databack['integrated_data']
        else:
            all_2dim_data = data['slices'] - databack['slices']
            all_int_data = data['integrated_data'] - \
                databack['integrated_data']
    else:
        all_winds_x = data['winds_x']
        all_winds_y = data['winds_y']
        all_2dim_data = data['slices']
        all_int_data = data['integrated_data']
        
    all_times = data['times']
    icut = data['icut']
    x_pos = data['x_pos']
    y_pos = data['y_pos']
    z_val = data['z_val']

    # Prepare the output filename

    if (args.tec):
        sVar = 'varTEC'
        varName = 'TEC'
    elif (args.on2):
        sVar = 'varON2'
        varName = 'ON2'
    else:
        sVar = 'var%02d' % args.var
        varName = header["vars"][args.var]
    if args.alog:
        sVar = sVar + 'L'
    
    filename = sVar + "_{:s}{:03d}".format(args.cut, icut)

    # ------------------------------------------------------------------
    # 1-d plot:
    if (args.timeplot):
        plot_value_vs_time(all_times, all_int_data,
                           varName, z_val, args,
                           args.timefile)
        if (args.timefile != "none"):
            write_value_vs_time(all_times, all_int_data,
                                varName, z_val, args,
                                args.timefile)
        exit()
    
    # ------------------------------------------------------------------
    # 2-d plots...
    
    # Does the user want to over plot locations?
    if (args.hiwind == 'none'):
        doPlotDataLoc = False
    else:
        doPlotDataLoc = True
        hiwind = read_hiwind_file(args.hiwind)
        DataLat = hiwind['lats']
        DataLon = (hiwind['lons'] + 360.0) % 360.0
        DataTime = hiwind['merid_time']
        DataTimeZ = hiwind['merid_time']
        DataZonal = hiwind['zonal_wind']
        DataMerid = hiwind['merid_wind']

    # Does the user want to add (a), (b), (c) labels?
    doLabel = False
    if (args.label):
        if (len(args.filelist) > 26):
            print("Can't label!  More than 26 files!!!")
        else:
            doLabel = True
            abc = "abcdefghijklmnopqrstuvwxyz"

    # If desired, take the log of the data
    if args.alog:
        all_2dim_data = np.log10(all_2dim_data)

    # scale data to eliminate exponentials, if needed    
    scaled_slice_data, factorString = rescale_data(all_2dim_data)
    all_2dim_data = scaled_slice_data

    if (args.tec):
        cbar_label = 'TEC'
    elif (args.on2):
        cbar_label = 'ON2'        
    else:
        cbar_label = header["vars"][args.var] + factorString
    
    # ----------------------------------------------------------
    # Define plotting limits    
    min_max_whole = get_min_max_data(all_2dim_data,
                                     y_pos, args, -1e32, 1e32)
    symmetric = min_max_whole['symmetric']
    cmap = min_max_whole['cmap']
    mini = min_max_whole['mini']
    maxi = min_max_whole['maxi']

    if args.winds:
        windMag = np.sqrt(all_winds_x**2 + all_winds_y**2)
        windScale = np.round(np.max(windMag))
        
    doPlotGlobal = True
    if args.cut == 'alt':

        min_max_north = get_min_max_data(all_2dim_data,
                                         y_pos, args, 40.0, 90.0)
        mask_north = min_max_north['mask']
        plot_north = min_max_north['doPlot']
        maxi_north = min_max_north['maxi']
        mini_north = min_max_north['mini']
        
        min_max_south = get_min_max_data(all_2dim_data,
                                         y_pos, args, -90.0, -40.0)
        mask_south = min_max_south['mask']
        plot_south = min_max_south['doPlot']
        maxi_south = min_max_south['maxi']
        mini_south = min_max_south['mini']

        if (np.min(y_pos) > -40.0):
            plot_south = False
        if (np.max(y_pos) < 40):
            plot_north = False

        if (args.north or args.south):
            doPlotGlobal = False
            if (args.north):
                plot_south = False
            else:
                plot_north = False
        if (noPole):
            plot_north = False
            plot_south = False
        

    # Define plot range
    minx, maxx = get_boundaries(x_pos, 2)
    miny, maxy = get_boundaries(y_pos, 2)

    if args.movie > 0:
        img_file_fmt = movie_routines.setup_movie_dir(filename)
    else:
        img_file_fmt = filename+'_{:}.'+args.ext

    dpi = 120
    # Create a plot for each time
    for itime, utime in enumerate(all_times):
        # Initialize the figure

        data2d = all_2dim_data[itime].transpose()
        if args.winds:
            xwind = np.array(all_winds_x[itime].transpose())
            ywind = np.array(all_winds_y[itime].transpose())

        outfile_add = ''
        if (doPlotGlobal):
            if (noPole):
                fig = plt.figure(constrained_layout=False, figsize=(10, 5.5),
                                 dpi = dpi)
            else:
                fig = plt.figure(constrained_layout=False, figsize=(10, 8.5),
                                 dpi = dpi)
            xSize = 10.0 * dpi
            ySize = 8.5 * dpi
            if (noPole):
                ax = fig.add_axes([0.07, 0.09, 0.97, 0.83])
            else:
                ax = fig.add_axes([0.07, 0.06, 0.97, 0.48])
                # Top Left Graph Northern Hemisphere
                if (plot_north):
                    ax2 = fig.add_axes([0.06, 0.55, 0.425, 0.43],
                                        projection='polar')
                # Top Right Graph Southern Hemisphere
                if (plot_south):
                    ax3 = fig.add_axes([0.535, 0.55, 0.425, 0.43],
                                       projection='polar')
        else:
            fig = plt.figure(constrained_layout=False, figsize=(5.4, 5))
            xSize = 5.4 * dpi
            ySize = 5.0 * dpi
            pos = [0.04, 0.04, 0.95, 0.85]
            if (plot_north):
                ax2 = fig.add_axes(pos, projection='polar')
                outfile_add = 'n'
            else:
                ax3 = fig.add_axes(pos, projection='polar')
                outfile_add = 's'

        if ((not args.tec) and (not args.on2)):
            title = "{:s}\n{:s}: {:.2f} {:s}".format(
                utime.strftime("%d %b %Y %H:%M:%S UT"), args.cut, z_val,
                'km' if args.cut == 'alt' else r'$^\circ$')
        else:
            title = utime.strftime("%d %b %Y %H:%M:%S UT")
        if args.winds:
            title = title + '; max wind: %.0f m/s' % windScale
        if (doLabel):
            title = '('+abc[itime]+') ' + title

        pwd = os.getcwd()
        filename = pwd + '/' + header['filename'][itime]
        if (doSubtrackBack):
            filename = filename + ' -- ' + args.backdir
            outfile_add = outfile_add + 'd'
        plt.text(xSize * 0.995, 10.0, filename, transform = None,
                 fontsize = 6,
                 verticalalignment='bottom',
                 horizontalalignment='center', rotation = 90)
                 

        if (doPlotGlobal):
            
            # Plot the global data set (square plot at bottom if three plots):

            dx = (x_pos[1] - x_pos[0])/2.0
            xp = np.append(x_pos - dx, x_pos[-1:]+dx)
            dy = (y_pos[1] - y_pos[0])/2.0
            yp = np.append(y_pos - dy, y_pos[-1]+dy)
            con = ax.pcolormesh(xp, yp, data2d,
                                vmin = mini, vmax = maxi, cmap = cmap,
                                shading = 'auto')
            ax.set_ylim([miny, maxy])
            ax.set_xlim([minx, maxx])

            # Set the labels and aspect ratio
            ax.set_title(title)
            ax.set_xlabel(r'Latitude ($^\circ$)' if args.cut == 'lon'
                          else r'Longitude ($^\circ$)')
            ax.set_ylabel(r'Latitude ($^\circ$)' if args.cut == 'alt'
                          else r'Altitude (km)')
            if args.cut == 'alt':
                ax.set_aspect(1.0)

            # Add the winds, if desired
            if args.winds:
                xwinds = xwind[::args.nstep, ::args.nstep]
                ywinds = ywind[::args.nstep, ::args.nstep]
                xp = np.array(x_pos)[::args.nstep]
                yp = np.array(y_pos)[::args.nstep]
                ax.quiver(xp, yp, xwinds, ywinds, scale = windScale * 10.0)

            # plot data locations:
            if (doPlotDataLoc):
                ax.plot(DataLon, DataLat, 'o', color = 'yellow', markersize=1)
                iTimeData = get_min_time_index(DataTime, utime)
                xp = [DataLon[iTimeData]]
                yp = [DataLat[iTimeData]]
                ax.plot(DataLon[iTimeData], DataLat[iTimeData],
                        'o', markersize = 10, color = 'yellow')
                if (args.winds):
                    yDataWind = [DataMerid[iTimeData]]
                    iTimeData = get_min_time_index(DataTimeZ, utime)
                    xDataWind = [DataZonal[iTimeData]]
                    ax.quiver(xp, yp, xDataWind, yDataWind,
                              scale = windScale * 10.0, color = 'yellow')
                            
            # Set the colorbar
            cbar = fig.colorbar(con, ax = ax, shrink = 0.75, pad = 0.02)
            cbar.set_label(cbar_label, rotation=90)

        title_temp = ''
        if (not doPlotGlobal):
            title_temp = title
            
        # If this is an altitude slice, add polar dials
        if (args.cut == 'alt') and \
           (plot_north or plot_south) and \
           (not noPole):

            if plot_north:
                plot_polar_region(fig, ax2, x_pos, y_pos, data2d,
                                  utime, mask_north, 'North',
                                  mini_north, maxi_north, cmap,
                                  title_temp, cbar_label)

                # overplot data locations:
                if (doPlotDataLoc):
                    plot_data_locations_polar(ax2, \
                                              DataLat, DataLon, DataTime, \
                                              utime, 'North', \
                                              data = hiwind, \
                                              doPlotWind = args.winds)
                
                # Add the winds, if desired
                if args.winds:
                    plot_winds_polar(ax2, x_pos, y_pos,
                                     all_winds_x[itime].transpose(),
                                     all_winds_y[itime].transpose(),
                                     utime, mask_north, 'North', args.nstep, windScale)
                
            if plot_south:
                plot_polar_region(fig, ax3, x_pos, y_pos, data2d,
                                  utime, mask_south, 'South',
                                  mini_south, maxi_south, cmap,
                                  title_temp, cbar_label)
                    
                # overplot data locations:
                if (doPlotDataLoc):
                    plot_data_locations_polar(ax3, DataLat, DataLon, DataTime,
                                              utime, 'South', \
                                              data = hiwind, \
                                              doPlotWind = args.winds)
                
                # Add the winds, if desired
                if args.winds:
                    plot_winds_polar(ax3, x_pos, y_pos,
                                     all_winds_x[itime].transpose(),
                                     all_winds_y[itime].transpose(),
                                     utime, mask_south, 'South', args.nstep, windScale)
                
        # Format the output filename
        if args.movie > 0:
            fmt_input = itime
        else:
            fmt_input = utime.strftime('%y%m%d_%H%M%S') + outfile_add
        outfile = img_file_fmt.format(fmt_input)

        # Save the output file
        print(" ==> Writing file : ", outfile)
        fig.savefig(outfile, dpi = dpi)
        plt.close(fig)

    # Create a movie, if desired
    if args.movie > 0:
        movie_routines.save_movie(filename, ext=args.ext,
                                  rate=args.rate)
        
    return


# Needed to run main script as the default executable from the command line
if __name__ == '__main__':
    plot_model_results()
