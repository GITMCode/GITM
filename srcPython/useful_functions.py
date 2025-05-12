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
from aetherpy.plot import data_prep, movie_routines

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

def fix_vars(vars):
    newvars = []
    for v in vars:
        nv = re.sub('!U', '', v)
        nv = re.sub('!N', '', nv)
        nv = re.sub('!D', '', nv)
        newvars.append(nv)

    return newvars

# ----------------------------------------------------------------------------
# Read in all of the model files:
# ----------------------------------------------------------------------------

def read_in_model_files(args, header):

    # Define the plotting inputs
    plot_vars = [0, 1, 2, args.var]

    # Update plotting variables to include the wind, if desired
    if args.winds:
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
            print('=> Reading file : ', filename)
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
            lons = np.degrees(data[0][:, 0, 0])  # Convert from rad to deg
            lats = np.degrees(data[1][0, :, 0])  # Convert from rad to deg
            # Find the desired index to cut along to get a 2D slice
            isgrid = False

            if (args.cut == 'lon'):
                pos = args.lon
            if (args.cut == 'lat'):
                pos = args.lat
            
            if (args.cut == 'alt'):
                pos = args.alt
                if (len(alts) == 1):
                    print("Only one alt found, setting alt pos = 0");
                    pos = 0
                    isgrid = True
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

