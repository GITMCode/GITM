#!/usr/bin/env python

import sys
import datetime as dt
import numpy as np
from omniweb import *
from swmf_imf import *
from create_fake_noaa_hpi_input import *
from supermag_download_ae import write_sme_file
import argparse
import matplotlib.pyplot as plt
import matplotlib.dates as dates

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args_omni():

    parser = argparse.ArgumentParser(description = 'Plot AMIE files')
    parser.add_argument('start', metavar = 'start', nargs = 1, \
                        help = 'start date as YYYYMMDD')
    parser.add_argument('end', metavar = 'end', nargs = 1, \
                        help = 'end date as YYYYMMDD')
    parser.add_argument('-swmf', \
                        help='output swmf style file (imfYYMMDD.dat)', \
                        action="store_true")
    parser.add_argument('-hp', \
                        help='output fake hemispheric power using AE', \
                        action="store_true")
    parser.add_argument('-ae', \
                        help='output AU, AL, AE in SME format', \
                        action="store_true")

    parser.add_argument('-hemi', \
                        help='Adjust hp by season? (must be run with -hp)', \
                        action="store_true")


    args = parser.parse_args()

    return args

def write_omni_file(lines, fileout):
    print("--> Writing file : ", fileout)
    with open(fileout, "w") as file:
        file.writelines("%s" % l for l in lines)

    
#------------------------------------------------------------------------------
#SCRIPT USE
#example command line input: python omniweb_read.py 20110620 20110623 -all

if __name__ == '__main__':

    args = parse_args_omni()

    #assuming first two args are the start/end dates, then info desired

    start = args.start
    if (not np.isscalar(start)):
        start = start[0]
    end = args.end
    if (not np.isscalar(end)):
        end = end[0]


    if args.hemi and not args.hp:
        raise ValueError(
            'To split HP by hemisphere, option to write Hemispheric '
            'Power must be enabled.')

    print("-> Downloading OMNI data using ", start, " -> ", end)
    results = download_omni_data(start, end, "-all")
    print("-> Parsing data")
    omniDirty = parse_omni_data(results)
    print("-> Cleaning data")
    data = clean_omni(omniDirty)

    if (args.swmf):
        fileout = data["times"][0].strftime('imf%Y%m%d.dat')
        message = "Data downloaded from OMNIWeb and processed by omniweb_read.py\n"
        print("-> Writing SWMF style output: ", fileout)
        write_swmf_imf_file(data, fileout, message)
    else:
        fileout = data["times"][0].strftime('omni_%Y%m%d.txt')
        write_omni_file(results, fileout)

    if (args.hp):
        fileout = data["times"][0].strftime('ae_power_%Y%m%d.dat')
        message = "Data downloaded from OMNIWeb and processed by create_fake_noaa_hpi_input\n"

        hp = calculate_hp_from_ae(np.array(data["ae"]))
        print("-> Writing hemispheric power type of file: ", fileout)
        write_derived_hp(data["times"], hp, output_filename = fileout,
                        add_seasonal_dependence=args.hemi)

    if (args.ae):
        message = "Data downloaded from OMNIWeb and written in SME format\n"
        print("-> Writing SME type of file")
        write_sme_file(data, message = message)

    print("-> Plotting data")
    plot_imf(data)
