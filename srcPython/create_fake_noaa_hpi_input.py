#!/usr/bin/env python

# Author: Dogacan S. Ozturk

import os
import argparse, sys
import datetime as dt
import numpy as np

try:
    from spacepy.pybats import kyoto
except ModuleNotFoundError:
    print('Spacepy not found. If using SME data, please pass it as a numpy array or a list to calculate_hp_from_ae.')

def dtime_to_doy(d_time):
    """
    Convert a pythonic datetime to day-of-year.

    Inputs
    -----
    d_time (datetime.datetime): Pythonic date_time.

    Outputs:
    -------
    int : Days of year.

    """

    return (d_time - dt.datetime(d_time.year, 1, 1)).days + 1

def calculate_hp_from_ae(ae):
    '''
    HP is in GW.
    AE is in nT.
    Formula taken from Wu et al, 2021. https://doi.org/10.1029/2020SW002629
    '''
    hp = 0.102*ae + 8.953
    return hp

def write_derived_hp(time_array, hp, 
                     output_filename = "empty", 
                     add_seasonal_dependence=False):

    if (output_filename == "empty"):
    
        savedir = './hp_from_ae'
        if not os.path.exists(savedir):
            os.mkdir(savedir)
        output_filename = os.path.join(savedir,
            'power_from_ae_{0:%Y%m%d}'.format(time_array[0]) + \
            '_to_{0:%Y%m%d}.txt'.format(time_array[-1]))

    output_file = open(output_filename, 'w')

    fmt_line = '{0:%Y-%m-%d} {0:%H:%M:%S} NOAA-17 (N)  6{1:7.2f}   0.75\n'
    ntimes = len(time_array)
    
    output_file.write(':Data_list: '+output_filename+'\n'.format(time_array[0]))
    output_file.write(':Created: {0:%a %b %d %H:%M:%S UTC %Y\n}'.format(dt.datetime.now()))
    output_file.write('# This file is created to replicate NOAA HPI files created before 2013.\n')
    output_file.write('# Please use with caution as Sat number, hemisphere, activity level, and normalizing factors are placeholders.\n')
    output_file.write('#\n')
    output_file.write('# Source: AE or SME index.\n')
    output_file.write('# Units: gigawatts\n\n')
    output_file.write('# Format:\n\n')
    output_file.write('# Each line is formatted as in this example:\n\n')
    output_file.write('# 2006-09-05 00:54:25 NOAA-16 (S)  7  29.67   0.82\n\n')
    output_file.write('# A19   Date and UT at the center of the polar pass as YYYY-MM-DD hh:mm:ss\n')
    output_file.write('# 1X    (Space)\n')
    output_file.write('# A7    NOAA POES Satellite number\n')
    output_file.write('# 1X    (Space)\n')
    output_file.write('# A3    (S) or (N) - hemisphere\n')
    output_file.write('# I3    Hemispheric Power Index (activity level)\n')
    output_file.write('# F7.2  Estimated Hemispheric Power in gigawatts\n')
    output_file.write('# F7.2  Normalizing factor\n\n')
    
    for i in range(ntimes):
        if add_seasonal_dependence:
            day_of_year = dtime_to_doy(time_array[i])

            # Ten percent correction factor (by *roughly* day of year)
            # Leap years would be off by 1 day, but that is OK.
            # From (https://doi.org/10.1029/2006GL028444)
            ten_percent_corr = .1 * np.cos(2 * np.pi * day_of_year/365)

            # Northern hemisphere
            line_north = fmt_line.format(time_array[i], 
                                        (1 + ten_percent_corr) * hp[i])
            output_file.write(line_north)

            # Southern hemisphere
            line_south = fmt_line.format(time_array[i], 
                                        (1 - ten_percent_corr) * hp[i])
            line_south = line_south.replace('(N)', '(S)') # tag w. sh
            output_file.write(line_south)

        else:
            output_file.write(fmt_line.format(time_array[i], hp[i]))
        
    output_file.close()

#__________________________________________________________________________#
#    To use the code and generate a fake HPI file modify the dates below.  #
#__________________________________________________________________________#
if __name__ == '__main__':


    # Set parameters here, or through command line.
    # Default is to use command line arguments, but if nothing is provided,
    #   these values will be used:

    t_start = dt.datetime(2018,8,1)
    t_end = dt.datetime(2018,9,1)

    add_seasonal_dependence = False
    # This will add a 10% seasonal correction factor to HP
    # Derived from Ridley, A. J. (2007), doi:10.1029/2006GL028444.

    out_filename = 'empty' # Defaults to hp_from_aepower_from_ae_[start]_to_[end].txt


    ########     ~~~~~~~~~~~~     ########
    # Parse command line arguments, if this script is called with any.
    # If there are none, the values set above are used to ensure compatibility with existing workflows.
    # There is no difference to functionality whether command line args are used or not.
    # They are all optional and will override the above values when set...

    parser = argparse.ArgumentParser()

    parser.add_argument('-start', type=str,
                        help='Start date, formatted as YYYYMMDD.'
                        ' Default is 20180801')
    parser.add_argument('-end', type=str, 
                        help='End date, formatted as YYYYMMDD.'
                        ' Default is 20180901')
    parser.add_argument('-out', type=str,
                        help='Name of output file. Defaults to '
                              'hp_from_ae/power_from_ae_[start]_to_[end].txt')
    parser.add_argument('--add_season', action='store_true', default=None,
                        help='Add a 10%% seasonal factor to hp? Default is False.'
                        ' Include this flag to enable.')

    args = parser.parse_args()
        
    if len(sys.argv) != 1: # Check if args provided, overwrites above if there are.
        # If an argument is not provided, use default value set above.
        t_start = dt.datetime.strptime(args.start, '%Y%m%d') if args.start else t_start
        t_end = dt.datetime.strptime(args.end, '%Y%m%d') if args.end else t_end
        out_filename = args.out if args.out else out_filename
        add_seasonal_dependence = args.add_season if args.add_season else add_seasonal_dependence

    ########    ~~~~~~~~~~~~    ########

    try:
        ae_data = kyoto.aefetch(t_start, t_end)
        ae = kyoto.KyotoAe(ae_data)
        hp = calculate_hp_from_ae(ae['ae'])
        write_derived_hp(ae['time'], hp, 
                         output_filename=out_filename,
                         add_seasonal_dependence=add_seasonal_dependence)

        if add_seasonal_dependence:
            print('\t-------------------------\n'
                  'Added a small seasonal correction factor to HP.\n'
                  'Make sure to enable DoSeparateHPI in GITM!\n'
                  '\t-------------------------\nDone.')
        
    except AttributeError:
        print('Spacepy Kyoto library not found. For SME derived fake HPI, use standalone functions in this file.')



