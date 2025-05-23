#!/usr/bin/env python3

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import argparse
from datetime import datetime
from datetime import timedelta
import os
from glob import glob
import numpy as np
import re

# running this without any parameters produces an f107_downloaded.txt file.
# this can be copied into ../srcData/f107.txt and committed.

# ----------------------------------------------------------------------
# There are two sites that seem to have this data:
#
# https://chain-new.chain-project.net/echaim_downloads/apf107.dat
#  - last three "columns" are F107 related.
#  - first of three is F107 adjusted
#  - second seems to be measured / observed
#  - third seems to be 81-day average
#
# ftp://ftpsedr.cls.fr/pub/previsol/solarflux/observation/radio_flux_adjusted_observation.txt
#  - columns are well labeled
#  - this code uses f10.7_c (observed values with interpolated data
#                            gaps and flare correction)
#  - above this, it says: "Daily averages adjusted to 1 AU"
# ----------------------------------------------------------------------

IsVerbose = True

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args_f107():

    parser = argparse.ArgumentParser(description = 'Download F107 data and process')

    parser.add_argument('-source', metavar = 'source',  default ='ftp', \
                        choices = ['ftp', 'https',
                                   'radioflux_temp.txt', 'apf107.dat'], 
                        help = 'ftp,https : which source would you like (including a file!)')
    parser.add_argument('-outfile', metavar = 'outfile', \
                        help = 'output file name',
                        default = 'f107_downloaded.txt')

    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------
#
# ----------------------------------------------------------------------

def write_f107_file(data, message = "none", fileout = "none"):

    if (fileout == "none"):
        ymd = data['times'][0].strftime('%Y%m%d')
        fileout = 'f107_' + ymd + '.txt'
    print(' -> Writing file ' + fileout)
        
    fp = open(fileout, 'wb')
    fp.write("#--------------------------------------------------\n".encode())
    fp.write("#>\n".encode())
    fp.write("#Element: adjusted\n".encode())
    fp.write("#Description: Adjusted daily solar radio flux\n".encode())
    fp.write("#Measure units: W/m^2/Hz\n".encode())
    fp.write("#Origin:  f107_download.py code\n".encode())
    if (message != "none"):
        m = "# " + message + "\n"
        fp.write(m.encode())
    fp.write("#\n".encode())
    fp.write("#Sampling: 1 day\n".encode())
    fp.write("#Missing value: 1.0E33\n".encode())
    fp.write("#>\n".encode())
    fp.write("#yyyy-MM-dd HH:mm value qualifier description\n".encode())

    for i, t in enumerate(data['times']):
        f107 = data['f107'][i]
        out = " %6.1f    ''   ''" % (f107)
        ymdhm = t.strftime('%Y-%m-%d %H:%M')
        line = ymdhm + out + "\n"
        fp.write(line.encode())

    fp.close()

    return fileout

# ----------------------------------------------------------------------
# do system command
# ----------------------------------------------------------------------

def run_command(command):
    if (IsVerbose):
        print("   -> Running Command : ")
        print("      ", command)
    os.system(command)
    return True

# ----------------------------------------------------------------------
# 
# ----------------------------------------------------------------------

def download_ap107_file():

    local_file = "apf107_temp.txt"
    command = "wget"
    command = command + " -O " + local_file 
    command = command + " https://chain-new.chain-project.net/echaim_downloads/apf107.dat"
    didWork = run_command(command)
    if (not didWork):
        local_file = "none"
    return local_file
    
# ----------------------------------------------------------------------
# 
# ----------------------------------------------------------------------

def read_ap107_file(file):

    if (file == "none"):
        local_file = download_ap107_file()
    else:
        local_file = file

    if (os.path.exists(local_file)):
        times = []
        f107 = []
        print(' -> Reading file : ', local_file)
        fpin = open(local_file, 'r')
        for line in fpin:
            iYear = int(line[1:3])
            if (iYear > 50):
                iYear += 1900
            else:
                iYear += 2000
            iMonth = int(line[4:6])
            iDay = int(line[7:9])
            iF107 = float(line[39:44])
            times.append(datetime(iYear, iMonth, iDay))
            f107.append(iF107)
        fpin.close()
    else:
        print("Can't seem to find file : " + local_file)
        print("download probably didn't work!")
        exit()

    data = {'times' : np.array(times),
            'f107' : np.array(f107)}
        
    return data
        
# ----------------------------------------------------------------------
# 
# ----------------------------------------------------------------------

def download_https_file():

    local_file = "radioflux_temp.txt"
    command = "wget"
    command = command + " -O " + local_file 
    command = command + " ftp://ftpsedr.cls.fr/pub/previsol/solarflux/observation/radio_flux_adjusted_observation.txt"
    didWork = run_command(command)
    if (not didWork):
        local_file = "none"
    return local_file
    
# ----------------------------------------------------------------------
# 
# ----------------------------------------------------------------------

def read_https_file(file):

    if (file == "none"):
        local_file = download_https_file()
    else:
        local_file = file

    if (os.path.exists(local_file)):
        times = []
        f107 = []
        print(' -> Reading file : ', local_file)
        fpin = open(local_file, 'r')

        # Read header:
        for line in fpin:
            m = re.match(r'.*year month.*',line)
            if m:
                break
            
        for line in fpin:
            cols = line.split()
            times.append(datetime(int(cols[0]),
                                  int(cols[1]),
                                  int(cols[2])))
            # This is the f10.7_c column
            f107.append(float(cols[13]))
        fpin.close()
    else:
        print("Can't seem to find file : " + local_file)
        print("download probably didn't work!")
        exit()

    data = {'times' : np.array(times),
            'f107' : np.array(f107)}
        
    return data
        
# ----------------------------------------------------------------------
# 
# ----------------------------------------------------------------------

def plot_f107(f107, file = "f107.png"):

    fig = plt.figure(figsize=(10, 8.5))
    ax = fig.add_subplot(111)
    ax.plot(f107['times'], f107['f107'])

    start = f107["times"][0].strftime('%b %d, %Y %H:%M')
    end = f107["times"][-1].strftime('%b %d, %Y %H:%M')
    ax.set_xlabel(start + ' to ' + end)
    title = 'F10.7 (Adjusted to 1 AU)'
    ax.set_ylabel(title)

    print('Writing file : ' + file)
    fig.savefig(file)

    return



#------------------------------------------------------------------------------
# main code:
#------------------------------------------------------------------------------

if __name__ == '__main__':
    
    args = parse_args_f107()

    if (args.source == "https"):
        data = read_https_file("none")
    if (args.source == "radioflux_temp.txt"):
        data = read_https_file(args.source)

    if (args.source == "ftp"):
        data = read_ap107_file("none")
    if (args.source == "apf107.dat"):
        data = read_ap107_file(args.source)
    write_f107_file(data, fileout = args.outfile)
    plot_f107(data)
