#!/usr/bin/env python3

import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as dates
import argparse
import os
import re
import juliandate as jd

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args_fism():

    parser = argparse.ArgumentParser(description = 'Download and process FISM data')
    parser.add_argument('-start', \
                        help = 'start date as YYYYMMDD', default = '0')
    parser.add_argument('-end', \
                        help = 'end date as YYYYMMDD', default = '0')
    parser.add_argument('-euvfile', \
                        help='EUV file that provides wavelength bins',
                        default = 'euv_59.csv')
    parser.add_argument('-fismfile', \
                        help='FISM file to convert',
                        default = 'none')
    parser.add_argument('-gitm', \
                        help='Write out GITM-style file',
                        action='store_true')
    parser.add_argument('-flare', \
                        help='Download and process FISM flare data',
                        action='store_true')

    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------
# Take string time YYYYMMDD.HHMM and convert to datetime
# ----------------------------------------------------------------------

def convert_ymdhm_to_dt(ymdhm):

    year = ymdhm[0:4]
    mo = ymdhm[4:6]
    da = ymdhm[6:8]

    if (len(ymdhm) >= 11):
        hr = ymdhm[9:11]
        if (len(ymdhm) >= 13):
            mi = ymdhm[11:13]
        else:
            mi = '00'
    else:
        hr = '00'
        mi = '00'
        
    time = dt.datetime(int(year), int(mo), int(da), int(hr), int(mi), 0)
    return time

# ----------------------------------------------------------------------
# Function to download FISM2 data.
# Need call to be in form:
# https://lasp.colorado.edu/lisird/latis/dap/fism_daily_hr.csv?&time>=2020-01-01T00:00:00.000Z&time<=2020-01-31T00:00:00.000Z
# ----------------------------------------------------------------------

def download_fism2(start, end, isFlare):

    sStart = start.isoformat()+'.000Z'
    sEnd = end.isoformat()+'.000Z'

    if (isFlare):
        site = 'https://lasp.colorado.edu/lisird/latis/dap/fism_flare_hr.csv'
    else:
        site = 'https://lasp.colorado.edu/lisird/latis/dap/fism_daily_hr.csv'

    url = site + '?&time>=' + sStart + '?&time<=' + sEnd
    
    ymdS = start.strftime('%Y%m%d')
    ymdE = end.strftime('%Y%m%d')
    filename = '.fism_raw_' + ymdS + '_to_' + ymdE + '.txt'

    command = 'curl "' + url + '" > ' + filename
    print("Running Command : ", command)
    os.system(command)
    
    return filename
    

#------------------------------------------------------------------------------
# Convert yearday (YYYYDDD) to date time
#   - center the day at 12 UT.
#------------------------------------------------------------------------------

def convert_time(yearday):
    year = int(yearday/1000.0)
    base = dt.datetime(year,1,1,12)
    doy = yearday - year * 1000.0
    time = base + dt.timedelta(days = doy-1)
    return time

#------------------------------------------------------------------------------
# Convert seconds since 1970 to date time
#   - center the day at 12 UT.
#------------------------------------------------------------------------------

def convert_time_seconds(seconds):
    base = dt.datetime(1970,1,1)
    time = base + dt.timedelta(seconds = seconds)
    return time

#------------------------------------------------------------------------------
# read euv file - this determines the wavelength bins
#------------------------------------------------------------------------------

def read_euv_csv_file(file):

    fpin = open(file, 'r')

    iFound = 0
    afac = []
    f74113 = []
    for line in fpin:
        aline = line.split(',')
        s = aline[-1].strip().split('.')[0]
        if (aline[0].strip() == "Short"):
            if (s.isnumeric()):
                short = np.asarray(aline[5:], dtype = float)
            else:
                short = np.asarray(aline[5:-1], dtype = float)
            iFound += 1
        if (aline[0].strip() == "Long"):
            if (s.isnumeric()):
                long = np.asarray(aline[5:], dtype = float)
            else:
                long = np.asarray(aline[5:-1], dtype = float)
        if (aline[0].strip() == "F74113"):
            if (s.isnumeric()):
                f74113 = np.asarray(aline[5:], dtype = float)
            else:
                f74113 = np.asarray(aline[5:-1], dtype = float)
            iFound += 1
        if (aline[0].strip() == "AFAC"):
            if (s.isnumeric()):
                afac = np.asarray(aline[5:], dtype = float)
            else:
                afac = np.asarray(aline[5:-1], dtype = float)
            iFound += 1
    # Save and convert from Angstroms to nm (FISM is in nm)
    wavelengths = {'short': short/10.0,
                   'long': long/10.0,
                   'afac': afac,
                   'f74113': f74113}
    return wavelengths


#------------------------------------------------------------------------------
# calculate the energy fluxes for euvac
#------------------------------------------------------------------------------

def calc_euvac(wavelengths, f107, f107a):

    euvac = []
    if (len(wavelengths['afac']) > 0):
        afac = wavelengths['afac']
        f74113 = wavelengths['f74113']
        f107bar = (f107 + f107a)/2
        slope = 1.0 + afac * (f107bar - 80.0)
        slope[slope < 0.8] = 0.8
        # in photons/cm2/s
        intensity = f74113 * slope * 1e9
        # convert to /m2
        intensity = intensity * 1e4

        # assume nm, convert to m:
        ave_wave = (wavelengths['long'] + wavelengths['short'])/2.0 * 1e-9
        # wavelength to energy:
        energy = 6.62607015e-34 * 2.99792e8 / ave_wave
        # convert from photons/s/m2 to W/m2:
        euvac = intensity * energy
        
    return euvac

#------------------------------------------------------------------------------
# read fism2 csv file
#------------------------------------------------------------------------------

def read_fism_csv_file(file):

    fpin = open(file, 'r')

    header = fpin.readline()
    vars = header.split(',')

    isSeconds = False

    m = re.match(r'.*seconds.*',vars[0])
    if m:
        isSeconds = True

    # Read in a 1d time array and 1d wavelength array
    # read in all of the irradiance and uncertainty data as
    # a large 1d array, then reshape it into a time vs wavelength 2d array
    # structure of the file is all wavelengths for one time, then
    # a new time and all wavelengths for that time.
    
    iRow = 0
    oldtime = 0.0
    nWaves = 0
    nTimes = 0
    times = []
    wavelengths = []
    irradiance = []
    uncertainty = []
    for line in fpin:
        aline = line.split(',')
        yearday = float(aline[0])
        if (yearday/1000 > 2100):
            isJulian = True
        else:
            isJulian = False
        irradiance.append(float(aline[2]))
        uncertainty.append(float(aline[3]))
        if (iRow == 0):
            oldtime = yearday

        if (yearday > oldtime):
            nWaves = 1
            if (isJulian):
                iTime = jd.to_gregorian(yearday)
                current = dt.datetime(iTime[0], iTime[1], iTime[2],
                                      iTime[3], iTime[4], iTime[5], iTime[6])
                times.append(current)
            else:
                if isSeconds:
                    times.append(convert_time_seconds(oldtime))
                else:
                    times.append(convert_time(oldtime))
            oldtime = yearday
            nTimes = nTimes + 1
        else:
            nWaves = nWaves + 1

        if (nTimes == 0):
            wavelengths.append(float(aline[1]))
            
        iRow = iRow+1

    if isSeconds:
        times.append(convert_time_seconds(oldtime))
    else:
        times.append(convert_time(oldtime))
    nTimes = nTimes + 1
    
    irradiance = np.array(irradiance).reshape((nTimes, nWaves))
    uncertainty = np.array(uncertainty).reshape((nTimes, nWaves))
    
    fism_data = {'vars': vars,
                 'time': times,
                 'nTimes': nTimes,
                 'wave': np.array(wavelengths),
                 'irr': irradiance,
                 'unc': uncertainty}
            
    return fism_data

#------------------------------------------------------------------------------
# 
#------------------------------------------------------------------------------

def find_index_from_time(data, timeToFind):
    tInDays = (timeToFind - data['times'][0]).total_seconds()/86400.0
    dist = np.abs(data['timesInDays'] - tInDays)
    iMin = np.argmin(dist)
    return iMin

#------------------------------------------------------------------------------
# 
#------------------------------------------------------------------------------

def read_ap107_file(local_file = "apf107_temp.txt"):

    times = []
    timesInDays = []
    f107 = []
    if (os.path.exists(local_file)):
        print(' -> Reading file : ', local_file)
        fpin = open(local_file, 'r')
        for iPt, line in enumerate(fpin):
            iYear = int(line[1:3])
            if (iYear > 50):
                iYear += 1900
            else:
                iYear += 2000
            iMonth = int(line[4:6])
            iDay = int(line[7:9])
            iF107 = float(line[39:44])
            times.append(dt.datetime(iYear, iMonth, iDay))
            if (iPt > 0):
                timesInDays.append((times[-1] - times[0]).total_seconds()/86400.0)
            else:
                timesInDays.append(0.0)
            f107.append(iF107)
        fpin.close()
    else:
        print('f107 read failed...')
        print("Can't seem to find file : " + local_file)
        print("download probably didn't work!")

    data = {'times' : np.array(times),
            'timesInDays' : np.array(timesInDays),
            'f107' : np.array(f107)}

    return data

#------------------------------------------------------------------------------
# rebin FISM data into new wavelength bins
#  some bins are single wavelength, and some span lots of wavelenghts
#------------------------------------------------------------------------------

def rebin_fism(fism_waves, fism_vals, wavelengths):

    shorts = wavelengths['short']
    longs = wavelengths['long']
    nWaves = len(shorts)
    new_irr = np.zeros(nWaves)
    ave_wav = np.zeros(nWaves)

    # first go through all of the wavelengths that are singular
    for iWave, short in enumerate(shorts):
        long = longs[iWave]
        if (long == short):
            d = np.abs(fism_waves - short)
            i = np.argmin(d)
            new_irr[iWave] = fism_vals[i] * \
                (fism_waves[i+1] - fism_waves[i])
            # zero out bin so we don't double count it.
            # fism_vals[i] = 0.0

    # then go through the ranges
    for iWave, short in enumerate(shorts):
        long = longs[iWave]
        ave_wav[iWave] = (short + long)/2.0
        if (long != short):
            d = np.abs(fism_waves - short)
            iStart = np.argmin(d)
            d = np.abs(fism_waves - long)
            iEnd = np.argmin(d)
            wave_int = 0.0
            for i in range(iStart+1, iEnd+1):
                new_irr[iWave] += fism_vals[i] * \
                    (fism_waves[i+1] - fism_waves[i])
                wave_int += (fism_waves[i+1] - fism_waves[i])
    return new_irr, ave_wav

#------------------------------------------------------------------------------
# main code:
#------------------------------------------------------------------------------

args = parse_args_fism()

if (args.fismfile == 'none'):

    if (args.start == '0'):
        print('Need to specify -start time! Use -h for help!')
        exit()
    
    start = convert_ymdhm_to_dt(args.start)

    if (args.flare):
        print('*** Downloading Flare data - Can only get 24 hours of data! ***')
        end = start + dt.timedelta(days = 1)
    else:
        if (args.end == '0'):
            print('Need to specify -end time! Use -h for help!')
            exit()
        end = convert_ymdhm_to_dt(args.end)

    fism_file = download_fism2(start, end, args.flare)

else:
    fism_file = args.fismfile

fism = read_fism_csv_file(fism_file)
wavelengths = read_euv_csv_file(args.euvfile)

nWaves = len(wavelengths['short'])

filetime = fism["time"][0].strftime('fism%Y%m')
filestart = filetime+'_nWaves_%03d' % nWaves

fig = plt.figure(figsize = (10,10))
ax = fig.add_subplot()

if (args.gitm):
    fileout = filestart + '_gitm.dat'
else:
    fileout = filestart + '.dat'
    
fp = open(fileout, 'wb')

if (args.gitm):
    fp.write('#START\n'.encode())

else:

    shortline = ' 0000 00 00 00 00 00 '
    for short in wavelengths['short']:
        shortline = shortline + "%8.1f" % short
    shortline = shortline + '\n'
    fp.write(shortline.encode())

    longline = ' 0000 00 00 00 00 00 '
    for long in wavelengths['long']:
        longline = longline + "%8.1f" % long
    longline = longline + '\n'
    fp.write(longline.encode())

for iTime, time in enumerate(fism['time']):
    new_irr, ave_wav = rebin_fism(fism['wave'], fism['irr'][iTime], wavelengths)

    if (args.gitm):
        ave_wav = np.flip(ave_wav)
        new_irr = np.flip(new_irr)
    ax.scatter(ave_wav, np.log10(new_irr))

    sTime = time.strftime(' %Y %m %d %H %M %S')
    sData = ' '
    for irr in new_irr:
        sData = sData + "%15.8e" % irr
    line = sTime + sData + '\n'
    fp.write(line.encode())

fp.close()

f107Data = read_ap107_file(local_file = "apf107_temp.txt")
if (len(f107Data['times']) > 0):
    iMid = int(len(fism['time'])/2)
    tMid = fism['time'][iMid]
    tLower = tMid - dt.timedelta(days = 40)
    tUpper = tMid + dt.timedelta(days = 40)

    iLower = find_index_from_time(f107Data, tLower)
    iMid = find_index_from_time(f107Data, tMid)
    iUpper = find_index_from_time(f107Data, tUpper)

    f107 = f107Data['f107'][iMid]
    f107a = np.mean(f107Data['f107'][iLower:iUpper+1])

    euvac = calc_euvac(wavelengths, f107, f107a)

    if (len(euvac) > 0):
        ave_wave = (wavelengths['long'] + wavelengths['short'])/2.0
        ax.plot(ave_wave, np.log10(euvac), label = 'EUVAC')
        ax.legend()

ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('log(W/m2)')

plotfile = filestart + '.png'
print('writing : ',plotfile)    
fig.savefig(plotfile)
plt.close()
