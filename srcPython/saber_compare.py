#!/usr/bin/env python

# Script for comparing GITM outputs to TIMED/SABER data.

# Top-level imports
import numpy as np
from astropy.utils.masked.function_helpers import datetime_as_string
from tensorboard.manager import start
import matplotlib.pyplot as plt
import sys, os, csv
from netCDF4 import Dataset
import xarray
from scipy.interpolate import InterpolatedUnivariateSpline, griddata
from datetime import datetime, timedelta
import pandas as pd
import aacgmv2
import scipy.stats
from netCDF4 import Dataset
import argparse
import matplotlib
matplotlib.use('Qt5Agg')
from scipy.interpolate import interp1d, griddata
from scipy.interpolate import RegularGridInterpolator
import scipy.integrate as integ
from tqdm import tqdm
import pickle
from scipy.signal import savgol_filter
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
import time
import pdb #, kaleido, math, apexpy
from scipy.stats import linregress
# import warnings
# warnings.filterwarnings("error")
import pathlib
import seaborn as sns

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description='Compare GITM NO profiles to SABER data.')
    parser.add_argument('gitm_loc', metavar='gitm_directory', nargs=1, help='location of GITM outputs')
    parser.add_argument('saber_loc', metavar='saber_directory', nargs=1,
                        help='location of SABER observations')
    parser.add_argument('plot_type', metavar='plotting_information', nargs=1, help="The type of plot to "
                                                                                  "produce. Can be 'peak' for a "
                                                                                  "timeseries of peak NO cooling, "
                                                                                  "'integrated' for a timeseries of "
                                                                                  "height-integrated NO cooling, or "
                                                                                  "'both' for generating both.")
    parser.add_argument('override', metavar='overwritting_cmd', nargs=1, help="Boolean argument that controls whether processed data is rewritten. Default value is False.", default=False)
    args = parser.parse_args()
    return args

# ----------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------
def clean_saber(nitric_oxide_profile, altitude_values, region=[100, 250]):
    "Clean SABER data corresponding to a NO volume emission rate altitude profile."
    good_inds_altitude = np.where(altitude_values != -999.)[0]
    altitude_subset = altitude_values[good_inds_altitude]
    nitric_oxide_subset = nitric_oxide_profile[good_inds_altitude]
    good_inds_nitric_oxide = np.where(nitric_oxide_subset > 0)[0]
    altitude_subset_2 = altitude_subset[good_inds_nitric_oxide]
    nitric_oxide_subset_2 = nitric_oxide_subset[good_inds_nitric_oxide]
    # Fit an interpolator:
    interp = InterpolatedUnivariateSpline(altitude_subset_2, nitric_oxide_subset_2, k=3)
    # Application of the interpolation:
    subset_length = len(np.where((altitude_values >= region[0]) & (altitude_values <= region[-1]))[0])
    altitude_range = np.linspace(region[0], region[1], num=subset_length)
    NO_values = interp(altitude_range)
    return NO_values, altitude_range

def gitm_file_mean_profile(filename):
    """
    Reads a single GITM NETCDF file (processed by Aaron Bukowski's code) and returns the mean altitude profile of NO
    emission between 100 and 250 km. DOES NOT use any interpolation to fill in bad values, due to GITM data being fully
    dense.
    Parameters
    ----------
    filename: str
        Name of the GITM NETCDF4 file to read in.
    Returns
    ----------
    NO_prof: array
        Mean NO emission rate values for the entire profile (index in lat and lon).
    mean_altitude_prof: array
        Altitude values for the profile.
    """
    # Read in the GITM file:
    info = Dataset(filename, "r")
    # Compute the averaged NO cooling rate altitude profile:
    mean_altitude_prof = np.asarray(info['alt'])
    NO_cooling = np.asarray(info['NOCooling'])
    # Convert the units of NO_prof from -K/s to W/m^3:
    Cp = np.asarray(info['Cp'])
    Rho = np.asarray(info['Rho'])
    NO_prof = -1 * np.multiply(np.multiply(NO_cooling, Cp), Rho)
    try:
        obs_time = datetime.strptime(info['time'].units[11:].replace(' ', 'T'), "%Y-%m-%dT%H:%M:%S") + timedelta(seconds=int(info['time'][:][0]/86400))
    except:
        obs_time = datetime.strptime(info['time'].units[11:].replace(' ', 'T'), "%Y-%m-%dT%H:%M:%S.%f") + timedelta(
            seconds=int(info['time'][:][0] / 86400))

    alt_inds = np.where((mean_altitude_prof >= 100) & (mean_altitude_prof <= 250))[0]

    return NO_prof[:,:, alt_inds], mean_altitude_prof[alt_inds], obs_time, np.asarray(info['lat']), np.asarray(info['lon'])

def saber_file_mean_profile(filename):
    """
    Reads a single SABER file and returns the mean altitude profile of NO emission between 100 and 250 km. Uses cubic
    interpolation to fill in bad values in the profile.
    Parameters
    ----------
    filename: str
        Name of the SABER NETCDF4 file to read in.
    Returns
    -------
    clean_mean_NO_prof: array
        Mean NO emission rate values for the entire profile.
    clean_mean_altitude_prof: array
        Altitude values for the profile.
    """
    # Read in the SABER file:
    info = Dataset(filename, "r")
    # Compute the averaged NO cooling rate altitude profile (both upper and lower regions):
    mean_altitude_profile = np.mean(info.variables['tpaltitude'], axis=0)
    mean_altitude_profile_top = np.mean(info.variables['tpaltitude_top'], axis=0)
    mean_NO_profile = np.mean(info.variables['NO_ver_unfilt'], axis=0)
    mean_NO_profile_top = np.mean(info.variables['NO_ver_top_unfilt'], axis=0)
    # Combine the profiles together:
    mean_alt_prof = np.flip(np.concatenate((mean_altitude_profile_top, mean_altitude_profile))) # We flip the result, since the altitudes in the file are arranged in descending order
    mean_NO_prof = np.flip(np.concatenate((mean_NO_profile_top, mean_NO_profile)))
    # Clean bad values (eliminate them and use cubic interpolation to fill them in); only consider the region between 100 and 250 km altitude (like Chand, et al. 2024):
    clean_mean_NO_prof, clean_mean_altitude_prof = clean_saber(mean_NO_prof, mean_alt_prof, region=[100,250])

    # Time values for the sounding measurements:
    dates = np.array(info.variables['date']) # Date YYYY DDD
    times = np.array(info.variables['time']) # Time in milliseconds since midnight.
    meanDate = dates[len(dates)//2]
    meanTime = np.nanmean(times)
    observation_time = datetime.strptime(str(meanDate)[:4] + "-" + str(meanDate)[4:], "%Y-%j") + timedelta(seconds=meanTime*0.001)

    # Extract Geolocation Information for the Sounding Measurements:
    mean_tplat = np.nanmean(info.variables['tplatitude'])
    mean_tplat_top = np.nanmean(info.variables['tplatitude_top'])
    mean_tplon = np.nanmean(info.variables['tplongitude'])
    mean_tplon_top = np.nanmean(info.variables['tplongitude_top'])
    obs_lat = np.mean([mean_tplat, mean_tplat_top])
    obs_lon = np.mean([mean_tplon, mean_tplon_top])

    return clean_mean_NO_prof, clean_mean_altitude_prof, observation_time, obs_lat, obs_lon

def saber_file_data(filename, **kwargs):
    """
    Return the (cleaned) NO data contained in a single SABER file.
    """
    # Read in the SABER file:
    info = Dataset(filename, "r")

    # Time values for the sounding measurements:
    dates = np.array(info.variables['date'])  # Date YYYY DDD
    times = np.array(info.variables['time'])  # Time in milliseconds since midnight.
    startDate = datetime.strptime(str(dates[0])[:4] + "-" + str(dates[0])[4:], "%Y-%j") + timedelta(seconds=times[0,0]*0.001)
    # Make a 2D array of datetimes for the observations:
    observation_times = np.zeros(times.shape, dtype=datetime)
    sst = np.zeros(times.shape, dtype=datetime)
    for i in range(times.shape[0]):
        # Rows
        for j in range(times.shape[1]):
            # Columns:
            if times[i,j] == -999.0:
                observation_times[i, j] = np.nan
                sst[i, j] = 0
            else:
                tV = datetime.strptime(str(dates[i])[:4] + "-" + str(dates[i])[4:], "%Y-%j") + timedelta(seconds=times[i,j]*0.001)
                dateRef = kwargs['dateRef']
                upperDateBoundary = dateRef + timedelta(hours=24)
                if tV > upperDateBoundary:
                    observation_times[i, j] = np.nan
                    sst[i, j] = 0
                else:
                    observation_times[i, j] = datetime.strptime(str(dates[i])[:4] + "-" + str(dates[i])[4:], "%Y-%j") + timedelta(seconds=times[i,j]*0.001)
                    sst[i, j] = (observation_times[i, j] - startDate).total_seconds()

    ssc_r = sst.ravel()
    smallDate = np.argmin(ssc_r)
    bigDate = np.argmax(ssc_r)
    smallDate_ = observation_times.ravel()[smallDate]
    bigDate_ = observation_times.ravel()[bigDate]

    # Grab the data:
    tp_alt = info.variables['tpaltitude']
    tp_alt_top = info.variables['tpaltitude_top']
    tp_lat = info.variables['tplatitude']
    tp_lat_top = info.variables['tplatitude_top']
    tp_lon = info.variables['tplongitude']
    tp_lon_top = info.variables['tplongitude_top']
    NO_data = info.variables['NO_ver_unfilt']
    NO_data_top = info.variables['NO_ver_top_unfilt']

    # Arrange the data:
    all_alts = np.hstack((np.array(tp_alt_top), np.array(tp_alt)))
    all_lats = np.hstack((np.array(tp_lat_top), np.array(tp_lat)))
    all_lons = np.hstack((np.array(tp_lon_top), np.array(tp_lon)))
    all_NO_data = np.hstack((np.array(NO_data_top), np.array(NO_data))) * 0.1 # Convert NO cooling from ergs/cm^3/sec to W/m^3

    # Clean the NO data:
    # all_NO_data_clean = interp2d(all_NO_data)

    # Replace bad values with NaNs:
    all_NO_data_na = all_NO_data.copy(); all_NO_data_na[all_NO_data_na == -9.9900002e+01] = np.nan
    all_alts_na = all_alts.copy(); all_alts_na[all_alts_na == -999.0] = np.nan
    all_lats_na = all_lats.copy(); all_lats_na[all_lats_na == -999.0] = np.nan
    all_lons_na = all_lons.copy(); all_lons_na[all_lons_na == -999.0] = np.nan

    # Mask the bad values:
    # all_NO_data_ma = np.ma.masked_where(all_NO_data == -999.0, all_NO_data)
    # all_alts_ma = np.ma.masked_where(all_NO_data == -999.0, all_alts)
    # all_lats_ma = np.ma.masked_where(all_NO_data == -999.0, all_lats)
    # all_lons_ma = np.ma.masked_where(all_NO_data == -999.0, all_lons)
    # all_observation_times_ma = np.ma.masked_where(observation_times == -999.0, observation_times)
    # all_NO_data_ma = np.ma.masked_where(all_NO_data_na == np.nan, all_NO_data_na)
    # all_alts_ma = np.ma.masked_where(all_NO_data_na == np.nan, all_alts_na)
    # all_lats_ma = np.ma.masked_where(all_NO_data_na == np.nan, all_lats_na)
    # all_lons_ma = np.ma.masked_where(all_NO_data_na == np.nan, all_lons_na)
    # all_observation_times_ma = np.ma.masked_where(observation_times == np.nan, observation_times)

    # Truncation step - SABER data often contain values for the next day; EXCLUDE those:
    # if kwargs:
    #     dateRef = kwargs['dateRef']
    #     upperDateBoundary = dateRef + timedelta(hours=24)
    #     ss_dR = (upperDateBoundary - dateRef).total_seconds()
    #     boolean_array = np.zeros_like(sst, dtype=bool)
    #     for row_index in range(sst.shape[0]):
    #         for col_index in range(sst.shape[1]):
    #             if sst[row_index, col_index] > ss_dR:
    #                 boolean_array[row_index, col_index] = False
    #             else:
    #                 boolean_array[row_index, col_index] = True

    # Cleaning data:
    all_NO_data_na_c = interp2d_join(all_NO_data_na, smooth=False)

    # Dealing with unusable data - if interp2d_join returns None, simply return None as well.
    if all_NO_data_na_c is None:
        pdb.set_trace()
        return
    else:
        all_alts_na_c = interp2d_join(all_alts_na)
        all_lats_na_c = interp2d_join(all_lats_na)
        all_lons_na_c = interp2d_join(all_lons_na)
        try:
            observation_times_c = interp2d_join(observation_times, dataType=datetime)
            if observation_times_c is None:
                observation_times_c = interp2d_join(observation_times, dataType=datetime, impute=True)
        except:
            raise Exception
        
        # Return the data:
        return all_alts_na_c, all_lats_na_c, all_lons_na_c, all_NO_data_na_c, observation_times_c, smallDate_, bigDate_

def interp2d_simple(arraylike, k=3, axis=-1):
    """
    Performs (naive) interpolation for 2D data, using the InterpolatedUnivariateSpline class from scipy.
    """
    if axis != -1:
        arraylike = arraylike.T

    # Make a copy of the input data:
    arr = np.zeros_like(arraylike)

    # Loop through each row (column) of the arraylike:
    for i in range(arraylike.shape[0]):
        current_line = arraylike[i, :]
        # If there are NaNs in the line, isolate only those values that are non-NaNs:
        good_indices = ~np.isnan(current_line)
        bad_indices = np.isnan(current_line)
        # Make no changes if ALL the values are good:
        if len(current_line[good_indices]) == len(current_line):
            arr[i, :] = current_line
        else:
            # Otherwise, proceed with interpolation:
            if len(current_line[good_indices]) <= k:
                # For lines that are ALL NaNs, do nothing (TODO: replace with bette strategy)
                arr[i, :] = current_line
            else:
                # Construct an x-axis at the valid points:
                all_locs = np.asarray([int(element) for element in np.linspace(0,len(current_line)-1, len(current_line))])
                good_locs = all_locs[good_indices]
                bad_locs = all_locs[bad_indices]
                # Fit a spline between the x-axis and the values at the valid points:
                spl = InterpolatedUnivariateSpline(good_locs, current_line[good_indices], k=k)
                # Using the spline, determine the missing data:
                out = spl(bad_locs)
                # Update the array:
                arr[i, :] = current_line
                j = 0
                for idx in bad_locs:
                    # If the interpolated value is larger than the peak cooling in the given profile, replace it with the
                    # MINIMUM cooling value in that given profile:
                    max_cooling = np.nanmax(arr[:, idx])
                    if out[j] >= max_cooling:
                        arr[i, idx] = np.nanmin(arr[:, idx])
                    else:
                        arr[i, idx] = out[j]
                    j += 1

    # Sanity Check - plot the before and after images on the same color scale:
    # g_vmin = np.nanmin([np.nanmin(arraylike), np.nanmin(arr)])
    # g_vmax = np.nanmax([np.nanmax(arraylike), np.nanmax(arr)])
    # fig, axs = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True)
    # axs[0].imshow(arraylike, aspect='auto', vmin=g_vmin, vmax=g_vmax, cmap='viridis')
    # pos = axs[1].imshow(arr, aspect='auto', vmin=g_vmin, vmax=g_vmax, cmap='viridis')
    # fig.colorbar(pos, ax=axs[1])

    return arr

def interp2d_join(arraylike, smooth=True, dataType=None, impute=False):
    # If dataType = datetime.datetime, convert everything to seconds:
    if dataType:
        new_array_like = np.zeros_like(arraylike, dtype=float)
        baseDatetime = arraylike[0, 0]
        for i in range(new_array_like.shape[0]):
            for j in range(new_array_like.shape[1]):
                try:
                    new_array_like[i, j] = (baseDatetime - arraylike[i, j]).total_seconds()
                except:
                    new_array_like[i, j] = np.nan
        arraylike = new_array_like

    if dataType and impute == True:
        # Fill in the array with interpolated data, uniformly between the smallest and largest time values:
        lowerBoundary = np.nanmin(arraylike)
        upperBoundary = np.nanmax(arraylike)
        imputedValues = np.linspace(lowerBoundary, upperBoundary, int(2*arraylike.size))
        arr = np.reshape(imputedValues, (arraylike.shape[0], int(2*arraylike.shape[1])))
        # Note: This approach assumes that the starting and ending date are VERY CLOSE (within minutes) of each other!
    else:
        # Make a copy of the input data:
        arr = np.zeros_like(arraylike)
        
        # STANDARD PROCESSING
        # Loop through each row of the arraylike:
        for i in range(arraylike.shape[0]):
            current_row = arraylike[i, :]
            # Clip all NaN values:
            current_row_clipped = current_row[~np.isnan(current_row)]
            # Smooth the data in a 50 km window:
            if smooth:
                # df = pd.Series(current_row_clipped)
                # c_y_ma = df.rolling(200, center=True).mean()
                try:
                    c_y_ma = smart_roll(current_row_clipped, 200)
                except:
                    # If this situation obtains, it is because the given row is ALL NaNs. We assume that there exist
                    # preceding rows not full of NaNs, and just linearly interpolate from those rows to the current one:
                    current_row_extrapolated = np.zeros_like(arr[0, :])
                    for j in range(arr.shape[1]):
                        xvals = np.linspace(1, len(arr[:i-1, j]), num=len(arr[:i-1, j]))
                        try:
                            result = linregress(xvals, arr[:i-1, j])
                        except:
                            # In this case, there are likely so many NaN values that the problem is likely intractible. We check if the number of NaN values exceeds half the number of values in the array. If so, processing this array is abandoned.
                            if arraylike[np.isnan(arraylike)].size >= int(arraylike.size/2):
                                print('More than half the values in this data are bad. Abandoning processing...')
                                return
                            else:
                                pdb.set_trace()
                        newResult =  result.intercept + result.slope*(xvals[-1] + 1)
                        current_row_extrapolated[j] = newResult
                    # Then we perform the desired operation:
                    c_y_ma = smart_roll(current_row_extrapolated, 200)
                i_d = c_y_ma
            else:
                i_d = current_row_clipped
            # Extend the size of the data to 800 datapoints, using interp/extrap:
            x = np.arange(len(i_d))
            try:
                interpolator = interp1d(
                    x, i_d,
                    kind="linear",
                    fill_value="extrapolate",
                    assume_sorted=True,
                )
            except:
                raise Exception
            # Perform the extrapolation:
            if dataType:
                x_new = np.linspace(0, 399, 400)
            else:
                x_new = np.linspace(0, 799, 800)
            new_row = interpolator(x_new)
            # Update the array:
            arr[i, :] = new_row

    if dataType:
        final_arraylike = np.zeros_like(arraylike, dtype=object)
        for i in range(final_arraylike.shape[0]):
            for j in range(final_arraylike.shape[1]):
                final_arraylike[i, j] = baseDatetime + timedelta(seconds=arr[i, j])
        arr = final_arraylike

    return arr

def smart_roll(u_d_arr, windowSize):
    """
    Computes rolling average (arithmetic mean) for some 1D array or list and a specified windowsize. Performs
    interpolation on the leading and trailing edge by default.
    """

    # Initial roll:
    df = pd.Series(u_d_arr)
    rolled = df.rolling(windowSize, center=True).mean()
    # Leading window:
    x = np.array([0, 1])
    leadingInterpolator = interp1d(
            x, u_d_arr[[0, windowSize]], kind="linear", assume_sorted=True
    )
    x_new = np.linspace(0, 1, windowSize)
    leadingInterpVals = leadingInterpolator(x_new)
    rolled[:windowSize] = leadingInterpVals
    # Trailing window:
    trailingInterpolator = interp1d(
        x, u_d_arr[[-windowSize, -1]], kind="linear", assume_sorted=True
    )
    trailingInterpVals = trailingInterpolator(x_new)
    rolled[-windowSize:] = trailingInterpVals

    return rolled

def movingAverage(arr, windowSize):
    """
    Compute a moving average along a 1D array with a given window size.
    """
    cumsum, moving_aves = ([0]
                           , [])
    for i, x in enumerate(arr, 1):
        cumsum.append(cumsum[i - 1] + x)
        if i >= windowSize:
            moving_ave = (cumsum[i] - cumsum[i - windowSize]) / windowSize
            # can do stuff with moving_ave here
            moving_aves.append(moving_ave)
    return moving_aves

def interp2d(arr, badVal=-999.0, smooth=True):
    """
    Wrapper for interp1d used to interpolate 2D data to fill NaNs. Proceeds one row at a time and interpolates values
    in each row; by default, extrapolates data outside the data range.
    """
    arr = arr.copy()
    if badVal:
        arr[arr==badVal] = np.nan

    def interp_nans(y, x=None):
        if x is None:
            x = np.arange(len(y))
        nans = np.isnan(y)
        interpolator = interp1d(
            x[~nans],
            y[~nans],
            kind="linear",
            fill_value="extrapolate",
            assume_sorted=True,
        )
        return interpolator(x)

    for i in range(arr.shape[0]):
        current_y = arr[i, :]
        if smooth:
            # Take moving average in a 50 km window
            df = pd.Series(current_y)
            c_y_ma = df.rolling(200, center=True).mean()
            i_d = c_y_ma
        else:
            i_d = current_y
        # Interpolate over the (smoothed) data
        result2 = interp_nans(i_d)
        # Replace the bad values with the good ones:
        current_y[np.isnan(current_y)] = result2[np.isnan(current_y)]
        arr[i, :] = current_y

    return arr

def gitmPlot(thermData, dateStr, latStr):
    """
    Helper function to generate a NO altitude profile from GITM 3DTHM data.
    :param thermData:
        3DTHM altitude profile data collected by the function 'gitm_saber_profile_comparison'.
    :param dateStr:
        A string describing the dates over which the profiles are being averaged.
    :param latStr: str
        A string describing what latitudinal region the data corresponds to.
    :return meanProfile_u: numpy.ndarray
        The mean altitude profile in the region under consideration.
    :return meanAltitudes: numpy.ndarray
        The mean altitudes to which the profile corresponds.
    """
    profs = np.asarray([element[0] for element in thermData])
    mean_prof = np.mean(profs, axis=0)
    # stddev = np.std(profs, axis=0)
    alts = thermData[0][1]

    # Compute 95%CI bounds:
    def conf_int(profile_data, confidence=0.95):
        lb = []
        ub = []
        for i in range(profile_data.shape[1]):
            a = 1.0 * np.array(profile_data[:, i])
            n = len(a)
            m, se = np.mean(a), scipy.stats.sem(a)
            h = se * scipy.stats.t.ppf((1 + confidence) / 2., n - 1)
            lb.append(m - h)
            ub.append(m + h)
        return np.asarray(lb), np.asarray(ub)

    profile_lb, profile_ub = conf_int(profs)

    # Convert the emission rates from ergs/cm^3/s to nW/m^3:
    meanProfile_u = mean_prof  #/ (1e-6) / (1e7) * -1e9
    profile_lb_u = profile_lb #/ (1e-6) / (1e7) * -1e9
    profile_ub_u = profile_ub #/ (1e-6) / (1e7) * -1e9
    plt.figure()
    plt.title(r'GITM NO Cooling (' + dateStr + '): ' + latStr)
    plt.fill_betweenx(alts, profile_lb_u, profile_ub_u, where=profile_ub_u > profile_lb_u, color='b',
                      alpha=0.3)
    plt.plot(meanProfile_u, alts, color='b')
    plt.xlabel('NO Cooling (W/m$^3$)')  # (ergs/cm$^3$/s)')
    plt.ylabel('Altitude (km)')
    plt.close()

    return meanProfile_u, profile_lb_u, profile_ub_u, alts

def saberPlot(soundingData, dateStr, latStr):
    """
    Helper function to generate a NO altitude profile from SABER sounding data.
    :param soundingData:
        Sounding data collected by the function 'gitm_saber_profile_comparison'.
    :param dateStr:
        A string describing the dates over which the profiles are being averaged.
    :param latStr: str
        A string describing what latitudinal region the data corresponds to.
    :return meanProfile_u: numpy.ndarray
        The mean altitude profile in the region under consideration.
    :return meanAltitudes: numpy.ndarray
        The mean altitudes to which the profile corresponds.
    """
    profs = [element[0] for element in soundingData]
    alts = [element[1] for element in soundingData]
    elementLength = min([len(element) for element in profs])
    profs_S = np.array([element[:elementLength] for element in profs])
    alts_S = np.array([element[:elementLength] for element in alts])
    meanProfile = np.mean(profs_S, axis=0)
    meanAltitudes = np.mean(alts_S, axis=0)
    profile_lb, profile_ub = conf_int(profs_S)
    # Convert the emission rates from ergs/cm^3/s to nW/m^3:
    meanProfile_u = meanProfile * 0.1 # / (1e-6) / (1e7) * 1e9
    profile_lb_u = profile_lb * 0.1 # / (1e-6) / (1e7) * 1e9
    profile_ub_u = profile_ub * 0.1 # / (1e-6) / (1e7) * 1e9
    #
    plt.figure()
    plt.title(r'SABER NO Cooling (' + dateStr + '): '+latStr)
    plt.fill_betweenx(meanAltitudes, profile_lb_u, profile_ub_u, where=profile_ub_u>profile_lb_u, color='b', alpha=0.3)
    plt.plot(meanProfile_u, meanAltitudes, color='b')
    plt.xlabel('NO Cooling (W/m$^3$)')  # (ergs/cm$^3$/s)')
    plt.ylabel('Altitude (km)')
    plt.close()
    return meanProfile_u, profile_lb_u,  profile_ub_u, meanAltitudes

def compPlot_ci(saberSoundingData, gitmSoundingData, dateStr, latStr):
    """
    Does the same thing as 'saberPlot' and 'gitmPlot', but takes BOTH SABER and GITM data and plots the mean profiles,
    along with confidence bands.
    """
    # SABER
    profs_saber = [element[0] for element in saberSoundingData]
    alts_saber = [element[1] for element in saberSoundingData]
    elementLength_saber = min([len(element) for element in profs_saber])
    profs_S = np.array([element[:elementLength_saber] for element in profs_saber])
    alts_S = np.array([element[:elementLength_saber] for element in alts_saber])
    meanProfile_saber = np.mean(profs_S, axis=0)
    meanAltitudes_saber = np.mean(alts_S, axis=0)
    profile_lb, profile_ub = conf_int(profs_S)
    # Convert the emission rates from ergs/cm^3/s to nW/m^3:
    meanProfile_u_saber = meanProfile_saber * 0.1  # / (1e-6) / (1e7) * 1e9
    profile_lb_u_saber = profile_lb * 0.1  # / (1e-6) / (1e7) * 1e9
    profile_ub_u_saber = profile_ub * 0.1  # / (1e-6) / (1e7) * 1e9

    # GITM
    profs_gitm = [element[0] for element in gitmSoundingData]
    alts_gitm = [element[1] for element in gitmSoundingData]
    elementLength_gitm = min([len(element) for element in profs_gitm])
    profs_G = np.array([element[:elementLength_gitm] for element in profs_gitm])
    alts_G = np.array([element[:elementLength_gitm] for element in alts_gitm])
    meanProfile_gitm = np.nanmean(profs_G, axis=0)
    meanAltitudes_gitm = np.nanmean(alts_G, axis=0)
    profile_lb_gitm, profile_ub_gitm = conf_int(profs_G)

    plt.figure()
    plt.title(r'NO Cooling (' + dateStr + '): ' + latStr)
    plt.fill_betweenx(meanAltitudes_saber, profile_lb_u_saber, profile_ub_u_saber, where=profile_ub_u_saber > profile_lb_u_saber, color='b',
                      alpha=0.3)
    plt.plot(meanProfile_u_saber, meanAltitudes_saber, color='b')
    #
    plt.fill_betweenx(meanAltitudes_gitm, profile_lb_gitm, profile_ub_gitm,
                      where=profile_ub_gitm > profile_lb_gitm, color='r',
                      alpha=0.3)
    plt.plot(meanProfile_gitm, meanAltitudes_gitm, color='r')
    plt.xlabel('NO Cooling (W/m$^3$)')  # (ergs/cm$^3$/s)')
    plt.ylabel('Altitude (km)')
    plt.close()
    # TODO: Save the figure...

    return

def compPlot(gitmProfile, saberProfile, dateStr, latStr, latLabel):
    """
    Helper function for taking outputs from gitmPlot and saberPlot and comparing the profiles together in the same
    figure. This function return nothing but rather saves a figure to a default location.
    :param gitmProfile: list
        Output from gitmPlot; the outputs from that function are arranged into a list: [meanProfile_u, profile_lb_u,
        profile_ub_u, alts]
    :param saberProfile: list
        Output from saberPlot; the outputs from that function are arranged into a list: [meanProfile_u, profile_lb_u,
        profile_ub_u, meanAltitudes]
    :param dateStr:
        A string describing the dates over which the profiles are being averaged.
    :param latStr: str
        A string describing what latitudinal region the data corresponds to.
    :param latLabel: str
        A label describing the latitudinal region (i.e. 'Polar', 'Auroral', 'Midlatitude', 'Equatorial')
    """
    plt.figure()
    plt.title('NO Cooling ('+dateStr+'): '+latStr)
    # SABER
    plt.fill_betweenx(saberProfile[-1], saberProfile[1], saberProfile[2], where=saberProfile[2]>saberProfile[1],
                      color='b', alpha=0.3)
    plt.plot(saberProfile[0], saberProfile[-1], color='b', label='SABER')
    # GITM
    plt.fill_betweenx(gitmProfile[-1], gitmProfile[1], gitmProfile[2], where=gitmProfile[2]>gitmProfile[1], color='r', alpha=0.3)
    plt.plot(gitmProfile[0], gitmProfile[-1], color='r', label='GITM')
    # Labels and Axes
    plt.legend(loc='best')
    plt.xlabel('NO Cooling (W/m$^3$)')
    plt.ylabel('Altitude (km)')
    # Save the figure:
    date_str = dateStr[:10]
    plt.savefig('gitmSaber_NO_compare_'+date_str+'_'+latLabel+'.png', dpi=300)
    plt.close()

    return

def gitm_saber_profile_comparison(gitm_data_dir, saveLoc):
    """
    Given the location of where GITM run data has been stored for an arbitrary GITM simulation, (1) download TIMED/SABER
    data during the same time period and (2) collect NO volume emission rates between GITM and SABER during the same
    time period. The data can then be broken down into comparisons by latitude region and is output to a file so that
    results between multiple runs can be compared.
    :param gitm_data_dir: str
        Location of where GITM simulation results are located (i.e. 'v10x20/UA/data')
    :param saveLoc: str
        The location where SABER data AND analysis output files will be placed.
    :return saber_NO: numpy.ndarray
        An array containing ALL the TIMED/SABER mean daily altitude profiles of NO emissions (per sounding).
        All soundings taken by SABER over the interval in question are contained in this matrix.
    :return saber_LOCS: numpy.ndarray
        An array of the same shape (minus the altitudinal dimension) as 'saber_NO' but with geolocation information
        for each sounding.
    :return saber_TIMES: numpy.ndarray
        An array with time values for the SABER soundings.
    :return gitm_NO: numpy.ndarray
        Contains the same NO emission profile information as 'saber_NO' but corresponding to GITM.
    :return gitm_LOCS: numpy.ndarray
        Contains the geolocation information for the NO emission profiles in GITM.
    :return gitm_TIMES: numpy.ndarray
        Contains the times corresponding to each GITM NO emission profile.
    """
    # Grab the names of all the files in the GITM data directory :
    filenames = os.listdir(gitm_data_dir)
    validFilenames = []
    datetime_strs = []
    datetime_stamps = []
    for file_ in filenames:
        if file_.endswith('.nc'):
            validFilenames.append(file_)
            datetime_strs.append( '20'+file_[1:3]+'-'+file_[3:5]+'-'+file_[5:7]+'T'+file_[8:10]+':'+file_[10:12]+':'+file_[12:14] )
            datetime_stamps.append( datetime.strptime(datetime_strs[-1], "%Y-%m-%dT%H:%M:%S" ) )
    # Arrange the files in chronological order:
    ordering = np.argsort(np.asarray(datetime_stamps))
    orderedValidFilenames = np.asarray(validFilenames)[ordering]
    orderedDatetimeStrs = np.asarray(datetime_strs)[ordering]
    orderedDatetimeStamps = np.asarray(datetime_stamps)[ordering]

    # From the list of filenames, grab the starting and ending dates:
    # startingFile = orderedValidFilenames[0]
    # endingFile = orderedValidFilenames[-1]
    startingTimeStr = orderedDatetimeStrs[0]
    endingTimeStr =orderedDatetimeStrs[-1]
    startingDatetime = orderedDatetimeStamps[0] # datetime.strptime('20' + startingTimeStr[1:3] + "-" + startingTimeStr[3:5] + "-" + startingTimeStr[5:7], "%Y-%m-%d") + timedelta(hours=int(startingTimeStr[8:10])) + timedelta(minutes=int(startingTimeStr[10:12]))
    endingDatetime = orderedDatetimeStamps[-1] # datetime.strptime('20' + endingTimeStr[1:3] + "-" + endingTimeStr[3:5] + "-" + endingTimeStr[5:7], "%Y-%m-%d") + timedelta(hours=int(endingTimeStr[8:10])) + timedelta(minutes=int(endingTimeStr[10:12]))
    if datetime.strftime(startingDatetime, "%Y-%m-%dT%H:%M:%S") == datetime.strftime(endingDatetime, "%Y-%m-%dT%H:%M:%S"):
        titleStr = datetime.strftime(startingDatetime, "%Y-%m-%dT%H:%M:%S")
    else:
        titleStr = datetime.strftime(startingDatetime, "%Y-%m-%dT%H:%M:%S") + ' to ' + datetime.strftime(endingDatetime, "%Y-%m-%dT%H:%M:%S")

    gitmData = []
    id = 0
    for gitm_file in orderedValidFilenames:
        current_NO_prof, current_altitude_prof, current_obs_time, current_obs_lat, current_obs_lon = gitm_file_mean_profile(
            gitm_data_dir + gitm_file)
        gitmData.append(
            [current_NO_prof, current_altitude_prof, current_obs_time, current_obs_lat, current_obs_lon, current_obs_lat,
             current_obs_lon])
        id += 1

    # Compute the GITM profiles in latitudinal chunks (using the boundaries from Brandt and Ridley 2022 - doi:10.1029/2021SW002922):
    print('Processing GITM data to obtain NO cooling profiles...')
    # 1: Polar Region |MLAT| >= 80
    polar_inds = np.where(gitmData[0][-2] >= 80)[0]
    polar_profs = []
    for gitmProfData in gitmData:
        polar_profs.append( [np.mean(np.mean(gitmProfData[0][:, polar_inds, :], axis=0), axis=0), current_altitude_prof] )
    latStrPol = '$|\mathrm{MLAT}| \geq 80^{\degree}$'
    meanPolarProfile_gitm, meanPolarProfile_gitm_lb, meanPolarProfile_gitm_ub, meanPolarAltitudes_gitm = gitmPlot(polar_profs, titleStr, latStrPol)
    writeCSV(meanPolarProfile_gitm, meanPolarProfile_gitm_lb, meanPolarProfile_gitm_ub, meanPolarAltitudes_gitm,
             source='GITM', label='Polar_' + titleStr, location=saveLoc[:-6]+'gitm/')

    # 2: Auroral Region 60 <= |MLAT| < 80
    auroral_inds = np.where((gitmData[0][-2] >= 60) & (gitmData[0][-2] < 80))[0]
    auroral_profs = []
    for gitmProfData in gitmData:
        auroral_profs.append([np.mean(np.mean(gitmProfData[0][:, auroral_inds, :], axis=0), axis=0), current_altitude_prof])
    latStrAur = '$60 \leq |\mathrm{MLAT}| < 80^{\degree}$'
    meanAuroralProfile_gitm, meanAuroralProfile_gitm_lb, meanAuroralProfile_gitm_ub, meanAuroralAltitudes_gitm = gitmPlot(
        auroral_profs, titleStr, latStrAur)
    writeCSV(meanAuroralProfile_gitm, meanAuroralProfile_gitm_lb, meanAuroralProfile_gitm_ub, meanAuroralAltitudes_gitm,
             source='GITM', label='Auroral_' + titleStr, location=saveLoc[:-6] + 'gitm/')

    # 3: Midlatitudes 30 <= |MLAT| < 60
    midlat_inds = np.where((gitmData[0][-2] >= 30) & (gitmData[0][-2] < 60))[0]
    midlat_profs = []
    for gitmProfData in gitmData:
        midlat_profs.append(
            [np.mean(np.mean(gitmProfData[0][:, midlat_inds, :], axis=0), axis=0), current_altitude_prof])
    latStrMid = '$30 \leq |\mathrm{MLAT}| < 60^{\degree}$'
    meanMidlatProfile_gitm, meanMidlatProfile_gitm_lb, meanMidlatProfile_gitm_ub, meanMidlatAltitudes_gitm = gitmPlot(
        midlat_profs, titleStr, latStrMid)
    writeCSV(meanMidlatProfile_gitm, meanMidlatProfile_gitm_lb, meanMidlatProfile_gitm_ub, meanMidlatAltitudes_gitm,
             source='GITM', label='Midlatitudes_' + titleStr, location=saveLoc[:-6] + 'gitm/')

    # 4: Low Latitudes / Equatorial region |MLAT| < 30
    eq_inds = np.where(gitmData[0][-2] < 30)[0]
    eq_profs = []
    for gitmProfData in gitmData:
        eq_profs.append(
            [np.mean(np.mean(gitmProfData[0][:, eq_inds, :], axis=0), axis=0), current_altitude_prof])
    latStrEq = '$|\mathrm{MLAT}| \leq 30^{\degree}$'
    meanEqProfile_gitm, meanEqProfile_gitm_lb, meanEqProfile_gitm_ub, meanEqAltitudes_gitm = gitmPlot(
        eq_profs, titleStr, latStrEq)
    writeCSV(meanEqProfile_gitm, meanEqProfile_gitm_lb, meanEqProfile_gitm_ub, meanEqAltitudes_gitm,
             source='GITM', label='Equatorial_' + titleStr, location=saveLoc[:-6] + 'gitm/')

    ####################################################################################################################

    # Download saber data between the starting and ending dates (inclusive):
    pdb.set_trace()
    here = pathlib.Path(__file__).parent.resolve()
    saber_download_file_str = '../srcPython/saber_download.py'
    saber_download_file = here.joinpath(saber_download_file_str)
    pdb.set_trace()
    cmd = 'python '+str(saber_download_file)+' '+str(startingDatetime)[:10].replace('-', '')+' '+str(endingDatetime)[:10].replace('-', '')+' '+saveLoc
    os.system(cmd)

    # Read in the SABER data and generate the SABER data matrix and all other associated information:
    dateList = np.asarray(pd.date_range(startingTimeStr[:-9], endingTimeStr[:-9], periods=(endingDatetime - startingDatetime).days+1).to_pydatetime())
    # np.asarray(pd.date_range('20' + startingTimeStr[1:3] + "-" + startingTimeStr[3:5] + "-" + startingTimeStr[5:7],
    #                          '20' + endingTimeStr[1:3] + "-" + endingTimeStr[3:5] + "-" + endingTimeStr[5:7], periods=(endingDatetime - startingDatetime).days+1).to_pydatetime())
    saberData = []
    for date in dateList:
        current_directory = saveLoc + str(date.year) + '/' + dateList[0].strftime('%j')
        current_saberfiles = os.listdir(current_directory)
        for saberfile in current_saberfiles:
            current_NO_prof, current_altitude_prof, current_obs_time, current_obs_lat, current_obs_lon = saber_file_mean_profile(current_directory+'/'+saberfile)
            out_lat, out_lon, out_r = aacgmv2.wrapper.convert_latlon(current_obs_lat, current_obs_lon, current_altitude_prof[0], current_obs_time) # TODO: DETERMINE IF AACGM IS ACTUALLY NEEDED
            saberData.append([current_NO_prof, current_altitude_prof, current_obs_time, current_obs_lat, current_obs_lon, out_lat, out_lon, out_r])

    # Compute the SABER profiles in latitudinal chunks (using the boundaries from Brandt and Ridley 2022 - doi:10.1029/2021SW002922):
    # 1: Polar Region |MLAT| >= 80
    polarSoundings = [element for element in saberData if abs(element[-3])>=80]

    # 2: Auroral Region 60 <= |MLAT| < 80
    auroralSoundings = [element for element in saberData if ((abs(element[-3]) >= 60) & (abs(element[-3]) < 80))]

    # 3: Midlatitudes 30 <= |MLAT| < 60
    midlatitudeSoundings = [element for element in saberData if ((abs(element[-3]) >= 30) & (abs(element[-3]) < 60))]

    # 4: Low Latitudes / Equatorial region |MLAT| < 30
    equatorialSoundings = [element for element in saberData if abs(element[-3]) < 30]

    # Plots the SABER profiles (with confidence bounds)!!
    # then, output the saber profile data (by latitudinal regions) individual .CSV files for each region:
    if len(polarSoundings) > 0:
        latStrPol = '$|\mathrm{MLAT}| \geq 80^{\degree}$'
        meanPolarProfile, meanPolarProfile_lb, meanPolarProfile_ub, meanPolarAltitudes = saberPlot(polarSoundings, titleStr, latStrPol)
        writeCSV(meanPolarProfile, meanPolarProfile_lb, meanPolarProfile_ub, meanPolarAltitudes, source='SABER', label='Polar_'+titleStr, location=saveLoc)
    if len(auroralSoundings) > 0:
        latStrAur = '$60 \leq |\mathrm{MLAT}| < 80^{\degree}$'
        meanAuroralProfile, meanAuroralProfile_lb, meanAuroralProfile_ub, meanAuroralAltitudes = saberPlot(auroralSoundings, titleStr, latStrAur)
        writeCSV(meanAuroralProfile, meanAuroralProfile_lb, meanAuroralProfile_ub, meanAuroralAltitudes,
                 source='SABER', label='Auroral_'+titleStr, location=saveLoc)
    if len(midlatitudeSoundings) > 0:
        latStrMid = '$30 \leq |\mathrm{MLAT}| < 60^{\degree}$'
        meanMidlatitudeProfile, meanMidlatitudeProfile_lb, meanMidlatitudeProfile_ub, meanMidlatitudeAltitudes = saberPlot(midlatitudeSoundings, titleStr, latStrMid)
        writeCSV(meanMidlatitudeProfile, meanMidlatitudeProfile_lb, meanMidlatitudeProfile_ub, meanMidlatitudeAltitudes,
                 source='SABER', label='Midlatitudes_'+titleStr, location=saveLoc)
    if len(equatorialSoundings) > 0:
        latStrEq = '$|\mathrm{MLAT}| \leq 30^{\degree}$'
        meanEqProfile, meanEqProfile_lb, meanEqProfile_ub, meanEqAltitudes = saberPlot(equatorialSoundings, titleStr, latStrEq)
        writeCSV(meanEqProfile, meanEqProfile_lb, meanEqProfile_ub, meanEqAltitudes, source='SABER', label='Equatorial_'+titleStr, location=saveLoc)

    ####################################################################################################################

    # Joint visualization; plots with GITM/SABER side-by-side or overlaid, in each latitudinal region (these plots are saved):
    if len(polarSoundings) > 0:
        compPlot(
            [meanPolarProfile_gitm, meanPolarProfile_gitm_lb, meanPolarProfile_gitm_ub, meanPolarAltitudes_gitm],
            [meanPolarProfile, meanPolarProfile_lb, meanPolarProfile_ub, meanPolarAltitudes],
            titleStr, latStrPol, 'Polar')
    if len(auroralSoundings) > 0:
        compPlot(
            [meanAuroralProfile_gitm, meanAuroralProfile_gitm_lb, meanAuroralProfile_gitm_ub, meanAuroralAltitudes_gitm],
            [meanAuroralProfile, meanAuroralProfile_lb, meanAuroralProfile_ub, meanAuroralAltitudes],
            titleStr, latStrAur, 'Auroral')
    if len(midlatitudeSoundings) > 0:
        compPlot([meanMidlatProfile_gitm, meanMidlatProfile_gitm_lb, meanMidlatProfile_gitm_ub, meanMidlatAltitudes_gitm],
                 [meanMidlatitudeProfile, meanMidlatitudeProfile_lb, meanMidlatitudeProfile_ub, meanMidlatitudeAltitudes], titleStr, latStrMid, 'Midlatitudes')
    if len(equatorialSoundings) > 0:
        compPlot([meanEqProfile_gitm, meanEqProfile_gitm_lb, meanEqProfile_gitm_ub, meanEqAltitudes_gitm],
                 [meanEqProfile, meanEqProfile_lb, meanEqProfile_ub, meanEqAltitudes], titleStr, latStrEq, 'Equatorial')

    return

def extract_feature(data, time_array, feature='peak', axis=0):
    """
    Given some 2D array, either extract the maximum value along a desired dimension of that array (and return the
    indices of those values), or integrate the quantities in the 2D array along the desired dimension and return the
    result (along with the index corresponding to the centroid of the data that is integrated along each dimension).
    Note that if feature = 'integrated', then 'data' should have two entries; the first is the data to be integrated
    and the second are the 'x-values' (grid points) which are used to determine bounds of integration.
    """
    stacked_time_array = np.hstack((time_array, time_array))
    if feature == 'integral':
        if len(data) != 2:
            raise ValueError("If determining the integral, 'Data' should have TWO elements!")
        else:
            if axis != 0:
                raise ValueError("If determining the integral, the axis must equal 0, and the second dimension of the"
                                 "the first element of 'data' should refer to altitude.")
            else:
                limits = [100, 250] # Defines the bounds of integration to between 100 and 250 km.
                values = data[0]
                alt_grid = data[1]
                # Subsetting
                mean_alt_profile = np.nanmean(alt_grid, axis=0)
                subset_inds = np.where( (mean_alt_profile >= limits[0]) & (mean_alt_profile <= limits[-1]) )[0]
                subset_values = values[:, subset_inds]
                subset_alt_grid = alt_grid[:, subset_inds]
                # Smooth the values in a 5 km window:
                cleaned_vals, cleaned_vals_low_CI, cleaned_vals_high_CI = smooth_array(subset_values, subset_alt_grid, resolution=5)
                # Integrate:
                vals = np.zeros(values.shape[0])
                times = np.zeros(values.shape[0], dtype=object)
                inds = np.full_like(vals, fill_value=subset_inds[len(subset_inds)//4])
                for i in range(len(vals)):
                    out = np.nansum(cleaned_vals[i, :]) * ( (subset_alt_grid[0, :][0] - subset_alt_grid[0, :][-1]) * 1e3 ) # W/m^2
                    vals[i] = out
                    times[i] = stacked_time_array[i, int(inds[i])]
    else:
        inds = np.zeros(data.shape[axis])
        vals = np.zeros(data.shape[axis])
        times = np.zeros(data.shape[axis], dtype=object)
        for i in range(data.shape[axis]):
            try:
                max_ind = np.nanargmax(data[i, :500])
            except:
                max_ind = 0
            inds[i] = 100 # 50 # max_ind
            vals[i] = data[i, max_ind]
            try:
                times[i] = stacked_time_array[i, max_ind]
            except:
                raise Exception

    return vals, inds, times

def gitm_saber_direct_compare(gitm_data_dir, saveLoc, contour=False, plot_type=None, override=False):
    """
    Does the same thing as gitm_saber_profile_comparison but instead extracts GITM data at locations identical to those
    corresponding to the SABER data. Uses 3D interpolation to perform the extrapolation if SABER locations are
    sufficiently far away from GITM locations.
    :param gitm_data_dir: str
        Location of where GITM simulation results are located (i.e. 'v10x20/out/')
    :param saveLoc: str
        The location where SABER data AND analysis output files will be placed.
    :param contour: bool
        Indicates whether altitude vs. time contour plots of NO are generated.
    :param plot_type: str
        Can either be 'peak' for which the peak cooling only is compared or 'integrated' for which the height-integrated
        cooling is compared. Default is None, for which comparisons across the entire altitude profiles of NO will be
        compared between GITM and SABER.
    :param override: bool
        Determines if pre-existing processed data is overriden. Default is False.
    """
    # Grab the names of all the files in the GITM data directory :
    filenames = os.listdir(gitm_data_dir)
    validFilenames = []
    datetime_strs = []
    datetime_stamps = []
    for file_ in filenames:
        if file_.endswith('.nc'):
            validFilenames.append(file_)
            datetime_strs.append('20' + file_[1:3] + '-' + file_[3:5] + '-' + file_[5:7] + 'T' + file_[8:10] + ':' + file_[
                                                                                                                 10:12] + ':' + file_[
                                                                                                                                12:14])
            datetime_stamps.append(datetime.strptime(datetime_strs[-1], "%Y-%m-%dT%H:%M:%S"))
    # Arrange the files in chronological order:
    ordering = np.argsort(np.asarray(datetime_stamps))
    orderedValidFilenames = np.asarray(validFilenames)[ordering]
    orderedDatetimeStrs = np.asarray(datetime_strs)[ordering]
    orderedDatetimeStamps = np.asarray(datetime_stamps)[ordering]
    
    # From the list of filenames, grab the starting and ending dates:
    # startingFile = orderedValidFilenames[0]
    # endingFile = orderedValidFilenames[-1]
    startingTimeStr = orderedDatetimeStrs[0]
    endingTimeStr = orderedDatetimeStrs[-1]
    startingDatetime = orderedDatetimeStamps[
        0]  # datetime.strptime('20' + startingTimeStr[1:3] + "-" + startingTimeStr[3:5] + "-" + startingTimeStr[5:7], "%Y-%m-%d") + timedelta(hours=int(startingTimeStr[8:10])) + timedelta(minutes=int(startingTimeStr[10:12]))
    endingDatetime = orderedDatetimeStamps[
        -1]  # datetime.strptime('20' + endingTimeStr[1:3] + "-" + endingTimeStr[3:5] + "-" + endingTimeStr[5:7], "%Y-%m-%d") + timedelta(hours=int(endingTimeStr[8:10])) + timedelta(minutes=int(endingTimeStr[10:12]))
    if datetime.strftime(startingDatetime, "%Y-%m-%dT%H:%M:%S") == datetime.strftime(endingDatetime,
                                                                                     "%Y-%m-%dT%H:%M:%S"):
        titleStr = datetime.strftime(startingDatetime, "%Y-%m-%dT%H:%M:%S")
    else:
        titleStr = datetime.strftime(startingDatetime, "%Y-%m-%dT%H:%M:%S") + ' to ' + datetime.strftime(endingDatetime,
                                                                                                         "%Y-%m-%dT%H:%M:%S")

    # Download saber data between the starting and ending dates (inclusive):
    here = pathlib.Path(__file__).parent.resolve()
    saber_download_file_str = '../srcPython/saber_download.py'
    saber_download_file = here.joinpath(saber_download_file_str)
    cmd = 'python '+str(saber_download_file)+' '+ str(startingDatetime)[:10].replace('-', '') + ' ' + str(endingDatetime)[
                                                                                            :10].replace('-',
                                                                                                         '') + ' ' + saveLoc
    
    os.system(cmd)

    # Read in the SABER data and generate the SABER data matrix and all other associated information:
    dateList = np.asarray(pd.date_range(startingTimeStr[:-9], endingTimeStr[:-9],
                                        periods=(endingDatetime - startingDatetime).days + 1).to_pydatetime())
    # np.asarray(pd.date_range('20' + startingTimeStr[1:3] + "-" + startingTimeStr[3:5] + "-" + startingTimeStr[5:7],
    #                          '20' + endingTimeStr[1:3] + "-" + endingTimeStr[3:5] + "-" + endingTimeStr[5:7], periods=(endingDatetime - startingDatetime).days+1).to_pydatetime())
    
    combined_NO_data = []
    print('Interpolating GITM data to the SABER locations...')
    gitm_time_starts = []
    gitm_time_ends = []
    goodFiles = []
    badFiles = []
    idx = 0
    for date in dateList:
        current_directory = saveLoc + str(date.year) + '/' + dateList[idx].strftime('%j')
        current_saberfiles = os.listdir(current_directory)
        # Order the files in 'current_saberfiles' chronologically:
        current_filename_tref = [int(element[10:17]) + int(element[18:23]) for element in current_saberfiles]
        
        ordered_current_saberfiles = np.array(current_saberfiles)[np.argsort(current_filename_tref)]
        # For each day, read in the SABER data:
        jdx = 0
        for saberfile in ordered_current_saberfiles:

            saber_file_output = saber_file_data(current_directory + '/' + saberfile, dateRef=date)

            if saber_file_output is None:
                badFiles.append(saberfile)
                pdb.set_trace()
                continue
            else:
                current_alts, current_lats, current_lons, current_NO_data, current_observation_times, lowestDate, highestDate = saber_file_output

            # Read in GITM data that is closest in time and space to the saber data; interpolate to the SABER location:
            gitm_time_start, gitm_ind_start = find_nearest(orderedDatetimeStamps, lowestDate)
            gitm_time_end, gitm_ind_end = find_nearest(orderedDatetimeStamps, highestDate)

            try:
                # TODO: Fix issue with NO GITM DATA AFTER 2011-08-07
                gitm_NO_data = gitm_NO_at_saber_loc_strict([current_alts, current_lats, current_lons, current_NO_data, current_observation_times],
                                                           dataLoc=gitm_data_dir, filenames=orderedValidFilenames ,fileInds=[gitm_ind_start, gitm_ind_end], startDate=gitm_time_start, endDate=gitm_time_end, override=override, fileFlag=[idx, jdx])
            except:
                pdb.set_trace()
            goodFiles.append(saberfile)
            
            if plot_type == 'peak':
                # SABER
                current_saber_peak_cooling, inds_saber_peak, saber_peak_times = extract_feature(current_NO_data,
                                                                              current_observation_times,
                                                                              feature=plot_type)
                # GITM
                current_gitm_peak_cooling, inds_gitm_peak, gitm_peak_times = extract_feature(gitm_NO_data,
                                                                            current_observation_times,
                                                                            feature=plot_type)
                combined_NO_data.append(
                    [current_alts, current_lats, current_lons, current_NO_data, gitm_NO_data, current_observation_times,
                     lowestDate, highestDate, current_saber_peak_cooling, inds_saber_peak, saber_peak_times,
                     current_gitm_peak_cooling, inds_gitm_peak, gitm_peak_times])
            elif plot_type == 'integral':
                # SABER
                current_saber_integrated_cooling, inds_saber_integrated, saber_integ_times = extract_feature([current_NO_data,
                                                                                           current_alts],
                                                                                          current_observation_times,
                                                                                          feature=plot_type)
                # GITM
                current_gitm_integrated_cooling, inds_gitm_integrated, gitm_integ_times = extract_feature(
                    [gitm_NO_data, current_alts], current_observation_times, feature=plot_type)
                combined_NO_data.append(
                    [current_alts, current_lats, current_lons, current_NO_data, gitm_NO_data, current_observation_times,
                     lowestDate, highestDate, current_saber_integrated_cooling, inds_saber_integrated,
                     saber_integ_times, current_gitm_integrated_cooling, inds_gitm_integrated, gitm_integ_times])
            elif plot_type == 'both':
                # SABER - Peak
                current_saber_peak_cooling, inds_saber_peak, saber_peak_times = extract_feature(current_NO_data,
                                                                              current_observation_times,
                                                                              feature='peak')
                # GITM - Peak
                current_gitm_peak_cooling, inds_gitm_peak, gitm_peak_times = extract_feature(gitm_NO_data, current_observation_times,
                                                                            feature='peak')
                # SABER - Integrated
                current_saber_integrated_cooling, inds_saber_integrated, saber_integ_times = extract_feature([current_NO_data,
                                                                                           current_alts],
                                                                                          current_observation_times,
                                                                                          feature='integral')
                # GITM - Integrated
                current_gitm_integrated_cooling, inds_gitm_integrated, gitm_integ_times = extract_feature(
                    [gitm_NO_data, current_alts], current_observation_times, feature='integral')
                #
                combined_NO_data.append(
                    [current_alts, current_lats, current_lons, current_NO_data, gitm_NO_data, current_observation_times,
                     lowestDate, highestDate, current_saber_peak_cooling, inds_saber_peak, saber_peak_times,
                     current_gitm_peak_cooling, inds_gitm_peak, gitm_peak_times, current_saber_integrated_cooling,
                     inds_saber_integrated, saber_integ_times, current_gitm_integrated_cooling, inds_gitm_integrated,
                     gitm_integ_times])
            else:
                combined_NO_data.append([current_alts, current_lats, current_lons, current_NO_data, gitm_NO_data, current_observation_times, lowestDate, highestDate])

            # gitm_time, gitm_ind = find_nearest(orderedDatetimeStamps, current_obs_time)
            # current_NO_prof_gitm, current_altitude_prof_gitm, current_lat_gitm, current_lon_gitm = (
            #     gitm_NO_at_saber_loc(dataLoc=gitm_data_dir, filenames=orderedValidFilenames, fileInd=gitm_ind,
            #                          lat=current_obs_lat, lon=current_obs_lon, alts=current_altitude_prof,
            #                          time=current_obs_time))

            # saberData.append([current_alts, current_lats, current_lons, current_NO_data, current_observation_times])
            # gitmData.append([current_alts_GITM, current_lats_GITM, current_lons_GITM, current_NO_data_GITM, current_observation_times_GITM])

            # Sanity check:
            # plt.figure()
            # plt.title('NO Cooling')
            # # SABER
            # plt.plot(current_NO_prof, current_altitude_prof, color='b', label='SABER')
            # # GITM
            # plt.plot(current_NO_prof_gitm, current_altitude_prof_gitm, color='r', label='GITM')
            # # Labels and Axes
            # plt.legend(loc='best')
            # plt.xlabel('NO Cooling (W/m$^3$)')
            # plt.ylabel('Altitude (km)')
            gitm_time_starts.append(gitm_time_start)
            gitm_time_ends.append(gitm_time_end)
            jdx += 1
        idx += 1
        print('--\nComplete for '+str(date)[:10]+'.\n--')
        time.sleep(0.5)
    
    if plot_type == 'both':
        pdb.set_trace()
        # Plot the Peak Cooling throughout the interval:
        all_peak_times_saber = np.concatenate([element[10] for element in combined_NO_data])
        all_peak_cooling_saber = np.concatenate([element[8] for element in combined_NO_data])
        all_peak_times_gitm = np.concatenate([element[13] for element in combined_NO_data])
        all_peak_cooling_gitm = np.concatenate([element[11] for element in combined_NO_data])
        plt.figure()
        sortInds_saber_peak = np.argsort(all_peak_times_saber)
        sortInds_gitm_peak = np.argsort(all_peak_times_gitm)
        plt.plot(all_peak_times_saber[sortInds_saber_peak], all_peak_cooling_saber[sortInds_saber_peak], label='SABER')
        plt.plot(all_peak_times_gitm[sortInds_gitm_peak], all_peak_cooling_gitm[sortInds_gitm_peak], label='GITM')
        plt.xlabel('Time')
        plt.ylabel('Peak NO Cooling (W/m$^3$)')
        plt.xlim([gitm_time_starts[0], gitm_time_ends[-1]])
        plt.xticks(rotation=45)
        plt.legend(loc='best')
        plt.tight_layout()
        plt.savefig(os.getcwd() + '/peak_cooling_' + str(gitm_time_starts[0]) + '_' + str(gitm_time_ends[-1]) + '.png', dpi=300)

        # Plot the height-integrated cooling throughout the interval:
        all_integ_times_saber = np.concatenate([element[16] for element in combined_NO_data])
        all_integ_cooling_saber = np.concatenate([element[14] for element in combined_NO_data])
        all_integ_times_gitm = np.concatenate([element[-1] for element in combined_NO_data])
        all_integ_cooling_gitm = np.concatenate([element[17] for element in combined_NO_data])
        plt.figure()
        sortInds_saber_integ = np.argsort(all_integ_times_saber)
        sortInds_gitm_integ = np.argsort(all_integ_times_gitm)
        plt.plot(all_integ_times_saber[sortInds_saber_integ], all_integ_cooling_saber[sortInds_saber_integ], label='SABER')
        plt.plot(all_integ_times_gitm[sortInds_gitm_integ], all_integ_cooling_gitm[sortInds_gitm_integ], label='GITM')
        plt.xlabel('Time')
        plt.ylabel('Height-Integrated NO Cooling (W/m$^2$)')
        plt.xlim([gitm_time_starts[0], gitm_time_ends[-1]])
        plt.xticks(rotation=45)
        plt.legend(loc='best')
        plt.tight_layout()
        plt.savefig(os.getcwd() + '/height_integrated_cooling_' + str(gitm_time_starts[0]) + '_' + str(gitm_time_ends[-1]) + '.png',
                    dpi=300)

        # TODO: Determine static x-limits for the two histograms below...
        # Histogram of residuals between GITM and SABER peak cooling:
        gitm_saber_residuals_peak_cooling = np.subtract(all_peak_cooling_gitm[sortInds_gitm_peak], all_peak_cooling_saber[sortInds_saber_peak])
        plt.figure()
        sns.histplot(gitm_saber_residuals_peak_cooling, bins=50, stat='density', label=r'$\mathrm{log}(\mu)$='+str(np.round(np.log10(np.nanmean(gitm_saber_residuals_peak_cooling)), 2))+
                                                                                       '\n$\mathrm{log}(\sigma)$='+str(np.round(np.log10(np.nanstd(gitm_saber_residuals_peak_cooling)), 2)) )
        sns.kdeplot(gitm_saber_residuals_peak_cooling, color='r')
        plt.xlabel('Residual (W/m$^3$)')
        plt.title('Distribution of Residuals of Peak NO Cooling')
        plt.legend(loc='best')
        plt.savefig(os.getcwd() + '/hist_peak_cooling_residuals.png', dpi=300)

        # Histogram of residuals between GITM and SABER integrated cooling:
        gitm_saber_residuals_int_cooling = np.subtract(all_integ_cooling_gitm[sortInds_gitm_integ],
                                                        all_integ_cooling_saber[sortInds_saber_integ])
        plt.figure()
        sns.histplot(gitm_saber_residuals_int_cooling, bins=50, stat='density',
                     label=r'$\mu$=' + str(np.round(np.nanmean(gitm_saber_residuals_int_cooling),4)) + '\n$\sigma$=' + str(np.round(np.nanstd(gitm_saber_residuals_int_cooling),4)))
        sns.kdeplot(gitm_saber_residuals_int_cooling, color='r')
        plt.xlabel('Residual (W/m$^2$)')
        plt.title('Distribution of Residuals of Height-Integrated Peak NO Cooling')
        plt.legend(loc='best')
        plt.savefig(os.getcwd() + '/hist_integrated_cooling_residuals.png', dpi=300)

        # Plot histograms for each and compute the Wasserstein Metric between them:
        from scipy.stats import wasserstein_distance
        from scipy.stats import gaussian_kde

        fig, axs = plt.subplots(nrows=2, ncols=1, sharex=False, figsize=(10, 8))
        # Top: Peak Cooling
        bins = 50
        binRange = (np.min([np.nanmin(all_peak_cooling_saber), np.nanmax(all_peak_cooling_gitm)]),
                    np.max([np.nanmin(all_peak_cooling_saber), np.nanmax(all_peak_cooling_gitm)]))
        kde_saber_peak = gaussian_kde(all_peak_cooling_saber)
        kde_gitm_peak = gaussian_kde(all_peak_cooling_gitm[~np.isnan(all_peak_cooling_gitm)])
        xx = np.linspace(binRange[0], binRange[-1], num=bins)
        axs[0].hist(all_peak_cooling_saber, bins=bins, range=binRange, label='SABER', edgecolor='black', density=True, alpha=0.8)
        axs[0].plot(xx, kde_saber_peak(xx), color='blue', linewidth=3)
        axs[0].hist(all_peak_cooling_gitm, bins=bins, range=binRange, label='GITM', alpha=0.5, edgecolor='black', density=True)
        axs[0].plot(xx, kde_gitm_peak(xx), color='darkorange', linewidth=3)
        axs[0].set_xlabel('Peak NO Cooling (W/m$^2$)')
        axs[0].set_ylabel('Probability Density')
        axs[0].set_title('Distributions of Peak Nitric Oxide Cooling')
        dist_peak = wasserstein_distance(all_peak_cooling_saber[~np.isnan(all_peak_cooling_saber)],
                                         all_peak_cooling_gitm[~np.isnan(all_peak_cooling_gitm)])
        axs[0].legend(loc='best', title=r'$W=${:0.3e}'.format(dist_peak))

        # Bottom: Integrated Cooling
        binRangeInteg = (np.min([np.nanmin(all_integ_cooling_saber), np.nanmax(all_integ_cooling_gitm)]),
                    np.max([np.nanmin(all_integ_cooling_saber), np.nanmax(all_integ_cooling_gitm)]))
        kde_saber_integ = gaussian_kde(all_integ_cooling_saber)
        kde_gitm_integ = gaussian_kde(all_integ_cooling_gitm[~np.isnan(all_integ_cooling_gitm)])
        xx = np.linspace(binRangeInteg[0], binRangeInteg[-1], num=bins)
        axs[1].hist(all_integ_cooling_saber, bins=bins, range=binRangeInteg, label='SABER', edgecolor='black', density=True,
                    alpha=0.8)
        axs[1].plot(xx, kde_saber_integ(xx), color='blue', linewidth=3)
        axs[1].hist(all_integ_cooling_gitm, bins=bins, range=binRangeInteg, label='GITM', alpha=0.5, edgecolor='black',
                    density=True)
        axs[1].plot(xx, kde_gitm_integ(xx), color='darkorange', linewidth=3)
        axs[1].set_xlabel('Height-integrated Nitric Oxide Cooling (W/m$^2$)')
        axs[1].set_ylabel('Probability Density')
        axs[1].set_title('Distributions of Integrated Nitric Oxide Cooling')
        dist_integrated = wasserstein_distance(all_integ_cooling_saber[~np.isnan(all_integ_cooling_saber)],
                                         all_integ_cooling_gitm[~np.isnan(all_integ_cooling_gitm)])
        axs[1].legend(loc='best', title=r'$W=${:0.3e}'.format(dist_integrated))

        # Minor adjustments and saving the figure:
        plt.tight_layout()
        plt.savefig(
            'NO_distributions_'+str(gitm_time_starts[0])+'_'+str(gitm_time_ends[0])+'.png', dpi=300)
        plt.close()

    # data_dist_figure = True
    # # pdb.set_trace()
    # if data_dist_figure:
    #     # Lat vs Local time plot for ALL data (AACGM coordinates)
    #     desired_shape = (sum([element[0].shape[0] for element in combined_NO_data]), sum([element[0].shape[-1] for element in combined_NO_data]))
    #     only_alts = [element[0] for element in combined_NO_data]
    #     only_lats = [element[1] for element in combined_NO_data]
    #     only_lons = [element[2] for element in combined_NO_data]
    #     only_times = [element[5] for element in combined_NO_data]
    #     alt_arr = np.concatenate(only_alts, axis=0)
    #     lat_arr = np.concatenate(only_lats, axis=0)
    #     lon_arr = np.concatenate(only_lons, axis=0)
    #     time_arr = np.concatenate(only_times, axis=0)
    #     mlats = np.zeros_like(alt_arr)
    #     mlons = np.zeros_like(alt_arr)
    #     mlts = np.zeros_like(alt_arr)
    #     for rowIndex in tqdm(range(alt_arr.shape[0])):
    #         for columnIndex in range(alt_arr.shape[1]):
    #             try:
    #                 m_coord = aacgmv2.wrapper.get_aacgm_coord(lat_arr[rowIndex, columnIndex],
    #                                                           lon_arr[rowIndex, columnIndex],
    #                                                           alt_arr[rowIndex, columnIndex],
    #                                                           time_arr[rowIndex, columnIndex])
    #             except:
    #                 try:
    #                     apex_out = apexpy.Apex(date=time_arr[rowIndex, columnIndex])
    #                     try:
    #                         mlat, mlt = apex_out.convert(lat=lat_arr[rowIndex, columnIndex], lon=lon_arr[rowIndex, columnIndex], source='geo', dest='mlt', datetime=time_arr[rowIndex, columnIndex])
    #                         mlon = apex_out.mlt2mlon(mlt=mlt, dtime=time_arr[rowIndex, columnIndex])
    #                         m_coord = [mlat, mlon, mlt]
    #                     except UserWarning:
    #                         m_coord = [np.nan, np.nan, np.nan]
    #                 except:
    #                     m_coord = [np.nan, np.nan, np.nan]
    #             mlats[rowIndex, columnIndex] = m_coord[0]
    #             mlons[rowIndex, columnIndex] = m_coord[1]
    #             mlts[rowIndex, columnIndex] = m_coord[2]
    #     mlat_1d = np.ravel(mlats)
    #     mlons_1d = np.ravel(mlons)
    #     mlts_1d = np.ravel(mlts)
    #     sorted_mlt_inds = np.argsort(mlts_1d)
    #     #
    #     plt.figure()
    #     plt.plot(mlts_1d[sorted_mlt_inds], mlat_1d[sorted_mlt_inds])
    #     plt.xlabel('MLT')
    #     plt.ylabel('MLAT')
    #     plt.title('SABER Data Coverage: '+str(gitm_time_starts[0]).replace(' ','T')+' to '+str(gitm_time_ends[-1]).replace(' ','T'))
    #     plt.show()
    #     # One orbit of SABER data:
    #     geographic_lats_1d = np.ravel(lat_arr)
    #     times_1d = np.ravel(time_arr)
    #     times_1d_scc = np.zeros_like(times_1d)
    #     for jdx in range(len(times_1d_scc)):
    #         if type(times_1d[jdx]) == float:
    #             times_1d_scc[jdx] = np.nan
    #         else:
    #             times_1d_scc[jdx] = (times_1d[jdx] - times_1d[0]).total_seconds()
    #     sorted_time_inds = np.argsort(times_1d_scc)
    #     sorted_geographic_lats_1d = geographic_lats_1d[sorted_time_inds]
    #     sorted_times_ssc_1d = times_1d_scc[sorted_time_inds]
    #     sorted_times_1d = times_1d[sorted_time_inds]
    #     one_orbit_time_inds = np.where((sorted_times_ssc_1d >= sorted_times_ssc_1d[0]) & (sorted_times_ssc_1d <= (sorted_times_ssc_1d[0] + 97*60)) )[0]
    #     sorted_mlts_1d = mlts_1d[sorted_time_inds]
    #     #
    #     plt.figure()
    #     plt.plot(sorted_mlts_1d[one_orbit_time_inds], sorted_geographic_lats_1d[one_orbit_time_inds])
    #     plt.xlabel('MLT (hours)')
    #     plt.ylabel('Geographic Latitude (degrees)')
    #     plt.title('SABER Data Coverage: ' + str(gitm_time_starts[0]).replace(' ', 'T') + ' to ' + str(
    #         gitm_time_starts[0] + timedelta(minutes=97)).replace(' ', 'T'))

    # pdb.set_trace()

    print('Interpolations complete. Generating comparison plots in each latitudinal boundary ('+str(gitm_time_starts[0])+' to '+str(gitm_time_ends[-1])+').\n')
    print('Out of '+str(len(goodFiles)+len(badFiles))+' total files, '+str(len(goodFiles))+' were good and '+str(len(badFiles))+' were bad and contained too many NaN values to use for analysis.')
    time.sleep(1)

    # Loop through the NO profiles and average them in 5 km sliding windows (with CIs):
    cleaned_NO_data = []
    for data in combined_NO_data:
        current_cleaned_SABER_NO_data, current_cleaned_SABER_NO_data_low_CI, current_cleaned_SABER_NO_data_high_CI = smooth_array(data[3], current_alts, resolution=5)
        current_cleaned_GITM_NO_data, current_cleaned_GITM_NO_data_low_CI, current_cleaned_GITM_NO_data_high_CI = smooth_array(data[4], current_alts, resolution=5)
        cleaned_NO_data.append([data[0], data[1], data[2], current_cleaned_SABER_NO_data,
                                current_cleaned_SABER_NO_data_low_CI, current_cleaned_SABER_NO_data_high_CI,
                                current_cleaned_GITM_NO_data, current_cleaned_GITM_NO_data_low_CI,
                                current_cleaned_GITM_NO_data_high_CI, data[5], data[6], data[7]])

    # Generate comparison plots in different latitudinal regions:
    # Polar Region |MLAT| >= 80:
    polInds = [np.where(element[1][:, 0] >= 80) for element in cleaned_NO_data]
    if len(set([len(element[0]) for element in polInds])) == 1:
        print('No polar data.\n')
    else:
        polData = []
        for i in range(len(polInds)):
            if len(polInds[i]) > 0:
                soundingInfo_pol = cleaned_NO_data[i]
                currentInds = polInds[i][0]
                polData.append([soundingInfo_pol[0][currentInds, :], soundingInfo_pol[1][currentInds, :],
                                soundingInfo_pol[2][currentInds, :],
                                soundingInfo_pol[3][currentInds, :], soundingInfo_pol[4][currentInds, :],
                                soundingInfo_pol[5][currentInds, :],
                                soundingInfo_pol[6][currentInds, :], soundingInfo_pol[7][currentInds, :],
                                soundingInfo_pol[8][currentInds, :],
                                soundingInfo_pol[9][currentInds, :], soundingInfo_pol[10], soundingInfo_pol[11]])
        complot(polData, label='Polar', timeStart=gitm_time_starts[0], timeEnd=gitm_time_ends[-1])
        if contour:
            gitm_saber_contour(polData, label='Auroral', timeStart=gitm_time_starts[0], timeEnd=gitm_time_ends[-1])
            gitm_saber_lat_contour(polData, alt=120, label='Polar', timeStart=gitm_time_starts[0],
                                   timeEnd=gitm_time_ends[-1])
            
    # Auroral Region 60 <= |MLAT| < 80:
    aurInds = [np.where((element[1][:, 0] >= 60) & (element[1][:, 0] < 80)) for element in cleaned_NO_data]
    if len(set([len(element[0]) for element in aurInds])) == 1:
        print('No auroral data.\n')
    else:
        aurData = []
        for i in range(len(aurInds)):
            if len(aurInds[i]) > 0:
                soundingInfo_aur = cleaned_NO_data[i]
                currentInds = aurInds[i][0]
                aurData.append([soundingInfo_aur[0][currentInds, :], soundingInfo_aur[1][currentInds, :],
                                soundingInfo_aur[2][currentInds, :],
                                soundingInfo_aur[3][currentInds, :], soundingInfo_aur[4][currentInds, :],
                                soundingInfo_aur[5][currentInds, :],
                                soundingInfo_aur[6][currentInds, :], soundingInfo_aur[7][currentInds, :],
                                soundingInfo_aur[8][currentInds, :],
                                soundingInfo_aur[9][currentInds, :], soundingInfo_aur[10], soundingInfo_aur[11]])
        complot(aurData, label='Auroral', timeStart=gitm_time_starts[0], timeEnd=gitm_time_ends[-1])
        if contour:
            gitm_saber_contour(aurData, label='Auroral', timeStart=gitm_time_starts[0], timeEnd=gitm_time_ends[-1])
            gitm_saber_lat_contour(aurData, alt=120, label='Auroral', timeStart=gitm_time_starts[0],
                                   timeEnd=gitm_time_ends[-1])

    # Mid-latitudes 30 <= |MLAT| < 60:
    midInds = [np.where((element[1][:, 0] >= 30) & (element[1][:, 0] < 60)) for element in cleaned_NO_data]
    if len(set([len(element[0]) for element in midInds])) == 1:
        print('No midlatitude data.\n')
    else:
        midData = []
        for i in range(len(midInds)):
            if len(midInds[i]) > 0:
                soundingInfo_mid = cleaned_NO_data[i]
                currentInds = midInds[i][0]
                midData.append( [ soundingInfo_mid[0][currentInds, :], soundingInfo_mid[1][currentInds, :], soundingInfo_mid[2][currentInds, :],
                                  soundingInfo_mid[3][currentInds, :], soundingInfo_mid[4][currentInds, :], soundingInfo_mid[5][currentInds, :],
                                  soundingInfo_mid[6][currentInds, :], soundingInfo_mid[7][currentInds, :], soundingInfo_mid[8][currentInds, :],
                                  soundingInfo_mid[9][currentInds, :], soundingInfo_mid[10], soundingInfo_mid[11] ] )
        complot(midData, label='Midlatitudes', timeStart=gitm_time_starts[0], timeEnd=gitm_time_ends[-1])
        if contour:
            gitm_saber_contour(midData, label='Midlatitudes', timeStart=gitm_time_starts[0], timeEnd=gitm_time_ends[-1])
            gitm_saber_lat_contour(midData, alt=120, label='Midlatitudes', timeStart=gitm_time_starts[0], timeEnd=gitm_time_ends[-1])

    # Low Latitudes/Equatorial Region |MLAT| < 30:
    eqInds = [np.where(element[1][:, 0] < 30) for element in cleaned_NO_data]
    if len(set([len(element[0]) for element in eqInds])) == 1:
        print('No equatorial data.\n')
    else:
        eqData = []
        for i in range(len(eqInds)):
            if len(eqInds[i]) > 0:
                soundingInfo_eq = cleaned_NO_data[i]
                currentInds = eqInds[i][0]
                try:
                    eqData.append(
                    [soundingInfo_eq[0][currentInds, :], soundingInfo_eq[1][currentInds, :], soundingInfo_eq[2][currentInds, :],
                     soundingInfo_eq[3][currentInds, :], soundingInfo_eq[4][currentInds, :], soundingInfo_eq[5][currentInds, :],
                     soundingInfo_eq[6][currentInds, :], soundingInfo_eq[7][currentInds, :], soundingInfo_eq[8][currentInds, :],
                     soundingInfo_eq[9][currentInds, :], soundingInfo_eq[10], soundingInfo_eq[11]])
                except:
                    try:
                        eqData.append(
                    [soundingInfo_eq[0][currentInds, :-1], soundingInfo_eq[1][currentInds, :-1], soundingInfo_eq[2][currentInds, :-1],
                     soundingInfo_eq[3][currentInds, :-1], soundingInfo_eq[4][currentInds, :-1], soundingInfo_eq[5][currentInds, :-1],
                     soundingInfo_eq[6][currentInds, :-1], soundingInfo_eq[7][currentInds, :-1], soundingInfo_eq[8][currentInds, :-1],
                     soundingInfo_eq[9][currentInds, :-1], soundingInfo_eq[10], soundingInfo_eq[11]])
                    except:
                        # Usually this occurs because of an IndexError in at indices 6 through 8. The current naive solution involves TRUNCATING the indices by one:
                        pdb.set_trace()
                        raise Exception
        complot(eqData, label='Equatorial', timeStart=gitm_time_starts[0], timeEnd=gitm_time_ends[-1])
        if contour:
            gitm_saber_contour(eqData, label='Equatorial', timeStart=gitm_time_starts[0], timeEnd=gitm_time_ends[-1])
            gitm_saber_lat_contour(eqData, alt=120, label='Equatorial', timeStart=gitm_time_starts[0],
                                   timeEnd=gitm_time_ends[-1])

    print('Comparisons complete!')
    return

def gitm_saber_lat_contour(info, alt, label, timeStart, timeEnd):
    """
    Helper function that takes data aggregated by gitm_saber_direct_compare and creates latitude vs. time contour plots
    at a selected altitude level by the user).
    """
    # Arrange the NO VER data (and time data) into (a) a giant 2D array:
    for i in range(len(info)):
        if i > 0:
            # Times...
            stacked_times = np.vstack((stacked_times, current_times))
            # Latitudes...
            stacked_latitudes = np.vstack((stacked_latitudes, current_latitudes))
            # SABER...
            stacked_saber_NO = np.vstack((stacked_saber_NO, current_saber_NO))
            # GITM
            stacked_gitm_NO = np.vstack((stacked_gitm_NO, current_gitm_NO))
        else:
            stacked_times = info[i][9]
            stacked_latitudes = info[i][8]
            stacked_saber_NO = info[i][3]
            stacked_gitm_NO = info[i][6]
        current_times = info[i][9]
        current_latitudes = info[i][1]
        current_saber_NO = info[i][3]
        current_gitm_NO = info[i][6]

    # With the NO VER array complete, extract the y-scale (altitudes):
    all_alts = []
    for i in range(len(info)):
        current_alts = np.nanmean(info[i][0], axis=0)
        all_alts.append(current_alts)
    alt_scale = np.nanmean(np.array(all_alts), axis=0)

    # Data sorting:
    stacked_seconds = np.zeros_like(stacked_times)
    for row in range(stacked_times.shape[0]):
        for col in range(stacked_times.shape[1]):
            if type(stacked_times[row, col]) != datetime:
                stacked_seconds[row, col] = np.nan
            else:
                stacked_seconds[row, col] = (stacked_times[row, col] - timeStart).total_seconds()
    stacked_seconds_double = np.hstack((stacked_seconds, stacked_seconds))
    stacked_times_double = np.hstack((stacked_times, stacked_times))
    
    # Extract data at the altitude of interest:
    altVal, altIdx = find_nearest(alt_scale, alt)
    gitm_NO_subset = stacked_gitm_NO[:, altIdx]
    saber_NO_subset = stacked_saber_NO[:, altIdx]
    latitudes_subset = stacked_latitudes[:, altIdx]
    times_subset = stacked_times_double[:, altIdx]
    t_dt = np.array([(element - times_subset[0]).total_seconds() for element in times_subset])
    # Make the Contour Plots:
    # plt.figure()
    # plt.tricontourf(t_dt, latitudes_subset, gitm_NO_subset)
    
    numel = 1000
    xvals = np.linspace(np.min(t_dt), np.max(t_dt), numel)
    yvals = np.linspace(np.min(latitudes_subset), np.max(latitudes_subset), numel)
    xx, yy = np.meshgrid(xvals, yvals)
    time_subset_upsampled = pd.date_range(str(times_subset[0]), str(times_subset[-1]), periods=numel).to_pydatetime()
    xx_dates, yy = np.meshgrid(time_subset_upsampled, yvals)
    gitm_grid = griddata((t_dt, latitudes_subset), gitm_NO_subset, (xx, yy), method='cubic')
    gitm_grid[gitm_grid < 0] = 0
    saber_grid = griddata((t_dt, latitudes_subset), saber_NO_subset, (xx, yy), method='cubic')
    saber_grid[saber_grid < 0] = 0
    perc_dev_grid = np.divide( 100*(np.abs(np.subtract(gitm_grid, saber_grid))), saber_grid)
    vmin_ = 0 # np.nanpercentile([np.nanmax(gitm_grid), np.nanmax(saber_grid)], 5)
    vmax_ = np.nanpercentile([np.nanmax(gitm_grid), np.nanmax(saber_grid)], 75)
    fig, axs = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(18,8))
    # Colorbar boundaries:
    numLevels = 25
    levels = np.linspace(vmin_, vmax_, numLevels)
    levels_r = np.linspace(np.nanpercentile(perc_dev_grid, 25), np.nanpercentile(perc_dev_grid, 75), numLevels)
    # Latitude cutoff if considering the midlatitudes:
    if label == 'Midlatitudes':
        lat_cutoff = 42.
    else:
        lat_cutoff = np.nanmin(yy)
    # SABER
    # im_saber = axs[0].contourf(xx_dates, yy, saber_grid, cmap='bwr', levels=levels, vmin=vmin_, vmax=vmax_)
    im_saber = axs[0].imshow(saber_grid, aspect='auto', cmap='bwr', interpolation='bicubic',
                             extent=[np.min(xx_dates), np.max(xx_dates), lat_cutoff, np.nanmax(yy)],
                             vmin=vmin_, vmax=vmax_, origin='lower')
    axs[0].set_title('SABER')
    axs[0].set_ylabel('Latitude')
    axs[0].set_xlabel('Time')
    cbar_saber = plt.colorbar(im_saber)
    cbar_saber.set_label('NO Cooling (W/m$^3$)')
    # GITM
    # im_gitm = axs[1].contourf(xx_dates, yy, gitm_grid, cmap='bwr', levels=levels, vmin=vmin_, vmax=vmax_)
    im_gitm = axs[1].imshow(gitm_grid, aspect='auto', cmap='bwr', interpolation='bicubic',
                             extent=[np.min(xx_dates), np.max(xx_dates), lat_cutoff, np.nanmax(yy)],
                            vmin=vmin_, vmax=vmax_, origin='lower')
    axs[1].set_title('GITM')
    axs[1].set_xlabel('Time')
    cbar_gitm = plt.colorbar(im_gitm)
    cbar_gitm.set_label('NO Cooling (W/m$^3$)')
    # Residual
    # im_residual = axs[2].contourf(xx_dates, yy, perc_dev_grid, levels=levels_r, cmap='Spectral', vmin=np.nanpercentile(perc_dev_grid, 25), vmax=np.nanpercentile(perc_dev_grid, 75))
    im_residual = axs[2].imshow(perc_dev_grid, aspect='auto', cmap='bwr', interpolation='bicubic',
                            extent=[np.min(xx_dates), np.max(xx_dates), lat_cutoff, np.nanmax(yy)],
                            vmin=-250., vmax=250., origin='lower')
                            # vmin=np.nanpercentile(perc_dev_grid, 5), vmax=np.nanpercentile(perc_dev_grid, 95),
                            # origin='lower')
    axs[2].set_title('Percent Difference')
    axs[2].set_xlabel('Time')
    cbar_residual = plt.colorbar(im_residual)
    cbar_residual.set_label('Percent Deviation from SABER')
    # Axes Labeling and Ticks:
    plt.setp(axs[0].xaxis.get_majorticklabels(), rotation=45)
    plt.setp(axs[1].xaxis.get_majorticklabels(), rotation=45)
    plt.setp(axs[2].xaxis.get_majorticklabels(), rotation=45)
    plt.suptitle('NO Cooling at '+str(np.round(altVal,2))+' km Altitude ('+str(timeStart)+' to '+str(timeEnd)+')')
    # Save the plot
    plt.savefig('saber_gitm_latitude_vs_time'+label+'_'+str(timeStart).replace(" ", "T").replace(":", "_")+'_to_'+str(timeEnd).replace(" ", "T").replace(":", "_")+'.png', dpi=300)

    return

def gitm_saber_contour(info, label, timeStart, timeEnd):
    """
    Helper function that takes data aggregated by gitm_saber_direct_compare and creates altitude vs. time contour plots
    (with identical scales for the NO Volume Emission Rate).
    """
    # Arrange the NO VER data (and time data) into (a) a giant 2D array:
    for i in range(len(info)):
        if i > 0:
            # Times...
            stacked_times = np.vstack((stacked_times, current_times))
            # SABER...
            stacked_saber_NO = np.vstack((stacked_saber_NO, current_saber_NO))
            # GITM
            stacked_gitm_NO = np.vstack((stacked_gitm_NO, current_gitm_NO))
        else:
            stacked_times = info[i][9]
            stacked_saber_NO = info[i][3]
            stacked_gitm_NO = info[i][6]
        current_times = info[i][9]
        current_saber_NO = info[i][3]
        current_gitm_NO = info[i][6]

    # VITAL: SORT the NO data BY TIME:
    stacked_seconds = np.zeros_like(stacked_times)
    for row in range(stacked_times.shape[0]):
        for col in range(stacked_times.shape[1]):
            if type(stacked_times[row, col]) != datetime:
                stacked_seconds[row, col] = np.nan
            else:
                stacked_seconds[row, col] = (stacked_times[row, col] - timeStart).total_seconds()
    stacked_seconds_double = np.hstack((stacked_seconds, stacked_seconds))
    stacked_times_double = np.hstack((stacked_times, stacked_times))
    sortInds = np.argsort(stacked_seconds_double, axis=0)
    sorted_times = np.take_along_axis(stacked_times_double, sortInds, axis=0)
    sorted_saber = np.take_along_axis(stacked_saber_NO, sortInds, axis=0)
    sorted_gitm = np.take_along_axis(stacked_gitm_NO, sortInds, axis=0)

    # TODO: FIX THE X VALUES!!!
    
    # With the NO VER array complete, extract the y-scale (altitudes):
    all_alts = []
    for i in range(len(info)):
        current_alts = np.nanmean(info[i][0], axis=0)
        all_alts.append(current_alts)
    alt_scale = np.nanmean(np.array(all_alts), axis=0)

    # Then extract the x-scale (time):
    t = np.linspace(pd.Timestamp(timeStart).value, pd.Timestamp(np.nanmax(sorted_times)).value, sorted_times.shape[0]) # pd.Timestamp(timeEnd).value, sorted_times.shape[0])
    t = pd.to_datetime(t)
    t_dt = np.array([element.to_pydatetime() for element in t])

    # Subset to between 100 and 250 km in altitude:
    alt_subset_inds = np.where((alt_scale >= 100) & (alt_scale <= 250))[0]
    alt_subset = alt_scale[alt_subset_inds]
    sorted_saber_subset = sorted_saber[:, alt_subset_inds]
    sorted_gitm_subset = sorted_gitm[:, alt_subset_inds]

    # Finally, generate a contour plot: (https://stackoverflow.com/questions/65899604/plotly-using-one-colorbar-for-multiple-contour-plots)
    vmin = np.log10(np.min([np.nanpercentile(sorted_saber_subset, 25), np.nanpercentile(sorted_gitm_subset, 25)]))
    vmax = np.log10(np.max([np.nanmax(sorted_saber_subset), np.nanmax(sorted_gitm_subset)]))
    fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(15,8))
    axs[0].imshow(np.log10(sorted_saber_subset.T), cmap='bwr', extent=[t_dt[0], t_dt[-1], 100, 250], vmin=vmin, vmax=vmax, aspect='auto', interpolation='bicubic')
    axs[0].tick_params(axis='x', labelrotation=45)
    axs[0].set_xlabel('Time')
    axs[0].set_ylabel('Altitude (km)')
    axs[0].set_title('SABER')
    axs[0].set_xlim([np.nanmin(sorted_times), timeEnd])
    out = axs[1].imshow(np.log10(sorted_gitm_subset.T), cmap='bwr', extent=[t_dt[0], t_dt[-1], 100, 250], vmin=vmin, vmax=vmax, aspect='auto', interpolation='bicubic')
    axs[1].tick_params(axis='x', labelrotation=45)
    axs[1].set_xlabel('Time')
    axs[1].set_title('GITM')
    axs[1].set_xlim([np.nanmin(sorted_times), timeEnd])
    cbar = fig.colorbar(out, ax=axs[1])
    cbar.ax.set_ylabel('Logarithm of NO Cooling (W/m$^3$)')
    fig.suptitle('Nitric Oxide Volume Emission Rate', fontsize=16)
    fig.tight_layout()
    fig.savefig('saber_gitm_contour_'+label+'_'+str(timeStart).replace(" ", "T").replace(":", "_")+'_to_'+str(timeEnd).replace(" ", "T").replace(":", "_")+'.png', dpi=300)
    plt.close()

    # Generate a contour plot of the DIFFERENCE between GITM and SABER:
    # gitm_saber_perc_diff = 100*(np.divide(np.abs(np.subtract(sorted_gitm_subset.T, sorted_saber_subset.T)), sorted_saber_subset.T))
    # Replace all infs with NaNs:
    # gitm_saber_perc_diff[gitm_saber_perc_diff ==np.inf] = np.nan
    # Distribution of Percent Difference
    # fig, ax_h = plt.subplots()
    # ax_h.hist(gitm_saber_perc_diff.ravel(), bins=75, range=(0, 400), density=True)
    # ax_h.set_xlabel('Percent Difference (%)')
    # ax_h.set_ylabel('Probability Density')
    # fig.tight_layout()
    # plt.show()
    
    # Contour plot of residuals
    gitm_residual = np.log10(np.subtract(sorted_gitm_subset.T, sorted_saber_subset.T))
    # gitm_perc_error = 100 * (np.divide(np.abs(np.subtract(sorted_gitm_subset.T, sorted_saber_subset.T)), sorted_saber_subset.T))
    limit = np.max([np.abs(np.nanmin(gitm_residual)), np.abs(np.nanmax(gitm_residual))])
    fig_d, ax = plt.subplots()
    out_d = ax.imshow(gitm_residual, cmap='bwr', extent=[t_dt[0], t_dt[-1], 100, 250], vmin=-limit, vmax=limit, aspect='auto', interpolation='bicubic')
    ax.tick_params(axis='x', labelrotation=45)
    ax.set_xlabel('Time')
    ax.set_ylabel('Altitude (km)')
    ax.set_title('GITM Deviation from SABER NO Cooling')
    cbar_d = plt.colorbar(out_d, ax=ax)
    cbar_d.ax.set_ylabel('Residual (GITM - SABER) (W/m$^3$)')
    fig.tight_layout()
    plt.savefig('saber_gitm_residual_'+label+'_'+str(timeStart).replace(" ", "T").replace(":", "_")+'_to_'+str(timeEnd).replace(" ", "T").replace(":", "_")+'.png', dpi=300)
    plt.close()
         
    # pio.renderers.default = "browser"
    # contours = dict(start=0, end=np.nanmax([np.nanmax(sorted_gitm_subset), np.nanmax(sorted_saber_subset)]))
    # figc = make_subplots(rows=1, cols=2, shared_xaxes=True, shared_yaxes=True, horizontal_spacing=0.05)
    # figc.add_trace(
    #     go.Contour(
    #         z=sorted_saber_subset.T,
    #         x=t_dt,  # horizontal axis
    #         y=alt_subset,  # vertical axis
    #         contours=contours
    #     ),
    #     row=1,
    #     col=1
    # )
    # figc.add_trace(
    #     go.Contour(
    #         z=sorted_gitm_subset.T,
    #         x=t_dt,  # horizontal axis
    #         y=alt_subset,  # vertical axis
    #         contours=contours,
    #         colorbar=dict(
    #             title=dict(
    #                 text='NO Cooling (W/m^3)',
    #                 side='right',
    #                 font=dict(
    #                     size=20,
    #                     family='Arial, sans-serif'
    #                 )
    #             )
    #         )
    #     ),
    #     row=1,
    #     col=2
    # )
    # figc.update_layout(
    #     title=dict(
    #         text="SABER and GITM NO Cooling ("+label+"):\n" + str(t_dt[0]) + " to " + str(t_dt[-1]),
    #         font=dict(
    #             size=25
    #         )
    #     ),
    #     title_x=0.5,
    #     yaxis_title="Altitude (km)",
    #     width=1800,
    #     height=1000
    # )
    # figc.update_yaxes(tickfont=dict(size=25))
    # show(figc)
    # # Save the contour plot:
    # pio.write_image(figc,
    #                 'saber_gitm_contour_'+label+
    #                 '_'+str(timeStart).replace(" ", "T").replace(":", "_")+
    #                 '_to_'+str(timeEnd).replace(" ", "T").replace(":", "_"), format='png')
    
    return

def show(fig):
    """
    Helper function to display Plotly images in a standalone window.
    Source: https://stackoverflow.com/questions/53570384/plotly-how-to-make-a-standalone-plot-in-a-window
    Takes as its sole argument a fig object generated by plotly.graph_objs
    """
    import io
    import plotly.io as pio
    from PIL import Image
    buf = io.BytesIO()
    pio.write_image(fig, buf)
    img = Image.open(buf)
    img.show()

def complot(input, saveLoc=os.getcwd(), label=None, timeStart=None, timeEnd=None):
    """
    Take processed NO data for SABER and GITM, and make a comparison plot (for a given latitudinal region). Takes as
    an optional argument the location with which to save a figure.
    input format:
    0 - Alts
    1 - Lats
    2 - Lons
    3 - Cleaned SABER NO data
    4 - Cleaned SABER NO data (low CI)
    5 - Cleaned SABER NO data (high CI)
    6 - Cleaned GITM NO data
    7 - Cleaned GITM NO data (low CI)
    8 - Cleaned GITM NO data (high CI)
    9 - Observation times
    10 - Start Date of Observations
    11 - End Date of Observations
    """
    averaged_alt_rulers = []
    averaged_profs_SABER = []
    averaged_profs_GITM = []
    numSoundings = 0
    for i in range(len(input)):
        averaged_alt_rulers.append(np.average(input[i][0], axis=0))
        averaged_profs_SABER.append(np.average(input[i][3], axis=0))
        averaged_profs_GITM.append(np.average(input[i][6], axis=0))
        numSoundings += input[i][0].shape[0]
    averaged_alt_rulers_arr = np.array(averaged_alt_rulers)
    averaged_profs_SABER_arr = np.array(averaged_profs_SABER)
    averaged_profs_GITM_arr = np.array(averaged_profs_GITM)

    # Averaging, with CIs:
    averaged_alt_ruler = np.average(averaged_alt_rulers_arr, axis=0)
    #
    average_prof_SABER = np.nanmean(averaged_profs_SABER_arr, axis=0)
    average_std_SABER = np.nanstd(averaged_profs_SABER_arr, axis=0)
    average_prof_SABER_upper_CI = np.add(average_prof_SABER, 1 * average_std_SABER)
    average_prof_SABER_lower_CI = np.subtract(average_prof_SABER, 1 * average_std_SABER)
    average_prof_SABER_lower_CI[average_prof_SABER_lower_CI < 0] = 0
    #
    average_prof_GITM = np.nanmean(averaged_profs_GITM_arr, axis=0)
    average_std_GITM = np.nanstd(averaged_profs_GITM_arr, axis=0)
    average_prof_GITM_upper_CI = np.add(average_prof_GITM, 1 * average_std_GITM)
    average_prof_GITM_lower_CI = np.subtract(average_prof_GITM, 1 * average_std_GITM)
    average_prof_GITM_lower_CI[average_prof_GITM_lower_CI < 0] = 0

    # Subset everything to between 100 and 250 km:
    relevant_alt_inds = np.where((averaged_alt_ruler >= 100) & (averaged_alt_ruler <= 250))[0]
    alt_subset = averaged_alt_ruler[relevant_alt_inds]
    average_prof_SABER_lower_CI_subset = average_prof_SABER_lower_CI[relevant_alt_inds]
    average_prof_SABER_upper_CI_subset = average_prof_SABER_upper_CI[relevant_alt_inds]
    average_prof_SABER_subset = average_prof_SABER[relevant_alt_inds]
    average_prof_GITM_lower_CI_subset = average_prof_GITM_lower_CI[relevant_alt_inds]
    average_prof_GITM_upper_CI_subset = average_prof_GITM_upper_CI[relevant_alt_inds]
    average_prof_GITM_subset = average_prof_GITM[relevant_alt_inds]

    # Plotting!
    plt.figure()
    #
    plt.fill_betweenx(alt_subset, average_prof_SABER_lower_CI_subset, average_prof_SABER_upper_CI_subset,
                     where=average_prof_SABER_lower_CI_subset<average_prof_SABER_upper_CI_subset, color='b', alpha=0.3)
    plt.plot(average_prof_SABER_subset, alt_subset, color='b', label='SABER')
    #
    plt.fill_betweenx(alt_subset, average_prof_GITM_lower_CI_subset, average_prof_GITM_upper_CI_subset,
                     where=average_prof_GITM_lower_CI_subset < average_prof_GITM_upper_CI_subset, color='r', alpha=0.3)
    plt.plot(average_prof_GITM_subset, alt_subset, color='r', label='GITM')
    plt.plot([], [], ' ', label=r'$\mathrm{log}(L^2)$: '+str(np.round(np.log10(euclidean(average_prof_GITM_subset, average_prof_SABER_subset)), 2)))
    #
    plt.xlabel('NO Cooling (W/m$^3$)')
    plt.ylabel('Altitude (km)')
    plt.legend(loc='best')
    timeStr = str(timeStart)+' to '+str(timeEnd)
    plt.title('SABER and GITM NO Cooling: '+label+' ('+str(numSoundings)+' soundings)\n'+timeStr)
    plt.tight_layout()
    # Save the plot...
    plt.savefig(saveLoc+'/NO_comparison_'+label+'_'+str(numSoundings)+'_soundings_'+timeStr.replace(" ","-").replace(":","_")+'.png', dpi=300)
    plt.close()

    return

def euclidean(d_1, d_2):
    """
    Compute the Euclidean distance between two arraylikes of equal length.
    :param d_1: arraylike
        Some univariate data.
    :param d_2: arraylike
        Some other univariate data (same length as d_1).
    :return dist: float
        The Euclidean distance between d_1 and d_2.
    """
    dist = np.sqrt( np.sum( np.power(np.subtract(d_1, d_2), 2)) )
    return dist

def smooth_array(array, ruler, resolution):
    """
    Given some 2D array, a 'ruler', and a 'resolution', smooth the array in the 'vertical direction' according to the
    desired resolution, which will be converted to a stepsize (in terms of elements) using the provided ruler.
    :param array: np.ndarray
        A 2D (nxm) array of data, where the m dimension is the 'vertical' direction.
    :param ruler: np.ndarray or list
        A 2D array (nxm), consonant with the m dimension of 'array' that gives the knots for which the vertical
        dimension of the 2D array is gridded. Used to compute a step size, when 'resolution' is combined with it. Values
        in this input have their own units, usually km.
    :param resolution: int
        A desired stepsize (or resolution) in the SAME UNITS as the elements of 'ruler'. Used to determine a
        unit-agnostic stepsize.
    """
    # Get unit-dependent steps in the vertical domain:
    vertSteps = []
    for i in range(ruler.shape[0]):
        vertSteps.append(np.diff(ruler[i, :]))
    vertStepsArray = np.array(vertSteps)
    # Get the average profile for unit-dependent steps in the vertical domain:
    ave_vert_step_profile = np.abs(np.nanmean(vertStepsArray, axis=1))
    # Extract the LARGEST unit-dependent step:
    biggest_step = np.max(ave_vert_step_profile) # km
    # Convert that unit-dependent step to an element-dependent stepsize:
    windowSize = int(np.round(resolution / biggest_step))

    # Now that the window size is obtained, smooth the data in the vertical domain (use a Savitzky-Golay Filter).
    # smoothedArray_sg = savgol_filter(array, window_length=windowSize, polyorder=3, axis=-1)

    smoothedArray, smoothedArray_low_CI, smoothedArray_high_CI = ma_ci(array, window_length=windowSize, axis=-1, cutoff=0, valreplace=0)

    # plt.figure();
    # plt.plot(array[0, :]), plt.plot(smoothedArray_sg[0, :]);
    # xvals = np.linspace(0, len(smoothedArray[0, :])-1, len(smoothedArray[0, :]))
    # plt.fill_between(xvals, smoothedArray_low_CI[0, :], smoothedArray_high_CI[0, :], where=smoothedArray_low_CI[0, :]<smoothedArray_high_CI[0, :], interpolate=True, color='tab:green', alpha=0.5)
    # plt.plot(smoothedArray[0, :], color='tab:green')

    return smoothedArray, smoothedArray_low_CI, smoothedArray_high_CI

def ma_ci(array, window_length, axis=-1, cutoff=0, valreplace=0):
    """
    Smooth an 2D array along a single axis, using a moving average. If supplied, replace bad values with a good value.
    """
    smoothedArray = np.zeros_like(array)
    smoothedArray_low_CI = np.zeros_like(array)
    smoothedArray_high_CI = np.zeros_like(array)
    for i in range(array.shape[0]):
        # Loop through each sounding and compute the sliding window along with 95%CI values:
        sounding_series = pd.Series(array[i, :])
        smoothedArray[i, :] = np.array(sounding_series.rolling(window=window_length, center=True, min_periods=1).mean())
        var = 2 * np.array(sounding_series.rolling(window=window_length, center=True, min_periods=1).std())
        smoothedArray_low_CI[i, :] = smoothedArray[i, :] - var
        smoothedArray_high_CI[i, :] = smoothedArray[i, :] + var
    # Replace bad values:
    smoothedArray[smoothedArray<cutoff] = valreplace
    smoothedArray_low_CI[smoothedArray_low_CI<cutoff] = valreplace
    smoothedArray_high_CI[smoothedArray_high_CI<cutoff] = valreplace
    return smoothedArray, smoothedArray_low_CI, smoothedArray_high_CI

def gitm_NO_at_saber_loc_strict(saberData, dataLoc, filenames, fileInds, startDate, endDate, override=False, **kwargs):
    """
    Does the same thing as 'gitm_NO_at_saber_loc' but recieves the output of 'saber_file_data' and finds GITM outputs
    at every single point of every single SABER sounding.
    """
    # String for filename to look for:
    timeStr = str(startDate).replace(" ", "T").replace(":", "_") + "_to_" + str(endDate).replace(" ", "T").replace(":",
                                                                                                                   "_")
    if kwargs:
        # This step prevents there from being duplicate .pkl files for time intervals so close that the starting and ending dates are actually identical.
        fileIndices = kwargs['fileFlag']
        GITM_NO_FILENAME = dataLoc+'gitm_NO_data_'+timeStr+'_'+str(fileIndices[0])+'_'+str(fileIndices[1])+'.pkl'
    else:
        GITM_NO_FILENAME = dataLoc+'gitm_NO_data_'+timeStr+".pkl"
    if os.path.isfile(GITM_NO_FILENAME) == False or override == True:
        # Unpack inputs
        current_alts, current_lats, current_lons, current_NO_data, current_observation_times = saberData
        stacked_observation_times = np.hstack((current_observation_times, current_observation_times))

        # Get list of ALL GITM files that will need to be considered:
        fileIndices = np.arange(fileInds[0]-1, fileInds[-1]+2, 1)
        gitm_data = []
        for ind in fileIndices:
            try:
                out = gitm_file_mean_profile(dataLoc + '/' + filenames[
                    ind])  # NO_prof[:,:, alt_inds], mean_altitude_prof[alt_inds], obs_time, np.asarray(info['lat']), np.asarray(info['lon'])
            except IndexError:
                # Just double-count the last file if nothing is beyond it...
                if ind == fileIndices[-1]:
                    out = gitm_file_mean_profile(dataLoc + '/' + filenames[ind - 1])
            gitm_data.append(out)

        # Extracting GITM variables:
        gitm_alts = [element[1] for element in gitm_data]
        gitm_lats = [element[3] for element in gitm_data]
        gitm_lons = [element[4] for element in gitm_data]
        gitm_NO_data_raw = np.asarray([element[0] for element in gitm_data])
        gitm_times = np.array([element[2] for element in gitm_data])
        gitm_times_s = [(element - gitm_times[0]).total_seconds() for element in gitm_times]

        # Loop through the SABER observation times and grab the closest GITM data for each observation; interpolate when
        # necessary:
        gitm_NO_data = np.zeros_like(current_NO_data)
        print('Extracting GITM data for SABER NETCDF4 file covering the time period: '+str(startDate) + ' to '+str(endDate)+'...')

        too_high = []
        # TODO: Include modifications for when data is outside lat/lon bounds...

        # Clean the observation times of NaNs:
        # if endDate >= datetime(2011, 8, 3, 23, 40):
        #     raise Exception
        # try:
        #     stacked_observation_times_c = clean_NaNs(stacked_observation_times, axis=-1)
        # except:
        #     raise Exception

        stacked_observation_times_c = clean_NaNs(stacked_observation_times, axis=-1)

        for i in tqdm(range(gitm_NO_data.shape[0])):
            # Rows
            for j in range(gitm_NO_data.shape[1]):
                # Columns
                currentTime = stacked_observation_times_c[i, j]
                if type(currentTime) == float:
                    gitm_NO_data[i, j] = np.nan
                    validData = False
                else:
                    # try:
                    # Isolate the GITM data that is relevant:
                    if currentTime < startDate:
                        currentTime = startDate
                    if currentTime > endDate:
                        currentTime = endDate
                    # if currentTime >= datetime(2011, 8, 7, 4):
                        # pdb.set_trace()
                    _, current_file_ind = find_nearest(gitm_times, currentTime)
                    current_gitm_subset = gitm_NO_data_raw[np.array([current_file_ind-1, current_file_ind, current_file_ind+1])]
                    gitm_times_subset = np.asarray(gitm_times)[np.array([current_file_ind-1, current_file_ind, current_file_ind+1])]
                    gitm_times_s_subset = np.asarray(gitm_times_s)[np.array([current_file_ind-1, current_file_ind, current_file_ind+1])]
                    validData = True

                    # Interpolate the GITM data:
                    interp = RegularGridInterpolator((np.array([0, 0.5, 1]), gitm_lons[0], gitm_lats[0], gitm_alts[0]), current_gitm_subset, bounds_error=False)

                    # Get the location of the datapoint to interpolate:
                    time_s = (currentTime - gitm_times[0]).total_seconds() / gitm_times_s_subset[-1] # Normalize the time value to be between 0 and 1.
                    current_point = np.array([time_s, current_lons[i,j], current_lats[i, j], current_alts[i, j]])
                    if (current_point[-1] >= 244) == True:
                        too_high.append([i,j,current_point])
                        current_point[-1] = np.nanmax(gitm_alts[0])
                        gitm_NO_data[i,j] = interp(current_point)
                    else:
                        out = interp(current_point)
                            # try:
                                # out = interp(current_point)
                            # except:
                                # pdb.set_trace()
                        # Save the interpolated datapoint:
                        gitm_NO_data[i, j] = out
                    # except:
                        # The presence of masked data indicates a bad value:
                        # pdb.set_trace()
                        # gitm_NO_data[i, j] = np.nan

        # Sanity checks (viewing profiles side-by-side):
        # for k in range(current_NO_data.shape[0]-90):
        #     plt.figure(); plt.plot(current_NO_data[k, :], current_alts[k, :], label='SABER');  plt.plot(gitm_NO_data[k, :], current_alts[k, :], label='GITM'); plt.legend(loc='best')

        # Final cleaning; any remaining NaNs are simply interpolated over (across profiles rather than vertically along
        # an individual profile):
        clean_data = True
        if clean_data:
            gitm_NO_data_cleaned = interp2d_simple(gitm_NO_data.T, k=3, axis=-1)
            gitm_NO_data = gitm_NO_data_cleaned.T

        # Sanity check - Looking at SABER data spatial coverage overtop GITM outputs.
        timeSelect = current_observation_times[45, 0]
        tVal, tDex = find_nearest(gitm_times, timeSelect)
        gitm_data_subset = gitm_NO_data_raw[tDex, :, :, :]
        aVal, aDex = find_nearest(gitm_alts[0], 125)
        gitm_data_at_125_km = gitm_data_subset[:, :, aDex] # SHAPE: Lon, Lat
        plt.figure()
        plt.imshow(gitm_data_at_125_km.T, extent=[np.min(gitm_lons[0]), np.max(gitm_lons[0]), np.min(gitm_lats[0]), np.max(gitm_lats[0])])
        for i in range(current_lons.shape[0]):
            plt.plot(current_lons[i, :], current_lats[i, :], color='r')
        plt.xlabel('Geographic Longitude')
        plt.ylabel('Geographic Latitude')
        plt.xlim([np.min(gitm_lons[0]), np.max(gitm_lons[0])])
        plt.ylim([np.min(gitm_lats[0]), np.max(gitm_lats[0])])
        plt.title('SABER Observations over GITM Data at 125 km: '+str(stacked_observation_times[0, 0])+' to '+str(stacked_observation_times[-1, 0]))
        plt.close()

        # CACHE the interpolated NO data (to the GITM directory in 'dataLoc'), so that if it already exists, it is just read-in:
        NO_dict = {
            'gitm_NO_data_raw': gitm_NO_data_raw,
            'gitm_NO_data': gitm_NO_data
        }
        
        pickle_save(NO_dict, GITM_NO_FILENAME)
    else:
        print('Loaded in GITM data for SABER NETCDF4 file covering the time period: ' + str(startDate) + ' to ' + str(
            endDate) + '...')
        NO_dict = pickle_load(GITM_NO_FILENAME)
        gitm_NO_data_raw = NO_dict['gitm_NO_data_raw']
        gitm_NO_data = NO_dict['gitm_NO_data']
    
    return gitm_NO_data

def clean_NaNs(data, axis=-1):
    new_data = np.zeros_like(data)
    for idx in range(data.shape[0]):
        row = data[idx,:]
        reference = firstNonNan(row)
        arr = []
        for element in row:
            if type(element) == float:
                arr.append(np.nan)
            else:
                arr.append((element - reference).total_seconds())
        #arr = np.array([(element - reference).total_seconds() for element in row])
        arr = np.asarray(arr)
        arr_clean = clean_arr(arr)
        new_data[idx, :] = np.array([reference + timedelta(seconds=element) for element in arr_clean])
        
    return new_data

def firstNonNan(listarr):
  for item in listarr:
    if type(item) != float:
      return item

def clean_arr(arr):
    x_vals = np.linspace(0, len(arr)-1, len(arr))
    try:
        interpolator = interp1d(x_vals[~np.isnan(arr)], arr[~np.isnan(arr)], kind='linear', bounds_error=False, fill_value='extrapolate')
    except:
        raise Exception
    clean_arr = interpolator(x_vals)
    return clean_arr

def gitm_NO_at_saber_loc(dataLoc, fileInd, filenames, lat, lon, alts, time):
    """
    Given the location of GITM data, a file to read in, and the desired lat and lon, extract a Nitric Oxide VER profile.
    Use interpolation in 3D to obtain the desired profile.
    :param dataLoc: str
        A string with the location of the GITM data.
    :param filenames: list of strings
        A list of strings of the filenames to get data from.
    :param fileInd: int
        The file number identifying the file to read in; provided by gitm_saber_direct_compare.
    :param lat: int
        The latitude at which to extract the data.
    :param lon: int
        The longitude at which to extract the data.
    :param alts: list or array
        The altitude values at which to extract the data.
    :param time: datetime
        The time at which to extract the data.
    """
    # Step 1: Read in the GITM data
    fileIndices = [fileInd-1, fileInd, fileInd+1]
    gitm_data = []
    for ind in fileIndices:
        try:
            out = gitm_file_mean_profile(dataLoc+'/'+filenames[ind]) # NO_prof[:,:, alt_inds], mean_altitude_prof[alt_inds], obs_time, np.asarray(info['lat']), np.asarray(info['lon'])
        except IndexError:
            # Just double-count the last file if nothing is beyond it...
            if ind == fileIndices[-1]:
                out = gitm_file_mean_profile(dataLoc + '/' + filenames[ind-1])
        gitm_data.append(out)
    gitm_alts = [element[1] for element in gitm_data]
    gitm_lats = [element[3] for element in gitm_data]
    gitm_lons = [element[4] for element in gitm_data]

    # Step 2: Interpolate the GITM data to the location and time desired
    NO_data = np.asarray([element[0] for element in gitm_data])
    gitm_times = [element[2] for element in gitm_data]
    gitm_times_s = [(element - gitm_times[0]).total_seconds() for element in gitm_times]
    time_s = gitm_times_s[-1] / (time - gitm_times[0]).total_seconds() # Normalize the time value to be between 0 and 3.

    # First, interpolate in time in order to get an idea of how densely to interpolate in space.
    # closest_time, closest_time_ind = find_nearest(gitm_times, time)
    # time_diff = np.abs(time - closest_time) # Round to nearest seconds
    # time_diff_pd = pd.Timedelta(time_diff)
    # res = time_diff_pd.round('s')
    # interp_times = np.asarray(pd.date_range(gitm_times[0], gitm_times[-1], freq=res).to_pydatetime())
    # final_time, final_time_ind = find_nearest(interp_times, time)
    # # Second, do the interpolation in space-time (in 4D):

    # Obtain interpolated data at each altitude:
    interp = RegularGridInterpolator((np.array([1, 2, 3]), gitm_lons[0], gitm_lats[0], gitm_alts[0]), NO_data, bounds_error=False)
    gitm_NO_interp = np.zeros_like(alts)
    if lon < 0:
        lon = 360 + lat
    for i in range(len(alts)):
        current_point = np.array([time_s, lon, lat, alts[i]])
        out = interp(current_point)
        gitm_NO_interp[i] = out

    # from scipy.interpolate import interpn
    # grid = (np.array([1, 2, 3]), gitm_lons[0], gitm_lats[0], gitm_alts[0])
    # point = np.array([time_s, lon, lat, alts[0]])
    # result = interpn(grid, NO_data, point, method='linear')

    return gitm_NO_interp, alts, lat, lon

def writeCSV(profile, profile_lb, profile_ub, altitudes, source, label=None, location=None):
    """
    Helper function to save processed GITM altitude profile data to a CSV file.
    :param profile: numpy.ndarray
        Mean Nitric Oxide altitude profile for SABER.
    :param profile_lb: numpy.ndarray
        5% CI lower boundary Nitric Oxide altitude profile for SABER.
    :param profile_ub: numpy.ndarray
        95% CI upper boundary Nitric Oxide altitude profile for SABER.
    :param altitudes: numpy.ndarray
        Altitude values in km, for the Nitric Oxide profile.
    :param source: str
        Either 'GITM' or 'SABER'.
    :param label: str
        Optional label for the latitude region the Nitric Oxide profile corresponds to.
    :param location: str
        The place where the CSV file should be saved. If None, defaults is the working directory.
    :return:
    """
    if location is None:
        location = os.getcwd() + '/'
    if source != 'SABER' and source != 'GITM':
        print(source)
        raise ValueError("Invalid value for argument 'source' - must be either 'SABER' or 'GITM'.")

    with open(location + source + '_NO_Profiles_'+label+'.csv', 'w') as f:
        writer = csv.writer(f, delimiter=',', lineterminator='\n')
        writer.writerow(['Mean_NO_Emission_(nW/m^3)', '5%CI_NO_Emission_(nW/m^3)', '95%CI_NO_Emission_(nW/m^3)', 'Altitude_(km)'])
        for i in range(len(profile)):
            writer.writerow([str(profile[i]), str(profile_lb[i]), str(profile_ub[i]), str(altitudes[i])])
    return

# def saberCSV(profile, profile_lb, profile_ub, altitudes, label=None, location=None):
#     """
#     Helper function to save processed SABER altiude profile data to a CSV file.
#     :param profile: numpy.ndarray
#         Mean Nitric Oxide altitude profile for SABER.
#     :param profile_lb: numpy.ndarray
#         5% CI lower boundary Nitric Oxide altitude profile for SABER.
#     :param profile_ub: numpy.ndarray
#         95% CI upper boundary Nitric Oxide altitude profile for SABER.
#     :param altitudes: numpy.ndarray
#         Altitude values in km, for the Nitric Oxide profile.
#     :param label: str
#         Optional label for the latitude region the Nitric Oxide profile corresponds to.
#     :param location: str
#         The place where the CSV file should be saved. If None, defaults is the working directory.
#     :return:
#     """
#     if location is None:
#         location = os.getcwd() + '/'
#     with open(location + 'SABER_NO_Profiles_'+label+'.csv', 'w') as saber_f:
#         writer = csv.writer(saber_f, delimiter=',', lineterminator='\n')
#         writer.writerow(['Mean_NO_Emission_(nW/m^3)', '5%CI_NO_Emission_(nW/m^3)', '95%CI_NO_Emission_(nW/m^3)', 'Altitude_(km)'])
#         for i in range(len(profile)):
#             writer.writerow([str(profile[i]), str(profile_lb[i]), str(profile_ub[i]), str(altitudes[i])])
#     return

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

# Compute 95%CI bounds:
def conf_int(profile_data, confidence=0.95):
    lb = []
    ub = []
    for i in range(profile_data.shape[1]):
        a = 1.0 * np.array(profile_data[:, i])
        n = len(a)
        m, se = np.nanmean(a), scipy.stats.sem(a, nan_policy='omit')
        h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
        lb.append(m-h)
        ub.append(m+h)
    return np.asarray(lb), np.asarray(ub)

def pickle_save(data, fname):
    with open(fname, 'wb') as handle:
        pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)

def pickle_load(fname):
    with open(fname, 'rb') as handle:
        b = pickle.load(handle)
    return b

# ----------------------------------------------------------------------
# Example command-line usage:
#
# python saber_compare.py /home/daabrand/Projects/TECHNO/Analysis2/20110805.Storm/v10x20/out/ /home/daabrand/Projects/TECHNO/Analysis2/20110805.Storm/v10x20/saber/ both False
#
args = parse_args()

gitm_data_loc = args.gitm_loc[0]
saber_data_loc = args.saber_loc[0]
plot_type = args.plot_type[0]
override = args.override[0]

gitm_saber_direct_compare(gitm_data_loc, saveLoc=saber_data_loc, contour=True, plot_type = plot_type, override = override)

# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# Execution (debugging)
# ----------------------------------------------------------------------
# if __name__ == '__main__':
    # # fname = '/home/daabrand/Projects/TECHNO/Analysis1/20110805.Storm/SABER/2011/217/SABER_L2A_2011217_52308_02.07.nc'
    # # NO_vals, altitude_vals, time_val, lat, lon = saber_file_mean_profile(fname)
    # #
    # # plt.figure()
    # # plt.plot(NO_vals, altitude_vals)
    # # plt.xlabel('NO Cooling (ergs/cm$^3$/s)')
    # # plt.ylabel('Altitude (km)')
    #
    # gitm_data_location = '/home/daabrand/Projects/TECHNO/Analysis2/20110805.Storm/v10x20/out/' # '/home/daabrand/Projects/TECHNO/Analysis2/20110805.Storm/v10x20/UA/data/'
    # saber_download_location = '/home/daabrand/Projects/TECHNO/Analysis2/20110805.Storm/v10x20/saber/'
    #
    # # gitm_saber_profile_comparison(gitm_data_location, saveLoc=saber_download_location)
    #
    # gitm_saber_direct_compare(gitm_data_location, saveLoc=saber_download_location, contour=True, plot_type = 'both', override = False)



