#  Takes in .txt file outputs from thermo_plot.py and compares the results of TWO runs.

# Top-level imports
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Qt5Agg')
import argparse
from datetime import datetime, timedelta
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.stats import linregress
import pdb

# ----------------------------------------------------------------------------------------------------------------------
# Global Plotting Settings:
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'large',
          'figure.figsize': (16, 8),
         'axes.labelsize': 'large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description='Compare the outputs of two GITM runs from thermo_plot.py.')
    parser.add_argument('var', nargs=1, \
                        help='The variable being plotted.')
    parser.add_argument('file1', nargs=1, \
                        help='Output .txt file for first GITM run.')
    parser.add_argument('file2', nargs=1, \
                        help='Output .txt file for second GITM run.')
    parser.add_argument('tag1', nargs=1, help='Label for first GITM run.')
    parser.add_argument('tag2', nargs=1, help='Label for second GITM run.')
    parser.add_argument('file3', nargs=1, help='Output .txt file for third GITM run.', default=None)
    parser.add_argument('tag3', nargs=1, help='Label for third GITM run.', default=None)
    parser.add_argument('loc', nargs='?', metavar='location',\
                        help='Location with which to save the combined figure.', default=os.getcwd())

    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------------------------------------------------------
def percErr(experimental, theoretical):
    return np.divide(np.subtract(experimental, theoretical), np.abs(theoretical)) * 100

def read_thermo(filename):
    """
    Read the output of thermo_goce.py text file, and return information about the contents.
    """
    times = []
    gitm_var = []
    print('Reading data from '+filename+'...')
    with open(filename, 'r') as File:
        lines = File.readlines()
        i = 0
        for line in lines:
            if i == 8:
                altitude = float(line.split()[0])
            if i >= 14:
                # Collect the observations by parsing each line at a time...
                splitted = line.split()
                times.append(
                    datetime(int(splitted[0]), int(splitted[1]), int(splitted[2]), int(splitted[3]), int(splitted[4]),
                             int(splitted[5])))
                gitm_var.append(float(splitted[6]))
            i += 1
    print('Complete!\n')
    return altitude, np.asarray(times), np.asarray(gitm_var)

def format_hour(hour_value):
    hour_str = str(hour_value)
    if len(hour_str) > 2:
        return f"{hour_str[:2]}:{hour_str[2:4]}:{hour_str[4:]}"
    return f"{hour_str.zfill(2)}:00:00"

def read_omni_txt_file(fname):
    """
    Performs the same function as 'readOMNI' below, but for a .txt file obtained with the wget command.
    :param fname: str
        The name of the file where the OMNI data is stored.
    :param numVars: int
        The number of variables in the data. Ensures that the data are parsed correctly.
    :return omniTimes: ndarray
        A 1D array of datetimes for the omni data.
    :return omniLabels: ndarray
        A 1D array of strings of each of the variables in the OMNI data.
    :return omniDataArray: ndarray
        A 2D array with all the OMNI data. The shape is nxm, where n is the number of time samples (each hour) and
        m is the number of variables collected at each time sample.
    """
    times = []
    data = []
    i = 0
    with open(fname, "r") as f:
        contents = f.readlines()
        numLines = len(contents)
        for line in contents:
            if i >= 9 and i < numLines-15:
                lineInfo = line.split()
                times.append( datetime(int(lineInfo[0]), 1, 1) + timedelta(days=int(lineInfo[1])-1) + timedelta(hours=int(lineInfo[2])) )
                data.append( float(lineInfo[-1]) )
            i += 1

    return times, data

def grab_omni_dst(dateStart, dateEnd):
    """
    Automatically obtain the Dst index from NASA OMNIWeb
    Parameters
    ----------
    dateStart: datetime
        The start date as a datetime object.
    dateEnd: datetime
        The end date as a dateteime object.
    Return:
    -------
    times: ndarray
        Datetimes values over the time interva.
    dst: ndarray
        The actual Dst values.
    """
    # Download the file.
    str1 = 'wget --post-data "activity=retrieve&res=hour&spacecraft=omni2&start_date='
    str2 = "&end_date="
    str3 = (
        "scale=Linear&ymin=&ymax=&view=0&charsize=&xstyle=0&ystyle=0&symbol=0&symsize=&linestyle=solid&table=0&"
        'imagex=640&imagey=480&color=&back=" https://omniweb.gsfc.nasa.gov/cgi/nx1.cgi -O'
    )
    startDate, endDate = dateStart.strftime("%Y-%m-%d"), (dateEnd + timedelta(days=1)).strftime("%Y-%m-%d") # Ceiling the end date for truncating later...
    fileToSave = 'omni_dst_'+startDate+'_'+endDate+'.txt'
    # override=True
    if os.path.isfile(fileToSave): # and override==False:
        print('Data already exists. Reading it in...')
    else:
        cmd = (
            str1
            + startDate.replace("-", "")
            + str2
            + endDate.replace("-", "")
            + "&vars=40&"
            + str3
            + " "
            + fileToSave
        )
        os.system(cmd)

    # Read the OMNI file and return the data:
    times, data = read_omni_txt_file(fileToSave)

    return np.asarray(times), np.asarray(data)

def getVarStr(varNum):
    masterList=[
        [0, 'Longitude'],
        [1, 'Latitude'],
        [2, 'Altitude'],
        [3, 'Rho'],
        [4, 'O(3P)'],
        [5, 'O2'],
        [6, 'N2'],
        [7, 'N(4S)'],
        [8, 'NO'],
        [9, 'He'],
        [10, 'N(2D)'],
        [11, 'N(2P)'],
        [12, 'H'],
        [13, 'CO2'],
        [14, 'O(1D)'],
        [15, 'Temperature'],
        [16, 'Vn(east)'],
        [17, 'Vn(north)'],
        [18, 'Vn(up)'],
        [19, 'Vn(up, O(3P))'],
        [20, 'Vn(up, O2)'],
        [21, 'Vn(up, N2)'],
        [22, 'Vn(up, N(4S))'],
        [23, 'Vn(up, NO)'],
        [24, 'Vn(up, He)'],
        [25, 'O_4SP_ +'],
        [26, 'NO +'],
        [27, 'O2 +'],
        [28, 'N2 +'],
        [29, 'N +'],
        [30, 'O(2D)+'],
        [31, 'O(2P)+'],
        [32, 'H +'],
        [33, 'He +'],
        [34, 'e -'],
        [35, 'eTemperature'],
        [36, 'iTemperature'],
        [37, 'Vi(east)'],
        [38, 'Vi(north)'],
        [39, 'Vi(up)']
    ]
    nums = np.array([element[0] for element in masterList])
    myInd = np.where(nums == int(varNum))[0][0]
    strs = np.array([element[-1] for element in masterList])
    return strs[myInd]

def combine_plot(var, file1, file2, tag1, tag2, file3=None, tag3=None, saveLoc=None):
    """
    Take two (or three) .txt files generated by thermo_goce.py (corresponding to different GITM runs over the same time period)
    and combine their results into a single figure. Generates and saves a .png to a user-defined location.
    Parameters
    ----------
    var: int
        The variable being plotted.
    file1: str
        The location/name of the first .txt file.
    file2:
        The location/name of the second .txt file.
    tag1: str
        A label for the first GITM run - used in the plot legend.
    tag2: str
        A label for the second GITM run - used in the plot legend.
    file3: str
    	The location/name of the third .txt file.
    tag3: str
    	A label for the third GITM run - used in the plot legend.
    saveLoc: str
        The full path (incl. the filename) with which to save the figure. Default is None, in which case the figure is
        saved to the current working directory
    Returns
    -------
    Nothing.
    """
    # Adjusting labels based on the variable being plotted:
    varString = getVarStr(var)

    # Read in the contents of the first file:
    altitude1, times1, gitm_var = read_thermo(file1)

    # Read in the contents of the second file:
    altitude2, times2, gitm_var2 = read_thermo(file2)
    
    if file3:
    	nfiles = 3
    	altitude3, times3, gitm_var3 = read_thermo(file3)
    else:
    	nfiles = 2

    # Make a plot of the time series for both runs - the plot includes the time series on the top and Dst over the sam
    # time period on the bottom. If there is a storm during the time period identified, it will be automatically
    # detected and its phases classified using the standards of Katus, et al. 2013 (doi:10.1029/2012JA017915)
    validInds = np.where((times1 >= times2[0]) & (times1 <= times2[-1]))[0]
    # Grab Dst during the time period:
    times_dst, dst = grab_omni_dst(times2[0], times2[-1])
    if np.min(dst) < -100:
        # Identify the beginning of the Main Phase and the Peak:
        peak_ind = np.argmin(dst)
        peak_time = times_dst[peak_ind]
        peak_dst = dst[peak_ind]
        #
        main_phase_start_ind = np.argmax(dst[np.where((times_dst >= peak_time - timedelta(hours=24)) & (times_dst <= peak_time))[0]])
        main_phase_start_time = times_dst[np.where((times_dst >= peak_time - timedelta(hours=24)) & (times_dst <= peak_time))[0]][main_phase_start_ind]
        true_main_phase_start_ind = np.where(times_dst == main_phase_start_time)[0][0]
        main_phase_start_dst = dst[true_main_phase_start_ind]

    #pdb.set_trace()
    fig, axs  = plt.subplots(2, 1, layout='constrained')
    axs[0].plot(times1[validInds], gitm_var[validInds], color='c', lw=4, label='GITM ('+tag1+')')
    axs[0].plot(times2, gitm_var2, color='m', lw=4, label='GITM ('+tag2+')')
    if nfiles == 3:
        axs[0].plot(times3, gitm_var3, color='forestgreen', lw=4, label='GITM ('+tag3+')')
    else:
        axs2 = axs[0].twinx()
        axs2.plot(times1[validInds], np.subtract(gitm_var[validInds], gitm_var2), color='b', linestyle='--')
        axs2.set_ylabel('Difference')
        axs2.tick_params(axis='y', labelcolor='b')
    axs[0].legend(loc='best')
    axs[0].set_xlabel('Time')
    axs[0].set_ylabel('Global Integral of '+varString)
    #
    validInds2 = np.where((times_dst >= times2[0]) & (times_dst <= times2[-1]))[0]
    axs[1].step(times_dst[validInds2], dst[validInds2], 'k-', lw=4)
    if np.min(dst) < -100:
        axs[0].axvline(x=main_phase_start_time, color='g', linestyle='-.', lw=4, label='Main Phase')
        axs[0].axvline(x=peak_time, color='r', linestyle='-.', lw=4, label='Peak')
        axs[1].axvline(x=main_phase_start_time, color='g', linestyle='-.', lw=4, label='Main Phase')
        axs[1].axvline(x=peak_time, color='r', linestyle='-.', lw=4, label='Peak')
        axs[1].legend(loc='best')
    axs[1].set_xlabel('Dst (nT)')
    # Save the figure:
    fig.suptitle('Evolution of Globally-Integrated '+varString+' at '+str(altitude1)+' km', fontsize=20, fontweight='bold')
    plt.savefig(saveLoc+'/gitm_time_series_comparison.png', dpi=300)

    if np.min(dst) < -100:
        # For the case that there is a storm, perform linear regression between the GITM variable and Dst:
        # 1 - Start to Peak:
        sp_inds = np.where((times_dst >= main_phase_start_time) & (times_dst <= peak_time))[0]
        dst_sp = dst[sp_inds]
        sp_len = len( sp_inds )
        downsampled_gitm_var = resample(gitm_var[validInds], sp_len)
        downsampled_gitm_var2 = resample(gitm_var2, sp_len)
        slope_1, intercept_1, r_value_1, p_value_1, std_err_1 = linregress(dst_sp, downsampled_gitm_var)
        slope_2, intercept_2, r_value_2, p_value_2, std_err_2 = linregress(dst_sp, downsampled_gitm_var2)
        if nfiles == 3:
            downsampled_gitm_var3 = resample(gitm_var3, sp_len)
            slope_3, intercept_3, r_value_3, p_value_3, std_err_3 = linregress(dst_sp, downsampled_gitm_var3)
        fig2, axs3 = plt.subplots(2, 1, layout='constrained')
        sample_dst = np.linspace(min(dst_sp)-20, max(dst_sp)+20, 10)
        axs3[0].scatter(dst_sp, downsampled_gitm_var, color='c', label='GITM ('+tag1+')')
        axs3[0].plot(sample_dst, intercept_1 + slope_1*sample_dst, 'c--', alpha=0.8, label='R$^2$='+str(np.round(r_value_1,2)))
        axs3[0].scatter(dst_sp, downsampled_gitm_var2, color='m', label='GITM ('+tag2+')')
        axs3[0].plot(sample_dst, intercept_2 + slope_2 * sample_dst, 'm--', alpha=0.8, label='R$^2$='+str(np.round(r_value_2,2)))
        if nfiles == 3:
            axs3[0].scatter(dst_sp, downsampled_gitm_var3, color='forestgreen', label='GITM ('+tag3+')')
            axs3[0].plot(sample_dst, intercept_3 + slope_3 * sample_dst, color='forestgreen', linestyle='--', alpha=0.8, label='R$^2$='+str(np.round(r_value_3,2)))
        axs3[0].set_xlabel('Dst (nT)')
        axs3[0].set_ylabel('Global Integral of ' + varString)
        axs3[0].legend(loc='best')
        axs3[0].set_title('Main Phase', fontweight='bold')
        axs3[0].set_xlim([sample_dst[0], sample_dst[-1]])
        if nfiles == 3:
            print('Main Phase Correlations:\n'
                  'GITM VAR 1 (Main Phase) - Slope: '+str(slope_1)+' Intercept: '+str(intercept_1)+' R-squared: '+str(r_value_1)+' P-value: '+str(p_value_1)+' Std Err: '+str(std_err_1)+'\n'
                  'GITM VAR 1 (Main Phase) - Slope: '+str(slope_2)+' Intercept: '+str(intercept_2)+' R-squared: '+str(r_value_2)+' P-value: '+str(p_value_2)+' Std Err: '+str(std_err_2)+'\n'
                  'GITM VAR 2 (Main Phase) - Slope: '+str(slope_3)+' Intercept: '+str(intercept_3)+' R-squared: '+str(r_value_3)+' P-value: '+str(p_value_3)+' Std Err: '+str(std_err_3)+'\n')
        else:
            print('Main Phase Correlations:\n'
                  'GITM VAR 1 (Main Phase) - Slope: '+str(slope_1)+' Intercept: '+str(intercept_1)+' R-squared: '+str(r_value_1)+' P-value: '+str(p_value_1)+' Std Err: '+str(std_err_1)+'\n'
                  'GITM VAR 1 (Main Phase) - Slope: '+str(slope_2)+' Intercept: '+str(intercept_2)+' R-squared: '+str(r_value_2)+' P-value: '+str(p_value_2)+' Std Err: '+str(std_err_2)+'\n')
        # 2 - Peak onwards:
        po_inds = np.where(times_dst >= peak_time)[0]
        dst_po = dst[po_inds]
        po_len = len(po_inds)
        downsampled_gitm_var_po = resample(gitm_var[po_inds], po_len)
        downsampled_gitm_var2_po = resample(gitm_var2, po_len)
        slope_1_po, intercept_1_po, r_value_1_po, p_value_1_po, std_err_1_po = linregress(dst_po, downsampled_gitm_var_po)
        slope_2_po, intercept_2_po, r_value_2_po, p_value_2_po, std_err_2_po = linregress(dst_po, downsampled_gitm_var2_po)
        if nfiles == 3:
            downsampled_gitm_var3_po = resample(gitm_var3, po_len)
            slope_3_po, intercept_3_po, r_value_3_po, p_value_3_po, std_err_3_po = linregress(dst_po, downsampled_gitm_var3_po)
        sample_dst_po = np.linspace(min(dst_po) - 20, max(dst_po) + 20, 10)
        axs3[1].scatter(dst_po, downsampled_gitm_var_po, color='c', label='GITM (' + tag1 + ')')
        axs3[1].plot(sample_dst_po, intercept_1_po + slope_1_po * sample_dst_po, 'c--', alpha=0.8,
                     label='R$^2$=' + str(np.round(r_value_1_po, 2)))
        axs3[1].scatter(dst_po, downsampled_gitm_var2_po, color='m', label='GITM (' + tag2 + ')')
        axs3[1].plot(sample_dst_po, intercept_2_po + slope_2_po * sample_dst_po, 'm--', alpha=0.8,
                     label='R$^2$=' + str(np.round(r_value_2_po, 2)))
        if nfiles == 3:
            axs3[1].scatter(dst_po, downsampled_gitm_var3_po, color='forestgreen', label='GITM ('+tag3+')')
            axs3[1].plot(sample_dst_po, intercept_3_po + slope_3_po * sample_dst_po, color='forestgreen', linestyle='--', alpha=0.8, label='R$^2$='+str(np.round(r_value_3_po, 3)))
        axs3[1].set_xlabel('Dst (nT)')
        axs3[1].set_ylabel('Global Integral of ' + varString)
        axs3[1].legend(loc='best')
        axs3[1].set_title('Recovery Phase', fontweight='bold')
        axs3[1].set_xlim([sample_dst_po[0], sample_dst_po[-1]])
        if nfiles == 3:
            print('Main Phase Correlations:\n'
              'GITM VAR 1 (Main Phase) - Slope: ' + str(slope_1_po) + ' Intercept: ' + str(
            intercept_1_po) + ' R-squared: ' + str(r_value_1_po) + ' P-value: ' + str(p_value_1_po) + ' Std Err: ' + str(
            std_err_1_po) + '\n'
                         'GITM VAR 1 (Main Phase) - Slope: ' + str(slope_2_po) + ' Intercept: ' + str(
            intercept_2_po) + ' R-squared: ' + str(r_value_2_po) + ' P-value: ' + str(p_value_2_po) + ' Std Err: ' + str(
            std_err_2_po) + '\n')
        else:
            print('Main Phase Correlations:\n'
                  'GITM VAR 1 (Main Phase) - Slope: ' + str(slope_1_po) + ' Intercept: ' + str(
            	intercept_1_po) + ' R-squared: ' + str(r_value_1_po) + ' P-value: ' + str(p_value_1_po) + ' Std Err: ' + str(
            	std_err_1_po) + '\n'
                         	'GITM VAR 1 (Main Phase) - Slope: ' + str(slope_2_po) + ' Intercept: ' + str(
            	intercept_2_po) + ' R-squared: ' + str(r_value_2_po) + ' P-value: ' + str(p_value_2_po) + ' Std Err: ' + str(
            	std_err_2_po) + 'GITM VAR 3 (Main Phase) - Slope: ' + str(slope_3_po) + ' Intercept: ' + str(
            	intercept_3_po) + ' R-squared: ' + str(r_value_3_po) + ' P-value: ' + str(p_value_3_po) + ' Std Err: ' + str(
            	std_err_3_po) + '\n')
        # Save the figure:
        plt.savefig(saveLoc + '/gitm_time_series_correlation_comparison.png', dpi=300)
    return 0

def resample(y, newLen):
    x = np.linspace(0, len(y)-1, len(y))
    f = InterpolatedUnivariateSpline(x, y)
    newX = np.linspace(0, len(y), newLen)
    newY = f(newX)
    return newY
#------------------------------------------------------------------------------
# SCRIPT USE
# example command line input: python thermo_compare.py 8 timefile_v1_112.txt timefile_vOvation_112.txt FTA OP

args = parse_args()

fileout = combine_plot(args.var[0], args.file1[0], args.file2[0], args.tag1[0], args.tag2[0], args.file3[0], args.tag3[0], args.loc)
