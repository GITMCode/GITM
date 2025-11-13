#!/usr/bin/env python3

import re
import datetime as dt
import numpy as np
import requests

# The parsing of Dst function was downloaded from spedas
# The downloading of the data started from the spedas code and then
#     was pretty much rewritten.

def parse_dst_html(html_text, year=None, month=None):
    """
    Parses the HTML content to extract relevant information.

    Parameters
    ----------
    html_text : str
        The HTML content to parse.
    year : int, optional
        The year to consider while parsing the HTML content. If not provided, all years are considered.
    month : int, optional
        The month to consider while parsing the HTML content. If not provided, all months are considered.

    Returns
    -------
    dict
        A dictionary containing the parsed information.
    """
    times = []
    data = []
    # remove all of the HTML before the table
    html_data = html_text[html_text.find("Hourly Equatorial Dst Values") :]
    # remove all of the HTML after the table
    html_data = html_data[: html_data.find("<!-- vvvvv S yyyymm_part3.html vvvvv -->")]
    html_lines = html_data.split("\n")
    data_strs = html_lines[5:]
    # loop over days
    for day_str in data_strs:
        # the first element of hourly_data is the day, the rest are the hourly Dst values
        hourly_data = re.findall(r'[-+]?\d+', day_str)
        ## if the data is not complete for a whole day (which is typically the case for real time data):
        if len(hourly_data[1:]) != 24:
            ## if the data is not completely missing for a whole day:
            if len(hourly_data[1:]) != 3:
                for idx, dst_value in enumerate(hourly_data[1:]):
                    ## The kyoto website uses a 4 digit format.
                    remainder = len(dst_value) % 4
                    ## if the remainder is not zero, it can be either the regular case '-23' or
                    ## the ill case '-159999'. index by 0:remainder gives the correct -23 or -15
                    if remainder > 0:
                        times.append(
                            dt.datetime(int(year),
                                        int(month),
                                        int(hourly_data[0]),
                                        int(idx),
                                        30,
                                        0))
                        data.append(float(dst_value[0:remainder]))
                    ## if the remainder is zero, it can be either the regular case '-1599999' or
                    ## the ill case '9999...9999' with the number of nine being the multiple of 4.
                    ## we further test if the first four digits are 9999. If not, we simply index by
                    ## [0:4], which gives -159 in the regular case. Else, we ignore missing data.
                    elif dst_value[0:4] != '9999':
                        times.append(
                            dt.datetime(int(year),
                                        int(month),
                                        int(hourly_data[0]),
                                        int(idx),
                                        30,
                                        0))
                        data.append(float(dst_value[0:4]))
        ## if the data is complete for a whole day.
        else:
            for idx, dst_value in enumerate(hourly_data[1:]):
                times.append(
                            dt.datetime(int(year),
                                        int(month),
                                        int(hourly_data[0]),
                                        int(idx),
                                        30,
                                        0))
                data.append(float(dst_value))

    return (times, data)


def get_dst(
    times = None,
    datatypes=["final", "provisional", "realtime"],
    time_clip=True,
    prefix="",
    suffix="",
    no_download=False,
    local_data_dir="",
    download_only=False,
    force_download=False,
):
    """
    Loads Dst index data from the Kyoto servers.

    Parameters
    ----------
    times : list of str, required
        Time range of interest with the format ['YYYY-MM-DD','YYYY-MM-DD'] or
        to specify more or less than a day ['YYYY-MM-DD/hh:mm:ss','YYYY-MM-DD/hh:mm:ss'].
    time_clip : bool, optional
        Time clip the variables to exactly the range specified in the trange keyword.
        Defaults to True.
    remote_data_dir : str, optional
        The remote directory from where to load the Dst index data.
        Defaults to "http://wdc.kugi.kyoto-u.ac.jp/".
    suffix : str, optional
        The tplot variable names will be given this suffix.
        By default, no suffix is added.
    force_download: bool
        Download file even if local version is more recent than server version
        Default: False

    Returns
    -------
    list of str
        List of tplot variables created.

    Notes
    -----
    There are three types of Dst data available: final, provisional, and realtime.
    Usually, only one type is available for a particular month.
    is function tries to download final data, if this is not available then
    it downloads provisional data, and if this is not available then it downloads
    realtime data.

    Examples
    --------
    >>> from pyspedas.projects.kyoto import dst
    >>> dst_data = dst(times=['2015-01-01', '2015-01-02'])
    >>> print(dst_data)
    kyoto_dst
    """

    vars = []  # list of tplot variables created

    if times is None:
        print("Keyword times is required to download data.")
        return vars
    if (len(times) > 1):
        if (times[0] >= times[1]):
            print("Invalid time range. End time must be greater than start time.")
            return vars

    if local_data_dir == "":
        local_data_dir = "./"
    if local_data_dir[-1] != "/":
        local_data_dir += "/"

    remote_data_dir = "https://wdc.kugi.kyoto-u.ac.jp/"

    datatypes = ["final", "provisional", "realtime"]
    
    remoteNames = []
    localNames = []
    allTimes = []
    allDst = []
    for time in times:
        ymd = time.split('-')
        filename = ymd[0] + ymd[1] + '/index.html'
        didWork = False
        iType = 0
        while ((not didWork) and (iType < len(datatypes))):
            remote = remote_data_dir + \
                'dst_' + datatypes[iType] + '/' + \
                filename
            print(remote, didWork, iType)
            result = requests.get(url = remote)
            if (len(result.text) > 1000):
                didWork = True
                t, d = parse_dst_html(result.text, year=int(ymd[0]), month=int(ymd[1]))
                if (len(allTimes) < 1):
                    allTimes = t
                    allDst = d
                else:
                    allTimes = np.concatenate(allTimes, t)
                    allDst = np.concatenate(allTimes, d)
            else:
                iType = iType + 1

    data = {'times': allTimes,
            'dst': allDst}

    return data
        
