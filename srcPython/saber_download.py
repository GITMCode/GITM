#!/usr/bin/env python

# Script for downloading SABER data for a specified time period.

# Top-level imports:
import sys, os
import numpy as np
import argparse
from datetime import datetime, timedelta
import requests
from bs4 import BeautifulSoup
import pooch

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description='Obtain TIMED/SABER data.')
    parser.add_argument('start', metavar='start', nargs=1, \
                        help='start date as YYYYMMDD')
    parser.add_argument('end', metavar='end', nargs=1, \
                        help='end date as YYYYMMDD')
    parser.add_argument('loc', metavar='location', nargs=1,\
                        help='directory to download data to')

    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------

def get_url_paths(url, ext='', params={}):
    response = requests.get(url, params=params)
    if response.ok:
        response_text = response.text
    else:
        return response.raise_for_status()
    soup = BeautifulSoup(response_text, 'html.parser')
    parent = [url + node.get('href') for node in soup.find_all('a') if node.get('href').endswith(ext)]
    return parent

def download_saber(dateStart, dateEnd, location):
    """
    Scrape the SABER website for Nitric Oxide volume emission data. If version 2.0 Level 2A data is available,
    preferentially grab that.

    Parameters
    ----------
    dateStart: str
        YYYYMMDD of the starting date to grab data.
    dateEnd: str
        YYYYMMDD of the ending date to grab data. This date is inclusive.
    location: str
        Place to save SABER data to.
    Returns
    -------
    Nothing. Just prints messages saying whether getting data was successful.
    """
    yearStart = int(dateStart[:4])
    yearEnd = int(dateEnd[:4])
    if not yearEnd == yearStart:
        raise ValueError('Ending year must equal starting year.')
    topUrl = 'https://data.gats-inc.com/saber/Version2_0/Level2A/'+str(yearStart)+'/'
    # Find the day numbers that you'll need to get data for:
    dateStartDatetime = datetime(yearStart, int(dateStart[4:6]), int(dateStart[6:]))
    dateEndDatetime = datetime(yearStart, int(dateEnd[4:6]), int(dateEnd[6:]))
    doy_start = dateStartDatetime.timetuple().tm_yday
    doy_end = dateEndDatetime.timetuple().tm_yday
    myDoys = np.linspace(doy_start, doy_end, num=doy_end-doy_start+1)
    myDoysInts = [int(element) for element in myDoys]
    # Loop over each individual day, and download ALL the files for that day. Arrange things so that subdirectories for
    # each day are made:
    for i in range(len(myDoysInts)):
        print('Downloading SABER files for day '+str(myDoysInts[i])+' of '+str(yearStart)+'...\n')
        currentURL = topUrl+str(myDoysInts[i])+'/'
        ext = 'nc'
        initial_result = get_url_paths(currentURL, ext)
        # Verify there are no duplicates
        result = np.unique(initial_result)
        # Check if there are data files corresponding to version 2; if so, ONLY grab those:
        combined = '_'.join(result)
        if '02.07.nc' in combined:
            newResult = []
            for fileName in result:
                if '02.07.nc' in fileName:
                    newResult.append(fileName)
        else:
            newResult = result
        fileNamesOnly = [element.split('/')[-1:][0] for element in newResult]
        # With the filenames obtained, grab all of them for this day and save them in the desired location, with a
        # sensible directory structure - do this with Pooch to avoid overwriting existing files:
        cacheLoc = location+'/'+str(yearStart)+'/'+str(myDoysInts[i])+'/'
        if os.path.isdir(cacheLoc) == False:
            os.makedirs(cacheLoc, exist_ok=True)
        for j in range(0,len(newResult)):
            fileStr = newResult[j]
            pooch.retrieve(
                url=fileStr,
                known_hash=None,
                fname=fileNamesOnly[j],
                path=cacheLoc,
                progressbar=True
            )
        print('SABER files for day '+str(myDoysInts[i])+' of '+str(yearStart)+' obtained!\n')
    if dateEnd == dateStart:
        print(
            'COMPLETE: SABER files for ' + dateStart + ' have been downloaded to: ' + location)
    else:
        print('COMPLETE: SABER files between '+dateStart+' and '+dateEnd+' inclusive have been downloaded to: '+location)

#------------------------------------------------------------------------------
# SCRIPT USE
# example command line input: python saber_download.py 20110805 20110808 /home/daabrand/Projects/TECHNO/Analysis1/20110805.Storm/SABER

args = parse_args()

#assuming first two args are the start/end dates, then info desired

start = args.start
if (not np.isscalar(start)):
    start = start[0]
end = args.end
if (not np.isscalar(end)):
    end = end[0]
loc = args.loc[0]

fileout = download_saber(start, end, loc)

# if __name__ == '__main__':
#     saveLoc = '/home/daabrand/Projects/TECHNO/Analysis2/20110805.Storm/v10x20/saber/'
#     dS = '20110803'
#     dE = '20110803'
#     download_saber(dS, dE, saveLoc)