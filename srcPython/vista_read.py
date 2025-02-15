#!/usr/bin/env python

from datetime import datetime, timedelta
import numpy as np
import h5py

def read_vista_file(filename):

    print('Reading VISTA file : ', filename)

    # Assume filename contains the date, since the time is not
    # stored in the file:
    
    iYr_ = filename.find('.hdf5') - 6
    year = 2000 + int(filename[iYr_ : iYr_ + 2])
    month = int(filename[iYr_ + 2 : iYr_ + 4])
    day = int(filename[iYr_ + 4 : iYr_ + 6])
    basetime = datetime(year, month, day, 0, 0, 0)
    print(' -> Date : ', basetime)

    f = h5py.File(filename, 'r')

    # vista data is in the following keys:
    vista_data = np.array(f['data']['VISTA'])
    f.close()
    
    nLats, nLons, nTimes = np.shape(vista_data)
    lats = np.arange(nLats) / (nLats - 1) * 180.0 - 90.0
    localtimes = np.arange(nLons) / (nLons - 1) * 24.0
    
    deltats = np.arange(nTimes) / nTimes * 86400.0
    times = []
    for deltat in deltats:
       times.append(basetime + timedelta(seconds = deltat)) 

    vista = {
        'tec' : vista_data,
        'lats' : lats,
        'localtimes' : localtimes,
        'times' : times}
    return vista



    
    
