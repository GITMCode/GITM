#!/bin/python3

import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import numpy as np
import h5py
import argparse
from gitm_routines import *


### Examples of usage:

# For a faceted plot of VISTA TEC, every 3 time-steps:
# python vista_gitm_plots.py --start_time 20130317 --end_time 20130317T02:00 --faceted_plots --plot_every_n 3


## To plot VISTA & GITM together:
# python vista_gitm_plots.py --start_time 20130317T13:00 --end_time 20130317T22:00 \
#                            --gitm_path /nobackup/Gitm/Runs_alb/advection_runs/20130317_on/ \
#                            --save_plots_to gitm_vs_vista_20130317


# Has some requirements, unfortunately. None too crazy.
# run:
# $ pip install numpy h5netcdf h5py xarray matplotlib pandas



def main(start_time, end_time, 
        vista_path,
        faceted_plots = True, 
        plot_every_n = 30, 
        gitm_path=None,
        save_plots_to='.'):

    if gitm_path is not None:
        faceted_plots = False

    # sanity check the dates...
    days = pd.date_range(start_time.date(), end_time.date(), freq='D')
    # Check if the last time is at 00:00:
    if end_time == pd.to_datetime(end_time.date()):
        days = days[:-1]

    vista_ds = load_vista_data(vista_path, days)
    vista_ds = vista_ds.where((vista_ds.time < end_time)
               & (vista_ds.time > start_time), drop=True)
    
    # then the faceted plots...
    if faceted_plots:
        vista_ds.isel(time=slice(0,-1,plot_every_n)).tec.plot(col='time', 
                                                              col_wrap=4,
                                                              robust=True,
                                                              cmap='plasma',
                                                              x='localtime',
                                                              subplot_kws={'ylabel': 'Latitude',
                                                                            'xlabel':'Local Time'},
                                                              cbar_kwargs={'label':'TEC (TECUnits)'})
        save_fig_path = os.path.join(save_plots_to, f'vista_{start_time.date()}-{start_time.time()}.png')
        plt.savefig(save_fig_path)
        plt.close()
        print('Saved plot:  ', save_fig_path)
        return
    
    else:
        if gitm_path:
            if '3DALL' not in gitm_path:
                gitm_path = os.path.join(gitm_path, '3DALL')
            headers = read_gitm_headers(gitm_path)
            if any([start_time < headers['time'][0], end_time > headers['time'][-1]]):
                raise ValueError('GITM data does not cover the full time range of the requested times')
            # Only use the GITM files closest to the VISTA outputs
            use_gitm_idxs = []
            for vista_time in vista_ds.time.values:
                use_gitm_idxs.append(np.argmin(np.abs((pd.to_datetime(headers['time'])
                                                      - vista_time).total_seconds())))

    ## Then go make the plots:

    if not os.path.exists(save_plots_to):
        print('creating directory: ', save_plots_to)
        os.makedirs(save_plots_to)

    for itime, vista_time in enumerate(vista_ds.time.values):
        print('Plotting ', vista_time)

        vista_time = pd.to_datetime(vista_time)

        vmin = None
        vmax = None

        if gitm_path:
            # read altitude & [e-]
            gitm_data = read_gitm_one_file(headers['filename'][use_gitm_idxs[itime]], [2, 34])
            # integrate the TEC
            tec_array = np.trapz(gitm_data[34], gitm_data[2])
            # get rid of ghost cells, convert to TECUnits
            tec_array = tec_array[2:-2, 2:-2].T / 1e16
            # now we need to shift the TEC from lon -> localtime
            tec_array = np.roll(tec_array, 
                                int((vista_time.hour) * 15/360 * (headers['nLons'])))

            #consistent colorbar limits btwn VISTA & GITM
            vmin = np.min(tec_array)
            vmax = np.max(tec_array)

            # Make figure, plot GITM:
            fig, ax = plt.subplots(3, 1, figsize=(10, 12), height_ratios=[1, 1, 0.1])
            ax[1].imshow(tec_array, cmap='plasma', vmin=vmin, vmax=vmax,
                        extent=[0, 24, -90, 90], aspect='auto')
            ax[1].set_title('GITM', fontsize='x-large')

            # These are not needed if GITM is not plotted:
            ax[0].set_title('VISTA', fontsize='x-large')
            ax[1].set_ylabel('Latitude')
            ax[1].set_xlabel('Local Time')


        else:
            fig, ax = plt.subplots(2, 1, figsize=(10, 7), height_ratios=[1, 0.08])

        # Plot VISTA:
        p = ax[0].imshow(vista_ds.isel(time=itime).tec.values, cmap='plasma',
                         vmin=vmin, vmax=vmax, extent=[0, 24, -90, 90], aspect='auto')
        
        fig.colorbar(p, cax = ax[-1], orientation='horizontal', label='TEC (TECUnits)',
                     extend='both')

        fig.suptitle('TEC at ' + vista_time.strftime('%Y-%m-%d %H:%M:%S'), fontsize='xx-large')
        ax[0].set_ylabel('Latitude')
        ax[0].set_xlabel('Local Time')

        save_fig_path = os.path.join(save_plots_to, f'vista_{str(itime).rjust(4,'0')}.png')
        plt.tight_layout()
        plt.savefig(save_fig_path)
        plt.close()
        print('Saved plot:  ', save_fig_path)

    return

   
def vista_grid(grid_path=None, vista_path=None):    
    # Load the grid
    if grid_path is None and vista_path is not None:
        grid_path = os.path.join('/',*vista_path.split('/')[:4], 'grid.hdf5')

    f = h5py.File(grid_path, "r")
    data = f["data"]
    
    lats = data['latitude']
    lts = data['local_time']
    
    lats = np.unique(lats)
    lts = np.unique(lts)
    
    f.close()

    return lats, lts


def load_vista_data(vista_path, days):

    lats, lts = vista_grid(vista_path=vista_path)
    
    # Load VISTA Data:
    vista_ds = []
    for day in days:
        print('loading ', day) 
    
        tmpds = xr.open_mfdataset(os.path.join(vista_path, 
                                               str(day.year),
                                               f"VISTA_{str(day.year)[2:]}{str(day.month).rjust(2, '0')}{str(day.day).rjust(2, '0')}.hdf5"),
                                  engine='h5netcdf', 
                                  group="/data", 
                                  phony_dims='access')
    
        ds_today = xr.Dataset()
        ds_today['lat'] = ('lat'), lats
        ds_today['localtime'] = ('localtime'), lts
        ds_today['time'] = ('time'), pd.date_range(day, 
                                                   day + pd.Timedelta('1 day'),
                                                   periods=tmpds.sizes['phony_dim_2']+1, 
                                                   inclusive='left')
        
        ds_today['tec'] = ('lat', 'localtime', 'time'), tmpds.VISTA.data
        vista_ds.append(ds_today)
    
    vista_ds = xr.concat(vista_ds, dim='time')

    return vista_ds


def parse_args():

    parser = argparse.ArgumentParser(description='Plot VISTA vs GITM TEC')

    parser.add_argument('--start_time', type=str, 
                        help='Start time in format YYYYMMDD[THH:MM], bracketed portion optional')
    parser.add_argument('--end_time', type=str,
                        help='End time in format YYYYMMDD[THH:MM], bracketed portion optional')
    parser.add_argument('--vista_path', type=str,
                        default = '/backup/Data/VISTA/VISTA_Database/HDF5/dtecrm/',
                        help='Path to VISTA data. Default is "/backup/Data/VISTA/VISTA_Database/HDF5/dtecrm/"')

    parser.add_argument('--faceted_plots', action='store_true',
                        help='If set, will create faceted plots, default is to make one file per time.'
                        ' (does not work w GITM TEC, only VISTA - will be False if gitm_path is set)')
    parser.add_argument('--plot_every_n', type=int, default=30,
                        help='If faceted_plots is set, will plot every n-th time step, default is 30')

    parser.add_argument('--gitm_path', type=str,
                        help='Path to GITM data, optional. If not set, will not plot GITM data')

    parser.add_argument('--save_plots_to', type=str, default='.',
                        help='Path to save plots to, default is to put them in the current directory')
    
    args = parser.parse_args()
    return args

if __name__ ==  '__main__':
    
    args = parse_args()

    args.start_time = pd.to_datetime(args.start_time)
    args.end_time = pd.to_datetime(args.end_time)
    print(f"Plotting from {args.start_time} to {args.end_time}")


    main(args.start_time, args.end_time,
        args.vista_path, args.faceted_plots, args.plot_every_n, args.gitm_path, args.save_plots_to)

