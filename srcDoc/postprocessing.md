# Postprocessing

GITM stores outputs at a cadence specified in the input `UAM.in` file (in #SAVEPLOTS). When
running at each `DtOutput`, each processor writes its own (block) file for each output type
(`3DALL`, `2DANC`, etc.), plus the first processor writes a header file. This creates
an extremely large number of files on long runs across many cores. To make the
number of files more manageable and easier to read, they are post-processed into
`.bin` files, one for each output type at each `DtOutput`.  This is different than other codes that write all of the output for all timesteps into a single file (so that file can be huge), or that writes different files for different states.  GITM's post-processor writes all of the requested variables to one file for each time-step.

When running `make rundir`, a file is placed within the `run/` folder called
`post_process.py` which is the entrypoint for postprocessing. This turns the raw
GITM outputs from this:

```bash
> ls UA/data
    3DALL_t021221_000000.b0001  3DALL_t021221_000000.header  3DALL_t021221_000500.b0004
    3DALL_t021221_000000.b0002  3DALL_t021221_000500.b0001   3DALL_t021221_000500.header
    3DALL_t021221_000000.b0003  3DALL_t021221_000500.b0002   log00000002.dat
    3DALL_t021221_000000.b0004  3DALL_t021221_000500.b0003   run_information.txt
```

Into this:

```bash
> ls UA/data
    3DALL_t021221_000000.bin  3DALL_t021221_000500.bin  log00000002.dat  run_information.txt
```

## post_process.py Usage

`post_process.py` supports burst processing (process existing files & exit), or
can watch for new outputs for a specified duration and process new outputs as
they are created. The output `.bin` files can be moved to a remote system,
tarred, and/or exist alongside the raw outputs.

This means that there are essentailly two ways of running the post-processor:
- Once the run is completed, run it and post process all of the files at once.
- Run it after starting the simulation, so it runs in parallel to the simulation.  This post processes the output files as the code is running, so that the files should be post processed soon after GITM is done running.

The help message for `post_process.py`:

```
./post_process.py -h
    usage: post_process.py [-h] [-remotefile REMOTEFILE] [-user USER] [-server SERVER] [-dir DIR] [-sleep SLEEP]
                        [-totaltime TOTALTIME] [-v] [-norm] [-tgz] [-nc] [--combine] [--runname RUNNAME]

    Post process and (optionally) move model results.
    - This functions similar to pGITM.py, but can copy files to
    a remote location, or postprocess files as they are created.

    options:
    -h, --help            show this help message and exit
    -remotefile REMOTEFILE
                            File that contains info. about remote system
    -user USER            remote user name (default none)
    -server SERVER        remote system name (default none)
    -dir DIR              remote directory to use (default none)
    -sleep SLEEP          how long to sleep between loops in seconds, (default 300)
    -totaltime TOTALTIME  specify how long to run in total in hours, (default 0 - only run once)
    -v, --verbose         Run with verbose
    -norm                 don't remove any files
    -tgz                  tar and zip raw GITM file instead of process
    -nc                   Postprocess to netCDF files instead on '.bin'?
    --combine             If processing to netCDF, we can combine each timestep to a 
                          single file per output type. (ex: 3DALL.nc, etc.). 
                          Will not work without -nc or if using remote.
    --runname RUNNAME     When combining, this is prepended to the output type in the
                          resulting files: '[runname]_3DALL.nc'. 
                          Default is no descriptor.
```

!!! note "Post-processing Speed" 

    The post_process.py code is relatively slow - it was written to be more robust and to be clean and not fast.  If extremely fast post processing is desired, the older fortran post processor (driven by a shell script) can be used. This is not the recommended way of doing things, but it can be much faster. The old post processor executable can be created by running `make POST` from GITM's root folder (`GITM/`). If this is done before `make rundir`, it will automatically be copied. Otherwise it must be linked to `run/` manually.

    The python post-processor will use `PostGITM.exe` if it is found in `run/`. This is not the default behavior as some systems limit which programs can be run from login nodes (but don't seem to limit python yet!).

One of the benefits of post_process.py is that it can post process the files and then scp them over to another system. In order to do this, you need to have keyless entry setup, so that files can be transferred without entering a password. The machine name, user name, and directory to move the run into have to be passed into the post_process.py script, which is done with a 'remotefile'.  If the python code finds this file, it will read it and start moving files. Specifically, the post processor:
- checks to see if the directory on the remote system exists. If not, it tries to create it.
- post processes any files that need to be post processed in UA/data.
- pushes the UA/data/*bin files over to the remote system.
- checks one-by-one to see if the files exist on the remote system. If the file exists on the remote system, it removes it from the UA/data directory on the compute node.
- pushes supplemental files (like log files and run files) to the remote system.
- if the processing time has not expired, it waits a give amount of time (5 minutes by default) and loops. If the processing time has been exceeded, it exists.

### Arguments

#### remotefile

Path to a file with information on the remote system to copy files to. This is
useful when running on a remote system with storage limitations such that files
can be copied externally faster than they are created. Within the `remotefile`,
three lines are checked and should contain:

```
username
hostname
/path/to/remote/directory
```

This information is used in a rsync command to copy the files, with the form
`rsync [output files] username@hostname:/path/to/remote/directory`.

It is not necessary to use a `remotefile` and it is also possible to specify
these three strings as individual arguments.

#### user

Username to use when copying files to a remote location.

#### server

Server/hostname to move `.bin` files to. As the script uses `rsync` to copy
these files, it is possible to use a hostname alias from `~/.ssh/config` as the
server.

#### dir

Path on remote server to move post-processed files to. If the final directory
does not yet exist, it will be created.

#### sleep

This sets the time, in seconds, to wait in between checking for new files. This
is only used if `-totaltime` is set. The default value is 300 seconds, or 5
minutes.

#### totaltime

`-totaltime` sets the maximum amount of time this script should run for, in
hours. The default value is 0 hours, so the script will run once and exit. The
script will continue running until `totaltime` has elapsed, even if the run has
finished, so care should be taken to estimate the actual time a run will take.

#### norm

`-norm` means no-rm, and if specified, will prevent the automatic removal of the
raw outputs.

#### tgz

`-tgz` tells the script to not post process the files at all, but to tar and zip them together. They can then be post-processed on a different machine. This is somewhat CPU intensive and slow. The remote feature will work with these files, so that they can be automatically pushed to another system. (This was created when Pleiades would not allow post processing on head nodes. After using this for a while, we made the post processor code do what the fortran code does, so that this option is mostly moot now.)

#### nc

`-nc` creates NetCDF files as outputs, instead of `.bin`. 
This requires [PyITMN](https://github.com/GITMCode/PyITM) be installed.

#### combine

`--combine` will append subsequent times to the first NetCDF file rather than creating a new file for each time. Some caveats: this can only be used with NetCDF, each output *type* is in a different file (3DALL.nc and 2DANC.nc are both created), and the Python routines in PyITM/GITM do not (yet) support this. This was designed for users of xarray who want to make a lot of custon plots quickly, as it works very well with Dask.

#### runname

`--runname` tells the post-processor to pre-pend the given name to the merged output files when using `--combine`. For example, `--runname example_run01` will create a file called `example_run01_3DALL.nc`. It can only be used when combining, and thus also only when making NetCDF file(s). 

---

## Plotting

For information on plotting, please see [the plotting description page](plotters.md).


## pGITM

If the run is creating files faster than they are postprocessed, do not panic.

GITM would create a `csh` script before the introduction of this Python script
which is still available. `pGITM` and `PostGITM.exe` use a Fortran backend
instead of post_process.py which is slow. The Fortran code should be much faster than
Python.

To access these, run `make POST` from GITM's root folder. If this is run before
`make rundir`, they will be placed in `run/` without manual intervention. If
`run` has already been created and moved or renamed, they must be copied
manually. `post_process.py` will prefer to use the Fortran post-processor and
only use Scipy if it is not found.

`pGITM` will only run one time and requires `csh` which doesn't seem to be available on systems by default anymore. The world has moved on.
