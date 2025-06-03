# Postprocessing

GITM stores outputs at a cadence specified in the input `UAM.in` file. When
running, there is one file created by each processor, for each output type
(`3DALL`, `2DANC`, etc.), plus one header file at each `DtOutput`. This creates
an extremely large number of files on long runs across many cores. To make the
number of files more manageable and easier to read, they are post-processed into
`.bin` files, one for each output type at each `DtOutput`.

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

The help message for `post_process.py`:

```
./post_process.py -h
    usage: post_process.py [-h] [-remotefile REMOTEFILE] [-user USER] [-server SERVER]
                        [-dir DIR] [-sleep SLEEP] [-totaltime TOTALTIME] [-q] [-v]
                        [-norm] [-tgz]

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
    -totaltime TOTALTIME  specify how long to run in total in hours, (default 0 - only run
    once)
    -q                    Run with verbose turned off
    -v                    Run with verbose
    -norm                 don't remove any files
    -tgz                  tar and zip raw GITM file instead of process
```

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

`-tgz` tells the script to tar and g-zip the outputs to save space. This is
somewhat CPU intensive.

---


## pGITM

If the run is creating files faster than they are postprocessed, do not panic.

GITM would create a `csh` script before the introduction of this Python script
which is still available. `pGITM` and `PostGITM.exe` use a Fortran backend
instead of post_process.py which uses Scipy. The Fortran code may be faster than
Python.

To access these, run `make POST` from GITM's root folder. If this is run before
`make rundir`, they will be placed in `run/` without manual intervention. If
`run` has already been created and moved or renamed, they must be copied
manually. `post_process.py` will prefer to use the Fortran post-processor and
only use Scipy if it is not found.

`pGITM` will only run one time and required `csh`.
