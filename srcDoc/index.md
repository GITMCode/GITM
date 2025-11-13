# Introduction

The Global Ionosphere Thermosphere Model (GITM) is a 3D model of the ionosphere and
thermosphere that has been applied to Earth, Venus, Mars, Jupiter, and Saturn's moon
Titan. The main repository for GITM is available on
[GitHub](https://github.com/GITMCode/GITM). The first paper describing GITM is
[here](https://www.sciencedirect.com/science/article/pii/S1364682606000071).

## Quick Start

If you have are confident downloading, compiling, and running codes, this is a good
place to start. A more detailed explanation can be found on the
[installation page](install.md).

Dependencies: fortran (assuming gfortran), mpi (mpich doesn't seem to work), make, perl
(we are old school!), and python 3+ for postprocessing & plotting. The version of
anything else should not matter as long as it is relatively recent.

Download, configure, and run:

```bash
git clone https://github.com/GITMCode/GITM.git
cd GITM
./Config.pl -install -compiler=gfortran -earth
make
make rundir
cd run
mpirun -np 4 ./GITM.exe
./post_process.py
```

!!! tip
    GITM now automatically detects the version of gfortran being used. If you notice
    `Config.pl` is incorrectly detecting the version of gfortran for your system, try
    specifying it manually.
  
    The compiler version only matters for gfortran, and the specific version is not
    important - only whether is it $`\pm 10`$. So use `compiler=gfortran` for 
    gfortran versions 9 and below, and `compiler=gfortran10` for versions 10 and up.

Then to produce a plot:
```bash
cd UA/data
python3 ../../../srcPython/gitm_plot_simple.py -alt=250 -var=3 *bin
```
This should produce two png files showing the density across the globe at about 250 km
altitude.  This is a relatively simple plotter that will make cuts in altitude,
latitude, and longitude. They are quite basic, but this should provide some framework
for reading in GITM files and visualizing the results.  More complex plotters exist,
and some of these are described in [this document](plotters.md).

## A Bit More Explanation

GITM is a fortran code that uses the message passing interface (MPI) to run on multiple
processors, hence the need for a fortran compiler and MPI.  If the code compiles
correctly, then the biggest issue is setting up the run.

GITM can create a run directory for you, so you don't have to worry about the executable
and input files and everything being in the right place - it does it for you with the
`make rundir` command.  This will create a directory called 'run'. You can move this
directory to where ever you want. Often the run directory is moved to a scratch disk
and renamed. For example:

```bash
mv run /scratchdisk/mydirectory/run.gitm.event01
ln -s /scratchdisk/mydirectory/run.gitm.event01 .
```

Then, within this directory, there is a file called `UAM.in`, which is what GITM reads
to set all of its parameters. 

## Basic Inputs

This is an extremely brief description of a few settings within `UAM.in` - to understand
more, consult other pages in the documentation! Within this file, there are several
things that need to be altered when running GITM:

- Start time, given as year, month, day, hour, minute, second. Each on its own line.
- End time, given as year, month, day, hour, minute, second.
- The grid, given as start and end latitude, start and end longitude, and the number of
  blocks (processors) to use in latitude and longitude. [See this grid description for
  more](internals/grid.md).
- Save plots, telling what type of files to output (`3DALL` by default) and how often
  (300s = 5 minutes). The restarts can be ignored unless you know you will need them.
- MHD_Indices - This is the the path to a file which is used as inputs to the
  high-latitude electrodynamics, specifying IMF and solar wind.  There are instructions
  in the UAM.in file on how to produce this with the correct format.
- SME_Indices - this is the AE file that some auroral precipitation models need to
  work. There are instructions in the UAM.in file on producing a file.
- Dynamo - Settings for the low-latitude dynamo. If you are running with anything better
  than 5 degrees by 5 degrees resolution, turn the dynamo on!

Once these things are set, you should be able to run GITM for the event of your
choosing.

We provide a code in srcPython called gitm_makerun.py that will take the default UAM.in
file, change the start/end times, download IMF and AE files, and flip switches that you
may want. 

For more information on plotting, please see [the plotting description page](plotters.md).

## The More Thorough Start

Please see the [installation page](install.md) page for more info on how to
download, configure, and install GITM.

