# Introduction

The Global Ionosphere Thermosphere Model (GITM) is a 3D model of the ionosphere and thermosphere that has been applied to Earth, Venus, Mars, Jupiter, and the moon Titan. The main repository for GITM can be found on [github](https://github.com/GITMCode/GITM). The first paper describing GITM is [here](https://www.sciencedirect.com/science/article/pii/S1364682606000071).

## Quick Start

If you have some expertise in downloading, compiling, and running codes, this is a good place to start.

Dependencies: fortran (assuming gfortran), mpi (mpich doesn't seem to work), make, perl (we are old school!), and python.

Download, configure, and run ('gfortran10' is for any gfortran version 10 or higher):
```bash
git clone https://github.com/GITMCode/GITM.git
cd GITM
./Config.pl -install -compiler=gfortran10 -earth
make
make rundir
cd run
mpirun -np 4 ./GITM.exe
./post_process.py
```

Then to produce a plot:
```bash
cd UA/data
python3 ../../../srcPython/gitm_plot_simple.py -alt=250 -var=3 *bin
```
This should produce two png files showing the density across the globe at about 250 km altitude.  This is a relatively simple plotter that will make cuts in altitude, latitude, and longitude. They are quite basic, but this should provide some framework for reading in GITM files and visualizing the results.  More complex plotters exist. Some of these are described in [this document](plotters.md).

## A Bit More Explaination

GITM is a fortran code that uses the message passing interface (MPI) to run on multiple processors, hence the need for a fortran compiler and MPI.  If the code compiles correctly, then the biggest issue is setting up the run.

GITM can create a run directory for you, so you don't have to worry about the executable and input files and everything being in the right place - it does it for you with the 'make rundir' command.  This will create a directory called 'run'. You can move this directory to where ever you want. Often the run directory is moved to a scratch disk and renamed. For example:
```bash
mv run /scratchdisk/mydirectory/run.gitm.event01
ln -s /scratchdisk/mydirectory/run.gitm.event01 .
```

Then, within this directory, there is a file called 'UAM.in', which is what GITM reads to set all of its parameters. This is an extremely brief description of settings within the file - to understand more, please read the manual! Within this file, there are several things that need to be altered when running GITM:

- Start time, given as year, month, day, hour, minute, second.
- End time, given as year, month, day, hour, minute, second.
- The grid, given as start and end latitude, start and end longitude, and the number of blocks (processors) to use in latitude and longitude. [See this grid description for more.](internals/grid.md).
- Save plots, telling what type of plots to output (3DALL by default) and how often (300s = 5 minutes). Ignore the restarts unless you know what you are doing.
- MHD_Indices - this is the IMF file that tells the high-latitude electrodynamics the IMF and solar wind.  There are instructions in the UAM.in file on producing a file.
- SME_Indices - this is the AE file that tells the auroral precipitation model how to work. There are instructions in the UAM.in file on producing a file.
- Dynamo - If you are running with anything better than 5 degrees by 5 degrees resolution, turn the dynamo on!

Once these things are set, you should be able to run GITM for the event of your choosing.

We provide a code in srcPython called gitm_makerun.py that will take the default UAM.in file, change the start/end times, download IMF and AE files, and flip switches that you may want. 

For more information on plotting, please see [the plotting description page](plotters.md).

## The More Thorough Start

Please see the [installation page](install.md) page for more info on how to
download, configure, and install GITM.

