# Installation

If a fortran compiler & MPI are already installed on your system, you may wish
to skip down to [Getting the Code](#getting-the-code).

If you are unsure whether or not things are installed & set up correctly, a good
way to check is to print the version number with `gfortran --version` (change to
`ifort` if using the Intel compiler).  The MPI implementation can be checked
with `mpirun --version`. If there is an error, something is not configured
correctly.

## Requirements

At a minimum, you need:

- a Fortran compiler (gfortran, ifort, ifx, etc.)
- MPI (~~mpich~~, openmpi, mvapich, etc.)
- GNU Make (`make`)

<!-- abbreviation definition -->
*[MPI]: Message Passing Interface

There is no difference in the outputs between different compilers, however some
compilers may produce slightly faster executables than others. For example,
using ifort on [Pleiades](https://www.nas.nasa.gov/hecc/resources/pleiades.html)
is faster than using gfortran (gcc) or aocc. 

As GITM can run on as many (or few) CPU cores as you wish, it is possible to run
GITM on a laptop or workstation. This is recommended for development, as the
turnaround for test runs will be much faster. 

### Linux Install Dependencies

On a Ubuntu-based linux distribution, the following commands will download GNU
Make, the GCC Fortran compiler & Open-MPI:

```bash
sudo apt-get install gfortran libopenmpi-dev make
```

Similarly, Fedora-based distributions can use:

```bash
sudo yum install gcc-gfortran openmpi-devel
```

### MacOS Install Dependencies

Installing gfortran on MacOS is most easily accomplished with a package manager
like [Homebrew](https://brew.sh/) or [MacPorts](https://www.macports.org/).
Homebrew users will need to run:

```bash
brew install gfortran open-mpi
```

And MacPorts users can run:

```bash
sudo port install gcc[??] open-mpi
```
> Make sure to specify a version number in place of `??` above! Use 
> `port search --name --glob 'gcc*'`, or `port search gfortran` to see
> available Ports.


### Recommended HPC Modules

Feel free to contribute to this list of recommended modules if you have
experience with GITM on other systems:

- Pleiades
    - `comp-intel`
    - `mpi-hpe`
- Derecho
    - `intel`
    - `cray-mpich`

Consult the documentation for each system to see which modules are available and
how to load them. The modules you use to compile GITM *may* also need to be
loaded when running, so it may be wise to set these to your defaults, or at
least take note of what was used so they can be loaded in the job script.

## Getting the code

GITM is hosted on GitHub. The `main` branch is stable and updated as important
features are added (in the `develop` branch). The `main` branch is default and
no additional steps are needed to use the latest, stable, version. Just clone
the GITM repository, and all dependencies like the
[share](https://github.com/SWMFsoftware/share),
[util](https://github.com/SWMFsoftware/util), and
[electrodynamics](https://github.com/GITMCode/Electrodynamics) libraries will be
downloaded during configuration.

```shell
git clone git@github.com:GITMCode/GITM.git
cd GITM
```

!!! note
    Replace `git@github.com:GITMCode/GITM.git` with
    `https://github.com/GITMCode/GITM.git` if you don't have Github ssh keys set
    up.

All of the following steps assume you have not changed out of the `GITM/`
directory. Most will error if run from another location, but this will not break
anything! Simply `cd` back to `GITM/` and try again. 

## Configuring & Compiling

This step configures the planet, compiler, and some paths GITM needs to
work properly. To configure GITM for Earth using the `gfortran` compiler, run:

```bash
./Config.pl -install -earth -compiler=gfortran
```

!!! caution
     <!--#TODO> </!--> 
    Prior GITM versions required specifying `compiler=gfortran10` when using
    `gfortran>=10.0.0`, however this is no longer necessary.

The full list of available configuration options can be found by running
`./Config.pl -h`. A useful flag while developing is `-debug` which will print a
traceback if an error is encountered, as well as performing some additional
validity checks when running.

This should have downloaded all necessary external libraries and created some
extra files necessary to compile. To compile, simply run:

```bash
make
```

!!! note
    Compilation can often take a while. This can be done in parallel with
    `make -j`, which will use all available cores. To limit the number of cores
    to 8, for example, use `make -j8`

If this runs without error, GITM is ready to be run!

## Running the Code

This will have created the executable `src/GITM.exe`. This is not where we want
to run it from, however. To create a directory where the GITM executable is run
from and organize some input files, run the command:

```bash
make rundir
```

Now there should be a new directory in the GITM folder called `run/`. This
folder can be moved, copied, renamed, etc. without issue. 

!!! Note 
    If a study requires multiple runs of GITM, it can be most time-effective to
    copy or move the folder that was just created to somewhere a lot of data
    can be stored. For example, one could now run
    `cp -R run/ /path/to/scratch/storage/thisproject/run01`, and repeat for as
    many runs are necessary.

It is now time to do some runs! The folder we just created contains a `UAM.in`
file which calls for 4 CPU cores and simulates 5 minutes in December of 2002.
This all can, of course, be adjusted later. To start this run, we first need to
`cd` into the run directory, then tell MPI to run `GITM.exe` on 4 cores.

```bash
cd run/
mpirun -np 4 ./GITM.exe
```

Wait a moment... And if no errors are reported then congratulations! You have
now run GITM! 

The next steps include exploring how to [postprocess](postprocessing.md) the 
outputs, and modifying the [input files](common_inputs.md).
