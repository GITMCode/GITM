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
- Python3 (most any version should work, but we have been testing with versions in the 3.10+ range)
- Perl (for configuring the code)

<!-- abbreviation definition -->
*[MPI]: Message Passing Interface

On linux systems (including Windows Subsystem for Linux), gfortran is often used.  This is the most robustly tested compiler for GITM.  One problem with gfortran is that the gcc-10 and above version don't place well with MPI for some reason.  A flag has to be specified to make them play nice, and therefore there are two different compiler options for gfortran (-compiler=gfortran for older versions and -compiler=-gfortran10 for version 10+). There is a very good chance that you have 10+.

There is no difference in the outputs between different compilers, however some
compilers may produce slightly faster executables than others. For example,
using ifort on [Pleiades](https://www.nas.nasa.gov/hecc/resources/pleiades.html)
is faster than using gfortran (gcc) or aocc. 

As GITM can run on as many (or few) CPU cores as you wish, it is possible to run
GITM on a laptop or workstation. This is recommended for development, as the
turnaround for test runs will be much faster. Most developers of GITM have 8-core machines and can run the default test problem of 4 processors with no problems.  Most modern computers are capable of this now.

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

On some computer systems, there are multiple versions of compilers, MPI, and other things like python available.  They deal with these things by making "modules".  This means that in order to compile the code, modules need to be loaded.  For example, on NASA's Pleiades computer, this line needs to be put into your .cshrc or .bashrc or whatever dot file you are using for your setup:
```bash
module load comp-intel mpi-hpe
```

The modules that need to be loaded on different systems will differ. Feel free to contribute to this list of recommended modules if you have experience with GITM on other systems:

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
downloaded during configuration. The following command is assuming that you will be using GITM and not developing GITM.  

```
git clone git@github.com:GITMCode/GITM.git
cd GITM
```

!!! tip
    Replace `git@github.com:GITMCode/GITM.git` with `https://github.com/GITMCode/GITM.git` if you don't have Github ssh keys set up. If you are going to be developing GITM, you may take the time to fork the repository (including all of the branches!), and then substitute your forked repository location in the git clone command. If you are NOT doing development, just use this command.

All of the following steps assume you have not changed out of the `GITM/`
directory. Most will error if run from another location, but this will not break
anything! Simply `cd` back to `GITM/` and try again. 

## Configuring & Compiling

This step configures the planet, compiler, and some paths GITM needs to
work properly. To configure GITM for Earth using the `gfortran` compiler, run:

```bash
./Config.pl -install -earth -compiler=gfortran10
```

!!! warning 
    <!--#TODO> </!--> The above example assumes you are using gfortran
     version 10+. If your gfortran version is <10 (you should upgrade your system!), use `compiler=gfortran`

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
    `make -j`, which will use all available cores on your computer. To limit the number of cores
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
folder can be moved, copied, renamed, etc. without issue. If a study requires multiple runs of GITM, it can be most time-effective to copy or move the folder that was just created to somewhere a lot of data can be stored. For example, one could now run:
```bash
mv run /path/to/scratch/storage/thisproject/run01
ln -s /path/to/scratch/storage/thisproject/run01 .
```
The 'ln' command will just create a link to the run directory in the GITM directory, so you don't lose it. If you want to make yet another run directory (you can make as many as you want!), you could then do:
```bash
make rundir
mv run /path/to/scratch/storage/thisproject/run02
ln -s /path/to/scratch/storage/thisproject/run02 .
```

It is now time to do some runs! The folder we just created contains a `UAM.in`
file which calls for 4 CPU cores and simulates 5 minutes in December of 2002.
This all can, of course, be adjusted later. To start this run, we first need to
`cd` into the run directory, then tell MPI to run `GITM.exe` on 4 cores.

```bash
cd run/
mpirun -np 4 ./GITM.exe
```

Twiddle your thumbs for a moment... And if no errors are reported then congratulations! You have
now run GITM! 

The next steps include exploring how to [postprocess](postprocessing.md) the 
outputs, and modifying the [input files](common_inputs.md).
