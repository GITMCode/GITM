# GITM

This is the home of the Global Ionosphere/Thermosphere Model (GITM).

GITM has been developed in fortran-90. It has been tested with gfortran
on Linux and MacOS as well as ifort on NASA's Pleiades computer.

The latest version of the PDF documentation can be accessed at:

[![Build Documentation](https://github.com/GITMCode/GITM/actions/workflows/build-docs.yml/badge.svg)](https://github.com/GITMCode/GITM/actions/workflows/build-docs.yml)

For a more complete documentation, see [GITM's Read the Docs Page](https://gitm.readthedocs.io).

---

<details>
<summary><b>GITM's default branch recently changed names. To access the latest features, please update your local refs.</b></summary>

From your local clone of this repository, run the following commands to update the name of the default branch.

```sh
git branch -m master main
git fetch origin
git branch -u origin/main main
git remote set-head origin -a
```

Optionally, run the following command to remove tracking references to the old branch name.

```sh
git remote prune origin
```

</details>

## Quick Start

> GITM needs MPI to work properly. This can be installed with your system's package manager, or loaded as a `module` on an HPC system.
> See [the installation guide](https://gitm.readthedocs.io/en/latest/quick_start/#requirements) for more information.

1\. Clone the repository and cd into the repo folder

```shell
git clone git@github.com:GITMCode/GITM.git
```

Substitute the URL with the https link from the "Code" button above if you do not have SSH keys set up.

2\. Go into the repo directory

```shell
cd GITM
```

3\. Configure the Fortran compiler and download external electrodynamics library
(gfortran versions 10+)

```shell
./Config.pl -install -earth -compiler=gfortran10
```

The biggest issue with the above command is that it assumes that you
have the gfortran (version 10+) compiler and things like mpif90 work
ok.  If you don't have gfortran and mpif90, then you need to get these
things for your computer.  If you have version 9 or before for gfortran,
you can do:

```shell
./Config.pl -install -earth -compiler=gfortran
```

In theory, Mars, Venus, Titan, and LV-426 should work.  These are in
various states of completion, so I wouldn't count on them being
perfect.

If running on Pleiades, you need to have these
in your start-up script (.cshrc, .bashrc, etc):

```
module load comp-intel
module load mpi-hpe
```

And you can use the below line below to configure the code:

```shell
./Config.pl -install -earth -compiler=ifort
```

4\. Make the binary:

```shell
make
```

5\. Creates a run directory that has all of the input files:

```shell
make rundir
```

6\. Go into the run directory:

```shell
cd run
```

7\. Run the code:

```shell
mpirun -np 4 ./GITM.exe
```

GITM reads in a file called `UAM.in`, which sets the configuration of
the simulation. The default `UAM.in` file has 2 lat blocks and 2 lon
blocks with 9 x 9 cells each, so the default resolution is 180 (deg
lat) / (2 \* 9) = 10 deg lat, by 360 (deg lon) / (2 \* 9) = 20 deg
lon. See below for how to set the resolution.

8\. Post process the output files by running:

```shell
./post_process.py
```

> This can be called from `run/` and will postprocess the files in UA/output.
> It has functionality to copy files to a remote location, monitor the
> output folder throughout a run, and more. Run `./post_process.py --help` to
> see available options.

> The legacy postprocessors are still available, but are not built by default. To build
> PostProcess.exe, run `make POST`. The csh script can be found at `src/pGITM`.

9\. Go into the output directory:

```shell
cd data
```

10\. Make some plots with an old plotter:

```shell
../../../srcPython/plot_model_results.py -var=3 -alt=120 3DALL_t021221_000500.bin
```

Then look at the png file that is created.  You can use a `-h` to see
how to run this code.

10b\. A more advanced plotter is available through aetherpy. This is a bit more
complicated, since you need to install aetherpy. If you don't use python
much, this is harder. Here is how to do this:

```shell
cd <directory where you started from>
git clone https://github.com/AetherModel/aetherpy
ls
```

(you should see 2 directories: GITM and aetherpy)

```shell
cd aetherpy
git checkout develop
```

(install the aetherpy libraries)

```shell
python setup.py develop --user
```

10c\. Test out the new plotter:

```shell
cd <directory where you started from>
cd GITM/run/UA/data
../../../srcPython/run_plot_model_results.py -var=34 -alt=300 3DALL_t021221_000500.bin
```

(look at the beautiful plot)

## Contributing

1. Please feel free to e-mail the development team to suggest ideas.

2. Please feel free to open an issue on github.  The development team
gets these issues and will review them.  We can then reach out to you
to figure out how to incorporate them.

3. Please feel free to fork this repository, make changes as you see
fit, and do a pull request.  Your suggested changes will be reviewed
and incorporated if they fit.

See [CONTRIBUTING.md](.github/CONTRIBUTING.md) for more details.

## External Codes

There are a number of external codes that are not developed by the GITM
team.  For example:

1. APEX - the magnetic coordinate system in the code.  Developed at
NCAR. The IGRF code that comes with it is also not developed at UM.

2. Many electrodynamics models (Weimer, Newell's Ovation Prime,
Mitchell's Ovation SME, others in the 
[GITMCode/Electrodynamics](https://github.com/GITMCode/Electrodynamics)
repository).

3. MSIS and IRI, which are in util/EMPIRICAL/srcUA. MSIS is used as a
lower BC at Earth and was developed at NRL. IRI is used to initialize
the code.

4. The horizonal wind model (HWM) is used as a lower BC at Earth and
was developed at NRL.

## Setting the Resolution

GITM uses a 2d domain decomposition with blocks. This means that in
each direction, you set the number of cells in each block in the
src/ModSize.f90 file. We almost always leave this as 9 cells in both
the latitude and longitude direction, since the math is easy with this
number. You can then set the resolution by asking for the number of
blocks you want in the `UAM.in` file.

If you wanted a grid that is (for example) 1 deg (lat) by 5 deg (lon),
you would need 20 blocks (9 x 20 = 180 cells) in latitude and 8 blocks
(9 x 8 = 72 cells) in longitude.  You need a total of 160 processors for this
simulation. You can change the `UAM.in` file for these number of cells.
To get a simple 5 deg x 5 deg resolution, you need 4 blocks in lat
and 8 blocks in longitude, or 32 processors.  If you have fewer
processors, you can change the parameters in `src/ModSize.f90` and adjust the
number of cells in each block to compensate. For example, you have 8
processors, so you can adjust `src/ModSize.f90` to have 18 cells in lat and
lon, then ask for 2 (lat) x 4 (lon) blocks to get 5 deg x 5 deg resolution.

See the 
[grid page](https://gitm.readthedocs.io/en/latest/grid/#horizontal-resolution)
in the documentation for more details and some examples.
