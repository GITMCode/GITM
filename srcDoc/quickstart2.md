# Quick Start (OLD)

!!! note "this is soooooo old"

    Will be dramatically thinned out shortly...

## Extracting the code from a tar file

Create a new and empty directory, and open the tar file you received,
e.g.:

    mkdir Gitm 
    cd Gitm 
    mv ../gitm.tgz . 
    tar -xvzf gitm.tgz 

## Checking out the code with CVS {#cvs.sec}

If CVS (Concurrent Versions System) is available on your computer and
you have an account on the CVS server machine herot.engin.umich.edu, you
can use CVS to install the current or a particular version of the code.
First of all have the following environment variables:

    setenv CVSROOT UserName@herot.engin.umich.edu:/CVS/FRAMEWORK
    setenv CVS_RSH ssh

where UserName is your user name on herot. Here it is assumed that you
use csh or tcsh. Also put these settings into your .cshrc file so it is
automatically executed at login. Alternatively, use

    CVSROOT=UserName@herot.engin.umich.edu:/CVS/FRAMEWORK 
    export CVSROOT 
    CVS_RSH=ssh 
    export CVS_RSH

under sh, ksh and bash shells, and also put these commands into your
.bashrc or .profile file so it is automatically executed at login.

Once the CVS environment variables are set, you can download the current
(HEAD) version of the GITM distribution with

    cvs checkout GITM2

If you want a particular version, use

    cvs checkout -r v2_0 GITM2

where v2_0 is the it tag associated with the version. To download bug
fixes or new features, the

    cvs update 

command can be used. See `man cvs` for more information.

A lot of times, you don't really want the `GITM2` directory to stay that
name, since you might download a couple different version (maybe one for
development and one for runs). Therefore, typically you will:

    mv GITM2 GITM2.Development

## Configuring and Making GITM

In order to compile GITM, you have to configure it first. The configure
script is inherited from the Space Weather Modeling Framework. There are
two primary reasons you need to do the configure: (1) put the right
`Makefile` in the right place, specifying the compiler and the version
of MPI that you will use to link the code; (2) put the right MPI header
in the right place. It also does some things like hard-codes the path of
the source code into the `Makefile`. Currently the configure script is
not capable of detecting the system and compilers available. Some
examples for commonly used set-ups are shown below. Make sure that your
.cshrc or .bashrc file is set up to detect the appropriate compilers
before attempting to install GITM.

Installing on Nyx:

    ./Config.pl -install -compiler=ifortmpif90 -earth

Installing with an Intel compiler and OpenMPI (such as Pleiades):

    ./Config.pl -install -compiler=ifort -earth

Installing a computer with gfortran and OpenMPI:

    ./Config.pl -install -compiler=gfortran -earth

Installing on a computer with gfortran and not using MPI:

    ./Config.pl -install -compiler=gfortran -earth -nompi

Sometimes people have a hard time with the ModUtilities.F90 file. If you
have errors with this file, try (for example):

    ./Config.pl -uninstall
    ./Config.pl -install -compiler=gfortran -earth -noflush

Don't forget, after configuring the Makefiles, you must still compile
the code!

    make
    make test_earth
    make install

## Running the Code

GITM requires a bunch of files to be in the right place in order to run.
Therefore, it is best to use the makefile to create a run directory:

    make rundir
    mv run myrun

where `myrun` can be whatever you want. I will use ` myrun` as an
example. You can actually put this directory where ever you want. On
many systems (such as nyx), there is a `nobackup` scratch disk that you
are supposed to use for runs, instead of your home directory. If you
need to ensure that your home directory doesn't use too much space,
moving the run directory onto a disk with more free space can solve the
problem:

    make rundir
    mv run /nobackup/myaccount/gitm/myrun
    ln -s /nobackup/myaccount/gitm/myrun .

This creates a shortcut to the `myrun` directory location on `nobackup`
in your GITM working directory. It allows you to treat the run directory
as if it were a local directory, but it isn't! It also means that you
don't have to compile and install GITM on the scratch disk, where
program storage may not be allowed.

Once you have created the run directory, you can run the default
simulation, by:

    cd myrun
    mpirun -np 4 GITM.exe

Or, if your system uses Mpiexec:

    cd myrun
    mpiexec ./GITM.exe

This, hopefully should run GITM for Earth for 5 minutes. If it doesn't
work, then you might have mpi set up incorrectly. The default is to
allow you to run 4 blocks per processor, and the default ` UAM.in` file
is set up for 4 blocks, so you could try just running GITM without mpi,
just to see if it works at all:

    ./GITM.exe

If that doesn't work, then it probably didn't compile correctly.
Hopefully, it just worked!

## Post Processing {#post_process.sec}

GITM, by default, produces one file per block per output. If you are
outputting often and you are running with many blocks, you can produce a
huge number of files. To post process all of these files, simply:

    cd UA
    ./pGITM

This merges all of the files for one time period, for one file type into
the same file. You can actually running this while the code is running,
since GITM doesn't use old files, unless you are using the APPENDFILE
option. As implied by the option's name, APPENDFILE opens an existing
file and appends the most recent data to it. This feature is typically
used only when running a satellite track though GITM. More information
on the APPENDFILE option is located in
Chapter [\[input.ch\]](#input.ch){reference-type="ref"
reference="input.ch"}
Section [\[def_out.sec\]](#def_out.sec){reference-type="ref"
reference="def_out.sec"}.

If you are NOT using satellites and NOT using APPENDFILE, then you are
free and clear to use pGITM as often as you want during a run. To avoid
deleting a file that GITM is currently writing to, it is recommended
that pGITM be run as part of a script in which there is a five minute
pause between executions.

Another useful script, when running GITM on another system, is given in
the block below. It occasionally executes an rsync between the computer
that you run GITM on and your home computer. This allows you to bring
over and evaluate the output files as they become available. To do this,
execute pGITM at a set cadence (60 seconds in the example) while GITM is
running and then rsync the remote and home directories (excluding all
unprocessed files). Finally, remove the processed and rsynced files to
prevent the remote directory from filling up. Be sure to replace
yourname@home.computer with your user and computer names! This is a very
simple, but very useful script.

    #!/bin/csh
    rm -f stop
    set LOC=$1

    if (-f remoteloc) set LOC=`cat remoteloc`

    while (!(-f stop))
      rsync -vrae ssh log.* UAM.* imf* yourname@home.computer:$LOC
      cd UA ; ./pGITM ; rsync --exclude '*.[bsh][0123ae][0123456789at]*' -vrae ssh d
    ata yourname@home.computer:$LOC ; cd ..
      sleep 60
    end

## The Code Won't Compile!!

I'm sorry. I tried to make this work on many different platforms, but
sometimes machines are very specific, and it just doesn't work out of
the box. Here are some ideas on how to quickly get this thing compiling.

### Can't find the right Makefile.whatever

If make does not work, then there is probably a problem with not finding
the FORTRAN 90 compiler. The platform and machine specific Makefiles are
in srcMake. If you type:

    uname 
    ls srcMake

If you don't see a file named something like Makefile.uname (where uname
is the output of the uname command), then you will have to build a
proper general Makefile.

You will need a little a little information about your computer, like
what the mpif90 compiler is called and where it is located. Take a look
at srcMake/Makefile.Linux, and try to figure out what all of the flags
are for your system. Then create a srcMake/Makefile.uname with the
correct information in it.

### The compiler doesn't recognize flag -x

You have an operating system that is recognized, but probably a
different compiler. In the srcMake/Makefile.uname file (where uname is
the output of the uname command), there is a line:

OSFLAGS = -w -dusty

You need to change this line to something more appropriate for your
compiler. Try deleting the flags and compile. If that doesn't work, you
will have to check the man pages of your compiler.

### src/ModHwm.90 doesn't compile

Certain versions of gfortran (4.6 and later) may give the following
error:

    src/ModHwm.f90:168.22:

            call HWMupdate(input,last,gfs,gfl,gfm,gvbar,gwbar,gbz,gbm,gzwght,glev,u
                          1
    Error: Dummy argument 'ebz' of procedure 'hwmupdate' at (1) has an attribute
     that requires an explicit interface for this procedure

    src/ModHwm.f90:168.22:

            call HWMupdate(input,last,gfs,gfl,gfm,gvbar,gwbar,gbz,gbm,gzwght,glev,u
                          1
    Error: Dummy argument 'ebz' of procedure 'hwmupdate' at (1) has an attribute
     that requires an explicit interface for this procedure

This is caused by the inputs in HWM. The latest incarnations of gfortran
don't allow optional inputs that are not declared. More information
about this can be found at:

    http://cosmocoffee.info/viewtopic.php?p=5136

A solution to this problem is currently being sought.
