# Quick Start

To get started, here is what you need:

## Requirements

at a minimum, you need:

- git
- a fortran compiler (gfortran, ifort, ifx, etc.)
- MPI (~~mpich~~, openmpi, mvapich, etc.)
- perl (you probably have this)

## Getting the code


do this:

```shell
git clone git@github.com:GITMCode/GITM.git
cd GITM
```

!!! note 
    Replace `git@github.com:GITMCode/GITM.git` with
    `https://github.com/GITMCode/GITM.git` if you don't have ssh set up!

## Configuring

Use:

```bash
./Config.pl -install [...]
```
