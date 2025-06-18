# Setting the Grid

Setting the grid resolution in GITM is not very complicated, but it does involve
some thought. There are a few variables that control this. 

In essence, GITM's resolution is dictated by the number of CPU cores it is run
on. The number of cores must equal the product of the number of blocks in
longitude and latitude, which are both set in `UAM.in`. The resolution in the
longitudinal and latitudinal directions are, by default, 9 times the number of
blocks in each direction.

## Important Variables

Each CPU core that GITM is run on is considered a **block**. Within each block,
the grid is further divided into individual grid **cells**. The number of cells
per block is set before compiling, and the number of blocks to model is decided
at runtime. The number of blocks per core cannot change from 1. How many cores
GITM is run on is specified in the `UAM.in` file, and the number of cells within
each block is specified in `src/ModSize.f90`.

### Blocks (UAM.in)

When GITM is first creating its grid, the region of interest (usually the entire
globe) if first divided into blocks. Each block is a region that is modeled on a
single CPU core. This is the most common variable to adjust when trying to
change the resolution.

These values are set in the `UAM.in` file and read at runtime. So the same
executable can be used for many different resolutions without needing to
recompile!

For example, the initial settings have 2 blocks in latitude and 2 in
longitude:

    #GRID
    2           lons
    2           lats
    -90.0       minimum latitude to model
    90.0        maximum latitude to model
    0.0		    longitude start to model (set to 0.0 for whole Earth)
    0.0         longitude end to model (set to 0.0 for whole Earth)

The number of CPU cores required is the product of the number of blocks in
longitude and latitude. Doubling the resolution in both latitude and longitude
requires four times as many cores since 
$`(2*nLons) \times (2*nLats) = 4\times (nLons*nLats)`$.

```math
nCPUs = (nBlocksLon) \times (nBlocksLat)
```


### Cells

In `ModSize.f90`, the following variables are defined:

```fortran
integer, parameter :: nLons = 9
integer, parameter :: nLats = 9
integer, parameter :: nAlts = 50

integer, parameter :: nBlocksMax = 1
```

The first three variables (`nLons`, `nLats` and `nAlts`) define the size of each
block. In the defaults shown above, there are 9 cells in latitude, 9 cells in
longitude and 50 altitude cells on each block.

The final variable (`nBlocksMax`) defines the maximum number of blocks you can
have on a single processor. This is deprecated and will cause issues when set to
something other than 1. It will be removed in a future release.

!!! note
    If you change any of these parameters, you will need to recompile the code.


## Horizontal Resolution

The number of *total* grid cells in the examples above will then be
$`(2\times9=)18`$ in both longitude & latitude, and 50 in altitude. Given that
there are $`360^\circ`$ in longitude and $`180^\circ`$ in latitude, the
resolution would be $`360^\circ/18 = 20.0^\circ`$ in longitude and
$`180^\circ/18 =10^\circ`$ in latitude. These settings are useful for debugging
and developing, but not for production runs.

For production runs, the resolution typically desired is $`4^\circ`$ longitude
by $`1^\circ`$ latitude. With the default values in `src/ModSize.f90`, this
requires 200 processing cores and the following values in the `UAM.in` file:

    #GRID
    10           lons
    20           lats
    -90.0        minimum latitude to model
    90.0         maximum latitude to model
    0.0		     longitude start to model (set to 0.0 for whole Earth)
    0.0          longitude end to model (set to 0.0 for whole Earth)

The longitudinal resolution ($`\Delta{\phi}`$) is set by:

```math
\Delta{\phi} = \frac{\phi_{end} - \phi_{start}}{nBlocksLon \times nCellsLon}
```

While the latitudinal resolution ($`\Delta{\theta}`$) is set by:

```math
\Delta{\theta} = \frac{\theta_{end} - \theta_{start}}{nBlocksLat \times nCellsLat}
```

---

Some recommended settings are listed below. An ideal setup would minimize the
unused cores on each node, while balancing runtime (increases with more cells
per block) and the number of files created (increases with more cores (blocks))

| Resolution <br> $`(^\circ Lon \times ^\circ Lat)`$ |   Total Cores   |  nBlocks <br> (nLons x nLats)  | nCells <br> (nLon x nLat) |
| :---  |      ----       |   :-----:   | :----: |
| 4 x 1 | 200             | 10 x 20     | 9 x 9  |
| 4 x 1 | 360 = (120 * 3) | 10 x  36    | 9 x 5  |
| 4 x 1 | 120             | 6 x 20      | 15 x 9 |
| 4 x 1 | 180             | 9 x 20      | 10 x 9 |
| 2 x 0.5 | 800           | 20 x 40     | 9 x 9  |

Experimentation may be necessary to find the parameters which work best on the
system GITM is run on. 

## Running 3D Over the Part of the Globe

GITM can be run over part of the globe - both in latitude and in longitude. It
can be run over a single polar region (by setting either the minimum or maximum
latitude to be greater (or less) than $`\pm 90^\circ`$). If this is selected,
message passing over the poles is implemented. If the pole is not selected, then
boundary conditions have to be set in `set_horizontal_bcs.f90`. By default, a
continuous gradient boundary condition is used on the densities and
temperatures, while a continuous value is used on the velocity. This is true in
both latitude and longitude. In longitude, message passing is implemented all of
the time, but the values are over-written by the boundary conditions if the
maximum and minimum longitude are not equal to each other.

## Altitudes

As defined in `src/ModEarth.f90`, each altitude block contains $`\frac{1}{3}`$
of a scale height, starting at 100 km and typically reaching up to ~500-700 km.
Increasing the number of altitude blocks will increase the altitude range, but
if this value is increased too much the model may become unstable and/or
inaccurate.

## Running in 1D

GITM can run in 1D mode, in which the call to advance_horizontal is not
completed. This means that GITM runs exactly the same way, but ignoring
all of the horizontal advection terms. You have to do two things to make
GITM run in 1D. First, in `ModSize.f90`:

    integer, parameter :: nLons = 1
    integer, parameter :: nLats = 1
    integer, parameter :: nAlts = 50

    integer, parameter :: nBlocksMax = 1

This tells the code that you only want one single latitude and longitude
location. To specify the exact location, in `UAM.in`:

    #GRID
    1           lons
    1           lats
    41.75       minimum latitude to model
    41.75       maximum latitude to model
    275.0       minimum longitude to model
    275.0       maximum longitude to model

This is pretty close to some place in Michigan. GITM will model this
exact point for as long as you specify. One thing to keep in mind with
running in 1D is that the Earth still rotates, so the spot will have a
strong day to night variation in temperature. In 3D, the winds decrease
some of the variability between day and night, but in 1D, this doesn't
happen. So, the results are not going to be perfect. But 1D is great for
debugging.
