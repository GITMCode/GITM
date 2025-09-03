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
    If you change any of these parameters, you have to recompile the code.


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

Some settings are listed below.  Machines like NASA's Pleiades have a certain number of cores (processors) per node, like Ivy has 20 cores/node. An ideal setup would minimize the
unused cores on each node, while balancing runtime (increases with more cells
per block) and the number of files created (increases with more cores (blocks)). What this practically means is that it is a good idea to have the total number of processors used being a multiple of the processors/node.  

| Resolution <br> $`(^\circ Lon \times ^\circ Lat)`$ |   Total Cores   |  nBlocks <br> (nLons x nLats)  | nCells <br> (nLon x nLat) |
| :---  |      ----       |   :-----:   | :----: |
| 4 x 1 | 200             | 10 x 20     | 9 x 9  |
| 4 x 1 | 360             | 10 x 36     | 9 x 5  |
| 4 x 1 | 120             | 6 x 20      | 15 x 9 |
| 4 x 1 | 180             | 9 x 20      | 10 x 9 |
| 4 x 1 | 10              | 2 x 5       | 45 x 36 |
| 2 x 0.5 | 800           | 20 x 40     | 9 x 9  |

Experimentation may be necessary to find the parameters which work best on the
system GITM is run on. You can always increase nCells in order to run on fewer processors - GITM will just take longer to run, but should be fine. For example, the second to last line shows that GITM can run at production levels on 10 processors, which relatively high-end laptops have currently and many modern desktop machines have.

## Running 3D Over Part of the Globe

GITM can be run over part of the globe - both in latitude and in longitude. It
can be run over a single polar region by setting either the minimum or maximum
latitude to be greater (or less) than $`\pm 90^\circ`$). If this is selected,
message passing over the poles is implemented. For example, this simulates the entire northern polar region at $`4^\circ`$ longitude and $`1^\circ`$ latitude resolution:
    #GRID
    8            lons
    4            lats
    45.0         minimum latitude to model
    90.0         maximum latitude to model
    0.0		     longitude start to model (set to 0.0 for whole Earth)
    0.0          longitude end to model (set to 0.0 for whole Earth)

If only part of the globe is selected, then boundary conditions have to be set in `set_horizontal_bcs.f90`. By default, a
continuous gradient boundary condition is used on the densities and
temperatures, while a continuous value is used on the velocity. This is true in
both latitude and longitude.

Running in a small region does not capture the global-scale dynamics, since the EUV heating on the dayside would not be specified correctly in the regions that are not modeled.  So, things like the global-scale neutral wind pattern will not be correct. This means that if the physics that you are desiring to simulate depends on things that are more global-scale, then it is not a good idea to run over an isolated region.  But, if the physics is much more localized and the time-scales are significantly shorter than a day, this is a fine feature to use. Example of good and bad cases would include:

- Good case: Exploring how an ionospheric travelling convection vortex interacts with the local thermosphere to drive strong vertical winds.
-Bad case: Exploring how to meridional neutral winds drive field-aligned flows over Michigan.

In order to mitigate the effects of this (but not completely eliminate the effects), we have implemented a feature that allows you to use global-scale simulation results as a horizontal boundary condition for a regional simulation.  In addition, GITM initializes the states with the global-scale simulation, so that you don't have to run the regional-scale simulation for a long period to get rid of the initial condition. In order to do this, you need to do the following:

- Run one simulation over the whole globe using some nominal resolution, outputting 3DALL files at a time-scale that is appropriate for capturing the global-scale dynamics (say every 15 minutes if there is nothing of geophysical interest occuring).  Remember that it is always a good idea to run 24 hours before an event in order to get rid of the initial condition!
- Post process the 3DALL files.
- Move the UA/data directory into some place like UA/GlobalScaleOutputs.
- Make a new UA/data directory.
- Modify the UAM.in file so that the grid captures the region of interest, and the start/end times align with the time of interest. Also, set the number of blocks so that you have the resolution you desire.
- Add the following line into the UAM.in file:
    #GITMBCS
    T
    UA/GlobalScaleOutputs

Then, when you run GITM, it will simulate only the region of interest.

## Altitudes

As defined in `src/ModEarth.f90`, the altitude spacing in GITM is 0.3  
scale heights. In the UAM.in file, the starting altitude is typically set to 100 km.  Given that in ModSize.f90, the number of altitudes is 50 (by default), GITM simulates 15 scale heights, typically reaching up to ~500-700 km.
Increasing the number of altitudes in ModSize.f90 will increase the altitude range, but
if this value is increased too much the model may become unstable and/or inaccurate.  It is also possible to reduce the starting altitude below 100 km, but decreasing it below the mesopause is not a great idea, since the physics of the mesopause (i.e., CO2 cooling) is not included in the current version of GITM.  Therefore, moving the lower boundary below about 95 km is not recommended.

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

This is some place in Michigan. GITM will model this
exact point for as long as you specify. One thing to keep in mind with
running in 1D is that the Earth still rotates, so the spot will have a
strong day to night variation in temperature and will not capture any horizontal dynamics at all, such as the horizontal neutral winds (which are primarily driven by pressure gradients). In 3D, the winds decrease
some of the variability between day and night, but in 1D, this doesn't
happen. So, the results are not going to be (even close to) perfect. But 1D is great for
debugging, since it is incredibly fast and many days can be run within a few minutes.
