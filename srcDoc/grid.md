# Setting the Grid

Setting the grid resolution in GITM is not very complicated, but it does
involve some thought. There are a few variables that control this. In
`ModSize.f90`, the following variables are defined:

    integer, parameter :: nLons = 9
    integer, parameter :: nLats = 9
    integer, parameter :: nAlts = 50

    integer, parameter :: nBlocksMax = 1

The first three variables (`nLons`, `nLats` and `nAlts`) define the size of a
single block. In the example above, there are 9 cells in latitude, 9 cells in
longitude and 50 cells in altitude. The latitude and longitude resolution in
GITM is defined by the cell numbers specified here and the block numbers
discussed in the following section. ~~The size of the altitude cells are
measured in scale heights instead of kilometers, as certain regions require
higher resolutions than others based on the dominating chemical and dynamical
processes~~. As defined in `src/ModEarth.,f90`, each altitude cell contains
$`\frac{1}{3}`$ of a scale height, starting at ~~80~~ 100 km and typically
reaching up to 500 km. Increasing the number of altitude blocks will increase
the altitude range, but if this value is increased too much the model becomes
unstable. ~~This altitude limit makes it problematic to model the equatorial
region during storms, where increased vertical plasma transport requires the
model consider altitudes of 2000 km and higher.~~

The final variable (`nBlocksMax`) defines the maximum number of blocks you can
have on a single processor. Most people run with one single block per processor,
as this is faster, and is necessary when running with the Dynamo option on. This
is a necessity because the Dynamo routine does not correctly transfer the
geomagnetically oriented data into the geographic coordinates used in the ghost
cells.

!!! note
    If you change any of these parameters you will need to recompile the code
    so that a new executable using the newly defined variables can be produced.

## Running 3D Over the Whole Globe

Once the number of cells is defined, then the number of blocks in
latitude and longitude need to be defined. This is done in the `UAM.in`
file, as shown in section [\[uam.sec\]](#uam.sec){reference-type="ref"
reference="uam.sec"}. For example, the initial settings have 8 blocks in
latitude and 8 in longitude:

    #GRID
    8           lons
    8           lats
    -90.0       minimum latitude to model
    90.0        maximum latitude to model
    0.0         minimum longitude to model
    0.0         maximum longitude to model

The number of cells in the simulation domain will then be 72 in
longitude, 72 in latitude, and 50 in altitude. Given that there are
$`360\deg`$ in longitude and $`180\deg`$ in latitude, the resolution would
be $`360\deg/72 = 5.0\deg`$ and $`180\deg/72 = 2.5\deg`$ in latitude. If one
block were put on each processor, 64 processors would be required.

If one desired to run without the Dynamo turned on, the problem could
also be run with multiple blocks per node. This grid fits quite nicely
on either four- or eight-core processors. However, on 12-core processors
this would not work very well at all. We have started to run simulations
of $`1\deg`$ in latitude and $`5\deg`$ in longitude, which can fit nicely on
a 12-core processor machine. For example, in `ModSize.f90`:

    integer, parameter :: nLons = 9
    integer, parameter :: nLats = 15

and in `UAM.in`:

    #GRID
    8       lons
    12      lats

This can then run on 96 cores, which is nicely divisible by 12.
Essentially, an infinite combination of cells per block and number of
blocks can be utilized. Typically, the number of blocks in latitude and
longitude are even numbers.

## Running 3D Over the Part of the Globe

GITM can be run over part of the globe - both in latitude and in
longitude. It can be run over a single polar region (by setting either
the minimum or maximum latitude to be greater (or less) than
$`\pm 90\deg`$). If this is selected, message passing over the poles is
implemented. If the pole is not selected, then boundary conditions have
to be set in `set_horizontal_bcs.f90`. By default, a continuous gradient
boundary condition is used on the densities and temperatures, while a
continuous value is used on the velocity. This is true in both latitude
and longitude. In longitude, message passing is implemented all of the
time, but the values are over-written by the boundary conditions if the
maximum and minimum longitude are not equal to each other.

The longitudinal resolution ($`\Delta{\phi}`$) is set by:

```math
\Delta{\phi} = \frac{\phi_{end} - \phi_{start}}{nBlocksLon \times nCellsLon}
```

while, the latitudinal resolution ($\Delta{\theta}$) is set by:

```math
\Delta{\theta} = \frac{\theta_{end} - \theta_{start}}{nBlocksLat \times nCellsLat}
```

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
happen. So, the results are going to be perfect. But, 1D is great for
debugging.
