# Overview

The Global Ionosphere Thermosphere Model (GITM) is a 3D model of the
upper atmosphere. It runs for Earth, Mars and Titan. A version is being
worked on for Jupiter. GITM solves for the coupled continuity, momentum
and energy equations of the neutrals and ions. For the ions, the time
rate of change of the velocity is ignored, so the steady-state ion flow
velocity is solved for. The ion temperature is a mixture of the electron
and neutral temperature.

The neutrals are solved for using the Navier Stokes Equations. The
continuity equation is solved for for each major species. One of the
problems with GITM that needs to be rectified is that there are no real
tracer species, so a species is either solved for completely or is not
at all. These species can still be included in the chemistry
calculation. There is only one horizontal velocity that is computed,
while there are vertical velocities for each of the major species. A
bulk vertical velocity is calculated as a mass weighted average. The
temperature is a bulk temperature.

## Source Terms

Chemistry is the only real source term for the continuity equation.
Typically, diffusion is added in the continuity equation to allow for
eddy diffusion, but this is not the case in GITM. What happens is that
the vertical velocities are solved for, then a friction term is applied
to that the velocities stay very close together in the eddy diffusion
part of the code. This way, the velocities can't differ too much from
each other. Diffusion is not needed, then.

For the horizontal momentum equation, there are the following sources:
(1) ion drag; (2) viscosity; and (3) gravity wave acceleration. For the
vertical velocity, the source terms are ion drag and friction between
the different neutral species.

For the neutral temperature, the following source terms are included:
(1) radiative cooling; (2) EUV heating; (3) auroral heating; (4) Joule
heating; (5) conduction; and (6) chemical heating. The biggest pain for
the temperature equation is the use of a normalized temperature. This
means that the `temperature` variable in GITM does not contain the
actual temperature, it contains the temperature multiplies by
Boltzmann's Constant divided by the mean mass. This turns out to be a
factor that is very similar to the specific heat, or roughly or order
1000. In order to get the actual temperature, the variable has to be
multiplied by `temp_unit`.

## Ghost Cells

GITM is a parallel code. It uses a 2D domain decomposition, with the
altitude domain being the only thing that is not broken up. Blocks of
latitude and longitude are used. These blocks are then distributed among
different processors. In order to communicate between the processors,
ghostcells are used. These are cells that essentially overlap with the
neighboring block. MPI (message passing interface) is then used to move
information from one block to another, filling in the ghostcells. The
code then loops from 1-N, where the flux is calculated at the boundaries
from the 0-1 boundary to the N-N+1 boundary. A second order scheme is
used to calculate the fluxes, along with a flux limiter. Therefore, two
ghost cells are needed.

In the vertical direction, ghost cells are also used to set boundary
conditions. The values in these cells are used to calculate the fluxes,
just as described above. Different types of boundary conditions
(constant values, constant fluxes, constant gradients, floating, zero
fluxes, etc) can be set by carefully choosing the right values in the
ghost cells.
