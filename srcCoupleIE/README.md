Two way IE-GITM coupling has been fully re-established within the SWMF. This README contains information regarding changes made to IE and GITM to achieve this goal, suggestions for optimal inputs for successful runs with GITM in the SWMF, and more(?).

# Changes to IE:
Use an auroral precipitation model (MAGNIT or IMP)
If using MAGNIT, ensure #TRACEIE is true in GM
If using IMP, ensure that 2 way coupling is enabled between IM and IE
The IMP model can intake and pass a full energy spectrum distribution of precipitation to GITM while calculating any acceleration due to a parallel potential drop. This option is the default when both IMP and GITM are in the simulation. 


# Changes to UA:
The primary changes exist in this folder (srcCoupleIE). GITM was previously updated to use an external electrodynamics library that provides auroras and potentials. The code in srcCoupleIE creates a class that mimics the class in the external IE library, but uses coupling information from IE as inputs, as opposed to empirical models. The class' primary purpose is to take in and store the most recent potential and auroral values from IE, then interpolate those values to GITM's grid when called.

> This is temporary & will, eventually, be merged with GITMCode/Electrodynamics into one library.
> When the new global potential solver is added it will go there as well.

Other changes to GITM include:
- Add energy spectrum output, 3DEDF, (probably temporary)
  - Had to change ED_N_Energies in ModSources.f90 to 54 so energies align with alt dimension & can use generic 3D outputter
  - The `z` dimension is the energy bin
- Change FAC quantity passed to IE from `DivJuAltMC` (in calc_electrodynamics) to `JParAlt` (in calc_thermoelectric_current)
- Change conductances passed to IE from field-line to altitude integrated quantities, add these to outputs
- Modify default dipole tilt & rotation and change how mag coords are handled without apex
- Ensure subsolar point rotates correctly
- Add `UseSpectrumAurora` which bypasses the aveE/eFlux -> spectrum calculations, this can be used in combination with non-spectrum aurora
- Do not get_potential when initializing, we do not have values from the SWMF yet
- Change calls to `gmres` in calc_electrodynamics to use the `nItersMax` and `MaxResidual` values from input file
- Remove unused `nMagLons`/`nMagLats` from `ModElecterodynamics.f90`; the number of points is set via the resolution

# PARAM additions to CON:
Add in 2 way coupling between UA and IE (add to COUPLEORDER as well)
set #COUPLETIME to true for UA

# PARAM additions to IE:
IE has the option to manually turn off spectral aurora, which will otherwise be defaulted to true when coupling with GITM (if the IMP model is being used). ! ~~This is not implemented yet~~ :(

Currently, turning off specific types of precipitation is not implemented when coupled with GITM. 

# PARAM additions to UA:

A baseline suggestion for UA run configuration can be found in the nightly test 3 PARAM file. This contains the conditions that were used when developing and testing the coupling, and seem to produce the best results. 

- If using the standard auroral precipitation outputs from IE, turn on all 4 booleans under #AURORATYPES. The #ELECTRODYNAMICS command should be set to swmf for both the aurora and potential sources.
- If using the spectral auroral precipitation outputs from IE, the #AURORATYPES command will have no effect. The #ELECTRODYNAMICS command should be set to spectrum for aurora and swmf for potential.
  - The broadband/wave aurora will be calculated from aveE/eFlux, overriding any user settings.
- The Dyanmo High Lat boundary in #DYNAMO should be set to 88.0 or higher (only tested with 88.0). This sets the magnetic grid on which values passed back to IE are calculated, and setting it at its default value (45) will cause unrealistic coupling/physics.