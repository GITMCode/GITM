# Electrodynamics

GITM uses [an external library](https://github.com/gitmcode/electrodynamics.git)
for high-latitude electrodynamics. This library is automatically cloned into
`ext/Electrodynamics` when running `Config.pl`, though an update is not
attempted.

GITM has the ability to use external auroral and/or potential models. These are
specified independently so one could, for instance, use AMIE (file-based) aurora
and Weimer potentials. 

## Specifying Electrodynamics Drivers

The Electrodynamics models are chosen in the `#ELECTRODYNAMICS` section of 
`UAM.in`. By default the values are set to:

```
#ELECTRODYNAMICS
zero        AuroralModel
60.0        DtAurora
zero        PotentialModel
60.0        DtPotential
```

This will cause warnings to be printed if running on Earth, as we often wish to
provide high-latitude electrodynamics drivers when modeling the Earth. However,
the run will complete.

If, for example, one wishes to perform a scientific run using a commonsense
configuration, the recommended settings are:

```
#ELECTRODYNAMICS
fta         AuroralModel
60.0        DtAurora
weimer05    PotentialModel
60.0        DtPotential
```

This will use FTA[^1] and Weimer05[^2] for the auroral and potential models,
respectively. 

[^1]: Wu, C., Ridley, A. J., DeJong, A. D., & Paxton, L. J. (2021). FTA: A Feature Tracking Empirical Model Of Auroral Precipitation. Space Weather, 19, e2020SW002629. <https://doi.org/10.1029/2020SW002629>.

[^2]: Weimer, D. R. (2005), Improved ionospheric electrodynamic models and application to calculating Joule heating rates, J. Geophys. Res., 110, A05306, <https://doi.org/10.1029/2004JA010884>.


The aurora and electric field model names are parsed in 
[`Electrodynamics/src/interpret_names.f90`](https://github.com/GITMCode/Electrodynamics/blob/main/src/interpret_names.f90).
See this file for the most up-to-date list of available modules and the acceptable names.

## Aurora

The following Auroral models are available:

- FTA
- FRE
- PEM
- OVATION
- AMIE

From the auroral module, GITM expects to receive Average Energy and Energy Flux,
for all of the selected auroral types. At the moment these all must be from the
same module, so one cannot use FTA for diffuse electron precipitation and AMIE
for monoenergetic electron precipitation.

### Aurora Types

Auroral types are specified in the `#AURORATYPES` section of `UAM.in`, and only
electron diffuse aurora are included by default:

```
#AURORATYPES
T         UseDiffuseAurora (logical)
F         UseMonoAurora (logical)
F         UseWaveAurora (logical)
F         UseIonAurora (logical)
```

Some notes on the different auroral types:

- `NormalizeAuroraToHP` is only recommended to be used in conjunction with FRE,
- `#AURORATYPES` are not supported by all auroral models. Presently, only
OVATION & MAGNIT (AMIE) can provide other than electron diffuse aurora.
- `AllowAurWODiffuse` was added for stability with OVATION-Prime; it restricts
mono/wave/ion aurora to only exist in locations which also contain electron
diffuse aurora. This can be set in [`#AURORAMODS`](../set_inputs.md#auroramods)

Internally, GITM represents Monoenergetic and Wave/broadband aurora with a
gaussian centered at the average energy. 

### Aurora Mods

The diffuse aurora can be represnted by
either a Maxwellian or Kappa distribution using the `#AURORAMODS` section of
`UAM.in`:

```
#AURORAMODS
F               NormalizeAuroraToHP     (logical)
1.0             AveEFactor    (real)
F               IsKappaAurora     (logical)
1.0             AuroraKappa    (real)
F               AllowAurWODiffuse (logical)
50.0            MaxAveEAurora    (real)
```



## Potentials

The following electric field models can be used:

- Weimer05
- Millstone-Hill
- Heppner Maynard
- AMIE

## Required Input Files

Each electrodynamics module has different inputs. At initialization, a
verification check will be performed where GITM ensures that all the required
input data are present and that the data file covers the entire simulation time
range. If an input file ends before the requested stop time for a run, errors
will be raised.

The check for valid data is located within 
[`Electrodynamics/src/indices_subroutines.f90`](https://github.com/GITMCode/Electrodynamics/blob/main/src/indices_subroutines.f90).
Some required input file-types are listed below:


| Model   |   IMF   |  AE   |  HP   | Kp    |
| :---    | :---:   | :---: | :---: | :---: | 
| Weimer  | Yes     | No    | No    | No    |
| FRE     | Yes     | No    | Yes   | No    |
| HepMay  | Yes     | No    | No    | Yes   |
| FTA     | No      | Yes   | No    | No    |
| PEM     | No      | No    | Yes   | No    |
| Ovation | Yes     | No    | No    | Yes   |

HP can be derived from AE, and is not necessarily required to be in a standalone
file. See [here](../common_inputs.md#sme-indices) for more details.

## File-based Electrodynamics

AMIE (Assimilative Mapping of Ionospheric Electrodynamics) is the name chosen
for the type of files which can be interpreted by Electrodynamics. A number of
[Python routines](https://github.com/GITMCode/Electrodynamics/tree/main/srcPython)
can be found in Electrodynamics which may be useful when generating these inputs.

AMIE files can be used, for instance, to input custom auroral patterns into
GITM. Once could generate an AMIE file with nominal diffuse electron aurora and
several intense monoenergetic beams at any number of locations.

## Running Electrodynamics Only

By using the [`#STATISTICALMODELSONLY`](../set_inputs.md#statisticalmodelsonly)
option in `UAM.in`, it is possible to run any configuration of Electrodynamics
models without GITM's physics, making the runs faster. 

By setting the desired output type to [`2DGEL`](../outputs.md#2dgel), and an
appropriate Dt for `#STATISTICALMODELSONLY` and `#OUTPUT`, GITM will read in the
necessary input files and output precipitation & potential patterns using the
specified electrodynamics modules. An example of this is located in
`srcTests/auto_test/UAM.in.04.ElectrodynamicsGeoCoords.test`.

Additionally, one can output data on a magnetic grid instead of geographic,
which is often desired when plotting outputs from the electrodynamics modules.
To do this, one must manually control the magnetic field configuration through
the use of `#APEX` and `#DIPOLE` in `UAM.in`.

By setting `#APEX` to F, GITM will use a tilted, offset dipole for the magnetic
field. The tilt and offset are normally set automatically, however with the use
of the `#DIPOLE` option, it is possible to force zero offset and tilt,
effectively aligning the geographic and magnetic poles. The output files will
then be in magnetic coordinates, on a magnetic grid, rather than geographic.

A complete example file for this can be found in
`srcTests/auto_test/UAM.in.05.ElectrodynamicsMagCoords.test`, where the
following sections are what differs this test from the previous:

```
#APEX 
F     Apex is turned off (so a dipole is used)

#DIPOLE
0.0         Magnetic Pole rotation
0.0         Magnetic pole tilt
0.0         x Dipole Center
0.0         y Dipole Center
0.0         z Dipole Center

```

