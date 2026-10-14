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
120.0       DtAurora
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
Names are not case sensitive. The accepted strings are:

| Type      | Model             | Accepted names                          |
| :---      | :---              | :---                                    |
| Aurora    | none              | `zero`                                  |
| Aurora    | FTA               | `fta`                                   |
| Aurora    | FRE               | `fre`, `ihp`, `hpi`                     |
| Aurora    | PEM               | `pem`                                   |
| Aurora    | Ovation Prime     | `ovation`, `ovationprime`               |
| Aurora    | AMIE file         | `amie`, `amiefile`, `swmffile`          |
| Potential | none              | `zero`                                  |
| Potential | Weimer05          | `weimer`, `weimer05`                    |
| Potential | Millstone Hill    | `millstone`                             |
| Potential | Heppner-Maynard   | `hepmay`, `heppnermaynard`              |
| Potential | AMIE file         | `amie`, `amiefile`, `swmffile`          |

An unrecognized name raises an error. See the file for the most up-to-date list.

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

- `NormalizeAuroraToHP` scales the aurora to the hemispheric power, and is only
  recommended in conjunction with FRE or FTA. GITM raises a warning if it is used with
  any other auroral model (and a second one if used with Ovation).
- `#AURORATYPES` are not supported by all auroral models. Presently, only
  OVATION & AMIE files can provide other than electron diffuse aurora (mono,
  wave and ion). Other models return a negligible placeholder for these types.
- `AllowAurWODiffuse` was added for stability with OVATION-Prime; it restricts
mono/wave/ion aurora to only exist in locations which also contain electron
diffuse aurora. This can be set in [`#AURORAMODS`](../set_inputs.md#auroramods)

Internally, GITM represents Monoenergetic and Wave/broadband aurora with a
gaussian centered at the average energy. 

### Aurora Mods

The diffuse aurora can be represented by
either a Maxwellian or Kappa distribution using the `#AURORAMODS` section of
`UAM.in`:

```
#AURORAMODS
F               NormalizeAuroraToHP     (logical)
1.0             AveEFactor    (real)
F               IsKappaAurora     (logical)
3.0             AuroraKappa    (real)
F               AllowAurWODiffuse (logical)
80.0            MaxAveEAurora    (real)
```

The values shown are the defaults.

### Reference Altitude

The aurora is not mapped at an arbitrary height. The altitude set by
[`#IEREFERENCEALT`](../set_inputs.md#iereferencealt) (default 120 km) is
matched to the nearest model level, which is printed at startup (with a high
enough debug level). The auroral patterns are taken from that level, and the
hemispheric power is integrated over the cells at that level. The dynamo does
not use this altitude: its magnetic grid is placed at the lowest model level
(`Altitude_GB(1,1,1,1)` in `calc_electrodynamics.f90`).



## Potentials

The following electric field models can be used:

- Weimer05
- Millstone-Hill
- Heppner Maynard
- AMIE

## Dynamo

The low- and mid-latitude electric field is solved by GITM's own dynamo, which
is driven by the neutral winds, with the high-latitude potential above as its
boundary condition. It is controlled by
[`#DYNAMO`](../set_inputs.md#dynamo) (`UseDynamo` defaults to F, with
`DynamoHighLatBoundary` 65.0) and [`#DYNAMOSOLVER`](../set_inputs.md#dynamosolver).
See `set_inputs.md` for the full list of options and defaults.

## Required Input Files

Each electrodynamics module has different inputs. At initialization, a
verification check will be performed where GITM ensures that all the required
input data are present and that the data file covers the entire simulation time
range. If an input file ends before the requested stop time for a run, errors
will be raised.

The check for valid data is located within 
[`Electrodynamics/src/indices_subroutines.f90`](https://github.com/GITMCode/Electrodynamics/blob/main/src/indices_subroutines.f90).
Some required input file-types are listed below:


| Model    | Needs                                  |
| :---     | :---                                   |
| Weimer05 | IMF By & Bz, solar-wind V and N        |
| HepMay   | IMF By & Bz, Kp                        |
| FTA      | AU and AL                              |
| FRE      | Hemispheric power (north and south)    |
| PEM      | Hemispheric power (north and south)    |
| Ovation  | IMF By & Bz, solar-wind V (not Kp)     |

Each input is supplied by the following `UAM.in` commands:

| Input          | Commands                                                                                  |
| :---           | :---                                                                                      |
| IMF & solar wind | [`#MHD_INDICES`](../set_inputs.md#mhd_indices) (time-varying) or [`#SOLARWIND`](../set_inputs.md#solarwind) (static) |
| AU & AL        | [`#SME_INDICES`](../set_inputs.md#sme_indices) or [`#SMESINGLE`](../set_inputs.md#smesingle)                          |
| Hemispheric power | [`#NOAAHPI_INDICES`](../set_inputs.md#noaahpi_indices), [`#HPI`](../set_inputs.md#hpi), or `#SME_INDICES` with `doUseAeForHp` = T |
| Kp             | [`#KP`](../set_inputs.md#kp), a single value only (no time-varying Kp file is read)      |

HP can be derived from AE, and is not necessarily required to be in a standalone
file. This only happens when the third value of `#SME_INDICES`
(`doUseAeForHp`) is T, in which case HP = 0.102 AE + 8.953 GW. See [here](../common_inputs.md#sme-indices) for more details.

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

`#APEX` defaults to F (`UseApex = .false.`), so GITM uses a tilted, offset dipole
for the magnetic field unless `#APEX T` is given, which enables the realistic
Apex field. The tilt and offset default to the IGRF 2020 epoch, however with the use
of the [`#DIPOLE`](../set_inputs.md#dipole) option, it is possible to force zero 
offset and tilt, effectively aligning the geographic and magnetic poles.
The output files will then be in magnetic
coordinates, on a magnetic grid, rather than geographic.

A complete example file for this can be found in
`srcTests/auto_test/UAM.in.05.ElectrodynamicsMagCoords.test`, where the
following sections are what differs this test from the previous:

```
#APEX 
F     Apex is turned off (so a dipole is used; this is the default)

#DIPOLE
0.0         Magnetic Pole rotation
0.0         Magnetic pole tilt
0.0         x Dipole Center
0.0         y Dipole Center
0.0         z Dipole Center

```

## Library Notes

The Electrodynamics library has had Weimer05 fixes: the Southern Hemisphere
tilt and IMF By handling were corrected, and the potential near the pole no
longer returns zero. Update `ext/Electrodynamics` if older results look odd.
