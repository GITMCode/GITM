# Outputs {#outputs.sec}

Now that you have managed to successfully complete a GITM run you've found yourself with a bunch of output files.
Output files are written to `UA/data` (unless [`#DIRECTORIES`](set_inputs.md#directories) says otherwise), along with a log file, `run_information.txt` & `output_schema.csv`.

Skip to [the output type section](#output-types-and-their-variables) if you are just looking for what to output.

## File formats

The format of these files depends on which [output backend\(`#OUTPUTBACKEND`\)](outputs/output_backends.md) was used during the run:

| Backend | Files for each output type, at each output time | Before you can read them |
|---|---|---|
| `legacy` (default) | one `.header` + one `.b####` block file per processor | [post-process](postprocessing.md) them into one `.bin` |
| `mpiio` | one `.header` + one raw `.bin` | [post-process](postprocessing.md) it into a standard `.bin` |
| `netcdf` | one `.nc` | nothing, they are complete as written |

Satellite files ([`#SATELLITES`](set_inputs.md#satellites)) are written in the legacy format for all backends.

For a given output type all three backends write the same variables, because each type is defined once in the code and every backend reads that definition.
[Adding new outputs](outputs/modifying_outputs.md) describes how that works.

## Units

All of the GITM outputs are in MKS units. 

!!! important "Horizontal coordinates: radians in `.bin`, degrees in netcdf"
    The one exception to "MKS units" is geographic longitude and latitude,
    and it differs **by backend**:

    | Backend | `Longitude` / `Latitude` |
    |---|---|
    | `legacy`, `mpiio` (`.bin`) | **radians** |
    | `netcdf` | **degrees** (`degrees_east` / `degrees_north`) |

    The variable tables below, and `output_schema.csv`, show the netcdf units
    for `Longitude`/`Latitude`.
    Only the magnetic coordinates `mlon`/`mlat` are in degrees on every backend.
    When reading a `.bin` file, convert with `numpy.degrees` yourself.

## The other files in `UA/data`

The log file provides useful information about the run, such as whether a restart was performed, which physical processes were used, and some variables which represent the state of the modeled system and its drivers.
This file can be very useful when sharing runs with other users, when revisiting an old run, or merely ensuring that GITM performed as expected.
Some example log files can be found in the `srcTests/auto_test/ref_solns/` directory.

The `run_information.txt` file has every setting used by GITM during a run.
Its format is similar to `UAM.in`, and every variable is written in the same order in which it is parsed by [`set_inputs.f90'](set_inputs.md).

The `output_schema.csv` file lists every output type and its variables, with units and a description, in the order they appear in the files.
`legacy` and `mpiio` files carry only the variable names, so this is where to look up their units.

## More information/FAQ

Here is where to look for more information or if you have questions:

- What is available? [see below](#output-types-and-their-variables)
- Output backends
  - [What is available & what should I use](outputs/output_backends.md)
  - [Details on how the `netcdf` backend handles restarts & appending](outputs/output_backends.md#netcdf)
- [Post-processing `legacy` and `mpiio` output](postprocessing.md)
- Modifying Outputs
  - [Should I create a new output type or add to an existing one?](outputs/modifying_outputs.md)
  - [Adding a temporary output inside a physics routine (USR outputs)](outputs/modifying_outputs.md#add-a-variable-to-a-usr-output)
  - [Adding a new variable to existing output types](outputs/modifying_outputs.md#modify-an-existing-output-type)
  - [I want to create a new output type](outputs/modifying_outputs.md#create-a-new-output-type)

## Output types and their variables

Request any of these in [`#SAVEPLOT`](set_inputs.md#saveplot).
`3DALL` is the usual choice.

| Type | Contents |
|---|---|
| `3DALL` | Neutral and ion state: densities, temperatures and velocities |
| `3DNEU` | The neutral half of `3DALL` |
| `3DION` | The ion half of `3DALL`, plus electric and magnetic fields, potential, currents and the ion-neutral collision frequency |
| `3DLST` | A shorter `3DALL`: mass density, the major species, winds, ion drifts and neutral temperature |
| `3DTHM` | Heating and cooling rates, and auroral ionization rates |
| `3DCHM` | Chemical heating, reaction by reaction |
| `3DMOH` | Terms of the horizontal momentum equation |
| `3DMOV` | Terms of the vertical momentum equation, species by species |
| `3DMAG` | Magnetic coordinates and the magnetic field |
| `3DEMI` | Emission at seven wavelengths |
| `3DHME` | `3DALL` plus heating and fields, written only inside the [`#SAVEHIMEPLOT`](set_inputs.md#savehimeplot) region |
| `3DGLO` | GLOW emissions. **Stops the run**, see the notes below |
| `3DUSR` | [User-defined variables](outputs/modifying_outputs.md#add-a-variable-to-a-usr-output) |
| `2DGEL` | Electrodynamics on the geographic grid: potential, conductances and auroral precipitation |
| `2DMEL` | Electrodynamics on the magnetic grid: conductances and the terms of the dynamo solver |
| `2DTEC` | Total electron content |
| `2DANC` | Total electron content and height-integrated heating and cooling |
| `2DHME` | Total electron content, written only inside the `#SAVEHIMEPLOT` region |
| `2DUSR` | User-defined variables |
| `1DALL` | The `3DALL` variables |
| `1DNEW` | `1DALL` plus local time, solar zenith angle and mixing ratios |
| `1DTHM` | Heating and cooling rates, and production and loss rates for each species |
| `1DCHM` | Chemical heating, reaction by reaction |
| `1DEMI` | Emission at seven wavelengths |
| `1DGLO` | GLOW emissions. **Stops the run**, see the notes below |
| `1DUSR` | User-defined variables |
| `0DALL` | The `3DALL` variables plus mixing ratios and a heat balance |
| `0DUSR` | User-defined variables |

The `1D` types are altitude profiles and the `0D` types single points.
A 1-D run (`nLons = nLats = 1`) writes a requested `3D` type as its `1D`
counterpart. [`#SATELLITES`](common_inputs.md#sat_aux.sec) writes `1DALL` along
each satellite track (`0DALL` when RCMR is running), or `1DUSR`/`0DUSR` if
asked. The `netcdf` backend does not write any of them; see
[Output backends](outputs/output_backends.md).

Notes:

- **Rates are zero at the first output.** Per-step diagnostics that GITM
  divides by the time step (3DTHM `Conduction`, `Molecular Conduction` and
  the chemical heating rate, 3DCHM, 3DMOH `Visc_*`, 2DANC `Chem_int`, and
  the probe types 0DALL, 1DTHM and 1DCHM) are rates. They are written as 0
  at the t = 0 output, where the time step is 0.
- **`3DGLO` and `1DGLO` stop the run.** GLOW computes emissions when `#GLOW`
  is on (`get_glow.f90` fills `vEmissionRate`), but the two output types were
  never wired to that array: their fills put placeholder zeros,
  indistinguishable from real measurements, so GITM stops if either is
  requested in `#SAVEPLOT`. They remain valid type codes so the fills can be
  pointed at `vEmissionRate` later; they produce nothing today and have no
  table below.
- **`2DMEL` is on a magnetic grid**, so its horizontal coordinates are
  `mlon`/`mlat` (in degrees), not `Longitude`/`Latitude`. Its
  `GeoLat`/`GeoLon` variables are the geographic position of each magnetic
  grid point (radians in `.bin`, degrees in netcdf).
- **USR types** list whatever the user code registered. The tables below show
  what GITM registers as shipped.

The tables below give each type's variables in the order they appear in the
file. `Variable` is the name in legacy/mpiio `.header` files; netcdf variable
names are the same after sanitization, with the original kept in each
variable's `gitm_name` attribute. The tables are generated from the
`output_schema.csv` that GITM itself writes, so they are updated
[with a script](outputs/modifying_outputs.md#update-the-documentation), not by
hand.

<!-- BEGIN GENERATED TABLES: edit with srcPython/generate_output_docs.py, not by hand -->

### 3DALL

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Latitude | degrees_north | Geographic latitude |
| Altitude | m | Altitude above surface |
| Rho | kg/m3 | Total neutral mass density |
| [O(3P)] | /m3 | Number density of O(3P) |
| [O2] | /m3 | Number density of O2 |
| [N2] | /m3 | Number density of N2 |
| [N(4S)] | /m3 | Number density of N(4S) |
| [NO] | /m3 | Number density of NO |
| [He] | /m3 | Number density of He |
| [CO2] | /m3 | Number density of CO2 |
| [N(2D)] | /m3 | Number density of N(2D) |
| [N(2P)] | /m3 | Number density of N(2P) |
| [H] | /m3 | Number density of H |
| [O(1D)] | /m3 | Number density of O(1D) |
| Temperature | K | Neutral temperature |
| Vn (east) | m/s | Eastward neutral wind |
| Vn (north) | m/s | Northward neutral wind |
| Vn (up) | m/s | Vertical neutral wind |
| Vn (up,O(3P)) | m/s | Vertical neutral wind for O(3P) |
| Vn (up,O2) | m/s | Vertical neutral wind for O2 |
| Vn (up,N2) | m/s | Vertical neutral wind for N2 |
| Vn (up,N(4S)) | m/s | Vertical neutral wind for N(4S) |
| Vn (up,NO) | m/s | Vertical neutral wind for NO |
| Vn (up,He) | m/s | Vertical neutral wind for He |
| Vn (up,CO2) | m/s | Vertical neutral wind for CO2 |
| [O(4S)+] | /m3 | Number density of O(4S)+ |
| [He+] | /m3 | Number density of He+ |
| [NO+] | /m3 | Number density of NO+ |
| [O2+] | /m3 | Number density of O2+ |
| [N2+] | /m3 | Number density of N2+ |
| [N+] | /m3 | Number density of N+ |
| [O(2D)+] | /m3 | Number density of O(2D)+ |
| [O(2P)+] | /m3 | Number density of O(2P)+ |
| [H+] | /m3 | Number density of H+ |
| [e-] | /m3 | Number density of e- |
| eTemperature | K | Electron temperature |
| iTemperature | K | Ion temperature |
| Vi (east) | m/s | Eastward ion drift |
| Vi (north) | m/s | Northward ion drift |
| Vi (up) | m/s | Vertical ion drift |

### 3DCHM

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Latitude | degrees_north | Geographic latitude |
| Altitude | m | Altitude above surface |
| N2+ + e | eV/m3/s | N2+ + e recombination heating rate |
| O2+ + e | eV/m3/s | O2+ + e recombination heating rate |
| N2+ + O | eV/m3/s | N2+ + O reaction heating rate |
| NO+ + e | eV/m3/s | NO+ + e recombination heating rate |
| N+ + O2 | eV/m3/s | N+ + O2 reaction heating rate |
| NO + N | eV/m3/s | NO + N reaction heating rate |
| O+ + O2 | eV/m3/s | O+ + O2 reaction heating rate |
| N + O2 | eV/m3/s | N + O2 reaction heating rate |
| O2+ + N | eV/m3/s | O2+ + N reaction heating rate |
| O2+ + NO | eV/m3/s | O2+ + NO reaction heating rate |
| O2+ + N2 | eV/m3/s | O2+ + N2 reaction heating rate |
| N2+ + O2 | eV/m3/s | N2+ + O2 reaction heating rate |
| N+ + O | eV/m3/s | N+ + O reaction heating rate |
| O+ + N2 | eV/m3/s | O+ + N2 reaction heating rate |
| O(1D) + N2 | eV/m3/s | O(1D) + N2 quench heating rate |
| O(1D) + O2 | eV/m3/s | O(1D) + O2 quench heating rate |
| O(1D) + O | eV/m3/s | O(1D) + O quench heating rate |
| O(1D) + e | eV/m3/s | O(1D) + e quench heating rate |
| N(2D) + O2 | eV/m3/s | N(2D) + O2 reaction heating rate |
| O+(2D) + e | eV/m3/s | O+(2D) + e recombination heating rate |
| N(2D) + O | eV/m3/s | N(2D) + O quench heating rate |
| N(2D) + e | eV/m3/s | N(2D) + e quench heating rate |
| O+(2D) + N2 | eV/m3/s | O+(2D) + N2 reaction heating rate |
| O+(2P) + e | eV/m3/s | O+(2P) + e recombination heating rate |
| O+(2P) + O | eV/m3/s | O+(2P) + O reaction heating rate |
| O+(2P) + N2 | eV/m3/s | O+(2P) + N2 reaction heating rate |
| Chemical Heating | eV/m3/s | Volumetric chemical heating rate |

### 3DEMI

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Latitude | degrees_north | Geographic latitude |
| Altitude | m | Altitude above surface |
| Emission 2470A |  | Emission at 2470A |
| Emission 7320A |  | Emission at 7320A |
| Emission 3726A |  | Emission at 3726A |
| Emission 5200A |  | Emission at 5200A |
| Emission 10400A |  | Emission at 10400A |
| Emission 6300A |  | Emission at 6300A |
| Emission 6364A |  | Emission at 6364A |

### 3DHME

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Latitude | degrees_north | Geographic latitude |
| Altitude | m | Altitude above surface |
| Rho | kg/m3 | Total neutral mass density |
| [O(3P)] | /m3 | Number density of O(3P) |
| [O2] | /m3 | Number density of O2 |
| [N2] | /m3 | Number density of N2 |
| [N(4S)] | /m3 | Number density of N(4S) |
| [NO] | /m3 | Number density of NO |
| [He] | /m3 | Number density of He |
| [CO2] | /m3 | Number density of CO2 |
| [N(2D)] | /m3 | Number density of N(2D) |
| [N(2P)] | /m3 | Number density of N(2P) |
| [H] | /m3 | Number density of H |
| [O(1D)] | /m3 | Number density of O(1D) |
| Temperature | K | Neutral temperature |
| Vn (east) | m/s | Eastward neutral wind |
| Vn (north) | m/s | Northward neutral wind |
| Vn (up) | m/s | Vertical neutral wind |
| Vn (up,O(3P)) | m/s | Vertical neutral wind for O(3P) |
| Vn (up,O2) | m/s | Vertical neutral wind for O2 |
| Vn (up,N2) | m/s | Vertical neutral wind for N2 |
| Vn (up,N(4S)) | m/s | Vertical neutral wind for N(4S) |
| Vn (up,NO) | m/s | Vertical neutral wind for NO |
| Vn (up,He) | m/s | Vertical neutral wind for He |
| Vn (up,CO2) | m/s | Vertical neutral wind for CO2 |
| [O(4S)+] | /m3 | Number density of O(4S)+ |
| [He+] | /m3 | Number density of He+ |
| [NO+] | /m3 | Number density of NO+ |
| [O2+] | /m3 | Number density of O2+ |
| [N2+] | /m3 | Number density of N2+ |
| [N+] | /m3 | Number density of N+ |
| [O(2D)+] | /m3 | Number density of O(2D)+ |
| [O(2P)+] | /m3 | Number density of O(2P)+ |
| [H+] | /m3 | Number density of H+ |
| [e-] | /m3 | Number density of e- |
| eTemperature | K | Electron temperature |
| iTemperature | K | Ion temperature |
| Vi (east) | m/s | Eastward ion drift |
| Vi (north) | m/s | Northward ion drift |
| Vi (up) | m/s | Vertical ion drift |
| PhotoElectron Heating | K | Photoelectron heating |
| Joule Heating | K | Joule heating |
| Specific Heat | J/kg/K | Specific heat Cp |
| Magnetic Latitude | deg | Magnetic latitude |
| Magnetic Longitude | deg | Magnetic longitude |
| B.F. East | T | Eastward magnetic field |
| B.F. North | T | Northward magnetic field |
| B.F. Vertical | T | Vertical magnetic field |
| B.F. Magnitude | T | Magnetic field magnitude |
| Potential | V | Electric potential |
| PotentialY | V | Potential Y component |
| E.F. East | V/m | Eastward electric field |
| E.F. North | V/m | Northward electric field |
| E.F. Vertical | V/m | Vertical electric field |

### 3DION

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Latitude | degrees_north | Geographic latitude |
| Altitude | m | Altitude above surface |
| [O(4S)+] | /m3 | Number density of O(4S)+ |
| [He+] | /m3 | Number density of He+ |
| [NO+] | /m3 | Number density of NO+ |
| [O2+] | /m3 | Number density of O2+ |
| [N2+] | /m3 | Number density of N2+ |
| [N+] | /m3 | Number density of N+ |
| [O(2D)+] | /m3 | Number density of O(2D)+ |
| [O(2P)+] | /m3 | Number density of O(2P)+ |
| [H+] | /m3 | Number density of H+ |
| [e-] | /m3 | Number density of e- |
| eTemperature | K | Electron temperature |
| iTemperature | K | Ion temperature |
| Vi (east) | m/s | Eastward ion drift |
| Vi (north) | m/s | Northward ion drift |
| Vi (up) | m/s | Vertical ion drift |
| Ed1 | V/m | Electric field component 1 (magnetic east) |
| Ed2 | V/m | Electric field component 2 (magnetic north) |
| Je1 | A/m2 | Current density component 1 (magnetic east) |
| Je2 | A/m2 | Current density component 2 (magnetic north) |
| Magnetic Latitude | deg | Magnetic latitude |
| Magnetic Longitude | deg | Magnetic longitude |
| B.F. East | T | Eastward magnetic field |
| B.F. North | T | Northward magnetic field |
| B.F. Vertical | T | Vertical magnetic field |
| B.F. Magnitude | T | Magnetic field magnitude |
| Potential | V | Electric potential |
| E.F. East | V/m | Eastward electric field |
| E.F. North | V/m | Northward electric field |
| E.F. Vertical | V/m | Vertical electric field |
| E.F. Magnitude | V/m | Electric field magnitude |
| IN Collision Freq | Hz | Ion-neutral collision frequency |
| IonPressGrad (east) | m/s2 | Ion Pressure gradient acceleration eastward |
| IonPressGrad (north) | m/s2 | Ion Pressure gradient acceleration northward |
| IonPressGrad (up) | m/s2 | Ion Pressure gradient acceleration upward |

### 3DLST

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Latitude | degrees_north | Geographic latitude |
| Altitude | m | Altitude above surface |
| Rho | kg/m3 | Total neutral mass density |
| [O(3P)] | /m3 | Number density of O(3P) |
| [O2] | /m3 | Number density of O2 |
| [N2] | /m3 | Number density of N2 |
| [NO] | /m3 | Number density of NO |
| Vn (east) | m/s | Eastward neutral wind |
| Vn (north) | m/s | Northward neutral wind |
| Vn (up) | m/s | Vertical neutral wind |
| [O(4S)+] | /m3 | Number density of O(4S)+ |
| [NO+] | /m3 | Number density of NO+ |
| [O2+] | /m3 | Number density of O2+ |
| [N2+] | /m3 | Number density of N2+ |
| [e-] | /m3 | Number density of e- |
| Vi (east) | m/s | Eastward ion drift |
| Vi (north) | m/s | Northward ion drift |
| Vi (up) | m/s | Vertical ion drift |
| Temperature | K | Neutral temperature |

### 3DMAG

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Latitude | degrees_north | Geographic latitude |
| Altitude | m | Altitude above surface |
| Magnetic Latitude | deg | Magnetic latitude |
| Magnetic Longitude | deg | Magnetic longitude |
| B.F. East | T | Eastward magnetic field |
| B.F. North | T | Northward magnetic field |
| B.F. Vertical | T | Vertical magnetic field |
| B.F. Magnitude | T | Magnetic field magnitude |

### 3DMOH

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Latitude | degrees_north | Geographic latitude |
| Altitude | m | Altitude above surface |
| Rho | kg/m3 | Total neutral mass density |
| Vn (east) | m/s | Eastward neutral wind |
| Vn (north) | m/s | Northward neutral wind |
| Visc_Ve_rshear | m/s2 | Viscous east shear accel |
| Visc_Vn_rshear | m/s2 | Viscous north shear accel |
| IonDrag_east | m/s2 | Ion drag east accel |
| IonDrag_north | m/s2 | Ion drag north accel |
| Horiz_Adv_Ve | m/s2 | Horiz advection east accel |
| Horiz_Adv_Vn | m/s2 | Horiz advection north accel |
| PressGrad (east) | m/s2 | Neutral Pressure gradient acceleration eastward |
| PressGrad (north) | m/s2 | Neutral Pressure gradient acceleration northward |
| Coriolis_east | m/s2 | Coriolis east accel |
| Coriolis_north | m/s2 | Coriolis north accel |
| Centrifugal_north | m/s2 | Centrifugal north accel |

### 3DMOV

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Latitude | degrees_north | Geographic latitude |
| Altitude | m | Altitude above surface |
| Vn (up,O(3P)) | m/s | Vertical wind for O(3P) |
| IonDrag_up_O(3P) | m/s2 | Ion drag upward for O(3P) |
| Horiz_Adv_Vn (up,O(3P)) | m/s2 | Horiz advection upward for O(3P) |
| Vn (up,O2) | m/s | Vertical wind for O2 |
| IonDrag_up_O2 | m/s2 | Ion drag upward for O2 |
| Horiz_Adv_Vn (up,O2) | m/s2 | Horiz advection upward for O2 |
| Vn (up,N2) | m/s | Vertical wind for N2 |
| IonDrag_up_N2 | m/s2 | Ion drag upward for N2 |
| Horiz_Adv_Vn (up,N2) | m/s2 | Horiz advection upward for N2 |
| Vn (up,N(4S)) | m/s | Vertical wind for N(4S) |
| IonDrag_up_N(4S) | m/s2 | Ion drag upward for N(4S) |
| Horiz_Adv_Vn (up,N(4S)) | m/s2 | Horiz advection upward for N(4S) |
| Vn (up,NO) | m/s | Vertical wind for NO |
| IonDrag_up_NO | m/s2 | Ion drag upward for NO |
| Horiz_Adv_Vn (up,NO) | m/s2 | Horiz advection upward for NO |
| Vn (up,He) | m/s | Vertical wind for He |
| IonDrag_up_He | m/s2 | Ion drag upward for He |
| Horiz_Adv_Vn (up,He) | m/s2 | Horiz advection upward for He |
| Vn (up,CO2) | m/s | Vertical wind for CO2 |
| IonDrag_up_CO2 | m/s2 | Ion drag upward for CO2 |
| Horiz_Adv_Vn (up,CO2) | m/s2 | Horiz advection upward for CO2 |
| Coriolis_up | m/s2 | Coriolis upward accel |
| Centrif_up | m/s2 | Centrifugal upward accel |
| EffectiveGravity | m/s2 | Effective gravity accel |

### 3DNEU

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Latitude | degrees_north | Geographic latitude |
| Altitude | m | Altitude above surface |
| Rho | kg/m3 | Total neutral mass density |
| [O(3P)] | /m3 | Number density of O(3P) |
| [O2] | /m3 | Number density of O2 |
| [N2] | /m3 | Number density of N2 |
| [N(4S)] | /m3 | Number density of N(4S) |
| [NO] | /m3 | Number density of NO |
| [He] | /m3 | Number density of He |
| [CO2] | /m3 | Number density of CO2 |
| [N(2D)] | /m3 | Number density of N(2D) |
| [N(2P)] | /m3 | Number density of N(2P) |
| [H] | /m3 | Number density of H |
| [O(1D)] | /m3 | Number density of O(1D) |
| Temperature | K | Neutral temperature |
| Vn (east) | m/s | Eastward neutral wind |
| Vn (north) | m/s | Northward neutral wind |
| Vn (up) | m/s | Vertical neutral wind |
| Vn (up,O(3P)) | m/s | Vertical neutral wind for O(3P) |
| Vn (up,O2) | m/s | Vertical neutral wind for O2 |
| Vn (up,N2) | m/s | Vertical neutral wind for N2 |
| Vn (up,N(4S)) | m/s | Vertical neutral wind for N(4S) |
| Vn (up,NO) | m/s | Vertical neutral wind for NO |
| Vn (up,He) | m/s | Vertical neutral wind for He |
| Vn (up,CO2) | m/s | Vertical neutral wind for CO2 |

### 3DTHM

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Latitude | degrees_north | Geographic latitude |
| Altitude | m | Altitude above surface |
| EUV Heating Rate | K/s | EUV heating rate |
| Conduction | K/s | Thermal conduction heating rate |
| Molecular Conduction | K/s | Molecular conduction heating rate |
| Chemical Heating Rate | K/s | Chemical heating rate |
| Joule Heating Rate | K/s | Joule heating rate |
| CO2 Cooling Rate | K/s | CO2 infrared cooling rate |
| NO Cooling Rate | K/s | NO infrared cooling rate |
| O Cooling Rate | K/s | O infrared cooling rate |
| Total Abs EUV | W/m2 | Total absorbed EUV flux |
| Cp | J/kg/K | Specific heat at constant pressure |
| Rho | kg/m3 | Total neutral mass density |
| E.F. Magnitude | V/m | Electric field magnitude |
| Sigma Ped | S/m | Pedersen conductivity |
| Ionization Rate O_3P | /m3/s | Auroral ionization rate of O(3P) |
| Ionization Rate O2 | /m3/s | Auroral ionization rate of O2 |
| Ionization Rate N2 | /m3/s | Auroral ionization rate of N2 |
| Ionization Rate He | /m3/s | Auroral ionization rate of He |

### 3DUSR

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Latitude | degrees_north | Geographic latitude |
| Altitude | m | Altitude above surface |
| Joule Heating | K/s | Joule Heating |
| JPara | A/m2 | JPara |

### 2DANC

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Latitude | degrees_north | Geographic latitude |
| Altitude | m | Altitude above surface |
| LT | hr | Local Time |
| SZA | rad | Solar zenith angle |
| TEC | TECU | Total electron content |
| JH_int | W/m2 | Alt-integrated Joule heating |
| HT_int | W/m2 | Alt-integrated Heat transfer |
| EUV_int | W/m2 | Alt-integrated EUV heating |
| PE_int | W/m2 | Alt-integrated PE heating |
| Chem_int | W/m2 | Alt-integrated Chem heating |
| RadCool_int | W/m2 | Alt-integrated Rad cooling |
| CO2Cool_int | W/m2 | Alt-integrated CO2 cooling |
| NOCool_int | W/m2 | Alt-integrated NO cooling |
| OCool_int | W/m2 | Alt-integrated O cooling |

### 2DGEL

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Latitude | degrees_north | Geographic latitude |
| Altitude | m | Altitude above surface |
| Potential | V | Electric potential |
| PedCond | S | Height-integrated Pedersen conductance |
| HalCond | S | Height-integrated Hall conductance |
| AveE | keV | Diffuse electron precipitation average energy |
| eFlux | ergs/cm2/s | Diffuse electron precipitation energy flux |
| AveE_W | keV | Wave-driven electron precipitation average energy |
| eFlux_W | ergs/cm2/s | Wave-driven electron precipitation energy flux |
| AveE_M | keV | Monoenergetic electron precipitation average energy |
| eFlux_M | ergs/cm2/s | Monoenergetic electron precipitation energy flux |
| AveE_I | keV | Ion precipitation average energy |
| eFlux_I | ergs/cm2/s | Ion precipitation energy flux |
| DivJuAlt | A/m2 | Height-integrated divergence of horizontal current |
| PedFLCond | S | Field-line integrated Pedersen conductance |
| HalFLCond | S | Field-line integrated Hall conductance |
| DivJuFL | A/m2 | Field-line divergence of current |
| FL_length | m | Field line length |

### 2DHME

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Latitude | degrees_north | Geographic latitude |
| Altitude | m | Altitude above surface |
| LT | hr | Local Time |
| TEC | TECU | Total electron content |

### 2DMEL

| Variable | Units | Description |
|---|---|---|
| mlon | degrees | Magnetic longitude |
| mlat | degrees | Magnetic latitude |
| Altitude | m | IE reference altitude (constant over the magnetic grid) |
| MLT | hr | Magnetic local time |
| GeoLat | degrees_north | Geographic latitude |
| GeoLon | degrees_east | Geographic longitude |
| PedCond | S | Height-integrated Pedersen conductance |
| HalCond | S | Height-integrated Hall conductance |
| DivJuAlt | A/m2 | Height-integrated divergence of horizontal current |
| FL_length | m | Field line length |
| Sigma_PP | S | Sigma PP conductance component |
| Sigma_LL | S | Sigma LL conductance component |
| Sigma_H | S | Hall conductance component |
| Sigma_C | S | Cowling conductance component |
| Sigma_PL | S | Sigma PL conductance component |
| Sigma_LP | S | Sigma LP conductance component |
| K_mphi |  | Dynamo K_mphi coefficient |
| K_mlambda |  | Dynamo K_mlambda coefficient |
| Solver_A |  | Dynamo solver coefficient A |
| Solver_B |  | Dynamo solver coefficient B |
| Solver_C |  | Dynamo solver coefficient C |
| Solver_D |  | Dynamo solver coefficient D |
| Solver_E |  | Dynamo solver coefficient E |
| Solver_S |  | Dynamo solver coefficient S |
| Phi_dynamo | V | Dynamo electric potential |
| Ed1new | V/m | Dynamo electric field component 1 |
| Ed2new | V/m | Dynamo electric field component 2 |
| Kphi |  | Dynamo Kphi coefficient |
| Klambda |  | Dynamo Klambda coefficient |

### 2DTEC

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Latitude | degrees_north | Geographic latitude |
| Altitude | m | Altitude above surface |
| SZA | rad | Solar zenith angle |
| TEC | TECU | Total electron content |

### 2DUSR

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Latitude | degrees_north | Geographic latitude |
| Altitude | m | Altitude above surface |
| Potential (kV) | kV | Potential (kV) |
| Average Energy (keV) | keV | Average Energy (keV) |
| Total Energy (ergs) | ergs/cm2/s | Total Energy (ergs) |
| Discrete Average Energy (keV) | keV | Discrete Average Energy (keV) |
| Discrete Total Energy (ergs) | ergs/cm2/s | Discrete Total Energy (ergs) |
| Wave Average Energy (keV) | keV | Wave Average Energy (keV) |
| Wave Total Energy (ergs) | ergs/cm2/s | Wave Total Energy (ergs) |

### 1DALL

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Latitude | degrees_north | Geographic latitude |
| Altitude | m | Altitude above surface |
| Rho | kg/m3 | Total neutral mass density |
| [O(3P)] | /m3 | [O(3P)] |
| [O2] | /m3 | [O2] |
| [N2] | /m3 | [N2] |
| [N(4S)] | /m3 | [N(4S)] |
| [NO] | /m3 | [NO] |
| [He] | /m3 | [He] |
| [CO2] | /m3 | [CO2] |
| [N(2D)] | /m3 | [N(2D)] |
| [N(2P)] | /m3 | [N(2P)] |
| [H] | /m3 | [H] |
| [O(1D)] | /m3 | [O(1D)] |
| Temperature | K | Neutral temperature |
| Vn (east) | m/s | Eastward neutral wind |
| Vn (north) | m/s | Northward neutral wind |
| Vn (up) | m/s | Vertical neutral wind |
| Vn (up,O(3P)) | m/s | Vn (up,O(3P)) |
| Vn (up,O2) | m/s | Vn (up,O2) |
| Vn (up,N2) | m/s | Vn (up,N2) |
| Vn (up,N(4S)) | m/s | Vn (up,N(4S)) |
| Vn (up,NO) | m/s | Vn (up,NO) |
| Vn (up,He) | m/s | Vn (up,He) |
| Vn (up,CO2) | m/s | Vn (up,CO2) |
| [O(4S)+] | /m3 | [O(4S)+] |
| [He+] | /m3 | [He+] |
| [NO+] | /m3 | [NO+] |
| [O2+] | /m3 | [O2+] |
| [N2+] | /m3 | [N2+] |
| [N+] | /m3 | [N+] |
| [O(2D)+] | /m3 | [O(2D)+] |
| [O(2P)+] | /m3 | [O(2P)+] |
| [H+] | /m3 | [H+] |
| [e-] | /m3 | [e-] |
| eTemperature | K | Electron temperature |
| iTemperature | K | Ion temperature |
| Vi (east) | m/s | Eastward ion drift |
| Vi (north) | m/s | Northward ion drift |
| Vi (up) | m/s | Vertical ion drift |

### 1DCHM

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Latitude | degrees_north | Geographic latitude |
| Altitude | m | Altitude above surface |
| N2+ + e | eV/m3/s | N2+ + e |
| O2+ + e | eV/m3/s | O2+ + e |
| N2+ + O | eV/m3/s | N2+ + O |
| NO+ + e | eV/m3/s | NO+ + e |
| N+ + O2 | eV/m3/s | N+ + O2 |
| NO + N | eV/m3/s | NO + N |
| O+ + O2 | eV/m3/s | O+ + O2 |
| N + O2 | eV/m3/s | N + O2 |
| O2+ + N | eV/m3/s | O2+ + N |
| O2+ + NO | eV/m3/s | O2+ + NO |
| O2+ + N2 | eV/m3/s | O2+ + N2 |
| N2+ + O2 | eV/m3/s | N2+ + O2 |
| N+ + O | eV/m3/s | N+ + O |
| O+ + N2 | eV/m3/s | O+ + N2 |
| O(1D) + N2 | eV/m3/s | O(1D) + N2 |
| O(1D) + O2 | eV/m3/s | O(1D) + O2 |
| O(1D) + O | eV/m3/s | O(1D) + O |
| O(1D) + e | eV/m3/s | O(1D) + e |
| N(2D) + O2 | eV/m3/s | N(2D) + O2 |
| O+(2D) + e | eV/m3/s | O+(2D) + e |
| N(2D) + O | eV/m3/s | N(2D) + O |
| N(2D) + e | eV/m3/s | N(2D) + e |
| O+(2D) + N2 | eV/m3/s | O+(2D) + N2 |
| O+(2P) + e | eV/m3/s | O+(2P) + e |
| O+(2P) + O | eV/m3/s | O+(2P) + O |
| O+(2P) + N2 | eV/m3/s | O+(2P) + N2 |
| Chemical Heating | eV/m3/s | Volumetric chemical heating rate |

### 1DEMI

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Latitude | degrees_north | Geographic latitude |
| Altitude | m | Altitude above surface |
| Emission 2470A |  | Emission 2470A |
| Emission 7320A |  | Emission 7320A |
| Emission 3726A |  | Emission 3726A |
| Emission 5200A |  | Emission 5200A |
| Emission 10400A |  | Emission 10400A |
| Emission 6300A |  | Emission 6300A |
| Emission 6364A |  | Emission 6364A |

### 1DNEW

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Local Time | hr | Local Time |
| Latitude | degrees_north | Geographic latitude |
| Solar Zenith Angle | rad | Solar Zenith Angle |
| Altitude | m | Altitude above surface |
| Rho | kg/m3 | Total neutral mass density |
| [O(3P)] | /m3 | [O(3P)] |
| [O2] | /m3 | [O2] |
| [N2] | /m3 | [N2] |
| [N(4S)] | /m3 | [N(4S)] |
| [NO] | /m3 | [NO] |
| [He] | /m3 | [He] |
| [CO2] | /m3 | [CO2] |
| [N(2D)] | /m3 | [N(2D)] |
| [N(2P)] | /m3 | [N(2P)] |
| [H] | /m3 | [H] |
| [O(1D)] | /m3 | [O(1D)] |
| Temperature | K | Neutral temperature |
| Vn (east) | m/s | Eastward neutral wind |
| Vn (north) | m/s | Northward neutral wind |
| Vn (up) | m/s | Vertical neutral wind |
| Vn (up,O(3P)) | m/s | Vn (up,O(3P)) |
| Vn (up,O2) | m/s | Vn (up,O2) |
| Vn (up,N2) | m/s | Vn (up,N2) |
| Vn (up,N(4S)) | m/s | Vn (up,N(4S)) |
| Vn (up,NO) | m/s | Vn (up,NO) |
| Vn (up,He) | m/s | Vn (up,He) |
| Vn (up,CO2) | m/s | Vn (up,CO2) |
| [O(4S)+] | /m3 | [O(4S)+] |
| [He+] | /m3 | [He+] |
| [NO+] | /m3 | [NO+] |
| [O2+] | /m3 | [O2+] |
| [N2+] | /m3 | [N2+] |
| [N+] | /m3 | [N+] |
| [O(2D)+] | /m3 | [O(2D)+] |
| [O(2P)+] | /m3 | [O(2P)+] |
| [H+] | /m3 | [H+] |
| [e-] | /m3 | [e-] |
| eTemperature | K | Electron temperature |
| iTemperature | K | Ion temperature |
| Vi (east) | m/s | Eastward ion drift |
| Vi (north) | m/s | Northward ion drift |
| Vi (up) | m/s | Vertical ion drift |
| O(3P) Mixing Ratio |  | O(3P) Mixing Ratio |
| O2 Mixing Ratio |  | O2 Mixing Ratio |
| N2 Mixing Ratio |  | N2 Mixing Ratio |
| N(4S) Mixing Ratio |  | N(4S) Mixing Ratio |
| NO Mixing Ratio |  | NO Mixing Ratio |
| He Mixing Ratio |  | He Mixing Ratio |
| CO2 Mixing Ratio |  | CO2 Mixing Ratio |

### 1DTHM

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Latitude | degrees_north | Geographic latitude |
| Altitude | m | Altitude above surface |
| EUV Heating (K/s) | K/s | EUV Heating (K/s) |
| Conduction (K/s) | K/s | Conduction (K/s) |
| Molecular Conduction (K/s) | K/s | Molecular Conduction (K/s) |
| Chemical Heating (K/s) | K/s | Chemical Heating (K/s) |
| Joule Heating (K/s) | K/s | Joule Heating (K/s) |
| CO2 Cooling (K/s) | K/s | CO2 Cooling (K/s) |
| NO Cooling (K/s) | K/s | NO Cooling (K/s) |
| O Cooling (K/s) | K/s | O Cooling (K/s) |
| Total Abs EUV | W/m2 | Total absorbed EUV flux |
| Production Rate O(3P) | /m3/s | Production Rate O(3P) |
| Production Rate O2 | /m3/s | Production Rate O2 |
| Production Rate N2 | /m3/s | Production Rate N2 |
| Production Rate N(4S) | /m3/s | Production Rate N(4S) |
| Production Rate NO | /m3/s | Production Rate NO |
| Production Rate He | /m3/s | Production Rate He |
| Production Rate CO2 | /m3/s | Production Rate CO2 |
| Production Rate N(2D) | /m3/s | Production Rate N(2D) |
| Production Rate N(2P) | /m3/s | Production Rate N(2P) |
| Production Rate H | /m3/s | Production Rate H |
| Production Rate O(1D) | /m3/s | Production Rate O(1D) |
| Loss Rate O(3P) | /m3/s | Loss Rate O(3P) |
| Loss Rate O2 | /m3/s | Loss Rate O2 |
| Loss Rate N2 | /m3/s | Loss Rate N2 |
| Loss Rate N(4S) | /m3/s | Loss Rate N(4S) |
| Loss Rate NO | /m3/s | Loss Rate NO |
| Loss Rate He | /m3/s | Loss Rate He |
| Loss Rate CO2 | /m3/s | Loss Rate CO2 |
| Loss Rate N(2D) | /m3/s | Loss Rate N(2D) |
| Loss Rate N(2P) | /m3/s | Loss Rate N(2P) |
| Loss Rate H | /m3/s | Loss Rate H |
| Loss Rate O(1D) | /m3/s | Loss Rate O(1D) |

### 1DUSR

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Latitude | degrees_north | Geographic latitude |
| Altitude | m | Altitude above surface |
| Electron Density | /m3 | Electron Density |

### 0DALL

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Latitude | degrees_north | Geographic latitude |
| Altitude | m | Altitude above surface |
| Rho | kg/m3 | Total neutral mass density |
| [O(3P)] | /m3 | [O(3P)] |
| [O2] | /m3 | [O2] |
| [N2] | /m3 | [N2] |
| [N(4S)] | /m3 | [N(4S)] |
| [NO] | /m3 | [NO] |
| [He] | /m3 | [He] |
| [CO2] | /m3 | [CO2] |
| [N(2D)] | /m3 | [N(2D)] |
| [N(2P)] | /m3 | [N(2P)] |
| [H] | /m3 | [H] |
| [O(1D)] | /m3 | [O(1D)] |
| Temperature | K | Neutral temperature |
| Vn (east) | m/s | Eastward neutral wind |
| Vn (north) | m/s | Northward neutral wind |
| Vn (up) | m/s | Vertical neutral wind |
| Vn (up,O(3P)) | m/s | Vn (up,O(3P)) |
| Vn (up,O2) | m/s | Vn (up,O2) |
| Vn (up,N2) | m/s | Vn (up,N2) |
| Vn (up,N(4S)) | m/s | Vn (up,N(4S)) |
| Vn (up,NO) | m/s | Vn (up,NO) |
| Vn (up,He) | m/s | Vn (up,He) |
| Vn (up,CO2) | m/s | Vn (up,CO2) |
| [O(4S)+] | /m3 | [O(4S)+] |
| [He+] | /m3 | [He+] |
| [NO+] | /m3 | [NO+] |
| [O2+] | /m3 | [O2+] |
| [N2+] | /m3 | [N2+] |
| [N+] | /m3 | [N+] |
| [O(2D)+] | /m3 | [O(2D)+] |
| [O(2P)+] | /m3 | [O(2P)+] |
| [H+] | /m3 | [H+] |
| [e-] | /m3 | [e-] |
| eTemperature | K | Electron temperature |
| iTemperature | K | Ion temperature |
| Vi (east) | m/s | Eastward ion drift |
| Vi (north) | m/s | Northward ion drift |
| Vi (up) | m/s | Vertical ion drift |
| O(3P) Mixing Ratio |  | O(3P) Mixing Ratio |
| O2 Mixing Ratio |  | O2 Mixing Ratio |
| N2 Mixing Ratio |  | N2 Mixing Ratio |
| N(4S) Mixing Ratio |  | N(4S) Mixing Ratio |
| NO Mixing Ratio |  | NO Mixing Ratio |
| He Mixing Ratio |  | He Mixing Ratio |
| CO2 Mixing Ratio |  | CO2 Mixing Ratio |
| RadCooling | K/s | RadCooling |
| EuvHeating | K/s | EuvHeating |
| Conduction | K/s | Thermal conduction heating rate |
| Heat Balance Total | K/s | Heat Balance Total |
| Heating Efficiency |  | Heating Efficiency |

### 0DUSR

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Latitude | degrees_north | Geographic latitude |
| Altitude | m | Altitude above surface |
| Electron Density | /m3 | Electron Density |
| Electron Temperature | K | Electron Temperature |
| Ion Temperature | K | Ion Temperature |

<!-- END GENERATED TABLES -->
