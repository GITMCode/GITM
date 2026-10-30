# Outputs {#outputs.sec}

Now that you have managed to successfully complete a GITM run you've
found yourself with a bunch of output files. All of the GITM output is
in MKS units. Output files are written to the directory given as the first
argument of `#DIRECTORIES` in `UAM.in` (default: `UA/data`). Satellite files
(`#SATELLITES`) always go to `UA/data`, whatever `#DIRECTORIES` says, and
are written in the legacy format whatever the backend.

!!! important "Horizontal coordinates: radians in `.bin`, degrees in netcdf"
    The one exception to "MKS units" is geographic longitude and latitude,
    and it differs **by backend**:

    | Backend | `Longitude` / `Latitude` |
    |---|---|
    | `legacy`, `mpiio` (`.bin`) | **radians**, exactly as the model holds them |
    | `netcdf` | **degrees** (`degrees_east` / `degrees_north`) |

    The variable tables below are generated from netcdf files, so the units
    they show for `Longitude`/`Latitude` are the netcdf ones. Reading a
    `.bin` file, convert with `numpy.degrees` yourself. This is what
    `pGITM.py` does when it writes netcdf from legacy/mpiio output, but a
    `.bin` file merged by `pGITM.py` stays in radians. 2DMEL's `GeoLat` and
    `GeoLon` follow the same rule (radians in `.bin`, degrees in netcdf).
    Only the magnetic coordinates `mlon`/`mlat` are in degrees on every
    backend.

The format of those files depends on which [output backend](outputs/output_backends.md)
was used during the run (`#OUTPUTBACKEND`; an unrecognized value falls
back to `legacy` with a warning). You will end up with a `log*.dat`
file and the output files for each type listed in `#SAVEPLOT` (see
[`#SAVEPLOT`](set_inputs.md#saveplot)). `netcdf` outputs are usable as
written. `legacy` writes per-block files that need
[postprocessing](postprocessing.md) to be merged, and `mpiio` writes one
raw `.bin` per output, which `pGITM.py` rewrites in the standard layout
(keeping the original as `.raw` unless it is removing files).

!!! warning "Don't mix the two netcdf dialects"
    `post_process.py`/`pGITM.py` can also convert legacy/mpiio output to
    netcdf, but the files it produces are **not** the same format the
    `netcdf` backend writes: variable names come from a different
    vocabulary (short names, no `gitm_name` attribute), the dimension
    order is reversed, edge ghost cells are kept on some types, and there
    are no CF attributes. Tools written for one dialect will not read the
    other. If you run the `netcdf` backend, use its files directly and
    skip the netcdf conversion in postprocessing; convert with pGITM only
    when starting from legacy/mpiio output.

The log file provides useful information about the run, such as whether a
restart was performed, which physical processes were used, and a list of some
variables which represent the state of the modeled system and its drivers. This
file can be very useful when sharing runs with other users, when revisiting an
old run, or merely ensuring that GITM performed as expected. Some example log
files can be found in the `srcTests/auto_test/ref_solns/` directory.

## Output system

GITM's output system is built around **containers** (`src/ModOutputContainer.f90`):
one per output type, holding the variable schema (names, units, shapes), the
data, and a pointer to the type's `fill_<type>` routine. For most types the
schema is derived from the fill itself. It runs once at init and each first
`put_field` defines its variable ("one-shot"). The rest declare their schema
explicitly in a `define_schema_<type>` routine instead (the authoritative
list is the `define_schema_*` calls in `init_output_containers`: the
magnetic-grid, HME, USR, and rank-0 probe types); the
`ModOutputContainer.f90` header describes both patterns. Either way the backends consume the container
directly, so headers, NetCDF metadata, and on-disk layout all come from one
schema and cannot drift from the data. None of the three backends
(`legacy`, `mpiio`, `netcdf`) has per-type code paths for the data itself.
(File headers and file opens do key on a few type codes, such as `2DMEL`,
`1DNEW`, and the `1D`/`0D` probe prefixes; see
[Adding new outputs](outputs/modifying_outputs.md).) Adding or modifying a
type touches `src/ModOutputProducers.f90` plus, for fixed variable names, a
metadata row in `src/ModOutputVarRegistry.f90`. See
[Adding new outputs](outputs/modifying_outputs.md) for recipes and
[Output backends](outputs/output_backends.md) for the backend descriptions.

## Output types and their variables

The tables below are generated from GITM output by
`srcPython/generate_output_docs.py`: variable names, units, and descriptions
are read straight from the netcdf files, not typed by hand. To regenerate
after changing an output type, run a sweep with the types you care about on
the `netcdf` backend, then splice the script's stdout into this page,
replacing the generated `### <TYPE>` sections (everything from the first one
down to the generation marker at the bottom of the file) and leaving the
hand-written text above them alone:

```bash
python srcPython/generate_output_docs.py <run>/UA/data
```

The `Variable` column is the GITM name as it appears in legacy/mpiio
`.header` files (netcdf variable names are the same after sanitization; the
raw name is kept in each variable's `gitm_name` attribute). Notes:

- **USR types** (`3DUSR`/`2DUSR`/`1DUSR`/`0DUSR`) list whatever the run's
  user code registered. The tables below show one example configuration.
- **Not listed** (absent from the generating run): `3DHME` (regional,
  3DALL-like variable set), the 1-D probe types (`1DALL`, `1DTHM`, `1DCHM`,
  `1DEMI`, `1DNEW`, `1DUSR`), and the 0-D satellite types (`0DALL`, `0DUSR`). The
  probe types halt on the `netcdf` backend and the generator reads only
  netcdf files, so they cannot have a generated table. For their variable
  list, read the `.header` file from a `legacy` or `mpiio` run.
- **Rates are zero at the first output.** Per-step diagnostics that GITM
  divides by the time step (3DTHM `Conduction`, `Molecular Conduction` and
  the chemical heating rate, 3DCHM, 3DMOH `Visc_*`, 2DANC `Chem_int`, and
  the probe types 0DALL, 1DTHM and 1DCHM) are rates. They are written as 0
  at the t = 0 output, where the time step is 0.
- **`3DGLO` and `1DGLO` cannot be generated at all.** GLOW computes
  emissions when `#GLOW` is on (`get_glow.f90` fills `vEmissionRate`), but
  the two output types were never wired to that array: their fills put
  placeholder zeros, indistinguishable from real measurements, so GITM
  stops if either is requested in `#SAVEPLOT`. They remain valid type
  codes so the fills can be pointed at `vEmissionRate` later; they produce
  nothing today.
- **`2DMEL` is on a magnetic grid**, so its horizontal coordinates are named
  `mlon`/`mlat` (in degrees), not `lon`/`lat`. Its `GeoLat`/`GeoLon`
  variables are the geographic position of each magnetic grid point (radians
  in `.bin`, degrees in netcdf).

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
| Ionization Rate O_3P | /s | Ionization rate for O(3P) |
| Ionization Rate O2 | /s | Ionization rate for O2 |
| Ionization Rate N2 | /s | Ionization rate for N2 |
| Ionization Rate He | /s | Ionization rate for He |

### 3DUSR

| Variable | Units | Description |
|---|---|---|
| Longitude | degrees_east | Geographic longitude |
| Latitude | degrees_north | Geographic latitude |
| Altitude | m | Altitude above surface |
| Joule Heating | K/s | Joule Heating |
| JPara | A/m2 | JPara |

<!-- generated by srcPython/generate_output_docs.py; types: 2DANC, 2DGEL, 2DHME, 2DMEL, 2DTEC, 2DUSR, 3DALL, 3DCHM, 3DEMI, 3DION, 3DLST, 3DMAG, 3DMOH, 3DMOV, 3DNEU, 3DTHM, 3DUSR -->
