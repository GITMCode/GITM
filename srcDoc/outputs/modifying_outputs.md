# Adding New Outputs

There are two ways to get new data into GITM's output files:

| You want to… | Use | Effort |
|---|---|---|
| Add a variable to an existing run's diagnostics | a **USR** output (`3DUSR`, `2DUSR`, `1DUSR`, `0DUSR`) | register the name once + a one-line `add_usr_output` call |
| Create a whole new output type with its own files | a **new output type** (e.g. `3DMYT`) | a schema/fill pair + one dispatch entry |

Both need a recompile. The first is by far the common case — start there.

[TOC]

---

## Add a variable to a USR output

The USR outputs exist so you can add a variable without writing any new output
machinery. Register the name once, store the data each step with a single
`add_usr_output` call, and recompile.

### 1. Register the name

In `init_usr_output_registry` (`src/user.f90`), add one line per variable. The first
argument is the dimensionality (3, 2, 1, or 0). Longitude/Latitude/Altitude are added
for you — don't list them.

```fortran
call register_usr_var(3, 'My Variable', 'units')
```

Registering here matters because a file's variable list is fixed when the file is first
written — the netcdf backend defines all variables at file creation and never adds more.
`init_usr_output_registry` runs before the first output, so anything registered here is
guaranteed to be in the schema.

### 2. Store the data

Wherever the variable is computed (or wherever you want it output from), hand it to
 `add_usr_output` **by name**.

```fortran
use ModOutputProducers, only: add_usr_output
call add_usr_output(MyVariable, 'My Variable', 'units', iBlock=iBlock)
```

`add_usr_output` picks the right storage from the array's rank (3D/2D/1D), strips ghost
cells automatically, and matches the variable by name — so registration order never
affects your output columns. 1D and 0D variables take no `iBlock`. The `Joule Heating`
output in `src/add_sources.f90` is a working example.

### 3. Turn it on and rebuild

Add the type to the `#SAVEPLOT` block in `UAM.in`, then `make`:

```
#SAVEPLOT
3600.0    (seconds between restart files)
1         (number of output types)
3DUSR     (output type)
60.0      (output interval, seconds)
```

See [set_inputs](../set_inputs.md) for the full `#SAVEPLOT` syntax.

!!! note "Lower-level alternative: write the array directly"
    `add_usr_output` is a wrapper over shared `UserData*D` arrays, which you can also
    assign yourself:

    ```fortran
    use ModUserGITM, only: UserData3D
    ! slot = registration order: the 1st register_usr_var(3,…) is slot 1, etc.
    UserData3D(1:nLons, 1:nLats, 1:nAlts, slot, iBlock) = MyVariable(1:nLons, 1:nLats, 1:nAlts)
    ```

    Only do this if you have a reason to manage slots by hand. The slot number must
    match the registration order — nothing checks it, and a mismatch silently swaps
    output columns. `add_usr_output` avoids that trap by matching on name.

!!! tip "Registering on the fly"
    If you call `add_usr_output` with a name that was never registered, it registers it
    on that first call. Handy for quick diagnostics — but that first call must happen
    **before the first output write**, or the variable won't make it into the file
    schema (and will be missing from netcdf output). Registering in step 1 is the
    reliable way; on-the-fly registration is best left for throwaway debugging.

!!! info "How many variables can I add?"
    `nUserOutputs` in `src/ModUser.f90` (default 100) caps how many variables each USR
    dimensionality can hold — it sizes the `UserData*D` arrays. If you hit the limit,
    raise it and recompile. That's the only knob.

---

## Create a new output type

Use this only when you want a genuinely new file type, not just more columns in a USR
file. You write two routines in `src/ModOutputProducers.f90` — a **schema** (lists the
variables) and a **fill** (copies data each dtOutput) — and add one entry to the
dispatch. **No need to modify the output backends:** every backend reads the
container directly.

The example below adds a 3-D geographic type called `3DMYT`.

### 1. Reserve a slot

Near the top of `src/ModOutputProducers.f90`, add an index constant and bump the count:

```fortran
integer, parameter :: iCont_3DMYT = 27
integer, parameter :: nContainers = 27   ! was 26
```

Then add `3DMYT` to `is_known_output_type` and `get_container_idx` so the rest of the
code recognises it.

### 2. Write the schema

List each variable once, in the order it should appear. `gridKind` and `nGhostCells`
describe the grid (see [Under the hood](#under-the-hood) for the full menu of options):

```fortran
subroutine define_schema_3dmyt(c)
  use ModGITM, only: nLons, nLats, nAlts
  type(OutputContainer), intent(inout) :: c
  integer :: nGC

  c%cType       = '3DMYT'
  c%gridKind    = GRID_GEO_3D     ! geographic, all ranks write
  c%nGhostCells = 2               ! 2 = include ghost cells, 0 = interior only
  nGC = c%nGhostCells

  call c%define_var('Longitude', units='degrees_east',  shape3=[nLons + 2*nGC, 1, 1], is_axis=.true.)
  call c%define_var('Latitude',  units='degrees_north', shape3=[nLats + 2*nGC, 1, 1], is_axis=.true.)
  call c%define_var('Altitude',  units='m', &
                    shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
  call c%define_var('My Variable', units='units', &
                    shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                    longName='Description of my variable')
end subroutine define_schema_3dmyt
```

Register it in `init_output_containers`:

```fortran
call define_schema_3dmyt(containers(iCont_3DMYT))
```

### 3. Write the fill

`prepare` decides whether this rank writes; then `put` each variable by name. Shapes
must match the schema, and each array must be `real(output_kind)`:

```fortran
subroutine fill_3dmyt(c, iBlock)
  use ModGITM
  use ModConst, only: cRadToDeg
  type(OutputContainer), intent(inout) :: c
  integer, intent(in) :: iBlock
  integer :: nGC

  call c%prepare(iBlock)
  if (.not. c%this_rank_writes) return
  nGC = c%nGhostCells

  call c%put('Longitude', real(Longitude(1-nGC:nLons+nGC, iBlock) * cRadToDeg, output_kind))
  call c%put('Latitude',  real(Latitude(1-nGC:nLats+nGC, iBlock) * cRadToDeg, output_kind))
  call c%put('Altitude',  real(Altitude_GB(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
  call c%put('My Variable', real(MyVar(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
end subroutine fill_3dmyt
```

If a name or shape doesn't match the schema, `put` calls `stop_gitm` — so a test run
catches the mistake immediately instead of writing a malformed file.

### 4. Hook into the dispatch

In `src/output_common.f90`, add a branch to the `if (cType == '...')` chain (around
line 265), following the existing entries:

```fortran
elseif (cType == '3DMYT') then
  call fill_3dmyt(containers(iCont_3DMYT), iBlock)
  if (associated(ActiveBackend%write_container)) &
    call ActiveBackend%write_container(containers(iCont_3DMYT), iBLK, iTimeArray)
  call containers(iCont_3DMYT)%reset()
```

### 5. Turn it on and rebuild

Add `3DMYT` to a `#SAVEPLOT` block in `UAM.in` and `make`. `pGITM.py` and all three
backends pick up the new type automatically — there is nothing else to change.

---

## Under the hood

You don't need any of this to add a variable. It explains *why* the code is shaped the
way it is, and details the options used above.

### Containers

Every output type is one **container** (`type(OutputContainer)`, in
`src/ModOutputContainer.f90`) holding its 5-character type code, a `gridKind`, and an
array of variables — each with a name, units, shape, and data. The containers live in
the `containers(:)` array in `src/ModOutputProducers.f90`, are built once
(`init_output_containers`), and are reused every timestep; only the data changes.

Each type pairs two routines: `define_schema_<type>` runs once to declare the
variables, and `fill_<type>` runs each timestep to copy data in. After the fill,
`output_common.f90` hands the container to the active backend's `write_container`.
Backends read the container directly, so none of them need per-type code.

### gridKind

`gridKind` tells the backends which ranks have data for a given timestep:

| Constant | Grid | Who writes |
|---|---|---|
| `GRID_GEO_3D` | geographic 3D (3DALL, 3DNEU, …) | all ranks |
| `GRID_GEO_2D` | geographic 2D / 1D / 0D (2DTEC, 2DGEL, …) | all ranks |
| `GRID_MAG_2D` | magnetic 2D (2DMEL) | `iProc == 0` only |
| `GRID_HIME` | regional (2DHME, 3DHME) | blocks inside the region (set externally) |

For `GRID_MAG_2D` and `GRID_HIME`, the netcdf backend switches to independent I/O to
avoid MPI collective deadlocks; see
[Output backends](output_backends.md#conditional-participation-types).

### Ghost cells

Each block carries **ghost cells** — a two-cell border copied from neighbouring blocks
(or boundary conditions) so stencil operations don't have to special-case the edges.
Interior cells are `1:nLons, 1:nLats, 1:nAlts`; with ghosts the range is `-1:nLons+2`,
and so on.

`c%nGhostCells` sets whether ghosts are written: `0` (interior only — all 2D, 1D, 0D,
magnetic, and USR types) or `2` (ghost-inclusive — most 3D geographic types). With
`nGhostCells = 2`, every variable's `shape3` is padded by `2*nGC` and the fill routine
`put`s ghost-inclusive slices; the netcdf backend strips them back off for geographic
grids. For a variable that is only valid on the interior, build a padded array with
clamped reads before `put`-ing it.

### `define_var` arguments

```fortran
call c%define_var(name, units, shape3, longName=…, is_axis=…)
```

- **`name`** — variable identifier. Appears in legacy/mpiio `.header` files and as the
  netCDF variable name; `put` looks variables up by exact match.
- **`units`** — unit string, written as the netCDF `units` attribute. Only the netcdf
  backend reads it, but it is not optional — pass `''` if you truly don't care.
- **`shape3`** — three dimension sizes, padding unused dims with `1` (e.g.
  `[nLons, nLats, 1]` for 2D). Drives allocation and `put` shape validation.
- **`longName`** *(optional)* — human-readable description; becomes the netCDF
  `long_name`. Defaults to `name`.
- **`is_axis`** *(optional)* — `.true.` for coordinate variables (Longitude, Latitude,
  …). These are 1D (`[N, 1, 1]`); netcdf writes them as coordinate variables, while
  legacy and mpiio replicate them across the grid.

Only netcdf reads `units` and `longName`; legacy and mpiio use just `name` and
`shape3`. By convention a container defines its axes first (Longitude, then Latitude),
then Altitude, then the physics variables.
