# Adding New Outputs

To add new outputs to GITM there are two options.

1. Add variables to one of the existing USR output types (`3DUSR`, `2DUSR`, `1DUSR`,
   or `0DUSR`).
1. Define an entirely new output type.

Both require changes to the Fortran source, so a recompile is needed before the new
outputs will appear.

[TOC]

---

## Background concepts

### Ghost cells

GITM's grid uses **ghost cells** — extra cells surrounding each block that are filled
with data from neighbouring blocks (or boundary conditions) to allow stencil operations
without special-casing the edges. The interior cells of a block are indexed
`1:nLons`, `1:nLats`, `1:nAlts`. Ghost cells extend that to `-1:nLons+2` in each
horizontal direction and `-1:nAlts+2` vertically.

Whether an output type includes ghost cells in the output file is set when the type is
registered. This is captured in the `nGhostCells` field:

| `nGhostCells` | Buffer size | Index loop |
|---|---|---|
| `2` | `(nLons+4) × (nLats+4) × (nAlts+4)` | `-1:nLons+2`, `-1:nLats+2`, `-1:nAlts+2` |
| `0` | `nLons × nLats` | `1:nLons`, `1:nLats` |

Most 3D types include ghost cells (`nGhostCells=2`) so that adjacent blocks share
boundary data. Most 2D types do not (`nGhostCells=0`).

Some physics variables are only computed on the interior (`1:nLons`, etc.) and do not
have meaningful values at ghost indices. When writing ghost cells for such a variable,
clamp the physics index to the interior range rather than reading out-of-bounds:

```fortran
iiLon = max(min(iLon, nLons), 1)   ! clamp to interior
buffer(iv, jx, jy, jz) = MyVar(iiLon, iiLat, iiAlt, iBlock)
```

The ghost cell slots in the output will then hold the value of the nearest interior
point, which is conventional for GITM outputs. See `gather_3dthm` in
`src/ModOutputGather.f90` for a complete example of this pattern.

### The output buffer

The gather routines fill a contiguous array called `buffer` with shape
`(nV, nX, nY, nZ)`, where:

- `nV` — number of output variables (must match the registry for this type)
- `nX`, `nY`, `nZ` — spatial dimensions, determined by `nGhostCells`

Variables are the **leading** (fastest-varying) dimension so that all values for one
spatial point sit next to each other in memory. `nZ=1` for 2D output types.

The buffer is 1-indexed: spatial index `iLon=-1` maps to `jx=1`, `iLon=0` maps to
`jx=2`, `iLon=1` maps to `jx=3`, and so on. In general:

```
jx = iLon + nGhostCells
jy = iLat + nGhostCells
jz = iAlt + nGhostCells
```

A clean pattern for filling the buffer sequentially is an `iv` counter (used throughout
`ModOutputGather.f90`):

```fortran
iv = 0
iv = iv + 1; buffer(iv, jx, jy, jz) = Longitude(iLon, iBlock)
iv = iv + 1; buffer(iv, jx, jy, jz) = Latitude(iLat, iBlock)
iv = iv + 1; buffer(iv, jx, jy, jz) = Altitude_GB(iLon, iLat, iAlt, iBlock)
iv = iv + 1; buffer(iv, jx, jy, jz) = MyVar(iLon, iLat, iAlt, iBlock)
```

This keeps the code in sync with the registry declaration order and is easy to extend.

### The registry and the `info` pointer

The registry is a fixed-size array of `OutputTypeInfo` descriptors (`RegisteredTypes`
in `ModOutputRegistry.f90`). `new_output_type()` allocates the next free slot and
returns a Fortran pointer, `info`, that points directly to it. After that call you
write variable metadata into it via `add_coord_vars(info)` and `add_var(info, ...)`.
No manual indexing is needed — `info` is just a convenient alias for the slot.

The complete set of fields on an `OutputTypeInfo` is:

| Field | Type | Meaning |
|---|---|---|
| `code` | `character(5)` | The 5-character type identifier used in `UAM.in` |
| `nDims` | `integer` | Spatial dimensionality: 3 (lon×lat×alt), 2 (lon×lat), 1 (alt), or 0 (point) |
| `nGhostCells` | `integer` | Ghost cells per side: 2 (3D types) or 0 (2D/1D/0D types) |
| `nVars` | `integer` | Number of variables registered so far (updated by `add_var`) |
| `isUserType` | `logical` | Marks USR types; variables are registered at runtime by `init_usr_output_registry` |
| `usesMagGrid` | `logical` | True if output is on magnetic (not geographic) coordinates. Only block 1 writes. (e.g., 2DMEL) |
| `isRegional` | `logical` | True if output is limited to a user-defined region. Only blocks within region write. (e.g., 2DHME/HIME) **For netcdf backend, both `usesMagGrid` and `isRegional` trigger independent I/O mode.** |
| `vars(:)` | `OutputVar` | Variable names and units, in output order |

---

## Modifying USR outputs

The USR output types (`3DUSR`, `2DUSR`, `1DUSR`, `0DUSR`) let you add any output
without touching the core output infrastructure. Three shared arrays in the
`src/user.f90` file hold the data:

| Array | Declared shape | Ghost cells? | Used by |
|---|---|---|---|
| `UserData3D` | `(-1:nLons+2, -1:nLats+2, -1:nAlts+2, nUserOutputs, nBlocksMax)` | Yes | `3DUSR` |
| `UserData2D` | `(-1:nLons+2, -1:nLats+2, 1, nUserOutputs, nBlocksMax)` | No[^gc2d] | `2DUSR` |
| `UserData1D` | `(1, 1, 1:nAlts, nUserOutputs)` | No | `1DUSR` and `0DUSR`[^0dusr] |

[^gc2d]: `UserData2D` is declared with ghost-cell bounds for consistency with other
    GITM arrays, but the gather routine only reads indices `1:nLons`, `1:nLats`. You
    only need to populate the interior.

[^0dusr]: `0DUSR` is a single-point output sampled along a satellite track. It reads
    from `UserData1D(1,1,1,1:nVarsUser0d-3)` — one value per variable at altitude
    index 1.

All arrays are initialised to `0.0` at startup. Ghost cell slots that you do not
populate will be zero in the output, which is acceptable for most uses.

To add a new variable to a USR output, follow these steps:

### 1. Populate the array

At the point in the physics code where your variable is available, copy it into the
appropriate `UserData` array. The **fourth index** selects which output column it will
occupy. Variable `n` in `UserData3D` (or `2D`) becomes output column `n+3` (after
Longitude, Latitude, and Altitude).

**3DUSR example** — outputting a 3D variable:

```fortran
use ModUserGITM, only: UserData3D

! Variable index 1 → output column 4
UserData3D(1:nLons, 1:nLats, 1:nAlts, 1, iBlock) = MyVariable(1:nLons, 1:nLats, 1:nAlts)
```

`3DUSR` includes ghost cells in the file, but you only need to fill the interior. The
ghost cell slots are already `0.0`. If your variable has meaningful ghost cell values
(e.g., it was filled during ghost cell exchanges), you may populate those too:

```fortran
UserData3D(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 1, iBlock) = &
  MyVariable(-1:nLons+2, -1:nLats+2, -1:nAlts+2)
```

**2DUSR example** — outputting a 2D (altitude-integrated or surface) variable:

```fortran
use ModUserGITM, only: UserData2D

! Variable index 1 → output column 4
UserData2D(1:nLons, 1:nLats, 1, 1, iBlock) = My2DVariable(1:nLons, 1:nLats)
```

The `2DUSR` gather only reads interior cells, so there is no need to fill ghost cells
regardless of how `UserData2D` is declared.

!!! tip "Where to put the assignment"
    The assignment should happen **after** the variable it copies from is computed and
    before the next output timestep. A common pattern is to place it at the end of
    whichever subroutine computes the variable (e.g., `calc_electrodynamics.f90`,
    `calc_all_fluxes_hydro.f90`), inside the block loop.

### 2. Register the variable

Open `src/user.f90` and find `init_usr_output_registry`. For your output type (`3DUSR`,
`2DUSR`, etc.), `idx` holds the index of that type in the registry. Call `add_var()` for
each new variable **in the same order as the array indices** you populated in step 1:

```fortran
! This snippet is from init_usr_output_registry in src/user.f90
idx = find_output_type('3DUSR')
if (idx > 0 .and. RegisteredTypes(idx)%nVars == 0) then
  call add_var(RegisteredTypes(idx), 'Longitude', 'rad')
  call add_var(RegisteredTypes(idx), 'Latitude', 'rad')
  call add_var(RegisteredTypes(idx), 'Altitude', 'm')
  call add_var(RegisteredTypes(idx), 'Joule Heating', 'K/s')
  call add_var(RegisteredTypes(idx), 'My Variable Name', 'units')  ! YOUR NEW VARIABLE
```

!!! danger "Variable order matters"
    The order of `add_var()` calls **must exactly match** the order of array indices
    you used when populating `UserData*D`. If you populate things out of order, the
    output columns will be scrambled. There is no error message.

The name you give here appears in the `.header` file and in NetCDF metadata. Use the
same naming style as existing entries in that section.

### 3. Update the variable count

In the same file, find the corresponding `set_nVarsUser3d` (or `2d`, `1d`, `0d`)
routine and update the count to reflect the new total. **Include the three coordinate
variables** (Longitude, Latitude, Altitude) in this count:

- 1 user variable → `nVarsUser3d = 4`
- 2 user variables → `nVarsUser3d = 5`
- and so on.

!!! warning
    If you skip this step or get the count wrong, GITM will output the wrong number of
    columns and post-processors will fail to read the file. The count **must** match
    the number of `add_var()` calls you made in step 2.

### 4. Enable the output type in UAM.in

Make sure the relevant type code (`3DUSR`, `2DUSR`, etc.) appears in `UAM.in` under
`#SAVEPLOT`. For example, to output `3DUSR` every 60 seconds:

```
#SAVEPLOT
3600.0    (Seconds between writing restart files)
1         (number of output types, adjust as necessary)
3DUSR     (output type)
60.0      (output interval in seconds)
```

See [set_inputs](../set_inputs.md) for full `#SAVEPLOT` syntax.

### 5. Recompile

```bash
make
```

!!! tip
    `nUserOutputs` in `src/ModUser.f90` sets the maximum number of user-defined
    variables across all USR types (default: 100). If you need more than 97 variables
    (100 minus the 3 coordinates), increase `nUserOutputs` and recompile.

---

## Creating new output types

To define a completely new output type, two files need to be changed. The backends
handle headers, file I/O, and format selection automatically once the type is
registered.

### 1. Register the type

In `src/ModOutputRegistry.f90`, write a `register_XXXX()` subroutine and call it from
`init_output_registry()`. Use `new_output_type()` to allocate a registry slot, then
`add_coord_vars()` and `add_var()` to declare each variable in output order:

```fortran
subroutine register_3dmytype()
  type(OutputTypeInfo), pointer :: info

  ! new_output_type arguments: code, nDims, nGhostCells, info
  !   code        -- 5-character identifier used in UAM.in
  !   nDims       -- 3 for lon×lat×alt, 2 for lon×lat, 1 for alt-only, 0 for point
  !   nGhostCells -- 2 to include ghost cells in output, 0 for interior only
  !   info        -- pointer to the new registry slot; use it directly below
  call new_output_type('3DMYT', 3, 2, info)
  ! Optional flags (if your type has special behavior):
  ! info%usesMagGrid = .true.   ! if output is on magnetic (not geographic) coordinates
  ! info%isRegional = .true.    ! if only blocks within a user-defined region write
  call add_coord_vars(info)            ! adds Longitude, Latitude, Altitude
  call add_var(info, 'My Variable', 'units')
  call add_var(info, 'Another Var',  'units')
end subroutine register_3dmytype
```

Then add the call inside `init_output_registry()`:

```fortran
call register_3dmytype()
```

The 5-character type code (`3DMYT` here) is what users put in `UAM.in`.

!!! info "Conditional participation flags"
    Set flags for output types where not all blocks write:
    
    - **`usesMagGrid = .true.`** — Output is on magnetic coordinates; only block 1 has the data. Example: 2DMEL
    - **`isRegional = .true.`** — Output is limited to a user-defined region; only blocks in that region write. Example: 2DHME (HIME region)
    
    When either flag is set, the **netcdf backend automatically switches to independent I/O mode** during writes to avoid MPI collective operation deadlocks. Other backends (legacy binary, mpiio) are unaffected but will respect the conditional participation rules.

### 2. Write the gather routine

In `src/ModOutputGather.f90`, write a subroutine that fills `buffer` from the model
state, and add a `case` for it in the `gather_output` dispatcher:

```fortran
case ('3DMYT'); call gather_3dmyt(iBlock, buffer, nV, nX, nY, nZ)
```

The `iv` counter pattern (used throughout `ModOutputGather.f90`) keeps variable
assignment in step with the registry declaration order and makes it easy to add or
remove variables without renumbering:

```fortran
subroutine gather_3dmyt(iBlock, buffer, nV, nX, nY, nZ)
  use ModGITM
  integer, intent(in) :: iBlock, nV, nX, nY, nZ
  real, intent(out)   :: buffer(nV, nX, nY, nZ)
  integer :: iAlt, iLat, iLon, jx, jy, jz, iv
  integer :: iiAlt, iiLat, iiLon   ! clamped interior indices (see below)

  do iAlt = -1, nAlts + 2
    jz = iAlt + 2
    iiAlt = max(min(iAlt, nAlts), 1)   ! clamp to interior
    do iLat = -1, nLats + 2
      jy = iLat + 2
      iiLat = max(min(iLat, nLats), 1)
      do iLon = -1, nLons + 2
        jx = iLon + 2
        iiLon = max(min(iLon, nLons), 1)
        iv = 0
        iv = iv + 1; buffer(iv, jx, jy, jz) = Longitude(iLon, iBlock)
        iv = iv + 1; buffer(iv, jx, jy, jz) = Latitude(iLat, iBlock)
        iv = iv + 1; buffer(iv, jx, jy, jz) = Altitude_GB(iLon, iLat, iAlt, iBlock)
        ! Physics variables: use clamped indices if the variable is only
        ! defined on the interior (1:nLons, 1:nLats, 1:nAlts).
        iv = iv + 1; buffer(iv, jx, jy, jz) = MyVar(iiLon, iiLat, iiAlt, iBlock)
        ! If the variable IS defined at ghost cells, use the un-clamped index:
        iv = iv + 1; buffer(iv, jx, jy, jz) = MyGhostedVar(iLon, iLat, iAlt, iBlock)
      enddo
    enddo
  enddo
  ! Optional but recommended: check the final iv matches the registry count.
  if (iv /= nV) call gather_error('gather_3dmyt', iv, nV)
end subroutine gather_3dmyt
```

!!! warning "Variable order matters"
  Take note of the order you added variables to the registry.
  The variables must be placed in the buffer with the same order.

  It is not required to fill the buffer in order (i.e. you can manually set `iV`),
  but filling the buffer in the same order that the variables were defined helps
  make sure things are consistent.

For a **2D output type** (registered with `nDims=2`, `nGhostCells=0`), there is no
altitude loop and the buffer size is `nLons × nLats`:

```fortran
subroutine gather_2dmyt(iBlock, buffer, nV, nX, nY, nZ)
  use ModGITM
  integer, intent(in) :: iBlock, nV, nX, nY, nZ
  real, intent(out)   :: buffer(nV, nX, nY, nZ)   ! nZ=1 for 2D types
  integer :: iLat, iLon, iv

  do iLat = 1, nLats
    do iLon = 1, nLons
      iv = 0
      iv = iv + 1; buffer(iv, iLon, iLat, 1) = Longitude(iLon, iBlock)
      iv = iv + 1; buffer(iv, iLon, iLat, 1) = Latitude(iLat, iBlock)
      iv = iv + 1; buffer(iv, iLon, iLat, 1) = Altitude_GB(iLon, iLat, 1, iBlock)
      iv = iv + 1; buffer(iv, iLon, iLat, 1) = My2DVar(iLon, iLat, iBlock)
    enddo
  enddo
  if (iv /= nV) call gather_error('gather_2dmyt', iv, nV)
end subroutine gather_2dmyt
```

### 3. Add to UAM.in

Add the new type code to `#SAVEPLOT` in `UAM.in`. No other changes to the input file
are needed — the registry metadata is used to generate headers automatically.

### 4. Recompile

```bash
make
```

!!! note
    You do not need to modify any backend code. The new output type will work with the
    `legacy`, `mpiio`, and `netcdf` backends without any further changes.
