# Adding New Outputs

There are two ways to get new data into GITM's output files, both needing a
recompile:

| You want to... | Use | Effort |
|---|---|---|
| Add a variable to an existing output | A **USR** output (`3DUSR`, `2DUSR`, `1DUSR`, `0DUSR`) | Register the name once, plus a one-line `add_usr_output` call |
| Create a whole new output file | A **new output type** | One fill routine plus one registration line |

Adding to a USR output is by far the common case. Start there.

[TOC]

---

## Add a variable to a USR output

The USR outputs exist so you can add a variable without writing any new output
machinery. Register the name once, store the data each step with a single
`add_usr_output` call, and recompile.

### 1. Register the name

In `init_usr_output_registry` (`src/user.f90`), add one line per variable. The first
argument is the dimensionality (3, 2, 1, or 0). Longitude/Latitude/Altitude are added
for you, so don't list them.

```fortran
call register_usr_var(3, 'My Variable', 'units')
```

Registering here matters because a file's variable list is fixed when the file is first
written. The netcdf backend defines all variables at file creation and never adds more.
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
cells automatically, and matches the variable by name, so registration order never
affects your output columns. 1D variables take no `iBlock`. The `Joule Heating`
output in `src/add_sources.f90` is a working example.

0D is the one exception: there is no `add_usr_output` path for scalars.
`fill_0dusr` reads the registered 0D variables directly from `UserData1D` at
the satellite's interpolation point, so a 0D variable only needs the
`register_usr_var(0, ...)` line and a `UserData1D` assignment.

You can also assign the shared `UserData*D` arrays directly (slot number =
registration order), but nothing checks the slot, and a mismatch silently
swaps output columns. `add_usr_output` avoids that trap by matching on name.
`nUserOutputs` in `src/ModUser.f90` (default 100) caps how many variables
each USR dimensionality can hold; raise it and recompile if you hit it.

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

!!! tip "Registering on the fly"
    If you call `add_usr_output` with a name that was never registered, it registers it
    on that first call. Handy for quick diagnostics, but that first call must happen
    **before the first output write**, or the variable won't make it into the file
    schema (and will be missing from netcdf output). Registering in step 1 is the
    reliable way; on-the-fly registration is best left for throwaway debugging.

---

## Modify an existing output type

The best-case scenario: one line of code plus one metadata row.

1. In the type's `fill_<type>` routine (`src/ModOutputProducers.f90`), add a
   `put_field` call where the variable belongs in the output order:

   ```fortran
   call c%put_field('My Variable', MyVar(:, :, :, iBlock))
   ```

2. Give the name metadata: either add a row to the table in
   `src/ModOutputVarRegistry.f90` (preferred for fixed names: units and
   long_name live in one shared place), or pass it inline for computed names:

   ```fortran
   call c%put_field('['//trim(cSpecies(i))//']', NDensityS(:, :, :, i, iBlock), units='/m3')
   ```

A `put_field` of a name with no registry row and no inline `units=` is a hard
`stop_gitm`: units cannot be guessed. Every variable must be put exactly once
per fill, unconditionally (see [Containers](#containers) under
[Under the hood](#under-the-hood) for why).

For the few types still using `define_schema_<type>` (magnetic-grid, HME,
USR, see below), also add a matching `define_var` line in schema order.

---

## Create a new output type

Use this only when you want a genuinely new file type, not just more columns in a USR
file. You write **one fill routine** in `src/ModOutputProducers.f90` and register it;
the schema is derived from the fill itself. **No need to modify the output backends:**
every backend reads the container directly.

The example below adds a 3-D geographic type called `3DMYT`.

### 1. Write the fill

Set the type identity, `prepare`, then `put_field`/`put_axis` each variable in
output order. Fixed names get their units/long_name from
`src/ModOutputVarRegistry.f90` rows; computed names pass `units=` inline:

```fortran
subroutine fill_3dmyt(c, iBlock)
  use ModGITM
  use ModConst, only: cRadToDeg
  type(OutputContainer), intent(inout) :: c
  integer, intent(in) :: iBlock

  c%cType = '3DMYT'
  c%gridKind = GRID_GEO_3D     ! geographic, all ranks write
  c%nGhostCells = 2            ! 2 = include ghost cells, 0 = interior only

  call c%prepare(iBlock)
  if (.not. c%this_rank_writes) return

  call c%put_axis('Longitude', Longitude(:, iBlock)*cRadToDeg)
  call c%put_axis('Latitude', Latitude(:, iBlock)*cRadToDeg)
  call c%put_field('Altitude', Altitude_GB(:, :, :, iBlock))
  call c%put_field('My Variable', MyVar(:, :, :, iBlock), &
                   units='units', longName='Description of my variable')
end subroutine fill_3dmyt
```

The first execution of the fill (a "priming" pass at init) defines each
variable on first touch and locks the schema. After that, a `put_field` of an
unknown name, a shape mismatch, or a missing/extra put is a `stop_gitm`, so a
test run catches the mistake immediately instead of writing a malformed file.
See [Containers](#containers) under [Under the hood](#under-the-hood) for the
full priming rule (every variable put exactly once, unconditionally).

### 2. Register the type

In `src/ModOutputProducers.f90`, add `'3DMYT'` to the `cContainerTypes`
array, and one line in `init_output_containers`:

```fortran
call set_fill('3DMYT', fill_3dmyt)
```

That's the whole wiring: the lookup functions and the dispatch in
`output_common.f90` are data-driven, and `assert_types_fully_wired` stops the
run at init if anything is missing.

### 3. Turn it on and rebuild

Add `3DMYT` to a `#SAVEPLOT` block in `UAM.in` and `make`. `pGITM.py` and all three
backends pick up the new type automatically. There is nothing else to change.

### Types with explicit schemas (mag-grid / HIME / USR)

A few types cannot be primed at init: magnetic-grid types read lazily
allocated arrays, and USR schemas depend on runtime registration. Those keep
the older two-routine pattern: a `define_schema_<type>` (called in
`init_output_containers`) declaring every variable with `define_var`, and a
fill using the plain `put` API against the declared schema. Only use this
pattern when priming is genuinely impossible; see `define_schema_2dmel` for
a worked example and the [Under the hood](#under-the-hood) section for
`define_var`'s arguments.

---

## Under the hood

You don't need any of this to add a variable. It explains *why* the code is shaped the
way it is, and details the options used above.

### Containers

Every output type is one **container** (`type(OutputContainer)`, in
`src/ModOutputContainer.f90`) holding its 5-character type code, a `gridKind`, an
array of variables (each with a name, units, shape, and data), and a pointer
to its fill routine. The containers live in the `containers(:)` array in
`src/ModOutputProducers.f90` (slot order = `cContainerTypes`), are built once
(`init_output_containers`), and are reused every timestep; only the data changes.

Most types are **one-shot**: their schema is materialized by running the fill
once at init and locking it (`finish_priming`). The remaining types declare
schemas explicitly via `define_schema_<type>`. Either way, `output_common.f90`
calls `containers(i)%fill(...)` each output step and hands the container to
the active backend's `write_container`; backends read the container directly,
so none of them need per-type code.

Priming is why every variable must be put exactly once per fill, unconditionally.
The first execution of the fill runs at init to derive the schema; after that,
`finish_fill`/`finish_priming` aborts on a later fill that puts an unknown name,
a shape mismatch, or a different count. Two things follow: don't gate a `put_field`
on anything that varies by rank or timestep (that would desync the schema across
ranks or steps), and since priming runs in every configuration, the fill must only
touch arrays that are always allocated.

### gridKind

`gridKind` tells the backends which ranks have data for a given timestep:

| Constant | Grid | Who writes |
|---|---|---|
| `GRID_GEO_3D` | geographic 3D (3DALL, 3DNEU, etc.) | all ranks |
| `GRID_GEO_2D` | geographic 2D / 1D / 0D (2DTEC, 2DGEL, etc.) | all ranks |
| `GRID_MAG_2D` | magnetic 2D (2DMEL) | `iProc == 0` only |
| `GRID_HIME` | regional (2DHME, 3DHME) | blocks inside the region (set externally) |

For `GRID_MAG_2D` and `GRID_HIME`, the netcdf backend switches to independent I/O to
avoid MPI collective deadlocks; see
[Output backends](output_backends.md#conditional-participation-types).

### Ghost cells

Each block carries **ghost cells**, a two-cell border copied from neighbouring blocks
(or boundary conditions) so stencil operations don't have to special-case the edges.
Interior cells are `1:nLons, 1:nLats, 1:nAlts`; with ghosts the range is `-1:nLons+2`,
and so on.

`c%nGhostCells` sets whether ghosts are written: `0` (interior only, all 2D, 1D, 0D,
magnetic, and USR types) or `2` (ghost-inclusive, most 3D geographic types). With
`nGhostCells = 2`, every variable's shape is padded by `2*nGC` and the fill routine
puts ghost-inclusive slices; the netcdf backend strips them back off for geographic
grids. For a variable that is only valid on the interior, build a padded array with
clamped reads before putting it.

### Metadata: the registry, `put_field`, and `define_var`

Variable metadata lives in two places by design. **Fixed names** get one row
in `src/ModOutputVarRegistry.f90`, shared by every type that outputs that
name, so units/long_name can't drift between types. **Computed names** (the
species/ion loops) pass `units=`/`longName=` at the put site, since the
registry can't enumerate them.

`put_field(name, arr[, units, longName])` defines the variable on first touch
(shape from the array, metadata from the registry unless overridden inline)
and stores the data; `put_axis` is the same but marks a 1-D coordinate
variable, which netcdf writes as a coordinate while legacy and mpiio
replicate it. For define-path types, `define_var(name, units, shape3[,
longName, is_axis])` declares each variable up front (`shape3` pads unused
dims with `1`, e.g. `[nLons, nLats, 1]` for 2D), and the fill's plain `put`
must then match name and shape exactly.

Only netcdf writes `units` and `long_name` to files; legacy and mpiio use just the
name and shape. By convention a container defines its axes first (Longitude, then
Latitude), then Altitude, then the physics variables.
