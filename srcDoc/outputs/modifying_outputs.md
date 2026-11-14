# Adding New Outputs

There are three ways to get new data into GITM's output files, all needing a
recompile:

| You want to... | Use | Effort |
|---|---|---|
| Add a variable to an existing output type | Edit its `fill_<type>` routine | One `put_field` call plus one registry row (see [Modify an existing output type](#modify-an-existing-output-type)) |
| Add a variable without changing any existing file's format | A **USR** output (`3DUSR`, `2DUSR`, `1DUSR`, `0DUSR`) | Register the name once, plus a one-line `add_usr_output` call (0D assigns `UserData1D` instead; see below) |
| Create a whole new output file | A **new output type** | One fill routine plus one registration line |

USR outputs are the low-friction path for personal diagnostics: they leave
the standard types' file formats (and everyone else's analysis scripts)
untouched. Adding to an existing type is just as easy, but it changes that
type's variable list for every user of the file.

[TOC]

## How the output system works

GITM's output system is built around **containers** (`src/ModOutputContainer.f90`):
one per output type, holding the variable schema (names, units, shapes), the
data, and a pointer to the type's `fill_<type>` routine.
For most types the schema is derived from the fill itself.
On initialization, each `put_field` defines its variable and schema ("one-shot").
The exceptions are the HME (regional), magnetic-grid (different block layout) and USR (not all variables are known at init) types.
They declare their schema explicitly in a `define_schema_<type>` routine (the authoritative list is the `define_schema_*` calls in `init_output_containers`).
The `ModOutputContainer.f90` header describes both patterns, if you need more information.

The backends then consume the container directly, so headers, NetCDF metadata, and on-disk layout all come from one schema and cannot drift from the data.
This prevents any of the backends (`legacy`, `mpiio`, `netcdf`) writing a specific output type or variable differently than the rest.

The recipes below need no more than this. [Under the hood](#under-the-hood) has the details.

!!! note "Coming from the old output system (the per-type `output_*` routines in `output_common.f90`)"
    If you last added a variable before the container refactor, four things
    changed:

    - The per-type header `write` lists are gone. A type's variable names,
      units, and long names come from its `fill_<type>` routine plus the
      shared table in `src/ModOutputVarRegistry.f90`. Headers and netcdf
      metadata are generated from that one schema, so they cannot disagree
      with the data.
    - Units are a separate field, not part of the name string (`'Rho'`
      with units `kg/m3`, not `'Rho (kg/m3)'`). `.header` variable lists
      show bare names.
    - The old conditional-write idiom (gating a variable on an
      `i*OutputList`-style flag) now aborts at startup unless the flag is
      fixed before init: the schema is derived by running each fill once at
      init ("priming"), and every later fill must put exactly the same
      variables. See [Containers](#containers).
    - `1DCMS` was removed, and type codes are matched exactly (5
      characters), not by substring. An old `UAM.in` requesting `1DCMS`
      or a truncated code stops at input parsing with `bad output type`.

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
    reliable way. Keep on-the-fly registration for throwaway debugging.

---

## Modify an existing output type

The best-case scenario: one line of code plus one metadata row.

1. In the type's `fill_<type>` routine (`src/ModOutputProducers.f90`), add a
   `put_field` call where the variable belongs in the output order:

   ```fortran
   call c%put_field('My Variable', MyVar(:, :, :, iBlock))
   ```

   Match the slicing of the puts already in that fill, which follows the
   type's `nGhostCells`: `fill_3dall` (`nGhostCells = 2`) puts whole
   ghosted arrays (`Altitude_GB(:, :, :, iBlock)`) while `fill_2dtec`
   (`nGhostCells = 0`) puts explicit interior slices
   (`Altitude_GB(1:nLons, 1:nLats, 1, iBlock)`). Position matters too: put
   order is the column order of legacy/mpiio files, so inserting mid-list
   renumbers every later column for scripts that index by position.
   Appending at the end is the safe default.

2. Give the name metadata: either add a row to the table in
   `src/ModOutputVarRegistry.f90` (preferred for fixed names: units and
   long_name live in one shared place), or pass it inline for computed names:

   ```fortran
   call c%put_field('['//trim(cSpecies(i))//']', NDensityS(:, :, :, i, iBlock), units='/m3')
   ```

   A registry row looks like this. The table is one long Fortran array
   constructor, so when you append, the previous last row needs a trailing
   `, &` and the closing `]` moves to your new row. It is assumed-shape,
   so there is no count to bump:

   ```fortran
   VarMeta('My Variable', 'units', 'One-line description')]
   ```

A `put_field` of a name with no registry row and no inline `units=` is a hard
`stop_gitm`: units cannot be guessed. Every variable must be put exactly once
per fill, unconditionally (see [Containers](#containers) under
[Under the hood](#under-the-hood) for why).

For the types that declare schemas with `define_schema_<type>` instead of
priming (see [below](#types-with-explicit-schemas)), also add a matching
`define_var` line in schema order.

### Changing an existing name or units

Renaming is safe to attempt: change the registry row and the put site
together, and if you miss one, the run stops at startup. Priming runs every
registered fill at init, so a put-site name with no metadata is caught on
the very next run, not months later on the one run that outputs that type.

Two things to check anyway:

- **A registry row is shared by every type that puts that name.** This
  keeps units and long_name from drifting between types, but editing the
  `'Rho'` row changes `3DALL`, `3DNEU`, `3DLST`, `3DMOH` and `3DTHM` at
  once, and a units-only edit produces no error.
  Grep the put sites before editing a shared row.
- Renames propagate to files and user scripts: the `.header` variable
  list, the netcdf variable name after sanitization (the original name
  survives in the `gitm_name` attribute described in [Outputs](../outputs.md)), and any
  analysis code keyed on the old string.

---

## Create a new output type

Use this only when you want a new file type, not just more columns in a USR
file. You write **one fill routine** in `src/ModOutputProducers.f90` and register it.
The schema is derived from the fill itself. **No need to modify the output backends:**
every backend reads the container directly for the data. (File headers and
file opens do check the type code in a few places: `2DMEL`, `1DNEW`, and the
`1D`/`0D` prefixes in the header writers in `ModOutputBackend.f90` and
`ModOutputMPIIO.f90`, and `2DMEL`/`HME` in the file open in `output_common.f90`.
A new magnetic-grid or regional type would get a wrong `BLOCKS` header or file
open until those checks are extended.)

The example below adds a 3D geographic type called `3DMYT`.

### 1. Write the fill

Set the type identity, `prepare`, then `put_field`/`put_axis` each variable in
output order. Fixed names get their units/long_name from
`src/ModOutputVarRegistry.f90` rows; computed names pass `units=` inline:

```fortran
subroutine fill_3dmyt(c, iBlock)
  use ModGITM
  type(OutputContainer), intent(inout) :: c
  integer, intent(in) :: iBlock

  c%cType = '3DMYT'
  c%gridKind = GRID_GEO_3D     ! geographic, all ranks write
  c%nGhostCells = 2            ! 2 = include ghost cells, 0 = interior only

  call c%prepare(iBlock)
  if (.not. c%this_rank_writes) return

  ! Longitude/Latitude go in as radians; the netcdf backend converts them to degrees
  call c%put_axis('Longitude', Longitude(:, iBlock))
  call c%put_axis('Latitude', Latitude(:, iBlock))
  call c%put_field('Altitude', Altitude_GB(:, :, :, iBlock))
  call c%put_field('My Variable', MyVar(:, :, :, iBlock), &
                   units='units', longName='Description of my variable')
end subroutine fill_3dmyt
```

The first execution of the fill (a "priming" pass at init) defines each
variable and locks the schema, so a bad put stops the first test run instead
of writing a malformed file. The rules are under [Containers](#containers).

### 2. Register the type

In `src/ModOutputProducers.f90`, add `'3DMYT'` to the `cContainerTypes`
array **and bump its declared size** (`cContainerTypes(28)` becomes `(29)`).
The constructor's length must match the declaration or the build fails
with a shape mismatch. Then add one line in `init_output_containers`:

```fortran
call set_fill('3DMYT', fill_3dmyt)
```

Type codes are exactly 5 characters and are matched exactly against
`#SAVEPLOT` entries. With the count bumped, that's the whole wiring: the
lookup functions and the dispatch in `output_common.f90` are data-driven,
and `assert_types_fully_wired` stops the run at init if anything else is
missing.

### 3. Enable the new type and rebuild

Add `3DMYT` to a `#SAVEPLOT` block in `UAM.in` and `make`. All three backends
pick up the new type automatically, and so does `pGITM.py` when it
post-processes `legacy` or `mpiio` output. Nothing else changes in
the code, but do add the new code to the `#SAVEPLOT` type list in
[set_inputs](../set_inputs.md) (and to the `#SAVEPLOT` help text in
`src/set_inputs.f90`), and [update the documentation](#update-the-documentation).

If the new type is a `1D*`/`0D*` probe type, the netcdf backend will refuse it
(see [Output backends](output_backends.md) for what each backend supports).

### Types with explicit schemas

A few types cannot be primed at init: magnetic-grid types read lazily
allocated arrays, USR schemas depend on runtime registration, and the
rank-0 probe types share the magnetic-grid participation rule. The
authoritative list is the set of `define_schema_*` calls in
`init_output_containers` (currently 3DHME, 2DHME, 1DGLO, 1DTHM, 1DCHM, 1DEMI,
2DMEL, and the four USR types). Those keep
the older two-routine pattern: a `define_schema_<type>` (called in
`init_output_containers`) declaring every variable with `define_var`, and a
fill using the plain `put` API against the declared schema. Only use this
pattern when priming is impossible. `define_schema_2dmel` is a worked
example, and [Under the hood](#under-the-hood) lists `define_var`'s
arguments.

---

## Update the documentation

The variable tables in [Outputs](../outputs.md) are generated, not typed.
At its first output GITM writes `UA/data/output_schema.csv`: one row
per variable of every output type, whatever the backend and whether or not
`#SAVEPLOT` asked for the type. A copy of that file is kept in
`srcDoc/outputs/output_schema.csv` and the tables are built from the copy.

After adding, renaming, or removing a variable or a type, run any case with
the new executable and then, from the top of the repository:

```bash
python srcPython/generate_output_docs.py run/UA/data/output_schema.csv
```

This replaces the copy and rewrites the tables in `outputs.md`. Commit both
files. A new type also needs one hand-written row in the type list at the top
of that section; the script tells you if the page never names the type.

The copy describes an Earth build with the `user.f90` in the repository.
Another planet, or user code that registers its own USR variables, writes a
different file, so only update the copy from a build of the code you are
committing.

`srcTests/auto_test/run_all_tests.sh` runs the same script with `--check`
after the tests and fails if the copy or the tables no longer match what the
code writes, so the GitHub tests catch a forgotten update.

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
so none of them need per-type code for the data (see the note on headers
under [Create a new output type](#create-a-new-output-type)).

Priming is why every variable must be put exactly once per fill, unconditionally.
The first execution of the fill runs at init to derive the schema; after that,
`finish_fill`/`finish_priming` aborts on a later fill that puts an unknown name,
a shape mismatch, or a different count. Two things follow: don't gate a `put_field`
on anything that varies by rank or timestep (that would desync the schema across
ranks or steps), and since priming runs in every configuration, the fill must only
touch arrays that are always allocated.

The precise form of the no-conditional-puts rule: a put may be gated only on
state that is identical on every writing rank and fixed before priming runs.
`fill_3dlst` is the sanctioned example. Its `i*OutputList` flags are set
once in `src/Earth.f90` before init and never change, so its conditional
puts resolve identically at priming and at every later fill. A flag that can
change mid-run, or differ across ranks, aborts in `finish_fill`. A container
holds at most `MaxContainerVars` (300) variables; the overflow message names
the parameter to raise (`src/ModOutputContainer.f90`).

### gridKind

`gridKind` tells the backends which ranks have data for a given timestep:

| Constant | Used by | Who writes |
|---|---|---|
| `GRID_GEO_3D` | geographic 3D (3DALL, 3DNEU, etc.) | all ranks |
| `GRID_GEO_2D` | geographic 2D (2DTEC, 2DGEL, etc.) and the all-rank probe types (1DALL, 0DALL, 1DNEW, 1DUSR, 0DUSR) | all ranks |
| `GRID_MAG_2D` | magnetic 2D (2DMEL), plus the probe types 1DGLO, 1DTHM, 1DCHM, 1DEMI, which use it purely as a "one rank writes" switch | `iProc == 0` only |
| `GRID_HIME` | regional (2DHME, 3DHME) | blocks inside the region (set externally) |

For 1D/0D types the choice is a participation rule, not a grid statement:
pick `GRID_MAG_2D` for a probe-like type that should write once from one
rank, `GRID_GEO_2D` only if every rank should contribute its own values.

For output types where not all blocks write, the `netcdf` backend automatically switches from
collective to **independent I/O mode** (`nfmpi_begin_indep_data` / `nfmpi_end_indep_data`):

- **`GRID_MAG_2D`** (e.g., 2DMEL): output is on the magnetic grid; only `iProc == 0`
  has the data and writes. `prepare()` sets `this_rank_writes` accordingly.
- **`GRID_HIME`** (e.g., 2DHME/3DHME): regional output; only blocks whose cells fall
  inside the user-defined HIME region write. Participation is set externally before
  `prepare()` is called.

The switch exists because a collective call that only some ranks reach is a
deadlock; independent mode lets the participating ranks write alone.

The `mpiio` backend handles conditional participation by striping output by global
block offset and leaving non-participating slots zero-filled (no separate per-block
file to mark NaN). For these two grid kinds it writes with independent
`MPI_File_write_at` rather than the collective `MPI_File_write_at_all`. A `0`
value in a regional `mpiio` output file is therefore expected for cells outside
the region. `validate_outputs.py` recognizes HME/HIME files and skips the
`altitude > 0` sanity check accordingly. The `legacy` backend skips file
creation on non-participating ranks.

### Ghost cells

Each block carries **ghost cells**, a two-cell border copied from neighboring blocks
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

### When it stops: the error messages

The output system fails loudly at startup rather than writing malformed
files. The messages you are most likely to meet:

| Message (abridged) | Cause | Fix |
|---|---|---|
| `put of unknown variable "X" in type Y after the schema was primed. A converted fill_<type> must put the same variables on every rank and every timestep; no conditional puts.` | A `put_field` gated on something that changed since init, or a name edited in the fill but not re-primed | Make the put unconditional (or gate only on pre-init-fixed state); rebuild after name changes |
| `variable "X" in type Y has no ModOutputVarRegistry entry; pass units= at the put site or register it` | A new or renamed name with no metadata | Add the registry row, or `units=` inline for computed names |
| `fill_X made N put(s) but the type has M variable(s); every variable must be put exactly once, unconditionally.` | A put skipped (conditional or early return) or duplicated relative to the primed schema | Every variable, exactly once, on every writing rank |
| `put shape mismatch for "X" ... defined as AxBxC, got DxExF` | The put's slice doesn't match the primed shape | Match the slicing of the fill's other puts (check `nGhostCells`) |
| `output container has no fill registered (forgot its set_fill line?): X` | A type added to `cContainerTypes` without its `set_fill` line | Add the line in `init_output_containers` |
| `ModOutputNetCDF: probe type X is not supported by the netcdf backend yet; use legacy or mpiio` | A `1D*`/`0D*` type requested under the netcdf backend, including 3D types renamed to `1D*` by a 1-D (`nLons=nLats=1`) run | Use the legacy or mpiio backend for probe types and 1-D runs |
