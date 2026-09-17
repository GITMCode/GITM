# Output Backends

GITM can write outputs using three different backends. All three consume the same
output [containers](../outputs.md#output-system) (defined in
`src/ModOutputProducers.f90`) and produce the same scientific data; they differ
only in how files are laid out on disk. All three are supported by `pGITM.py`.
Which one you use depends on your system and workflow.

Two current restrictions to know up front:

- The `netcdf` backend does not support the probe-style `1D*`/`0D*` types
  (including 3D types renamed to `1D*` by a 1-D `nLons=nLats=1` run) and
  stops with an error naming the type; use `legacy` or `mpiio` for those.
- The `.nc` files this backend writes are a different dialect from the ones
  `pGITM.py` produces when converting legacy/mpiio output — see the warning
  in [Outputs](../outputs.md). Don't mix the two in one analysis.

To select a backend, add this block to `UAM.in`:

```
#OUTPUTBACKEND
mpiio
```

Without this block, GITM uses the **legacy** backend to ensure backwards compatibility
with existing runs and tools.

[TOC]

---

## legacy

The legacy backend is the original GITM output format. Each MPI rank writes its own
block file independently, producing:

- One ASCII `.header` file (written by rank 0)
- One binary `.b####` file per MPI rank (e.g., `.b0001`, `.b0002`, ...)

These per-block files must be merged into a single `.bin` file before they can be
analysed. `post_process.py` handles this step automatically. See
[Postprocessing](../postprocessing.md) for details.

!!! warning "Lustre and parallel filesystems"
    On Lustre-based HPC systems, having hundreds of MPI ranks each create small files
    simultaneously can saturate the metadata server and slow the run significantly.
    If this is a concern, consider switching to the `mpiio` backend.

---

## mpiio

The MPI-IO backend writes all blocks into a single shared file using MPI collective
I/O (`MPI_File_open` / `MPI_File_write_at`). Each rank writes its block(s)
at a fixed byte offset, so all ranks write simultaneously without
serialization. The result is:

- One ASCII `.header` file (same format as legacy, with an added `FORMAT/mpiio` line)
- One `.bin` file containing all blocks in a single raw-binary file (no Fortran
  record markers)

`pGITM.py` detects the `FORMAT/mpiio` marker in the header and reads accordingly, so
no special post-processing step is required to merge files. Running `post_process.py`
is still recommended if you are using the standard Python analysis tools.

This backend is generally preferred on Lustre-based systems: it produces far fewer
files per timestep and benefits from large sequential writes rather than many small
independent ones.

---

## netcdf

The NetCDF backend uses PnetCDF to write a single self-describing file per
output type per timestep. Each MPI rank writes its interior (non-ghost) grid slices
directly into the correct position in the global array, and ghost cells are stripped
automatically.

The result is a single `.nc` file per output type per timestep. No post-processing or
merging is required. Files can be read directly with xarray, NCO, PyITM, or any
NetCDF-aware tool.

### Multi-time files, run names, and restarts

By default the NetCDF backend creates one file per output type per timestep
(e.g., `3DALL_t230101_000000.nc`, `3DALL_t230101_000500.nc`, ...). `#NETCDFOPTIONS`
controls this with two sequential flags (the first required, the second
optional):

```
#NETCDFOPTIONS
single          netCDFAppendMethod ('single'/'daily'/'none'=default)
storm1          netCDFFileName (optional run name)
```

`netCDFAppendMethod` selects the file layout:

- `none` (default): one file per output type per timestep.
- `daily`: one file per output type per day, e.g. `3DALL_2002-12-21.nc`.
- `single`: one file per output type for the entire run, e.g. `3DALL.nc`.

`daily` and `single` files use an unlimited time dimension; data variables are
shaped `(time, [z,] lat, lon)` (time first) and read natively by xarray
(on a magnetic grid the horizontal dimensions are named `mlat`/`mlon`
instead, and carry plain `degrees`):

```python
import xarray as xr
ds = xr.open_dataset("3DALL.nc", decode_times=True)
```

`netCDFFileName`, when set, is inserted into every netcdf filename after the
type code (`3DALL_storm1.nc`, `3DALL_storm1_t021221_000000.nc`), in every
append mode. This lets several scenario runs share one directory without
colliding.

A restarted `single`/`daily` run with a run name appends a
`restart_yymmdd_hhmmss` stamp (the restart's simulation time) to the name.
Each leg of the run writes its own file and prior legs are never reopened,
so a crash can only affect the current leg. Consecutive legs share the
restart epoch (it is re-emitted), so analysis tools should drop duplicate
times, keeping the later leg.

!!! warning "Overwrite guard"
    In `single`/`daily` mode, a run that is a restart or carries a run name
    refuses to overwrite an existing target file and stops with instructions
    (move/delete the file, change `netCDFFileName`, or use a fresh run
    directory). An unnamed, non-restart run keeps the historical behavior of
    overwriting. An unnamed restart always trips the guard, so set
    `netCDFFileName` to continue a run across restarts. `none` mode and the
    `legacy`/`mpiio` backends are unguarded: colliding epochs are recomputed
    and overwritten in place. (The recompute is not guaranteed bit-identical:
    restarted runs are known to diverge slightly, likely via the dynamo /
    high-latitude electric fields.)

`#NETCDFOPTIONS` only affects the `netcdf` backend; it is ignored for
`legacy` and `mpiio`.

### Conditional participation types

For output types where not all blocks write, the NetCDF backend automatically switches from
collective to **independent I/O mode** (`nfmpi_begin_indep_data` / `nfmpi_end_indep_data`).
Participation is driven by the container's `gridKind` (see
`src/ModOutputContainer.f90`):

- **`GRID_MAG_2D`** (e.g., 2DMEL): output is on the magnetic grid; only `iProc == 0`
  has the data and writes. `prepare()` sets `this_rank_writes` accordingly.
- **`GRID_HIME`** (e.g., 2DHME/3DHME): regional output; only blocks whose cells fall
  inside the user-defined HIME region write. Participation is set externally before
  `prepare()` is called.

The switch exists because a collective call that only some ranks reach is a
deadlock; independent mode lets the participating ranks write alone.

The `mpiio` backend handles conditional participation by striping output by global
block offset and leaving non-participating slots zero-filled (no separate per-block
file to mark NaN). A `0` value in a regional `mpiio` output file is therefore
expected for cells outside the region. `validate_outputs.py` recognises HME/HIME
files and skips the `altitude > 0` sanity check accordingly. The `legacy` backend
simply skips file creation on non-participating ranks.

!!! note "Build requirements"
    PnetCDF must be installed and `pnetcdf-config` on `PATH` when compiling;
    the build system detects it and defines `-DHavePNetCDF`. Without it the
    NetCDF module compiles as a stub and GITM stops if `netcdf` is requested
    at runtime. NetCDF is entirely optional: `legacy` and `mpiio` need
    nothing beyond MPI.

---
