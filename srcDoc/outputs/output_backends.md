# Output Backends

GITM can write outputs using three different backends. All three produce the
same  data. They differ in how files are laid out on disk and in
the few limits listed in the table below.
Which one you "should" use depends on your system and workflow.

To select a backend, add this block to `UAM.in`:

```
#OUTPUTBACKEND
mpiio
```

Without this block, GITM uses the `legacy` backend, for backwards compatibility
with existing runs and tools.

| | `legacy` | `mpiio` | `netcdf` |
|---|---|---|---|
| Files for each output | `.header` + one `.b####` per MPI rank | `.header` + one raw `.bin` | one `.nc` |
| [Post-processing](../postprocessing.md) | required | required | none |
| `1D*`/`0D*` output types | yes | yes | no, GITM stops |
| Needs to compile | MPI | MPI | MPI and PnetCDF |

Satellite (`#SATELLITES`) files are written in the legacy format whatever the
backend.

[TOC]

---

## legacy

The legacy backend is the original GITM output format. Each MPI rank writes its own
block file independently, producing:

- One ASCII `.header` file (written by rank 0)
- One binary `.b####` file per MPI rank (e.g., `.b0001`, `.b0002`, ...)

These per-block files must be merged into a single `.bin` file before they can be
analyzed. `post_process.py` handles this step automatically. See
[Post-processing](../postprocessing.md) for details.

!!! warning "Lustre and parallel filesystems"
    On Lustre-based HPC systems, having hundreds of MPI ranks each create small files
    simultaneously can saturate the metadata server and slow the run significantly.
    If this is a concern, consider switching to the `mpiio` backend.

---

## mpiio

The MPI-IO backend writes all blocks into a single shared file using MPI collective
I/O (`MPI_File_open` / `MPI_File_write_at_all`). Each rank writes its block(s)
at a fixed byte offset, so all ranks write simultaneously without
serialization. (The magnetic-grid and regional types, where only some ranks
have data, are written with independent calls instead.) The result is:

- One ASCII `.header` file (same format as legacy, with an added `FORMAT/mpiio` line)
- One `.bin` file containing all blocks in a single raw-binary file (no Fortran
  record markers)

This raw `.bin` is not the layout the analysis tools read, so `mpiio` output
still needs [post-processing](../postprocessing.md). `pGITM.py` detects the
`FORMAT/mpiio` marker in the header and rewrites the file as a standard `.bin`,
identical to post-processed `legacy` output (with `-norm` the original is kept
as `.raw`). What `mpiio` saves is the per-rank file count, not the
post-processing step.

This backend is preferred on Lustre-based systems: it produces far fewer
files per timestep and benefits from large sequential writes rather than many small
independent ones.

A `0` value in a regional (`2DHME`/`3DHME`) `mpiio` file is expected for cells
outside the region: the file has a slot for every block, and the blocks that
do not write leave theirs zero-filled.

---

## netcdf

The `netcdf` backend uses PnetCDF to write a single self-describing file per
output type per timestep. Each MPI rank writes its interior (non-ghost) grid slices
directly into the correct position in the global array, and ghost cells are stripped
automatically.

The result is a single `.nc` file per output type per timestep. No post-processing or
merging is required. Files can be read directly with xarray, NCO, PyITM, or any
NetCDF-aware tool.

!!! note "Build requirements"
    PnetCDF must be installed and `pnetcdf-config` on `PATH` when compiling;
    the build system detects it and defines `-DHavePNetCDF`. Without it the
    NetCDF module compiles as a stub and GITM stops if `netcdf` is requested
    at runtime. NetCDF is entirely optional: `legacy` and `mpiio` need
    nothing beyond MPI.

Two limits:

- The `netcdf` backend does not support the probe-style `1D*`/`0D*` types
  (including 3D types renamed to `1D*` by a 1D `nLons=nLats=1` run) and
  stops with an error naming the type. Use `legacy` or `mpiio` for those.
- The `.nc` files this backend writes are a different dialect from the ones
  `pGITM.py` produces when converting legacy/mpiio output (see the warning
  under [`-nc`](../postprocessing.md#nc)). Don't mix the two in one analysis.

### Multi-time files, run names, and restarts

By default the `netcdf` backend creates one file per output type per timestep
(e.g., `3DALL_t021221_000000.nc`, `3DALL_t021221_000500.nc`, ...). `#NETCDFOPTIONS`
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

`daily` and `single` files use an unlimited time dimension. Data variables are
shaped `(time, [z,] lat, lon)` in C/numpy/xarray order (time is the last axis
in Fortran order). On a magnetic grid the horizontal dimensions are named
`mlat`/`mlon` instead, and carry plain `degrees`. xarray reads the files
natively:

```python
import xarray as xr
ds = xr.open_dataset("3DALL.nc", decode_times=True)
```

`netCDFFileName`, when set, is inserted into every netcdf filename after the
type code (`3DALL_storm1.nc`, `3DALL_storm1_t021221_000000.nc`).
This lets several scenario runs share one directory without colliding.

Restarts in `single`/`daily` mode depend on whether the run is named:

- **Unnamed**: the restarted run continues the previous leg's file, so a
  run split across restarts still ends up with one time series per file
  (one per day for `daily`). The first output of the new leg overwrites the
  last record if it is for the same time (a restart can re-emit the epoch
  it starts from), and is appended after it otherwise. GITM stops instead
  if the file has records *after* that time (an earlier run went past the
  restart point, and netcdf cannot remove those records), or if the file's grid
  or variables do not match the run. Move the file aside to start a new
  one.
- **Named**: a `restart_yymmdd_hhmmss` stamp (the restart's simulation
  time) is appended to the name. Each leg writes its own file and prior
  legs are never reopened, so a crash can only affect the current leg.
  Consecutive legs can share the restart epoch, so analysis tools should
  drop duplicate times, keeping the later leg.

!!! warning "Overwrite guard"
    In `single`/`daily` mode, a run with a run name refuses to overwrite an
    existing target file and stops with instructions (move/delete the file,
    change `netCDFFileName`, or use a fresh run directory). This catches
    re-running the same leg and branching a scenario in place. An unnamed,
    non-restart run keeps the historical behavior of overwriting. `none`
    mode and the `legacy`/`mpiio` backends are unguarded: colliding epochs
    are recomputed and overwritten in place. (The recompute is not
    guaranteed bit-identical: restarted runs are known to diverge slightly,
    likely via the dynamo / high-latitude electric fields.)

`#NETCDFOPTIONS` only affects the `netcdf` backend; it is ignored for
`legacy` and `mpiio`.

---
