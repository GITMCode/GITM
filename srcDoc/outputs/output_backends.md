# Output Backends

GITM can write outputs using three different backends. All three produce the same data
and are supported by `pGITM.py`. Which one you use depends on your system and workflow.

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
analysed. `post_process.py` handles this step automatically — see
[Postprocessing](../postprocessing.md) for details.

!!! warning "Lustre and parallel filesystems"
    On Lustre-based HPC systems, having hundreds of MPI ranks each create small files
    simultaneously can saturate the metadata server and slow the run significantly.
    If this is a concern, consider switching to the `mpiio` backend.

---

## mpiio

The MPI-IO backend writes all blocks into a single shared file using MPI collective
I/O[^collectiveio] (`MPI_File_open` / `MPI_File_write_at`). Each rank writes its block(s)
at a fixed byte offset[^byteoffset], so all ranks write simultaneously without
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

The NetCDF backend uses PnetCDF[^pnetcdf] to write a single self-describing file per
output type per timestep. Each MPI rank writes its interior grid slices[^interiorslices]
directly into the correct position in the global array, and ghost cells are stripped
automatically.

The result is a single `.nc` file per output type per timestep — no post-processing or
merging is required. Files can be read directly with xarray, NCO, PyITM, or any
NetCDF-aware tool.

### Conditional participation types

For output types where not all blocks write, the NetCDF backend automatically switches from
collective to **independent I/O mode** (`nfmpi_begin_indep_data` / `nfmpi_end_indep_data`).
This applies to:

- **Magnetic grid types** (`usesMagGrid=.true.`, e.g., 2DMEL) — only block 1 has the data
- **Regional types** (`isRegional=.true.`, e.g., 2DHME) — only blocks within the user region write

Independent mode allows only certain blocks/processors to write without causing MPI collective
operation deadlocks. The transition happens transparently during file open/close; no special
configuration is needed.

!!! note "Build requirements"
    PnetCDF must be installed and `pnetcdf-config` must be on your `PATH` when
    compiling. The build system detects it automatically and defines `-DHavePNetCDF`.
    Without this flag, the NetCDF module compiles as a stub and GITM will stop if
    `netcdf` is requested at runtime.

    NetCDF is entirely optional. The `legacy` and `mpiio` backends have no additional
    external dependencies beyond MPI.

---

[^collectiveio]: Collective I/O — an MPI operation where all ranks coordinate file writes
    simultaneously rather than independently, reducing metadata server load.

[^byteoffset]: Byte offset — the starting position (in bytes) within a file where data
    should be written. Each rank writes at a pre-computed offset so blocks do not overlap.

[^interiorslices]: Interior grid slices — the portion of each rank's grid excluding ghost
    cells (indices 1:nLons, 1:nLats, [1:nAlts for 3D]). This interior data is written to
    the global array at the correct position.

[^pnetcdf]: PnetCDF (Parallel NetCDF) — a parallel-capable version of NetCDF that allows
    multiple MPI ranks to write a single file simultaneously without serialization.
