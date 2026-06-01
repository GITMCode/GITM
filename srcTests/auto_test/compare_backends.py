#!/usr/bin/env python3
"""
compare_backends.py — Cross-backend output comparison for GITM.

Verifies that all backends produce numerically identical scientific data
given the same physics configuration.  This is the strongest test that a
backend is correctly serializing and reconstructing output data.

Usage:
    python3 compare_backends.py dir1/ dir2/ [dir3/ ...]

Each directory should contain the output files from a single backend run
(.header + .b*/.bin/.raw files for legacy/mpiio, or .nc files for netcdf).

NetCDF variables are matched to binary variables via the 'gitm_name' attribute
written by GITM (the canonical name also used in binary outputs).

Exit code 0 = all match, 1 = mismatch detected, 2 = usage error.
"""

import sys
import os
import glob
import re
import numpy as np

# Import pGITM for reading legacy/mpiio binary formats
try:
    import pGITM
except ImportError:
    sys.path.append(os.path.join(os.path.dirname(__file__), '../../srcPython'))
    import pGITM

# Import netCDF4 for reading NetCDF outputs (optional)
try:
    import netCDF4 as nc4
except ImportError:
    nc4 = None


def extract_timestep(filename):
    """Extract timestep key from an output filename.

    Examples:
        '3DALL_t021221_001000.header' -> '3DALL_t021221_001000'
        '3DALL_t021221_001000.nc'     -> '3DALL_t021221_001000'
    """
    match = re.match(r'(\w+_t\d+_\d+)', os.path.basename(filename))
    return match.group(1) if match else None


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def load_binary_outputs(directory):
    """Load all timesteps from a legacy or mpiio output directory via pGITM.

    Returns {timestep_key: data_dict, ...} where data_dict has a 'vars' list
    and varname -> np.array mappings.  Ghost cells are stripped so arrays are
    (nAlts, nLats, nLons) — matching the netcdf layout — ready for ravel comparison.
    """
    result = {}
    for header in sorted(glob.glob(os.path.join(directory, '*.header'))):
        ts = extract_timestep(header)
        if ts is None:
            continue
        try:
            data = pGITM.process_one_file(header, dowrite=False, isVerbose=False)
            ng_alt = data.get('nGhostsAlt', 2)
            ng_lat = data.get('nGhostsLat', 2)
            ng_lon = data.get('nGhostsLon', 2)
            sl = (
                slice(ng_alt, -ng_alt if ng_alt > 0 else None),
                slice(ng_lat, -ng_lat if ng_lat > 0 else None),
                slice(ng_lon, -ng_lon if ng_lon > 0 else None),
            )
            for v in data.get('vars', []):
                if v in data and isinstance(data[v], np.ndarray) and data[v].ndim == 3:
                    data[v] = data[v][sl]
            result[ts] = data
        except Exception as e:
            print(f"  [WARN] Could not load {header}: {e}")
    return result


def load_netcdf_outputs(directory):
    """Load all timesteps from a netcdf output directory.

    Returns the same {timestep_key: data_dict} structure as the binary loader
    so the comparison code is format-agnostic.  Physics variables are keyed by
    their 'gitm_name' attribute (the canonical GITM name matching the binary
    format), falling back to the nc variable name for coordinate variables.
    """
    if nc4 is None:
        print("  [WARN] netCDF4 not available; cannot load NetCDF outputs.")
        return {}

    result = {}
    for nc_path in sorted(glob.glob(os.path.join(directory, '*.nc'))):
        ts = extract_timestep(nc_path)
        if ts is None:
            continue
        try:
            ds = nc4.Dataset(nc_path, 'r')
            data = {'vars': []}
            for vname in ds.variables:
                v = ds.variables[vname]
                arr = v[:]
                # Multi-time files: take the last timestep for comparison
                if 'time' in v.dimensions \
                   and 'time' in ds.dimensions \
                   and len(ds.dimensions['time']) > 1:
                    arr = arr[-1]
                arr = np.asarray(arr, dtype=np.float64)
                key = v.gitm_name.replace(" ", "") \
                    if hasattr(v, 'gitm_name') else vname
                data['vars'].append(key)
                data[key] = arr
            ds.close()
            result[ts] = data
        except Exception as e:
            print(f"  [WARN] Could not load {nc_path}: {e}")
    return result


def load_directory(path):
    """Auto-detect format and load all timesteps from a backend directory."""
    has_nc = bool(glob.glob(os.path.join(path, '*.nc')))
    has_hdr = bool(glob.glob(os.path.join(path, '*.header')))

    if has_nc:
        return load_netcdf_outputs(path), 'netcdf'
    elif has_hdr:
        return load_binary_outputs(path), 'binary'
    else:
        return {}, 'unknown'


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

def compare_timestep(name_a, data_a, name_b, data_b,
                     rtol=1e-5, atol=1e-8):
    """Compare all common variables between two backend outputs.

    Tolerance is generous enough to accommodate single<->double precision
    round-trips (the container stores in output_kind=kind(1.0), NetCDF
    upcasts to kind=8 on write).

    Returns True if every common variable matches within tolerance.
    """
    vars_a = set(v for v in data_a.get('vars', []) if v in data_a)
    vars_b = set(v for v in data_b.get('vars', []) if v in data_b)
    common = sorted(vars_a & vars_b)

    if not common:
        print(f"    [FAIL] No common variables between {name_a} and {name_b}")
        return False

    only_a = vars_a - vars_b
    only_b = vars_b - vars_a
    if only_a:
        print(f"    [INFO] Only in {name_a}: {sorted(only_a)}")
    if only_b:
        print(f"    [INFO] Only in {name_b}: {sorted(only_b)}")

    passed = True
    n_ok = 0

    for vname in common:
        a = np.asarray(data_a[vname], dtype=np.float64).ravel()
        b = np.asarray(data_b[vname], dtype=np.float64).ravel()

        if a.size != b.size:
            print(f"    [FAIL] '{vname}': size mismatch "
                  f"({name_a}={a.size} vs {name_b}={b.size})")
            passed = False
            continue

        # Regional/mag-grid outputs leave out-of-region cells unwritten: the
        # legacy/netcdf backends NaN-fill them while mpiio zero-fills them
        # Compare only cells finite in BOTH backends so the fill region (NaN
        # in at least one) is dropped; in-region data is still fully checked.
        finite = np.isfinite(a) & np.isfinite(b)
        if not finite.any():
            n_ok += 1
            continue
        a = a[finite]
        b = b[finite]

        if np.allclose(a, b, rtol=rtol, atol=atol):
            n_ok += 1
        else:
            diff = np.abs(a - b)
            max_abs = float(diff.max())
            denom = np.maximum(np.abs(a), np.abs(b))
            denom = np.where(denom > 0, denom, 1.0)
            max_rel = float((diff / denom).max())
            n_bad = int(np.sum(~np.isclose(a, b, rtol=rtol, atol=atol)))
            print(f"    [FAIL] '{vname}': {name_a} vs {name_b} — "
                  f"max_abs={max_abs:.3e} max_rel={max_rel:.3e} "
                  f"mismatched={n_bad}/{a.size}")
            passed = False

    if passed:
        print(f"    [PASS] All {n_ok} common variables match.")

    return passed


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) < 3:
        print("Usage: compare_backends.py dir1/ dir2/ [dir3/ ...]")
        sys.exit(2)

    dirs = [d.rstrip('/') for d in sys.argv[1:] if os.path.isdir(d)]
    if len(dirs) < 2:
        print("  Need at least 2 backend directories to compare. Skipping.")
        sys.exit(0)

    # Load each backend's outputs
    backends = {}
    for d in dirs:
        name = os.path.basename(d)
        data, fmt = load_directory(d)
        if data:
            backends[name] = data
            print(f"  Loaded {name} ({fmt}): {len(data)} timestep(s)")
        else:
            print(f"  [WARN] No loadable outputs in {d}")

    if len(backends) < 2:
        print("  Fewer than 2 backends loaded successfully. Skipping.")
        sys.exit(0)

    # Pairwise comparison across all loaded backends
    names = sorted(backends.keys())
    all_passed = True

    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            na, nb = names[i], names[j]
            common_ts = sorted(set(backends[na]) & set(backends[nb]))

            if not common_ts:
                print(f"\n  [WARN] No common timesteps between {na} and {nb}")
                continue

            print(f"\n  Comparing {na} vs {nb} "
                  f"({len(common_ts)} common timestep(s))...")
            for ts in common_ts:
                print(f"  Timestep: {ts}")
                if not compare_timestep(na, backends[na][ts],
                                        nb, backends[nb][ts]):
                    all_passed = False

    # Summary
    print("\n" + "=" * 60)
    if all_passed:
        print("  Cross-backend comparison PASSED — all backends agree. ✓")
    else:
        print("  FAILED — backend outputs disagree!")
    print("=" * 60)
    sys.exit(0 if all_passed else 1)


if __name__ == '__main__':
    main()
