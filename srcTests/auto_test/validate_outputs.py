#!/usr/bin/env python3
"""
validate_outputs.py — Physical sanity checks for GITM scientific outputs.

Checks absolute physical constraints: coordinate units, sign, non-zero where
required.  Cross-backend numerical agreement is handled by compare_backends.py.

Usage:
    ./validate_outputs.py <directory_or_file> [--verbose]
"""

import sys
import os
import glob
import argparse
import numpy as np

try:
    import pGITM
except ImportError:
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../srcPython')))
    import pGITM

try:
    import netCDF4 as nc
except ImportError:
    nc = None

verbose = False


def log_fail(label, detail):
    print(f"  [FAIL] {label}: {detail}")


def log_pass(label, detail):
    if verbose:
        print(f"  [PASS] {label}: {detail}")


def log_check(label, condition, detail):
    if condition:
        log_pass(label, detail)
    else:
        log_fail(label, detail)
    return condition


def validate_variable(name, vals):
    """Apply name-driven physical rules.  Returns True if all checks pass."""
    passed = True
    name_lower = name.lower()

    flat = np.asarray(vals, dtype=np.float64).ravel()
    flat = flat[np.isfinite(flat)]
    if len(flat) == 0:
        return True

    # --- Coordinate unit check: lat/lon must be in degrees, not radians ---
    if 'lat' in name_lower or 'lon' in name_lower:
        max_abs = float(np.abs(flat).max())
        passed = log_check(
            f"'{name}' in degrees (not radians)",
            max_abs > 3.15,
            f"max_abs={max_abs:.3f}"
        ) and passed

    # --- Latitude bounds [-130, 130] (accommodates ghost cells) ---
    if 'lat' in name_lower:
        vmin, vmax = float(flat.min()), float(flat.max())
        passed = log_check(
            f"'{name}' latitude bounds [-130, 130]",
            vmin >= -130.0 and vmax <= 130.0,
            f"min={vmin:.3f} max={vmax:.3f}"
        ) and passed

    # --- Longitude bounds [-370, 370] ---
    elif 'lon' in name_lower:
        vmin, vmax = float(flat.min()), float(flat.max())
        passed = log_check(
            f"'{name}' longitude bounds [-370, 370]",
            vmin >= -370.0 and vmax <= 370.0,
            f"min={vmin:.3f} max={vmax:.3f}"
        ) and passed

    # --- Altitude: exact name match only (avoids AltInt* variables) ---
    elif name_lower == 'altitude':
        vmin = float(flat.min())
        passed = log_check(
            f"'{name}' > 0",
            vmin > 0.0,
            f"min={vmin:.3f}"
        ) and passed

    # --- Temperature must be positive ---
    elif name_lower in ('temperature', 'etemperature', 'itemperature', 'tn') \
            or (name_lower.startswith('temp') and 'rate' not in name_lower):
        vmin = float(flat.min())
        passed = log_check(
            f"'{name}' > 0 K",
            vmin > 0.0,
            f"min={vmin:.3f}"
        ) and passed

    # --- Density/pressure: non-negative and not entirely zero ---
    elif any(t in name_lower for t in ('rho', 'density')) \
            or name_lower in ('press', 'pressure'):
        vmin = float(flat.min())
        passed = log_check(
            f"'{name}' >= 0",
            vmin >= 0.0,
            f"min={vmin:.3g}"
        ) and passed
        passed = log_check(
            f"'{name}' non-zero",
            int(np.count_nonzero(flat)) > 0,
            f"{int(np.count_nonzero(flat))} non-zero values"
        ) and passed

    elif verbose and np.count_nonzero(flat) == 0:
        print(f"  [INFO] '{name}' is all-zero (placeholder/inactive field).")

    return passed


def check_netcdf_file(nc_path):
    if verbose:
        print(f"\nChecking NetCDF: {os.path.basename(nc_path)}")
    if nc is None:
        log_fail("netCDF4 available", "package not installed")
        return False

    try:
        ds = nc.Dataset(nc_path, 'r')
    except Exception as e:
        log_fail(f"open {os.path.basename(nc_path)}", str(e))
        return False

    passed = True
    variables = list(ds.variables.keys())
    dims = ds.dimensions

    has_geo = "GeoLat" in variables and "GeoLon" in variables
    has_std = "Longitude" in variables and "Latitude" in variables
    if not log_check("Coordinate variables present", has_geo or has_std,
                     str([v for v in variables if v in ('GeoLat', 'GeoLon', 'Longitude', 'Latitude')])):
        passed = False

    is_append = 'time' in dims and len(dims['time']) > 1

    for varname in variables:
        var = ds.variables[varname]
        data = var[:]
        if is_append and 'time' in var.dimensions:
            t_last = len(dims['time']) - 1
            if not validate_variable(varname, data[t_last]):
                passed = False
        else:
            if not validate_variable(varname, data):
                passed = False

    ds.close()
    return passed


def check_header_file(header_path):
    if verbose:
        print(f"\nChecking binary: {os.path.basename(header_path)}")
    try:
        data = pGITM.process_one_file(header_path, dowrite=False, isVerbose=False)
    except Exception as e:
        log_fail(f"load {os.path.basename(header_path)}", str(e))
        return False

    passed = True
    nX = data.get('nLonsTotal', 0)
    nY = data.get('nLatsTotal', 0)
    nZ = data.get('nAltsTotal', 0)
    if not log_check("Grid dims > 0", nX > 0 and nY > 0 and nZ > 0,
                     f"nLons={nX} nLats={nY} nAlts={nZ}"):
        passed = False

    for varname in data.get('vars', []):
        if varname not in data:
            log_fail(varname, "not found in processed data")
            passed = False
            continue
        if not validate_variable(varname, data[varname]):
            passed = False

    return passed


def main():
    global verbose

    parser = argparse.ArgumentParser()
    parser.add_argument("target", help="Directory or single output file")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Print PASS lines (default: failures only)")
    args = parser.parse_args()
    verbose = args.verbose

    if not os.path.exists(args.target):
        print(f"ERROR: {args.target} does not exist")
        sys.exit(2)

    files_to_check = []
    if os.path.isdir(args.target):
        files_to_check.extend(sorted(glob.glob(os.path.join(args.target, '*.header'))))
        files_to_check.extend(sorted(glob.glob(os.path.join(args.target, '*.nc'))))
    else:
        files_to_check.append(args.target)

    if not files_to_check:
        print("  [FAIL] No output files (*.nc or *.header) found.")
        sys.exit(1)

    overall = True
    n_files = len(files_to_check)
    n_failed = 0

    for f in files_to_check:
        ok = False
        if f.endswith('.header'):
            ok = check_header_file(f)
        elif f.endswith('.nc'):
            ok = check_netcdf_file(f)
        if not ok:
            if not verbose:
                print(f"  (from {os.path.basename(f)})")
            n_failed += 1
            overall = False

    print("\n" + "=" * 60)
    if overall:
        print(f"  All {n_files} output file(s) passed validation. ✓")
    else:
        print(f"  {n_failed}/{n_files} file(s) FAILED validation.")
    print("=" * 60)
    sys.exit(0 if overall else 1)


if __name__ == "__main__":
    main()
