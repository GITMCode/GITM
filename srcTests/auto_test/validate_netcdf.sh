#!/bin/bash
# validate_netcdf.sh
#
# Validate container-based output backends (netcdf and mpiio).
# Runs GITM once per provided UAM.in file, saves outputs to a timestamped
# directory, and optionally checks structure/data.
#
# Usage (from srcTests/auto_test/run/ or any rundir with GITM.exe):
#   ../validate_netcdf.sh [options] UAM.in.file1 [UAM.in.file2 ...]
#
# Options:
#   -c, --check       Check outputs after each run.
#   -n, --np N        MPI rank count (default: 4).
#   --oversubscribe   Pass --oversubscribe to mpirun.
#   -h, --help        Show this help.
#
# Defaults (no UAM files given):
#   Runs tests 12 (mpiio), 13 (netcdf single-time), and 14 (netcdf append).
#
# Output directories are created under UA/data/ in the run dir:
#   UA/data/<uam_basename>_v<N>/
#
# Checks performed (with --check):
#   mpiio backend (.bin files):
#     - 2DTEC, 2DGEL: header nVars, coord ranges in degrees, non-zero physics.
#     - 3DNEU, 3DALL, 3DION, 3DCHM: header nVars, 3D dims, coord ranges.
#     - 2DMEL: binary file presence (structure validated separately).
#   netcdf backend (.nc files):
#     - ncdump -h on every 2DMEL*.nc: dims, coord vars, and all var names.
#     - GeoLat/GeoLon in degrees (not radians), non-zero values.
#     - For append files: time dimension > 1 record.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_UAM_12="${SCRIPT_DIR}/UAM.in.12.BACKEND_TEST.mpiio.test"
DEFAULT_UAM_13="${SCRIPT_DIR}/UAM.in.13.BACKEND_TEST.netcdf.test"
DEFAULT_UAM_14="${SCRIPT_DIR}/UAM.in.14.BACKEND_TEST.netcdf.appendSingle.test"

NP=4
OVERSUBSCRIBE=""
DO_CHECK=false
UAM_FILES=()

usage() {
  sed -n '2,/^$/p' "$0" | grep '^#' | sed 's/^# \?//'
  exit 0
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)       usage ;;
    -c|--check)      DO_CHECK=true; shift ;;
    -n|--np)         NP="$2"; shift 2 ;;
    --oversubscribe) OVERSUBSCRIBE="--oversubscribe"; shift ;;
    -*)              echo "Unknown option: $1"; usage ;;
    *)               UAM_FILES+=("$1"); shift ;;
  esac
done

if [[ ${#UAM_FILES[@]} -eq 0 ]]; then
  UAM_FILES=("$DEFAULT_UAM_12" "$DEFAULT_UAM_13" "$DEFAULT_UAM_14")
fi

# Must be run from a rundir containing GITM.exe
if [[ ! -f ./GITM.exe ]]; then
  echo "ERROR: GITM.exe not found in $(pwd)"
  echo "Run this script from a rundir (e.g. srcTests/auto_test/run/ or run/)."
  exit 1
fi

# -------------------------------------------------------------------------
# Helper: next available version suffix (avoids clobbering old results)
# -------------------------------------------------------------------------
next_outdir() {
  local base="UA/data/$1"
  local n=1
  while [[ -d "${base}_v${n}" ]]; do
    n=$((n+1))
  done
  echo "${base}_v${n}"
}

# -------------------------------------------------------------------------
# Helper: describe a UAM file's output mode from filename
# -------------------------------------------------------------------------
describe_uam() {
  local f
  f="$(basename "$1")"
  if [[ "$f" == *"appendSingle"* || "$f" == *"append"* ]]; then
    echo "netcdf-append"
  elif [[ "$f" == *"netcdf"* ]]; then
    echo "netcdf-singletime"
  elif [[ "$f" == *"mpiio"* ]]; then
    echo "mpiio"
  else
    echo "$(basename "$1")"
  fi
}

# -------------------------------------------------------------------------
# Helper: next-version output dir name from UAM filename
# -------------------------------------------------------------------------
outdir_name_for() {
  local f
  f="$(basename "$1" .test)"
  f="${f#UAM.in.}"      # strip leading UAM.in.
  f="${f//[. ]/_}"      # dots/spaces → underscores
  next_outdir "$f"
}

# -------------------------------------------------------------------------
# run_one: run GITM with one UAM file, save outputs, optionally check
# -------------------------------------------------------------------------
run_one() {
  local uam_path="$1"
  local label
  label="$(describe_uam "$uam_path")"
  local outdir
  outdir="$(outdir_name_for "$uam_path")"

  if [[ ! -f "$uam_path" ]]; then
    echo "ERROR: UAM file not found: $uam_path"
    return 1
  fi

  echo ""
  echo "================================================================"
  echo " Running: $label"
  echo " UAM:     $uam_path"
  echo " Output:  $outdir"
  echo "================================================================"

  # For append mode, we run GITM twice so the file accumulates two records.
  local is_append=false
  if [[ "$label" == "netcdf-append" ]]; then
    is_append=true
  fi

  ln -sf "$uam_path" UAM.in
  rm -f GITM.DONE UA/data/*.nc UA/data/*.bin UA/data/*.header

  echo ">> Pass 1 ..."
  mpirun -np "$NP" $OVERSUBSCRIBE ./GITM.exe
  if [[ ! -f GITM.DONE ]]; then
    echo "ERROR: GITM did not finish (no GITM.DONE) on pass 1."
    return 1
  fi
  rm -f GITM.DONE

  if $is_append; then
    echo ">> Pass 2 (append mode: second run to accumulate time records) ..."
    mpirun -np "$NP" $OVERSUBSCRIBE ./GITM.exe
    if [[ ! -f GITM.DONE ]]; then
      echo "ERROR: GITM did not finish (no GITM.DONE) on pass 2."
      return 1
    fi
    rm -f GITM.DONE
  fi

  # Collect all container output files (.nc, .bin, .header)
  mkdir -p "$outdir"
  cp UA/data/*.nc  "$outdir/" 2>/dev/null || true
  cp UA/data/*.bin "$outdir/" 2>/dev/null || true
  cp UA/data/*.header "$outdir/" 2>/dev/null || true

  local nc_count bin_count
  nc_count=$(ls "$outdir"/*.nc  2>/dev/null | wc -l)
  bin_count=$(ls "$outdir"/*.bin 2>/dev/null | wc -l)
  echo ""
  echo "Saved to $outdir: ${nc_count} .nc file(s), ${bin_count} .bin file(s)"

  if $DO_CHECK; then
    check_outputs "$outdir" "$is_append" "$label"
  fi
}

# -------------------------------------------------------------------------
# check_mpiio_outputs: validate .bin/.header files for mpiio backend
# -------------------------------------------------------------------------
check_mpiio_outputs() {
  local outdir="$1"
  local PYTHON=/home/aaron/miniconda3/envs/py3/bin/python
  local CHECK_SCRIPT
  CHECK_SCRIPT="$(dirname "${BASH_SOURCE[0]}")/check_geo2d_bin.py"
  local pass=true

  echo ""
  echo "--- mpiio checks: $outdir ---"

  # 2DMEL: just confirm .bin exists (full validation done in the netcdf tests)
  local mel_bins
  mel_bins=$(ls "$outdir"/2DMEL*.bin 2>/dev/null | wc -l)
  echo ""
  if [[ "$mel_bins" -gt 0 ]]; then
    echo "  [PASS] 2DMEL .bin files present: $mel_bins"
  else
    echo "  [WARN] No 2DMEL*.bin files found (may not be in output list)"
  fi

  # 2DTEC, 2DGEL, 3DNEU, 3DALL, 3DION, 3DCHM: full value checks via Python
  local any_container=false
  for ctype in 2DTEC 2DGEL 3DNEU 3DALL 3DION 3DCHM; do
    local count
    count=$(ls "$outdir"/${ctype}_*.bin 2>/dev/null | wc -l)
    if [[ "$count" -gt 0 ]]; then
      any_container=true
      echo ""
      echo "  [${ctype}] ${count} .bin file(s) found — running value checks"
      if "$PYTHON" "$CHECK_SCRIPT" "$outdir" --type "$ctype"; then
        : # pass
      else
        pass=false
      fi
    else
      echo ""
      echo "  [WARN] No ${ctype}_*.bin files found in $outdir"
    fi
  done

  if ! $any_container; then
    echo ""
    echo "  FAIL: no container .bin files found — container output may not be working"
    pass=false
  fi

  echo ""
  if $pass; then
    echo "  >>> All mpiio checks PASSED for $outdir <<<"
  else
    echo "  >>> One or more mpiio checks FAILED for $outdir <<<"
    return 1
  fi
}

# -------------------------------------------------------------------------
# check_outputs: dispatch to mpiio or netcdf checks based on what's in outdir
# -------------------------------------------------------------------------
check_outputs() {
  local outdir="$1"
  local is_append="$2"
  local label="${3:-}"

  echo ""
  echo "--- Checking $outdir ---"

  # Dispatch based on backend: mpiio produces .bin; netcdf produces .nc.
  local has_bin has_nc
  has_bin=$(ls "$outdir"/*.bin 2>/dev/null | wc -l)
  has_nc=$(ls "$outdir"/*.nc  2>/dev/null | wc -l)

  if [[ "$has_bin" -gt 0 ]]; then
    check_mpiio_outputs "$outdir"
    return $?
  fi

  if [[ "$has_nc" -eq 0 ]]; then
    echo "  FAIL: no .nc or .bin files found in $outdir"
    return 1
  fi

  # --- netcdf path: existing 2DMEL checks ---
  local pass=true
  local PYTHON=/home/aaron/miniconda3/envs/py3/bin/python
  local CHECK_SCRIPT
  CHECK_SCRIPT="$(dirname "${BASH_SOURCE[0]}")/check_2dmel_nc.py"

  local any_nc
  any_nc=$(ls "$outdir"/2DMEL*.nc 2>/dev/null | head -1)
  if [[ -z "$any_nc" ]]; then
    echo "  FAIL: no 2DMEL*.nc files found in $outdir"
    return 1
  fi

  for nc_file in "$outdir"/2DMEL*.nc; do
    echo ""
    echo "  File: $nc_file"

    # --- Structure check: ncdump -h ---
    if command -v ncdump &>/dev/null; then
      echo ""
      echo "  [ncdump -h]"
      ncdump -h "$nc_file" | sed 's/^/    /'
    else
      echo "  WARNING: ncdump not found, skipping structure check."
    fi

    # --- Value checks via Python/xarray ---
    echo ""
    echo "  [Python/xarray value checks]"
    local append_flag=""
    if [[ "$is_append" == "true" ]]; then
      append_flag="--append"
    fi
    if "$PYTHON" "$CHECK_SCRIPT" "$nc_file" $append_flag; then
      : # pass — check_2dmel_nc.py already printed per-check results
    else
      pass=false
    fi

    # --- Append mode time-dimension check (belt-and-suspenders) ---
    if [[ "$is_append" == "true" ]] && command -v ncdump &>/dev/null; then
      echo ""
      echo "  [ncdump time dimension cross-check]"
      local time_count
      time_count=$(ncdump -h "$nc_file" | grep -E 'time\s*=\s*UNLIMITED' | \
        grep -oP '\(\d+ currently\)' | grep -oP '\d+' || echo "0")
      printf "    time records: %s\n" "$time_count"
      if [[ -n "$time_count" && "$time_count" -gt 1 ]]; then
        echo "    PASS: multiple time records confirmed"
      else
        echo "    FAIL: expected >1 time records for append mode, got: $time_count"
        pass=false
      fi
    fi
  done

  echo ""
  if $pass; then
    echo "  >>> All checks PASSED for $outdir <<<"
  else
    echo "  >>> One or more checks FAILED for $outdir <<<"
    return 1
  fi
}

# -------------------------------------------------------------------------
# Main loop
# -------------------------------------------------------------------------
overall_pass=true
for uam in "${UAM_FILES[@]}"; do
  if ! run_one "$uam"; then
    overall_pass=false
  fi
done

echo ""
echo "================================================================"
if $overall_pass; then
  echo " All runs completed successfully."
else
  echo " One or more runs FAILED."
  exit 1
fi
if ! $DO_CHECK; then
  echo " Re-run with --check to validate output structure and values."
  echo " (mpiio: checks 2D+3D container binaries; netcdf: checks 2DMEL structure/values)"
fi
echo "================================================================"
