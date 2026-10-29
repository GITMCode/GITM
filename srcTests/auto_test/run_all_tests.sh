#!/bin/bash

# Exit immediately if an unhandled pipeline/command fails
set -e

get_help(){

  printf "
------------------------------------------------------------------------------------

> This script compiles GITM and runs every test case in srcTests/auto_test/
  matching UAM.*.test (alphabetically).
  - The same tests run automatically on GitHub.
  - run/ is deleted and recreated on every invocation; copy anything you want
    to keep out of it first.

> When adding a new feature, it is recommended to create a test.
  - First create a UAM.in file which uses the new feature.
  - Then save a reference solution for it with:
    > ./run_all_tests.sh --save_solution -o [your uam file]
  - Name the file UAM.in.##.[a_few_words].test (no spaces). Numbers need not be
    sequential, but sequential numbers make outputs easier to compare.

> Output backend tests have BACKEND_TEST in the name, plus the backend
  (UAM.in.##.BACKEND_TEST.[legacy|mpiio|netcdf].*.test). On top of the log diff,
  these also run:
  - validate_outputs.py on the run's outputs (physical sanity checks).
  - compare_backends.py across backends, after all tests have run.
  - For appendSingle tests, a check that each netcdf file accumulated time
    records. A RESTART leg (e.g. 14.1) reopens the files its first leg
    (14.0) wrote, so run the pair together, not the RESTART leg alone.

  netcdf failures are reported as optional, but still fail the run.

------------------------------------------------------------------------------------
Usage:
------

> First, from root of repository, install GITM with './Config.pl [...]'
> Next, 'cd srcTests/auto_test/'

- To check that the code compiles and all tests pass:
    > ./run_all_tests.sh

- To re-run a single test after a failure:
    > ./run_all_tests.sh --only UAM.in.00.default.test

- To run only the output backend tests:
    > ./run_all_tests.sh --backends_only

- To force a clean re-compile:
    > ./run_all_tests.sh -c

------------------------------------------------------------------------------------
Arguments:
----------

        [none]                    Compile, then run all tests with defaults.
        -h, --help                Show this information.
        -c, --clean               Run 'make clean' before compiling.
        -u, --uninstall           Run 'make distclean', then reinstall with
                                    './Config.pl -install -earth -compiler=gfortran -debug'.
                                    Useful when dependencies change, or when a lot of
                                    changes have been made (by you or git).
        -o, --only FILE           Only run this UAM test file.
        --backends_only           Only run the BACKEND_TEST tests.
        --nocompare               Run the tests without diffing logs against the
                                    reference solutions or comparing backends.
                                    validate_outputs.py and the append check still run.
        --save_solution           Rewrite (or create) the reference solution(s).
                                    - Similar to the SWMF's '-BLESS'.
                                    - Pauses briefly first, so you can Ctrl+C if
                                      this was set by mistake.
                                    - Skips all output checks for that run.
                                    - If a solution was overwritten by mistake, use
                                      'git restore' to put it back.
        --oversubscribe           Pass '--oversubscribe' to mpirun.
                                    Needed on GitHub's CI servers, but not for *most*
                                    users. MPI libraries implement this differently,
                                    so it is off by default.

"
}

get_backend_from_filename(){
  # Extract backend from filename: UAM.in.##.description.BACKEND.anything.test
  # Returns: legacy, mpiio, netcdf, or "unknown" if not found
  local filename="$1"

  if [[ "$filename" =~ \.(legacy|mpiio|netcdf)\. ]]; then
    echo "${BASH_REMATCH[1]}"
  else
    echo "unknown"
  fi
}

is_optional_backend(){
  # Returns 0 (true) if backend is optional, 1 (false) if required
  local backend="$1"

  if [[ "$backend" == "netcdf" ]]; then
    return 0  # netcdf is optional
  else
    return 1  # legacy and mpiio are required
  fi
}

warnsavesolution(){
  # refactored to clean up arg parsing
  # This is only called if --savesolutions is set

  printf "
!!!!!!!!!!!!!!!!    WARNING    !!!!!!!!!!!!!!!!!

Using --save_solution will overwrite the reference logfiles

Only use this if you know what you are doing.

If you just created a new test, run with -o/--only
to save a solution for the test you created.

Cancel with Ctrl+C

"
  sleep 4

  return

}



checkoutputs(){
  # refactored for cleanliness
  # Returns 0 on success, 1 on failure

  # 1. Handle solution saving for global log files
  if [ $do_save = true ]; then
    cp UA/data/log*.dat ../ref_solns/log.$test_uam
    echo
    echo
    echo " Created ref solution: log." $test_uam
    return 0
  fi

  # 2. Diff the global log files against the reference
  if [ $do_compare = true ]; then
    echo
    # Capture diff status without triggering set -e
    local diff_status=0
    ../../../share/Scripts/DiffNum.py -r=5e-5 -a=1e-2 -t ../ref_solns/log.$test_uam UA/data/log*.dat || diff_status=$?

    if [ $diff_status != 0 ]; then
      echo
      echo " ============    ERROR!!!!!    ============"
      echo "  Output differs from reference solution."
      echo "  Something has gone terribly wrong!!"
      rm -f GITM.DONE
      return 1
    fi
  fi

  # 3. Physical value/sanity checks (BACKEND_TEST runs only)
  if [[ "$test_uam" == *"BACKEND_TEST"* ]] && [ -x ../validate_outputs.py ]; then
    echo ">> Running output validation checks..."
    local python_status=0
    ../validate_outputs.py "UA/data/" || python_status=$?
    if [ $python_status != 0 ]; then
      rm -f GITM.DONE
      return 1
    fi
  fi

  # 4. Single-file netcdf append: every file must have accumulated records
  if [[ "$test_uam" == *"appendSingle"* ]]; then
    if ! check_append_records; then
      rm -f GITM.DONE
      return 1
    fi
  fi

  return 0
}

nc_time_records(){
  # Print the length of the time dimension of netcdf file $1
  python3 -c "import netCDF4, sys; print(len(netCDF4.Dataset(sys.argv[1]).dimensions['time']))" "$1"
}

check_append_records(){
  # The first leg (14.0) writes twice in one run, so each file needs >= 2
  # records. The RESTART leg (14.1) reopens the files 14.0 left in
  # $append_carry, so each must end up with more records than it started with.
  # Returns 0 on success, 1 on failure
  echo ">> Checking netcdf append time records..."
  local ok=true
  local nc_files=(UA/data/*.nc)
  if [ ! -f "${nc_files[0]}" ]; then
    echo "  [FAIL] no .nc files in UA/data/"
    return 1
  fi
  for f in "${nc_files[@]}"; do
    local name=$(basename "$f")
    local n=$(nc_time_records "$f")
    if [[ "$test_uam" == *"RESTART"* ]]; then
      if [ ! -f "$append_carry/$name" ]; then
        echo "  [FAIL] $name: no file carried over from the first leg"
        ok=false
        continue
      fi
      local n0=$(nc_time_records "$append_carry/$name")
      if [ "$n" -gt "$n0" ]; then
        echo "  [PASS] $name: $n0 -> $n time records"
      else
        echo "  [FAIL] $name: $n0 -> $n time records (restart did not append)"
        ok=false
      fi
    else
      if [ "$n" -ge 2 ]; then
        echo "  [PASS] $name: $n time records"
      else
        echo "  [FAIL] $name: $n time record(s), expected >= 2"
        ok=false
      fi
    fi
  done
  $ok
}

run_a_test(){
  # Run a single test and return status:
  #   0 = success
  #   1 = failure (required backend)
  #   2 = failure (optional backend)

  printf "\n\n>> Testing with $test_uam ...\n"
  # Copy UAM
  ln -sf $test_uam UAM.in
  
  # Remove output files left over from previous runs
  rm -f GITM.DONE UA/data/*.nc UA/data/*.bin UA/data/*.header UA/data/*.dat UA/data/*.b[0-9][0-9][0-9][0-9]

  # An append RESTART leg must reopen the files its first leg wrote. A first
  # leg clears any stale carry-over, so a failed 14.0 can't feed an old file
  # to 14.1.
  if [[ "$test_uam" == *"appendSingle"* ]]; then
    if [[ "$test_uam" != *"RESTART"* ]]; then
      rm -rf "$append_carry"
    elif [ -d "$append_carry" ]; then
      cp "$append_carry"/*.nc UA/data/
    else
      echo " WARNING: no first-leg output in $append_carry; the append check will fail"
    fi
  fi

  # Run GITM, capture status without triggering set -e
  local mpi_status=0
  mpirun -np 4 $oversubscribe ./GITM.exe || mpi_status=$?

  # this will either save, or diff, the output log files.
  local checkout_status=0
  if [ -f GITM.DONE ] && [ $mpi_status = 0 ]; then
    checkoutputs || checkout_status=$?
  else
    checkout_status=1
  fi

  if [ -f GITM.DONE ] && [ $checkout_status = 0 ]; then
      printf "\n\n>>> $test_uam ran successfully! <<< \n\n"
      mv $test_uam $test_uam.success
      mv UA/data/log*.dat UA/data/log_$test_uam.success
      rm -f GITM.DONE

      # Keep a first append leg's files for its RESTART leg; drop them after
      if [[ "$test_uam" == *"appendSingle"* ]]; then
        if [[ "$test_uam" == *"RESTART"* ]]; then
          rm -rf "$append_carry"
        else
          rm -rf "$append_carry"
          mkdir -p "$append_carry"
          cp UA/data/*.nc "$append_carry"/
        fi
      fi

      # Passed: clean up
      if [ $do_save = false ]; then
        # For backend tests, stage outputs for cross-backend comparison.
        # Only the first run per backend is kept (e.g. netcdf test 13 is
        # staged, the append variants 14.x are not).
        local backend=$(get_backend_from_filename "$test_uam")
        local stage_dir="UA/data/backend_compare/${backend}"
        if [[ "$test_uam" == *"BACKEND_TEST"* ]] && [[ "$backend" != "unknown" ]] && [ ! -d "$stage_dir" ]; then
          mkdir -p "$stage_dir"
          cp UA/data/*.nc UA/data/*.bin UA/data/*.header UA/data/*.b[0-9][0-9][0-9][0-9] "$stage_dir/" 2>/dev/null || true
        fi
        rm -f UA/data/*.nc UA/data/*.bin UA/data/*.header UA/data/*.b[0-9][0-9][0-9][0-9]
      fi
      return 0
  else
      printf "\n\n>>> $test_uam   UNSUCCESSFUL! <<< \n\n"
      
      # Failed: keep the outputs for debugging
      mkdir -p UA/data/failed_outputs_$test_uam
      mv UA/data/*.nc UA/data/*.bin UA/data/*.header UA/data/*.dat UA/data/*.b[0-9][0-9][0-9][0-9] UA/data/failed_outputs_$test_uam/ 2>/dev/null || true
      rm -f GITM.DONE
      # paths relative to srcTests/auto_test/, where this script is run from
      if [[ $do_compare = true && -f ../ref_solns/log.$test_uam ]]; then
        printf "  To compare (from srcTests/auto_test/):\n"
        printf "    ../../share/Scripts/DiffNum.py -r=5e-5 -a=1e-2 ref_solns/log.%s run/UA/data/failed_outputs_%s/log*.dat\n\n" \
          "$test_uam" "$test_uam"
      fi

      # Determine if this is optional or required backend
      local backend=$(get_backend_from_filename "$test_uam")
      if is_optional_backend "$backend"; then
        return 2  # optional backend failure
      else
        return 1  # required backend failure
      fi
  fi

}

print_summary(){
  # Print test results summary
  local passed_count=${#passed[@]}
  local failed_required_count=${#failed_required[@]}
  local failed_optional_count=${#failed_optional[@]}
  local total_count=$((passed_count + failed_required_count + failed_optional_count))

  echo ""
  echo "============================================================"
  echo "TEST SUMMARY"
  echo "============================================================"
  echo "Total tests: $total_count"
  echo "Passed: $passed_count"

  if [ $failed_required_count -gt 0 ]; then
    echo "Failed (required): $failed_required_count"
    for test in "${failed_required[@]}"; do
      echo "   - $test"
    done
  fi

  if [ $failed_optional_count -gt 0 ]; then
    echo "Failed (optional): $failed_optional_count"
    for test in "${failed_optional[@]}"; do
      echo "   - $test (netcdf - optional feature)"
    done
  fi
  echo "============================================================"
  echo ""
}

do_tests(){
    # Configure & compile
    cd ../../
    
    if [ $clean = true ]; then
      make clean
    elif [ $distclean = true ]; then
      make distclean
    fi

    if [ $config = true ]; then
        ./Config.pl -install -earth -compiler=gfortran -debug
    fi

    # Capture make status without triggering set -e
    make -j || {
      echo
      echo "Could not compile!"
      echo "The tests have failed. Exiting."
      exit 1
    }

    # Since we can't tell when (or how) srcTest/auto_test/run was made, replace it
    # 'make rundir' defaults to creating 'run', if we set RUNDIR, it creates that instead.
    rm -rf srcTests/auto_test/run
    make rundir RUNDIR=srcTests/auto_test/run

    # Copy the test files into run/
    cd srcTests/auto_test/
    rm -f run/UAM*
    cp UAM.*.test run/

    # begin running:
    cd run/

    # Initialize result tracking arrays
    passed=()
    failed_required=()
    failed_optional=()

    if [ $onlyone = false ]; then
      local test_glob="UAM.*.test"
      if [ $backends_only = true ]; then test_glob="UAM.*BACKEND_TEST*.test"; fi
      for test_uam in $test_glob; do
          run_a_test
          local status=$?

          if [ $status = 0 ]; then
            passed+=("$test_uam")
          elif [ $status = 1 ]; then
            failed_required+=("$test_uam")
          elif [ $status = 2 ]; then
            failed_optional+=("$test_uam")
          fi
      done

      # Cross-backend data comparison (runs after all individual backend tests)
      if [ $do_save = false ] && [ $do_compare = true ]; then
        local compare_dir="UA/data/backend_compare"
        if [ -d "$compare_dir" ]; then
          local ndirs=$(find "$compare_dir" -mindepth 1 -maxdepth 1 -type d | wc -l)
          if [ "$ndirs" -ge 2 ]; then
            echo ""
            echo ">> Cross-backend data comparison..."
            local compare_status=0
            python3 ../compare_backends.py "$compare_dir"/*/ || compare_status=$?
            if [ $compare_status != 0 ]; then
              failed_required+=("CROSS_BACKEND_COMPARISON")
              # Keep staged per-backend outputs so the failure can be
              # inspected / compare_backends.py re-run by hand.
              echo "   Outputs kept for inspection at: $compare_dir"
            else
              passed+=("CROSS_BACKEND_COMPARISON")
              rm -rf "$compare_dir"
            fi
          else
            rm -rf "$compare_dir"
          fi
        fi
      fi

      # Print summary
      print_summary

      # Exit with error code if ANY tests failed (required or optional)
      if [ ${#failed_required[@]} -gt 0 ] || [ ${#failed_optional[@]} -gt 0 ]; then
        exit 1
      else
        exit 0
      fi
    else # if we are only running one UAM file
      test_uam=$onlyone
      run_a_test
      local status=$?

      if [ $status = 0 ]; then
        exit 0
      else
        # Any test failure (required or optional) fails the workflow
        if [ $status = 2 ]; then
          # Optional backend failure - still fails workflow, but report as warning
          echo ""
          echo "WARNING: Test failed: $test_uam (netcdf - optional feature)"
        fi
        exit 1
      fi
    fi
}



## -------------------------- ##
## ---- Main Entry Point ---- ##
## -------------------------- ##

debug=""
clean=false
distclean=false
config=false
onlyone=false
oversubscribe=""
backends_only=false

do_save=false
do_compare=true

# Where test 14.0's netcdf files wait for the 14.1 RESTART leg (run dir-relative)
append_carry="UA/data/append_carry"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      get_help
      exit 0
      ;;

    -c|--clean)
      echo "Forcing a 'make clean' before compiling!"
      clean=true
      shift
      ;;

    -u|--uninstall)
      echo "Uninstalling with 'make distclean' before running!"
      distclean=true
      config=true
      shift
      ;;

    --nocompare)
      do_compare=false
      echo "Not checking outputs!"
      shift
      ;;

    --save_solution)
      # warnsavesolution is called after args are all parsed, before running.
      do_save=true
      shift      
      ;;

    -o|--only)
      if [ -e "$2" ]; then
        echo "Only running test: $2" 
        onlyone=$2
      else
        echo "ERROR!! Testfile $2 does not appear to exist!"
        exit 1
      fi
      shift 2
      ;;

    --backends_only)
      backends_only=true
      shift
      ;;

    --oversubscribe)
      echo "Oversubscribing! Using 4 threads."
      oversubscribe="--oversubscribe"
      shift
      ;;

    *)
      get_help
      echo "> Unrecognized argument: $1"
      if [ -e $1 ]; then echo "  Run with '-o $1' to test one file"; fi
      exit 1

  # end arg parsing
  esac 
done

# get mad about overwriting ref solutions
if [ $do_save = true ]; then warnsavesolution; fi

# wait a sec to show users that settings are being used
sleep 2

do_tests
