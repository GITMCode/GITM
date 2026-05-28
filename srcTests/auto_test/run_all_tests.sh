#!/bin/bash

# Exit immediately if an unhandled pipeline/command fails
set -e

get_help(){

  printf "
------------------------------------------------------------------------------------

> This script will automatically compile GITM and
  run the test cases located within srcTests/auto_test/
  - All UAM files matching the pattern UAM.*.test
    are automatically tested on GitHub

> When adding a new feature, it is recommended to create a test.
  - To do this, first create a UAM.in file which uses the new feature.
  - Then, run this script with:
    > ./run_all_tests.sh -d -c --save_solution -o [your uam file]
  - Add a few words to the file name to explain the test. See below for more.
  - Please try to keep the numbers increasing sequentially, if possible.

Additional notes about the test may be added. Name must conform to:
            UAM.in.##.*.test
and no spaces can be added. Numbers do not need to be increasing, but will be useful
when comparing outputs. It is recommended to follow the pattern:
            UAM.in.##.[anything_you_want].test

------------------------------------------------------------------------------------
Usage:
------

> First, from root of repository, install GITM with './Config.pl [...]'
> Next, 'cd srcTests/auto_test/'

- To just check that the code compiles and all tests run, use:
    > ./run_all_tests.sh

- If a test fails and you have made edits, you don't need to re-run every test.
  (substitute for another test):
    > ./run_all_tests.sh --only UAM.in.00.default.test

- For a sanity check, you can force the code to be re-compiled with the (-c/--clean) 
  flag:
    > ./run_all_tests.sh -c

------------------------------------------------------------------------------------
Arguments:
----------

        [none]                    run automatically with defaults.
        -h, --help                see this information
        -c, --clean               Run a 'make clean' before make-ing?
        -u, --uninstall           Force a 'make distclean' before running?
                                    This effectively uninstalls GITM before running. 
                                    Useful when dependencies change, or when a lot of
                                    changes have been made (by you or git). 
                                    * Will only work with gfortran10 & in debug mode. *
        -o, --only                Only run this UAM test file.
                                    Useful if a higher number test has failed.
                                    By default, all tests are run alphabetically.
        --nocompare               Default configuration diff's the log file.
                                    Use this if you want to only do a run & not
                                    check the output.
        --save_solution           Rewrite (or create) reference solution?
                                    - Similar to the SWMF's '-BLESS'.
                                    - Will take a moment to rewrite a solution in case
                                      this is set erroreously.
                                    - If a test was overwritten by mistake, use
                                      'git restore' to put it back.
        --oversubscribe           Run 'mpirun' with the '--oversubscribe' flag?
                                    This is necessary on GitHub's CI/CD servers, but
                                    is not be necessary for *most* users. As some MPI
                                    libraries implement this differently, this is
                                    off by default.

"
}

get_backend_from_filename(){
  # Extract backend from filename: UAM.in.##.description.BACKEND.test
  # Returns: legacy, mpiio, netcdf, or "unknown" if not found
  local filename="$1"

  if [[ "$filename" =~ \.(legacy|mpiio|netcdf)\.test$ ]]; then
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
!!!!!!!!!!!!!!!!!    WARNING    !!!!!!!!!!!!!!!!!

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

  if [ $do_save = true ]; then
    cp UA/data/log*.dat ../ref_solns/log.$test_uam
    echo
    echo
    echo " Created ref solution: log." $test_uam
    return 0
  fi

  if [ $do_compare = true ]; then
    echo
    # Capture diff status without triggering set -e
    local diff_status=0
    ../../../share/Scripts/DiffNum.pl -r=5e-5 -a=1e-1 -t ../ref_solns/log.$test_uam UA/data/log*.dat || diff_status=$?

    if [ $diff_status = 0 ]; then
      # test was a success. no differences found.
      return 0

    else
      echo
      echo " ============    ERROR!!!!!    ============"
      echo "  Output differs from reference solution."
      echo "  Something has gone terribly wrong!!"
      rm GITM.DONE
      return 1
    fi
  fi
  return 0
}

run_a_test(){
  # Run a single test and return status:
  #   0 = success
  #   1 = failure (required backend)
  #   2 = failure (optional backend)

  printf "\n\n>> Testing with $test_uam ...\n"
  # Copy UAM
  ln -sf $test_uam UAM.in
  rm -f GITM.DONE

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
      return 0
  else
      printf "\n\n>>> $test_uam   UNSUCCESSFUL! <<< \n\n"
      mv UA/data/log*.dat log_$test_uam.fail 2>/dev/null

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
  echo "✓ Passed: $passed_count"

  if [ $failed_required_count -gt 0 ]; then
    echo "> Failed (REQUIRED): $failed_required_count"
    for test in "${failed_required[@]}"; do
      echo "   - $test"
    done
  fi

  if [ $failed_optional_count -gt 0 ]; then
    echo "> Failed (optional): $failed_optional_count"
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
      for test_uam in UAM.*.test; do
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
          echo "⚠ Test failed: $test_uam (netcdf - optional feature)"
        fi
        exit 1
      fi
    fi
}



## --------------------------- ##

debug=""
clean=false
distclean=false
config=false
onlyone=false
oversubscribe=""

do_save=false
do_compare=true

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      get_help
      exit 1
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
