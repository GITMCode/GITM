#!/bin/bash


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
  exit 1
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

  if [ $do_save = true ]; then
    cp UA/data/log*.dat ../ref_solns/log.$test_uam
    echo
    echo
    echo " Created ref solution: log." $test_uam 

    return;
  fi

  if [ $do_compare = true ]; then
    echo
    ../../../share/Scripts/DiffNum.pl -r=5e-5 -a=1e-1 -t ../ref_solns/log.$test_uam UA/data/log*.dat
    
    if [ $? = 0 ]; then
      # test was a success. no differences found.
      return 
    
    else
      echo
      echo " ============    ERROR!!!!!    ============"
      echo "  Output differs from reference solution."
      echo "  Something has gone terribly wrong!!"
      rm GITM.DONE
      return
    fi

  fi

}

run_a_test(){

  printf "\n\n>> Testing with $test_uam ...\n"
  # Copy UAM 
  ln -sf $test_uam UAM.in
  rm -f GITM.DONE

  # Run GITM, stop if error.
  mpirun -np 4 $oversubscribe ./GITM.exe

  # this will either save, or diff, the output log files.
  # (double sanity-check): only run if GITM.DONE exists & above exited with no error code.
  if [[ -f GITM.DONE && $? = 0 ]]; then
    checkoutputs
  fi

  if [ -f GITM.DONE ]; then
      printf "\n\n>>> $test_uam ran successfully! <<< \n\n"
      mv $test_uam $test_uam.success 
      mv UA/data/log*.dat UA/data/log_$test_uam.success
      rm -f GITM.DONE
  else
      printf "\n\n>>> $test_uam   UNSUCCESSFUL! <<< \n\n EXITING\n\n"
      mv UA/data/log*.dat log_$test_uam.fail
      exit 1
  fi

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
        ./Config.pl -install -earth -compiler=gfortran10 $debug
    fi

    make -j

    if [ $? != 0 ]; then
      echo
      echo "Could not compile!"
      echo "The tests have failed. Exiting."
      exit 1
    fi

    # make GITM/run to srcTest/auto_test/run
    # - Old rundir's can cause unintended issues.
    # - Rename GITM/run it if it exists & make a new rundir
    if [[ -d run ]]; then
      echo
      echo " ==========================================================="
      echo " Found a run directory at GITM's root folder!!"
      echo " Moving its contents to run_madebytestscript/, just in case"
      echo " You have 10 seconds to cancel ..."
      echo 
      echo " > ls" $(pwd)"/run"
      ls run/
      sleep 10

      mv run run_madebytestscript
    fi

    # Cnsecutive test runs won't have to sleep for 10 seconds
    # Since we can't tell when (or how) srcTest/auto_test/run was made, replace it
    make rundir
    rm -rf srcTests/auto_test/run
    mv run srcTests/auto_test/

    # Copy the test files into run/
    cd srcTests/auto_test/
    rm -f run/UAM*
    cp UAM.*.test run/

    # begin running:
    cd run/
    if [ $onlyone = false ]; then
      for test_uam in UAM.*.test; do
          run_a_test
      done
      exit 0
    else # if we are only running one UAM file
      test_uam=$onlyone
      run_a_test
      exit 0
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

    -d|--debug)
      echo "Using -debug"
      debug="-debug"
      shift
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
      echo "Unrecognized argument: $1"
      if [ -e $1 ]; then echo "Run with '-o $1' to test one file"; fi
      exit 1

  # end arg parsing
  esac 
done

# get mad about overwriting ref solutions
if [ $do_save = true ]; then warnsavesolution; fi

# wait a sec to show users that settings are being used
sleep 1.5

do_tests
done
