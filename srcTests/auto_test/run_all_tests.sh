#!/bin/bash


get_help(){

  printf "
------------------------------------------------------------------------------------

> This script will automatically compile GITM
  run all tests scripts located within srcTests/auto_test/

> New UAM.in files should be placed within this folder,
  with the naming convention UAM.in.__.test
  Information regarding the test configuration may be added, 
  though no spaces can be used.

Additional notes about the test may be added. Name must conform to: 
            UAM.*.##.test
and no spaces can be added. Numbers do not need to be increasing, but will be useful 
when comparing outputs.

------------------------------------------------------------------------------------
Usage:
------

- To just check that the code compiles and all tests run, use:
    > ./run_all_tests.sh

- If a test fails and you have made edits, you don't need to re-config:
    > ./run_all_tests.sh --skip_config

- For a sanity test, you can force re-configuring & re-compiling everything in 
  debug mode with:
    > ./run_all_tests.sh -c -d

- When making large changes that could affect outputs, it may be useful to
  compare with the latest release. 
  First, commit your changes then 'git checkout main'. Save the outputs from 
  the latest release with:
    > ./run_all_tests.sh --save_to testoutputs_default
  Then, checkout your branch 'git checkout branch_name' and compare the results with:
    > ./run_all_tests.sh --compare_with testoutputs_default

------------------------------------------------------------------------------------
Arguments:
----------

        [none]                    run automatically with defaults.
        -h, --help                see this information
        -d, --debug               Configure & compile GITM in -debug
        -c, --clean               run a 'make clean' before make-ing?
        --skip_config             skip running Config.pl?
        --check_outputs           Check against reference solutions?
        --save_outputs            Save outputs? Useful to run BEFORE making
                                    changes that could affect outputs.
        --only_test [file]        To only run a single test, specify its path here.
                                    Does not support multiple test files.

"
  return
}


do_tests(){

    # go back to root of repo
    cd ../../ 
    
    if [ $config = true ]; then
        ./Config.pl -install -earth -compiler=gfortran10 $debug
    fi
    
    if [ $clean = true ]; then
      make clean
    fi

    # always compile in case any changes were made.
    make
    
    if [ ! -f run/GITM.exe ]; then
      # only make rundir if it does not already exist
      make rundir
      cp -fr run srcTests/auto_test/
    fi

    # move back into folder w this script.
    cd srcTests/auto_test/

    # check if we're only running one test
    if [ ! $only_test_one -eq false ]; then
      rm run/UAM.*.test
      cp $only_test_one run/
    else
      # Copy all test files into run/
      cp UAM* run/
    fi

    # begin running:
    cd run/
    for test_uam in UAM.*.test; do
        printf "\n\n>> Testing with $test_uam ...\n"
        # Copy UAM 
        cp $test_uam UAM.in

        # Run GITM, stop if error.
        mpirun -np 4 ./GITM.exe
        if [ $? -eq 0 ]; then
            printf "\n\n>>> $test_uam ran successfully! <<< \n\n"
            mv $test_uam $test_uam.success
        else
            printf "\n\n>>> $test_uam   UNSUCCESSFUL! <<< \n\n EXITING\n\n"
            exit 1
        fi

        if [ $do_save = true ]; then
          cp data/log00000002.dat ../ref_soln_logs/log.$test_uam
        fi

        if [ $do_compare = true ]; then
          diff_answer=$(diff -y ../ref_soln_logs/log.$test_uam  data/log00000002.dat)
          if [ $? -eq 0 ]; then
            printf "\n\n>>> $test_uam diff'ed successfully! <<< \n\n"
          else
            echo
            echo
            diff ../ref_soln_logs/log.$test_uam  data/log00000002.dat
            echo
            echo
            printf "\n\n>>> $test_uam has differences. UNSUCCESSFUL! <<< \n\n EXITING\n\n"

            exit 1
          fi
        fi


    done
    exit 0
}



## --------------------------- ##

debug=""
clean=false
config=true

do_save=false
do_compare=false
only_test_one=false

while [[ $# -gt 0 ]]; do
  case "$1" in
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

    --skip_config)
      echo "skipping config!"
      config=false
      shift
      ;;

    --check_outputs)
      do_compare=true
      shift
      ;;

    --save_outputs)
      do_save=true
      shift
      ;;

    --only_test)
      echo "Only running test file: " $2
      only_test_one=$2
      shift 2
      ;;

    -h|--help)
      get_help
      ;;

    *)
      get_help
      echo "------------------------------------------------------------------------------------"
      echo
      echo "==> Argument '$1' not recognized."
      echo "See above for help."
      exit 1

  esac

done


do_tests
done