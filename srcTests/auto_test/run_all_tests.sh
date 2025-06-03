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
        -o, --only                Only run this UAM test file.
                                    Useful if a higher number test has failed.
                                    By default, all tests are run alphabetically.
        --skip_config             skip running Config.pl?
        --compare_with [path]     Path to the run directory which has outputs
                                    from all tests (not implemented yet)
        --save_to [path]          Save outputs? Useful to compare when making
                                    changes that could affect outputs. (not implemented yet)

"
  exit 1

}

run_a_test(){

  printf "\n\n>> Testing with $test_uam ...\n"
  # Copy UAM 
  ln -sf $test_uam UAM.in
  rm -f GITM.DONE

  # Run GITM, stop if error.
  mpirun -np 4 --oversubscribe ./GITM.exe
  if [ -f GITM.DONE ]; then
      printf "\n\n>>> $test_uam ran successfully! <<< \n\n"
      mv $test_uam $test_uam.success && rm -f GITM.DONE
  else
      printf "\n\n>>> $test_uam   UNSUCCESSFUL! <<< \n\n EXITING\n\n"
      exit 1
  fi

}

do_tests(){
    # setup run directory
    
    cd ../../ 
    
    if [ $config = true ]; then
        ./Config.pl -install -earth -compiler=gfortran10 $debug
    fi
    
    if [ $clean = true ]; then
      make clean
    fi

    make
    
    if [ ! -f run/GITM.exe ]; then
        # only make rundir if it does not already exist
        make rundir
    fi
    cp -fr run srcTests/auto_test/

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
config=true
onlyone=false

do_save=false
do_compare=false

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
    --skip_config)
      echo "skipping config!"
      config=false
      shift
      ;;
    --compare_with)
      if [[ -d "$2" ]]; then
        echo "Comparing with" $2
        compare_dir=$2
      else
        echo "ERROR: --compare_with directory $2 not found!"
        exit 1
      fi
      shift 2
      ;;
    --save_to)
      if [[ -d "$2" ]]; then
        echo "--save_to directory $2 already exists! Waiting 5 seconds then overwriting."
        echo "   cancel with 'Ctrl C'"
        sleep 5
        compare_dir=$2
      else
        echo "ERROR: --compare_with directory $2 not found!"
        exit 1
      fi
      shift 2
      ;;
    -o|--only)
      if [[ -e "$2" ]]; then
        echo "Only running test: $2" 
        onlyone=$2
      else
        echo "ERROR!! Testfile $2 does not appear to exist!"
        exit 1
      fi
      shift 2
      ;;
    *)
      echo "Unrecognized argument: $1"
      if [[ -e $1 ]]; then echo "Run with '-o $1' to test one file"; fi
      exit 1


  esac

done


do_tests
done