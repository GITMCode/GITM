#!/bin/bash


get_help(){
    echo \> This script will automatically compile GITM
    echo \& run all tests scripts located within srcTests/auto_test/
    echo
    echo \> New UAM.in files should be placed within this folder,
    echo  with the naming convention UAM.in.__.test
    echo  Information regarding the test configuration may be added, 
    echo  though no spaces can be used.
    echo
    echo Additional notes about the test may be added. Name must conform to: 
    echo "            UAM.\*.test"
    echo and no spaces can be added
    echo
    echo Usage:
    echo
    echo "      [none]         run automatically"
    echo "      -h, --help     see this information"
    echo "      -d, --debug    Configure & compile GITM in -debug"
    echo "      -c, --clean    run a make clean before make-ing?"
    echo "      --skip_config  skip running Config.pl?"
    echo

    exit 1

}


do_tests(){
    # setup run directory
    
    cd ../../ 
    
    if [ $clean -lt 1 ]; then
        ./Config.pl -install -earth -compiler=gfortran10 $debug
    fi
    
    if [ $clean -gt 0 ]; then
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
    cp UAM* run/

    # begin running:
    cd run/
    for test_uam in UAM.*.test; do
        printf "\n\n>> Testing with $test_uam ...\n"
        # Copy UAM (not the first one though)
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
    done
    exit 0
}



## --------------------------- ##

debug=""
clean=0
config=0

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
      echo "cleaning before compiling!"
      clean=1
      shift
      ;;
    --skip_config)
      echo "skipping config!"
      config=1
      shift
      ;;
  esac
#   shift
done


do_tests
done