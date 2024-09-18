#!/bin/bash

# This script will automatically compile GITM
# & run all tests scripts located within srcTests/auto_test/

# New UAM.in files should be placed within this folder,
# with the naming convention UAM.in.##.test

# Additional notes about the test may be added. Name must conform to: 
# > UAM.*.test
# but the beginning characters "UAM.in."* must not be changed,
# and no spaces can be added

## --------------------------- ##

# setup run directory
cd ../../ 
./Config.pl -install -earth -compiler=gfortran10
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



