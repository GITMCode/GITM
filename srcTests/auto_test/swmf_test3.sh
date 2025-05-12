#!/bin/bash


get_help(){

  printf "
------------------------------------------------------------------------------------

> This script will automatically run SWMF test3 (UA coupling)
> Just input a directory where to put the files, and the rest will be automatic.
> You can either move this script & it will 'git clone' GITM, or if it is run from
  GITM/srcTests/auto_test, it will copy the local version of GITM to test.

------------------------------------------------------------------------------------
> Example workflow/usage examples (All using the -debug flag):
--------------------------------------------------------------

- You just made changes to your local copy of GITM and want to check they work
  before pushing. This is the most common use case.
  # Go to this folder, run & put things into 'test3':
  > pwd
    /home/[you]/Documents/GITM
  > cd srcTests/auto_test/
  > ./swmf_test3.sh -p run_test3 -O 0

- You made changes on a GITM branch 'couplefix', which have been pushed to
  GITMCode/GITM, and *have not run this before*.
  # Assuming standalone GITM and SWMF will go into ~/Documents:
  > pwd
     /home/[you]/Documents/GITM
  > ls ../
     [Make sure there is NO SWMF here!]

  # Clone the SWMF, switch to the correct GITM branch, just check if it compiles (-O 0):
  >  ./srcTests/auto_test/swmf_test3.sh -d -b couplefix -O 0 -p ../

- The step above did not compile. You made changes to ~/Documents/GITM and want to test 
  before comitting. We do not need to pull or checkout any branches. 
  > pwd
     /home/[you]/Documents/GITM
  
  # Copy your changes to the SWMF:
  > cp -r * ../SWMF/UA/GITM/

  # Run compilation test, do not pull SWMF or checkout a different GITM branch:
  >  ./srcTests/auto_test/swmf_test3.sh -d --skip_config -c -O 0 ../SWMF

- The SWMF did compile! Now run test3 with -O3 (so that it gets done faster):
  > pwd
     /home/[you]/Documents/GITM

  > ./srcTests/auto_test/swmf_test3.sh -d -c ../SWMF

  # If the 'diff' files are empty, good job! If not, something in the outputs has changed.
    If the diff files are not printed, the test did not complete. See the run log with:

  > tail -n 50 ../SWMF/run_test/runlog
  
  (and if the test failed in the restart phase)

  > tail -n 50 ../SWMF/run_test/runlog_restart

------------------------------------------------------------------------------------
Arguments:
----------

        -h, --help                See this information
        -d, --debug               Run installation with -debug enabled?
        -c, --clean               Run a 'make clean' before make-ing?
        --skip_config             Skip running Config.pl?
        --pull                    Run a 'gitall pull', updating the SWMF and GITM repos?
        (-b, --branch)  NAME      Branch name to pull from (on GITMCode/GITM.git).
                                  - If left empty, 'develop' is used. Branch is only
                                    used if --pull is used as well, otherwise GITM
                                    is left in place!
        (-O/--omax)  ##           Optimization level & just_compile. Use a space, please.
                                  - Example:  '-O 0' for fastest compilation.
                                    If this is used, we will not run test3. Just compile.
                                  - By default, -O3 is used, in most cases. Set this to '-O 3' to
                                    'just_compile' and exit, not running test3.
        
        [-p/--path]     PATH      Path to clone SWMF to and run from. Required!
                                  - This can be an existing path.
                                  - -p/--path flags optional, but if no arg is used, 
                                    PATH is assumed to be the last argument.
                                  - If the path provided is not set with -p/--path, and the
                                    path does not exist,the script will print help & exit.


This is tested, but may not be free from bugs. Commit & push, or backup, before running.
Reach out with questions/comments/problems.


tldr: you probably just want to make sure SWMF compiles. do this:
  > cd srcTests/auto_test/
  > ./swmf_test3.sh -p run_test3 -O 0

"

}

do_tests(){
    # setup run directory

    # First, check to see if we're running from GITM's test folder
    if [ -d ../../../GITM ]; then
      # If so, get GITM's absolute path for copying later
      gitm_path=$(pwd)/../../../GITM/
    fi

    # get into the SWMF folder
    cd $testdir
    if [[ -d SWMF/ ]]; then 
      cd SWMF;
    elif [[ -d ../SWMF/ ]]; then
      cd ../SWMF  #Do nothing 
    else
      git clone --depth 1 git@github.com:SWMFSoftware/SWMF.git
      cd SWMF
    fi

    # Get a specific version of GITM, if desired
    if [[ $gitm_path != "false" ]]; then # local gitm
      rsync -a --exclude='srcTests' $gitm_path UA/GITM ;
      # "uninstall"
      cd UA/GITM
      make distclean
      cd ../../

    else # pull a branch
      if [ dopull == true ]; then
        echo "Checking out GITM branch $branch"
        sleep 2
        cd UA/GITM
        git checkout --force $branch
        cd ../../

        echo "Pulling all SWMF dependencies:"
        sleep 2
        gitall pull
      fi
    fi

    # Go thru options:
    if [ $config = true ]; then
        if [ -f Makefile.def ]; then
            # Only clean if running config, don't uninstall.
            make clean
        fi
        ./Config.pl -install $debug $OFlags -v=UA/GITM
    fi
    
    if [ $clean = true ]; then
      make clean
    fi

    # run test3
    if [ $just_compile == true ]; then
      make -j test3_compile
    else
      make -j test3
    fi

    if [ $? -eq 0 ]; then
        printf "\n\n>>> Success! <<< \n\n"
    else
        printf "\n\n>>> FAIL!! See above for more information. <<< \n\n EXITING\n\n"
        printf "If the only error is differing outputs, you may be in the clear."
        exit 1
    fi
    
    exit 0
}



## --------------------------- ##

debug=""
branch="develop"
dopull=false
clean=false
config=true

OFlags=""
just_compile=false

testdir=""

# not set by user; just used by this script
path_set_explicit=false
gitm_path=false

while [[ $# -gt 0 ]]; do
  echo $1
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
      echo "Skipping config!"
      config=false
      shift
      ;;
    -O|--omax)
      echo "Using optimization level -O$2"
      OFlags="-O$2"
      just_compile=true
      shift
      ;;
    -b|--branch)
      echo "Using GITM branch named $2"
      branch=$2
      shift
      ;;
    -p|--path)
      echo "Using the provided path: $2"
      testdir=$2
      path_set_explicit=true
      shift
      ;;
  esac
  shift
done

# Check if the directory is explicitly set
if [[ "$path_set_explicit" != true ]]; then
  testdir=$1
fi

## Check if directory exists
if [[ -d $testdir ]]; then
  echo "Found directory: " $testdir
else
    if [ $path_set_explicit == true ]; then
       echo "Making directory: " testdir
       mkdir -p $testdir
    else
      get_help
      echo
      echo "=>> Provided directory: " $testdir " does not exist!"
      echo "    Refusing to make unless specified by --path."
      echo "Exiting!"
      exit 1
    fi
fi


echo "Waiting 3 seconds then starting."
printf "
Your inputs:
debug= $debug
branch= $branch
dopull= $dopull
clean= $clean
config= $config
OFlags= $OFlags
just_compile= $just_compile
testdir= $testdir

"
sleep 3

do_tests
done
