#!/bin/bash


get_help(){

  printf "
------------------------------------------------------------------------------------

> This script will automatically run SWMF test3 (UA coupling)
> Just input a directory where to put the files, and the rest will be automatic.


------------------------------------------------------------------------------------
> Example workflows (All using the -debug flag):
-------------------

- You just made changes on a GITM branch "couplefix", which have been pushed to
  GITMCode/GITM, and *have not run this before*.
  - Assuming standalone GITM and SWMF will go into ~/Documents:
  > pwd
     /home/[you]/Documents/GITM
  > ls ../
     [Make sure there is NO SWMF here!]

  - Clone the SWMF, switch to the correct GITM branch, just compile:
  >  ./srcTests/auto_test/swmf_test3.sh -d -b couplefix -O 0 -p ../

- The step above did not compile. You made changes to ~/Documents/GITM and want to test 
  before comitting. We do not need to pull or checkout any branches. 
  > pwd
     /home/[you]/Documents/GITM
  
  - Copy your changes to the SWMF:
   > cp -r * ../SWMF/UA/GITM/

  - Run compilation test, do not pull SWMF or checkout a different GITM branch:
  >  ./srcTests/auto_test/swmf_test3.sh -d --skip_config -c -O 0 ../SWMF

- The SWMF did compile! Now run test3 with -O3 (so that it gets done faster):
  > pwd
     /home/[you]/Documents/GITM

  > ./srcTests/auto_test/swmf_test3.sh -d -c ../SWMF

  - If the 'diff' files are empty, good job! If not, something in GITM has changed.
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
                                  - If left empty, "develop" is used. Branch is only
                                    used it --pull is used as well, otherwise GITM
                                    is left in place!
        (-O/--omax)     #         Optimization level. Use a space, please.
                                  - Example:  '-O 0' for fastest compilation.
                                    If this is used, we will not run test3. Just compile.
        
        [-p/--path]     PATH      Path to clone SWMF to and run from. Required!
                                  - This can be an existing path.
                                  - -p/--path flags optional, but if no arg is used, 
                                    PATH is assumed to be the last argument.


(This is tested, but may not be free from bugs. Commit & push before running.)
(Reach out with questions/comments/problems.)

"

exit 1
}

do_tests(){
    # setup run directory
    
    cd $testdir

    if [ -d SWMF/ ]; then cd SWMF; fi

    # Check if we're inside swmf folder with git status:
    git status --porcelain
    if [ $? -eq 0 ]; then
      # We're probably within SWMF... Don't need to clone.
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

    else
      # Probably not inside an existing SWMF folder...
      git clone --depth 1 git@github.com:SWMFSoftware/SWMF.git
      cd SWMF
    fi


    
    if [ $config = true ]; then
        if [ -f Makefile.def ]; then
            # Just make clean if running config, don't uninstall.
            make clean
        fi

        ./Config.pl -install $debug $OFlags -v=UA/GITM
    fi
    


    if [ $clean = true ]; then
      make clean
    fi

    if [ just_compile = true ]; then
      make -j 
    else
      make -j test3
    fi

    if [ $? -eq 0 ]; then
        printf "\n\n>>> Success! <<< \n\n"
        mv $test_uam $test_uam.success
    else
        printf "\n\n>>> FAIL!! See above for more information. <<< \n\n EXITING\n\n"
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


while [[ $# -gt 1 ]]; do
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
      echo "Using max. optimization level -O$2"
      OFlags="-O$2"
      just_compile=true
      shift 2
      ;;
    -b|--branch)
      echo "Using GITM branch names $2"
      branch=$2
      shift 2
      ;;
    -p|--path)
      echo "Using GITM branch names $2"
      branch=$2
      shift 2
      ;;
  esac

done

## Get directory. remaining argument if not set already
if [ $# -gt 0 ]; then
  if [[ -d $1 ]]; then
    echo "Found directory $1"
    testdir=$1
  else
    echo "Directory $1 not found! Exiting."
    echo "Found remaining $# arguments"
    echo "There could be a probelm with your arguments given!"
    echo "Try again please. Bye."
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
OFlags= $OFLAGS
just_compile= $just_compile
testdir= $testdir

"
sleep 3

do_tests
done