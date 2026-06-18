#!/bin/bash

# This file writes some 'version' info about GITM to a file
# The file is located in src/.version and is compiled into GITM.
# This info is then written out to run_information.txt.

## Function we call if .git folder exists, and $(git status) has changed from last time
writegitversion(){
    # This is the float-type version placed in output files
    # Format is (date of last commit).(# of files changed from HEAD)
    # first line defines the format, second line puts in the  date of last commit
    # Third line calls git status, counts the number of different files, and strips that of whitespace
    printf 'character(25), parameter :: GitmVersion = & \n "%s.%s"\n' \
        "$(git log -1 --date=format:'%Y%m%d' --pretty='format:%ad')" \
        "$(git status --porcelain | grep -v '^??' | wc -l | awk '{$1=$1};1')" > src/.version


    # The "GitmVersionFull" variable has info about the last commit hash in addition to above
    printf '\ncharacter(50), parameter :: GitmVersionFull = "%s_%s.%s"\n\n' \
        "$(git rev-parse --abbrev-ref HEAD)" \
        "$(git log -1 --date=format:'%Y%m%d' --pretty='format:%h-%ad')" \
        "$(git status --porcelain | grep -v '^??' | wc -l | awk '{$1=$1};1')" >> src/.version

    # Same as above, but for Electrodynamics
    printf '\ncharacter(50), parameter :: ElectrodynamicsVersionFull = "%s_%s.%s"\n\n' \
        "$(git -C ext/Electrodynamics rev-parse --abbrev-ref HEAD)" \
        "$(git -C ext/Electrodynamics log -1 --date=format:'%Y%m%d' --pretty='format:%h-%ad')" \
        "$(git -C ext/Electrodynamics status --porcelain | grep -v '^??' | wc -l | awk '{$1=$1};1')" \
        >> src/.version


    # Here, for completeness, all of the changes files are listed.
    printf 'character(*), parameter :: DifferentFilesGitm = ""' >> src/.version

    files=$(git status --porcelain | awk '{print $NF}')

    # Storing multi-line strings in Fortran is a chore
    # format is var = "some text" // NEW_LINE("A") // & ...
    # NEW_LINE("A") is the newline character, where A is any text (used for decoding)
    # This implementation seemed line the least amount of extra code & conditionals
    # We parse thru the different files (from git status) and for each write the
    # continuation on the prev line, this file name, and a newline character.
    # This allows the first line to be a linrbreak & the last line to not have anything
    for file in $files; do
        printf ' // &\n\"'$file'\" // NEW_LINE(\"A\")' >> src/.version
    done
    echo '' >> src/.version
    echo '' >> src/.version

    printf 'character(*), parameter :: DifferentFilesElectrodynamics = ""' >> src/.version

    if [ -d ext/Electrodynamics ]; then
        files=$(git -C ext/Electrodynamics status --porcelain | awk '{print $NF}')
        for file in $files; do
            printf ' // &\n\"'$file'\" // NEW_LINE(\"A\")' >> src/.version
            # echo  "$file,&" >> src/.version
        done
    fi
    echo '' >> src/.version


}

################
##    MAIN    ##
################

if [ -d .git ]; then # Check if .git/ directory exists

    # check to see if git status has changed, saves time compiling if not
    nFilesDiff=$(git status --porcelain | wc -l)
    if [ -f gitstatus.txt ]; then
        nFilesPrev=$(wc -l < gitstatus.txt) 
    else # Our first time here...
        nFilesPrev=0
    fi

    if [ $nFilesDiff ==  $nFilesPrev ] && [ -f src/.version ]; then
        echo "Version has not changed since last compile: $(head -n 2 src/.version | tail -n 1)"
        exit 0
    else
        writegitversion
        git status --porcelain > gitstatus.txt
        echo "Writing GITM version to file: $(head -n 2 src/.version | tail -n 1)"

    fi

else # if .git/ does not exist, this is probably a release
    cp version.def src/.version
    echo "Wrote version file for GITM: $(head -n 2 src/.version | tail -n 1)"
fi

