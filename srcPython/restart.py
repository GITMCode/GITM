#!/usr/bin/env python

import os
import argparse

IsVerbose = False

#-----------------------------------------------------------------------------
# get arguments from the user
#-----------------------------------------------------------------------------

def get_args():
    parser = argparse.ArgumentParser(
        description = 'mv restartOUT and set up restartIN')
    parser.add_argument('name', \
                        help = 'name to add to restart directory to save -> restart.name')

    parser.add_argument('-q', '--quiet',
                        help = 'Run withOUT verbose',
                        action = 'store_true')
    
    args = parser.parse_args()
    return args

# ----------------------------------------------------------------------
# do system command
# ----------------------------------------------------------------------

def run_command(command):
    if (IsVerbose):
        print("   -> Running Command : ")
        print("      ", command)
    os.system(command)
    return True

# ----------------------------------------------------------------------
# add cd into the UA directory
# ----------------------------------------------------------------------

def add_cd(command):
    outCommand = 'cd UA ; ' + command + ' ; cd - '
    return outCommand

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Main Code!
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

# Get the input arguments
args = get_args()

name = args.name
IsVerbose = not args.quiet

rOut = 'restartOUT'
rIn = 'restartIN'

# Need to see if we are in the UA directory or in the run directory
inUA = False
isFound = False

if (os.path.exists('UA/' + rOut)):
    inUA = False
    isFound = True
else:
    if (os.path.exists(rOut)):
        inUA = True
        isFound = True

if (not isFound):
    print('*** Error ***')
    print('Cant find the directory : ' + rOut)
    print('Make sure you are running this from the run or UA directory')
    print('Can not continue!')
    exit()

dir = 'restart.' + name
if (not inUA):
    dir = 'UA/' + dir

if (os.path.exists(dir)):
    print('*** Error ***')
    print('output path already exists!')
    print('Can not continue!')
    exit()
    
command = 'mv ' + rOut + ' restart.' + name
if (not inUA):
    command = add_cd(command)
run_command(command)

command = 'mkdir ' + rOut
if (not inUA):
    command = add_cd(command)
run_command(command)

command = 'rm -f ' + rIn 
if (not inUA):
    command = add_cd(command)
run_command(command)

command = 'ln -s  restart.' + name + ' ' + rIn
if (not inUA):
    command = add_cd(command)
run_command(command)


