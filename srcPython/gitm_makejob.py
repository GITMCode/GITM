#!/usr/bin/env python3

import argparse
import re
from datetime import datetime
from datetime import timedelta
import os

IsVerbose = False

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args():

    parser = argparse.ArgumentParser(description = 'Prepare GITM job file')

    parser.add_argument('-input', \
                        help = 'Use this file as the baseline jobfile', \
                        default = 'job_file_ivy')
    
    parser.add_argument('-gitm', \
                        help = 'Base directory to find GITM', \
                        default = 'none')
    
    parser.add_argument('-name', \
                        help = 'job name', \
                        default = 'boring_name')
    
    parser.add_argument('-wall', \
                        help = 'walltime (just the hours portion)', \
                        default = '24')
    
    parser.add_argument('-machine', \
                        help = 'machine to run on (ivy)', \
                        default = 'ivy')
    
    parser.add_argument('-gid', \
                        help = 'group id for charging run', \
                        default = 'none')
    
    parser.add_argument('-procspernode', \
                        help = 'processors per node', \
                        default = '20')
    
    parser.add_argument('-exe', \
                        help = 'executable to run (GITM.exe)', \
                        default = 'GITM.exe')
    
    parser.add_argument('-log', \
                        help = 'logfile to output to (log.txt)', \
                        default = 'log.txt')
    
    parser.add_argument('-cpus', \
                        help = 'number of cpus to use (the total for the job)', \
                        default = '100')
    
    parser.add_argument('-output', \
                        help = 'Output filename', \
                        default = 'job')
    
    parser.add_argument('-restart', \
                        help='uncomment the cp UAM.in.Restart line', \
                        action="store_true")
    
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

#-----------------------------------------------------------------------------
# Reads a job file into a dictionary
#-----------------------------------------------------------------------------

def read_job(file):

    fpin = open(file, 'r')
    lines = fpin.readlines()
    fpin.close()
    
    job = []

    for line in lines:
        line = line.strip('\n')
        job.append(line)

    return job

# ----------------------------------------------------------------------
#
# ----------------------------------------------------------------------

def write_job(file, job):
    print('--> Writing file: ', file)
    fp = open(file, 'w')
    for line in job:
        fp.write(line + "\n")
    fp.close()
    command = 'chmod a+x ' + file
    didWork = False
    didWork = run_command(command)
    return didWork

# ----------------------------------------------------------------------
# modify job script with inputs
# ----------------------------------------------------------------------

def modify_job(inLines, \
               name = 'test_job_name', \
               walltime = '02',
               gid = 'unknown',
               cpus = '100',
               nProcsPerNode = 20,
               machine = 'ivy',
               exe = 'GITM.exe',
               log = 'log.txt',
               restart = False):

    outLines = []

    nodes = int(float(cpus) / nProcsPerNode)
    if (nodes * nProcsPerNode < int(cpus)):
        nodes += 1
    
    for line in inLines:
        m = re.match(r'#PBS -N (.*)', line)
        if m:
            oldname = m.group(1)
            line = '#PBS -N ' + name
        m = re.match(r'.*walltime.*', line)
        if m:
            line = '#PBS -lwalltime=' + walltime + ':00:00'
        m = re.match(r'.*group_list.*', line)
        if m:
            line = '#PBS -W group_list=' + gid
        m = re.match(r'.*ncpus.*', line)
        if m:
            line = '#PBS -lselect=%d:ncpus=%d' % \
                (int(nodes), int(nProcsPerNode))
            line = line + ':model=' + machine
        m = re.match(r'.*mpiexec.*', line)
        if m:
            line = 'mpiexec -np ' + cpus + ' ' + exe + ' > ' + log
        if (restart):
            m = re.match(r'#(.*UAM.in.Restart.*)', line)
            if m:
                line = m.group(1)
            m = re.match(r'.*UAM.in.Start.*', line)
            if m:
                line = '#' + line
        outLines.append(line)

    return outLines

# ----------------------------------------------------------------------
# main code:
# ----------------------------------------------------------------------

args = parse_args()

gitmDir = args.gitm
if (os.path.exists(gitmDir)):
    useGitmDir = True
    gitmDir = gitmDir + '/'
else:
    useGitmDir = False

# ----------------------------------------
# read baseline UAM.in file:

file = args.input

inputFound = False
if (os.path.exists(file)):
    inputFound = True
else:
    if (useGitmDir):
        file = gitmDir + 'GITM_Input_Files/Pleiades/'+file
        if (os.path.exists(file)):
            inputFound = True
if (inputFound):
    print('--> Reading file : ', file)
    job = read_job(file)
else:    
    print('Can not find input file ', file, '!')
    print('  Please run with -h if you need help!')
    exit()

newjob = modify_job(job, \
                    name = args.name, \
                    walltime = args.wall, \
                    gid = args.gid,
                    cpus = args.cpus, \
                    nProcsPerNode = int(args.procspernode), \
                    machine = args.machine, \
                    exe = args.exe, \
                    log = args.log, \
                    restart = args.restart)

didWork = write_job(args.output, newjob)



