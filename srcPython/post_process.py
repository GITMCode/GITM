#!/usr/bin/env python3

import argparse
import os
from glob import glob
import time
from pGITM import *
from datetime import datetime

IsVerbose = False
DoRm = True

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args_post():

    parser = argparse.ArgumentParser(
        description = "Post process and (optionally) move model results.\n"+
        "- This functions similar to pGITM.py, but can copy files to\n"+
        "  a remote location, or postprocess files as they are created.",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-remotefile',
                        help = 'File that contains info. about remote system',
                        default = 'remote')

    parser.add_argument('-user',
                        help = 'remote user name (default none)',
                        default = 'none')
    
    parser.add_argument('-server',
                        help = 'remote system name (default none)',
                        default = 'none')
    
    parser.add_argument('-dir',
                        help = 'remote directory to use (default none)',
                        default = 'none')
    
    parser.add_argument('-sleep',
                        help = 'how long to sleep between loops in seconds, (default 300)',
                        default = 300, type = int)

    parser.add_argument('-totaltime',
                        help = 'specify how long to run in total in hours, (default 0 - only run once)',
                        default = 0, type = int)
    
    parser.add_argument('-v', '--verbose',
                        help = 'Run with verbose',
                        action = 'store_true')
    
    parser.add_argument('-norm',
                        help = "don't remove any files",
                        action = 'store_true')
    
    parser.add_argument('-tgz',
                        help = "tar and zip raw GITM file instead of process",
                        action = 'store_true')

    parser.add_argument('-nc', action = 'store_true',
                        help = "Postprocess to netCDF files instead on '.bin'?")
    
    parser.add_argument('--combine', action='store_true',
                        help="If processing to netCDF, we can combine each timestep to a"
                        " single file per output type. (ex: 3DALL.nc, etc.). "
                        "Will not work without -nc or if using remote.")
    
    parser.add_argument('--runname', type=str, default='', help=
                        "When combining, this is prepended to the output type in the "
                        "resulting files: '[runname]_3DALL.nc'. Default is no descriptor.")


    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------
# parse remote file
#   - remote file has the form:
# username
# remote server
# remote directory
# ----------------------------------------------------------------------

def parse_remote_file(file):

    if (IsVerbose):
        print('Reading file ', file)
    
    fpin = open(file, 'r')
    user = fpin.readline()
    server = fpin.readline()
    dir = fpin.readline()
    fpin.close()

    remote = {'user': user.strip(),
              'server': server.strip(),
              'dir': dir.strip()}
    return remote

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
# Check inputs:
# ----------------------------------------------------------------------

def check_inputs(user, server, dir):
    
    IsRemote = True
    if ((len(user) == 0) or (user == 'none')):
        if (IsVerbose):
            print("Can't parse user information")
        IsRemote = False
    if ((len(server) == 0) or (server == 'none')):
        if (IsVerbose):
            print("Can't parse server information")
        IsRemote = False
    if ((len(dir) == 0) or (dir == 'none')):
        if (IsVerbose):
            print("Can't parse dir information")
        IsRemote = False

    return IsRemote

# ----------------------------------------------------------------------
# Check the test file to see if it is good
# ----------------------------------------------------------------------

def check_for_stop_file(stopFile = 'stop', delete = False):
    if (os.path.exists(stopFile)):
        DoesExist = True
        print("  --> Stop file exists...")
        if (delete):
            print("  --> Removing Stop file")
            command = "/bin/rm -f " + stopFile
            DidWork = run_command(command)
            # Recheck to see if the stop file still exists...
            DoesExist = os.path.exists(stopFile)
    else:
        DoesExist = False
    return DoesExist

# ----------------------------------------------------------------------
# Check the test file to see if it is good
# ----------------------------------------------------------------------

def parse_test_file(file):

    if (os.path.exists(file)):
        IsGood = True
        if (IsVerbose):
            print('   --> Reading test file ', file)
        fpin = open(file, 'r')
        for line in fpin:
            if (line.find('No such') >= 0):
                IsGood = False
            if (line.find('Connection reset') >= 0):
                IsGood = False
        fpin.close()
    else:
        IsGood = False

    return IsGood
    
# ----------------------------------------------------------------------
# Checks to see if remote file or directory exists
# ----------------------------------------------------------------------

def test_if_remote_exists(user, server, dir):

    DidWork = True
    
    remote_command = user + "@" + server + " 'ls " + dir + "'"
    # check to see if the remote directory exists:

    if (IsVerbose):
        print('Checking to see if remote file or directory exists')
    command = 'ssh ' + remote_command + ' >& .test_file'
    DidWork = run_command(command)
    DidWork = parse_test_file('.test_file')

    if (DidWork):
        if (IsVerbose):
            print('   --> Remote directory (or file) exists!')
    else:
        print('--> Remote directory (or file) does NOT exist!')
        print('    Need to make this directory!')
        
    return DidWork

# ----------------------------------------------------------------------
# Make a remote directory
# ----------------------------------------------------------------------

def make_remote_dir(user, server, dir):

    DidWork = True
    
    remote_command = user + "@" + server + " 'mkdir " + dir + "'"
    # check to see if the remote directory exists:

    print('Making remote directory : ', dir)
    command = 'ssh ' + remote_command + ' >& .test_file'
    DidWork = run_command(command)
    DidWork = parse_test_file('.mkdir_command')

    return DidWork

# ----------------------------------------------------------------------
# Push log and data files to remote server
# ----------------------------------------------------------------------

def transfer_log_files(user, server, dir):

    remote = user + '@' + server + ':' + dir
    files = 'job* log.* UAM* *.txt *.dat*'
    outfile = '.output_rsync_log'
    command = 'rsync -vrae ssh ' + files + ' ' + remote + ' >& ' + outfile
    DidWork = run_command(command)

    return DidWork
    
# ----------------------------------------------------------------------
# Determine base filenames from header file
# ----------------------------------------------------------------------

def determine_base(headerFile, dir):

    start = len(dir)
    end = len(headerFile) - len('.header')
    
    baseFile = headerFile[start:end]

    return baseFile
    
# ----------------------------------------------------------------------
# tar and zip raw GITM files
# ----------------------------------------------------------------------

def tar_and_zip_gitm():

    print('-> Tar and Zipping files!')
    data_here = 'UA/data'
    
    # get header file list:
    headerFiles = sorted(glob(data_here + '/*.header'))

    for headerFile in headerFiles:
        baseFile = determine_base(headerFile, data_here + '/')
        print('--> Processing header file : ', baseFile)

        tarFile = baseFile + '.tgz'
        command = 'cd ' + data_here + ' ; rm -f ' + tarFile + ' ; cd ../..'
        DidWork = run_command(command)
        
        command = 'cd ' + data_here + ' ; tar -cvzf ' + tarFile + ' ' + baseFile + '.* ; cd ../..'
        DidWork = run_command(command)

        # Remove raw files:
        command = 'cd ' + data_here + ' ; rm -f ' + baseFile + '.b[0-9][0-9][0-9][0-9] ; cd ../..'
        DidWork = run_command(command)
        command = 'cd ' + data_here + ' ; rm -f ' + baseFile + '.header ; cd ../..'
        DidWork = run_command(command)
        command = 'cd ' + data_here + ' ; rm -f ' + baseFile + '.sat ; cd ../..'
        DidWork = run_command(command)

        # There are times when you want to break out of this loop,
        # so allow user to break loop if stop file exists:
        stopCheck = check_for_stop_file()
        if (stopCheck):
            print("  --> Stopping!")
            break

    DidWork = True
        
    return DidWork

## ----------------------------------------------------------------------
## post process GITM files
## ----------------------------------------------------------------------
#
#def post_process_gitm():
#
#    command = \
#        'cd UA ; ' + \
#        'chmod a+rx data ; '
#    if (IsVerbose):
#        command = command + './pGITM ; '
#    else:
#        command = command + './pGITM > .post_process ; '
#    command = command + 'chmod a+r data/*.bin ; cd ..'
#    DidWork = run_command(command)
#
#    return DidWork

# ----------------------------------------------------------------------
# Transfer file, check if it made it, then delete local (if requested)
# ----------------------------------------------------------------------

def transfer(filelist, user, server, dir, DoRemove):

    DidWork = True
    remote = user + '@' + server + ':' + dir

    files = ''
    outfile = '.output_rsync_log'

    if (len(filelist) > 0):
        for file in filelist:
            chmod = 'chmod a+r ' + file
            DidWork = run_command(chmod)
            files = files + ' ' + file
            
        rsync = 'rsync -rav ' + files + ' ' + remote
        if (not IsVerbose):
            rsync = rsync + ' >> ' + outfile + ' 2>&1'
        DidWork = run_command(rsync)
        
    if (DoRemove):
        for file in filelist:
            sep = file.split('/')
            test_file = sep[-1]
            DidTransfer = test_if_remote_exists(user,
                                                server,
                                                dir + '/' + test_file)
            if (DidTransfer):
                if (IsVerbose):
                    print('   --> Remote file (' + test_file + ') exists!' +
                          '  Deleting local!')
                if (DoRm):
                    command = '/bin/rm -f ' + file
                    DidWork = run_command(command)
                    # Systems reject ssh commands if too many happen in
                    # too short of time, so sleep in between the commands. 
                    time.sleep(5)
            else:
                if (IsVerbose):
                    print('Remote file (' + test_file + ') does not exist!')

    return DidWork

# ----------------------------------------------------------------------
# Transfer a bunch of different file types to remote machine:
# ----------------------------------------------------------------------

def transfer_model_output_files(user, server, dir):

    data_here = 'UA/data'
    remote = user + '@' + server + ':' + dir

    print("Transfering files...")
    
    # run info file:
    file = 'run_information.txt'
    filelist = sorted(glob(data_here + '/' + file))
    DoRemove = True
    DidWork = transfer(filelist, user, server, dir, DoRemove)

    file = 'log*.dat'
    filelist = sorted(glob(data_here + '/' + file))
    DoRemove = False
    DidWork = transfer(filelist, user, server, dir, DoRemove)

    # data files - remove by default
    DoRemove = True
    file = '3D*.bin'
    filelist = sorted(glob(data_here + '/' + file))
    DidWork = transfer(filelist, user, server, dir, DoRemove)

    file = '2D*.bin'
    filelist = sorted(glob(data_here + '/' + file))
    DidWork = transfer(filelist, user, server, dir, DoRemove)

    file = '1D*.bin'
    filelist = sorted(glob(data_here + '/' + file))
    DidWork = transfer(filelist, user, server, dir, DoRemove)

    file = '0D*.bin'
    filelist = sorted(glob(data_here + '/' + file))
    DidWork = transfer(filelist, user, server, dir, DoRemove)

    # tar and zip - remove by default
    DoRemove = True
    file = '3D*.tgz'
    filelist = sorted(glob(data_here + '/' + file))
    DidWork = transfer(filelist, user, server, dir, DoRemove)

    file = '2D*.tgz'
    filelist = sorted(glob(data_here + '/' + file))
    DidWork = transfer(filelist, user, server, dir, DoRemove)

    file = '1D*.tgz'
    filelist = sorted(glob(data_here + '/' + file))
    DidWork = transfer(filelist, user, server, dir, DoRemove)

    file = '0D*.tgz'
    filelist = sorted(glob(data_here + '/' + file))
    DidWork = transfer(filelist, user, server, dir, DoRemove)

    return DidWork

# ----------------------------------------------------------------------
# Post process and then transfer files once:
# ----------------------------------------------------------------------

def do_loop(doTarZip, user, server, dir, IsRemote,
            write_nc=False, combine=False, runname='', # For writing outputs to netCDF
            ):

    DidWork = True

    # 1 - check to see if the remote directory exists:
    if (IsRemote):
        DidWork = test_if_remote_exists(user, server, dir)
        if (not DidWork):
            # if the directory doesn't exist, try to make it:
            DidWork = make_remote_dir(user, server, dir)
            # then check again to see if it exists:
            DidWork = test_if_remote_exists(user, server, dir)
            if (not DidWork):
                return DidWork

    # 2 - push log files
    if (IsRemote and DidWork):
        DidWork = transfer_log_files(user, server, dir)
        if (not DidWork):
            return DidWork

    # 3 - Post process GITM files:
    print('Post Processing GITM files...')
    if (doTarZip):
        DidWork = tar_and_zip_gitm()
    else:
        # Default should be UA/data
        processDir = 'UA/data' 
        if (not os.path.exists(processDir)):
            # Maybe we are in the UA directory?
            processDir = 'data'
            if (not os.path.exists(processDir)):
                # Maybe we are already in the data directory???
                processDir = '.'
        DidWork = post_process_gitm(processDir, DoRm, isVerbose = IsVerbose,
                                    write_nc=write_nc, combine=combine, runname=runname)

    # 4 - Check if remote data directory exists, make it if it doesn't:
    data_remote = '/data'
    if (IsRemote and DidWork):
        dir = dir + data_remote
        DidWork = test_if_remote_exists(user, server, dir)
        if (not DidWork):
            DidWork = make_remote_dir(user, server, dir)
            DidWork = test_if_remote_exists(user, server, dir)
            if (not DidWork):
                return DidWork

    # 5 - transfer output file:
    if (IsRemote and DidWork):
        DidWork = transfer_model_output_files(user, server, dir)
    
    return DidWork

# ----------------------------------------------------------------------
# load remote file
# ----------------------------------------------------------------------

def load_remote_file(args):

    remoteFile = args.remotefile
    
    # figure out remote system information:
    if (os.path.exists(remoteFile)):
        print('Found file: ', remoteFile)
        remote = parse_remote_file(remoteFile)
        user = remote['user']
        server = remote['server']
        dir = remote['dir']
    else:
        user = args.user
        server = args.server
        dir = args.dir

    # Check remote system inputs:
    IsRemote = check_inputs(user, server, dir)
    return IsRemote, user, server, dir
    
# ----------------------------------------------------------------------
# Main Code
# ----------------------------------------------------------------------

if __name__ == '__main__':  # main code block

    args = parse_args_post()

    # make local variables for arguments:
    doTarZip = args.tgz
    IsVerbose = args.verbose
    if (args.norm):
        DoRm = False

    # Sanity check netCDF-related arguments...
    if args.nc or args.combine or (args.runname != ''):
        if not args.nc:
            raise ValueError(f"Combine and Runname can only be used with NetCDF files!")
        if (args.runname != '') and not args.combine:
            print(f"Using runname '{args.runname}' and combining")
            args.combine = True
        # For compatibility, make sure we can import PyITM if we're using netcdf files
        try:
            import pyitm
        except ImportError:
            raise ImportError("\n>> PyITM is required for NetCDF post-processing!\n"
                              " It is available at https://github.com/GITMCode/PyITM.git")

    # Check if stop file exists:
    check_for_stop_file(delete = True)

    # Very simple main loop - run forever
    DidWork = True
    startTime = datetime.now()
    while DidWork:
        # load remote file every iteration so we can change it if needed:
        IsRemote, user, server, dir = load_remote_file(args)
        if (IsVerbose):
            print('Move files to remote system? ', IsRemote)
            print('Process into netCDF files? ', args.nc)
        if IsRemote and args.nc and args.combine:
            raise ValueError(
                "The remote transfer & netcdf combine options cannot be used together")

        DidWork = do_loop(doTarZip, user, server, dir, IsRemote,
                          args.nc, args.combine, args.runname)
        if (DidWork):
            # Check if stop file exists:
            stopCheck = check_for_stop_file()
            if (stopCheck):
                print("  --> Stopping due to stop file existing!")
                # want to break out of loop, so set loop breaking condition:
                DidWork = False
            else:
                currentTime = datetime.now()
                dt = ((currentTime - startTime).total_seconds())/3600.0
                if (dt > args.totaltime):
                    if args.totaltime == 0:
                        # Different exit message for non-continuous runs 
                        print(" -> All done!")
                    else:
                        print(" -> Stopping due to totaltime exceeded!")
                    # want to break out of loop, so set loop breaking condition:
                    DidWork = False

        if (DidWork):
            # everything ok, go to sleep for a bit
            print('Sleeping ... ', args.sleep, ' sec.')
            time.sleep(args.sleep)
