#!/usr/bin/env python3
'''
Handle all unconcatenated output fragments 
produced during a GITM simulation.

This script works robustly and securely on a range of
unix-like environments with Python available
(i.e., any modern computer.)
'''

import os
from glob import glob
from subprocess import run
import argparse
import re
import datetime
import numpy as np

from scipy.io import FortranFile
from struct import unpack

# This should work on most systems:
endChar='<'

# ----------------------------------------------------------------------------
# 
# ----------------------------------------------------------------------------

def get_args_pgitm():
    parser = argparse.ArgumentParser(
        description = 'Post-process GITM model results')
    parser.add_argument('-v',  \
                        action='store_true', default = False, \
                        help = 'set verbose to true')
    parser.add_argument('-norm',  \
                        action='store_true', default = False, \
                        help = 'DO NOT REMOVE processed files')
    parser.add_argument('-dir',
                        default = 'UA/data',
                        help = 'Directory to find unprocessed files')
    parser.add_argument('-header',
                        default = 'none',
                        help = 'Can process just one file if we want')
    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------------
# 
# ----------------------------------------------------------------------------

def write_gitm_file(file, data, isVerbose = False):

    if (isVerbose):
        print(' -> Writing Fortran Binary file : ', file)
    fp = FortranFile(file, 'w')

    val = np.array([data["version"]], dtype = np.float64)
    fp.write_record(val)
    iVals = [data["nLonsTotal"], data["nLatsTotal"], data["nAltsTotal"]]
    fp.write_record(np.array(iVals, dtype = np.int32))

    # Write out variables!
    val = np.array([data["nVars"]], dtype = np.int32)
    fp.write_record(val)
    for var in data["vars"]:
        varpad = var.ljust(40).encode('utf-8')
        fp.write_record(varpad)

    # Write out time:
    ymdhmsm = np.array([data["time"].year, \
                        data["time"].month, \
                        data["time"].day, \
                        data["time"].hour, \
                        data["time"].minute, \
                        data["time"].second, \
                        int(data["time"].microsecond/1000.0)], \
                      dtype=np.int32)
    fp.write_record(ymdhmsm)
            
    for iVar in np.arange(data["nVars"]):
        v = data["vars"][iVar]
        vals = np.array(data[v], dtype = np.float64)
        fp.write_record(vals)
    
    fp.close()
    

# ----------------------------------------------------------------------------
# 
# ----------------------------------------------------------------------------

def read_record_float(fp, nPoints):
    recLen = unpack(endChar+'l',fp.read(4))[0]
    values = np.array(unpack(endChar+'%id'%(nPoints),fp.read(recLen)))
    endLen = np.array(unpack(endChar+'l',fp.read(4))[0])
    return values
    
# ----------------------------------------------------------------------------
# We want to eliminate all interior ghostcells, so we need to figure out
# whether we are on interior or exterior block, and then get the range of
# indices to be correct (i.e., include left ghosts or right ghosts or neither)
# ----------------------------------------------------------------------------

def determine_bounds(iBlock, nBlocks, nPts, nGCs):
    
    if (iBlock == 0):
        iWholeStart = 0
        iStart = 0
    else:
        iWholeStart = nGCs + iBlock * nPts
        iStart = nGCs
    if (iBlock == nBlocks - 1):
        iWholeEnd = nGCs + (iBlock + 1) * nPts + nGCs
        iEnd = nGCs + nPts + nGCs
    else:
        iWholeEnd = nGCs + (iBlock + 1) * nPts
        iEnd = nGCs + nPts
    return iWholeStart, iStart, iWholeEnd, iEnd

# ----------------------------------------------------------------------------
# 
# ----------------------------------------------------------------------------

def read_one_binary_file(file, nVars, nLons, nLats, nAlts, isVerbose = False):
    # Store in reversed order, due to the whole fortran / python order diff
    data = np.zeros((nVars, nAlts, nLats, nLons))
    if (isVerbose):
        print('   --> Reading file : ', file)
    fpin = open(file, 'rb')
    for iAlt in range(nAlts):
        for iLat in range(nLats):
            for iLon in range(nLons):
                data[:, iAlt, iLat, iLon] = read_record_float(fpin, nVars)
    fpin.close()
    return data

# ----------------------------------------------------------------------------
# 
# ----------------------------------------------------------------------------

def delete_files(header, isVerbose = False):

    if (isVerbose):
        print(' -> Deleting files associated with header : ', header)

    
# ----------------------------------------------------------------------------
# 
# ----------------------------------------------------------------------------

def read_header(headerFile, isVerbose = False):

    print(' --> Reading file : ', headerFile)
    fpin = open(headerFile, 'r')
    lines = fpin.readlines()
    fpin.close()

    nGhostsLat = 2
    nGhostsLon = 2
    nGhostsAlt = 2
    useGhostCells = True

    nBlocksAlt = 1
    nBlocksLat = 1
    nBlocksLon = 1

    nVars = 0
    nAlts = 0
    
    isDone = False
    iLine = 0
    nLines = len(lines)

    def get_int(line):
        l7 = line[0:7]
        smashed = l7.replace(" ", "")
        return int(smashed)
    
    headerInfo = {}
    while (not isDone):
        line = lines[iLine].strip('\n')
        smashed = line.replace(" ", "")
        line = line.replace("  ", "")

        m = re.match(r'BLOCKS', smashed)
        if m:
            nBlocksAlt = get_int(lines[iLine + 1])
            nBlocksLat = get_int(lines[iLine + 2])
            nBlocksLon = get_int(lines[iLine + 3])
            iLine += 3

        m = re.match(r'NGHOSTCELLS', smashed)
        if m:
            nGhostsAlt = get_int(lines[iLine + 1])
            nGhostsLat = get_int(lines[iLine + 2])
            nGhostsLon = get_int(lines[iLine + 3])
            iLine += 3
            
        m = re.match(r'NOGHOSTCELLS', smashed)
        if m:
            nGhostsLat = 0
            nGhostsLon = 0
            nGhostsAlt = 0
            useGhostCells = False

        m = re.match(r'NUMERICAL', smashed)
        if m:
            nVars = get_int(lines[iLine + 1])
            headerInfo['nVars'] = nVars
            nAlts = get_int(lines[iLine + 2])
            nLats = get_int(lines[iLine + 3])
            nLons = get_int(lines[iLine + 4])
            iLine += 4

        m = re.match(r'VARIABLE', smashed)
        if m:
            headerInfo['vars'] = []
            for iV in range(headerInfo['nVars']):
                v = lines[iLine + 1 + iV][7:]
                v = v.strip('\n')
                v = v.replace(" ", "")
                nv = re.sub('!U', '', v)
                nv = re.sub('!N', '', nv)
                nv = re.sub('!D', '', nv)
                if (isVerbose):
                    print('   --> Var ', iV, ': ', nv)
                headerInfo['vars'].append(nv)
            
        m = re.match(r'TIME', smashed)
        if m:
            iYear = get_int(lines[iLine + 1])
            iMonth = get_int(lines[iLine + 2])
            iDay = get_int(lines[iLine + 3])
            iHour = get_int(lines[iLine + 4])
            iMinute = get_int(lines[iLine + 5])
            iSecond = get_int(lines[iLine + 6])
            iMilli = get_int(lines[iLine + 7])
            headerInfo['time'] = \
                datetime.datetime(iYear, iMonth, iDay, iHour, iMinute, iSecond, iMilli)
            iLine += 4
            
        m = re.match(r'VERSION', smashed)
        if m:
            line = lines[iLine + 1]
            smashed = line.replace(" ", "")
            headerInfo['version'] = float(smashed)

        m = re.match(r'END', smashed)
        if m:
            isDone = True
        iLine += 1
        if (iLine >= nLines):
            isDone = True

    headerInfo['useGhostCells'] = useGhostCells
    headerInfo['nLons'] = nLons - 2 * nGhostsLon
    headerInfo['nBlocksLon'] = nBlocksLon
    headerInfo['nGhostsLon'] = nGhostsLon

    headerInfo['nLats'] = nLats - 2 * nGhostsLat
    headerInfo['nBlocksLat'] = nBlocksLat
    headerInfo['nGhostsLat'] = nGhostsLat
    
    headerInfo['nAlts'] = nAlts - 2 * nGhostsAlt
    headerInfo['nBlocksAlt'] = nBlocksAlt
    headerInfo['nGhostsAlt'] = nGhostsAlt
    
    return headerInfo

# ----------------------------------------------------------------------------
# 
# ----------------------------------------------------------------------------

def remove_files(header, isVerbose = False):

    if (isVerbose):
        print('  --> removing files associated with : ', header)
    headerInfo = read_header(header, isVerbose = isVerbose)
    iBlock = 1
    for iAltBlock in range(headerInfo['nBlocksAlt']):
        for iLatBlock in range(headerInfo['nBlocksLat']):
            for iLonBlock in range(headerInfo['nBlocksLon']):
                cBlock = '.b%4.4i' % iBlock
                file = header.replace('.header', cBlock)
                command = 'rm -f ' + file
                if (isVerbose):
                    print('     --> running command : ' + command)
                os.system(command)
                iBlock += 1
    command = 'rm -f ' + header
    if (isVerbose):
        print('     --> running command : ' + command)
    os.system(command)
    return

# ----------------------------------------------------------------------------
# 
# ----------------------------------------------------------------------------

def process_one_file(header, isVerbose = False):

    fileOut = header.replace('.header','.bin')

    print('Processing : ', header)
    headerInfo = read_header(header, isVerbose = isVerbose)

    fpin = open(header, 'r')
    lines = fpin.readlines()
    fpin.close()

    data = {}
    # Move all of the header info into the data dictionary:
    for key in headerInfo.keys():
        if (isVerbose):
            print('   -> found key in header : ', key)
            print('      value : ', headerInfo[key])
        data[key] = headerInfo[key]

    data['nLonsTotal'] = data['nBlocksLon'] * data['nLons'] + 2 * data['nGhostsLon']
    data['nLatsTotal'] = data['nBlocksLat'] * data['nLats'] + 2 * data['nGhostsLat']
    data['nAltsTotal'] = data['nBlocksAlt'] * data['nAlts'] + 2 * data['nGhostsAlt']
    
    for iVar in np.arange(data['nVars']):
        v = data['vars'][iVar]
        # Make the data array in reverse order so that it will write out correctly:
        data[v] = np.zeros((data['nAltsTotal'], data['nLatsTotal'], data['nLonsTotal']))

    # Now read in all of the binary files and put the data where it belong:
    iBlock = 1
    for iAltBlock in range(data['nBlocksAlt']):
        for iLatBlock in range(data['nBlocksLat']):
            for iLonBlock in range(data['nBlocksLon']):
                cBlock = '.b%4.4i' % iBlock
                file = header.replace('.header', cBlock)
                oneFile = read_one_binary_file(file, \
                                               data['nVars'], \
                                               data['nLons'] + 2 * data['nGhostsLon'], \
                                               data['nLats'] + 2 * data['nGhostsLat'], \
                                               data['nAlts'] + 2 * data['nGhostsAlt'],
                                               isVerbose = isVerbose)

                iWholeLonStart, iLonStart, iWholeLonEnd, iLonEnd = \
                    determine_bounds(iLonBlock, \
                                     data['nBlocksLon'], data['nLons'], data['nGhostsLon'])
                iWholeLatStart, iLatStart, iWholeLatEnd, iLatEnd = \
                    determine_bounds(iLatBlock, \
                                     data['nBlocksLat'], data['nLats'], data['nGhostsLat'])
                iWholeAltStart, iAltStart, iWholeAltEnd, iAltEnd = \
                    determine_bounds(iAltBlock, \
                                     data['nBlocksAlt'], data['nAlts'], data['nGhostsAlt'])
                iLonOff = iLonStart - iWholeLonStart
                iLatOff = iLatStart - iWholeLatStart
                iAltOff = iAltStart - iWholeAltStart
                for iVar in np.arange(data["nVars"]):
                    v = data["vars"][iVar]
                    data[v][iWholeAltStart:iWholeAltEnd, \
                            iWholeLatStart:iWholeLatEnd, \
                            iWholeLonStart:iWholeLonEnd] = \
                                oneFile[iVar][iAltStart:iAltEnd, \
                                              iLatStart:iLatEnd, \
                                              iLonStart:iLonEnd]
                iBlock += 1
        
    write_gitm_file(fileOut, data, isVerbose = isVerbose)
    return data
    

def post_process_gitm(dir, doRemove, isVerbose = False):

    # Get into the correct directory for processing
    processDir = dir
    if (not os.path.exists(processDir)):
        print('Processing directory doesnt exist: ', processDir)
        print('Make sure to put in a valid -dir= ')
        exit()
    # Move into the UA/data directory
    if (isVerbose):
        print('Moving into processing directory : ', processDir)
    cwd = os.getcwd()
    os.chdir(processDir)

    # Get list of header files to process:
    headers = sorted(glob('*.header'))
    if (isVerbose):
        print('Found %i files to process' % len(headers))
    
    # Process each single header:
    for head in headers:
        # Attempt to process file:
        if (isVerbose):
            print(f' --> Processing {head}')

        data = process_one_file(head, isVerbose = isVerbose)
        if (doRemove):
            remove_files(head, isVerbose = isVerbose)

    if (isVerbose):
        print('Moving into old directory : ', cwd)
    os.chdir(cwd)
            
    return True
    
# ----------------------------------------------------------------------------
# Main Code
# ----------------------------------------------------------------------------

if __name__ == '__main__':

    # Get the input arguments
    args = get_args_pgitm()

    doPrint = args.v
    doRemove = not args.norm
    
    if (args.header != 'none'):
        data = process_one_file(args.header, isVerbose = doPrint)
        if (doRemove):
            remove_files(args.header, isVerbose = doPrint)
        exit()

    post_process_gitm(args.dir, doRemove, isVerbose = doPrint)

    
    
