#!/usr/bin/env python3

# renameless 20251111 - added additonal debug prints, deactivated all except print when removing SINEX due to zero user stations

import subprocess, logging
from math import sqrt
from pathlib import Path
import shutil
import geodepy.gnss
from ngca_path import NGCA_DIR, NGCA_APREF_DIR
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-j", "--juris", required=True, help="Jurisdiction to process")

args = parser.parse_args()

# Directory
JURIS_DIR = NGCA_DIR / args.juris.lower()
SINEX_DIR = NGCA_DIR / args.juris.lower() / 'sinexFiles'
RINEXANTLS_DIR = NGCA_DIR / args.juris.lower() / 'rinexantls'

# Create a list of the APREF stations that are used as constraints, that is,
# those that have more than 2 years of data
go = 0
constraints = []

# Find the APREF solution file
apref_sols = list(NGCA_APREF_DIR.glob("apref*.snx"))

# If not exactly one solution file, raise error
if len(apref_sols) == 0:
    raise FileNotFoundError(f"No apref*.snx file found in {apref_sols}")

if len(apref_sols) > 1:
    raise RuntimeError(f"Expected one apref.snx file, found {len(apref_sols)}: {apref_sols}")  

aprefSol = apref_sols[0]  

# Create a list of APREF stations that are constraints (In solution file)
f = open(aprefSol)
for line in f:
    if line[0:8] == '-SITE/ID':
        break 
    if go and line[0:5] != '*CODE':
        constraints.append(line[1:5])
    if line[0:8] == '+SITE/ID':
        go = 1 
f.close()

# Find most recent NGCA directory
dirs = sorted(JURIS_DIR.glob('20??????'))
lates_dir = dirs[-1]

# Creating log file in most recent NGCA directory
logFile = lates_dir / 'selectRef.log'
tempLogFile = lates_dir / 'temp.log'
logging.basicConfig(filename=logFile, level=logging.INFO, filemode='w')

# Loop over the SINEX files
for snx in SINEX_DIR.glob('*.AUS'):
    print(snx)

    # Create name of new sinex file
    ngcaVer = snx.with_name(snx.name.replace('AUS', args.juris.upper() + '.NGCA'))

    # Find name of related rinexantls file
    rinexantls = RINEXANTLS_DIR / snx.name.replace('.SNX.ITRF2014.AUS', '_ls')
    
    sub = []

    # Adds name of stations in rinexantls to sub list. This gives a list of the NGCA stations for the job
    for line in open(rinexantls):
        sub.append(line[0:4])
    
    # Get station coordinates from sinex file
    proc = subprocess.Popen(["rdsinex", "-X", snx], stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    snxLines = proc.communicate()
    snxLines=snxLines[0].decode("utf-8")
    snxLines = snxLines.split('\n')
    
    line = {}
    ref = []
    refJuris = []
    
    # Loop over each station in sinex file and add to dictionary with the station name as the key 
    for snxLine in snxLines:

        if snxLine != '':
            cols = snxLine.split()
            line[cols[0]] = snxLine

            if cols[0] not in sub:
                ref.append(cols[0]) # If station is not in the submitted list, add to ref list

    sumX = 0
    sumY = 0
    sumZ = 0
    # Find the number of submitted stations
    numStn = len(sub)

    stnsRemoved = []
    
    # Loop over each submitted station
    for stn in sub:
        try:
            # Get the X, Y, Z coordinates of the station from the sinex file and add to the sum of X, Y, Z
            cols = line[stn].split()

            sumX = sumX + float(cols[4])
            sumY = sumY + float(cols[5])
            sumZ = sumZ + float(cols[6])

        # If submitted station is not in the sinex file, log the error and add to removed stations list
        except KeyError:
            logging.info('Station ' + stn + ' in ' + snx +
                         ' was submitted but failed processing')

            numStn -= 1
            stnsRemoved.append(stn)

    # If there are still submitted station in sinex file, calculate the centre of submitted stations
    if numStn != 0:
        cenX = sumX / len(sub)
        cenY = sumY / len(sub)
        cenZ = sumZ / len(sub)
    
        statDist = []
        
        # Loop through the reference stations and calculate distance from centre of submitted stations
        for stn in ref:
            if stn in constraints:
                data = []
                data.append(stn)
                cols = line[stn].split()
                dist = sqrt( (cenX - float(cols[4]))**2 + 
                             (cenY - float(cols[5]))**2 +
                             (cenZ - float(cols[6]))**2 )
                data.append(dist)
                statDist.append(data)

        # Sort by distance
        statDist = sorted(statDist, key=lambda x: x[1])

        # Remove the two cloest reference stations so they will stay in the sinex file.
        ref.remove(statDist[0][0])
        ref.remove(statDist[1][0])

        # Remove the remaining reference stations from the sinex file
        geodepy.gnss.remove_stns_sinex(snx, ref)
        Path('output.snx').rename(ngcaVer)
        
        # If any stations were not found in the sinex file, remove them from the rinexantls file
        if stnsRemoved:
            lines = [line for line in rinexantls.read_text().splitlines(True) if line[0:4] not in stnsRemoved]
            rinexantls.write_text(''.join(lines))
        
    # If there are no submitted stations in the sinex file, delete the sinex and rinexantls files and log        
    else:
        fullSNX = snx.with_suffix('')

        fullSNX.unlink()
        snx.unlink()
        rinexantls.unlink()

        print(f"Adding log - Sinex deleted, no stations left: {snx}")
        logging.info(snx + 
                ' was deleted as it had no submitted stations left.')

logging.shutdown()
shutil.copyfile(logFile, tempLogFile) # I am not sure yet why this tmp log file exists

if logFile.stat().st_size == 0:
    logFile.unlink()
