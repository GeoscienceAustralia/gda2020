#!/usr/bin/env python3

# This script removes stations that are outside the extent of GDA2020 from an 
# NGCA SINEX file
#
# Requirements: 
#	- GMT program select
#	- GA program rdsinex
#
# Usage: exciseStationsNGCA.pl <files>
#
# Note: wildcards are accepted

import os
import glob
import subprocess
import geodepy.gnss


# Copy over gda2020.dat
os.system('cp ~/apref/workDir/gda2020.dat ../sinexFiles')

# Create a list of APREF stations to exclude
exclude = ['MAC1', 'CEDU', 'CA19', 'HOBA', 'CCPL', 'IGWG'];

# Loop over the input SINEX files
os.chdir('../sinexFiles')
for f in glob.glob('*.SNX'):

    # Use subprocess to run rdsinex and GMT and create a list of stations 
    # outside GDA2020
    stns = set()
    rdsnx = subprocess.Popen(["rdsinex","-P",f], stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL)
    gmt = subprocess.Popen(["gmt", "select", "-Ef","-Fgda2020.dat", "-If"],
                           stdin=rdsnx.stdout, stdout=subprocess.PIPE)
    all_details_b = gmt.stdout.read()
    all_details = all_details_b.decode()
    stns_details = all_details.split('\n')
    for stn_details in stns_details:
        if stn_details != '':
            details = stn_details.split()
            stns.add(details[6])
    
    # Check to see if any of the excluded stations are in the SINEX file and add
    # them to the stations to be removed
    rdsnx = subprocess.Popen(["rdsinex","-P",f], stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL)
    all_details_b = rdsnx.stdout.read()
    all_details = all_details_b.decode()
    stns_details = all_details.split('\n')
    for stn_details in stns_details:
        if stn_details != '':
            details = stn_details.split()
            if details[6] in exclude:
                stns.add(details[6])

    # Remove stations and rename output file
    geodepy.gnss.remove_stns_sinex(f, stns)
    ofile = f + '.AUS'
    os.rename('output.snx', ofile)

# Tidy up
os.remove('gda2020.dat')
print('*** EXCLUDING ***')
print(exclude)
