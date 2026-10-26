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

import subprocess
import geodepy.gnss
from ngca_path import NGCA_DIR
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()

parser.add_argument("-j", "--juris", required=True, help="Jurisdiction to process")

args = parser.parse_args()


# Copy over gda2020.dat
GDA2020_DAT_PATH = NGCA_DIR / 'auxiliary_files' / 'gda2020.dat'
SINEX_DIR = NGCA_DIR / args.juris.lower() / 'sinexFiles'

# Create a list of APREF stations to exclude
exclude = ['MAC1', 'CEDU', 'CA19', 'HOBA', 'CCPL', 'IGWD'];

# Loop over the input SINEX files
for f in SINEX_DIR.glob('*.SNX.ITRF2014'):

    # Use subprocess to run rdsinex and GMT and create a list of stations 
    # outside GDA2020
    stns = set()
    rdsnx = subprocess.Popen(["rdsinex","-P",f], stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL)
    gmt = subprocess.Popen(["gmt", "select", "-Ef","-F" + str(GDA2020_DAT_PATH), "-If"],
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
    ofile = f.with_name(f.name + ".AUS")
    Path("output.snx").rename(ofile)

print('*** EXCLUDING ***')
print(exclude)
