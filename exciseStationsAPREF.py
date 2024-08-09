#!/usr/bin/env python3

# This script removes stations that are outside the extent of GDA2020 from
# the APREF cumulative solution 
#
# Requirements: 
#	- GMT program select
#	- GA program rdsinex
#
# Usage: exciseStationsAPREF.py

import os
import glob
import sys
import subprocess
import geodepy.gnss 


# Get the path to the home directory and move to the working directory
home = os.path.expanduser("~")
os.chdir(home + '/apref/workDir/')

# Check that the file containing the GDA2020 polygons is in the work directory
if not os.path.exists('gda2020.dat'):
    sys.exit('The file of GDA2020 polygons is not in the working directory')

# Get the input file name
for f in glob.glob('AUS0OPSSNX_*_00U_SOL.SNX'):
    ifile = f
try:
    ifile
except NameError:
    sys.exit('There is no APREF SINEX file in the working directory')

# Use subprocess to run rdsinex and GMT and create a list of stations inside
# GDA2020
stns = set()
rdsnx = subprocess.Popen(["rdsinex","-P",ifile], stdout=subprocess.PIPE,
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

# Remove stations and rename output file
geodepy.gnss.remove_stns_sinex(ifile, stns)
ofile = ifile + '.AUS'
os.rename('output.snx', ofile)
