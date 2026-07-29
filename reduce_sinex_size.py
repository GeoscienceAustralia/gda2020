#!/usr/bin/env python3

# This script reduces the size of the SINEX file by removing the
# velocity parameters (which are not needed for this workflow).
# If velocity parameters are needed, then change the function call
# to remove_matrixzeros_sinex() insetad of remove_velocity_sinex() below.
# This will remove zero lines only from +SOLUTION/MATRIX_ESTIMATE block,
# Whilst reatining the velocity parameters.
#
# NOTE: this file does overwrite the input file, but a copy is made,
#       with the extension tag "VEL"
#
#
# Usage: reduce_sinex_size.py

import os
import glob
import geodepy.gnss
import sys
from shutil import copy2

# Get the path to the home directory and move to the working directory
home = os.path.expanduser("~")
os.chdir(home + '/apref/workDir/')

# Get the input file name
files = glob.glob('apref*.snx')

if not files:
    sys.exit('There is no APREF SINEX file by that name in the working directory')

if len(files)> 1:
    sys.exit('There are multiple APREF SINEX files in the working directory') 

ifile = files[0]  

# Make copy of original file
copy2(ifile, f"{ifile}.VEL")

# Remove velocity parameters and rename output file
try:
    geodepy.gnss.remove_velocity_sinex(ifile)
except: 
    raise NameError("Velocities could not be removed from sinex file")

ofile = ifile
os.rename('output.snx', ofile)
