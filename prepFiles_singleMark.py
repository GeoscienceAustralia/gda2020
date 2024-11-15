#!/usr/bin/env python3

"""This script prepares all the RINEX files and antenna information files for
automatic processing
"""

import re
import os
import datetime

import sys


# Compile the regular expressions
#p1 = re.compile('^\w{8}\.\d{2}o$', re.I)

# Create the directory for the individual RINEX antenna information files
os.mkdir('rinexantls')

# Create single-mark clusters from RinexAntLs.txt
for line in open('RinexAntLs.txt'):

    # Form cluster name
    rnxantls = line[0:8] + '_' + line[9:11]  + '_ls'
    
    # Make cluster file
    fout = open(rnxantls, 'w')
    fout.write(line)
    fout.close()

# Move antenna information files to rinexantls/
os.system('mv *_ls rinexantls')

