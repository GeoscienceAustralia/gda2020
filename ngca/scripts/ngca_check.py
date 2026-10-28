#!/usr/bin/env python3

"""Identify problems arising from NGCA processing

Go through the .adj files produced from getSigma0.pl and pull out stations
that have large reduiuals or large uncertainties on their coordinates

checkNgca.py <residual> <h uncert> <v uncert>

residual: baseline compnents with residuals larger than this will be flagged
            (in mm)
h uncert: stations that have horizontal uncertainties larger than this will
            be flagged (in mm)
v uncert: stations that have vertical uncertainties larger than this will be
            flagged (in mm)

Revised the column/array numbers for added E/N/z information in adjustment scripts.            
"""

from glob import glob
import argparse
from pathlib import Path

# Create an ArgumentParser object
parser = argparse.ArgumentParser(
    description='Identify problems arising from NGCA processing.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Add the arguments and parse the command line
parser.add_argument('residual', metavar='R', type=int,
                    help='flag residuals larger than this (in mm)')
parser.add_argument('horizontal', metavar='H', type=int,
                    help='flag hori. uncertainties larger than this (in mm)')
parser.add_argument('vertical', metavar='V', type=int,
                    help='flag vert. uncertainties larger than this (in mm)')
parser.add_argument('--directory', type=Path, default=Path.cwd(),
                    help='direcotry containing adjustment files')
args = parser.parse_args()

# Loop over the .adj files
for file in args.directory.glob('*.adj'):
    badR = False
    badU = False
    outLinesR = []
    outLinesU = []

    # Loop over the lines in the .adj file
    with file.open() as f:
        for line in f:

            # Check the baseline residuals
            if line[0:2] == 'X ':
                resid = line[105:117]
                resid = '{:.1f}'.format(float(resid.strip()) * 1000)
                
                if abs(float(resid)) > args.residual:
                    badR = True
                    stat1 = line[2:22].strip()
                    stat2 = line[22:42].strip()
                    comp = line[65:66]
                    outLinesR.append(
                        f"{stat1:20s}{stat2:20s} {comp}   {resid}"
                    )

            # Check the coordinate uncertainties
            if line[20:23] == 'FFF':
                east = line[158:170]
                east = '{:.1f}'.format(float(east.strip()) * 1000)
                
                if float(east) > args.horizontal:
                    badU = True
                    stat = line[0:20].strip()
                    outLinesU.append(
                        f"{stat:40s} E   {east}"
                    )

                north = line[170:180]
                north = '{:.1f}'.format(float(north.strip()) * 1000)
                
                if float(north) > args.horizontal:
                    badU = True
                    stat = line[0:20].strip()
                    outLinesU.append(
                        f"{stat:40s} N   {north}"
                    )
                
                up = line[180:190]
                up = '{:.1f}'.format(float(up.strip()) * 1000)
                
                if float(up) > args.vertical:
                    badU = True
                    stat = line[0:20].strip()
                    outLinesU.append(
                        f"{stat:40s} U   {up}"
                    )

    # Output the results, if there are any
    if badR or badU:
        name = file.name.replace(".simult.adj", "")
        print(f"### {name} ###")
        
        if badR:
            for line in outLinesR:
                print(line)
        if badU:
            for line in outLinesU:
                print(line)

