#!/usr/bin/env python3

""" This script generates the file typeb.dat that specifies the Type B
    uncertainties that should be applied to each station in the given
    SINEX file.

    RVS stations: 3mm, 3mm, 6mm (enu)
    Non-RVS stations: 6mm, 6mm, 12mm (enu)

    If an RVS station has multiple segments then only the last is treated as an
    RVS station, the other segments are treated like non-RVS stations
"""    

import sys


# Read in the list of RVS stations
rvsStations = []
for line in open('RVS_GDA2020.txt'):
    rvsStations.append(line[0:4])
rvsStations = rvsStations[1:]

# Loop over the lines in the SINEX file and get the solution lines
solLines = []
go = False
for line in open(sys.argv[1]):
    line = line.strip()
    if line[0:15] == '*Code PT SOLN T':
        go = True
    elif line[0:16] == '-SOLUTION/EPOCHS':
        break
    elif go:
        solLines.append(line)

# Loop over the solution lines and determine the last segment for each station
lastSegment = {}
for solLine in solLines:
    try:
        lastSegment[solLine[0:4]]
        if lastSegment[solLine[0:4]] < solLine[10:12]:
                lastSegment[solLine[0:4]] = solLine[10:12]
    except KeyError:
        lastSegment[solLine[0:4]] = solLine[10:12]

# Open the output file and loop over the solution lines and write out the Type
# B information
fout = open('typeb.dat', 'w')
for solLine in solLines:
    stat = solLine[0:4]
    pt = solLine[10:12]
    if stat in rvsStations and pt == lastSegment[stat]:
        fout.write('0.003 0.003 0.006 ' + stat + ' A 50127M001 ' + pt + '\n')
    else:
        fout.write('0.006 0.006 0.012 ' + stat + ' A 50127M001 ' + pt + '\n')
