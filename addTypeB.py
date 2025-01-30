#!/usr/bin/env python3

"""This script adds Type B uncertainties to those given in the .apu file.
"""

from sys import exit
from glob import glob
from geodepy.convert import dec2hp, hp2dec
from numpy import matrix
from geodepy.statistics import vcv_cart2local, error_ellipse, circ_hz_pu
from math import sqrt
from pathlib import Path

# Determine the files to use 
apuFiles = glob('*.apu')
if (len(apuFiles) == 1):
    apuFile = apuFiles[0]
elif (len(apuFiles) == 0):
    exit('\nThere is no apu file to work on\n')
else:
    print('\nThere are multiple apu files:')
    i = 0
    for apuFile in apuFiles:
        i += 1
        print('\t' + str(i) + '\t' + apuFile)
    fileNum = input('Type the number of the file you want to check: ')
    if int(fileNum) < 1 or int(fileNum) > len(apuFiles):
        exit('Invalid response. Select a number between 1 and ' +
            str(len(apuFiles)))
    apuFile = apuFiles[int(fileNum) - 1]

# Set the Type B uncertainties
rvsE = 0.003
rvsN = 0.003
rvsU = 0.006
nonRvsE = 0.006
nonRvsN = 0.006
nonRvsU = 0.012

# Read in the RVS stations
rvsStations = []
for line in open(f'{Path.home()}/rvs_stations.dat'):
    rvsStations.append(line.rstrip())

# Open output file
fout = open(apuFile + '.typeB', 'w')

# Read in the apu file
apuLines = []
i = 0
with open(apuFile) as f:
    for line in f:
        if line[:9] == 'Station  ':
            j = i + 2
        apuLines.append(line.rstrip())
        i += 1

# Print out the header info
for line in apuLines[:j]:
    fout.write(line + '\n')

# Loop over the .apu file and read in the uncertainty info
stations = []
hpLat = {}
hpLon = {}
lat = {}
lon = {}
hPU = {}
vPU = {}
semiMajor = {}
semiMinor = {}
orient = {}
xLine = {}
xVar = {}
xyCoVar = {}
xzCoVar = {}
yLine = {}
yVar = {}
yzCoVar = {}
zLine = {}
zVar = {}
for line in apuLines[j:]:
    cols = line.split()
    numCols = len(cols)
    if numCols == 2:
        yLine[station] = line
        yVar[station] = float(line[131:150].strip())
        yzCoVar[station] = float(line[150:].strip())
    elif numCols == 1:
        zLine[station] = line
        zVar[station] = float(line[150:].strip())
    else:
        station = line[:20].rstrip()
        stations.append(station)
        hpLat[station] = float(line[23:36])
        hpLon[station] = float(line[38:51])
        lat[station] = hp2dec(hpLat[station])
        lon[station] = hp2dec(hpLon[station])
        hPU[station] = float(line[51:62].strip())
        vPU[station] = float(line[62:73].strip())
        semiMajor[station] = float(line[73:86].strip())
        semiMinor[station] = float(line[86:99].strip())
        orient[station] = float(line[99:112].strip())
        xLine[station] = line[112:]
        xVar[station] = float(line[112:131].strip())
        xyCoVar[station] = float(line[131:150].strip())
        xzCoVar[station] = float(line[150:].strip())

# Create the full Cartesian VCV from the upper triangular
vcv_cart = {}
for stat in stations:
    vcv_cart[stat] = matrix([[xVar[stat], xyCoVar[stat], xzCoVar[stat]],
                             [xyCoVar[stat], yVar[stat], yzCoVar[stat]],
                             [xzCoVar[stat], yzCoVar[stat], zVar[stat]]
                            ])

# Loop over all the stations
for stat in stations:
    # Transform the XYZ VCV to ENU
    vcv_local = vcv_cart2local(vcv_cart[stat], lat[stat], lon[stat])
    
    # Add the Type B uncertainty
    if stat in rvsStations:
        vcv_local[0, 0] += rvsE**2
        vcv_local[1, 1] += rvsN**2
        vcv_local[2, 2] += rvsU**2
    else:
        vcv_local[0, 0] += nonRvsE**2
        vcv_local[1, 1] += nonRvsN**2
        vcv_local[2, 2] += nonRvsU**2
    
    # Calculate the semi-major axis, semi-minor axis and orientation, and
    # convert the orientation from deciaml degrees to HP notation
    a, b, orientation = error_ellipse(vcv_local)
    orientation = dec2hp(orientation)

    # Calculate the PUs
    hz_pu = circ_hz_pu(a, b)
    vt_pu = 1.96 * sqrt(vcv_local[2, 2])

    # Output the uncertainties
    line = '{:20}{:>16.9f}{:>15.9f}{:11.4f}{:11.4f}{:13.4f}{:13.4f}{:13.4f}'. \
            format(stat, hpLat[stat], hpLon[stat], hz_pu, vt_pu, a, b,
                    orientation)
    line += xLine[stat]
    fout.write(line + '\n')
    fout.write(yLine[stat] + '\n')
    fout.write(zLine[stat] + '\n')

#AD: Check this: print('What about the xyz and adj files?')
