#!/usr/bin/env python3

"""
This is a modified version of createBLs.py v2.10. It is designed to work on
SNXEPO.SNX.NONCON.AUS, the SINEX file of APREF stations that have less than
2 years of data. It is in GDA2020 and has Type B uncertainties added

NAME:
    createBLsAPREF.py
PURPOSE:
    Form baselines from the stations in the APREF SINEX file and create DynaML
    formatted files.
EXPLANATION:
    The code takes SNXEPO.SNX.NONCON.AUS and returns a pair of DynaML
    formatted station measurement files for input into DynAdjust
USAGE:
    createBLsAPREF.py YYYYMMDD
INPUT:
    SNXEPO.SNX.NONCON.AUS
OUTPUT:
    aprefYYYYMMDD_stn.xml
    aprefYYYYMMDD_msr.xml
HISTORY:
    0.01    2019-03-25  Craig Harrison
            - Written
"""
import sys
import os
import datetime
import re
from glob import glob
from numpy import matrix, zeros, copy


# Set some variables
refFrame = 'GDA2020'
epoch = '01.01.2020'
inputFile = 'SNXEPO.SNX.NONCON.AUS'
aprefVer = sys.argv[1]

# Get the starting epoch of the solution epoch data from inputFile
epochs = {}
go = False
for line in open(inputFile):
    if '-SOLUTION/EPOCHS' in line:
        break
    if go and '*Code' not in line:
        stat = line[1:5] + line[11:13].strip()
        epochs[stat] = line[16:18] + line[19:22]
        if int(epochs[stat]) < 94000:
            epochs[stat] = '20' + epochs[stat]
        else:
            epochs[stat] = '19' + epochs[stat]
    if '+SOLUTION/EPOCHS' in line:
        go = True

# Read in the discontinuity information for renaming purposes
yearDOY = {}
go = False
for line in open('soln.snx'):
    if '-SOLUTION/DISCONTINUITY' in line:
        break
    if go and line[0:1] != '*':
        cols = line.rstrip().split()
        if cols[6] == 'P' and cols[4] != cols[5]:
            segment = cols[0] + cols[2]
            try: 
                epochs[segment]
                yearDOY[segment] = line[16:18] + line[19:22]
                if int(yearDOY[segment]) < 94000:
                    yearDOY[segment] = '20' + yearDOY[segment]
                else:
                    yearDOY[segment] = '19' + yearDOY[segment]
                if yearDOY[segment] == '2000000':
                    yearDOY[segment] = '1900001'
            except KeyError:
                pass
    if '+SOLUTION/DISCONTINUITY' in line:
        go = True
go = False

# Determine which stations have discontinuities
discontStats = set()
for key, value in yearDOY.items():
    stat = key[0:4]
    discontStats.add(stat)

# Open the output files
stn = open('apref' + aprefVer + '_stn.xml', 'w')
msr = open('apref' + aprefVer + '_msr.xml', 'w')

# Write headers
stn.write('<?xml version="1.0"?>\n')
stn.write('<DnaXmlFormat type="Station File" referenceframe="GDA2020" ' +
            'epoch="01.01.2020" ' +
            'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ' +
            'xsi:noNamespaceSchemaLocation="DynaML.xsd">\n')

msr.write('<?xml version="1.0"?>\n')
msr.write('<DnaXmlFormat type="Measurement File" referenceframe="GDA2020" ' +
            'epoch="01.01.2020" ' +
            'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ' +
            'xsi:noNamespaceSchemaLocation="DynaML.xsd">\n')

# Open the SINEX file and read in all lines
snxFile = open(inputFile)
lines = snxFile.readlines()

# Create lists to hold the site ID, station coordinate estimate, and the VCV 
# matrix lines
estimateLines = []
matrixLines = []
goE = 0
goM = 0
for line in lines:
    if re.match('\+SOLUTION/ESTIMATE', line):
        goE = 1
    if re.match('\+SOLUTION/MATRIX_ESTIMATE', line):
        goM = 1
    if goE:
        if not re.match('\+|\*|\-', line):
            estimateLines.append(line)
    if goM:
        if not re.match('\+|\*|\-', line):
            (a, b) = divmod(int(line[:6].lstrip()) - 1, 3)
            if a % 2 == 0:
                matrixLines.append(line)
    if re.match('\-SOLUTION/ESTIMATE', line):
        goE = 0
    if re.match('\-SOLUTION/MATRIX_ESTIMATE', line):
        goM = 0

# Create a list of dictionaries to hold the station names and their coordinates
data = []
estimateLines.reverse()
while estimateLines:
    col = estimateLines.pop().rstrip().split()
    if col[1][0:3] != 'VEL':
        source = {}
        if col[2].upper() in discontStats:
            segment = col[2] + col[4]
            stat = col[2].upper() + '_' + yearDOY[segment] 
        else:
            stat = col[2].upper()
        source['site'] = stat
        source['x'] = float(col[8])
        col = estimateLines.pop().rstrip().split()
        source['y'] = float(col[8])
        col = estimateLines.pop().rstrip().split()
        source['z'] = float(col[8])
        data.append(source)

# Create the variance-covariance matrix. In the SINEX file it is given as a
# upper triangular matrix, but it is immediately converted to a lower
# trianguler matrix
vcvXVL = matrix(zeros((6*len(data), 6*len(data))))
for line in matrixLines:
    cols = line.rstrip().split()
    row = int(cols[0]) - 1
    col = int(cols[1]) - 1
    for i in range(2, len(cols)):
        vcvXVL[col+i-2, row] = float(cols[i])
vcvL = matrix(zeros((3*len(data), 3*len(data))))
for i in range(len(vcvXVL)):
    (a, b) = divmod(i, 6)
    row = i - 3 * a
    for j in range(len(vcvXVL)):
        (a, b) = divmod(j, 6)
        col = j - 3 * a
        if col < len(vcvL) and row < len(vcvL):
            vcvL[row, col] = float(vcvXVL[i, j])
vcvU = copy(vcvL.transpose())
vcv = vcvU + vcvL

# Create the design matrix
desMatrix = matrix(zeros((3*(len(data)-1), 3*len(data))))
for i in range(len(data)-1):
    desMatrix[3*i, 0] = -1
    desMatrix[3*i+1, 1] = -1
    desMatrix[3*i+2, 2] = -1
    desMatrix[3*i, 3*(i+1)] = 1
    desMatrix[3*i+1, 3*(i+1)+1] = 1
    desMatrix[3*i+2, 3*(i+1)+2] = 1
                                                                
# Create the matrix of observed antenna positions
coords = matrix(zeros((3*len(data), 1)))
for i in range(len(data)):
    coords[3*i, 0] = data[i]['x']
    coords[3*i+1, 0] = data[i]['y']
    coords[3*i+2, 0] = data[i]['z']

# Calculate the deltas and the corresponding VCV matrix
deltas = desMatrix * coords
delVCV = desMatrix * vcv * desMatrix.transpose()

# Loop over the sites and write the station data to the output XML file
for i in range(len(data)):
    stn.write('\t<DnaStation>\n')
    stn.write('\t\t<Name>%s</Name>\n'%(data[i]['site']))
    stn.write('\t\t<Constraints>FFF</Constraints>\n')
    stn.write('\t\t<Type>XYZ</Type>\n')
    stn.write('\t\t<StationCoord>\n')
    stn.write('\t\t\t<Name>%s</Name>\n'%(data[i]['site']))
    stn.write('\t\t\t<XAxis>%20.14e</XAxis>\n'%(data[i]['x']))
    stn.write('\t\t\t<YAxis>%20.14e</YAxis>\n'%(data[i]['y']))
    stn.write('\t\t\t<Height>%20.14e</Height>\n'%(data[i]['z']))
    stn.write('\t\t\t<HemisphereZone></HemisphereZone>\n')
    stn.write('\t\t</StationCoord>\n')
    stn.write('\t\t<Description></Description>\n')
    stn.write('\t</DnaStation>\n')

# Create an array of the stations minus the core station
coreStat = data[0]['site']
nonCoreStns = []
for i in range(1, len(data)):
    nonCoreStns.append(data[i]['site'])
numCovar = len(nonCoreStns) - 1

# Loop over the non-core stations and write the measurement data to the output
# XML file
msr.write('\t<!--Type X GNSS baseline cluster (full correlations)-->\n')
msr.write('\t<DnaMeasurement>\n')
msr.write('\t\t<Type>X</Type>\n')
msr.write('\t\t<Ignore/>\n')
msr.write('\t\t<ReferenceFrame>%s</ReferenceFrame>\n'%(refFrame))
msr.write('\t\t<Epoch>%s</Epoch>\n'%(epoch))
msr.write('\t\t<Vscale>1.000</Vscale>\n')
msr.write('\t\t<Pscale>1.000</Pscale>\n')
msr.write('\t\t<Lscale>1.000</Lscale>\n')
msr.write('\t\t<Hscale>1.000</Hscale>\n')
msr.write('\t\t<Total>%s</Total>\n'%(len(nonCoreStns)))
for i in range(len(nonCoreStns)):
    msr.write('\t\t<First>%s</First>\n'%(coreStat))
    msr.write('\t\t<Second>%s</Second>\n'%(nonCoreStns[i]))
    msr.write('\t\t<GPSBaseline>\n')
    msr.write('\t\t\t<X>%20.14e</X>\n'%(deltas[3*i, 0]))
    msr.write('\t\t\t<Y>%20.14e</Y>\n'%(deltas[3*i+1, 0]))
    msr.write('\t\t\t<Z>%20.14e</Z>\n'%(deltas[3*i+2, 0]))
    msr.write('\t\t\t<SigmaXX>%20.14e</SigmaXX>\n'%(delVCV[3*i, 3*i]))
    msr.write('\t\t\t<SigmaXY>%20.14e</SigmaXY>\n'%(delVCV[3*i+1, 3*i]))
    msr.write('\t\t\t<SigmaXZ>%20.14e</SigmaXZ>\n'%(delVCV[3*i+2, 3*i]))
    msr.write('\t\t\t<SigmaYY>%20.14e</SigmaYY>\n'%(delVCV[3*i+1, 3*i+1]))
    msr.write('\t\t\t<SigmaYZ>%20.14e</SigmaYZ>\n'%(delVCV[3*i+2, 3*i+1]))
    msr.write('\t\t\t<SigmaZZ>%20.14e</SigmaZZ>\n'%(delVCV[3*i+2, 3*i+2]))
    for j in range(numCovar):
        msr.write('\t\t\t<GPSCovariance>\n')
        msr.write('\t\t\t\t<m11>%20.14e</m11>\n'% \
            (delVCV[3*(i+1)+3*j, 3*i]))
        msr.write('\t\t\t\t<m12>%20.14e</m12>\n'% \
            (delVCV[3*(i+1)+3*j+1, 3*i]))
        msr.write('\t\t\t\t<m13>%20.14e</m13>\n'% \
            (delVCV[3*(i+1)+3*j+2, 3*i]))
        msr.write('\t\t\t\t<m21>%20.14e</m21>\n'% \
            (delVCV[3*(i+1)+3*j, 3*i+1]))
        msr.write('\t\t\t\t<m22>%20.14e</m22>\n'% \
            (delVCV[3*(i+1)+3*j+1, 3*i+1]))
        msr.write('\t\t\t\t<m23>%20.14e</m23>\n'% \
            (delVCV[3*(i+1)+3*j+2, 3*i+1]))
        msr.write('\t\t\t\t<m31>%20.14e</m31>\n'% \
            (delVCV[3*(i+1)+3*j, 3*i+2]))
        msr.write('\t\t\t\t<m32>%20.14e</m32>\n'% \
            (delVCV[3*(i+1)+3*j+1, 3*i+2]))
        msr.write('\t\t\t\t<m33>%20.14e</m33>\n'% \
            (delVCV[3*(i+1)+3*j+2, 3*i+2]))
        msr.write('\t\t\t</GPSCovariance>\n')
    numCovar -= 1
    msr.write('\t\t</GPSBaseline>\n')
msr.write('\t\t<Source>%s</Source>\n'%os.path.basename(inputFile))
msr.write('\t</DnaMeasurement>\n')
stn.write('</DnaXmlFormat>\n')
msr.write('</DnaXmlFormat>\n')
