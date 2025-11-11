#!/usr/bin/env python3

# renameless 20251111 - added additonal debug prints, deactivated all except print when removing SINEX due to zero user stations

import sys, subprocess, os, logging
from glob import glob
from math import sqrt
from pathlib import Path
import shutil
import geodepy.gnss

# Move to the SINEX directory
d = os.getcwd() 
j = d.split('/')[-2]
os.chdir(f'{Path.home()}/ngca/' + j + '/sinexFiles')

# Create a list of the APREF stations that are used as constraints, that is,
# those that have more than 2 years of data
go = 0
constraints = []
for f in glob(f'{Path.home()}/apref/apref*.snx'):
    aprefSol = f
f = open(aprefSol)
for line in f:
    if line[0:8] == '-SITE/ID':
        break 
    if go and line[0:5] != '*CODE':
        constraints.append(line[1:5])
    if line[0:8] == '+SITE/ID':
        go = 1 
f.close()

# Set up logging
for f in glob('../20??????'):
    pass
logFile = f + '/selectRef.log'
tempLogFile = f + '/temp.log'
if os.path.isfile(logFile):
    os.remove(logFile)
logging.basicConfig(filename = logFile, level = logging.INFO)

# Loop over the SINEX files
for snx in glob('*.AUS'):
    print(snx)
    # os.system('printf "looping on snx: "' + snx)
    # os.system('printf "\n"')
    ngcaVer = snx.replace('AUS', j.upper() + '.NGCA')
    rinexantls = snx.replace('.SNX.AUS', '_ls')
    rinexantls = '../rinexantls/' + rinexantls
    sub = []

    for line in open(rinexantls):
        sub.append(line[0:4])
        # os.system('printf "looping on rinexantls line : "' + line)
        # os.system('printf "\n"')
        debug_printarray = []
        #trying this - want all 4 elements
        debug_printarray = ','.join(line[0:4])
        # os.system('printf "line elements 0-4 added to sub: "' + debug_printarray)
        # os.system('printf "\n"')
    proc = subprocess.Popen(["rdsinex", "-X", snx], stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    snxLines = proc.communicate()
    snxLines=snxLines[0].decode("utf-8")
    snxLines = snxLines.split('\n')
    line = {}
    ref = []
    refJuris = []
    for snxLine in snxLines:
        # os.system('printf "looping on snxLine : "' + snxLine)
        # os.system('printf "\n"')
        if snxLine != '':
            cols = snxLine.split()
            line[cols[0]] = snxLine
            # os.system('printf "Line cols element 0 : "' + line[cols[0]])
            # os.system('printf "\n"')
            if cols[0] not in sub:
                ref.append(cols[0])
                # os.system('printf "cols element 0 not in sub, appending to ref : "' + cols[0])
                # os.system('printf "\n"')
    sumX = 0
    sumY = 0
    sumZ = 0
    numStn = len(sub)
    # os.system('printf "numStn value : "' + str(numStn))
    # os.system('printf "\n"')
    stnsRemoved = []
    for stn in sub:
        # os.system('printf "looping on stn in sub - stn value : "' + str(stn))
        # os.system('printf "\n"')
        try:
            line[stn]

            # os.system('printf "trying line of stn, value: "' + str(line[stn]))
            # os.system('printf "\n"')
            cols = line[stn].split()
            debug_printarray = []
            debug_printarray = ','.join(cols)
            # os.system('printf "trying cols over line of stn, value: "' + str(debug_printarray))
            # os.system('printf "\n"')
            sumX = sumX + float(cols[4])
            # os.system('printf "trying sumX value : "' + str(sumX))
            # os.system('printf "\n"')
            sumY = sumY + float(cols[5])
            # os.system('printf "trying sumY value : "' + str(sumY))
            # os.system('printf "\n"')
            sumZ = sumZ + float(cols[6])
            # os.system('printf "trying sumZ value : "' + str(sumZ))
            # os.system('printf "\n"')
        except KeyError:
            logging.info('Station ' + stn + ' in ' + snx +
                         ' was submitted but failed processing')
            # os.system('printf "Exception : Item added to log - Station stn in snx submitted but failed processing "')
            # os.system('printf "\n"')
            # os.system('printf "Exception : stn "' + stn)
            # os.system('printf "\n"')
            # os.system('printf "Exception : snx "' + stn)
            # os.system('printf "\n"')

            numStn -= 1
            stnsRemoved.append(stn)
            # os.system('printf "stnRemoved item appended : "' + stn)
            # os.system('printf "\n"')
    if numStn != 0:
        # os.system('printf "numSTN is nonzero"')
        # os.system('printf "\n"')
        cenX = sumX / len(sub)
        cenY = sumY / len(sub)
        cenZ = sumZ / len(sub)
        # os.system('printf "centre X value : "' + str(cenX))
        # os.system('printf "\n"')
        # os.system('printf "centre Y value : "' + str(cenY))
        # os.system('printf "\n"')
        # os.system('printf "centre Z value : "' + str(cenY))
        # os.system('printf "\n"')
    
        statDist = []
        for stn in ref:
            if stn in constraints:
                data = []
                data.append(stn)
                # os.system('printf "looping on constraint stn: "' + stn)
                # os.system('printf "\n"')
                cols = line[stn].split()
                dist = sqrt( (cenX - float(cols[4]))**2 + 
                             (cenY - float(cols[5]))**2 +
                             (cenZ - float(cols[6]))**2 )
                data.append(dist)
                # os.system('printf "calculated dist: "' + str(dist))
                # os.system('printf "\n"')
                statDist.append(data)

        os.system('printf "sorted statDist: "' + str(debug_printarray))
        os.system('printf "\n"')
        statDist = sorted(statDist, key=lambda x: x[1])

        # os.system('printf "removing ref element: "' + str(statDist[0][0]))
        # os.system('printf "\n"')
        ref.remove(statDist[0][0])
        # os.system('printf "removing ref element: "' + str(statDist[1][0]))
        # os.system('printf "\n"')
        ref.remove(statDist[1][0])
        # os.system('printf "Calling geodepy remove stations sinex "')
        # os.system('printf "\n"')
        # os.system('printf "Geodepy passed snx: "' + snx)
        # os.system('printf "\n"')
        debug_printarray = []
        debug_printarray = ','.join(ref)
        # os.system('printf "Geodepy passed ref: "' + str(debug_printarray))
        # os.system('printf "\n"')
        geodepy.gnss.remove_stns_sinex(snx, ref)
        ngcaSnx = snx.replace('AUS', j.upper() + '.NGCA')
        os.rename('output.snx', ngcaSnx)
        if stnsRemoved:
            os.rename(rinexantls, 'tmp')
            f = open(rinexantls, 'w')
            for line in open('tmp'):
                if line[0:4] not in stnsRemoved:
                    # os.system('printf "Element not in stnsRemoved: "' + str(line[0:4]))
                    # os.system('printf "\n"')
                    f.write(line)
                    # os.system('printf "line added to f "' + line)
                    # os.system('printf "\n"')
            f.close()
            os.remove('tmp')
    else:
        # os.system('printf "else case - numstn is zero "')
        # os.system('printf "\n"')
        fullSNX = snx.replace('.AUS', '')
        # os.system('printf "Deleting fullSNX: "' + fullSNX)
        # os.system('printf "\n"')
        os.remove(fullSNX)
        # os.system('printf "Deleting snx: "' + snx)
        # os.system('printf "\n"')
        os.remove(snx)
        # os.system('printf "Deleting rinexantls: "' + snx)
        # os.system('printf "\n"')
        os.remove(rinexantls)
        os.system('printf "Adding log - Sinex deleted, no stations left: "' + snx)
        os.system('printf "\n"')
        logging.info(snx + 
                ' was deleted as it had no submitted stations left.')
logging.shutdown()
shutil.copyfile(logFile, tempLogFile)
if os.stat(logFile).st_size == 0:
    os.remove(logFile)
