#!/usr/bin/env python3

import sys, subprocess, os, logging
from glob import glob
from math import sqrt
import shutil
import geodepy.gnss


# Move to the SINEX directory
d = os.getcwd() 
j = d.split('/')[-2]
os.chdir('/home/fedora/ngca/' + j + '/sinexFiles')

# Create a list of the APREF stations that are used as constraints, that is,
# those that have more than 2 years of data
go = 0
constraints = []
for f in glob('/home/fedora/apref/apref*.snx'):
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
    ngcaVer = snx.replace('AUS', j.upper() + '.NGCA')
    rinexantls = snx.replace('.SNX.AUS', '_ls')
    rinexantls = '../rinexantls/' + rinexantls
    sub = []
    for line in open(rinexantls):
        sub.append(line[0:4])
    proc = subprocess.Popen(["rdsinex", "-X", snx], stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    snxLines = proc.communicate()
    snxLines=snxLines[0].decode("utf-8")
    snxLines = snxLines.split('\n')
    line = {}
    ref = []
    refJuris = []
    for snxLine in snxLines:
        if snxLine != '':
            cols = snxLine.split()
            line[cols[0]] = snxLine
            if cols[0] not in sub:
                ref.append(cols[0])
    sumX = 0
    sumY = 0
    sumZ = 0
    numStn = len(sub)
    stnsRemoved = []
    for stn in sub:
        try:
            line[stn]
            cols = line[stn].split()
            sumX = sumX + float(cols[4])
            sumY = sumY + float(cols[5])
            sumZ = sumZ + float(cols[6])
        except KeyError:
            logging.info('Station ' + stn + ' in ' + snx +
                         ' was submitted but failed processing')
            numStn -= 1
            stnsRemoved.append(stn)
    if numStn != 0:
        cenX = sumX / len(sub)
        cenY = sumY / len(sub)
        cenZ = sumZ / len(sub)
    
        statDist = []
        for stn in ref:
            if stn in constraints:
                data = []
                data.append(stn)
                cols = line[stn].split()
                dist = sqrt( (cenX - float(cols[4]))**2 + 
                             (cenY - float(cols[5]))**2 +
                             (cenZ - float(cols[6]))**2 )
                data.append(dist)
                statDist.append(data)
        statDist = sorted(statDist, key=lambda x: x[1])
        ref.remove(statDist[0][0])
        ref.remove(statDist[1][0])
        geodepy.gnss.remove_stns_sinex(snx, ref)
        ngcaSnx = snx.replace('AUS', j.upper() + '.NGCA')
        os.rename('output.snx', ngcaSnx)
        if stnsRemoved:
            os.rename(rinexantls, 'tmp')
            f = open(rinexantls, 'w')
            for line in open('tmp'):
                if line[0:4] not in stnsRemoved:
                    f.write(line)
            f.close()
            os.remove('tmp')
    else:
        fullSNX = snx.replace('.AUS', '')
        os.remove(fullSNX)
        os.remove(snx)
        os.remove(rinexantls)
        logging.info(snx + 
                ' was deleted as it had no submitted stations left.')
logging.shutdown()
shutil.copyfile(logFile, tempLogFile)
if os.stat(logFile).st_size == 0:
    os.remove(logFile)
