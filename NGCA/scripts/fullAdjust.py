#!/usr/bin/env python3

# This script automates the fully constrained adjustment performed after NGCA scaling

# renameless 20251111 - include E N coordinates via stn-coord-types and retain .xyz files for use with DynaDiff

import sys
import os
import glob
import re
from pathlib import Path

p1 = re.compile(r'\d{8}$')
d = os.getcwd()
jur = d.split('/')[-2]
jur = jur.upper()
os.chdir('../')
dirs = os.listdir('.')
dirs.sort()
for dir in dirs:
    if p1.match(dir):
        ngcaVer = dir
os.chdir(d)
aprefPath = glob.glob(f'{Path.home()}/apref/apref*.snx')
aprefFile = aprefPath[-1].split('/')[-1]
discontsFile = aprefFile.replace('apref', 'disconts')

os.system('dnaimport -n ' + jur + '_RENAME ' + jur + '_NGCA_' +
        ngcaVer + '_stn.xml ' + jur + '_NGCA_' + ngcaVer +
        '_msr.xml -r GDA2020 --flag-unused-stations --remove-ignored-msr '
        '--export-xml')
os.system('cp ' + jur + '_RENAMEmsr.xml ' + jur + '_GDA2020_' + ngcaVer +
        '_msr.xml')
os.system('dnaimport -n ' + jur + '_GDA2020 ' + aprefFile +
        ' ' + jur + '_RENAMEstn.xml ' + jur +
        '_RENAMEmsr.xml -r GDA2020 --discontinuity-file ' + discontsFile)
os.system('dnasegment ' + jur + '_GDA2020 --min 800 --max 800')
os.system('dnaadjust ' + jur + '_GDA2020 --phased --max-iter 30 ' 
        '--stn-coord-types "ENzPLHhXYZ" '
        '--output-adj-msr --sort-adj-msr-field 7 --export-xml-stn '
        '--export-xml-msr')
os.system('dnaimport -n ' + jur + '_GDA2020_' + ngcaVer + 
        ' ' + jur + '_GDA2020.phased.adj.stn.xml ' + jur + '_GDA2020_' +
        ngcaVer + '_msr.xml -r GDA2020 --flag-unused-stations --export-xml')
os.system('rm ' + jur + '_RENAME*')
os.system('cp ' + jur + '_GDA2020.phased.adj temp')
os.system('cp ' + jur + '_GDA2020.phased.xyz temp2')
os.system('rm ' + jur + '_GDA2020[.-]*')
os.system('mv temp ' + jur + '_GDA2020.phased.adj')
os.system('mv temp2 ' + jur + '_GDA2020.phased.xyz')
os.system('rm ' + jur + '_GDA2020_' + ngcaVer + '.*')
os.system('mv ' + jur + '_GDA2020_' + ngcaVer + 'stn.xml ' + jur + '_GDA2020_' +
        ngcaVer + '.adj.xml')
os.system('rm ' + jur + '_GDA2020_' + ngcaVer + 'msr.xml')
