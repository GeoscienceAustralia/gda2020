#!/usr/bin/env python3

#changed checkngca version to v2 - accounting for ENz info in v2.
import sys
import os
import re

p1 = re.compile(r'\d{8}$')
d = os.getcwd()
jur = d.split('/')[-2]
os.chdir('../')
dirs = os.listdir('.')
dirs.sort()
for dir in dirs:
    if p1.match(dir):
        date = dir
os.chdir(d)
os.system('checkNgca.py 20 20 40 > checkNgca.log')
os.system('mkdir all')
os.system('mv *all*.xml all')
network = jur.upper() + '_NGCA_' + date
os.system('dnaimport -n ' + network + 
        ' *stn.xml *msr.xml -r ITRF2014 --export-xml')
os.system('dnareftran ' + network + ' -r GDA2020 --export-xml')
os.system('mv ' + network + 'stn.xml ' + network + '_EoO_stn.xml')
os.system('mv ' + network + 'msr.xml ' + network + '_EoO_msr.xml')
os.system('mv ' + network + '.GDA2020.1.1.2020stn.xml ' + network +
    '_stn.xml')
os.system('mv ' + network + '.GDA2020.1.1.2020msr.xml ' + network +
    '_msr.xml')
os.system('rm ' + network + '.*')
