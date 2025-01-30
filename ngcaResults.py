#!/usr/bin/env python3

import sys
import os
import shutil
import glob
from pathlib import Path

d = os.getcwd()
jur = d.split('/')[-2]
os.chdir('../')
dirs = glob.glob('20*')
dirs.sort()
date = dirs[-1]
os.chdir(d)
campaign = jur.upper() + '_NGCA_' + date
os.chdir(f'{Path.home()}/ngca/sent')
try:
    os.mkdir(campaign)
except OSError:
    shutil.rmtree(campaign)
    os.mkdir(campaign)
os.chdir(campaign)
os.mkdir('dynaML')
os.chdir('dynaML')
os.system('cp ~/ngca/' + jur + '/baselines/*.xml .')
os.system('rm -f *GDA2020*')
os.system('mv *NGCA* ../')
os.chdir('../')
os.mkdir('adjustments')
os.chdir('adjustments')
os.system('cp ~/ngca/' + jur + '/baselines/*.adj .')
os.chdir('../')
os.mkdir('feedback')
os.chdir('feedback')
for item in ['sigma0.dat', 'checkNgca.log']:
    os.system('cp ~/ngca/' + jur + '/baselines/' + item + ' .')
for item in ['log_4_submit.txt', 'selectRef.log', 'verifySub.log']:
    target = '../../../' + jur + '/' + date + '/' + item
    sys_cmd = 'cp ../../../' + jur + '/' + date + '/' + item + ' .'
    if os.path.isfile(target):
        os.system(sys_cmd)
os.chdir('../')
os.mkdir('other')
os.chdir('other')
os.system('cp /home/fedora/ngca/' + jur + '/' + date + '/nameChanges.dat ./' )
os.system('cp /home/fedora/transTables/' + jur + 'TransTable*.csv ./' )
os.chdir('../')
sys_cmd = 'cp -r ../../' + jur + '/sinexFiles ./snx'
os.system(sys_cmd)
os.system('rm snx/*.SNX')
sys_cmd = 'cp -r ../../' + jur + '/rinexantls ./clusters'
os.system(sys_cmd)
os.chdir('../')
zip_cmd = 'zip -r ' + campaign + '.zip ' + campaign
os.system(zip_cmd)
os.system('rm -rf ' + campaign)
cmd = 'aws s3 rm s3://gda2020-ngca/ngcaResults/ --quiet --recursive --exclude "*" --include "*' + jur.upper() + '*"'
print(cmd)
os.system(cmd)
cmd = 'aws s3 cp ' + campaign  + '.zip s3://gda2020-ngca/ngcaResults/ --quiet'
print(cmd)
os.system(cmd)
