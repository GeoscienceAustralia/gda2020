#!/usr/bin/env python3

import os, sys
import glob
import shutil


# Create a directory for the solutions to be copied to
if os.path.exists('solutions'):
    shutil.rmtree('solutions')
os.mkdir('solutions')

# Loop over all jobs
with open('log_4_submit.txt') as f:
    for line in f:

        # Get job number and cluster name
        col = line.split('--')
        job = col[0].split(' ')[0]
        cluster = col[2].split('_ls')[0]
        cluster = cluster.lstrip()
        
        # Copy over the SINEX file
        cmd = 'mv ../SNX/' + job + '.SNX solutions/' + cluster + '.SNX'
        try:
            os.system(cmd)
        except FileNotFoundError:
            print('No solution could be found for job ' + job +
                  ' (' + cluster + ')')

os.system('cp *_ls/*_ls solutions')
