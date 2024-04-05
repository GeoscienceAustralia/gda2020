#!/usr/bin/env python3

"""Update the Vscale in the NGCA baseline files based on the values in
sigma0.dat
"""

import glob
import os
import shutil


# Read in the sigma0 values
sigma = {}
for line in open('sigma0.dat'):
    (clus, x, y, sigma0) = line.split()
    clus = clus.replace('_all.simult.adj:Rigorous', '')
    if float(sigma0) < 1:
        sigma[clus] = '1.000'
    else:
        sigma[clus] = sigma0

# Loop over the msr files
for file in glob.glob('*msr.xml'):
    os.rename(file, 'tmp')
    clus = file.replace('_msr.xml', '')
    clus = clus.replace('_all', '')
    try:
        sigma[clus]
        fOut = open(file, 'w')
        for line in open('tmp'):
            line = line.replace('<Vscale>1.000', '<Vscale>' + str(sigma[clus]))
            fOut.write(line)
    except KeyError:
        print('No sigma0 for clus ' + clus)
        shutil.move( 'tmp', file)

# Clean up
os.remove('tmp')
