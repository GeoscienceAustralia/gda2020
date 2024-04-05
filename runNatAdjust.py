#!/usr/bin/env python3

"""This script runs DynAdjust to do the national adjustment and various pieces
of QA
"""

import argparse
from datetime import date
import os
from glob import glob


# Set the parameter defaults and intialise some variables
min_inner = 3200
max_block = 3200
geoidFile = 'AUSGeoid2020_20180201.gsb'

# Create an ArgumentParser object
parser = argparse.ArgumentParser(
    description='Run DynAdjust to do the NADJ and various pieces of QA.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Add the arguments and parse the command line
parser.add_argument('-a', '--apu', action='store_true',
                    help='Create an APU file that contains the covariances')
parser.add_argument('-g', '--gen', action='store_true',
                    help='Generate bash script and exit')
parser.add_argument('-i', '--iter', type=int, 
                    help='Set the maximum number of iterations')
parser.add_argument('-m', '--mode', choices=['simult', 'phased', 'staged'
                                             'multi'], 
                    help='Select the adjutsment mode')
parser.add_argument('-n', '--nohup', action='store_true',
                    help='Run the adjustment in the background')
parser.add_argument('-qa', action='store_true',
                    help='Perform a 1-iteration adjustment for QA purposes')
parser.add_argument('-s', '--search', choices=['dup', 'near'],
                    help='Perform a station search')
parser.add_argument('-t', '--thresh', type=float,
                    help='The threshold for an iteration to converge')
args = parser.parse_args()

# Get today's date and create the epoch and network name
epoch = date.today()
epoch = epoch.strftime('%Y%m%d')
network = 'gda2020_' + epoch
if args.search == 'dup':
    network += '.dup'
elif args.search == 'near':
    network += '.near'

# Create the bash script
fout = open(network, 'w')
fout.write('#!/bin/bash\n')

# If doing a station search
if args.search:
    fout.write('dnaimport -n ' + network + ' ')
    files = glob('aprefRename_???.xml')
    files.sort(reverse=True)
    for file in files:
        fout.write(file + ' ')
    for file in glob('apref????????_stn.xml'):
        fout.write(file + ' ')
    for file in glob('apref????????_msr.xml'):
        fout.write(file + ' ')
    files = glob('stn/*.xml')
    files.sort()
    for file in files:
        fout.write(file + ' ')
    files = glob('msr/*.xml')
    files.sort()
    for file in files:
        fout.write(file + ' ')
    if args.search == 'dup':
        fout.write('--prefer-single-x-as-g --flag-unused --ignore-similar-msr ')
        fout.write('--remove-ignored-msr -r GDA2020\n')
    elif args.search == 'near':
        fout.write('--search-nearby-stn --prefer-single-x-as-g --flag-unused ')
        fout.write('--ignore-similar-msr --remove-ignored-msr -r GDA2020\n')
    fout.close()
    if not args.gen: 
        cmd = 'bash ' + network
        os.system(cmd)
    for file_pattern in ['*.aml', '*.asl', '*.bms', '*.bst', '*.dbid',
                         '*.dms', '*.dnaproj', '*.imp', '*.map']: 
        for file in glob(file_pattern):
            os.remove(file)

# If running an adjustment
else:

    # Create the dnaimport command
    fout.write('dnaimport -n ' + network + ' ')
    if os.path.isfile('apriori.xml'):
        fout.write('apriori.xml ')
    for file in glob('apref????????.snx'):
        fout.write(file + ' ')
    for file in glob('apref????????_stn.xml'):
        fout.write(file + ' ')
    for file in glob('apref????????_msr.xml'):
        fout.write(file + ' ')
    files = glob('stn/*.xml')
    files.sort()
    for file in files:
        fout.write(file + ' ')
    files = glob('msr/*.xml')
    files.sort()
    for file in files:
        fout.write(file + ' ')
    fout.write('--prefer-single-x-as-g --flag-unused -r GDA2020 ')
    fout.write('--remove-ignored-msr --discontinuity-file ')
    for file in glob('disconts????????.snx'):
        fout.write(file)
    if args.qa:
        fout.write('\n')
    else: 
        fout.write(' --output-msr-to-stn --export-xml-files\n')

    # Transform the data to GDA2020 before making the geoid corrections
    fout.write('dnareftran -n ' + network + ' -r GDA2020\n')

    # Create the dnageoid command
    fout.write('dnageoid -n ' + network + ' -g ' + geoidFile + ' --convert ')
    fout.write('--verbose-level 1\n')

    # Create the dnasegment command
    fout.write('dnasegment ' + network + ' --min ' + str(min_inner) +
               ' --max ' + str(max_block) + '\n')

    # Create the dnaadjust command
    fout.write('dnaadjust ' + network + ' ')
    if args.mode == 'simult':
        mode = '--simultaneous-adjustment '
    elif args.mode == 'phased':
        mode = '--phased-adjustment '
    elif args.mode == 'staged':
        mode = '--staged --create-stage-files '
    elif args.mode == 'multi':
        mode = '--multi-thread '
    
    # If doing QA set the number of iterations to one
    if args.qa:
        args.iter = 1
    if args.iter:
        fout.write('--max-iterations ' + str(args.iter) + ' ')
    if args.thresh:
        fout.write('--iteration-threshold ' + str(args.thresh) + ' ')
    if args.apu:
        fout.write('--output-all-covariances ')

    # Set the output. Only output the adjusted measurements if doing QA
    if not args.qa:
        fout.write('--output-pos-uncertainty --output-corrections-file ')
        fout.write('--angular-msr-type 1 --stn-coord-types "ENzPLHhXYZ" ')
        fout.write('--export-xml-stn-file ')
    fout.write('--output-adj-msr\n')

    # Rename the output stn and msr files
    fout.write('mv ' + network + 'stn.xml ' + network + '_stn.xml\n')
    fout.write('mv ' + network + 'msr.xml ' + network + '_msr.xml\n')

    # Delete the large unnecessary files and move the rest to the adjustments
    # folder
    fout.write('rm -f *.mtx *.bms *.bst\n')
    fout.write('mkdir adjustments/gda2020_' + epoch + '\n')
    fout.write('cp *.adj.stn.xml apriori.xml\n')
    fout.write('mv *' + epoch + '* adjustments/gda2020_' + epoch + '\n')
    fout.close()
    if not args.gen:
        if args.nohup:
            cmd = 'nohup bash ' + network + '&'
        else:
            cmd = 'bash ' + network
        os.system(cmd)
    if args.qa:
        for file_pattern in ['*.aml', '*.asl', '*.bms', '*.bst', '*.dbid',
                             '*.dnaproj', '*.dst', '*.imp', '*.map', '*.mtx',
                             '*.seg', '*.xyz']:
            for file in glob(file_pattern):
                os.remove(file)
        if not args.gen:
            cmd = 'vi ' + network + '.phased-stage.adj'
            os.system(cmd)
