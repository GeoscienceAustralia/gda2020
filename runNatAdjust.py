"""This script runs DynAdjust to do the national adjustment and various pieces
of QA
"""

import argparse
from datetime import date
import os
from glob import glob
import sys

def abspaths(pattern, **kwargs):
    return [os.path.abspath(f) for f in sorted(glob(pattern), **kwargs)]

# Set the parameter defaults and initialise some variables
min_inner = 3200
max_block = 3200
geoidFile = 'AUSGeoid2020_20180201.gsb'
DYNADJUST_DIR="/home/ubuntu/dynadjust-linux-static"

apref_snx = abspaths('APREF*/apref*.snx', reverse=True)
apref_stn = abspaths('APREF*/apref????????_stn.xml')
apref_msr = abspaths('APREF*/apref????????_msr.xml')
disconts = abspaths('APREF*/disconts*.snx')
stn_files = abspaths('stn/*.xml')
msr_files = abspaths('msr/*.xml')

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
parser.add_argument('-m', '--mode', choices=['simult', 'phased',
                    'staged', 'multi'], default='staged',
                    help='Select the adjustment mode')
parser.add_argument('-qa', action='store_true',
                    help='Perform a 1-iteration adjustment for QA purposes')
parser.add_argument('-s', '--search', choices=['dup', 'near', 'both'],
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
bash_file = network + '.sh'

# Create the bash script
fout = open(bash_file, 'w')

fout.write('#!/bin/bash\n')
fout.write(f'export PATH="{DYNADJUST_DIR}:${{PATH}}"\n')
fout.write(f'export LD_LIBRARY_PATH="{DYNADJUST_DIR}:${{LD_LIBRARY_PATH}}"\n')

# If doing a station search
if args.search:
    import_parts = ['dnaimport -n ' + network]
    for f in abspaths('APREF*/aprefRename_???.xml', reverse=True):
        import_parts.append(f)
    for f in apref_stn:
        import_parts.append(f)
    for f in apref_msr:
        import_parts.append(f)
    for f in stn_files:
        import_parts.append(f)
    for f in msr_files:
        import_parts.append(f)
    if args.search == 'dup':
        import_parts.append('--prefer-single-x-as-g --flag-unused '
                            '--ignore-similar-msr --remove-ignored-msr -r GDA2020')
        import_parts.append('\nperl checkDST.pl')
    elif args.search == 'near':
        import_parts.append('--search-nearby-stn --prefer-single-x-as-g '
                            '--flag-unused --ignore-similar-msr '
                            '--remove-ignored-msr -r GDA2020')
        import_parts.append('\nperl filterNearStns.pl')
    elif args.search == 'both':
        import_parts.append('--search-nearby-stn --prefer-single-x-as-g '
                            '--flag-unused --ignore-similar-msr '
                            '--remove-ignored-msr -r GDA2020')
        import_parts.append('\nperl checkDST.pl')
        import_parts.append('\nperl filterNearStns.pl')
    commands = ' '.join(import_parts)

    fout.write(commands + '\n')

    if not args.gen:
        os.chmod(bash_file, 0o755)
        bash_file = "./" + bash_file
        print(bash_file)
        fout.close()
        os.system(bash_file)
    for file_pattern in ['*.aml', '*.asl', '*.bms', '*.bst', '*.dbid',
                         '*.dms', '*.dnaproj', '*.imp', '*.map', 'notUsedIgnore.dat']: 
        for file in glob(file_pattern):
            os.remove(file)

# If running an adjustment
else:

    # Create the import command
    fout.write('dnaimport -n ' + network + ' ')
    if os.path.isfile('apriori.xml'):
        fout.write('apriori.xml ')
    for file in glob('APREF*/apref????????.snx'):
        fout.write(file + ' ')
    for file in glob('APREF*/apref????????_stn.xml'):
        fout.write(file + ' ')
    for file in glob('APREF*/apref????????_msr.xml'):
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
    for file in glob('APREF*/disconts????????.snx'):
        fout.write(file)
    if args.qa:
        fout.write('\n')
    else: 
        fout.write(' --output-msr-to-stn --export-xml-files\n')

    # Transform the data to GDA2020 before making the geoid corrections
    fout.write('dnareftran -n ' + network + ' -r GDA2020\n')

    # Create the geoid command
    fout.write('dnageoid -n ' + network + ' -g ' + geoidFile + ' --convert ')
    fout.write('--verbose-level 1\n')

    # Create the segment command
    fout.write('dnasegment ' + network + ' --min ' + str(min_inner) +
               ' --max ' + str(max_block) + '\n')

    # Create the adjust command
    fout.write('dnaadjust ' + network + ' ')
    if args.mode == 'simult':
        mode = '--simultaneous-adjustment '
    elif args.mode == 'phased':
        mode = '--phased-adjustment '
    elif args.mode == 'staged':
        mode = '--staged --create-stage-files --stage-path /mnt/data --max-threads 64 '
    elif args.mode == 'multi':
        mode = '--multi-thread '
    fout.write(mode)
    
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
    fout.write('--output-adj-msr --export-xml-stn-file\n')

    # Rename the output stn and msr files
    if not args.qa:
        fout.write('mv ' + network + 'stn.xml ' + network + '_stn.xml\n')
        fout.write('mv ' + network + 'msr.xml ' + network + '_msr.xml\n')

        # Delete the large unnecessary files and move the rest to the
        # adjustments folder
        fout.write('rm /mnt/data/* *.mtx *.bms *.bst\n')
        fout.write('mkdir adjustments/gda2020_' + epoch + '\n')
        fout.write('cp *.adj.stn.xml apriori.xml\n')
        fout.write('mv *' + epoch + '* adjustments/gda2020_' + epoch + '\n')
    fout.close()
    if not args.gen:
        os.chmod(bash_file, 0o755)
        bash_file = "./" + bash_file
        os.system(bash_file)
    if args.qa:
        for file_pattern in ['*.aml', '*.asl', '*.bat', '*.bms', '*.bst',
                             '*.dbid', '*.dnaproj', '*.dst', '*.imp', '*.map',
                             '*.mtx', '*.rft', '*.seg', '*.xyz']:
            for file in glob(file_pattern):
                os.remove(file)
