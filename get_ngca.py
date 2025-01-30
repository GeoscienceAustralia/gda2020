#!/usr/bin/env python3

""" This script downloads the NGCA from the AWS S3 bucket gda2020-ngca
"""

# Import modules
import argparse
import re
import pysftp
import sys
import os
import shutil
import glob
from ftplib import FTP
from datetime import datetime
from pathlib import Path


# Create an ArgumentParser object
parser = argparse.ArgumentParser(
    description='Download the NGCA from the AWS S3 bucket gda2020-ngca.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Add the arguments and parse the command line
parser.add_argument('-j', type=str, metavar='JURIS', default='all',
                    help='Download only the data from JURIS')
parser.add_argument('-v', action='version', version='%(prog)s v1.0')
args = parser.parse_args()
if args.j == 'all':
    juris = ['act', 'tas', 'sa', 'vic', 'nt', 'wa', 'qld', 'nsw']
else:
    juris = [args.j]

# Compile regular expressions
p1 = re.compile('\w{8}\.\d{2}o$', re.I)

# Move to the NGCA directory
print('* Moving to the NGCA directory')
os.chdir(f'{Path.home()}/ngca/')

# Get today's date in the format YYYYMMDD
today = str(datetime.today())
archiveDate = today[0:10].replace('-','')

# Loop over the jurisdictions
for jur in juris:
    print('* Processing ' + jur.upper())

    # Create the archive directory and move to it
    os.mkdir(jur + '/' + archiveDate)
    os.chdir(jur + '/' + archiveDate)

    # Download the files
    print('* Downloading files')
    os.system('aws s3 cp s3://gda2020-ngca/ngca/' + jur.lower() + \
            ' . --quiet --recursive --include "*"')

    # Deleting old data
    print('* Deleting old SINEX files and RINEX antenna information files')
    os.chdir('../')
    for cluster in glob.glob('rinexantls/*'):
        os.remove(cluster)
    for snxFile in glob.glob('sinexFiles/*') :
        os.remove(snxFile)

    # Move back up to main directory
    os.chdir('../')

