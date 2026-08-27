#!/usr/bin/env python3

""" This script downloads the NGCA from the AWS S3 bucket gda2020-ngca
"""

# Import modules
import argparse
import re
from datetime import datetime
from ngca_path import NGCA_DIR
import subprocess


# Create an ArgumentParser object
parser = argparse.ArgumentParser(
    description='Download the NGCA from the AWS S3 bucket gda2020-ngca.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Add the arguments and parse the command line
parser.add_argument('-j', type=str, metavar='JURIS', default='all', nargs="+",
                    help='Download only the data from JURIS')
parser.add_argument('-v', action='version', version='%(prog)s v1.0')
args = parser.parse_args()
if args.j == 'all':
    juris = ['act', 'tas', 'sa', 'vic', 'nt', 'wa', 'qld', 'nsw']
else:
    juris = args.j

# Compile regular expressions
p1 = re.compile(r'\w{8}\.\d{2}o$', re.I)

# Get today's date in the format YYYYMMDD
today = str(datetime.today())
archiveDate = today[0:10].replace('-','')
print(archiveDate)

# Loop over the jurisdictions
for jur in juris:
    print('* Processing ' + jur.upper())

    # Create the archive directory
    archive_dir = NGCA_DIR / jur / archiveDate
    archive_dir.mkdir(parents=True, exist_ok=True)

    # Download the files
    print('* Downloading files')

    subprocess.run(["aws", "s3", "cp", 
                    f"s3://gda2020-ngca/ngca/{jur.lower()}/" ,
                    archive_dir, "--quiet", "--recursive"], check=True
                   )

    # Deleting old data
    print('* Deleting old SINEX files and RINEX antenna information files')
        
    for folder in ["rinexantls", "sinexFiles"]:
        ngca_fold_path = NGCA_DIR / jur / folder

        for item in ngca_fold_path.iterdir():
            item.unlink()
