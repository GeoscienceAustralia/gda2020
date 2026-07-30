#!/usr/bin/env python3
  
"""This script updates the translation tables in transTables/ with those in the
S3 bucket
"""

import  os
import glob
from pathlib import Path

# Define the jurisdictions
jurisdictions = ['nsw', 'act', 'vic', 'tas', 'sa', 'wa', 'nt', 'qld']

# Move to the translation tables directory and delete the contents 
os.chdir(f'{Path.home()}/transTables')
for f in os.listdir('.'):
    os.remove(f)

# Copy the files directory from the S3 bucket
os.system('aws s3 cp s3://gda2020-ngca/files/ . --quiet --recursive --exclude "*" --include "*Trans*"')

# Loop over the csv files
for f in os.listdir('.'):
    
# Determine the jurisdiction
    if f[:3].lower() in jurisdictions:
        jurisdiction = f[:3].lower()
    elif f[:2].lower() in jurisdictions:
        jurisdiction = f[:2].lower()
    else:
        print('File ' + f + ' has an unknown jurisdiction')
        unlink(f)

    # Determine the date
    date = f.split('.')[0]
    date = date[-8:]

    # Create the new file name 
    f_new = jurisdiction + 'TransTable' + date + '.csv'
    os.rename(f, f_new)
