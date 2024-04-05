#!/usr/bin/env python3

"""This script copies the files directory from the AWS S3 bucket, gda2020-ngca,
making sure the file names correctly follow the convention, places them in the
correct directory on NADJ
"""

import re
import os


# Set up regular expressions
p1 = re.compile('.csv$', re.I)
p2 = re.compile('.ignore$', re.I)
p3 = re.compile('.near$', re.I)
p4 = re.compile('.renaming$', re.I)

# Define the jurisdictions
jurisdictions = ['nsw', 'act', 'vic', 'tas', 'sa', 'wa', 'nt', 'qld']

# Create a temporary directory and move there
os.mkdir('aws_files')
os.chdir('aws_files')

# Copy the files directory from the S3 bucket
os.system('aws s3 cp s3://gda2020-ngca/files/ . --quiet --recursive --include "*"')

# Loop over the files
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
    if p1.search(f):
        f_new = jurisdiction + 'TransTable' + date + '.csv'
    elif p2.search(f):
        f_new = jurisdiction + '_' + date + '.ignore'
    elif p3.search(f):
        f_new = jurisdiction + '_' + date + '.near'
    elif p4.search(f):
        f_new = jurisdiction + '_' + date + '.renaming'
    else:
        print('File ' + f + ' is of an unknown type')
        unlink(f)
    os.rename(f, f_new)

# Remove all the old files
os.system('rm /home/fedora/transTables/*.csv')
os.system('rm /home/fedora/renaming/*.ignore')
os.system('rm /home/fedora/nearStns/*.near')
os.system('rm /home/fedora/renaming/*.renaming')

# Move the files to their correct directory
os.system ('mv *.csv /home/fedora/transTables/')
os.system ('mv *.renaming /home/fedora/renaming/')
os.system ('mv *.ignore /home/fedora/renaming/')
os.system ('mv *.near /home/fedora/nearStns/')

# Remove empty directory
os.chdir('/home/fedora/')
os.system('rm -r aws_files/')
